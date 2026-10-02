#!/usr/bin/env python3
"""Deterministic self-contained ZIP; packaging performs zero science evaluations.

The complete checkpoint payload, including recorded outputs and freeze receipt,
is protected. Only the literal package bookkeeping names and Python caches are
excluded. A fresh ZIP replay is made outside this payload after construction.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import stat
import zipfile

from ledger_integrity import (CHECKPOINT_NAME, EXCLUDED, ZIP_NAME, ZIP_SHA_NAME,
                              digest, file_sha, fresh_external, payload_files, read_json, require,
                              safe_file, safe_relative, verify_frozen, write_json)


def construct_zip(root, names, destination):
    with zipfile.ZipFile(destination, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9,
                         allowZip64=True) as archive:
        for name in sorted(names):
            info = zipfile.ZipInfo(CHECKPOINT_NAME + '/' + name, date_time=(2026, 10, 2, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            with archive.open(info, 'w', force_zip64=True) as member, safe_file(root, name).open('rb') as source:
                shutil.copyfileobj(source, member, 1024 * 1024)


def verify_zip(path, manifest_sha256):
    with zipfile.ZipFile(path) as archive:
        members = archive.infolist()
        names = [item.filename for item in members]
        require(len(names) == len(set(names)), 'Duplicate ZIP member')
        relative = {}
        for item in members:
            require(item.filename.startswith(CHECKPOINT_NAME + '/'), 'ZIP checkpoint prefix differs')
            name = item.filename[len(CHECKPOINT_NAME) + 1:]
            safe_relative(name)
            require(not item.is_dir() and stat.S_IFMT(item.external_attr >> 16) == stat.S_IFREG,
                    'ZIP member must be a regular file: ' + name)
            relative[name] = item
        require('MANIFEST.json' in relative, 'ZIP package manifest missing')
        manifest_raw = archive.read(relative['MANIFEST.json'])
        require(digest(manifest_raw) == manifest_sha256, 'ZIP manifest SHA differs')
        manifest = json.loads(manifest_raw)
        require(manifest.get('schema_version') == 1 and manifest.get('excluded') == sorted(EXCLUDED),
                'ZIP manifest schema or exclusions differ')
        require('FREEZE_RECEIPT.json' in manifest.get('files', {}), 'ZIP must protect freeze receipt')
        require(set(relative) == set(manifest['files']) | {'MANIFEST.json'}, 'ZIP payload membership differs')
        for name, pin in manifest['files'].items():
            safe_relative(name)
            import hashlib
            hasher = hashlib.sha256()
            with archive.open(relative[name]) as member:
                for block in iter(lambda: member.read(1024 * 1024), b''):
                    hasher.update(block)
            require(hasher.hexdigest() == pin, 'ZIP payload checksum differs: ' + name)
    return {'zip_members': len(members), 'package_manifest_sha256': manifest_sha256}


def extract_verified(path, destination, manifest_sha256):
    """Path-safe fresh extraction for a subsequent external replay."""
    verify_zip(path, manifest_sha256)
    destination = Path(destination).absolute()
    require(not destination.exists(), 'Fresh extraction destination required')
    destination.mkdir(parents=True, exist_ok=False)
    with zipfile.ZipFile(path) as archive:
        for item in archive.infolist():
            target = destination.joinpath(*safe_relative(item.filename))
            require(target.resolve().is_relative_to(destination.resolve()), 'ZIP extraction escapes destination')
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(item) as source, target.open('xb') as output:
                shutil.copyfileobj(source, output, 1024 * 1024)
    return destination / CHECKPOINT_NAME


def build(root, output, registration_sha256, freeze_commit):
    root = Path(root).absolute()
    output = fresh_external(root, output)
    require(root.name == CHECKPOINT_NAME, 'Unexpected checkpoint directory name')
    inputs = verify_frozen(root, registration_sha256, freeze_commit)
    files = payload_files(root)
    require('FREEZE_RECEIPT.json' in files, 'Package must protect freeze receipt')
    manifest = {'schema_version': 1, 'files': files, 'excluded': sorted(EXCLUDED),
                'scope': 'Complete source-free ledger checkpoint, including frozen inputs and recorded outputs; packaging only.'}
    raw = (json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    manifest_path = root / 'MANIFEST.json'
    if manifest_path.exists():
        require(manifest_path.read_bytes() == raw, 'Existing package manifest differs; no overwrite')
    else:
        with manifest_path.open('xb') as destination:
            destination.write(raw)
    output.mkdir(parents=True, exist_ok=False)
    zip_path, repeat = output / ZIP_NAME, output / (ZIP_NAME + '.repeat')
    construct_zip(root, [*files, 'MANIFEST.json'], zip_path)
    construct_zip(root, [*files, 'MANIFEST.json'], repeat)
    require(file_sha(zip_path) == file_sha(repeat) and zip_path.stat().st_size == repeat.stat().st_size,
            'Repeated deterministic ZIP differed')
    repeat.unlink()
    checks = verify_zip(zip_path, digest(raw))
    current = verify_frozen(root, registration_sha256, freeze_commit)
    require({k: v for k, v in current.items() if k != 'package_manifest_sha256'} ==
            {k: v for k, v in inputs.items() if k != 'package_manifest_sha256'},
            'Frozen inputs or receipt changed during packaging')
    require(payload_files(root) == files and manifest_path.read_bytes() == raw,
            'Payload changed during packaging')
    zip_pin = file_sha(zip_path)
    (output / ZIP_SHA_NAME).write_text(zip_pin + '  ' + ZIP_NAME + '\n')
    receipt = {'status': 'PASS_PACKAGING_ONLY', 'science_evaluated': False,
               'repeat_build_identical': True, 'frozen_inputs_unchanged': True,
               'zip_path': str(zip_path), 'zip_sha256': zip_pin, 'zip_bytes': zip_path.stat().st_size,
               'payload_files': len(files), **current, **checks}
    write_json(output / 'PACKAGE_RECEIPT.json', receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--freeze-commit', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--extract-to', type=Path)
    args = parser.parse_args()
    receipt = build(args.checkpoint_root, args.output_dir, args.registration_sha256, args.freeze_commit)
    if args.extract_to is not None:
        extract_verified(receipt['zip_path'], args.extract_to, receipt['package_manifest_sha256'])
    print(json.dumps(receipt, sort_keys=True))


if __name__ == '__main__':
    main()
