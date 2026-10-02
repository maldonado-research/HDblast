#!/usr/bin/env python3
"""Build deterministic ZIP and <=64MiB GitHub parts; never invoke science.

FREEZE_RECEIPT is part of the payload manifest. Logical ZIP, its SHA, the part
ledger and the fixed part namespace are excluded to avoid receipt/hash cycles.
The unchanged replay driver understands the literal manifest exclusion list.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

from reassemble_package import (CHECKPOINT_NAME, ZIP_NAME, ZIP_SHA_NAME,
    PARTS_MANIFEST_NAME, PART_BYTES, MAX_PARTS, PART_NAMES, EXCLUDED_FILES,
    CHUNK_BYTES, digest, file_sha, require, safe_file, valid_sha, verify_zip)

REQUIRED_RECEIPT_FIELDS = ('public_freeze_commit', 'registration_sha256', 'independent_manifest_sha256')


def verify_inputs(root):
    registration_raw = safe_file(root, 'FULL_REGISTRATION.json').read_bytes()
    registration = json.loads(registration_raw)
    receipt_raw = safe_file(root, 'FREEZE_RECEIPT.json').read_bytes()
    receipt = json.loads(receipt_raw)
    require(all(field in receipt for field in REQUIRED_RECEIPT_FIELDS), 'Incomplete post-freeze receipt')
    commit = receipt['public_freeze_commit']
    require(isinstance(commit, str) and len(commit) == 40 and all(c in '0123456789abcdef' for c in commit),
            'Invalid public freeze commit')
    require(receipt['registration_sha256'] == digest(registration_raw), 'Receipt registration hash mismatch')
    require(isinstance(registration.get('files'), dict) and registration['files'], 'Missing registered file table')
    for name, pin in registration['files'].items():
        require(name not in EXCLUDED_FILES and valid_sha(pin), 'Invalid or bookkeeping registered input: ' + name)
        require(file_sha(safe_file(root, name)) == pin, 'Frozen input changed: ' + name)
    raw = safe_file(root, 'independent/MANIFEST.json').read_bytes()
    require(digest(raw) == receipt['independent_manifest_sha256'] == registration['independent_manifest_sha256'],
            'Independent manifest pin mismatch')
    require(registration['files'].get('independent/MANIFEST.json') == digest(raw),
            'Independent manifest omitted from registration')
    entries = json.loads(raw)['files']
    if isinstance(entries, dict):
        entries = [{'path': name, 'sha256': pin} for name, pin in entries.items()]
    require(isinstance(entries, list) and entries, 'Missing independent source manifest')
    seen = set()
    for item in entries:
        relative = item['path']; name = 'independent/' + relative
        require(relative not in seen, 'Duplicate independent source path'); seen.add(relative)
        require(file_sha(safe_file(root, name)) == item['sha256'], 'Independent frozen input changed: ' + relative)
        require(registration['files'].get(name) == item['sha256'], 'Independent input differs from registration: ' + relative)
    return {'public_freeze_commit': commit, 'registration_sha256': digest(registration_raw),
            'independent_manifest_sha256': digest(raw), 'freeze_receipt_sha256': digest(receipt_raw)}


def payload_files(root):
    files = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'Symlink forbidden: ' + str(path))
        if path.is_file():
            name = path.relative_to(root).as_posix()
            require('__pycache__' not in path.parts and path.suffix != '.pyc', 'Cache forbidden: ' + name)
            if name not in EXCLUDED_FILES:
                files[name] = file_sha(safe_file(root, name))
    return files


def construct_zip(root, names, destination):
    with zipfile.ZipFile(destination, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9,
                         allowZip64=True) as archive:
        for name in names:
            info = zipfile.ZipInfo(CHECKPOINT_NAME + '/' + name, date_time=(2026, 10, 2, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED; info.create_system = 3
            info.external_attr = 0o100644 << 16
            with archive.open(info, 'w', force_zip64=True) as member, safe_file(root, name).open('rb') as source:
                shutil.copyfileobj(source, member, CHUNK_BYTES)


def equal_files(first, second):
    with first.open('rb') as a, second.open('rb') as b:
        while True:
            left, right = a.read(CHUNK_BYTES), b.read(CHUNK_BYTES)
            if left != right:
                return False
            if not left:
                return True


def build(root, external_dir=None):
    root = Path(root).resolve()
    require(root.name == CHECKPOINT_NAME, 'Unexpected checkpoint name')
    inputs = verify_inputs(root)
    files = payload_files(root)
    require('FREEZE_RECEIPT.json' in files, 'Package must protect freeze receipt')
    manifest = {'schema_version': 2, 'files': files, 'excluded': sorted(EXCLUDED_FILES),
                'freeze_receipt_sha256': inputs['freeze_receipt_sha256'],
                'part_bytes': PART_BYTES, 'maximum_parts': MAX_PARTS,
                'scope': 'Exact curated metric payload including freeze receipt; only explicit package bookkeeping excluded.'}
    manifest_raw = (json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
    (root/'MANIFEST.json').write_bytes(manifest_raw)
    external = (Path(external_dir).resolve() if external_dir is not None else
                Path(tempfile.mkdtemp(prefix='hdblast-metric-memory-package-')).resolve())
    require(not external.is_relative_to(root), 'Full ZIP staging must be external to checkpoint')
    external.mkdir(parents=True, exist_ok=True)
    zip_path, second = external/ZIP_NAME, external/(ZIP_NAME+'.repeat')
    require(not zip_path.exists() and not second.exists(), 'Fresh external ZIP destination required')
    names = sorted([*files, 'MANIFEST.json'])
    construct_zip(root, names, zip_path); construct_zip(root, names, second)
    require(equal_files(zip_path, second), 'Repeated deterministic archive construction differed')
    second.unlink()
    checks = verify_zip(zip_path, digest(manifest_raw), inputs['public_freeze_commit'])
    require(verify_inputs(root) == inputs and payload_files(root) == files,
            'Frozen inputs, receipt or payload changed during packaging')
    require((root/'MANIFEST.json').read_bytes() == manifest_raw, 'Package manifest changed during packaging')
    size, pin = zip_path.stat().st_size, file_sha(zip_path)
    require(size <= PART_BYTES*MAX_PARTS, 'Logical ZIP exceeds fixed 64-part namespace')
    for name in PART_NAMES:
        path = root/name
        if path.exists():
            require(path.is_file() and not path.is_symlink(), 'Invalid previous bookkeeping part')
            path.unlink()
    tracked_zip = root/ZIP_NAME
    if tracked_zip.exists():
        require(tracked_zip.is_file() and not tracked_zip.is_symlink(), 'Invalid previous bookkeeping ZIP')
        tracked_zip.unlink()
    parts = []
    if size <= PART_BYTES:
        shutil.copyfile(zip_path, tracked_zip)
        parts = [{'path': ZIP_NAME, 'bytes': size, 'sha256': pin}]
        storage = 'single'
    else:
        storage = 'split'
        with zip_path.open('rb') as source:
            remaining = size
            for name in PART_NAMES:
                if not remaining:
                    break
                length = min(PART_BYTES, remaining)
                with (root/name).open('xb') as destination:
                    left = length
                    while left:
                        block = source.read(min(CHUNK_BYTES, left))
                        require(block, 'Unexpected EOF while splitting ZIP')
                        destination.write(block); left -= len(block)
                parts.append({'path': name, 'bytes': length, 'sha256': file_sha(root/name)})
                remaining -= length
            require(remaining == 0 and not source.read(1), 'Split ZIP length mismatch')
    (root/ZIP_SHA_NAME).write_text(pin+'  '+ZIP_NAME+'\n')
    ledger = {'schema_version': 1, 'checkpoint': CHECKPOINT_NAME, 'zip_name': ZIP_NAME,
              'zip_sha256': pin, 'zip_bytes': size, 'storage': storage,
              'part_bytes': PART_BYTES, 'parts': parts, 'package_manifest_sha256': digest(manifest_raw),
              **inputs, 'scope': 'Packaging integrity only; scientific experiment status is unchanged.'}
    (root/PARTS_MANIFEST_NAME).write_text(json.dumps(ledger, indent=2, sort_keys=True)+'\n')
    require(payload_files(root) == files and verify_inputs(root) == inputs, 'Final packaging changed protected payload')
    return {'status': 'PASS_PACKAGING_ONLY', 'science_evaluated': False,
            'payload_files': len(files), 'zip_members': len(names), 'zip_bytes': size, 'zip_sha256': pin,
            'zip_path': str(zip_path), 'tracked_storage': storage, 'tracked_zip_present': tracked_zip.exists(),
            'tracked_parts': parts, 'part_bytes': PART_BYTES, 'repeat_build_identical': True,
            'frozen_inputs_unchanged': True, **checks, **inputs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--external-dir', type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.checkpoint, args.external_dir), sort_keys=True))


if __name__ == '__main__':
    main()
