#!/usr/bin/env python3
"""Reassemble and verify a metric memory-followup ZIP outside its parts directory.

This standalone standard-library helper never invokes science or repository code.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

CHECKPOINT_NAME = 'HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP'
ZIP_NAME = CHECKPOINT_NAME + '.zip'
ZIP_SHA_NAME = ZIP_NAME + '.sha256'
PARTS_MANIFEST_NAME = 'PACKAGE_PARTS.json'
PART_BYTES = 64 * 1024 * 1024
MAX_PARTS = 64
PART_NAMES = tuple(ZIP_NAME + '.part' + str(n).zfill(4) for n in range(1, MAX_PARTS+1))
EXCLUDED_FILES = frozenset(('MANIFEST.json', ZIP_NAME, ZIP_SHA_NAME,
                            PARTS_MANIFEST_NAME, *PART_NAMES))
CHUNK_BYTES = 1024 * 1024


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def hash_stream(stream):
    value = hashlib.sha256()
    for block in iter(lambda: stream.read(CHUNK_BYTES), b''):
        value.update(block)
    return value.hexdigest()


def file_sha(path):
    with Path(path).open('rb') as stream:
        return hash_stream(stream)


def valid_sha(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def safe_file(root, name):
    root = root.resolve()
    relative = Path(name)
    require(isinstance(name, str) and name and not relative.is_absolute()
            and '..' not in relative.parts, 'Unsafe payload path: ' + str(name))
    target = root / relative
    require(target.resolve().is_relative_to(root), 'Payload path escapes checkpoint: ' + name)
    require(target.is_file() and not target.is_symlink(), 'Missing or symlink file: ' + name)
    for parent in target.parents:
        if parent == root:
            break
        require(not parent.is_symlink(), 'Symlink directory: ' + name)
    return target


def verify_zip(path, expected_manifest_sha=None, expected_freeze=None):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names) == len(set(names)), 'Duplicate ZIP member')
        manifest_name = CHECKPOINT_NAME + '/MANIFEST.json'
        require(manifest_name in names, 'ZIP package manifest missing')
        raw_manifest = archive.read(manifest_name)
        if expected_manifest_sha is not None:
            require(digest(raw_manifest) == expected_manifest_sha, 'ZIP manifest hash differs from part ledger')
        manifest = json.loads(raw_manifest)
        require(isinstance(manifest.get('files'), dict) and manifest['files'], 'ZIP payload file table missing')
        require(set(manifest.get('excluded', [])) == EXCLUDED_FILES, 'ZIP exclusion constants differ')
        require('FREEZE_RECEIPT.json' in manifest['files'], 'ZIP manifest must hash freeze receipt')
        expected = [CHECKPOINT_NAME + '/' + name for name in sorted([*manifest['files'], 'MANIFEST.json'])]
        require(names == expected, 'ZIP membership/order differs from manifest')
        for name, pin in manifest['files'].items():
            relative = Path(name)
            require(not relative.is_absolute() and '..' not in relative.parts and
                    name not in EXCLUDED_FILES and valid_sha(pin), 'Unsafe or invalid ZIP payload entry')
            with archive.open(CHECKPOINT_NAME + '/' + name) as member:
                require(hash_stream(member) == pin, 'ZIP member bytes differ: ' + name)
        require(archive.testzip() is None, 'ZIP CRC verification failed')
        receipt_raw = archive.read(CHECKPOINT_NAME + '/FREEZE_RECEIPT.json')
        receipt = json.loads(receipt_raw)
        registration_raw = archive.read(CHECKPOINT_NAME + '/FULL_REGISTRATION.json')
        registration = json.loads(registration_raw)
        require(receipt['registration_sha256'] == digest(registration_raw), 'ZIP freeze registration binding differs')
        require(manifest['freeze_receipt_sha256'] == digest(receipt_raw), 'ZIP receipt digest differs')
        if expected_freeze is not None:
            require(receipt['public_freeze_commit'] == expected_freeze, 'ZIP public freeze differs from part ledger')
        require(isinstance(registration.get('files'), dict) and registration['files'], 'ZIP registration missing')
        for name, pin in registration['files'].items():
            require(manifest['files'].get(name) == pin, 'ZIP registration/payload mismatch: ' + name)
        independent_name = 'independent/MANIFEST.json'
        independent_raw = archive.read(CHECKPOINT_NAME + '/' + independent_name)
        require(digest(independent_raw) == receipt['independent_manifest_sha256'] ==
                registration['independent_manifest_sha256'] == registration['files'].get(independent_name),
                'ZIP independent manifest binding differs')
    return {'zip_members': len(names), 'payload_files': len(manifest['files']),
            'package_manifest_sha256': digest(raw_manifest), 'all_members_equal_manifest_hashes': True,
            'archive_crc_verified': True, 'freeze_receipt_included_and_hashed': True}


def reassemble(ledger_path, output, expected_sha=None):
    ledger_path = Path(ledger_path)
    require(ledger_path.is_file() and not ledger_path.is_symlink(), 'Regular part ledger required')
    ledger_path, output = ledger_path.resolve(), Path(output).resolve()
    root = ledger_path.parent
    require(ledger_path.name == PARTS_MANIFEST_NAME and not ledger_path.is_symlink(), 'Canonical part ledger required')
    ledger = json.loads(ledger_path.read_text())
    require(ledger.get('checkpoint') == CHECKPOINT_NAME and ledger.get('zip_name') == ZIP_NAME,
            'Unexpected checkpoint or logical ZIP')
    require(ledger.get('part_bytes') == PART_BYTES and valid_sha(ledger.get('zip_sha256')),
            'Part size or ZIP hash contract differs')
    require(expected_sha is None or expected_sha == ledger['zip_sha256'], 'Requested trusted ZIP SHA differs')
    require(not output.exists() and not output.is_relative_to(root), 'Fresh reconstructed ZIP must be external')
    parts = ledger.get('parts')
    require(isinstance(parts, list) and 0 < len(parts) <= MAX_PARTS, 'Missing or excessive ZIP parts')
    if ledger.get('storage') == 'single':
        require(len(parts) == 1 and parts[0]['path'] == ZIP_NAME and ledger['zip_bytes'] <= PART_BYTES,
                'Invalid single-file package')
    else:
        require(ledger.get('storage') == 'split' and ledger['zip_bytes'] > PART_BYTES,
                'Invalid split package')
        require([p['path'] for p in parts] == list(PART_NAMES[:len(parts)]), 'Part order/names differ')
    require(sum(p['bytes'] for p in parts) == ledger['zip_bytes'], 'Part lengths differ from logical ZIP')
    output.parent.mkdir(parents=True, exist_ok=True)
    total = hashlib.sha256()
    with output.open('xb') as destination:
        for index, part in enumerate(parts):
            require(valid_sha(part.get('sha256')) and isinstance(part.get('bytes'), int) and
                    0 < part['bytes'] <= PART_BYTES, 'Invalid part digest/length')
            if ledger['storage'] == 'split' and index < len(parts)-1:
                require(part['bytes'] == PART_BYTES, 'Non-final split part must be exactly 64MiB')
            source = safe_file(root, part['path'])
            require(source.stat().st_size == part['bytes'], 'Part file length differs')
            item_hash = hashlib.sha256()
            with source.open('rb') as stream:
                for block in iter(lambda: stream.read(CHUNK_BYTES), b''):
                    item_hash.update(block); total.update(block); destination.write(block)
            require(item_hash.hexdigest() == part['sha256'], 'Part SHA mismatch: ' + part['path'])
    require(output.stat().st_size == ledger['zip_bytes'] and total.hexdigest() == ledger['zip_sha256'],
            'Reconstructed ZIP length/hash differs')
    checks = verify_zip(output, ledger['package_manifest_sha256'], ledger['public_freeze_commit'])
    return {'status': 'PASS_REASSEMBLY_ONLY', 'science_evaluated': False, 'zip_path': str(output),
            'zip_sha256': total.hexdigest(), 'zip_bytes': output.stat().st_size,
            'parts_verified': len(parts), **checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parts-manifest', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expected-sha256')
    args = parser.parse_args()
    print(json.dumps(reassemble(args.parts_manifest, args.output, args.expected_sha256), sort_keys=True))


if __name__ == '__main__':
    main()
