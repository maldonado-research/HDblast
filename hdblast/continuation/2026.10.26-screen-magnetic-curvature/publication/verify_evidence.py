#!/usr/bin/env python3
"""Check the complete local evidence inventory and retained byte hashes, read-only."""
import hashlib
import json
from pathlib import Path


def need(value, label):
    if not value:
        raise ValueError(label)


def main():
    root = Path(__file__).resolve().parent
    manifest_path = root / 'EVIDENCE_MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    need(manifest['schema'] == 'HDBLAST_SCREEN_MAGNETIC_PUBLICATION_EVIDENCE_MANIFEST_V1',
         'EXPECTED_MANIFEST_SCHEMA')
    rows = manifest['files']
    paths = [row['path'] for row in rows]
    need(len(paths) == len(set(paths)), 'UNIQUE_PATHS')
    need(manifest['files_including_manifest'] == len(paths) + 1, 'DECLARED_FILE_COUNT')
    for row in rows:
        path = Path(row['path'])
        need(not path.is_absolute() and '..' not in path.parts
             and path.as_posix() == row['path'] and row['path'] != manifest_path.name,
             'SAFE_RELATIVE_PATH')
        target = root / path
        need(target.is_file() and not target.is_symlink(), 'REGULAR_FILE:' + str(path))
        data = target.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'EXACT_BYTES:' + str(path))
    actual = sorted(p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file())
    need(actual == sorted(paths + [manifest_path.name]), 'CLOSED_FILE_INVENTORY')
    print(json.dumps({'status':'PASS','files_including_manifest':len(actual),
                      'manifest_sha256':hashlib.sha256(manifest_path.read_bytes()).hexdigest()}))


if __name__ == '__main__':
    main()
