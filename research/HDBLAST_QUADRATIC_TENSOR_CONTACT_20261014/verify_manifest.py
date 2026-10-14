"""Verify the frozen source capsule using only the Python standard library."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re


def verify(root, expected):
    manifest = root / 'MANIFEST.json'
    if not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('A complete external manifest SHA256 is required')
    if manifest.is_symlink() or hashlib.sha256(manifest.read_bytes()).hexdigest() != expected:
        raise ValueError('Manifest SHA256 does not match the external pin')
    data = json.loads(manifest.read_text())
    if data['capsule'] != 'HDBLAST_QUADRATIC_TENSOR_CONTACT_20261014':
        raise ValueError('Unexpected capsule identity')
    seen = set()
    for entry in data['files']:
        name = entry['path']
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or str(relative) != name or name in seen:
            raise ValueError('Noncanonical or duplicate manifest path: ' + name)
        seen.add(name)
        path = root / name
        if any(part.is_symlink() for part in [path, *path.parents]) or not path.is_file():
            raise ValueError('Missing or symlinked capsule file: ' + name)
        if not path.resolve().is_relative_to(root.resolve()):
            raise ValueError('File outside capsule: ' + name)
        content = path.read_bytes()
        if len(content) != entry['bytes'] or hashlib.sha256(content).hexdigest() != entry['sha256']:
            raise ValueError('Capsule bytes differ: ' + name)
    if not seen:
        raise ValueError('Empty manifest')
    return {'status': 'PASS', 'files_verified': len(seen), 'manifest_sha256': expected}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected', required=True)
    args = parser.parse_args()
    print(json.dumps(verify(Path(__file__).resolve().parent, args.expected), indent=2))


if __name__ == '__main__':
    main()
