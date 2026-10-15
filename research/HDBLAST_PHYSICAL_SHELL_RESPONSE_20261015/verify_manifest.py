"""Authenticate the full physical-shell-response capsule using the standard library."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re

CAPSULE = 'HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015'


def verify(root, expected):
    root = Path(root).absolute()
    if any(p.is_symlink() for p in [root, *root.parents]):
        raise ValueError('Real capsule path required')
    root = root.resolve()
    if not re.fullmatch(r'[0-9a-f]{64}', expected):
        raise ValueError('Complete external manifest SHA256 required')
    manifest = root / 'MANIFEST.json'
    if manifest.is_symlink() or hashlib.sha256(manifest.read_bytes()).hexdigest() != expected:
        raise ValueError('External manifest pin differs')
    data = json.loads(manifest.read_text())
    if data['capsule'] != CAPSULE:
        raise ValueError('Capsule identity differs')
    seen = set()
    for row in data['files']:
        name = row['path']
        relative = PurePosixPath(name)
        if (not name or '\\' in name or relative.is_absolute() or '..' in relative.parts
                or relative.as_posix() != name or name in seen or name == 'MANIFEST.json'):
            raise ValueError('Noncanonical, reserved or duplicate source path')
        seen.add(name)
        path = root / name
        if any(p.is_symlink() for p in [path, *path.parents]) or not path.is_file():
            raise ValueError('Real regular source required: ' + name)
        if not path.resolve().is_relative_to(root):
            raise ValueError('Source outside capsule')
        raw = path.read_bytes()
        if len(raw) != row['bytes'] or hashlib.sha256(raw).hexdigest() != row['sha256']:
            raise ValueError('Frozen bytes differ: ' + name)
    actual = set()
    for path in root.rglob('*'):
        if path.is_symlink():
            raise ValueError('Symlink in frozen tree')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    if not seen or actual != seen | {'MANIFEST.json'}:
        raise ValueError('Frozen tree contains missing or extra files')
    return {'status': 'PASS', 'files_verified': len(seen), 'manifest_sha256': expected}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected', required=True)
    args = parser.parse_args()
    print(json.dumps(verify(Path(__file__).absolute().parent, args.expected), indent=2))
