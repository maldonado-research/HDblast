#!/usr/bin/env python3
"""Build and verify the portable checkpoint without altering frozen inputs."""
import hashlib
import io
import json
from pathlib import Path
import zipfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    root = Path(__file__).resolve().parents[1]
    registration = json.loads((root / 'FULL_REGISTRATION.json').read_text())
    for name, digest in registration['files'].items():
        require(sha((root / name).read_bytes()) == digest, 'Frozen input changed: ' + name)
    archive_name = root.name + '.zip'
    excluded = {'MANIFEST.json', archive_name, archive_name + '.sha256'}
    files = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'Symlink forbidden: ' + str(path))
        if path.is_file():
            name = path.relative_to(root).as_posix()
            require('__pycache__' not in path.parts and path.suffix != '.pyc', 'Cache forbidden')
            if name not in excluded:
                files[name] = sha(path.read_bytes())
    manifest = {'schema_version': 1, 'files': files,
                'excluded': sorted(excluded),
                'scope': 'Exact complete curated payload; package bookkeeping excluded.'}
    (root / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    names = sorted([*files, 'MANIFEST.json'])

    def construct():
        output = io.BytesIO()
        with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name in names:
                item = zipfile.ZipInfo(root.name + '/' + name, date_time=(2026, 10, 2, 0, 0, 0))
                item.compress_type = zipfile.ZIP_DEFLATED
                item.create_system = 3
                item.external_attr = 0o100644 << 16
                archive.writestr(item, (root / name).read_bytes(), compresslevel=9)
        return output.getvalue()

    raw = construct()
    require(raw == construct(), 'Repeated deterministic archive construction differed')
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        require(archive.namelist() == [root.name + '/' + name for name in names], 'Archive membership')
        for name in names:
            require(archive.read(root.name + '/' + name) == (root / name).read_bytes(), 'Archive bytes: ' + name)
    (root / archive_name).write_bytes(raw)
    (root / (archive_name + '.sha256')).write_text(sha(raw) + '  ' + archive_name + '\n')
    print(json.dumps({'status': 'PASS', 'payload_files': len(files), 'package_files': len(files) + 3,
                      'zip_bytes': len(raw), 'zip_sha256': sha(raw),
                      'repeat_build_identical': True, 'all_members_equal_source_bytes': True}))


if __name__ == '__main__':
    main()
