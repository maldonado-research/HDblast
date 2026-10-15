"""Copy an assembled source tree into a new authenticated capsule and fixed-layout ZIP.

The source tree is read only. The output directory must not already exist.
No publication, scientific validation, dependency execution or Git action occurs.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import stat
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(name, raw):
    return {'path': name, 'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def scan(source):
    files = {}
    for path in sorted(source.rglob('*')):
        name = path.relative_to(source).as_posix()
        require(not path.is_symlink(), 'Symlink forbidden: ' + name)
        require(path.is_dir() or stat.S_ISREG(path.stat().st_mode), 'Special file forbidden: ' + name)
        require(not set(path.relative_to(source).parts) & {'.git', '__pycache__'},
                'Repository metadata or bytecode cache forbidden: ' + name)
        if path.is_file():
            require('\\' not in name and path.suffix not in {'.pyc', '.pyo'}, 'Noncanonical source: ' + name)
            require(name != 'MANIFEST.json', 'Source is already frozen; assemble a fresh source tree')
            files[name] = path.read_bytes()
    require(bool(files), 'Empty source tree')
    return files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--capsule', required=True)
    parser.add_argument('--created-utc', help='Fixed ISO8601 UTC timestamp for byte-repeatable freezes')
    args = parser.parse_args()
    source, output = args.source.absolute(), args.output_dir.absolute()
    for path in (source, output):
        require(not any(p.is_symlink() for p in [path, *path.parents]), 'Real input/output paths required')
    require(source.is_dir(), 'Source directory missing')
    require(not output.exists(), 'Output must be a new directory')
    require(not output.resolve().is_relative_to(source.resolve()), 'Output must be outside source')
    require(re.fullmatch(r'HDBLAST_[A-Z0-9_]+', args.capsule) is not None, 'Canonical capsule name required')
    created = args.created_utc or datetime.now(timezone.utc).isoformat()
    parsed = datetime.fromisoformat(created.replace('Z', '+00:00'))
    require(parsed.utcoffset() is not None and parsed.utcoffset().total_seconds() == 0, 'UTC timestamp required')
    files = scan(source)
    manifest = {'schema': 'hdblast-frozen-source-v1', 'capsule': args.capsule,
                'created_utc': created,
                'scope': 'Complete source-byte inventory; no mathematical claim is certified by this manifest.',
                'files': [record(n, b) for n, b in files.items()]}
    manifest_raw = (json.dumps(manifest, indent=2) + '\n').encode()
    frozen = dict(files, **{'MANIFEST.json': manifest_raw})
    output.mkdir(parents=True)
    capsule = output / args.capsule
    for name, raw in frozen.items():
        target = capsule / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
    archive = output / (args.capsule + '.zip')
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for name, raw in sorted(frozen.items()):
            info = zipfile.ZipInfo(args.capsule + '/' + name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, raw, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(archive) as zf:
        require(len(zf.namelist()) == len(frozen), 'ZIP member count differs')
        require(set(zf.namelist()) == {args.capsule + '/' + n for n in frozen}, 'ZIP inventory differs')
        for name, raw in frozen.items():
            require(zf.read(args.capsule + '/' + name) == raw == (capsule / name).read_bytes(),
                    'Source/ZIP byte mismatch: ' + name)
    require(scan(source) == files, 'Source changed during freeze')
    receipt = {'schema': 'hdblast-source-freeze-v1', 'status': 'PASS',
               'created_utc': created, 'capsule': args.capsule,
               'source_files': len(files), 'archive_members': len(frozen),
               'manifest': record(args.capsule + '/MANIFEST.json', manifest_raw),
               'archive': record(archive.name, archive.read_bytes()),
               'scope': 'Byte-preserving source freeze and ZIP readback only; not scientific validation or publication.'}
    (output / 'SOURCE_FREEZE_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
