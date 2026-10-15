"""Assemble a byte-pinned, approved source roster into a fresh directory."""
import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import re

CAPSULE = 'HDBLAST_COORDINATE_SHELL_SOURCES_20261015'


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def canonical(name):
    require(isinstance(name, str) and name and '\\' not in name, 'Relative string path required')
    path = PurePosixPath(name)
    require(path.parts and not path.is_absolute() and '..' not in path.parts and path.as_posix() == name,
            'Canonical relative path required')
    require(not ({part.lower() for part in path.parts} &
                 {'internal', 'private', '_internal', '_private', '.git', '__pycache__'}),
            'Private, repository or cache path forbidden')
    require(path.suffix not in {'.pyc', '.pyo'}, 'Bytecode forbidden')
    return name


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--source-root', type=Path, default=task)
    parser.add_argument('--output-dir', type=Path, default=task/'assembled-source')
    parser.add_argument('--receipt', type=Path, default=task/'publication/SOURCE_ASSEMBLY_RECEIPT.json')
    args = parser.parse_args()
    config_raw = args.config.read_bytes()
    config = json.loads(config_raw)
    require(config['schema'] == 'hdblast-approved-source-assembly-v1'
            and config['status'] == 'APPROVED_FOR_FROZEN_SOURCE'
            and config['capsule'] == CAPSULE, 'Approved current source assembly required')
    root, destination = args.source_root.absolute(), args.output_dir.absolute()
    for path in (root, destination):
        require(not any(p.is_symlink() for p in [path, *path.parents]), 'Real source/destination required')
    require(root.is_dir() and not destination.exists(), 'Existing source and new assembly directory required')
    require(not args.receipt.resolve().is_relative_to(destination.resolve()), 'Assembly receipt must remain outside source')
    require(not args.receipt.exists() and not args.receipt.is_symlink(), 'New assembly receipt required')
    files = {}
    for row in config['files']:
        name, origin = canonical(row['path']), canonical(row.get('source', row['path']))
        require(name != 'MANIFEST.json' and name not in files, 'Unique nonreserved output path required')
        path = root/origin
        require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents])
                and path.resolve().is_relative_to(root.resolve()), 'Regular reviewed input required: '+origin)
        raw = path.read_bytes()
        require(type(row['bytes']) is int and len(raw) == row['bytes']
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None
                and sha256(raw).hexdigest() == row['sha256'], 'Reviewed source bytes differ: '+origin)
        files[name] = (origin, raw)
    require(files and {'reproduce.py', 'REPLAY_ROSTER.json', 'verify_manifest.py', 'requirements.txt'} <= set(files),
            'Nonempty complete runnable source roster required')
    require(not any(a != b and b.startswith(a+'/') for a in files for b in files), 'Overlapping file/directory paths')
    destination.mkdir(parents=True)
    for name, (origin, raw) in files.items():
        require((root/origin).read_bytes() == raw, 'Source changed during assembly')
        target = destination/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        require(target.read_bytes() == raw, 'Assembly write/readback differs')
    for origin, raw in files.values():
        require((root/origin).read_bytes() == raw, 'Source changed before assembly completed')
    records = [{'path': name, 'source': origin, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()}
               for name, (origin, raw) in sorted(files.items())]
    receipt = {'status': 'PASS', 'capsule': CAPSULE, 'source_files': len(records),
               'approved_config': {'bytes': len(config_raw), 'sha256': sha256(config_raw).hexdigest()},
               'files': records,
               'scope': 'Explicit approved byte-pinned source assembly only; no scientific result, replay or publication is certified.'}
    args.receipt.write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'source_files': len(records)}))


if __name__ == '__main__':
    main()
