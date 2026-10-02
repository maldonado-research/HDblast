#!/usr/bin/env python3
"""Verify faithful public input copies without importing or evaluating them."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    require(not args.output.exists(), 'Refusing to overwrite previous evidence')
    require(not args.output.resolve().is_relative_to(root), 'Use an external output path')
    manifest_path = root / 'INHERITED_REFERENCE_MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    records = []
    require(manifest['schema_version'] == 1 and manifest['files'], 'Malformed source manifest')
    for name, entry in manifest['files'].items():
        relative = Path(name)
        path = (root / relative).resolve()
        require(not relative.is_absolute() and '..' not in relative.parts,
                'Inherited input path escape')
        require(path.is_relative_to(root) and not (root / relative).is_symlink(),
                'Inherited input must be a regular in-package file')
        require(entry['original_public_commit'] == 'b807a4d549a40bc78b66122cb5719987091d2ce7',
                'Inherited public source commit differs')
        require(name == 'reference_inputs/' + entry['original_repository_path'],
                'Relocated public source path differs')
        raw = path.read_bytes()
        require(len(raw) == entry['bytes'], 'Inherited byte length differs: ' + name)
        require(hashlib.sha256(raw).hexdigest() == entry['sha256'],
                'Inherited bytes differ: ' + name)
        records.append({'path': name, 'sha256': entry['sha256'], 'status': 'PASS'})
    result = {'status': 'PASS', 'faithful_public_inputs': len(records),
              'checks': records, 'python_optimization': sys.flags.optimize,
              'manifest_sha256': hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'physical_evaluations': 0,
              'scope': 'Byte identities only; no source imports, roots, modes or response evaluation.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': 'PASS', 'faithful_public_inputs': len(records)}))


if __name__ == '__main__':
    main()
