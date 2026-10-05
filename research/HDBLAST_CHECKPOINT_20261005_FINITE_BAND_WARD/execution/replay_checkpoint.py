"""Authenticate a frozen mathematical checkpoint and replay its exact proofs."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys


def require(value, message):
    if not value:
        raise RuntimeError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def verify(root, expected):
    manifest = root / 'MANIFEST.json'
    require(not manifest.is_symlink() and manifest.is_file(), 'Manifest is not a regular file')
    raw = manifest.read_bytes()
    require(sha(raw) == expected, 'Manifest hash differs from the external pin')
    ledger = json.loads(raw)
    expected_paths = set()
    for item in ledger['files']:
        name = item['path']
        path = PurePosixPath(name)
        require(not path.is_absolute() and all(x not in ('', '.', '..') for x in name.split('/')),
                'Unsafe manifest path')
        require(name not in expected_paths and name != 'MANIFEST.json', 'Duplicate/self manifest entry')
        expected_paths.add(name)
        file = root / path
        require(not file.is_symlink() and file.is_file(), 'Payload is not a regular file: ' + name)
        body = file.read_bytes()
        require(len(body) == item['bytes'] and sha(body) == item['sha256'], 'Payload bytes differ: ' + name)
    actual = set()
    for file in root.rglob('*'):
        require(not file.is_symlink(), 'Symlink in checkpoint')
        if file.is_file():
            actual.add(file.relative_to(root).as_posix())
        else:
            require(file.is_dir(), 'Special file in checkpoint')
    require(actual == expected_paths | {'MANIFEST.json'}, 'Complete payload membership differs')
    return len(actual)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--checkpoint', type=Path, required=True)
    ap.add_argument('--expected-manifest-sha256', required=True)
    ap.add_argument('--output-directory', type=Path, required=True)
    args = ap.parse_args()
    root, output = args.checkpoint.resolve(), args.output_directory.resolve()
    require(root.is_dir() and not args.checkpoint.is_symlink(), 'Invalid checkpoint directory')
    require(output != root and root not in output.parents and output not in root.parents,
            'Output must be outside and disjoint from checkpoint')
    require(not output.exists(), 'Fresh output directory required')
    payload_files = verify(root, args.expected_manifest_sha256)
    import sympy
    require(sympy.__version__ == '1.14.0', 'Replay requires pinned SymPy 1.14.0')
    output.mkdir(parents=True)
    routes = [
        ('primary', 'verify_finite_band.py', 'EXACT_CHECKS_FINAL_NORMAL.json', 'EXACT_CHECKS_FINAL_OPTIMIZED.json', ()),
        ('theory', 'theory/verify_finite_band_theory.py', 'theory/SYMBOLIC_RELEASE_PHASE_NORMAL.json',
         'theory/SYMBOLIC_RELEASE_PHASE_OPTIMIZED.json', ('python_version', 'python_optimization')),
        ('independent', 'independent/verify_independent_ward.py', 'independent/EXACT_WARD_FINAL_NORMAL.json',
         'independent/EXACT_WARD_FINAL_OPTIMIZED.json', ('python_optimized',)),
    ]
    results = []
    for label, program, normal, optimized, runtime_fields in routes:
        for mode, reference in (('normal', normal), ('optimized', optimized)):
            destination = output / (label + '_' + mode + '.json')
            command = [sys.executable, '-B'] + (['-O'] if mode == 'optimized' else [])
            command += [str(root / program), '--output', str(destination)]
            execution = subprocess.run(command, capture_output=True, timeout=120)
            (output / (label + '_' + mode + '.stdout.log')).write_bytes(execution.stdout)
            (output / (label + '_' + mode + '.stderr.log')).write_bytes(execution.stderr)
            require(execution.returncode == 0, 'Verifier failed: ' + label + '/' + mode)
            produced_raw, reference_raw = destination.read_bytes(), (root / reference).read_bytes()
            produced, expected = json.loads(produced_raw), json.loads(reference_raw)
            for key in runtime_fields:
                produced.pop(key, None)
                expected.pop(key, None)
            require(produced == expected, 'Scientific receipt differs: ' + label + '/' + mode)
            if not runtime_fields:
                require(produced_raw == reference_raw, 'Primary exact receipt byte reproduction differs')
            results.append({'route': label, 'mode': mode, 'scientific_receipt_equal': True,
                            'byte_equal': produced_raw == reference_raw,
                            'ignored_runtime_fields': list(runtime_fields),
                            'produced_sha256': sha(produced_raw), 'reference_sha256': sha(reference_raw)})
    require(verify(root, args.expected_manifest_sha256) == payload_files, 'Payload changed during replay')
    receipt = {'status': 'PASS_FRESH_EXACT_FINITE_BAND_REPLAY', 'manifest_sha256': args.expected_manifest_sha256,
               'payload_files': payload_files, 'routes': results, 'all_payload_bytes_unchanged': True,
               'python_version': sys.version.split()[0], 'sympy_version': sympy.__version__,
               'source_callbacks': 0, 'quantum_arrays_decoded': 0, 'numerical_trajectory_cases': 0,
               'full_twelve_case_pressure_contact_certificate': 'UNRESOLVED'}
    (output / 'REPLAY_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': receipt['status'], 'payload_files': payload_files, 'verifier_executions': len(results)}))


if __name__ == '__main__':
    main()
