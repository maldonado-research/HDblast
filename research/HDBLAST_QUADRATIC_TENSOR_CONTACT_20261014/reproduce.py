"""Replay the six frozen independent tensor checkers in fresh processes."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time


STEPS = [
    ('bulk/verify_quadratic_tt.py', 'bulk/CHECKS_NORMAL.json'),
    ('bulk/verify_tensor_operator.py', 'bulk/TENSOR_OPERATOR_NORMAL.json'),
    ('bulk/verify_sphere_projection.py', 'bulk/SPHERE_NORMAL.json'),
    ('bulk/verify_regular_amplitude.py', 'bulk/REGULAR_AMPLITUDE_NORMAL.json'),
    ('independent/check_full_quadratic_tt.py', 'independent/full_quadratic_tt_results.json'),
    ('boundary/verify_tensor_contact.py', 'boundary/tensor_contact_normal.json'),
]


def scientific(value):
    return {k: v for k, v in value.items()
            if k not in ('python_optimization_level', 'elapsed_seconds')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('tensor_integrity', root / 'verify_manifest.py')
    integrity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(integrity)
    authentication = integrity.verify(root, args.manifest_sha256)
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('Fresh empty output directory required')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for optimized in (False, True):
        mode = 'optimized' if optimized else 'normal'
        for script_name, saved_name in STEPS:
            work = output / mode / Path(script_name).parent / Path(script_name).stem
            work.mkdir(parents=True)
            script = work / Path(script_name).name
            shutil.copy2(root / script_name, script)
            cmd = [sys.executable, '-B'] + (['-O'] if optimized else [])
            cmd += [str(script), '--output', 'fresh.json']
            started = time.monotonic()
            with (work / 'execution.log').open('w') as log:
                result = subprocess.run(cmd, cwd=work, stdout=log, stderr=subprocess.STDOUT)
            fresh = work / 'fresh.json'
            if result.returncode != 0 or not fresh.is_file():
                raise RuntimeError('Checker failed or fresh output missing: ' + script_name)
            actual = json.loads(fresh.read_text())
            expected = json.loads((root / saved_name).read_text())
            if actual.get('status') not in ('passed', 'pass', 'PASS') or scientific(actual) != scientific(expected):
                raise RuntimeError('Scientific output differs: ' + script_name)
            rows.append({'mode': mode, 'script': script_name, 'exit_code': result.returncode,
                         'scientific_output_matches': True,
                         'wall_seconds': round(time.monotonic() - started, 3),
                         'sha256': hashlib.sha256(fresh.read_bytes()).hexdigest()})
        for name, flag, reason in [
            ('bulk/verify_quadratic_tt.py', '--inject-failure', 'injected_matter_only_claim'),
            ('boundary/verify_tensor_contact.py', '--inject-error', 'full radial-time canonical conjugation'),
        ]:
            work = output / mode / Path(name).parent / Path(name).stem
            cmd = [sys.executable, '-B'] + (['-O'] if optimized else [])
            cmd += [str(work / Path(name).name), flag]
            result = subprocess.run(cmd, cwd=work, text=True, capture_output=True)
            diagnostic = result.stdout + result.stderr
            (work / 'corrupted-execution.log').write_text(diagnostic)
            if result.returncode != 1 or reason not in diagnostic:
                raise RuntimeError('Actual corruption control did not fail as expected: ' + name)
            rows.append({'mode': mode, 'script': name, 'control': flag,
                         'expected_failure_exit_code': result.returncode,
                         'expected_reason_detected': True})
    receipt = {'status': 'PASS', 'authentication': authentication,
               'steps': rows, 'python': sys.version,
               'scope': 'Fresh exact algebra and operator replay; no numerical spacetime or full coupled contact proof.'}
    (output / 'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'successful_checker_runs': 12,
                      'actual_corruption_controls': 4, 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
