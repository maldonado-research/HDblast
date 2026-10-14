"""Fresh replay of six independent full-contact coefficient verifiers."""
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
    ('astra-bulk/verify_scalar_bulk.py', 'astra-bulk/normal.json', ['--model-root', 'INPUTS']),
    ('sol-bulk/checker/direct_metric_checker.py', 'sol-bulk/checker/normal.json', ['--include-controls']),
    ('sol-bulk/verify_angular_selection.py', 'sol-bulk/ANGULAR_NORMAL.json', []),
    ('boundary/verify_scalar_boundary.py', 'boundary/verification-normal.json', []),
    ('boundary/principal/verify_principal_symbol.py', 'boundary/principal/NORMAL.json', []),
    ('harmonic-referee/verify_harmonic_contact.py', 'harmonic-referee/normal.json', []),
]
CONTROLS = [
    ('astra-bulk/verify_scalar_bulk.py', ['--inject-corruption'], 'direct quadratic Ricci 01'),
    ('sol-bulk/checker/direct_metric_checker.py', ['--control', 'omit-d2-factor'], 'differences'),
    ('sol-bulk/verify_angular_selection.py', ['--wrong-formula', 'wrong-l4-coefficient'], 'stress-SVT-decomposition'),
    ('boundary/verify_scalar_boundary.py', ['--inject-error'], 'quadratic_extrinsic_curvature_factorial'),
    ('boundary/principal/verify_principal_symbol.py', ['--inject-error'], 'normal harmonic trace-reversal identity after Israel'),
    ('harmonic-referee/verify_harmonic_contact.py', ['--inject-corruption', 'subsidiary'], 'subsidiary_coordinate_Robin_0'),
    ('harmonic-referee/verify_harmonic_contact.py', ['--inject-corruption', 'second-source'], 'second_boundary_relative_lapse_source'),
]


def scientific(data):
    return {k: v for k, v in data.items() if k not in ('optimized', 'python_optimized')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location('contact_integrity', root / 'verify_manifest.py')
    integrity = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(integrity)
    authentication = integrity.verify(root, args.manifest_sha256)
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise RuntimeError('New or empty output directory required')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for optimized in (False, True):
        mode = 'optimized' if optimized else 'normal'
        command = [sys.executable, '-B'] + (['-O'] if optimized else [])
        workdirs = {}
        for name, saved, extras in STEPS:
            work = output / mode / Path(name).parent / Path(name).stem
            work.mkdir(parents=True)
            script = work / Path(name).name
            shutil.copy2(root / name, script)
            workdirs[name] = work
            extras = [str(root / 'inputs') if value == 'INPUTS' else value for value in extras]
            started = time.monotonic()
            with (work / 'execution.log').open('w') as log:
                result = subprocess.run(command + [str(script), '--output', 'fresh.json'] + extras,
                                        cwd=work, stdout=log, stderr=subprocess.STDOUT)
            fresh = work / 'fresh.json'
            if result.returncode != 0 or not fresh.is_file():
                raise RuntimeError('Fresh checker failed: ' + name)
            actual = json.loads(fresh.read_text())
            expected = json.loads((root / saved).read_text())
            if actual.get('status') not in ('PASS', 'passed', 'pass') or scientific(actual) != scientific(expected):
                raise RuntimeError('Scientific output differs: ' + name)
            rows.append({'mode': mode, 'script': name, 'exit_code': 0,
                         'scientific_output_matches': True,
                         'wall_seconds': round(time.monotonic() - started, 3),
                         'output_sha256': hashlib.sha256(fresh.read_bytes()).hexdigest()})
        for index, (name, flags, reason) in enumerate(CONTROLS):
            work = workdirs[name]
            result = subprocess.run(command + [str(work / Path(name).name), '--output',
                                     'corrupt-' + str(index) + '.json'] + flags,
                                    cwd=work, text=True, capture_output=True)
            diagnostic = result.stdout + result.stderr
            (work / ('corrupt-' + str(index) + '.log')).write_text(diagnostic)
            if result.returncode != 1 or reason not in diagnostic:
                raise RuntimeError('Actual corruption did not reject for expected reason: ' + name)
            rows.append({'mode': mode, 'script': name, 'flags': flags,
                         'expected_failure_exit_code': 1, 'expected_reason_detected': True})
    receipt = {'status': 'PASS', 'authentication': authentication, 'steps': rows,
               'python': sys.version,
               'scope': 'Fresh exact identity and formula replay; no numerical spacetime or computer-certified analytic existence proof.'}
    (output / 'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'successful_checker_runs': 12,
                      'actual_corruption_controls': 14, 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
