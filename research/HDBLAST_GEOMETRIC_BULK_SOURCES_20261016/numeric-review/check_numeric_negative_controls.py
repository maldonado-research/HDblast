#!/usr/bin/env python3
"""Direct normal/optimized semantic rejection processes for numeric review."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
NUMERIC = ROOT.parent/'numeric'


def need(ok, text):
    if not ok:
        raise ValueError(text)


def process(script, arguments, expected, output, optimized):
    command = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])+[str(script), *arguments, '--output', str(output)]
    child = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    need(child.returncode != 0, 'negative control passed '+str(arguments))
    need(expected in child.stderr, 'wrong negative control reason '+str(arguments)+' '+child.stderr[-800:])
    need(not output.exists(), 'failed control wrote a certificate')
    return {'optimized': optimized, 'script': script.name, 'arguments': arguments,
            'exit_nonzero': True, 'expected_stderr': expected, 'output_absent': True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--include-source', action='store_true')
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    records = []
    with tempfile.TemporaryDirectory(prefix='geometric-numeric-controls-', dir='/tmp') as temp:
        temporary = Path(temp)
        for optimized in (False, True):
            for mutation, reason in (
                ('wrong_factorial', 'independent ordinary quotient derivative'),
                ('wrong_scalar_shift', 'independent coordinate-to-frame ordinary derivative s'),
                ('drop_radial_rho', 'independent coordinate-to-frame ordinary derivative A'),
                ('drop_radial_conversion', 'independent coordinate-to-frame ordinary derivative Ayy')):
                output = temporary/f'own-{mutation}-{optimized}.json'
                records.append(process(ROOT/'check_independent_actual_jets.py', ['--inject-error', mutation], reason, output, optimized))
            # Isolated corrupted inputs cannot affect producer or parents.
            for mutation in ('pinned_input_changed', 'certificate_review_mismatch'):
                target = temporary/f'{mutation}-{optimized}'
                target.mkdir()
                pins = json.loads((NUMERIC/'PREDECESSOR_PINS.json').read_text())
                paths = set(pins)|{'actual_first_jets.py', 'PREDECESSOR_PINS.json', 'exterior_error_caps.py'}
                for path in paths:
                    dest = target/path
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(NUMERIC/path, dest)
                if mutation == 'pinned_input_changed':
                    path = target/'inputs/BUMP_DERIVATIVE_BOUNDS.json'
                    path.write_bytes(path.read_bytes()+b'\n')
                    expected = 'predecessor authentication: inputs/BUMP_DERIVATIVE_BOUNDS.json'
                else:
                    path = target/'inputs/APPROXIMATION_ERRORS_NORMAL.json'
                    report = json.loads(path.read_text())
                    report['master_error_budget_by_scaled_total_order']['actual_minus_comparison'][3]['exact_upper'] = '0'
                    path.write_text(json.dumps(report, indent=2)+'\n')
                    raw = path.read_bytes()
                    pins['inputs/APPROXIMATION_ERRORS_NORMAL.json']['sha256'] = hashlib.sha256(raw).hexdigest()
                    pins['inputs/APPROXIMATION_ERRORS_NORMAL.json']['bytes'] = len(raw)
                    (target/'PREDECESSOR_PINS.json').write_text(json.dumps(pins, indent=2)+'\n')
                    expected = 'review binds inherited error certificate'
                output = temporary/f'provenance-{mutation}-{optimized}.json'
                records.append(process(target/'actual_first_jets.py', [], expected, output, optimized))
            if args.include_source:
                for mutation, reason in (
                    ('comparator_only', 'actual solution derivative errors retained rather than comparator substitution'),
                    ('halve_hessian', 'actual full-cell Einstein Hessian negative witness and second-derivative normalization'),
                    ('reverse_metric_forcing', 'registered reduced metric forcing sign and scale'),
                    ('drop_corner_cell', 'complete clipped causal triangle box cover including diagonal cells'),
                    ('omit_shell_normalization', 'physical shell unit-normal affine Hessian and first-row elimination')):
                    output = temporary/f'source-{mutation}-{optimized}.json'
                    records.append(process(NUMERIC/'evaluate_actual_sources.py', ['--inject-error', mutation], reason, output, optimized))
    files = [ROOT/'check_independent_actual_jets.py', NUMERIC/'actual_first_jets.py', NUMERIC/'evaluate_actual_sources.py', ROOT.parent/'projection/source_coefficients.py', ROOT.parent/'projection/GEOMETRIC_COEFFICIENT_MAPS.json']
    out = {'status': 'PASS_INDEPENDENT_NUMERIC_DIRECT_REJECTIONS', 'direct_child_processes': len(records),
           'all_rejections_require': 'nonzero process, mutation-specific stderr, absent output file',
           'records': records, 'input_sha256': {str(path.relative_to(ROOT.parent)): hashlib.sha256(path.read_bytes()).hexdigest() for path in files}}
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(out['status'])


if __name__ == '__main__':
    main()
