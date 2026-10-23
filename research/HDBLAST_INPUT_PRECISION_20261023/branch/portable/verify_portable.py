#!/usr/bin/env python3
"""Relocate and check the portable replay, including expected-failure controls."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

def require(test, message):
    if not test:
        raise ValueError(message)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path.cwd()/'PORTABLE_CHECKS.json')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    cases = []
    with tempfile.TemporaryDirectory(prefix='hdblast-portable-replay-') as dirname:
        relocated = Path(dirname)
        for filename in ('replay_branch.py', 'NUMERICAL_REFERENCE.json'):
            shutil.copyfile(here/filename, relocated/filename)
        arbitrary_cwd = relocated/'unrelated-working-directory'
        arbitrary_cwd.mkdir()
        executable = str(relocated/'replay_branch.py')
        reports = []
        for optimized in (False, True):
            command = [sys.executable, '-I', '-S', '-B'] + (['-O'] if optimized else [])
            report = relocated/('optimized.json' if optimized else 'normal.json')
            run = subprocess.run(command+[executable, '--output', str(report)],
                                 cwd=arbitrary_cwd, text=True, capture_output=True)
            require(run.returncode == 0, 'positive replay failed: '+run.stdout+run.stderr)
            data = json.loads(report.read_text())
            require(data['status'] == 'PASS_CONDITIONAL_RETAINED_BRANCH_RATIONAL_REPLAY',
                    'positive replay status')
            reports.append(report.read_bytes())
            cases.append({'mode': 'optimized' if optimized else 'normal', 'exit': run.returncode,
                          'status': data['status']})
        require(reports[0] == reports[1], 'normal/optimized byte identity')
        failures = {'coefficient': 'degree4 residual cancellation',
                    'damping': 'secant damping bound', 'rounding': 'outward beta rounding'}
        for optimized in (False, True):
            command = [sys.executable, '-I', '-S', '-B'] + (['-O'] if optimized else [])
            for control, reason in failures.items():
                run = subprocess.run(command+[executable, '--negative-control', control],
                                     cwd=arbitrary_cwd, text=True, capture_output=True)
                require(run.returncode != 0, 'negative control unexpectedly passed: '+control)
                failure = json.loads(run.stdout)
                require(failure == {'status': 'FAIL', 'reason': reason},
                        'negative control failed at wrong stage: '+control)
                cases.append({'mode': 'optimized' if optimized else 'normal',
                              'negative_control': control, 'exit': run.returncode,
                              'failure_reason': reason})
        reference_path = relocated/'NUMERICAL_REFERENCE.json'
        corrupted = json.loads(reference_path.read_text())
        corrupted['beta_radius'] = '1/100000000000000000000'
        reference_path.write_text(json.dumps(corrupted)+'\n')
        for optimized in (False, True):
            command = [sys.executable, '-I', '-S', '-B'] + (['-O'] if optimized else [])
            run = subprocess.run(command+[executable], cwd=arbitrary_cwd,
                                 text=True, capture_output=True)
            expected = {'status': 'FAIL', 'reason': 'exact portable numerical reference comparison'}
            require(run.returncode != 0 and json.loads(run.stdout) == expected,
                    'corrupted reference detection')
            cases.append({'mode': 'optimized' if optimized else 'normal',
                          'negative_control': 'corrupted_reference', 'exit': run.returncode,
                          'failure_reason': expected['reason']})
        result = {'status': 'PASS_PORTABLE_RELOCATION_AND_NEGATIVE_CONTROLS',
                  'stdlib_only': True, 'runtime_inputs': ['replay_branch.py', 'NUMERICAL_REFERENCE.json'],
                  'relocated_execution': True, 'arbitrary_working_directory': True,
                  'child_python_flags': ['-I', '-S', '-B'],
                  'normal_optimized_byte_identical': True,
                  'positive_report_sha256': hashlib.sha256(reports[0]).hexdigest(), 'cases': cases}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('status', 'relocated_execution',
                                                 'normal_optimized_byte_identical')}, sort_keys=True))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
