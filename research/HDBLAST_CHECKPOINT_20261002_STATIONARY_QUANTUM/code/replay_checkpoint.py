#!/usr/bin/env python3
"""Replay the registered roots, original validator failure, repair, and independent checks."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

FREEZE = 'b4f77f5826e0a603400513832aa53b0673df7533'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not sys.flags.optimize, 'Run the replay driver with assertions enabled')
    package = Path(__file__).resolve().parents[1]
    repo, output = args.repo.resolve(), args.output.resolve()
    require(repo != output and repo not in output.parents,
            'Choose an unused output directory outside the repository')
    manifest = read(package / 'MANIFEST.json')
    actual = {p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()
              and p.relative_to(package).as_posix() not in manifest['excluded']}
    require(actual == set(manifest['files']), 'Checkpoint membership changed')
    for name, expected in manifest['files'].items():
        require(digest(package / name) == expected, 'Changed payload: ' + name)
    for name, expected in read(package / 'FULL_REGISTRATION.json')['files'].items():
        require(digest(package / name) == expected, 'Changed prospective source: ' + name)
    protocol = read(package / 'PREREGISTRATION.json')
    for name, expected in protocol['source_sha256'].items():
        require(digest(repo / name) == expected, 'Changed inherited input: ' + name)
    repair = read(package / 'amendments/VALIDATOR_IMPLEMENTATION_AMENDMENT.json')
    old = (package / 'code/validate_stationary.py').read_bytes()
    new = (package / 'code/validate_stationary_v2.py').read_bytes()
    require(old.replace(b"condition=row['jacobian_condition']",
                        b"condition_number=row['jacobian_condition']") == new,
            'Validator amendment is not the declared one-keyword repair')
    require(digest(package / 'code/validate_stationary_v2.py') == repair['replacement_validator_sha256'],
            'Changed amended validator')
    output.mkdir(parents=True, exist_ok=False)
    work = output / 'checkpoint'
    shutil.copytree(package, work)
    fresh, logs = output / 'fresh', output / 'logs'
    fresh.mkdir()
    logs.mkdir()
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               MPLBACKEND='Agg')
    registration = work / 'independent/REGISTRATION.md'
    ind_args = ['--registration', str(registration), '--registration-sha256', digest(registration)]
    commands = [
        ('theory', ['theory/verify_stationary_model.py', '--output', str(fresh / 'THEORY.json')], 0),
        ('theory_optimized', ['-O', 'theory/verify_stationary_model.py', '--output', str(fresh / 'THEORY_OPTIMIZED.json')], 0),
        ('primary', ['code/stationary_shell.py', '--repo', str(repo), '--public-freeze', FREEZE,
                     '--output', str(fresh / 'primary')], 0),
        ('original_validator', ['code/validate_stationary.py', '--repo', str(repo),
                                '--runs', str(fresh / 'primary/results.json'),
                                '--output', str(fresh / 'ORIGINAL_CHECKS.json')], 1),
        ('amended_validator', ['code/validate_stationary_v2.py', '--repo', str(repo),
                               '--runs', str(fresh / 'primary/results.json'),
                               '--output', str(fresh / 'CHECKS_V2.json')], 0),
        ('amended_validator_optimized', ['-O', 'code/validate_stationary_v2.py', '--repo', str(repo),
                                         '--runs', str(fresh / 'primary/results.json'),
                                         '--output', str(fresh / 'CHECKS_V2_OPTIMIZED.json')], 0),
        ('independent_roots', ['independent/independent_stationary.py', *ind_args,
                               '--amplitudes=0,100,1000000', '--slopes=-1,1',
                               '--reference=0.00011848029588657912', '--eta-ref=-0.00008405268308670538',
                               '--output', str(fresh / 'ROOTS.json')], 0),
        ('independent_proper_time', ['independent/check_proper_time.py', *ind_args,
                                     '--roots', str(fresh / 'ROOTS.json'),
                                     '--output', str(fresh / 'PROPER_TIME.json')], 0),
        ('independent_audit', ['independent/audit_primary_v2.py', '--checkpoint', str(work),
                               '--primary', str(fresh / 'primary/results.json'),
                               '--roots', str(fresh / 'ROOTS.json'), '--output', str(fresh / 'AUDIT.json')], 0),
        ('summary', ['code/summarize_stationary.py', '--runs', str(fresh / 'primary/results.json'),
                      '--checks', str(fresh / 'CHECKS_V2.json'), '--output', str(fresh / 'summary')], 0),
        ('figure', ['code/plot_stationary.py', '--summary', str(fresh / 'summary/SUMMARY.json'),
                     '--output', str(fresh / 'figures')], 0),
    ]
    receipts = []
    for name, command, expected in commands:
        with (logs / (name + '.log')).open('w') as stream:
            run = subprocess.run([sys.executable, *command], cwd=work, env=env,
                                 stdout=stream, stderr=subprocess.STDOUT)
        receipts.append({'check': name, 'command': [sys.executable, *command],
                         'exit_code': run.returncode, 'expected_exit_code': expected})
        (output / 'EXECUTION.json').write_text(json.dumps(receipts, indent=2) + '\n')
        require(run.returncode == expected, f'{name}: unexpected exit; see {logs / (name + ".log")}')
        print(name + (' preserved expected FAIL' if expected else ' PASS'), flush=True)
    original = read(fresh / 'ORIGINAL_CHECKS.json')
    failed = [g for g in original['gates'] if g['pass'] is False]
    require(original['status'] == 'FAIL' and len(failed) == 48
            and all(g['name'].endswith('/validation_exception')
                    and "multiple values for argument 'condition'" in g['traceback'] for g in failed),
            'Original validator did not reproduce precisely the documented failure')
    checks = read(fresh / 'CHECKS_V2.json')
    optimized = read(fresh / 'CHECKS_V2_OPTIMIZED.json')
    require(checks['status'] == optimized['status'] == 'PASS'
            and len(checks['gates']) == 1584 and checks['gates'] == optimized['gates']
            and checks['run_count'] == 48 and not checks['failures']
            and all(g['pass'] is True for g in checks['gates']), 'Amended checks incomplete')
    primary = read(fresh / 'primary/results.json')
    require(primary['execution_status'] == 'COMPLETE' and len(primary['runs']) == 48
            and len(primary['negative_controls']) == 3 and not primary['failures'], 'Root coverage incomplete')
    for name in ['THEORY.json', 'THEORY_OPTIMIZED.json']:
        result = read(fresh / name)
        require(result['status'] == 'PASS' and len(result['checks']) == 51
                and len(result['negative_controls']) == 20, 'Symbolic counts incomplete')
    roots = read(fresh / 'ROOTS.json')
    require(roots['status'] == 'PASS' and len(roots['runs']) == 12
            and len(roots['refinement_gates']) == 6 and not roots['failures'], 'Independent roots incomplete')
    proper = read(fresh / 'PROPER_TIME.json')
    require(proper['status'] == 'PASS' and len(proper['results']) == 8 and len(proper['gates']) == 48
            and all(g['passed'] is True for g in proper['gates']), 'Proper-time checks incomplete')
    audit = read(fresh / 'AUDIT.json')
    require(audit['status'] == 'PASS' and len(audit['gates']) == 564
            and audit['primary_rows_checked'] == 48 and audit['independent_endpoints_checked'] == 6
            and all(g['passed'] is True for g in audit['gates']), 'Independent audit incomplete')
    receipt = {'status': 'PASS', 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'original_validator_status': 'FAIL preserved', 'amended_validator_status': 'PASS',
               'scope': 'registered stationary quantum-action model; sampled hierarchy reported separately',
               'dynamical_evolution_executed': False, 'eft_certified': False,
               'manifest_sha256': digest(package / 'MANIFEST.json'), 'commands': receipts,
               'counts': {'primary_roots': 48, 'wrong_model_controls': 3, 'amended_gates': 1584,
                          'preserved_original_exceptions': 48, 'symbolic_checks': 51,
                          'symbolic_controls': 20, 'independent_roots': 12,
                          'proper_time_evaluations': 8, 'proper_time_gates': 48,
                          'independent_audit_gates': 564}}
    (output / 'VALIDATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'counts': receipt['counts'], 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
