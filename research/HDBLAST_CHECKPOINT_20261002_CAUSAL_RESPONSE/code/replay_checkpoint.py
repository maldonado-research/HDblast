#!/usr/bin/env python3
"""Verify and replay the self-contained registered fixed-geometry calibration."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

FREEZE = 'a01d17e015f070be2ad2018745e877148fe2f4e4'
REGISTRATION = 'f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6'
INDEPENDENT = '5bf042c7fe686e23dc469031eef81f0e682d56b4db4ff4c0381cf0d8a2ecee22'
PACKAGE = 'HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not sys.flags.optimize, 'Run the driver with optimization disabled')
    require(sys.version_info[:2] == (3, 12), 'Use Python 3.12')
    package = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    require(output != package and package not in output.parents,
            'Choose an unused output directory outside the source checkpoint')
    require(not output.exists(), 'Refusing to overwrite an existing replay')
    manifest = read(package / 'MANIFEST.json')
    excluded = {'MANIFEST.json', PACKAGE + '.zip', PACKAGE + '.zip.sha256'}
    require(set(manifest['excluded']) == excluded, 'Unexpected manifest exclusions')
    actual = {p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()
              and p.relative_to(package).as_posix() not in excluded}
    require(actual == set(manifest['files']), 'Checkpoint membership changed')
    for name, expected in manifest['files'].items():
        require(digest(package / name) == expected, 'Changed payload: ' + name)
    require(digest(package / 'FULL_REGISTRATION.json') == REGISTRATION, 'Registration changed')
    frozen = read(package / 'FULL_REGISTRATION.json')['files']
    require(len(frozen) == 32, 'Prospective coverage changed')
    for name, expected in frozen.items():
        require(digest(package / name) == expected, 'Changed frozen input: ' + name)
    require(digest(package / 'independent/MANIFEST.json') == INDEPENDENT,
            'Independent prospective manifest changed')
    for item in read(package / 'PREREGISTRATION.json')['inherited_files'].values():
        require(digest(package / item['local_copy']) == item['sha256'], 'Inherited local copy changed')
    output.mkdir(parents=True, exist_ok=False)
    work = output / 'checkpoint'
    shutil.copytree(package, work)
    fresh, logs = output / 'fresh', output / 'logs'
    fresh.mkdir()
    logs.mkdir()
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0', MPLBACKEND='Agg',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    primary, modes = fresh / 'primary/results.json', fresh / 'modes/results.json'
    validation = ['code/validate_causal.py', '--primary', str(primary), '--modes', str(modes),
                  '--registration-sha256', REGISTRATION]
    commands = [
        ('theory', ['theory/verify_causal_theory.py', '--output', str(fresh / 'THEORY.json')]),
        ('theory_optimized', ['-O', 'theory/verify_causal_theory.py', '--output', str(fresh / 'THEORY_OPTIMIZED.json')]),
        ('primary', ['code/causal_primary.py', '--output', str(fresh / 'primary'),
                     '--registration-sha256', REGISTRATION, '--public-freeze-commit', FREEZE]),
        ('independent_modes', ['independent/forced_modes.py', '--output-dir', str(fresh / 'modes'),
                               '--freeze-commit', FREEZE, '--manifest-sha256', INDEPENDENT]),
        ('validation', [*validation, '--output', str(fresh / 'CHECKS.json')]),
        ('validation_optimized', ['-O', *validation, '--output', str(fresh / 'CHECKS_OPTIMIZED.json')]),
        ('archive_audit', ['independent/audit_archived_modes.py', '--checkpoint', str(work),
                           '--primary', str(primary), '--modes', str(modes),
                           '--checks', str(fresh / 'CHECKS.json'), '--output', str(fresh / 'ARCHIVE_AUDIT.json')]),
        ('summary', ['postrun/summarize_causal.py', '--checkpoint', str(work), '--primary', str(primary),
                     '--modes', str(modes), '--checks', str(fresh / 'CHECKS.json'),
                     '--checks-optimized', str(fresh / 'CHECKS_OPTIMIZED.json'), '--output', str(fresh / 'summary')]),
        ('figures', ['postrun/plot_causal.py', '--summary', str(fresh / 'summary/SUMMARY.json'),
                     '--output', str(fresh / 'figures')]),
    ]
    receipts = []
    for name, command in commands:
        timed_out = False
        with (logs / (name + '.log')).open('w') as stream:
            try:
                result = subprocess.run([sys.executable, *command], cwd=work, env=env,
                                        stdout=stream, stderr=subprocess.STDOUT, timeout=1000)
                exit_code = result.returncode
            except subprocess.TimeoutExpired:
                timed_out, exit_code = True, None
                stream.write('\nReplay command exceeded its 1000-second wrapper limit.\n')
        receipts.append({'check': name, 'command': [sys.executable, *command],
                         'exit_code': exit_code, 'expected_exit_code': 0, 'timed_out': timed_out})
        (output / 'EXECUTION.json').write_text(json.dumps(receipts, indent=2) + '\n')
        require(exit_code == 0 and not timed_out,
                name + ': unexpected exit or timeout; see ' + str(logs / (name + '.log')))
        print(name + ' PASS', flush=True)
    for name, expected in frozen.items():
        require(digest(work / name) == digest(package / name) == expected,
                'Frozen input changed during replay: ' + name)
    theory, optimized_theory = read(fresh / 'THEORY.json'), read(fresh / 'THEORY_OPTIMIZED.json')
    require(theory == optimized_theory and theory['passed'] is True
            and theory['identity_count'] == 46 and theory['mutation_count'] == 14
            and all(c['passed'] for c in theory['checks'])
            and all(c['detected'] for c in theory['mutations']), 'Exact theory coverage failed')
    checks, optimized = read(fresh / 'CHECKS.json'), read(fresh / 'CHECKS_OPTIMIZED.json')
    require(checks['status'] == optimized['status'] == 'passed'
            and {k:v for k,v in checks.items() if k != 'evidence'}
                == {k:v for k,v in optimized.items() if k != 'evidence'}
            and len(checks['comparisons']) == 36 and len(checks['wrong_formula_controls']) == 17
            and all(c['rejected'] for c in checks['wrong_formula_controls']), 'Registered validation incomplete')
    require(len(read(primary)['rows']) == len(read(modes)['rows']) == 12
            and len(read(modes)['runs']) == 4, 'Producer coverage incomplete')
    audit = read(fresh / 'ARCHIVE_AUDIT.json')
    require(audit['status'] == 'passed'
            and audit['counts']['reconstructed_finite_K_responses'] == 72
            and audit['counts']['mode_observation_snapshots'] == 24
            and audit['counts']['panel_moment_checks'] == 98304, 'Archive audit incomplete')
    figure_files = [p for p in (fresh / 'figures').iterdir() if p.suffix in {'.svg', '.png', '.pdf'}]
    require(len(figure_files) == 9, 'Figure coverage incomplete')
    fresh_hashes = {p.relative_to(fresh).as_posix(): digest(p) for p in sorted(fresh.rglob('*')) if p.is_file()}
    receipt = {'status': 'PASS', 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'public_freeze_commit': FREEZE, 'registration_sha256': REGISTRATION,
               'manifest_sha256': digest(package / 'MANIFEST.json'), 'commands': receipts,
               'scope': 'two registered homogeneous scalar-variance responses on fixed de Sitter geometry',
               'stress_response_executed': False, 'coupled_geometry_executed': False,
               'heating_established': False, 'stability_established': False,
               'total_numerical_error_certified': False,
               'counts': {'commands': 9, 'symbolic_identities': 46, 'symbolic_mutations': 14,
                          'response_points': 12, 'finite_K_comparisons': 36, 'wrong_formula_controls': 17,
                          'independent_mode_runs': 4, 'archive_integrals': 72, 'archive_snapshots': 24,
                          'panel_moment_checks': 98304, 'figure_files': 9},
               'fresh_output_sha256': fresh_hashes}
    (output / 'VALIDATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'counts': receipt['counts']}, indent=2))


if __name__ == '__main__':
    main()
