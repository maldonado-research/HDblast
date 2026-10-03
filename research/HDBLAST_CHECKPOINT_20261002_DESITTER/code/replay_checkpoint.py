#!/usr/bin/env python3
"""Verify payload and input pins, then execute fresh sources and sensitivities."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


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
    require(not sys.flags.optimize, 'Run the replay driver without -O or PYTHONOPTIMIZE')
    package = Path(__file__).resolve().parents[1]
    repo, output = args.repo.resolve(), args.output.resolve()
    require(repo != output and repo not in output.parents,
            'Use a fresh output directory outside the repository')
    manifest = read(package / 'MANIFEST.json')
    excluded = set(manifest['excluded'])
    actual = {p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file()
              and p.relative_to(package).as_posix() not in excluded}
    require(actual == set(manifest['files']), 'Checkpoint membership differs from manifest')
    for name, expected in manifest['files'].items():
        require(digest(package / name) == expected, 'Changed payload: ' + name)
    registration = read(package / 'PREREGISTRATION.json')
    inherited = dict(registration['source_sha256'])
    for name, pin in read(package / 'sensitivity/INPUT_PINS.json').items():
        if name == 'registration':
            require(digest(package / 'sensitivity/REGISTRATION.md') == pin['sha256'],
                    'Changed sensitivity registration')
        else:
            require(name not in inherited or inherited[name] == pin['sha256'],
                    'Conflicting inherited pins')
            inherited[name] = pin['sha256']
    for name, expected in inherited.items():
        require(digest(repo / name) == expected, 'Changed inherited input: ' + name)
    for name, expected in registration['local_sha256'].items():
        require(digest(package / name) == expected, 'Changed primary registration input: ' + name)
    output.mkdir(parents=True, exist_ok=False)
    work = output / 'checkpoint'
    shutil.copytree(package, work)
    logs, fresh = output / 'logs', output / 'fresh'
    logs.mkdir()
    fresh.mkdir()
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    commands = [
        ('theory', ['theory/verify_desitter_action.py', '--output', str(fresh / 'THEORY.json')]),
        ('theory_optimized', ['-O', 'theory/verify_desitter_action.py', '--output', str(fresh / 'THEORY_OPTIMIZED.json')]),
        ('primary', ['code/desitter_sources.py', '--repo', str(repo), '--output', str(fresh / 'primary')]),
        ('primary_validator', ['code/validate_sources.py', str(fresh / 'primary')]),
        ('proper_time', ['independent/independent_proper_time.py', '--output', str(fresh / 'PROPER_TIME.json')]),
        ('exact_modes', ['independent/exact_mode_bridge.py', '--output', str(fresh / 'EXACT_MODES.json')]),
        ('exact_modes_optimized', ['-O', 'independent/exact_mode_bridge.py', '--output', str(fresh / 'EXACT_MODES_OPTIMIZED.json')]),
        ('independent_audit', ['independent/audit_independent_sources.py', '--producer', str(fresh / 'primary/results.json'), '--proper', str(fresh / 'PROPER_TIME.json'), '--output', str(fresh / 'AUDIT.json')]),
        ('summary', ['code/summarize_sources.py', '--primary', str(fresh / 'primary/results.json'), '--independent', str(fresh / 'PROPER_TIME.json'), '--output', str(fresh / 'summary')]),
        ('sensitivity', ['sensitivity/compute_sensitivity.py', '--repo', str(repo), '--output', str(fresh / 'RUNS.json')]),
        ('sensitivity_validator', ['sensitivity/verify_sensitivity.py', '--runs', str(fresh / 'RUNS.json'), '--output', str(fresh / 'SENSITIVITY.json')]),
    ]
    receipts = []
    for name, command in commands:
        with (logs / (name + '.log')).open('w') as stream:
            run = subprocess.run([sys.executable, *command], cwd=work, env=env,
                                 stdout=stream, stderr=subprocess.STDOUT)
        receipts.append({'check': name, 'command': [sys.executable, *command],
                         'exit_code': run.returncode})
        (output / 'EXECUTION.json').write_text(json.dumps(receipts, indent=2) + '\n')
        require(run.returncode == 0, f'{name} failed; see {logs / (name + ".log")}')
        print(name + ' PASS', flush=True)
    theory = read(fresh / 'THEORY.json')
    require(theory == read(fresh / 'THEORY_OPTIMIZED.json') and theory['assertions'] == 37
            and theory['negative_controls_detected'] == 6, 'Incomplete symbolic checks')
    primary = read(fresh / 'primary/results.json')
    require(primary['status'] == 'PASS' and len(primary['points']) == 4
            and len(primary['gates']) == 57 and len(primary['negative_controls']) == 3,
            'Incomplete primary checks')
    proper = read(fresh / 'PROPER_TIME.json')
    require(proper['status'] == 'PASS' and len(proper['results']) == 8
            and len(proper['refinement_gates']) == 12, 'Incomplete proper-time checks')
    exact = read(fresh / 'EXACT_MODES.json')
    require(exact == read(fresh / 'EXACT_MODES_OPTIMIZED.json') and exact['status'] == 'PASS'
            and exact['pressure_inferred_from_trace'] is False
            and [exact['results'][q]['integral_times_pi_squared'] for q in ['rho','p','Q']]
            == ['11/960', '-11/960', '1/12'], 'Incomplete exact pressure check')
    audit = read(fresh / 'AUDIT.json')
    require(audit['status'] == 'PASS' and len(audit['cross_method_gates']) == 12
            and len(audit['post_run_proper_time_trace_gates']) == 8
            and len(audit['validator_negative_controls']) == 5, 'Incomplete independent audit')
    sensitivity = read(fresh / 'SENSITIVITY.json')
    require(sensitivity['status'] == 'PASS' and sensitivity['run_count'] == 19
            and len(sensitivity['gates']) == 11
            and all(v is True for v in sensitivity['gates'].values())
            and len(sensitivity['negative_controls']) == 3 and not sensitivity['failures'],
            'Incomplete stationary sensitivity checks')
    receipt = {'status': 'PASS', 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'scope': 'fresh mathematical de Sitter sources and classical differential response',
               'coupled_evolution_executed': False, 'total_error_enclosure': False,
               'commands': receipts, 'manifest_sha256': digest(package / 'MANIFEST.json'),
               'inherited_input_count': len(inherited),
               'counts': {'primary_points': 4, 'primary_gates': 57, 'primary_negative_controls': 3,
                          'symbolic_assertions': 37, 'symbolic_negative_controls': 6,
                          'proper_time_evaluations': 8, 'proper_time_refinement_gates': 12,
                          'cross_method_gates': 12, 'post_run_trace_gates': 8,
                          'validator_mutations_rejected': 5, 'sensitivity_runs': 19,
                          'sensitivity_gates': 11, 'sensitivity_negative_controls': 3}}
    (output / 'VALIDATION.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'output': str(output), 'counts': receipt['counts']}, indent=2))


if __name__ == '__main__':
    main()
