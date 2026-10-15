"""Authenticated fresh finite checks; not computer certification of the analytic proof."""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import time

CAPSULE = 'HDBLAST_PHYSICAL_SHELL_RESPONSE_20261015'
CHECKERS = [
    {'name': 'observable', 'script': 'check_physical_observable.py',
     'golden': 'OBSERVABLE_CHECKS_NORMAL.json', 'groups': 7, 'switch': '--fault',
     'controls': [
         ('drop_embedding', 'physical scalar first variation from an explicit moving pullback'),
         ('wrong_gauge_sign', 'bulk normal gauge shift cancels the shell embedding shift'),
         ('half_hessian', 'second derivative of squared physical response includes its factor two'),
         ('wrong_angular_norm', 'actual unit S3 quadrupole norm from Hopf-coordinate integration'),
         ('missing_proper_time', 'one physical shell and proper-time integral give the rho_b power four'),
         ('second_jet_leak', 'second derivative of squared physical response includes its factor two'),
         ('double_shell', 'one physical shell and proper-time integral give the rho_b power four'),
     ]},
    {'name': 'astra', 'script': 'astra/check_physical_response.py',
     'golden': 'astra/CHECK_NORMAL.json', 'groups': 17, 'switch': '--corrupt',
     'controls': [
         ('reverse_spatial_drift', 'registered_master_in_inward_depth'),
         ('reverse_observation', 'physical_scalar_observation'),
         ('drop_normal_lapse', 'complete_scalar_junction_including_normal_variation'),
         ('drop_embedding_shift', 'physical_pullback_gauge_cancellation'),
         ('insert_measure_variation', 'exact_observable_second_variation'),
         ('double_count_common_shell', 'one_common_shell_observable_normalization'),
     ]},
    {'name': 'sol', 'script': 'sol/check_response_algebra.py',
     'golden': 'sol/normal-positive.stdout.json', 'groups': 17, 'switch': '--inject-error',
     'controls': [
         ('potential', 'full_drift_conjugation'),
         ('shell-sign', 'shell_conjugation_sign_and_mass'),
         ('ambient-velocity', 'ambient_to_global_initial_velocity'),
     ]},
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def authenticate(root, expected):
    """Bootstrap the entire manifest before importing or executing local code."""
    root = Path(root).absolute()
    require(not any(p.is_symlink() for p in [root, *root.parents]), 'Real capsule path required')
    require(re.fullmatch('[0-9a-f]{64}', expected) is not None, 'Complete external manifest SHA256 required')
    manifest = root / 'MANIFEST.json'
    require(not manifest.is_symlink() and hashlib.sha256(manifest.read_bytes()).hexdigest() == expected,
            'External manifest pin differs')
    data = json.loads(manifest.read_text())
    require(data['capsule'] == CAPSULE, 'Capsule identity differs')
    seen = set()
    for row in data['files']:
        name = row['path']
        relative = PurePosixPath(name)
        require(name and '\\' not in name and not relative.is_absolute() and '..' not in relative.parts
                and relative.as_posix() == name and name not in seen and name != 'MANIFEST.json',
                'Noncanonical, reserved or duplicate source path')
        seen.add(name)
        path = root / name
        require(not any(p.is_symlink() for p in [path, *path.parents]) and path.is_file()
                and path.resolve().is_relative_to(root.resolve()), 'Real source inside capsule required: ' + name)
        raw = path.read_bytes()
        require(len(raw) == row['bytes'] and hashlib.sha256(raw).hexdigest() == row['sha256'],
                'Frozen bytes differ: ' + name)
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Symlink in frozen tree')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    require(seen and actual == seen | {'MANIFEST.json'}, 'Frozen tree contains missing or extra files')
    return {'status': 'PASS', 'files_verified': len(seen), 'manifest_sha256': expected}


def record(path, relative_to=None):
    raw = path.read_bytes()
    return {'path': path.relative_to(relative_to).as_posix() if relative_to else str(path),
            'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def scientific(name, data):
    return {k: v for k, v in data.items() if k != 'python_optimized'} if name == 'observable' else data


def execute(root, output, checker, mode, control=None):
    work = output / mode / checker['name'] / (control[0] if control else 'positive')
    work.mkdir(parents=True)
    source = root / checker['script']
    script = work / source.name
    shutil.copy2(source, script)
    require(source.read_bytes() == script.read_bytes(), 'Fresh executable bytes differ')
    flags = ['-I', '-B'] + (['-O'] if mode == 'optimized' else [])
    arguments = [] if checker['name'] == 'sol' else ['--output', 'fresh.json']
    if control:
        arguments += [checker['switch'], control[0]]
    started = time.monotonic()
    with (work / 'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(script), *arguments], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT, timeout=180)
    diagnostic = (work / 'execution.log').read_text()
    fresh = work / 'fresh.json'
    if checker['name'] == 'sol':
        data = json.loads(diagnostic)
        fresh.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')
        output_origin = 'Complete captured checker JSON; positive stdout or failure stderr retained verbatim in execution.log'
    else:
        require(fresh.is_file(), 'Checker did not write its JSON result')
        data = json.loads(fresh.read_text())
        output_origin = 'Checker-written JSON result'
    if checker['name'] == 'observable':
        require(data['python_optimized'] is (mode == 'optimized'), 'Python optimization mode mismatch')
    row = {'kind': 'external_corruption' if control else 'positive', 'checker': checker['name'],
           'mode': mode, 'script': checker['script'], 'source_sha256': record(source)['sha256'],
           'flags': flags, 'arguments': arguments, 'exit_code': result.returncode,
           'wall_seconds': round(time.monotonic() - started, 3),
           'copied_script': record(script, output), 'log': record(work / 'execution.log', output),
           'output': record(fresh, output), 'output_origin': output_origin}
    if control:
        require(result.returncode == 1 and data.get('status', '').lower() == 'fail',
                'Injected formula did not report FAIL with exit code one')
        if checker['name'] == 'observable':
            failures = [check['check'] for check in data['checks'] if check['status'] == 'FAIL']
            require(data['fault'] == control[0] and failures == [control[1]],
                    'Observable fault rejected for an unexpected named reason')
        elif checker['name'] == 'astra':
            failures = [check['name'] for check in data['checks'] if check['passed'] is False]
            require(data['corruption'] == control[0] and failures == [control[1]],
                    'Astra corruption rejected for an unexpected named reason')
        else:
            require(isinstance(data.get('error'), str)
                    and data['error'].startswith(control[1] + ': nonzero residual '),
                    'Sol injection rejected for an unexpected named reason')
        row.update({'control': control[0], 'expected_reason': control[1],
                    'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        require(result.returncode == 0 and data.get('status', '').lower() == 'pass',
                'Fresh positive checker did not report success')
        golden = root / checker['golden']
        expected = json.loads(golden.read_text())
        require(scientific(checker['name'], data) == scientific(checker['name'], expected),
                'Scientific output differs from frozen reference: ' + checker['name'])
        if checker['name'] == 'observable':
            groups = data['finite_positive_check_groups']
            require(groups == len(data['checks']) and data['fault'] == 'none'
                    and all(c['status'] == 'PASS' for c in data['checks']), 'Observable positive inventory differs')
        elif checker['name'] == 'astra':
            groups = len(data['checks'])
            require(data['corruption'] == 'none' and all(c['passed'] is True for c in data['checks']),
                    'Astra positive inventory differs')
        else:
            groups = data['positive_identity_count']
            require(groups == len(data['identities']), 'Sol positive inventory differs')
        require(groups == checker['groups'], 'Finite positive group count differs')
        row.update({'golden': record(golden, root), 'scientific_output_matches': True,
                    'finite_positive_check_groups': groups,
                    'ignored_top_level_fields': ['python_optimized'] if checker['name'] == 'observable' else []})
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).absolute().parent
    authentication = authenticate(root, args.manifest_sha256)
    # Local verifier execution occurs only after its bytes and the full tree authenticate.
    completed = subprocess.run([sys.executable, '-I', '-S', '-B', str(root / 'verify_manifest.py'),
                                '--expected', args.manifest_sha256], capture_output=True, text=True, timeout=30)
    require(completed.returncode == 0 and json.loads(completed.stdout) == authentication,
            'Independent manifest verifier disagrees with bootstrap authentication')
    requested = args.output_dir.absolute()
    require(not any(p.is_symlink() for p in [requested, *requested.parents]), 'Real output path required')
    output = requested.resolve()
    require(not output.is_relative_to(root.resolve()), 'Output must be outside the frozen capsule')
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())),
            'New or empty output directory required')
    require(importlib.metadata.version('sympy') == '1.14.0', 'Pinned sympy 1.14.0 required')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for mode in ['normal', 'optimized']:
        for checker in CHECKERS:
            rows.append(execute(root, output, checker, mode))
        for checker in CHECKERS:
            for control in checker['controls']:
                rows.append(execute(root, output, checker, mode, control))
    positives = sum(row['kind'] == 'positive' for row in rows)
    external = sum(row['kind'] == 'external_corruption' for row in rows)
    groups = {mode: sum(row.get('finite_positive_check_groups', 0) for row in rows if row['mode'] == mode)
              for mode in ['normal', 'optimized']}
    require(positives == 6 and external == 32 and groups == {'normal': 41, 'optimized': 41},
            'Incomplete fresh replay inventory')
    require(authenticate(root, args.manifest_sha256) == authentication, 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-physical-shell-response-fresh-replay-v1', 'status': 'PASS',
               'capsule': CAPSULE, 'authentication': authentication,
               'runner': record(Path(__file__).absolute(), root),
               'integrity_verifier': record(root / 'verify_manifest.py', root),
               'requirements': record(root / 'requirements.txt', root),
               'python': {'version': sys.version, 'executable': sys.executable}, 'sympy_version': '1.14.0',
               'successful_checker_runs': positives, 'external_corruption_processes_rejected': external,
               'embedded_wrong_formula_controls_rejected': 0,
               'finite_positive_check_groups_each_mode': 41, 'finite_positive_check_groups_by_mode': groups,
               'scope': 'Exact finite algebra, pullback and normalization checks only; not certification of analytic activation, existence or nonlinear remainder bounds. External controls are actual failing formula processes. Only observable python_optimized is ignored in scientific comparison and is independently checked against each process mode.',
               'steps': rows}
    (output / 'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'successful_checker_runs': positives,
                      'external_corruption_processes_rejected': external,
                      'finite_positive_check_groups_each_mode': 41,
                      'embedded_wrong_formula_controls_rejected': 0, 'output': str(output)}, indent=2))


if __name__ == '__main__':
    main()
