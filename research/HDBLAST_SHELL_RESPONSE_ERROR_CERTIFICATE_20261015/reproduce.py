"""Authenticated fresh arithmetic and algebra checks, not analytic proof certification."""
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

CHECKERS = [
    {'name': 'certificate', 'script': 'check_response_certificate.py',
     'dependencies': ['certify_response.py'], 'golden': 'RESPONSE_CERTIFICATE_NORMAL.json',
     'output': None, 'groups': 25, 'embedded_controls': 0},
    {'name': 'archive', 'script': 'archive/derive_shell_box.py',
     'dependencies': [
         'archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/theory/JACOBI_ENDPOINT_MONOTONICITY.md',
         'archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/bounds/RESCALED_BOUND_CERTIFICATE.json'],
     'golden': 'archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json',
     'output': 'archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json', 'groups': 14, 'embedded_controls': 0},
    {'name': 'sol', 'script': 'sol/algebra/check_shell_response_algebra.py',
     'dependencies': [], 'golden': 'sol/algebra/symbolic_check_results.json',
     'output': 'sol/algebra/symbolic_check_results.json', 'groups': 35, 'embedded_controls': 9},
    {'name': 'astra', 'script': 'astra/check_residual_response.py',
     'dependencies': [], 'golden': 'astra/ASTRA_REPLAY_RECEIPT.json',
     'output': None, 'groups': 17, 'embedded_controls': 0},
    {'name': 'referee', 'script': 'referee/independent_coarse_response_audit.py',
     'dependencies': ['archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json'],
     'golden': 'referee/INDEPENDENT_COARSE_RESPONSE_AUDIT.json',
     'output': 'referee/INDEPENDENT_COARSE_RESPONSE_AUDIT.json', 'groups': 23, 'embedded_controls': 0},
    {'name': 'sol_crossreview', 'script': 'sol/crossreview/independent_quantitative_review.py',
     'dependencies': [
         'SHELL_RESPONSE_ERROR_CERTIFICATE.md', 'certify_response.py', 'check_response_certificate.py',
         'archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md',
         'archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json', 'archive/derive_shell_box.py',
         'archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/theory/JACOBI_ENDPOINT_MONOTONICITY.md',
         'archive/sources/git/research/HDBLAST_SCALAR_RESPONSE_MONOTONICITY_20261011/bounds/RESCALED_BOUND_CERTIFICATE.json',
         'sol/RESIDUAL_TO_PHYSICAL_RESPONSE.md'],
     'golden': 'sol/crossreview/INDEPENDENT_QUANTITATIVE_REVIEW.json',
     'output': None, 'groups': 35, 'embedded_controls': 0},
]

ASTRA_CORRUPTIONS = [
    {'name': 'boundary-source', 'reason': 'boundary flux and inertia cancellation'},
    {'name': 'moving-flux', 'reason': 'shrinking right boundary flux'},
    {'name': 'physical-prefactor', 'reason': 'physical coefficient normalization'},
]

# Each source fragment must occur exactly once in the authenticated module.
# The independent checker is never changed. A failed process must emit JSON
# with the required named reason; an exception or unrelated failure is rejected.
CORRUPTIONS = [
    {'name': 'sqrt_inward_upper', 'reason': 'sqrt_outward_2',
     'before': 'else F(integer + 1, scale)', 'after': 'else F(integer, scale)'},
    {'name': 'wrong_inertia_weight', 'reason': 'positive_inertia_forcing',
     'before': "d['G_abs_sup'] ** 2 / d['mu_lower']", 'after': "d['G_abs_sup'] ** 2 * d['mu_lower']"},
    {'name': 'omit_physical_shift', 'reason': 'physical_shift_trace',
     'before': "1 / d['mu_lower'] + d['c_abs'] ** 2 / d['nu']", 'after': "1 / d['mu_lower']"},
    {'name': 'reverse_norm_triangle', 'reason': 'norm_triangle_lower',
     'before': "sqrt_bounds(d['trial_D_sq_lower'])[0] - delta", 'after': "sqrt_bounds(d['trial_D_sq_lower'])[0] + delta"},
    {'name': 'wrong_response_prefactor', 'reason': 'proper_time_prefactor',
     'before': 'prefactor_lower = 3 * pil ** 2', 'after': 'prefactor_lower = 2 * pil ** 2'},
    {'name': 'erase_exponential_growth', 'reason': 'exp_positive_bound',
     'before': 'return (1 - value / count) ** (-count)', 'after': 'return F(1)'},
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def authenticate(root, expected):
    """Full stdlib bootstrap, before any local module import or execution."""
    root = Path(root).absolute()
    require(not any(p.is_symlink() for p in [root, *root.parents]), 'Real capsule path required')
    require(re.fullmatch('[0-9a-f]{64}', expected) is not None, 'Complete external manifest SHA256 required')
    manifest = root / 'MANIFEST.json'
    require(manifest.is_file() and not manifest.is_symlink(), 'Real regular manifest required')
    raw = manifest.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, 'External manifest pin differs')
    data = json.loads(raw)
    require(data['capsule'] == root.name, 'Capsule directory identity differs')
    seen = set()
    for row in data['files']:
        name = row['path']
        relative = PurePosixPath(name)
        require(name and '\\' not in name and not relative.is_absolute() and '..' not in relative.parts
                and relative.as_posix() == name and name not in seen and name != 'MANIFEST.json',
                'Noncanonical, reserved or duplicate source path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid source record')
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
        else:
            require(path.is_dir(), 'Special file in frozen tree')
    require(seen and actual == seen | {'MANIFEST.json'}, 'Frozen tree contains missing or extra files')
    return {'status': 'PASS', 'files_verified': len(seen), 'manifest_sha256': expected}


def record(path, base):
    raw = path.read_bytes()
    return {'path': path.relative_to(base).as_posix(), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}


def positive_inventory(checker, data):
    if checker['name'] == 'certificate':
        require(data['status'] == 'PASS' and not data['failures'], 'Certificate reported failures')
        require(len(data['groups']) == len(set(data['groups'])) == data['finite_check_groups'],
                'Certificate group inventory differs')
        return data['finite_check_groups'], 0
    if checker['name'] == 'archive':
        require(data['status'] == 'PASS_EXACT_RATIONAL_CONSEQUENCES_OF_RETAINED_ANALYTIC_THEOREM',
                'Archive certificate did not pass')
        collar = data['collar']
        require(collar['status'] == 'PASS_EXACT_FIRST_EXIT_TUBE_AND_COEFFICIENT_BOUNDS'
                and len(collar['first_exit_inclusions']) == 4
                and len(collar['gates']) == 10 and all(v is True for v in collar['gates'].values()),
                'Archive inclusion/gate inventory differs')
        return 14, 0
    if checker['name'] == 'astra':
        require(data['status'] == 'pass' and data['checks'] == len(data['names'])
                == len(set(data['names'])), 'Astra identity inventory differs')
        return data['checks'], 0
    if checker['name'] == 'referee':
        require(data['status'] == 'PASS_INDEPENDENT_COARSE_RATIONAL_CONSEQUENCES'
                and all(value is True for value in data['checks'].values()), 'Referee checks differ')
        return len(data['checks']), 0
    if checker['name'] == 'sol_crossreview':
        require(data['status'] == 'PASS_CONDITIONAL_QUANTITATIVE_CROSS_REVIEW'
                and data['check_count'] == len(data['checks']) == len(set(data['checks']))
                and data['production_module_imported'] is False
                and data['physical_waveform_computed'] is False
                and data['analytic_hypotheses_certified_by_program'] is False,
                'Independent Sol cross-review inventory or scope differs')
        return data['check_count'], 0
    require(data['positive_checks_passed'] == len(data['positive_checks'])
            and all(c['residual'] == '0' for c in data['positive_checks']), 'Sol identities differ')
    require(data['negative_controls_rejected'] == len(data['negative_controls'])
            and all(c['exact_witness_defect'] != '0' for c in data['negative_controls']),
            'Sol embedded wrong-formula witnesses differ')
    return data['positive_checks_passed'], data['negative_controls_rejected']


def execute(root, output, checker, mode, corruption=None):
    work = output / mode / checker['name'] / (corruption['name'] if corruption else 'positive')
    work.mkdir(parents=True)
    copied = []
    for name in [checker['script'], *checker['dependencies']]:
        target = work / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / name, target)
        require(target.read_bytes() == (root / name).read_bytes(), 'Fresh copy differs: ' + name)
        copied.append({'source': record(root / name, root), 'copy': record(target, output)})
    mutation = None
    if corruption and checker['name'] == 'certificate':
        target = work / 'certify_response.py'
        before = target.read_text()
        require(before.count(corruption['before']) == 1, 'Corruption must match exactly one source fragment')
        after = before.replace(corruption['before'], corruption['after'], 1)
        require(before != after, 'Corruption must change formula bytes')
        target.write_text(after)
        mutation = {**corruption, 'source_module': record(root / 'certify_response.py', root),
                    'mutated_module': record(target, output),
                    'independent_checker_unchanged': True}
    flags = ['-I', '-B'] + (['-O'] if mode == 'optimized' else [])
    arguments = ['--inject-error', corruption['name']] if corruption and checker['name'] == 'astra' else []
    script = work / checker['script']
    started = time.monotonic()
    with (work / 'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(script), *arguments], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT, timeout=180)
    if checker['name'] == 'astra' and corruption:
        diagnostic = (work / 'execution.log').read_text()
        prefix = 'CHECK_FAILED: ' + corruption['reason'] + ': nonzero discrepancy '
        require(result.returncode == 1 and diagnostic.startswith(prefix)
                and diagnostic.endswith('\n') and diagnostic.count('\n') == 1
                and diagnostic[len(prefix):].strip() not in ('', '0'),
                'Astra fault did not report the exact named nonzero discrepancy')
        data = {'status': 'FAIL', 'failures': [corruption['reason']], 'diagnostic': diagnostic,
                'result_origin': 'Runner parses the complete one-line checker failure; original text retained in execution.log'}
    elif checker['output'] is None:
        data = json.loads((work / 'execution.log').read_text())
    else:
        result_path = work / checker['output']
        require(result_path.is_file(), 'Checker result missing')
        data = json.loads(result_path.read_text())
    fresh = work / 'fresh.json'
    fresh.write_text(json.dumps(data, indent=2) + '\n')
    row = {'kind': 'external_corruption' if corruption else 'positive', 'checker': checker['name'],
           'mode': mode, 'script': checker['script'], 'flags': flags, 'arguments': arguments,
           'exit_code': result.returncode,
           'wall_seconds': round(time.monotonic() - started, 3), 'copied_inputs': copied,
           'copied_inputs_snapshot': 'Initial byte-identical copies before any separately recorded formula mutation',
           'log': record(work / 'execution.log', output), 'output': record(fresh, output),
           'output_origin': ('Runner-parsed complete CHECK_FAILED line' if checker['name'] == 'astra' and corruption
                             else checker['output'] or 'complete stdout JSON, retained verbatim in execution.log')}
    require(script.read_bytes() == (root / checker['script']).read_bytes(), 'Independent checker changed')
    if corruption:
        require(result.returncode == 1 and data.get('status') == 'FAIL'
                and corruption['reason'] in data.get('failures', []),
                'Fault did not report the expected named failure with exit one')
        row.update({'mutation': mutation, 'control': corruption['name'],
                    'expected_reason': corruption['reason'], 'reported_failures': data['failures'],
                    'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        require(result.returncode == 0, 'Positive checker exited unsuccessfully')
        golden = root / checker['golden']
        expected = json.loads(golden.read_text())
        if checker['name'] == 'astra':
            positives = [r for r in expected['runs'] if 'kind' not in r
                         and r['optimized'] is False and r['mutation'] is None]
            require(len(positives) == 1 and positives[0]['exit_code'] == 0, 'Astra golden selection differs')
            expected = json.loads(positives[0]['stdout'])
        require(data == expected, 'Fresh scientific JSON differs from frozen golden')
        groups, embedded = positive_inventory(checker, data)
        require(groups == checker['groups'] and embedded == checker['embedded_controls'],
                'Finite check inventory differs')
        row.update({'golden': record(golden, root), 'scientific_output_matches': True,
                    'finite_positive_check_groups': groups, 'embedded_wrong_formula_controls': embedded,
                    'ignored_top_level_fields': []})
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).absolute().parent
    authentication = authenticate(root, args.manifest_sha256)
    verified = subprocess.run([sys.executable, '-I', '-S', '-B', str(root / 'verify_manifest.py'),
                               '--expected', args.manifest_sha256], capture_output=True, text=True, timeout=30)
    require(verified.returncode == 0 and json.loads(verified.stdout) == authentication,
            'Independent verifier disagrees with bootstrap')
    requested = args.output_dir.absolute()
    require(not any(p.is_symlink() for p in [requested, *requested.parents]), 'Real output path required')
    output = requested.resolve()
    require(not output.is_relative_to(root.resolve()), 'Output must be outside capsule')
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())), 'New or empty output required')
    require(importlib.metadata.version('sympy') == '1.14.0', 'Pinned SymPy 1.14.0 required')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for mode in ['normal', 'optimized']:
        for checker in CHECKERS:
            rows.append(execute(root, output, checker, mode))
        for corruption in CORRUPTIONS:
            rows.append(execute(root, output, CHECKERS[0], mode, corruption))
        astra = next(c for c in CHECKERS if c['name'] == 'astra')
        for corruption in ASTRA_CORRUPTIONS:
            rows.append(execute(root, output, astra, mode, corruption))
    positives = sum(r['kind'] == 'positive' for r in rows)
    external = sum(r['kind'] == 'external_corruption' for r in rows)
    embedded = sum(r.get('embedded_wrong_formula_controls', 0) for r in rows)
    groups = {mode: sum(r.get('finite_positive_check_groups', 0) for r in rows if r['mode'] == mode)
              for mode in ['normal', 'optimized']}
    total_groups = sum(c['groups'] for c in CHECKERS)
    require(positives == 2 * len(CHECKERS) and external == 2 * (len(CORRUPTIONS) + len(ASTRA_CORRUPTIONS))
            and embedded == 2 * sum(c['embedded_controls'] for c in CHECKERS)
            and groups == {'normal': total_groups, 'optimized': total_groups}, 'Incomplete replay inventory')
    require(authenticate(root, args.manifest_sha256) == authentication, 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-validated-response-fresh-replay-v1', 'status': 'PASS',
               'capsule': root.name, 'authentication': authentication,
               'runner': record(Path(__file__).absolute(), root),
               'integrity_verifier': record(root / 'verify_manifest.py', root),
               'requirements': record(root / 'requirements.txt', root),
               'python': {'version': sys.version, 'executable': sys.executable}, 'sympy_version': '1.14.0',
               'successful_checker_runs': positives, 'external_corruption_processes_rejected': external,
               'embedded_wrong_formula_controls_rejected': embedded,
               'finite_positive_check_groups_each_mode': total_groups,
               'finite_positive_check_groups_by_mode': groups,
               'scope': 'Fresh exact finite arithmetic/algebra checks and named wrong-formula processes only. Analytic energy/trace estimates, retained background existence and nonlinear remainder statements remain mathematical hypotheses and arguments; no waveform solution or formal proof certification. Complete scientific JSON compared with no ignored fields.',
               'steps': rows}
    (output / 'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'successful_checker_runs',
                     'external_corruption_processes_rejected', 'embedded_wrong_formula_controls_rejected',
                     'finite_positive_check_groups_each_mode']}, indent=2))


if __name__ == '__main__':
    main()
