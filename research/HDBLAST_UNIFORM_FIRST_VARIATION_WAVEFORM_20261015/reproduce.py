"""Authenticated fresh interval calculations, separate from analytic PDE proof."""
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

CAPSULE = 'HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015'
GAIN = 'archive/PHYSICAL_GAIN_ENCLOSURE.json'
BACKGROUND = 'inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json'
CHECKERS = [
    {'name': 'root', 'script': 'check_waveform.py',
     'dependencies': ['waveform_intervals.py', 'compute_waveform.py', GAIN],
     'golden': 'CHECK_WAVEFORM_NORMAL.json', 'output': None, 'groups': 44, 'embedded': 0},
    {'name': 'generator', 'script': 'compute_waveform.py',
     'dependencies': ['waveform_intervals.py', GAIN],
     'arguments': ['--output-dir', 'generated'],
     'golden': 'computed/WAVEFORM_ENCLOSURE.json', 'output': 'generated/WAVEFORM_ENCLOSURE.json',
     'additional_outputs': [('generated/WAVEFORM_CELLS.csv', 'computed/WAVEFORM_CELLS.csv')],
     'groups': 0, 'embedded': 0},
    {'name': 'gain', 'script': 'archive/enclose_physical_gain.py',
     'dependencies': ['archive/inputs/' + n for n in [
         'OCT11_JACOBI_ENDPOINT_MONOTONICITY.md', 'OCT11_RESCALED_BOUND_CERTIFICATE.json',
         'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json', 'REFERENCE_ONLY_prior_derive_shell_box.py',
         'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md',
         'PRIOR_SHELL_RESPONSE_ERROR_CERTIFICATE.md']],
     'golden': GAIN, 'output': GAIN, 'groups': 0, 'embedded': 0},
    {'name': 'sol', 'script': 'sol/check_uniform_waveform.py',
     'dependencies': ['sol/UNIFORM_WAVEFORM_CERTIFICATE.md',
         'inputs/previous/SHELL_RESPONSE_ERROR_CERTIFICATE.md',
         'inputs/previous/sol/RESIDUAL_TO_PHYSICAL_RESPONSE.md',
         'inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md', BACKGROUND],
     'golden': 'sol/UNIFORM_WAVEFORM_CHECKS.json',
     'optimized_golden': 'sol/UNIFORM_WAVEFORM_CHECKS_OPTIMIZED.json',
     'output': 'sol/UNIFORM_WAVEFORM_CHECKS.json', 'groups': 34, 'embedded': 0},
    {'name': 'astra', 'script': 'astra/check_uniform_pointwise.py', 'dependencies': [],
     'golden': 'astra/UNIFORM_POINTWISE_NORMAL.json', 'output': None, 'groups': 14, 'embedded': 15},
    {'name': 'trace', 'script': 'referee/independent_uniform_trace_audit.py',
     'dependencies': [BACKGROUND], 'arguments': ['--background', BACKGROUND],
     'golden': 'referee/INDEPENDENT_UNIFORM_TRACE_AUDIT_NORMAL.json',
     'output': None, 'groups': 20, 'embedded': 0},
    {'name': 'csv', 'script': 'referee/independent_csv_audit.py',
     'dependencies': ['computed/WAVEFORM_ENCLOSURE.json', 'computed/WAVEFORM_CELLS.csv',
                      'waveform_intervals.py', 'compute_waveform.py', GAIN],
     'golden': 'referee/INDEPENDENT_CSV_AND_NORM_AUDIT.json',
     'output': None, 'groups': 23, 'embedded': 0},
    {'name': 'peak', 'script': 'astra/peak_analysis/verify_peaks.py', 'dependencies': [],
     'golden': 'astra/peak_analysis/certificate.json', 'output': None, 'groups': 0, 'embedded': 4},
]
# Only the copied production module is changed. The checker remains byte-identical.
CORRUPTIONS = [
    {'name': 'round_upper_inward', 'module': 'waveform_intervals.py', 'reason': 'negative_ceiling',
     'before': 'return -((-n) // d)', 'after': 'return n // d'},
    {'name': 'erase_zero_crossing', 'module': 'waveform_intervals.py', 'reason': 'square_through_zero',
     'before': 'low = 0 if a[0] <= 0 <= a[1] else min(a[0] ** 2, a[1] ** 2)',
     'after': 'low = min(a[0] ** 2, a[1] ** 2)'},
    {'name': 'discard_exponential_tail', 'module': 'waveform_intervals.py', 'reason': 'exp_bracket_128',
     'before': 'return 0, 1 << (BITS - exponent)', 'after': 'return 0, 0'},
    {'name': 'reverse_bump_derivative', 'module': 'waveform_intervals.py', 'reason': 'rising_derivative_positive',
     'before': 'derivative = mul(interval(3 - 2 * right, 3 - 2 * left), (elo, ehi))',
     'after': 'derivative = neg(mul(interval(3 - 2 * right, 3 - 2 * left), (elo, ehi)))'},
    {'name': 'omit_majorant_coupling', 'module': 'waveform_intervals.py', 'reason': 'majorant_coupled_rows',
     'before': 'determinant = (1 - r) * (1 - s) - T * r / m',
     'after': 'determinant = (1 - r) * (1 - s)'},
    {'name': 'omit_uniform_clock_term', 'module': 'waveform_intervals.py', 'reason': 'uniform_output_sources',
     'before': 'eta = ell * (M * T * A + (1 + gamma * T) * V)',
     'after': 'eta = ell * (M * T * A + V)'},
    {'name': 'omit_tail_inertia_term', 'module': 'waveform_intervals.py', 'reason': 'tail_output_sources',
     'before': 'V_tail = ell / m ** 2 + ell / m + r * A / m + s * V',
     'after': 'V_tail = ell / m + r * A / m + s * V'},
    {'name': 'wrong_coefficient_scale', 'module': 'compute_waveform.py', 'reason': 'coefficient_scaling',
     'before': 'C=(K[0]*true_integral[0]/ell,K[1]*true_integral[1]/ell)',
     'after': 'C=(K[0]*true_integral[0]*ell,K[1]*true_integral[1]/ell)'},
]

SOL_CORRUPTIONS = [
    {'name': 'drop-determinant-coupling', 'reason': 'first majorant equation'},
    {'name': 'flip-boundary-damping', 'reason': 'boundary jet correct sign identity'},
    {'name': 'use-full-C', 'reason': 'wrong gamma equals full C rejected'},
    {'name': 'omit-velocity-input', 'reason': 'dropped velocity input rejected at peak'},
    {'name': 'delete-boundary-memory', 'reason': 'deleted boundary memory rejected'},
]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def authenticate(root, expected):
    """Stdlib bootstrap, before any capsule module import or execution."""
    root = Path(root).absolute()
    require(not any(p.is_symlink() for p in [root, *root.parents]), 'Real capsule path required')
    require(re.fullmatch('[0-9a-f]{64}', expected) is not None, 'Complete external manifest SHA256 required')
    manifest = root / 'MANIFEST.json'
    require(manifest.is_file() and not manifest.is_symlink(), 'Real regular manifest required')
    raw = manifest.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, 'External manifest pin differs')
    data = json.loads(raw)
    require(data['capsule'] == root.name == CAPSULE, 'Capsule directory identity differs')
    records = {}
    for row in data['files']:
        name = row['path']
        relative = PurePosixPath(name)
        require(name and '\\' not in name and not relative.is_absolute() and '..' not in relative.parts
                and relative.as_posix() == name and name not in records and name != 'MANIFEST.json',
                'Noncanonical, reserved or duplicate source path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid source record')
        records[name] = row
        read_source(root, records, name)
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Symlink in frozen tree')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            require(path.is_dir(), 'Special file in frozen tree')
    require(records and actual == set(records) | {'MANIFEST.json'}, 'Frozen tree contains missing or extra files')
    return {'status': 'PASS', 'files_verified': len(records), 'manifest_sha256': expected}, records


def read_source(root, inventory, name):
    require(name in inventory, 'Missing declared source dependency: ' + name)
    path = root / name
    require(not any(p.is_symlink() for p in [path, *path.parents]) and path.is_file()
            and path.resolve().is_relative_to(root.resolve()), 'Real source inside capsule required: ' + name)
    raw = path.read_bytes()
    item = inventory[name]
    require(len(raw) == item['bytes'] and hashlib.sha256(raw).hexdigest() == item['sha256'],
            'Frozen bytes differ: ' + name)
    return raw


def record(path, base):
    raw = path.read_bytes()
    return {'path': path.relative_to(base).as_posix(), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}


def positive_inventory(checker, data, mode):
    name = checker['name']
    if name in ('root', 'sol', 'astra'):
        require(data['status'] == 'PASS' and len(data['checks']) == len(set(data['checks'])),
                'Unique successful check names required')
        groups = len(data['checks'])
        if name == 'root':
            require(groups == data['finite_check_groups'], 'Root group count differs')
        if name == 'sol':
            require(groups == data['check_count'] and data['python_optimization'] == (mode == 'optimized'),
                    'Sol group count or optimization flag differs')
        controls = data.get('wrong_formula_controls_rejected', [])
        require(len(controls) == len(set(controls)), 'Repeated embedded control')
        return groups, len(controls)
    if name in ('trace', 'csv'):
        expected = ('PASS_INDEPENDENT_UNIFORM_CHARACTERISTIC_TRACE_BOUND' if name == 'trace'
                    else 'PASS_INDEPENDENT_CSV_AND_NORM_AUDIT')
        require(data['status'] == expected and all(v is True for v in data['checks'].values())
                and len(data['checks']) == data['finite_check_groups'], 'Independent audit failed')
        return len(data['checks']), 0
    if name == 'generator':
        require(data['status'] == 'PASS' and data['quadrature']['support_cells'] == 4096
                and data['cell_data']['rows'] == 4098, 'Generator grid or status differs')
        return 0, 0
    if name == 'gain':
        require(data['status'] == 'PASS_EXACT_RATIONAL_PHYSICAL_GAIN_ENCLOSURE'
                and data['archived_programs_executed'] is False
                and len(data['polynomial_identities_verified']) == 3, 'Gain enclosure inventory differs')
        return 0, 0
    require(name == 'peak' and data['status'] == 'PASS' and len(data['enclosures']) == 14
            and len(data['negative_controls']) == 4
            and len({c['name'] for c in data['negative_controls']}) == 4
            and all(c['rejected'] is True for c in data['negative_controls']), 'Peak enclosure inventory differs')
    return 0, 4


def execute(root, inventory, output, checker, mode, corruption=None):
    work = output / mode / checker['name'] / (corruption['name'] if corruption else 'positive')
    work.mkdir(parents=True)
    copied = []
    for name in [checker['script'], *checker['dependencies']]:
        raw = read_source(root, inventory, name)
        target = work / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        require(target.read_bytes() == raw, 'Fresh copy differs: ' + name)
        copied.append({'source': inventory[name], 'copy': record(target, output)})
    mutation = None
    if corruption and checker['name'] == 'root':
        target = work / corruption['module']
        before = target.read_text()
        require(before.count(corruption['before']) == 1, 'Mutation must match exactly one source fragment')
        after = before.replace(corruption['before'], corruption['after'], 1)
        require(before != after, 'Mutation must change formula bytes')
        target.write_text(after)
        mutation = {**corruption, 'source_module': inventory[corruption['module']],
                    'mutated_module': record(target, output), 'independent_checker_unchanged': True}
    flags = ['-I', '-B'] + (['-O'] if mode == 'optimized' else [])
    arguments = (['--inject-error', corruption['name']] if corruption and checker['name'] == 'sol'
                 else checker.get('arguments', []))
    result_path = work / checker['output'] if checker['output'] else work / 'execution.log'
    if checker['output']:
        require(not result_path.exists(), 'Stale output present before execution')
    started = time.monotonic()
    with (work / 'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(work / checker['script']), *arguments],
                                cwd=work, stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if corruption and checker['name'] == 'sol':
        diagnostic = (work / 'execution.log').read_text()
        reason = 'RuntimeError: failed check: ' + corruption['reason']
        require(result.returncode == 1 and not result_path.exists()
                and diagnostic.startswith('Traceback (most recent call last):\n')
                and diagnostic.rstrip().splitlines()[-1] == reason
                and 'PASS:' not in diagnostic,
                'Sol fault did not reject with its exact named failure and no passing output')
        data = {'status': 'FAIL', 'failure': corruption['reason'], 'exception': 'RuntimeError',
                'result_origin': 'Runner parses final named RuntimeError; complete traceback retained in log',
                'passing_output_absent': True}
        raw = (json.dumps(data, indent=2) + '\n').encode()
    else:
        require(result_path.is_file(), 'Fresh result missing')
        raw = result_path.read_bytes()
        data = json.loads(raw)
    fresh = work / 'fresh.json'
    fresh.write_bytes(raw)
    row = {'kind': 'external_corruption' if corruption else 'positive', 'checker': checker['name'],
           'mode': mode, 'script': checker['script'], 'flags': flags, 'arguments': arguments,
           'exit_code': result.returncode, 'wall_seconds': round(time.monotonic() - started, 3),
           'copied_inputs': copied, 'copied_inputs_snapshot': 'Authenticated initial copies before the recorded formula mutation',
           'log': record(work / 'execution.log', output), 'output': record(fresh, output),
           'additional_outputs': [], 'output_origin': ('Runner-parsed named RuntimeError traceback'
            if corruption and checker['name'] == 'sol' else checker['output'] or 'complete stdout JSON')}
    require((work / checker['script']).read_bytes() == read_source(root, inventory, checker['script']),
            'Independent checker script changed')
    if corruption:
        require(result.returncode == 1 and data.get('status') == 'FAIL'
                and data.get('failure') == corruption['reason'] and data.get('exception') == 'RuntimeError',
                'Fault did not report its exact named assertion with exit one: ' + corruption['name'])
        row.update({'mutation': mutation, 'control_mechanism': 'copied module formula replacement' if mutation else 'checker CLI formula injection',
                    'control': corruption['name'],
                    'expected_reason': corruption['reason'], 'reported_failure': data['failure'],
                    'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        require(result.returncode == 0, 'Positive checker exited unsuccessfully: ' + checker['name'])
        golden_name = checker.get('optimized_golden', checker['golden']) if mode == 'optimized' else checker['golden']
        golden_raw = read_source(root, inventory, golden_name)
        require(data == json.loads(golden_raw), 'Full fresh JSON differs from frozen golden: ' + checker['name'])
        if checker['name'] in ('generator', 'gain'):
            require(raw == golden_raw, 'Complete generated JSON bytes differ')
        groups, embedded = positive_inventory(checker, data, mode)
        require(groups == checker['groups'] and embedded == checker['embedded'], 'Finite inventory differs')
        for path, golden in checker.get('additional_outputs', []):
            fresh_path = work / path
            require(fresh_path.is_file() and fresh_path.read_bytes() == read_source(root, inventory, golden),
                    'Complete generated data bytes differ: ' + path)
            row['additional_outputs'].append({'output': record(fresh_path, output),
                                              'golden': inventory[golden], 'complete_bytes_equal': True})
        row.update({'golden': inventory[golden_name], 'scientific_output_matches': True,
                    'complete_generated_json_bytes_equal': checker['name'] in ('generator', 'gain'),
                    'finite_positive_check_groups': groups, 'embedded_wrong_formula_controls': embedded,
                    'ignored_golden_fields': []})
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).absolute().parent
    authentication, inventory = authenticate(root, args.manifest_sha256)
    verified = subprocess.run([sys.executable, '-I', '-S', '-B', str(root / 'verify_manifest.py'),
                               '--expected', args.manifest_sha256], capture_output=True, text=True, timeout=30)
    require(verified.returncode == 0 and json.loads(verified.stdout) == authentication,
            'Independent manifest verifier disagrees with bootstrap')
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
            rows.append(execute(root, inventory, output, checker, mode))
        for corruption in CORRUPTIONS:
            rows.append(execute(root, inventory, output, CHECKERS[0], mode, corruption))
        sol = next(c for c in CHECKERS if c['name'] == 'sol')
        for corruption in SOL_CORRUPTIONS:
            rows.append(execute(root, inventory, output, sol, mode, corruption))
    agreements = []
    for checker in CHECKERS:
        pair = [r for r in rows if r['kind'] == 'positive' and r['checker'] == checker['name']]
        require(len(pair) == 2 and [r['mode'] for r in pair] == ['normal', 'optimized'], 'Mode pair differs')
        data = [json.loads((output / r['output']['path']).read_bytes()) for r in pair]
        ignored = ['python_optimization'] if checker['name'] == 'sol' else []
        for item in data:
            for field in ignored:
                del item[field]
        require(data[0] == data[1], 'Normal/optimized scientific outputs differ')
        agreements.append({'checker': checker['name'], 'scientific_json_equal': True,
                           'only_normalized_fields': ignored})
    positives = sum(r['kind'] == 'positive' for r in rows)
    external = sum(r['kind'] == 'external_corruption' for r in rows)
    embedded = sum(r.get('embedded_wrong_formula_controls', 0) for r in rows)
    groups = {m: sum(r.get('finite_positive_check_groups', 0) for r in rows if r['mode'] == m)
              for m in ['normal', 'optimized']}
    require(positives == 16 and external == 26 and embedded == 38
            and groups == {'normal': 135, 'optimized': 135}, 'Incomplete fresh replay inventory')
    require(authenticate(root, args.manifest_sha256) == (authentication, inventory), 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-uniform-waveform-fresh-replay-v1', 'status': 'PASS',
        'capsule': root.name, 'authentication': authentication, 'runner': inventory['reproduce.py'],
        'integrity_verifier': inventory['verify_manifest.py'], 'requirements': inventory['requirements.txt'],
        'python': {'version': sys.version, 'executable': sys.executable}, 'sympy_version': '1.14.0',
        'successful_checker_runs': positives, 'external_corruption_processes_rejected': external,
        'embedded_wrong_formula_controls_rejected': embedded,
        'finite_positive_check_groups_each_mode': 135, 'finite_positive_check_groups_by_mode': groups,
        'separate_generated_enclosures_each_mode': {'continuous_csv_cells': 4098,
            'gain_polynomial_identities': 3, 'peak_enclosures': 14},
        'mode_agreements': agreements,
        'scope': 'Fresh finite interval/arithmetic checks and conditional first-variation enclosures only. Generator JSON and CSV bytes reproduce in full. These executions do not prove the retained analytic PDE/exact-family hypotheses, compute an exact nonlinear waveform, set effective nonlinear amplitude or observational scales, or provide external peer review. All golden JSON fields compared; only the independently checked Sol optimization flag is normalized for comparison between modes.',
        'steps': rows}
    (output / 'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'successful_checker_runs',
        'external_corruption_processes_rejected', 'embedded_wrong_formula_controls_rejected',
        'finite_positive_check_groups_each_mode']}, indent=2))


if __name__ == '__main__':
    main()
