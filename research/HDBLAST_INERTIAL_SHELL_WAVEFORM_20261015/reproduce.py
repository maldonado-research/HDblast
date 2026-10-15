"""Authenticated fresh inertial-waveform calculations, separate from analytic PDE proof."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import time

CAPSULE = 'HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015'
SHELL = 'inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json'
BACKGROUND = 'archive/BACKGROUND_DERIVATIVE_BOUNDS.json'
GAIN = 'inputs/previous/archive/PHYSICAL_GAIN_ENCLOSURE.json'
CHECKERS = [
    {'name': 'root', 'script': 'check_inertial_waveform.py',
     'dependencies': ['inertial_error_bounds.py', 'compute_inertial_waveform.py', 'sol/inertial_comparator.py',
        BACKGROUND, GAIN, SHELL, 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json', 'computed/INERTIAL_WAVEFORM_CELLS.csv'],
     'arguments': ['--data-dir', 'computed', '--output', 'checked.json'],
     'golden': 'computed/CHECK_INERTIAL_WAVEFORM.json', 'output': 'checked.json', 'groups': 31, 'embedded': 0},
    {'name': 'cells', 'script': 'referee/check_independent_cells.py',
     'dependencies': ['inertial_error_bounds.py', 'compute_inertial_waveform.py', 'sol/inertial_comparator.py',
        BACKGROUND, GAIN, SHELL, 'referee/INDEPENDENT_INERTIAL_BOUND_NORMAL.json',
        'computed/INERTIAL_WAVEFORM_ENCLOSURE.json', 'computed/INERTIAL_WAVEFORM_CELLS.csv'],
     'arguments': ['--data-dir', 'computed', '--output', 'checked.json'],
     'golden': 'referee/INDEPENDENT_CELLS_NORMAL.json',
     'optimized_golden': 'referee/INDEPENDENT_CELLS_OPTIMIZED.json', 'optimization_field': 'python_optimization',
     'output': 'checked.json', 'groups': 18, 'embedded': 0},
    {'name': 'generator', 'script': 'compute_inertial_waveform.py',
     'dependencies': ['inertial_error_bounds.py', 'sol/inertial_comparator.py', BACKGROUND, GAIN, SHELL],
     'arguments': ['--output-dir', 'generated'],
     'golden': 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json', 'output': 'generated/INERTIAL_WAVEFORM_ENCLOSURE.json',
     'additional_outputs': [('generated/INERTIAL_WAVEFORM_CELLS.csv', 'computed/INERTIAL_WAVEFORM_CELLS.csv')],
     'groups': 0, 'embedded': 0, 'byte_equal': True},
    {'name': 'sol', 'script': 'sol/inertial_comparator.py', 'dependencies': [SHELL],
     'arguments': ['--source-root', '.'], 'golden': 'sol/COMPARATOR_NORMAL.json',
     'optimized_golden': 'sol/COMPARATOR_OPTIMIZED.json', 'optimization_field': 'python_optimization',
     'output': None, 'groups': 13, 'embedded': 0},
    {'name': 'astra', 'script': 'astra/check_inertial_error.py',
     'dependencies': [SHELL, 'inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md'],
     'arguments': ['--source-root', '.'], 'golden': 'astra/INERTIAL_ERROR_NORMAL.json',
     'optimized_golden': 'astra/INERTIAL_ERROR_OPTIMIZED.json', 'optimization_field': 'optimized_mode',
     'output': None, 'groups': 23, 'embedded': 5},
    {'name': 'referee', 'script': 'referee/check_independent_inertial_bound.py',
     'dependencies': [BACKGROUND, 'archive/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'],
     'arguments': ['--source-root', '.'], 'golden': 'referee/INDEPENDENT_INERTIAL_BOUND_NORMAL.json',
     'optimized_golden': 'referee/INDEPENDENT_INERTIAL_BOUND_OPTIMIZED.json', 'optimization_field': 'python_optimization',
     'output': None, 'groups': 21, 'embedded': 0},
    {'name': 'archive', 'script': 'archive/bound_background_derivatives.py',
     'dependencies': ['archive/inputs/'+n for n in ['OCT11_JACOBI_ENDPOINT_MONOTONICITY.md',
         'OCT11_RESCALED_BOUND_CERTIFICATE.json', 'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json']],
     'golden': BACKGROUND, 'output': BACKGROUND, 'groups': 0, 'embedded': 0, 'byte_equal': True},
    {'name': 'identities', 'script': 'archive/check_background_derivative_identities.py', 'dependencies': [],
     'golden': 'archive/BACKGROUND_DERIVATIVE_IDENTITY_CHECK.json',
     'output': 'archive/BACKGROUND_DERIVATIVE_IDENTITY_CHECK.json', 'groups': 6, 'embedded': 4, 'byte_equal': True},
]
CONTROLS = [
    {'checker': 'root', 'name': 'wrong_inverse_width_coefficient', 'reason': 'correct inverse-width norm Jacobian',
     'failure_format': 'fail_line', 'target': 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json',
     'json_path': ['quadrature', 'C_interval'], 'transform': 'multiply_by_ell_squared'},
    {'checker': 'root', 'name': 'omit_incoming_velocity', 'reason': 'specified nonzero incoming velocity retained',
     'failure_format': 'fail_line', 'target': 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json',
     'json_path': ['state', 'initial_velocity'], 'value': '0'},
    {'checker': 'root', 'name': 'claim_zero_uniqueness', 'reason': 'only existence claimed for two zero intervals',
     'failure_format': 'fail_line', 'target': 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json',
     'json_path': ['zero_uniqueness_claimed'], 'value': True},
    {'checker': 'root', 'name': 'shorten_signed_tail_domain', 'reason': 'strict nonzero signed tail theorem',
     'failure_format': 'fail_line', 'target': 'computed/INERTIAL_WAVEFORM_ENCLOSURE.json',
     'json_path': ['tail_theorem', 'domain'], 'value': [2, 2.9]},
    {'checker': 'root', 'name': 'alter_complete_csv', 'reason': 'complete CSV authenticated',
     'failure_format': 'fail_line', 'target': 'computed/INERTIAL_WAVEFORM_CELLS.csv', 'append': '\n'},
    {'checker': 'sol', 'name': 'reverse_kernel_initial', 'reason': 'kernel derivative strictly positive on complete interval',
     'arguments': ['--source-root', '.', '--mass-cells', '64', '--fault', 'reverse_kernel_initial'], 'failure_format': 'fail_line'},
    {'checker': 'sol', 'name': 'reverse_kernel_sign', 'reason': 'positive postpulse integrand on full continuous rectangle',
     'arguments': ['--source-root', '.', '--mass-cells', '64', '--fault', 'reverse_kernel_sign'], 'failure_format': 'fail_line'},
    {'checker': 'sol', 'name': 'omit_local_inertia', 'reason': 'center retains negative local inertia correction',
     'arguments': ['--source-root', '.', '--mass-cells', '64', '--fault', 'omit_local_inertia'], 'failure_format': 'fail_line'},
    {'checker': 'sol', 'name': 'reverse_q_derivative', 'reason': 'positive incoming lobe',
     'arguments': ['--source-root', '.', '--mass-cells', '64', '--fault', 'reverse_q_derivative'], 'failure_format': 'fail_line'},
    {'checker': 'sol', 'name': 'alter_authenticated_shell', 'reason': 'shell input authentication failed',
     'arguments': ['--source-root', '.', '--mass-cells', '64'], 'failure_format': 'fail_line',
     'target': SHELL, 'append': '\n'},
    {'checker': 'archive', 'name': 'alter_retained_coefficient_certificate',
     'reason': 'Input pin differs: OCT11_RESCALED_BOUND_CERTIFICATE.json', 'failure_format': 'traceback',
     'exception': 'ValueError', 'target': 'archive/inputs/OCT11_RESCALED_BOUND_CERTIFICATE.json', 'append': '\n'},
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
    if 'optimization_field' in checker:
        field = checker['optimization_field']
        expected = bool(mode == 'optimized') if field == 'optimized_mode' else int(mode == 'optimized')
        require(type(data[field]) is type(expected) and data[field] == expected, 'Optimization metadata differs: '+name)
    statuses = {'root': 'PASS_INERTIAL_SHELL_WAVEFORM',
        'cells': 'PASS_INDEPENDENT_CONTINUOUS_INERTIAL_CELLS_AND_NORM',
        'sol': 'PASS_EXACT_FREE_INERTIAL_COMPARATOR',
        'astra': 'PASS_CONDITIONAL_INERTIAL_COMPARATOR_AND_POSITIVE_TAIL',
        'referee': 'PASS_INDEPENDENT_CONDITIONAL_INERTIAL_RESIDUAL_BOUND',
        'generator': 'PASS', 'archive': 'PASS_EXACT_RATIONAL_WIDER_COLLAR_AND_DERIVATIVE_BOUNDS',
        'identities': 'PASS_FORMAL_COEFFICIENT_IDENTITIES_ONLY'}
    require(data['status'] == statuses[name], 'Positive status differs: '+name)
    if name in ('root', 'cells', 'sol', 'astra', 'referee'):
        checks = data['checks']
        require(len(checks) == len(set(checks)) == checker['groups'], 'Distinct named check inventory differs')
        require(name == 'astra' or data['check_count'] == len(checks), 'Declared check count differs')
        embedded = data.get('negative_controls_rejected', [])
        require(len(embedded) == len(set(embedded)) == checker['embedded'], 'Embedded control inventory differs')
        if name == 'sol':
            require(data['mass_mesh_cells'] == 4096, 'Full comparator mass grid required')
        if name in ('root', 'cells'):
            require(data['csv_rows' if name == 'root' else 'rows_reconstructed'] == 4353, 'Complete continuous grid required')
        return len(checks), len(embedded)
    if name == 'identities':
        require(len(data['checks']) == 6 and all(v is True for v in data['checks'].values())
                and len(data['wrong_formula_controls']) == 4 and all(v is True for v in data['wrong_formula_controls'].values())
                and data['sympy_version'] == '1.14.0', 'Formal identity inventory differs')
        return 6, 4
    if name == 'generator':
        require(data['cell_data']['support_cells'] == 4096 and data['cell_data']['tail_cells'] == 256
                and data['cell_data']['rows'] == 4353 and data['zero_uniqueness_claimed'] is False,
                'Generator grid or existence-only scope differs')
    if name == 'archive':
        require(data['first_exit']['status'] == 'PASS' and len(data['first_exit']['component_margins']) == 4
                and data['shell_correlation']['scalar_panels'] == 128, 'Derivative producer inventory differs')
    return 0, 0


def execute(root, inventory, output, checker, mode, control=None):
    work = output / mode / checker['name'] / (control['name'] if control else 'positive')
    work.mkdir(parents=True)
    copied = []
    for name in [checker['script'], *checker['dependencies']]:
        raw = read_source(root, inventory, name)
        target = work / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        require(target.read_bytes() == raw, 'Fresh copied bytes differ: '+name)
        copied.append({'source': inventory[name], 'copy': record(target, output)})
    mutation = None
    if control and 'target' in control:
        target = work / control['target']
        original = target.read_bytes()
        if 'append' in control:
            mutated = original+control['append'].encode()
        elif 'json_path' in control:
            document = json.loads(original)
            parent = document
            for key in control['json_path'][:-1]:
                parent = parent[key]
            key = control['json_path'][-1]
            if control.get('transform') == 'multiply_by_ell_squared':
                parent[key] = [str(F(v)*F(document['state']['ell'])**2) for v in parent[key]]
            else:
                parent[key] = control['value']
            mutated = (json.dumps(document, indent=2)+'\n').encode()
        else:
            require(original.decode().count(control['before']) == 1, 'Formula replacement must be unique')
            mutated = original.decode().replace(control['before'], control['after'], 1).encode()
        require(mutated != original, 'Corruption must alter copied bytes')
        target.write_bytes(mutated)
        mutation = {'source': inventory[control['target']], 'mutated': record(target, output),
                    'recipe': control, 'independent_checker_unchanged': True}
    arguments = control.get('arguments', checker.get('arguments', [])) if control else checker.get('arguments', [])
    flags = ['-I', '-B']+(['-O'] if mode == 'optimized' else [])
    result_path = work / checker['output'] if checker['output'] else work/'execution.log'
    if checker['output']:
        require(not result_path.exists(), 'Stale output present before execution')
    started = time.monotonic()
    with (work/'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(work/checker['script']), *arguments], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if control:
        diagnostic = (work/'execution.log').read_text()
        require(result.returncode == 1, 'Corruption did not fail with exit one: '+control['name'])
        require(not checker['output'] or not result_path.exists(), 'A corrupted process emitted a passing output file')
        if control['failure_format'] == 'fail_line':
            require(diagnostic == 'FAIL: '+control['reason']+'\n', 'Unrelated or ambiguous failure diagnostic')
        else:
            require(control['failure_format'] == 'traceback'
                    and diagnostic.startswith('Traceback (most recent call last):\n')
                    and diagnostic.rstrip().splitlines()[-1] == control['exception']+': '+control['reason'],
                    'Traceback did not end with the exact named expected failure')
        data = {'status': 'FAIL', 'failure': control['reason'], 'exit_code': result.returncode,
                'passing_output_absent': True, 'complete_diagnostic': diagnostic,
                'result_origin': 'Runner-parsed exact named failure; original diagnostic retained in execution.log'}
        raw = (json.dumps(data, indent=2)+'\n').encode()
    else:
        require(result.returncode == 0 and result_path.is_file(), 'Successful fresh result missing: '+checker['name'])
        raw = result_path.read_bytes()
        data = json.loads(raw)
    fresh = work/'fresh.json'
    fresh.write_bytes(raw)
    row = {'kind': 'external_corruption' if control else 'positive', 'checker': checker['name'], 'mode': mode,
           'script': checker['script'], 'flags': flags, 'arguments': arguments, 'exit_code': result.returncode,
           'wall_seconds': round(time.monotonic()-started, 3), 'copied_inputs': copied,
           'copied_inputs_snapshot': 'Authenticated initial copies before any separately recorded corruption',
           'log': record(work/'execution.log', output), 'output': record(fresh, output), 'additional_outputs': [],
           'output_origin': 'Exact named failure diagnostic' if control else checker['output'] or 'complete stdout JSON'}
    require((work/checker['script']).read_bytes() == read_source(root, inventory, checker['script']),
            'Independent executed checker script changed')
    if control:
        row.update({'mutation': mutation, 'control': control['name'],
            'control_mechanism': 'copied input/certificate bytes' if mutation else 'production formula CLI injection',
            'expected_reason': control['reason'], 'reported_failure': data['failure'],
            'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        golden_name = checker.get('optimized_golden', checker['golden']) if mode == 'optimized' else checker['golden']
        golden_raw = read_source(root, inventory, golden_name)
        require(data == json.loads(golden_raw), 'Complete scientific JSON differs from frozen golden: '+checker['name'])
        if checker.get('byte_equal'):
            require(raw == golden_raw, 'Complete generated JSON bytes differ: '+checker['name'])
        groups, embedded = positive_inventory(checker, data, mode)
        require(groups == checker['groups'] and embedded == checker['embedded'], 'Positive finite inventory differs')
        for path, golden in checker.get('additional_outputs', []):
            fresh_path = work/path
            require(fresh_path.is_file() and fresh_path.read_bytes() == read_source(root, inventory, golden),
                    'Complete generated data bytes differ: '+path)
            row['additional_outputs'].append({'output': record(fresh_path, output), 'golden': inventory[golden],
                                              'complete_bytes_equal': True})
        row.update({'golden': inventory[golden_name], 'scientific_output_matches': True,
                    'complete_generated_json_bytes_equal': checker.get('byte_equal', False),
                    'finite_positive_check_groups': groups, 'embedded_wrong_formula_controls': embedded,
                    'ignored_golden_fields': []})
    return row


def roster_counts():
    return {'successful_checker_runs': 2*len(CHECKERS),
            'external_corruption_processes_rejected': 2*len(CONTROLS),
            'embedded_wrong_formula_controls_rejected': 2*sum(c['embedded'] for c in CHECKERS),
            'finite_positive_check_groups_each_mode': sum(c['groups'] for c in CHECKERS)}


def validate_declared_roster(declared):
    counts = roster_counts()
    expected_positive = [{'name': c['name'], 'script': c['script'],
                          'finite_positive_groups_each_mode': c['groups'],
                          'embedded_controls_each_mode': c['embedded']} for c in CHECKERS]
    expected_controls = [{'checker': c['checker'], 'name': c['name'], 'expected_reason': c['reason'],
                          'mechanism': 'copied input/certificate corruption' if 'target' in c
                                       else 'production formula CLI injection'} for c in CONTROLS]
    require(declared['schema'] == 'hdblast-replay-roster-v1' and declared['status'] == 'APPROVED_FOR_FROZEN_REPLAY'
            and declared['capsule'] == CAPSULE and declared['modes'] == ['normal', 'optimized']
            and declared['counts'] == counts and declared['positive_processes'] == expected_positive
            and declared['external_controls'] == expected_controls,
            'Frozen reviewed roster differs from actual runner')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).absolute().parent
    authentication, inventory = authenticate(root, args.manifest_sha256)
    verified = subprocess.run([sys.executable, '-I', '-S', '-B', str(root/'verify_manifest.py'),
                               '--expected', args.manifest_sha256], capture_output=True, text=True, timeout=30)
    require(verified.returncode == 0 and json.loads(verified.stdout) == authentication,
            'Independent source verifier disagrees with bootstrap')
    declared = json.loads(read_source(root, inventory, 'REPLAY_ROSTER.json'))
    counts = validate_declared_roster(declared)
    requested = args.output_dir.absolute()
    require(not any(p.is_symlink() for p in [requested, *requested.parents]), 'Real output path required')
    output = requested.resolve()
    require(not output.is_relative_to(root.resolve()), 'Output must be outside frozen source')
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())), 'New or empty output required')
    require(importlib.metadata.version('sympy') == '1.14.0', 'Pinned SymPy 1.14.0 required')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    checkers_by_name = {c['name']: c for c in CHECKERS}
    require(len(checkers_by_name) == len(CHECKERS), 'Duplicate checker identity')
    for mode in declared['modes']:
        for checker in CHECKERS:
            rows.append(execute(root, inventory, output, checker, mode))
        for control in CONTROLS:
            rows.append(execute(root, inventory, output, checkers_by_name[control['checker']], mode, control))
    agreements = []
    for checker in CHECKERS:
        pair = [r for r in rows if r['kind'] == 'positive' and r['checker'] == checker['name']]
        require(len(pair) == 2 and [r['mode'] for r in pair] == declared['modes'], 'Complete mode pair required')
        data = [json.loads((output/r['output']['path']).read_bytes()) for r in pair]
        normalized = [checker['optimization_field']] if 'optimization_field' in checker else []
        for item in data:
            for key in normalized:
                del item[key]
        require(data[0] == data[1], 'Normal and optimized scientific outputs differ')
        agreements.append({'checker': checker['name'], 'scientific_json_equal': True, 'only_normalized_fields': normalized})
    actual = {'successful_checker_runs': sum(r['kind'] == 'positive' for r in rows),
        'external_corruption_processes_rejected': sum(r['kind'] == 'external_corruption' for r in rows),
        'embedded_wrong_formula_controls_rejected': sum(r.get('embedded_wrong_formula_controls', 0) for r in rows)}
    require(all(actual[k] == counts[k] for k in actual), 'Incomplete process inventory')
    groups = {m: sum(r.get('finite_positive_check_groups', 0) for r in rows if r['mode'] == m)
              for m in declared['modes']}
    require(all(v == counts['finite_positive_check_groups_each_mode'] for v in groups.values()), 'Per-mode group count differs')
    require(authenticate(root, args.manifest_sha256) == (authentication, inventory), 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-inertial-waveform-fresh-replay-v1', 'status': 'PASS', 'capsule': CAPSULE,
        'authentication': authentication, 'runner': inventory['reproduce.py'],
        'integrity_verifier': inventory['verify_manifest.py'], 'requirements': inventory['requirements.txt'],
        'roster': inventory['REPLAY_ROSTER.json'], 'python': {'version': sys.version, 'executable': sys.executable},
        'sympy_version': '1.14.0', **counts, 'finite_positive_check_groups_by_mode': groups,
        'separate_generated_enclosures_each_mode': {'continuous_csv_cells': 4353, 'support_cells': 4096,
            'tail_cells': 256, 'archive_collar_components': 4, 'archive_correlated_shell_panels': 128},
        'mode_agreements': agreements, 'steps': rows,
        'scope': 'Fresh finite arithmetic and conditional inertial first-variation enclosures only. Complete generated JSON and CSV reproduce byte for byte. All golden JSON fields are compared; only individually verified execution-mode metadata is normalized across modes. Executions do not prove analytic PDE/exact-family hypotheses, compute a finite-amplitude nonlinear waveform, set an effective nonlinear amplitude or observational scale, or replace external peer review.'}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', **counts}, indent=2))


if __name__ == '__main__':
    main()
