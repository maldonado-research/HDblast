"""Authenticated fresh prescribed quadratic-forcing calculations; the full U2 response remains unresolved."""
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

# BEGIN EMBEDDED APPROVED CONFIGURATION
CAPSULE = 'HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
PYTHON_VERSION = '3.12.14'
DEPENDENCY_VERSIONS = {}
CONFIG_SHA256 = None
CHECKERS = []
CONTROLS = []
# END EMBEDDED APPROVED CONFIGURATION

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


def read_generated(path, work):
    require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents])
            and path.resolve().is_relative_to(work.resolve()),
            'Fresh generated output must be a real regular file inside its workspace')
    return path.read_bytes()


def record(path, base):
    raw = read_generated(path, base)
    return {'path': path.relative_to(base).as_posix(), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}


def json_equal(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(json_equal(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(json_equal(a, b) for a, b in zip(left, right))
    return left == right


def at_path(data, path):
    for key in path:
        data = data[key]
    return data


def inventory_names(data, specification):
    if specification is None:
        return []
    value = at_path(data, specification['path'])
    if specification['kind'] == 'names':
        require(isinstance(value, list) and all(isinstance(x, str) and x for x in value), 'Named list required')
        names = value
    else:
        require(specification['kind'] == 'true_map' and isinstance(value, dict)
                and all(isinstance(k, str) and k and v is True for k, v in value.items()), 'Passing named boolean map required')
        names = list(value)
    require(len(names) == len(set(names)), 'Duplicate named group/control')
    return names


def positive_inventory(checker, data, mode):
    require(at_path(data, checker['status_path']) == checker['status'], 'Positive status differs: '+checker['name'])
    for field in checker['mode_fields']:
        value = at_path(data, field['path'])
        expected = field[mode]
        require(type(value) is type(expected) and value == expected, 'Execution-mode metadata differs')
    checks = inventory_names(data, checker['checks'])
    embedded = inventory_names(data, checker['embedded_checks'])
    require(checks == checker['check_names'] and embedded == checker['embedded_control_names'],
            'Exact approved group/control identities differ')
    for item in checker['required_fields']:
        value = at_path(data, item['path'])
        require(json_equal(value, item['equals']), 'Required positive output field differs')
    if checker.get('check_count_path') is not None:
        require(type(at_path(data, checker['check_count_path'])) is int
                and at_path(data, checker['check_count_path']) == len(checks), 'Reported group count differs')
    if checker.get('embedded_count_path') is not None:
        require(type(at_path(data, checker['embedded_count_path'])) is int
                and at_path(data, checker['embedded_count_path']) == len(embedded), 'Reported embedded count differs')
    return len(checks), len(embedded)


def validate_failure_diagnostic(control, diagnostic):
    if control['failure_format'] == 'fail_line':
        require(diagnostic == 'FAIL: '+control['reason']+'\n', 'Unrelated or ambiguous failure diagnostic')
    elif control['failure_format'] == 'json':
        failure = json.loads(diagnostic)
        require(at_path(failure, control.get('status_path', ['status'])) == control.get('status', 'FAIL')
                and at_path(failure, control['failure_path']) == control['reason'],
                'JSON failure does not contain the exact expected named reason')
    else:
        require(control['failure_format'] == 'traceback'
                and diagnostic.startswith(control.get('diagnostic_prefix', '')+'Traceback (most recent call last):\n')
                and diagnostic.endswith(control['exception']+': '+control['reason']+'\n')
                and diagnostic.splitlines()[-1] == control['exception']+': '+control['reason'],
                'Traceback did not end with the exact named expected failure')


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
        require(control['target'] != checker['script']
                and control['target'] in checker['dependencies'], 'Corruption must target a declared copied dependency, preserving the checker')
        original = target.read_bytes()
        if 'append' in control:
            mutated = original+control['append'].encode()
        elif 'json_path' in control:
            document = json.loads(original)
            parent = document
            for key in control['json_path'][:-1]:
                parent = parent[key]
            key = control['json_path'][-1]
            if 'factor' in control:
                value = parent[key]
                parent[key] = [str(F(v)*F(control['factor'])) for v in value] if isinstance(value, list) else str(F(value)*F(control['factor']))
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
    for path, _ in checker.get('additional_outputs', []):
        require(not (work/path).exists(), 'Stale additional generated output present')
    started = time.monotonic()
    with (work/'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(work/checker['script']), *arguments], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if control:
        diagnostic = read_generated(work/'execution.log', work).decode()
        require(result.returncode == 1, 'Corruption did not fail with exit one: '+control['name'])
        require(not checker['output'] or not (result_path.exists() or result_path.is_symlink()), 'A corrupted process emitted a passing output file')
        validate_failure_diagnostic(control, diagnostic)
        data = {'status': 'FAIL', 'failure': control['reason'], 'exit_code': result.returncode,
                'passing_output_absent': True, 'complete_diagnostic': diagnostic,
                'result_origin': 'Runner-parsed exact named failure; original diagnostic retained in execution.log'}
        raw = (json.dumps(data, indent=2)+'\n').encode()
    else:
        require(result.returncode == 0 and result_path.is_file(), 'Successful fresh result missing: '+checker['name'])
        raw = read_generated(result_path, work)
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
            'control_mechanism': control['mechanism'],
            'expected_reason': control['reason'], 'reported_failure': data['failure'],
            'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        golden_name = checker.get('optimized_golden', checker['golden']) if mode == 'optimized' else checker['golden']
        golden_raw = read_source(root, inventory, golden_name)
        require(json_equal(data, json.loads(golden_raw)), 'Complete scientific JSON differs from frozen golden: '+checker['name'])
        if checker.get('byte_equal'):
            require(raw == golden_raw, 'Complete generated JSON bytes differ: '+checker['name'])
        groups, embedded = positive_inventory(checker, data, mode)
        require(groups == checker['groups'] and embedded == checker['embedded'], 'Positive finite inventory differs')
        for path, golden in checker.get('additional_outputs', []):
            fresh_path = work/path
            require(read_generated(fresh_path, work) == read_source(root, inventory, golden),
                    'Complete generated data bytes differ: '+path)
            row['additional_outputs'].append({'output': record(fresh_path, output), 'golden': inventory[golden],
                                              'complete_bytes_equal': True})
        row.update({'golden': inventory[golden_name], 'scientific_output_matches': True,
                    'complete_generated_json_bytes_equal': checker.get('byte_equal', False),
                    'finite_positive_check_groups': groups, 'embedded_wrong_formula_controls': embedded,
                    'finite_positive_check_names': checker['check_names'],
                    'embedded_control_names': checker['embedded_control_names'],
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
                          'embedded_controls_each_mode': c['embedded'],
                          'finite_positive_check_names': c['check_names'],
                          'embedded_control_names': c['embedded_control_names']} for c in CHECKERS]
    expected_controls = [{'checker': c['checker'], 'name': c['name'], 'expected_reason': c['reason'],
                          'mechanism': c['mechanism']} for c in CONTROLS]
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
    require(CONFIG_SHA256 is not None and CHECKERS, 'An approved embedded runner configuration is required')
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
    require(sys.version.split()[0] == PYTHON_VERSION, 'Tested Python version differs')
    dependency_versions = {name: importlib.metadata.version(name) for name in DEPENDENCY_VERSIONS}
    require(dependency_versions == DEPENDENCY_VERSIONS, 'Pinned dependency version differs')
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
        normalized = [field['path'] for field in checker['mode_fields']]
        for item in data:
            for path in normalized:
                parent = item
                for key in path[:-1]:
                    parent = parent[key]
                del parent[path[-1]]
        require(json_equal(data[0], data[1]), 'Normal and optimized scientific outputs differ')
        agreements.append({'checker': checker['name'], 'scientific_json_equal': True, 'only_normalized_fields': normalized})
    actual = {'successful_checker_runs': sum(r['kind'] == 'positive' for r in rows),
        'external_corruption_processes_rejected': sum(r['kind'] == 'external_corruption' for r in rows),
        'embedded_wrong_formula_controls_rejected': sum(r.get('embedded_wrong_formula_controls', 0) for r in rows)}
    require(all(actual[k] == counts[k] for k in actual), 'Incomplete process inventory')
    groups = {m: sum(r.get('finite_positive_check_groups', 0) for r in rows if r['mode'] == m)
              for m in declared['modes']}
    require(all(v == counts['finite_positive_check_groups_each_mode'] for v in groups.values()), 'Per-mode group count differs')
    require(authenticate(root, args.manifest_sha256) == (authentication, inventory), 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-quadratic-shell-forcing-fresh-replay-v1', 'status': 'PASS', 'capsule': CAPSULE,
        'authentication': authentication, 'runner': inventory['reproduce.py'],
        'integrity_verifier': inventory['verify_manifest.py'], 'requirements': inventory['requirements.txt'],
        'roster': inventory['REPLAY_ROSTER.json'], 'python': {'version': sys.version, 'executable': sys.executable},
        'dependency_versions': dependency_versions, 'embedded_configuration_sha256': CONFIG_SHA256,
        **counts, 'finite_positive_check_groups_by_mode': groups,
        'mode_agreements': agreements, 'steps': rows,
        'scope': 'Fresh finite arithmetic for prescribed quadratic material shell forcing and known constitutive stress terms only. The full second-order U2 response remains unresolved quantitatively. Complete generated JSON and CSV reproduce byte for byte. All golden JSON fields are compared; only individually verified execution-mode metadata is normalized across modes. Executions do not prove analytic PDE/exact-family hypotheses, compute a finite-amplitude nonlinear waveform, set an effective nonlinear amplitude or observational scale, or replace external peer review.'}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', **counts}, indent=2))


if __name__ == '__main__':
    main()
