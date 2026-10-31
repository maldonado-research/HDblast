#!/usr/bin/env python3
"""Authenticate and replay the conditional exact-family TT transfer checkpoint.

The closed source authentication is a provenance check. The replay checks finite
algebra, inequalities, and rejection behavior; it does not reprove analytic PDE
or weighted constraint theorems.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def need(value, witness):
    if not value:
        raise RuntimeError(witness)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def regular(path):
    return (path.is_file() and not path.is_symlink()
            and all(not parent.is_symlink() for parent in path.absolute().parents))


def relative(value):
    need(isinstance(value, str) and bool(value), 'NONEMPTY_RELATIVE_PATH')
    path = PurePosixPath(value)
    need(not path.is_absolute() and '..' not in path.parts
         and path.as_posix() == value and value != '.', 'SAFE_RELATIVE_PATH:' + value)
    return Path(*path.parts)


def authenticate(root, expected=None):
    need(root.is_dir() and not root.is_symlink(), 'CAPSULE_ROOT:' + root.name)
    manifest = root / 'MANIFEST.json'
    need(regular(manifest), 'MANIFEST_REGULAR:' + root.name)
    digest = sha(manifest)
    if expected is not None:
        need(digest == expected, 'EXTERNAL_MANIFEST_PIN:' + root.name)
    document = read(manifest)
    need(isinstance(document['files'], list), 'MANIFEST_FILES:' + root.name)
    listed = {'MANIFEST.json'}
    for row in document['files']:
        rel = relative(row['path'])
        key = rel.as_posix()
        need(key not in listed, 'UNIQUE_MANIFEST_MEMBER:' + key)
        listed.add(key)
        path = root / rel
        need(regular(path) and type(row['bytes']) is int and row['bytes'] >= 0
             and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
             'CLOSED_MEMBER:' + root.name + '/' + key)
    actual = set()
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'CLOSED_SYMLINK:' + root.name)
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            need(path.is_dir(), 'CLOSED_SPECIAL_FILE:' + root.name)
    need(actual == listed, 'CLOSED_INVENTORY:' + root.name)
    return digest


def sources(repo, plan):
    result = {name: authenticate(repo / 'research' / relative(name), pin)
              for name, pin in sorted(plan['source_capsules'].items())}
    for name, row in sorted(plan['selected_source_files'].items()):
        path = repo / relative(name)
        need(regular(path) and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
             'SELECTED_SOURCE:' + name)
    return result


def resolve(value, bindings):
    return value.format_map(bindings)


def dig(document, path):
    for key in path.split('.'):
        document = document[key]
    return document


def controls(document, contract, mode=None):
    for field, expected in contract.get('fields', {}).items():
        need(dig(document, field) == expected, 'CONTROL_REQUIRED_FIELD:' + field)
    rows = dig(document, contract['rows'])
    if 'filter_field' in contract:
        rows = [row for row in rows if row.get(contract['filter_field']) in contract['filter_values']]
    need(len(rows) == contract['count'], 'COMPLETE_CONTROL_COUNT')
    pairs = {}
    for row in rows:
        code = row[contract['returncode']]
        need(type(code) is int and code != 0, 'ACTUAL_NONZERO_CONTROL_EXIT')
        if 'expected_returncode' in contract:
            need(code == contract['expected_returncode'], 'EXPECTED_CONTROL_EXIT')
        for field, expected in contract.get('row_fields', {}).items():
            need(dig(row, field) == expected, 'CONTROL_ROW_REQUIRED_FIELD:' + field)
        for field in contract.get('true', []):
            need(dig(row, field) is True, 'CONTROL_TRUE:' + field)
        for field in contract.get('false', []):
            need(dig(row, field) is False, 'CONTROL_FALSE:' + field)
        for field in contract.get('nonempty', []):
            need(isinstance(dig(row, field), str) and bool(dig(row, field).strip()),
                 'ACTUAL_CONTROL_WITNESS:' + field)
        if 'mode' in contract:
            pairs.setdefault(row[contract['case']], []).append(row[contract['mode']])
    if 'mode' in contract:
        expected_modes = contract['mode_values']
        need(all(len(values) == 2 and set(values) == set(expected_modes)
                 for values in pairs.values()), 'CONTROL_BOTH_PYTHON_MODES')
    categories = {}
    for row in rows:
        label = contract.get('category_by_case', {}).get(row.get(contract.get('case')),
                    row.get(contract.get('category_field'), contract.get('category', 'unclassified')))
        categories[label] = categories.get(label, 0) + 1
    return {'rejected_subprocesses': len(rows), 'categories': categories}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--jobs', type=int, choices=(1, 2, 3, 4), default=3)
    args = parser.parse_args()
    repo = args.repo_root.resolve(strict=True)
    cap = HERE.parent
    out = args.output_dir.absolute()
    need(repo.is_dir() and not args.repo_root.is_symlink(), 'REPOSITORY_ROOT')
    need(not args.output_dir.exists() and not args.output_dir.is_symlink(), 'FRESH_OUTPUT_DIRECTORY')
    need(all(not parent.is_symlink() for parent in out.parents), 'OUTPUT_ANCESTOR_SYMLINK')
    out = out.resolve()
    need(repo != out and repo not in out.parents and cap != out and cap not in out.parents,
         'OUTPUT_OUTSIDE_FROZEN_SOURCE')
    own = authenticate(cap, args.manifest_sha256)
    plan = read(HERE / 'REPLAY_PLAN.json')
    need(plan['schema'] == 'HDBLAST_EXACT_FAMILY_TT_ROOT_REPLAY_V1', 'PLAN_SCHEMA')
    before = sources(repo, plan)
    # No project scientific code or third-party module has been imported above.
    env = dict(os.environ)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    command = [sys.executable, '-I', '-B', '-c',
               'import json,sympy,mpmath; print(json.dumps({"sympy":sympy.__version__,"mpmath":mpmath.__version__},sort_keys=True))']
    process = subprocess.run(command, capture_output=True, text=True, timeout=30, env=env)
    need(process.returncode == 0, 'PINNED_SYMBOLIC_REQUIREMENTS_MISSING')
    versions = json.loads(process.stdout)
    need(versions == plan['requirements'], 'PINNED_SYMBOLIC_REQUIREMENTS_VERSION')
    out.mkdir(parents=True)
    (out / 'logs').mkdir()
    bindings = {'CAP': str(cap), 'SOURCE': str(repo), 'OUT': str(out)}
    timings = {}
    started = time.monotonic()

    def execute(group):
        canonical, modes, rejected, extras, local_times = {}, {}, {}, {}, {}
        for step in group['steps']:
            argv = [resolve(value, bindings) for value in step['argv']]
            flags = ['-I', '-B'] if step.get('requires_sympy') else ['-I', '-S', '-B']
            command = [sys.executable, *flags, *(['-O'] if step.get('optimized') else []), *argv]
            for key in ('output', 'receipt'):
                if key in step and not (key == 'receipt' and step['kind'] == 'retained_checkpoint'):
                    Path(resolve(step[key], bindings)).parent.mkdir(parents=True, exist_ok=True)
            clock = time.monotonic()
            result = subprocess.run(command, capture_output=True, timeout=step.get('timeout_seconds', 300),
                                    env=env, cwd=repo)
            local_times[step['name']] = round(time.monotonic() - clock, 6)
            log = out / 'logs' / (group['name'] + '-' + step['name'])
            log.with_suffix('.stdout').write_bytes(result.stdout)
            log.with_suffix('.stderr').write_bytes(result.stderr)
            need(result.returncode == 0, 'REPLAY_STEP:' + group['name'] + '/' + step['name'] + ':'
                 + result.stderr.decode(errors='replace')[-3000:])
            if step['kind'] == 'canonical':
                produced = Path(resolve(step['output'], bindings))
                if step.get('stdout_is_report'):
                    need(not produced.exists(), 'FRESH_STDOUT_REPORT')
                    produced.write_bytes(result.stdout)
                expected = cap / relative(step['reference'])
                need(sha(expected) == step['expected_sha256'], 'REFERENCE_PIN:' + step['name'])
                need(regular(produced) and produced.read_bytes() == expected.read_bytes(),
                     'CANONICAL_BYTE_IDENTITY:' + step['name'])
                key = step['canonical_key']
                digest = sha(produced)
                need(key not in canonical or canonical[key] == digest, 'NORMAL_OPTIMIZED_IDENTITY')
                canonical[key] = digest
                modes.setdefault(key, []).append(bool(step.get('optimized')))
                document = read(produced)
                for field, expected_value in step.get('require_fields', {}).items():
                    need(dig(document, field) == expected_value, 'CANONICAL_REQUIRED_FIELD:' + field)
                if 'controls' in step:
                    rejected[step['name']] = controls(document, step['controls'], bool(step.get('optimized')))
                for field, count in step.get('inprocess_formula_controls', {}).items():
                    value = dig(document, field)
                    need((len(value) if isinstance(value, (list, dict)) else value) == count,
                         'INPROCESS_FORMULA_CONTROLS:' + field)
                    extras[step['name'] + ':' + field] = count
            elif step['kind'] == 'controls':
                receipt = Path(resolve(step['receipt'], bindings))
                rejected[step['name']] = controls(read(receipt), step['contract'])
                if 'reference' in step:
                    expected = cap / relative(step['reference'])
                    need(sha(expected) == step['expected_sha256'] and receipt.read_bytes() == expected.read_bytes(),
                         'CONTROL_RECEIPT_BYTE_IDENTITY:' + step['name'])
            elif step['kind'] == 'retained_checkpoint':
                document = read(Path(resolve(step['receipt'], bindings)))
                need(document['status'] == 'PASS_COMPLETE_POSTCONTACT_SCREEN_CHECKPOINT', 'RETAINED_STATUS')
                need(document['capsule_manifest_sha256'] == step['manifest_sha256'], 'RETAINED_MANIFEST_PIN')
                need(document['canonical_reports_reproduced_normal_and_optimized'] == 9
                     and document['real_negative_subprocesses'] == 68, 'RETAINED_COMPLETE_ROSTER')
                need(document['all_frozen_sources_authenticated_before_scientific_imports_and_after'] is True
                     and document['scientific_sources_modified'] is False, 'RETAINED_SOURCE_IMMUTABILITY')
                extras['retained_checkpoint'] = {
                    'canonical_report_pairs': 9, 'real_rejected_subprocesses': 68,
                    'inprocess_report_data_controls_distinct': document['report_only_controls_distinct'],
                    'inprocess_report_data_mode_evaluations': document['report_only_controls_mode_evaluations'],
                    'focused_code_checks_per_mode': document['focused_code_checks_per_mode'],
                    'capsule_manifest_sha256': document['capsule_manifest_sha256']}
            else:
                raise RuntimeError('UNKNOWN_REPLAY_STEP_KIND')
        need(all(len(values) == 2 and set(values) == {False, True} for values in modes.values()),
             'EACH_NEW_CANONICAL_REPORT_BOTH_MODES')
        summary = {'group': group['name'], 'canonical_sha256': canonical,
                   'real_rejected_subprocesses': rejected, 'distinct_finite_check_evaluations': extras}
        print(json.dumps({'completed_group': group['name'], 'new_report_pairs': len(canonical),
                          'new_rejected_subprocesses': sum(item['rejected_subprocesses'] for item in rejected.values())}),
              flush=True)
        return summary, local_times

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        executed = list(pool.map(execute, plan['groups']))
    results = [item[0] for item in executed]
    timings.update({group['group']: timing for group, timing in executed})
    need(authenticate(cap, own) == own and sources(repo, plan) == before, 'ALL_FROZEN_SOURCES_UNCHANGED')
    total = sum(sum(item['rejected_subprocesses'] for item in group['real_rejected_subprocesses'].values())
                for group in results)
    reports = sum(len(group['canonical_sha256']) for group in results)
    retained = [group['distinct_finite_check_evaluations']['retained_checkpoint'] for group in results
                if 'retained_checkpoint' in group['distinct_finite_check_evaluations']]
    need(len(retained) == 1, 'RETAINED_CHECKPOINT_EXECUTED_ONCE')
    need(total == plan['expected_new_rejected_subprocesses']
         and reports == plan['expected_new_canonical_reports'], 'COMPLETE_NEW_REPLAY_ROSTER')
    summary = {'schema': 'HDBLAST_EXACT_FAMILY_TT_ROOT_VERIFICATION_V1',
               'status': 'PASS_COMPLETE_EXACT_FAMILY_TT_TRANSFER_CHECKPOINT',
               'capsule_manifest_sha256': own, 'source_capsule_manifest_sha256': before,
               'all_frozen_sources_authenticated_before_scientific_imports_and_after': True,
               'scientific_sources_modified': False, 'symbolic_versions': versions,
               'new_canonical_report_pairs': reports, 'new_real_rejected_subprocesses': total,
               'retained_checkpoint': retained[0], 'groups': results, 'scope': plan['scope']}
    (out / 'VERIFICATION.json').write_text(json.dumps(summary, sort_keys=True, indent=2) + '\n')
    # Timing/logs are runtime evidence outside the authenticated publication.
    elapsed = round(time.monotonic() - started, 6)
    (out / 'RUNTIME.json').write_text(json.dumps({'step_seconds': timings, 'total_seconds': elapsed}, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': summary['status'], 'new_report_pairs': reports,
                      'new_rejected_subprocesses': total, 'retained_report_pairs': 9,
                      'retained_rejected_subprocesses': 68, 'runtime_seconds': elapsed}), flush=True)


if __name__ == '__main__':
    main()
