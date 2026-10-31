#!/usr/bin/env python3
"""Authenticate and replay the complete proper-clock screen checkpoint.

Closed source authentication is provenance, not a replay of inherited PDE
theorems or prior checkpoints. Only the new proper-clock computations run.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def need(ok, witness):
    if not ok:
        raise ValueError(witness)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def regular(path):
    return path.is_file() and not path.is_symlink() and all(not x.is_symlink() for x in path.absolute().parents)


def relative(value):
    need(isinstance(value, str) and bool(value), 'NONEMPTY_RELATIVE_PATH')
    p = PurePosixPath(value)
    need(not p.is_absolute() and '..' not in p.parts and p.as_posix() == value and value != '.', 'SAFE_RELATIVE_PATH')
    return Path(*p.parts)


def authenticate(root, pin=None):
    need(root.is_dir() and not root.is_symlink(), 'CAPSULE_ROOT:'+root.name)
    manifest = root/'MANIFEST.json'
    need(regular(manifest), 'MANIFEST_REGULAR:'+root.name)
    actual_pin = sha(manifest)
    if pin is not None:
        need(actual_pin == pin, 'EXTERNAL_MANIFEST_PIN:'+root.name)
    rows = read(manifest)['files']
    if isinstance(rows, dict):
        rows = [{'path': rel, **row} for rel, row in rows.items()]
    need(isinstance(rows, list), 'MANIFEST_ROW_LIST')
    seen = {'MANIFEST.json'}
    for row in rows:
        rel = relative(row['path'])
        need(str(rel) not in seen, 'DUPLICATE_MANIFEST_MEMBER')
        seen.add(str(rel))
        p = root/rel
        need(regular(p) and type(row['bytes']) is int and row['bytes'] >= 0 and
             p.stat().st_size == row['bytes'] and sha(p) == row['sha256'],
             'CLOSED_MEMBER:'+root.name+'/'+str(rel))
    actual = set()
    for p in root.rglob('*'):
        need(not p.is_symlink(), 'CLOSED_SYMLINK:'+root.name)
        if p.is_file():
            actual.add(str(p.relative_to(root)))
        else:
            need(p.is_dir(), 'CLOSED_SPECIAL_FILE:'+root.name)
    need(actual == seen, 'CLOSED_INVENTORY:'+root.name)
    return {'manifest_sha256': actual_pin, 'members_including_manifest': len(seen)}


def sources(repo, cap, plan):
    result = {name: authenticate(repo/'research'/relative(name), pin)
              for name, pin in sorted(plan['source_capsules'].items())}
    for rel, row in plan['selected_source_files'].items():
        p = repo/relative(rel)
        need(regular(p) and p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], 'SELECTED_OLD_SOURCE:'+rel)
    for rel, row in plan['critical_new_files'].items():
        p = cap/relative(rel)
        need(regular(p) and p.stat().st_size == row['bytes'] and sha(p) == row['sha256'], 'SELECTED_NEW_SOURCE:'+rel)
    return result


def dig(value, path):
    for key in path.split('.'):
        value = value[key]
    return value


def rejection_contract(document, contract):
    for key, expected in contract.get('fields', {}).items():
        need(dig(document, key) == expected, 'CONTROL_REQUIRED_FIELD:'+key)
    rows = dig(document, contract['rows'])
    need(len(rows) == contract['count'], 'COMPLETE_CONTROL_COUNT')
    categories = {}
    cases = {}
    for row in rows:
        code = row[contract['returncode']]
        need(type(code) is int and code != 0, 'ACTUAL_NONZERO_CONTROL_EXIT')
        for key in contract.get('true', []):
            need(dig(row, key) is True, 'CONTROL_TRUE:'+key)
        for key in contract.get('nonempty', []):
            need(isinstance(dig(row, key), str) and bool(dig(row, key).strip()), 'CONTROL_WITNESS:'+key)
        if 'observed_witness' in row:
            need(row['observed_witness'] == row['expected_witness'], 'EXACT_CONTROL_WITNESS')
        if 'case' in contract and 'mode' in contract:
            cases.setdefault(row[contract['case']], []).append(row[contract['mode']])
        case_key = row.get(contract.get('case', 'name'))
        category = contract.get('category_by_case', {}).get(case_key,
                    row.get(contract.get('category_field'), contract.get('category', 'formula_or_input_integrity')))
        categories[category] = categories.get(category, 0)+1
    if cases:
        need(all(len(modes) == 2 and set(modes) == {False, True} for modes in cases.values()), 'REJECTION_BOTH_PYTHON_MODES')
    return {'actual_rejected_subprocesses': len(rows), 'categories': categories}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--jobs', type=int, choices=(1, 2, 3, 4), default=3)
    parser.add_argument('--authenticate-only', action='store_true')
    args = parser.parse_args()
    repo, cap = args.repo_root.resolve(strict=True), HERE.parent
    out = args.output_dir.absolute()
    need(not args.output_dir.exists() and not args.output_dir.is_symlink(), 'FRESH_OUTPUT_DIRECTORY')
    need(all(not p.is_symlink() for p in out.parents), 'OUTPUT_ANCESTOR_SYMLINK')
    out = out.resolve()
    need(repo != out and repo not in out.parents and cap != out and cap not in out.parents, 'OUTPUT_OUTSIDE_SOURCE')
    own = authenticate(cap, args.manifest_sha256)
    plan = read(HERE/'REPLAY_PLAN.json')
    need(plan['schema'] == 'HDBLAST_PROPER_CLOCK_ROOT_REPLAY_V1', 'REPLAY_PLAN_SCHEMA')
    before = sources(repo, cap, plan)
    if args.authenticate_only:
        print(json.dumps({'status': 'PASS_PROPER_CLOCK_CLOSED_SOURCE_AUTHENTICATION', 'old_capsules': len(before)}, sort_keys=True))
        return
    # No project scientific module or third-party library was imported above.
    need(sys.version.split()[0] == plan['python_version'], 'PINNED_PYTHON_VERSION')
    env = dict(os.environ)
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    command = [sys.executable, '-I', '-B', '-c',
               'import json,sympy,mpmath,numpy;print(json.dumps({"sympy":sympy.__version__,"mpmath":mpmath.__version__,"numpy":numpy.__version__},sort_keys=True))']
    version_proc = subprocess.run(command, capture_output=True, text=True, env=env, timeout=30)
    need(version_proc.returncode == 0, 'PINNED_REQUIREMENTS_MISSING')
    versions = json.loads(version_proc.stdout)
    need(versions == plan['requirements'], 'PINNED_REQUIREMENTS_VERSION')
    out.mkdir(parents=True)
    (out/'logs').mkdir()
    bindings = {'CAP': str(cap), 'REPO': str(repo), 'OUT': str(out)}

    def resolve(v):
        return v.format_map(bindings)

    def execute(group):
        canonical, modes, rejected, diagnostics = {}, {}, {}, {}
        for step in group['steps']:
            flags = ['-I', '-B'] if step.get('third_party') else ['-I', '-S', '-B']
            argv = [sys.executable, *flags, *(['-O'] if step.get('optimized') else []),
                    *[resolve(x) for x in step['argv']]]
            if 'output' in step and step.get('create_output_parent', True):
                Path(resolve(step['output'])).parent.mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(argv, capture_output=True, env=env, cwd=repo,
                                  timeout=step.get('timeout_seconds', 600))
            log = out/'logs'/(group['name']+'-'+step['name'])
            log.with_suffix('.stdout').write_bytes(proc.stdout)
            log.with_suffix('.stderr').write_bytes(proc.stderr)
            need(proc.returncode == 0, 'REPLAY_STEP:'+group['name']+'/'+step['name']+':'+proc.stderr.decode(errors='replace')[-2000:])
            if step['kind'] == 'canonical':
                produced = Path(resolve(step['output']))
                if step.get('stdout_is_report'):
                    need(not produced.exists(), 'FRESH_STDOUT_REPORT')
                    produced.write_bytes(proc.stdout)
                expected = cap/relative(step['reference'])
                need(regular(produced) and sha(expected) == step['expected_sha256'] and
                     produced.read_bytes() == expected.read_bytes(), 'CANONICAL_BYTE_IDENTITY:'+step['name'])
                key = step['canonical_key']
                pin = sha(produced)
                need(key not in canonical or canonical[key] == pin, 'NORMAL_OPTIMIZED_CANONICAL_IDENTITY')
                canonical[key] = pin
                modes.setdefault(key, []).append(bool(step.get('optimized')))
                document = read(produced)
                for keypath, value in step.get('fields', {}).items():
                    need(dig(document, keypath) == value, 'CANONICAL_REQUIRED_FIELD:'+keypath)
                for keypath, count in step.get('inprocess_controls', {}).items():
                    value = dig(document, keypath)
                    need((len(value) if isinstance(value, (list, dict)) else value) == count, 'INPROCESS_CONTROL_ROSTER:'+keypath)
                    diagnostics[step['name']+':'+keypath] = count
            elif step['kind'] == 'controls':
                document = read(Path(resolve(step['output'])))
                rejected[step['name']] = rejection_contract(document, step['contract'])
                for keypath, count in step.get('positive_counts', {}).items():
                    need(len(dig(document, keypath)) == count, 'POSITIVE_CHECK_COUNT')
                    diagnostics[step['name']+':'+keypath] = count
                for keypath, count in step.get('scalar_check_counts', {}).items():
                    need(dig(document, keypath) == count, 'SEPARATE_SCALAR_CHECK_COUNT')
                    diagnostics[step['name']+':'+keypath] = count
            elif step['kind'] == 'success':
                need(step['stdout_contains'].encode() in proc.stdout, 'POSITIVE_STDOUT_WITNESS')
            else:
                raise ValueError('UNKNOWN_REPLAY_STEP_KIND')
            print(json.dumps({'stage': 'completed', 'group': group['name'], 'step': step['name']}, sort_keys=True), flush=True)
        need(all(len(v) == 2 and set(v) == {False, True} for v in modes.values()), 'EVERY_NEW_REPORT_BOTH_MODES')
        return {'canonical': canonical, 'rejections': rejected, 'separately_counted_checks': diagnostics}

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(execute, plan['groups']))
    count = sum(len(r['canonical']) for r in results)
    rejected = sum(v['actual_rejected_subprocesses'] for r in results for v in r['rejections'].values())
    need(count == plan['expected_new_canonical_report_pairs'], 'COMPLETE_NEW_REPORT_ROSTER')
    need(rejected == plan['expected_rejected_subprocesses'], 'COMPLETE_REJECTED_SUBPROCESS_ROSTER')
    need(authenticate(cap, own['manifest_sha256']) == own and sources(repo, cap, plan) == before, 'ALL_SOURCES_UNCHANGED_AFTER_REPLAY')
    result = {'schema': 'HDBLAST_PROPER_CLOCK_ROOT_REPLAY_RECEIPT_V1',
              'status': 'PASS_COMPLETE_PROPER_CLOCK_SCREEN_CHECKPOINT',
              'capsule_manifest_sha256': own['manifest_sha256'],
              'canonical_report_pairs_reproduced': count, 'actual_rejected_subprocesses': rejected,
              'all_sources_authenticated_before_scientific_imports_and_after': True,
              'inherited_closed_capsules_authenticated': before,
              'inherited_checkpoint_proof_replay_executed': False,
              'inherited_source_authentication_is_proof_replay': False,
              'new_source_campaign': {'actual_first_jet_calls_per_numeric_mode': 65, 'TT_new_source_oracle_calls': 0},
              'requirements': versions, 'python_version': plan['python_version'],
              'groups': {g['name']: r for g, r in zip(plan['groups'], results)},
              'scope': plan['scope']}
    (out/'VERIFICATION.json').write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'canonical_report_pairs': count,
                      'actual_rejected_subprocesses': rejected}, sort_keys=True))


if __name__ == '__main__':
    main()
