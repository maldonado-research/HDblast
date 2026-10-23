#!/usr/bin/env python3
"""Replay the selected boundary certificate from authenticated read-only input.

The authenticated replay plan names each scientific role and its raw goldens.
Each role runs once normally and once with -O, using only the standard library.
Retained upstream campaigns are dependencies and are not executed again.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    result = {}
    for name, value in pairs:
        need(name not in result, 'DUPLICATE_JSON_KEY:' + name)
        result[name] = value
    return result


def decode(raw):
    return json.loads(raw, object_pairs_hook=unique_object)


def canonical(name):
    pure = PurePosixPath(name)
    need(name == pure.as_posix() and not pure.is_absolute() and name not in ('', '.')
         and '..' not in pure.parts and '\\' not in name, 'NONCANONICAL_PLAN_PATH')
    return name


def authenticate(expected):
    """Bootstrap closure before importing even our local auth/helper modules."""
    manifest_path = ROOT / 'MANIFEST.json'
    need(manifest_path.is_file() and not manifest_path.is_symlink(), 'SOURCE_MANIFEST_NOT_REGULAR')
    raw = manifest_path.read_bytes()
    need(sha(raw) == expected, 'SOURCE_MANIFEST_EXTERNAL_PIN')
    manifest = decode(raw)
    need(manifest['capsule'] == 'HDBLAST_BOUNDARY_RESPONSE_20261023', 'SOURCE_CAPSULE')
    expected_files = {}
    for row in manifest['files']:
        name = canonical(row['path'])
        need(name != 'MANIFEST.json' and name not in expected_files, 'SOURCE_DUPLICATE_MEMBER')
        expected_files[name] = row
    actual = set()
    for path in [ROOT] + sorted(ROOT.rglob('*')):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), 'SOURCE_NONREGULAR_MEMBER')
        need(not path.stat().st_mode & 0o222, 'SOURCE_WRITABLE_MEMBER')
        if path.is_file():
            name = path.relative_to(ROOT).as_posix(); actual.add(name)
            if name != 'MANIFEST.json':
                need(name in expected_files, 'SOURCE_UNEXPECTED_MEMBER:' + name)
                row = expected_files[name]; data = path.read_bytes()
                need(len(data) == row['bytes'] and sha(data) == row['sha256'], 'SOURCE_MEMBER_PIN:' + name)
    need(actual == set(expected_files) | {'MANIFEST.json'}, 'SOURCE_CLOSURE')


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    need(spec is not None and spec.loader is not None, 'HELPER_MODULE_LOADER')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def get_field(value, keys):
    for key in keys:
        value = value[key]
    return value


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode, 'USE_PYTHON_I_S_B')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest-sha256', required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    need(len(args.manifest_sha256) == 64 and all(c in '0123456789abcdef' for c in args.manifest_sha256),
         'INVALID_MANIFEST_HASH')
    output = args.output_dir.resolve()
    need(args.output_dir.is_absolute() and not output.exists() and not output.is_relative_to(ROOT),
         'FRESH_EXTERNAL_OUTPUT_DIRECTORY_REQUIRED')
    authenticate(args.manifest_sha256)
    auth = load('selected_source_auth', ROOT / 'source_auth.py')
    before = auth.authenticate(args.manifest_sha256)
    dependency = load('selected_dependency_chain', ROOT / 'dependency_chain.py')
    output.mkdir(parents=True)
    previous, old, chain_record = dependency.extract_chain(
        ROOT / 'inputs/FORCED_RESPONSE_DEPENDENCY.zip', output / 'retained-dependencies')
    old_before = dependency.inventory(old, dependency.GEOMETRIC)
    previous_before = dependency.inventory(previous, dependency.PREVIOUS)
    plan = decode((ROOT / 'REPLAY_PLAN.json').read_bytes())
    need(plan['capsule'] == before['capsule'] and isinstance(plan['jobs'], list) and plan['jobs'], 'REPLAY_PLAN')
    roles = set(); records = []; fresh = {}; controls = {}; joined = {}; join_records = []
    for job in plan['jobs']:
        role = canonical(job['role'])
        need('/' not in role and role not in roles, 'DUPLICATE_REPLAY_ROLE')
        roles.add(role)
        program = canonical(job['program'])
        mode_outputs = []
        for optimized in (False, True):
            mode = 'optimized' if optimized else 'normal'
            folder = output / (role + '_' + mode)
            folder.mkdir()
            replacements = {
                'SOURCE': str(ROOT), 'OLD': str(old), 'PREVIOUS': str(previous),
                'FOLLOWON': str(ROOT / 'inputs/followon'), 'JOB': str(folder),
                'MANIFEST': args.manifest_sha256,
            }
            for prior_role, paths in fresh.items():
                if optimized in paths:
                    replacements['FRESH_' + prior_role.upper()] = str(paths[optimized])
            replacements.update({name: str(path) for name, path in joined.items()})
            argv = []
            for argument in job['args']:
                need(type(argument) is str, 'REPLAY_ARGUMENT_TYPE')
                for key, value in replacements.items():
                    argument = argument.replace('{' + key + '}', value)
                need('{' not in argument and '}' not in argument, 'UNRESOLVED_REPLAY_ARGUMENT')
                argv.append(argument)
            command = [sys.executable] + (['-O'] if optimized else []) + [
                '-I', '-S', '-B', str(ROOT / program)] + argv
            completed = subprocess.run(command, cwd=output, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       timeout=600, check=False)
            (output / (role + '_' + mode + '.stdout.txt')).write_text(completed.stdout)
            (output / (role + '_' + mode + '.stderr.txt')).write_text(completed.stderr)
            need(completed.returncode == 0, 'REPLAY_CHILD_FAILED:' + role + ':' + mode + ':' + completed.stderr[-2400:])
            actual_files = set()
            for path in folder.rglob('*'):
                need(not path.is_symlink() and (path.is_dir() or path.is_file()), 'REPLAY_NONREGULAR_OUTPUT')
                if path.is_file():
                    actual_files.add(path.relative_to(folder).as_posix())
            outputs = job['outputs']
            need(actual_files == {canonical(row['fresh']) for row in outputs}, 'REPLAY_OUTPUT_CLOSURE:' + role)
            artifacts = []; raws = []
            for row in outputs:
                fresh_name = canonical(row['fresh'])
                golden_name = canonical(row[mode + '_golden'])
                path = folder / fresh_name
                raw = path.read_bytes(); saved = (ROOT / golden_name).read_bytes()
                need(raw == saved, 'REPLAY_RAW_GOLDEN_MISMATCH:' + role + ':' + fresh_name)
                if row['typed_json']:
                    need(decode(raw) == decode(saved), 'REPLAY_TYPED_GOLDEN_MISMATCH')
                raws.append(raw)
                artifacts.append({'output': path.relative_to(output).as_posix(), 'golden': golden_name,
                                  'bytes': len(raw), 'sha256': sha(raw), 'raw_golden_equal': True,
                                  'typed_json_golden_equal': row['typed_json']})
            mode_outputs.append(raws)
            primary = folder / canonical(job['primary_output'])
            need(primary.relative_to(folder).as_posix() in actual_files, 'PRIMARY_OUTPUT_MISSING')
            fresh.setdefault(role, {})[optimized] = primary
            for group in job.get('control_groups', []):
                selected = get_field(decode(primary.read_bytes()), group['field'])
                if type(selected) is list:
                    need(selected and all(type(row) is dict and row.get('rejected') is True for row in selected),
                         'TARGETED_REJECTION_OR_SYNTHETIC_CONTROL_STATUS')
                elif type(selected) is dict:
                    need(selected and all(type(name) is str and value is True for name, value in selected.items()),
                         'SPECIFIC_FORMULA_DISTINCTION_STATUS')
                else:
                    raise ValueError('CONTROL_GROUP_MUST_BE_REJECTION_LIST_OR_DISTINCTION_DICT')
                count = len(selected)
                label = role + ':' + group['name']
                controls.setdefault(label, {})[mode] = count
            records.append({'role': role, 'program': program, 'optimized': optimized,
                            'exit_code': completed.returncode, 'artifacts': artifacts})
        need(mode_outputs[0] == mode_outputs[1], 'REPLAY_EXECUTION_MODE_MISMATCH:' + role)
        if 'join_pair' in job:
            pair = job['join_pair']
            directory = canonical(pair['directory'])
            need('/' not in directory and pair['token'] not in joined, 'INVALID_FRESH_PAIR_JOIN')
            target_dir = output / directory
            target_dir.mkdir()
            artifacts = []
            for optimized in (False, True):
                mode = 'optimized' if optimized else 'normal'
                name = canonical(pair[mode])
                need('/' not in name, 'FRESH_PAIR_FILENAME')
                target = target_dir / name
                raw = fresh[role][optimized].read_bytes()
                with target.open('xb') as handle:
                    handle.write(raw)
                target.chmod(0o444)
                artifacts.append({'output': target.relative_to(output).as_posix(),
                                  'source': fresh[role][optimized].relative_to(output).as_posix(),
                                  'bytes': len(raw), 'sha256': sha(raw), 'mode': 0o444})
            target_dir.chmod(0o555)
            joined[pair['token']] = target_dir
            join_records.append({'role': role, 'token': pair['token'], 'directory': directory,
                                 'artifacts': artifacts, 'directory_mode': 0o555,
                                 'fresh_mode_pair_equal_and_read_only': True})
    need(dependency.inventory(old, dependency.GEOMETRIC) == old_before and
         dependency.inventory(previous, dependency.PREVIOUS) == previous_before, 'DEPENDENCY_CHANGED_DURING_REPLAY')
    need(auth.authenticate(args.manifest_sha256) == before, 'SOURCE_CHANGED_DURING_REPLAY')
    for join in join_records:
        folder = output / join['directory']
        need(not folder.is_symlink() and folder.stat().st_mode & 0o777 == join['directory_mode'], 'FRESH_PAIR_DIRECTORY_CHANGED')
        need({p.name for p in folder.iterdir()} == {Path(row['output']).name for row in join['artifacts']},
             'FRESH_PAIR_CLOSURE_CHANGED')
        for row in join['artifacts']:
            path = output / row['output']
            need(path.is_file() and not path.is_symlink() and path.stat().st_mode & 0o777 == row['mode']
                 and len(path.read_bytes()) == row['bytes'] and sha(path.read_bytes()) == row['sha256'],
                 'FRESH_PAIR_MEMBER_CHANGED')
    for row in controls.values():
        need(row['normal'] == row['optimized'], 'CONTROL_MODE_COUNT_MISMATCH')
    receipt = {'status': 'PASS_AUTHENTICATED_BOUNDARY_RESPONSE_REPRODUCTION',
               'capsule': before['capsule'], 'manifest_sha256': args.manifest_sha256,
               'files_including_manifest': before['files_including_manifest'],
               'top_level_children': len(records), 'roles': [row['role'] for row in plan['jobs']],
               'processes_per_role': {role: sum(row['role'] == role for row in records) for role in sorted(roles)},
               'targeted_in_process_controls': controls,
               'retained_dependencies': chain_record, 'children': records, 'fresh_report_joins': join_records,
               'source_read_only': True, 'source_unchanged': True,
               'retained_dependencies_read_only_and_unchanged': True,
               'upstream_scientific_producers_replayed': False,
               'scope': plan['scope']}
    target = output / 'REPRODUCTION_RECEIPT.json'
    with target.open('x') as handle:
        handle.write(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': receipt['status'], 'top_level_children': len(records),
                      'receipt_sha256': sha(target.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        sys.exit('BOUNDARY_RESPONSE_REPRODUCTION_FAILED:' + str(error))
