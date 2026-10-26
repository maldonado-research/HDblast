#!/usr/bin/env python3
"""Authenticate the complete screen magnetic curvature capsule, then replay it.

Only the Python standard library is used by this entry point.  The manifest
and archive SHA-256 values must be obtained outside the capsule.  Every run in
the authenticated REPLAY_PLAN.json uses isolated Python, without site loading
or bytecode writes.  Frozen upstream campaigns are authenticated, not rerun.
"""
import sys
sys.dont_write_bytecode = True

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import tempfile
import zipfile

CAPSULE_NAME = 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026'
MANIFEST_SCHEMA = 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_MANIFEST_V1'
PLAN_SCHEMA = 'HDBLAST_SCREEN_MAGNETIC_REPLAY_PLAN_V1'
UPSTREAM = {
    'exterior': ('HDBLAST_EXTERIOR_TENSOR_JETS_20261026', '0bc16f65598ae4977a7b8ab70a66c13ec85d844c871d6fee7739c4ec472d1c70'),
    'current': ('HDBLAST_COUPLED_MEAN_PILOT_20261025', 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de'),
    'I': ('HDBLAST_INPUT_PRECISION_20261023', 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61'),
    'P': ('HDBLAST_COMMON_GAUGE_FORCING_20261016', '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6'),
    'N': ('HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44'),
    'G': ('HDBLAST_GEOMETRIC_BULK_SOURCES_20261016', '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),
    'C': ('HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'),
    'T': ('HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023', 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'),
    'F': ('HDBLAST_FULL_ANGULAR_RESPONSE_20261026', 'a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac'),
    'precision35': ('HDBLAST_SHELL_FEEDBACK_PRECISION_20261025', '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4'),
}
ID = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def valid_sha(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def member_path(value):
    need(isinstance(value, str) and value and '\\' not in value and '\x00' not in value,
         'MANIFEST_MEMBER_PATH')
    path = PurePosixPath(value)
    need(not path.is_absolute() and path.as_posix() == value and ':' not in value
         and all(part not in ('.', '..') for part in path.parts)
         and value not in ('.', '..'), 'MANIFEST_MEMBER_PATH')
    return path


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, 'JSON_DUPLICATE_KEY')
        result[key] = value
    return result


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique_object)


def no_symlink_components(path, reason):
    absolute = Path(os.path.abspath(path))
    need(all(not part.is_symlink() for part in (absolute,) + tuple(absolute.parents)), reason)
    return absolute


def authenticate(root, pin, own=False):
    """Read and authenticate every file, without importing project Python."""
    need(valid_sha(pin), 'MANIFEST_PIN_FORMAT')
    root = no_symlink_components(root, 'ROOT_SYMLINK_OR_TYPE')
    need(root.is_dir(), 'ROOT_SYMLINK_OR_TYPE')
    manifest = root / 'MANIFEST.json'
    need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST_PIN')
    data = read_json(manifest)
    need(isinstance(data, dict) and isinstance(data.get('files'), list), 'MANIFEST_SCHEMA')
    if own:
        need(data.get('schema') == MANIFEST_SCHEMA, 'MANIFEST_SCHEMA')
    expected = {'MANIFEST.json'}
    total = manifest.stat().st_size
    for row in data['files']:
        need(isinstance(row, dict), 'MANIFEST_ENTRY')
        relative = row.get('path')
        member_path(relative)
        need(relative not in expected, 'MANIFEST_MEMBER_PATH')
        expected.add(relative)
        need(type(row.get('bytes')) is int and row['bytes'] >= 0 and valid_sha(row.get('sha256')),
             'MANIFEST_ENTRY_FIELDS')
        path = root / relative
        need(path.is_file() and not path.is_symlink(), 'MANIFEST_MEMBER_TYPE')
        need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'MANIFEST_MEMBER_BYTES')
        total += row['bytes']
    actual = set()
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TREE_SYMLINK')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            need(path.is_dir(), 'TREE_SPECIAL_MEMBER')
    need(actual == expected, 'MANIFEST_COMPLETE_CLOSURE')
    need(sha(manifest) == pin, 'MANIFEST_CHANGED_DURING_AUTHENTICATION')
    return {'manifest_sha256': pin, 'schema': data.get('schema', 'legacy_manifest_without_schema'),
            'files_including_manifest': len(expected), 'bytes_including_manifest': total}


def authenticate_upstream(repo):
    return {key: authenticate(repo / 'research' / name, pin)
            for key, (name, pin) in UPSTREAM.items()}


def archive_check(root, archive, pin):
    need(valid_sha(pin), 'ARCHIVE_PIN_FORMAT')
    archive = no_symlink_components(archive, 'ARCHIVE_SYMLINK')
    need(archive.is_file() and sha(archive) == pin, 'ARCHIVE_PIN')
    prefix = CAPSULE_NAME + '/'
    expected = {prefix + path.relative_to(root).as_posix(): path
                for path in root.rglob('*') if path.is_file()}
    with zipfile.ZipFile(archive) as packed:
        names = packed.namelist()
        need(len(names) == len(set(names)) and set(names) == set(expected), 'ARCHIVE_EXACT_MEMBERS')
        for entry in packed.infolist():
            member_path(entry.filename)
            mode = entry.external_attr >> 16
            need(not entry.is_dir() and not stat.S_ISLNK(mode)
                 and (stat.S_IFMT(mode) in (0, stat.S_IFREG)), 'ARCHIVE_MEMBER_TYPE')
            need(not entry.flag_bits & 1, 'ARCHIVE_ENCRYPTED_MEMBER')
            path = expected[entry.filename]
            need(entry.file_size == path.stat().st_size
                 and packed.read(entry) == path.read_bytes(), 'ARCHIVE_MEMBER_BYTES')
        need(packed.testzip() is None, 'ARCHIVE_CRC')
    need(sha(archive) == pin, 'ARCHIVE_CHANGED_DURING_AUTHENTICATION')
    return {'sha256': pin, 'members': len(expected), 'crc_pass': True,
            'exact_member_set': True, 'all_member_bytes_equal': True}


def guard_output(root, repo, output, archive):
    output = no_symlink_components(output, 'OUTPUT_DIRECTORY_SYMLINK')
    need(not output.exists(), 'FRESH_EXTERNAL_OUTPUT_DIRECTORY')
    need(output.parent.is_dir(), 'OUTPUT_PARENT_REQUIRED')
    resolved = output.resolve()
    protected = [root.resolve()] + [(repo / 'research' / name).resolve()
                                   for name, _ in UPSTREAM.values()]
    need(all(resolved != path and path not in resolved.parents for path in protected),
         'OUTPUT_INSIDE_PROTECTED_CAPSULE')
    need(resolved != archive.resolve(), 'OUTPUT_INPUT_COLLISION')
    return resolved


def output_path(output, relative, create_parents=False):
    member_path(relative)
    path = output / relative
    no_symlink_components(path, 'OUTPUT_SYMLINK')
    need(output.resolve() in path.resolve().parents, 'OUTPUT_ESCAPES_REPLAY_DIRECTORY')
    if create_parents:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


def isolated_command(script, optimized=False):
    return [sys.executable, '-I', '-S', '-B'] + (['-O'] if optimized else []) + [str(script)]


def output_arguments(arguments):
    """Recognize both argparse spellings of the supported output options."""
    result = []
    for index, argument in enumerate(arguments):
        for option in ('--output', '--output-dir'):
            if argument == option:
                need(index + 1 < len(arguments), 'REPLAY_PLAN_OUTPUT_ARGUMENT')
                result.append((option, Path(arguments[index + 1])))
            elif argument.startswith(option + '='):
                need(bool(argument[len(option) + 1:]), 'REPLAY_PLAN_OUTPUT_ARGUMENT')
                result.append((option, Path(argument[len(option) + 1:])))
    return result


def compact(value):
    print(json.dumps(value, sort_keys=True, separators=(',', ':')), flush=True)


def plan_check(root, output, repo):
    plan = read_json(root / 'REPLAY_PLAN.json')
    need(isinstance(plan, dict) and plan.get('schema') == PLAN_SCHEMA, 'REPLAY_PLAN_SCHEMA')
    runs = plan.get('runs')
    need(isinstance(runs, list) and bool(runs), 'REPLAY_PLAN_RUNS')
    tokens = {'output': str(output), 'repo': str(repo), 'capsule': str(root)}
    for key, relative in plan.get('pins', {}).items():
        need(isinstance(key, str) and ID.fullmatch(key) and key not in tokens
             and key + '_sha' not in tokens, 'REPLAY_PLAN_PIN_NAME')
        member_path(relative)
        source = root / relative
        need(source.is_file() and not source.is_symlink(), 'REPLAY_PLAN_PIN_MEMBER')
        tokens[key], tokens[key + '_sha'] = str(source), sha(source)
    ids, declared_outputs = set(), set()
    reserved_outputs = {'REPRODUCTION_RECEIPT.json', 'REPRODUCTION_FAILURE.json'}
    for run in runs:
        if isinstance(run, dict) and isinstance(run.get('id'), str):
            reserved_outputs.update(('logs/' + run['id'] + '.stdout.log',
                                     'logs/' + run['id'] + '.stderr.log'))
    for run in runs:
        need(isinstance(run, dict), 'REPLAY_PLAN_RUN')
        identity = run.get('id')
        need(isinstance(identity, str) and ID.fullmatch(identity) and identity not in ids,
             'REPLAY_PLAN_RUN_ID')
        ids.add(identity)
        script = run.get('script')
        member_path(script)
        need(script.endswith('.py') and (root / script).is_file(), 'REPLAY_PLAN_SCRIPT')
        need(run.get('mode') in ('normal', 'optimized'), 'REPLAY_PLAN_MODE')
        argv = run.get('argv', [])
        need(isinstance(argv, list) and all(isinstance(value, str) for value in argv), 'REPLAY_PLAN_ARGV')
        try:
            expanded = [value.format_map(tokens) for value in argv]
        except (KeyError, ValueError) as error:
            raise ValueError('REPLAY_PLAN_PLACEHOLDER:' + str(error)) from error
        # Output options must name fresh replay paths, even when a typo in an
        # authenticated plan would otherwise target a source/input file.
        for option, destination in output_arguments(expanded):
            need(destination.is_absolute() and output in destination.resolve().parents,
                 'REPLAY_PLAN_OUTPUT_COLLISION')
            no_symlink_components(destination, 'OUTPUT_SYMLINK')
        goldens = run.get('goldens', [])
        need(isinstance(goldens, list), 'REPLAY_PLAN_GOLDENS')
        reason = run.get('reject_reason')
        absent = run.get('absent', [])
        if reason is not None:
            need(isinstance(reason, str) and bool(reason) and not goldens
                 and isinstance(absent, list) and bool(absent), 'REPLAY_PLAN_NEGATIVE_CONTROL')
        else:
            need(bool(goldens), 'REPLAY_PLAN_POSITIVE_NEEDS_GOLDEN')
        for relative in absent:
            output_path(output, relative)
        produced = []
        for golden in goldens:
            need(isinstance(golden, dict), 'REPLAY_PLAN_GOLDEN')
            relative, reference = golden.get('output'), golden.get('reference')
            output_path(output, relative)
            need(relative not in reserved_outputs, 'REPLAY_PLAN_RESERVED_OUTPUT')
            member_path(reference)
            need((root / reference).is_file(), 'REPLAY_PLAN_REFERENCE')
            if 'sha256' in golden:
                need(golden['sha256'] == sha(root / reference), 'REPLAY_PLAN_GOLDEN_PIN')
            produced.append(relative)
        if 'stdout' in run:
            need(run['stdout'] in produced and reason is None, 'REPLAY_PLAN_STDOUT_GOLDEN')
        need(len(produced) == len(set(produced)), 'REPLAY_PLAN_DUPLICATE_OUTPUT')
        need(not declared_outputs.intersection(produced), 'REPLAY_PLAN_DUPLICATE_OUTPUT')
        declared_outputs.update(produced)
    counts = plan.get('expected_counts')
    need(isinstance(counts, dict) and all(type(counts.get(key)) is int and counts[key] >= 0
         for key in ('positive_runs', 'negative_runs', 'golden_matches')), 'REPLAY_PLAN_EXPECTED_COUNTS')
    expected = {'positive_runs': sum('reject_reason' not in run for run in runs),
                'negative_runs': sum('reject_reason' in run for run in runs),
                'golden_matches': sum(len(run.get('goldens', [])) for run in runs)}
    need(counts == expected, 'REPLAY_PLAN_COUNT_MISMATCH')
    for group in plan.get('parity_groups', []):
        need(isinstance(group, list) and len(group) >= 2 and len(group) == len(set(group))
             and all(identity in ids for identity in group), 'REPLAY_PLAN_PARITY_GROUP')
    roles = plan.get('program_roles')
    need(isinstance(roles, dict), 'REPLAY_PLAN_PROGRAM_ROLES')
    all_programs = {path.relative_to(root).as_posix() for path in root.rglob('*.py')}
    direct = {run['script'] for run in runs}
    need(set(roles) == all_programs - direct - {'reproduce.py'}, 'REPLAY_PLAN_PROGRAM_COVERAGE')
    for relative, row in roles.items():
        need(isinstance(row, dict) and row.get('role') in ('dependency', 'provenance')
             and isinstance(row.get('description'), str) and bool(row['description']),
             'REPLAY_PLAN_PROGRAM_ROLE')
        if row['role'] == 'dependency':
            need(isinstance(row.get('used_by'), list) and bool(row['used_by'])
                 and all(actor in direct for actor in row['used_by']), 'REPLAY_PLAN_DEPENDENCY_COVERAGE')
    return plan, tokens


class Replay:
    def __init__(self, root, repo, output, plan, tokens):
        self.root, self.repo, self.output = root, repo, output
        self.plan, self.tokens = plan, tokens
        self.runs, self.matches = [], []
        (output / 'logs').mkdir()

    def run(self, run):
        identity = run['id']
        optimized = run['mode'] == 'optimized'
        command = isolated_command(self.root / run['script'], optimized)
        command += [argument.format_map(self.tokens) for argument in run.get('argv', [])]
        stdout = output_path(self.output, run.get('stdout', 'logs/' + identity + '.stdout.log'), True)
        stderr = output_path(self.output, 'logs/' + identity + '.stderr.log', True)
        for row in run.get('goldens', []):
            produced = output_path(self.output, row['output'])
            need(not produced.exists() and not produced.is_symlink(), 'REPLAY_OUTPUT_EXISTS:' + identity)
        for option, destination in output_arguments(command):
            no_symlink_components(destination, 'OUTPUT_SYMLINK')
            need(not destination.exists(), 'REPLAY_OUTPUT_EXISTS:' + identity)
            destination.parent.mkdir(parents=True, exist_ok=True)
        for relative in run.get('absent', []):
            absent = output_path(self.output, relative, True)
            need(not absent.exists() and not absent.is_symlink(), 'NEGATIVE_OUTPUT_ALREADY_EXISTS:' + identity)
        compact({'phase': 'start', 'run': identity, 'mode': run['mode']})
        with stdout.open('xb') as out_stream, stderr.open('xb') as err_stream:
            result = subprocess.run(command, stdout=out_stream, stderr=err_stream,
                                    cwd=self.output, check=False)
        reason = run.get('reject_reason')
        if reason is None:
            need(result.returncode == 0, 'REPLAY_FAILED:' + identity + ':' + stderr.read_text(errors='replace')[-5000:])
        else:
            # Require the terminal exception marker, not a stray substring in
            # progress output, to establish the intended rejection gate.
            terminal = stderr.read_text(errors='replace').strip().splitlines()
            observed = terminal[-1] if terminal else ''
            need(result.returncode != 0 and observed in (reason, 'ValueError: ' + reason),
                 'WRONG_NEGATIVE_CONTROL:' + identity + ':' + observed)
            need(all(not output_path(self.output, relative).exists()
                     and not output_path(self.output, relative).is_symlink()
                     for relative in run['absent']), 'NEGATIVE_CONTROL_WROTE_OUTPUT:' + identity)
        record = {'id': identity, 'script': run['script'], 'mode': run['mode'],
                  'flags': ['-I', '-S', '-B'] + (['-O'] if optimized else []),
                  'exit_code': result.returncode,
                  'stdout': stdout.relative_to(self.output).as_posix(), 'stdout_sha256': sha(stdout),
                  'stderr': stderr.relative_to(self.output).as_posix(), 'stderr_sha256': sha(stderr)}
        if reason is not None:
            record.update(expected_rejection=reason, reason_verified=True, output_absent=True)
        self.runs.append(record)
        for row in run.get('goldens', []):
            produced, reference = self.output / row['output'], self.root / row['reference']
            need(produced.is_file() and not produced.is_symlink(), 'GOLDEN_OUTPUT_MISSING:' + identity)
            need(produced.read_bytes() == reference.read_bytes(), 'GOLDEN_BYTES:' + identity + ':' + row['reference'])
            match = {'run': identity, 'output': row['output'], 'reference': row['reference'],
                     'sha256': sha(produced), 'bytes': produced.stat().st_size, 'byte_identical': True}
            self.matches.append(match)
            compact(dict(phase='golden', **match))
        compact({'phase': 'pass', 'run': identity, 'expected_rejection': reason})

    def finish(self):
        for group in self.plan.get('parity_groups', []):
            outputs = [[self.output / row['output'] for run in self.plan['runs'] if run['id'] == identity
                        for row in run.get('goldens', [])] for identity in group]
            need(bool(outputs[0]) and all(len(rows) == len(outputs[0]) for rows in outputs), 'MODE_PARITY_OUTPUT_COUNT')
            need(all(left.read_bytes() == right.read_bytes() for rows in outputs[1:]
                     for left, right in zip(outputs[0], rows)), 'MODE_PARITY:' + ','.join(group))
        counts = {'positive_runs': sum('expected_rejection' not in row for row in self.runs),
                  'negative_runs': sum('expected_rejection' in row for row in self.runs),
                  'golden_matches': len(self.matches)}
        need(counts == self.plan['expected_counts'], 'REPLAY_FINAL_COUNTS')
        return counts


def integrity_controls(root, output):
    """Two copied source payloads, inert under a complete fixture manifest.

    Every rejection is performed in a fresh isolated authentication subprocess.
    The copied primary and dependency sources are never imported or executed.
    No control invokes the full replay entry point recursively.
    """
    programs = sorted(path for path in root.rglob('*.py') if path.name != 'reproduce.py')
    need(len(programs) >= 2, 'INTEGRITY_CONTROLS_NEED_TWO_SOURCES')
    rows = []
    with tempfile.TemporaryDirectory(prefix='integrity-fixture-', dir=output) as temporary:
        fixture = Path(temporary) / 'source'
        fixture.mkdir()
        members = [fixture / 'primary.py', fixture / 'dependency.py']
        for path, source in zip(members, programs[:2]):
            path.write_bytes(source.read_bytes())
        manifest = fixture / 'MANIFEST.json'
        manifest.write_text(json.dumps({'schema': MANIFEST_SCHEMA, 'files': [
            {'path': path.name, 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in members
        ]}, indent=2, sort_keys=True) + '\n', encoding='utf-8')
        original = manifest.read_bytes()
        pin = sha(manifest)
        checker = root / 'reproduce.py'

        def check(identity, reason=None, supplied_pin=None, optimized=False):
            command = isolated_command(checker, optimized) + ['--authenticate-tree', str(fixture),
                '--manifest-sha256', supplied_pin or pin]
            proc = subprocess.run(command, capture_output=True, text=True, cwd=output, check=False)
            terminal = proc.stderr.strip().splitlines()
            observed = terminal[-1] if terminal else ''
            if reason is None:
                need(proc.returncode == 0, 'FIXTURE_AUTHENTICATION_FAILED:' + observed)
            else:
                need(proc.returncode != 0 and observed == 'ValueError: ' + reason,
                     'WRONG_INTEGRITY_CONTROL:' + identity + ':' + observed)
                rows.append({'id': identity, 'mode': 'optimized' if optimized else 'normal',
                             'expected_rejection': reason, 'exit_code': proc.returncode,
                             'copied_project_sources_executed': False, 'reason_verified': True})
                compact({'phase': 'integrity-control', 'run': identity, 'reason': reason, 'pass': True})

        for optimized in (False, True):
            suffix = '-optimized' if optimized else '-normal'
            check('fixture-baseline' + suffix, optimized=optimized)
            check('wrong-manifest-pin' + suffix, 'MANIFEST_PIN', '0' * 64, optimized)
            for index, member in enumerate(members):
                saved = member.read_bytes()
                member.write_bytes(saved + b'\n# integrity control\n')
                check(('changed-primary-source' if index == 0 else 'changed-dependency-source') + suffix,
                      'MANIFEST_MEMBER_BYTES', optimized=optimized)
                member.write_bytes(saved)
            extra = fixture / 'UNDECLARED'
            extra.write_bytes(b'undeclared\n')
            check('undeclared-member' + suffix, 'MANIFEST_COMPLETE_CLOSURE', optimized=optimized)
            extra.unlink()
            saved = members[0].read_bytes()
            members[0].unlink()
            members[0].symlink_to(programs[0])
            check('source-symlink' + suffix, 'MANIFEST_MEMBER_TYPE', optimized=optimized)
            members[0].unlink()
            members[0].write_bytes(saved)
            for identity, mutate in (
                ('duplicate-member', lambda data: data['files'].append(data['files'][0].copy())),
                ('parent-traversal', lambda data: data['files'][0].update(path='../outside')),
                ('absolute-member', lambda data: data['files'][0].update(path='/outside')),
            ):
                data = json.loads(original)
                mutate(data)
                manifest.write_text(json.dumps(data), encoding='utf-8')
                check(identity + suffix, 'MANIFEST_MEMBER_PATH', sha(manifest), optimized)
                manifest.write_bytes(original)
            check('fixture-restored' + suffix, optimized=optimized)
    return rows


def archive_controls(root, archive, pin):
    """Reject incorrect external pins and correctly repinned byte changes."""
    need(sha(archive) == pin, 'ARCHIVE_PIN')
    rows = []
    for identity, expected_pin in (('wrong-archive-pin', '0' * 64),):
        try:
            archive_check(root, archive, expected_pin)
        except ValueError as error:
            need(str(error) == 'ARCHIVE_PIN', 'WRONG_ARCHIVE_CONTROL:' + str(error))
            rows.append({'id': identity, 'expected_rejection': 'ARCHIVE_PIN', 'reason_verified': True})
        else:
            raise ValueError('ARCHIVE_CONTROL_ACCEPTED:' + identity)
    # The second ZIP retains the member roster and valid CRCs but changes one
    # file.  Its freshly computed pin prevents an external-hash failure from
    # masking the all-member byte comparison being tested.
    with tempfile.TemporaryDirectory(prefix='archive-fixture-') as temporary:
        changed = Path(temporary) / 'changed.zip'
        with zipfile.ZipFile(archive) as source, zipfile.ZipFile(changed, 'x', compression=zipfile.ZIP_STORED) as target:
            for index, entry in enumerate(source.infolist()):
                payload = source.read(entry)
                if index == 0:
                    payload += b'!'
                target.writestr(entry, payload)
        try:
            archive_check(root, changed, sha(changed))
        except ValueError as error:
            need(str(error) == 'ARCHIVE_MEMBER_BYTES', 'WRONG_ARCHIVE_CONTROL:' + str(error))
            rows.append({'id': 'changed-archive-member', 'expected_rejection': 'ARCHIVE_MEMBER_BYTES',
                         'reason_verified': True})
        else:
            raise ValueError('ARCHIVE_CONTROL_ACCEPTED:changed-archive-member')
        incomplete = Path(temporary) / 'missing-members.zip'
        with zipfile.ZipFile(incomplete, 'x', compression=zipfile.ZIP_STORED):
            pass
        try:
            archive_check(root, incomplete, sha(incomplete))
        except ValueError as error:
            need(str(error) == 'ARCHIVE_EXACT_MEMBERS', 'WRONG_ARCHIVE_CONTROL:' + str(error))
            rows.append({'id': 'missing-archive-members', 'expected_rejection': 'ARCHIVE_EXACT_MEMBERS',
                         'reason_verified': True})
        else:
            raise ValueError('ARCHIVE_CONTROL_ACCEPTED:missing-archive-members')
    for row in rows:
        compact(dict(phase='integrity-control', **row, **{'pass': True}))
    return rows


def program_index(root, plan):
    result = []
    for path in sorted(root.rglob('*.py')):
        relative = path.relative_to(root).as_posix()
        row = {'path': relative, 'sha256': sha(path)}
        direct = [run['id'] for run in plan['runs'] if run['script'] == relative]
        if relative == 'reproduce.py':
            row.update(role='entry', description='Executing standard-library authentication and replay entry point.')
        elif direct:
            row.update(role='direct', run_ids=direct)
        else:
            row.update(plan['program_roles'][relative])
        result.append(row)
    return result


def output_index(output):
    result = []
    for path in sorted(output.rglob('*')):
        need(not path.is_symlink(), 'REPLAY_OUTPUT_SYMLINK')
        if path.is_file() and path.name not in ('REPRODUCTION_RECEIPT.json', 'REPRODUCTION_FAILURE.json'):
            result.append({'path': path.relative_to(output).as_posix(),
                           'bytes': path.stat().st_size, 'sha256': sha(path)})
        elif not path.is_file():
            need(path.is_dir(), 'REPLAY_OUTPUT_SPECIAL_MEMBER')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--repo-root', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--archive-sha256')
    parser.add_argument('--verify-only', action='store_true', help='Authenticate all files/ZIP/upstream roots; run no project code.')
    parser.add_argument('--authenticate-tree', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         'ENTRY_REQUIRES_ISOLATED_NO_SITE_NO_BYTECODE_FLAGS')
    entry = no_symlink_components(__file__, 'ENTRY_SOURCE_SYMLINK')
    entry_pin = sha(entry)
    if args.authenticate_tree is not None:
        need(args.repo_root is None and args.output_dir is None and args.archive is None
             and args.archive_sha256 is None and not args.verify_only, 'AUTHENTICATE_TREE_ARGUMENTS')
        result = authenticate(args.authenticate_tree, args.manifest_sha256, own=True)
        need(sha(entry) == entry_pin, 'ENTRY_SOURCE_CHANGED')
        compact({'status': 'PASS_TREE_AUTHENTICATION', 'tree': result, 'project_execution': False})
        return
    need(args.repo_root is not None and args.archive is not None and args.archive_sha256 is not None,
         'REPO_ARCHIVE_AND_EXTERNAL_ARCHIVE_PIN_REQUIRED')
    need(args.verify_only or args.output_dir is not None, 'FRESH_EXTERNAL_OUTPUT_DIRECTORY_REQUIRED')
    root = entry.parent
    need(root.name == CAPSULE_NAME, 'CAPSULE_ROOT_NAME')
    repo = no_symlink_components(args.repo_root, 'REPO_ROOT_SYMLINK')
    need(repo.is_dir(), 'REPO_ROOT_TYPE')
    output = guard_output(root, repo, args.output_dir, args.archive) if args.output_dir is not None else None
    before = authenticate(root, args.manifest_sha256, own=True)
    upstream_before = authenticate_upstream(repo)
    archive_before = archive_check(root, args.archive, args.archive_sha256)
    need(sha(entry) == entry_pin, 'ENTRY_SOURCE_CHANGED')
    if args.verify_only:
        after = authenticate(root, args.manifest_sha256, own=True)
        upstream_after = authenticate_upstream(repo)
        archive_after = archive_check(root, args.archive, args.archive_sha256)
        need(before == after and upstream_before == upstream_after and archive_before == archive_after
             and sha(entry) == entry_pin, 'SOURCE_CHANGED_DURING_VERIFICATION')
        compact({'status': 'PASS_COMPLETE_AUTHENTICATION', 'capsule': before,
                 'upstream_capsules': len(upstream_before), 'archive': archive_before,
                 'project_execution': False, 'source_unchanged': True})
        return
    plan, tokens = plan_check(root, output, repo)
    coverage = program_index(root, plan)
    output.mkdir()
    replay = Replay(root, repo, output, plan, tokens)
    try:
        source_controls = integrity_controls(root, output)
        zip_controls = archive_controls(root, args.archive, args.archive_sha256)
        # The original closure and archive are rechecked after the control
        # subprocesses and before the first scientific producer is launched.
        need(authenticate(root, args.manifest_sha256, own=True) == before
             and authenticate_upstream(repo) == upstream_before
             and archive_check(root, args.archive, args.archive_sha256) == archive_before
             and sha(entry) == entry_pin, 'SOURCE_CHANGED_BEFORE_PROJECT_EXECUTION')
        for run in plan['runs']:
            replay.run(run)
        counts = replay.finish()
        after = authenticate(root, args.manifest_sha256, own=True)
        upstream_after = authenticate_upstream(repo)
        archive_after = archive_check(root, args.archive, args.archive_sha256)
        need(before == after and upstream_before == upstream_after and archive_before == archive_after
             and sha(entry) == entry_pin, 'SOURCE_CHANGED_DURING_REPLAY')
        receipt = {'schema': 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_REPRODUCTION_V1', 'status': 'PASS',
            'python_version': sys.version, 'entry_mode': 'optimized' if sys.flags.optimize else 'normal',
            'entry_sha256': entry_pin, 'replay_plan_sha256': sha(root / 'REPLAY_PLAN.json'),
            'complete_capsule_before': before, 'complete_capsule_after': after,
            'complete_ten_upstream_roots_before': upstream_before,
            'complete_ten_upstream_roots_after': upstream_after,
            'all_authentication_before_any_project_execution': True,
            'archive_before': archive_before, 'archive_after': archive_after,
            'source_integrity_controls': source_controls, 'archive_integrity_controls': zip_controls,
            'direct_runs': replay.runs, 'golden_matches': replay.matches, 'counts': counts,
            'pass_ids': [row['id'] for row in replay.runs], 'program_index': coverage,
            'output_file_index': output_index(output), 'source_unchanged_after_replay': True,
            'runtime_dependencies': 'Python standard library only; no SymPy provenance generator is executed.',
            'scope': plan.get('scope', 'Only this new capsule is replayed; upstream full campaigns are not rerun.'),
            'limits': plan.get('limits', [])}
        with (output / 'REPRODUCTION_RECEIPT.json').open('x', encoding='utf-8') as stream:
            json.dump(receipt, stream, indent=2, sort_keys=True)
            stream.write('\n')
        compact({'status': 'PASS', **counts, 'source_integrity_controls': len(source_controls),
                 'archive_integrity_controls': len(zip_controls), 'upstream_capsules': len(UPSTREAM),
                 'manifest_sha256': args.manifest_sha256, 'archive_sha256': args.archive_sha256,
                 'receipt_sha256': sha(output / 'REPRODUCTION_RECEIPT.json')})
    except BaseException as error:
        # Failure receipts never replace the original exception or claim that
        # authentication succeeded after a failed scientific/control run.
        failure = {'status': 'FAIL', 'error': type(error).__name__ + ': ' + str(error),
                   'completed_runs': replay.runs, 'completed_golden_matches': replay.matches}
        with (output / 'REPRODUCTION_FAILURE.json').open('x', encoding='utf-8') as stream:
            json.dump(failure, stream, indent=2, sort_keys=True)
            stream.write('\n')
        raise


if __name__ == '__main__':
    main()
