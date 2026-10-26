#!/usr/bin/env python3
"""Authenticate and replay only the new bounded full23 angular capsule.

Run the packaged copy at the capsule root. Scientific goldens must reproduce
byte for byte in both Python modes. All outputs and mutation copies remain
outside the source capsule and all eight frozen upstream roots.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile
import threading
import zipfile

sys.dont_write_bytecode = True
UPSTREAM = {
    'current': ('HDBLAST_COUPLED_MEAN_PILOT_20261025', 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de'),
    'I': ('HDBLAST_INPUT_PRECISION_20261023', 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61'),
    'P': ('HDBLAST_COMMON_GAUGE_FORCING_20261016', '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6'),
    'N': ('HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44'),
    'G': ('HDBLAST_GEOMETRIC_BULK_SOURCES_20261016', '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),
    'C': ('HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'),
    'T': ('HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023', 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'),
    'precision35': ('HDBLAST_SHELL_FEEDBACK_PRECISION_20261025', '1a5ee8aa22953f949da000ba03da3b5a0872e4467a0fbaec9f2d3fe7fba11ec4'),
}
PRIMARY_INITIAL = ('INITIAL_REPORT.json', 'INITIAL_PROFILES.json')
PRIMARY_RESPONSE = ('RESPONSE_REPORT.json', 'RESPONSE_PROFILES.json', 'EVALUATED_CANDIDATE.json')
MODES = (('normal', False), ('optimized', True))


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def digest_string(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def member_path(value):
    need(isinstance(value, str) and value and '\\' not in value and '\x00' not in value,
         'MANIFEST_MEMBER_PATH')
    path = PurePosixPath(value)
    need(value == path.as_posix() and not path.is_absolute()
         and value not in ('.', '..') and '..' not in path.parts and '.' not in path.parts,
         'MANIFEST_MEMBER_PATH')
    return path


def authenticate(root, pin, require_schema=True):
    """Complete file closure, checked without importing any project module."""
    need(digest_string(pin), 'MANIFEST_PIN_FORMAT')
    need(root.is_dir() and not root.is_symlink(), 'ROOT_SYMLINK_OR_TYPE')
    manifest = root/'MANIFEST.json'
    need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST_PIN')
    data = json.loads(manifest.read_text())
    need(isinstance(data, dict) and isinstance(data.get('files'), list), 'MANIFEST_SCHEMA')
    if require_schema:
        need(isinstance(data.get('schema'), str) and data['schema'], 'MANIFEST_SCHEMA')
    seen, total_bytes = {'MANIFEST.json'}, manifest.stat().st_size
    for row in data['files']:
        need(isinstance(row, dict), 'MANIFEST_ENTRY')
        relative = row.get('path')
        member_path(relative)
        need(relative not in seen, 'MANIFEST_MEMBER_PATH')
        seen.add(relative)
        need(type(row.get('bytes')) is int and row['bytes'] >= 0
             and digest_string(row.get('sha256')), 'MANIFEST_ENTRY_FIELDS')
        path = root/relative
        need(path.is_file() and not path.is_symlink(), 'MANIFEST_MEMBER_TYPE')
        need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'MANIFEST_MEMBER_BYTES')
        total_bytes += row['bytes']
    actual = {'MANIFEST.json'}
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TREE_SYMLINK')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            need(path.is_dir(), 'TREE_SPECIAL_MEMBER')
    need(actual == seen, 'MANIFEST_COMPLETE_CLOSURE')
    return {'manifest_sha256': pin, 'schema': data.get('schema', 'legacy_manifest_without_schema'),
            'files_including_manifest': len(seen), 'declared_bytes_including_manifest': total_bytes}


def authenticate_upstream(repo):
    return {key: authenticate(repo/'research'/name, pin, require_schema=False)
            for key, (name, pin) in UPSTREAM.items()}


def archive_check(root, archive, pin):
    need(archive.is_file() and digest_string(pin) and sha(archive) == pin, 'ARCHIVE_PIN')
    with zipfile.ZipFile(archive) as packed:
        need(packed.testzip() is None, 'ARCHIVE_CRC')
        names = packed.namelist()
        prefix = root.name+'/'
        expected = {prefix+path.relative_to(root).as_posix()
                    for path in root.rglob('*') if path.is_file()}
        need(len(names) == len(set(names)) and set(names) == expected, 'ARCHIVE_EXACT_MEMBERS')
        for name in names:
            member_path(name)
            need(packed.read(name) == (root/name[len(prefix):]).read_bytes(), 'ARCHIVE_MEMBER_BYTES')
    return {'sha256': pin, 'members': len(names), 'crc_pass': True, 'all_source_bytes_equal': True}


def corruption_controls(root, pin, output):
    """Seven precise failures on a disposable copy; no project execution."""
    with tempfile.TemporaryDirectory(prefix='capsule-corruption-', dir=output) as temporary:
        copy = Path(temporary)/'copy'
        shutil.copytree(root, copy)
        manifest = copy/'MANIFEST.json'
        original = manifest.read_bytes()
        rows = json.loads(original)['files']
        need(bool(rows), 'CORRUPTION_CONTROL_NEEDS_MEMBER')
        target = copy/rows[0]['path']
        saved = target.read_bytes()
        passed = []

        def reject(name, marker, supplied_pin=pin):
            try:
                authenticate(copy, supplied_pin)
            except ValueError as error:
                need(str(error) == marker, 'WRONG_CORRUPTION_REASON:'+name+':'+str(error))
                passed.append({'control': name, 'rejected': True, 'reason': marker})
            else:
                raise ValueError('CORRUPTION_ACCEPTED:'+name)

        reject('wrong_external_manifest_pin', 'MANIFEST_PIN', '0'*64)
        target.write_bytes(saved+b'!')
        reject('changed_member_bytes', 'MANIFEST_MEMBER_BYTES')
        target.write_bytes(saved)
        extra = copy/'_UNDECLARED_REPRODUCTION_CONTROL_'
        need(not extra.exists(), 'CORRUPTION_EXTRA_COLLISION')
        extra.write_bytes(b'undeclared')
        reject('undeclared_member', 'MANIFEST_COMPLETE_CLOSURE')
        extra.unlink()
        target.unlink()
        target.symlink_to(root/rows[0]['path'])
        reject('member_symlink', 'MANIFEST_MEMBER_TYPE')
        target.unlink()
        target.write_bytes(saved)
        for name, alteration in (
            ('duplicate_manifest_member', lambda data: data['files'].append(data['files'][0].copy())),
            ('parent_traversal', lambda data: data['files'][0].update(path='../outside')),
            ('absolute_member_path', lambda data: data['files'][0].update(path='/outside')),
        ):
            modified = json.loads(original)
            alteration(modified)
            manifest.write_text(json.dumps(modified))
            reject(name, 'MANIFEST_MEMBER_PATH', sha(manifest))
            manifest.write_bytes(original)
        need(authenticate(copy, pin) == authenticate(root, pin), 'RESTORED_CORRUPTION_COPY')
        return passed


class Replay:
    def __init__(self, root, repo, output):
        self.root, self.repo, self.output = root, repo, output
        self.runs, self.matches, self.nested = [], [], []
        self.lock = threading.Lock()
        (output/'logs').mkdir()

    def run(self, identity, script, arguments, optimized=False):
        path = self.root/script
        need(path.is_file(), 'MISSING_REPLAY_PROGRAM:'+script)
        command = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])
        command += [str(path)]+[str(value) for value in arguments]
        log = self.output/'logs'/(identity+'.log')
        with self.lock:
            print(json.dumps({'phase': 'start', 'run': identity, 'script': script,
                              'optimized': optimized}), flush=True)
        with log.open('x') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
        need(result.returncode == 0, 'REPLAY_FAILED:'+identity+':'+log.read_text()[-6000:])
        record = {'id': identity, 'script': script, 'optimized': optimized,
                  'exit_code': result.returncode, 'log': log.relative_to(self.output).as_posix(),
                  'log_sha256': sha(log)}
        with self.lock:
            self.runs.append(record)
            print(json.dumps({'phase': 'pass', 'run': identity}), flush=True)
        return record

    def golden(self, produced, reference, actor):
        source = self.root/reference
        need(produced.is_file() and source.is_file(), 'GOLDEN_MEMBER_MISSING:'+reference)
        need(produced.read_bytes() == source.read_bytes(), 'GOLDEN_BYTES:'+reference+':'+actor)
        self.matches.append({'actor': actor, 'output': produced.relative_to(self.output).as_posix(),
                             'golden': reference, 'sha256': sha(produced), 'bytes': produced.stat().st_size,
                             'byte_identical': True})

    def pair(self, identity, script, builder):
        # The two outputs are independent and have disjoint fresh destinations.
        with ThreadPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(self.run, identity+'-'+mode, script, builder(mode), optimized)
                       for mode, optimized in MODES]
            return [future.result() for future in futures]

    def produce(self):
        fresh = {}
        for family in ('initial', 'response'):
            for mode, _ in MODES:
                fresh[(family, mode)] = self.output/(family+'-'+mode)
        self.pair('initial', 'initial/full23_initial.py',
                  lambda mode: ['--repo-root', self.repo, '--output-dir', fresh[('initial', mode)]])
        for mode, _ in MODES:
            for name in PRIMARY_INITIAL:
                self.golden(fresh[('initial', mode)]/name, 'initial/results/'+name, 'initial-'+mode)
        self.pair('native-response', 'compute/full23_response.py',
                  lambda mode: ['--repo-root', self.repo, '--output-dir', fresh[('response', mode)]])
        for mode, _ in MODES:
            for name in PRIMARY_RESPONSE:
                self.golden(fresh[('response', mode)]/name, 'compute/response-normal/'+name, 'native-response-'+mode)
        tt = self.output/'tensor'
        tt.mkdir()
        self.pair('TT-total', 'tensor/tt_response.py',
                  lambda mode: ['--repo-root', self.repo, '--output', tt/('TT_TOTAL_'+mode+'.json')])
        for mode, _ in MODES:
            self.golden(tt/('TT_TOTAL_'+mode+'.json'),
                        'tensor/TT_TOTAL_RESPONSE_'+mode.upper()+'.json', 'TT-total-'+mode)
        return fresh

    def wrappers(self):
        controls = self.output/'structure-controls.json'
        self.run('structure-controls', 'structure/check_structure_controls.py',
                 ['--repo-root', self.repo, '--output', controls])
        self.golden(controls, 'structure/STRUCTURE_CONTROLS.json', 'structure-controls')
        result = json.loads(controls.read_text())
        need(result['normal_optimized_byte_identical'] and result['retained_outputs_equal_to_replay']
             and result['status'] == 'PASS_ALL_CONTROLS' and len(result['controls']) == 4,
             'STRUCTURE_WRAPPER_COVERAGE')
        # The authenticated wrapper itself checks both fresh temporary outputs
        # byte for byte against both archived scientific reports before deletion.
        for mode, _ in MODES:
            reference = 'structure/ANGULAR_STRUCTURE_'+mode.upper()+'.json'
            need(sha(self.root/reference) == result['output_sha256'], 'STRUCTURE_NESTED_GOLDEN')
            self.matches.append({'actor': 'structure-controls:'+mode,
                                 'golden': reference, 'sha256': result['output_sha256'],
                                 'byte_identical': True, 'comparison': 'fresh temporary producer bytes checked by authenticated control wrapper'})
        self.nested.append({'wrapper': 'structure/check_structure_controls.py',
                            'program': 'structure/exact_angular_structure.py',
                            'positive_modes': ['normal', 'optimized'], 'positive_runs': 2,
                            'negative_modes': ['optimized'], 'negative_runs': 4,
                            'controls': result['controls']})
        for identity, script, directory, report_prefix, archived, extra in (
            ('shell-controls', 'shell/verify_shell_audit.py', 'shell-checks', 'FULL23_SHELL',
             'shell/verified', []),
            ('kernel-controls', 'shell/kernel-review/verify_kernel_review.py', 'kernel-checks', 'ANGULAR_RIDGE',
             'shell/kernel-review/verified-final', ['--work-root', self.root]),
        ):
            destination = self.output/directory
            self.run(identity, script, ['--repo-root', self.repo]+extra+['--output-dir', destination])
            receipt = json.loads((destination/'EXECUTION_RECEIPT.json').read_text())
            need(receipt['status'] == 'PASS' and receipt['normal_optimized_byte_identical']
                 and len(receipt['positive_runs']) == 2 and len(receipt['negative_controls']) == 8,
                 'WRAPPER_BOTH_MODES_AND_MUTANTS:'+script)
            for mode, _ in MODES:
                name = report_prefix+'_'+mode.upper()+'.json'
                self.golden(destination/name, archived+'/'+name, identity+':'+mode)
            producer = ('shell/check_full23_shell.py' if identity == 'shell-controls'
                        else 'shell/kernel-review/check_angular_ridge.py')
            self.nested.append({'wrapper': script, 'program': producer,
                                'positive_modes': ['normal', 'optimized'], 'positive_runs': 2,
                                'negative_modes': ['normal', 'optimized'], 'negative_runs': 8,
                                'controls': receipt['negative_controls']})

    def reviews(self, fresh):
        reviews = self.output/'reviews'
        reviews.mkdir()
        for mode, optimized in MODES:
            folder = reviews/mode
            folder.mkdir()
            jobs = [
                ('initial-exact', 'initial/review/check_exact_adapter.py',
                 ['--repo-root', self.repo], 'EXACT_ADAPTER_REVIEW.json', 'initial/review/EXACT_ADAPTER_REVIEW.json'),
                ('initial-numeric', 'initial/review/check_numeric_report.py',
                 ['--repo-root', self.repo, '--results-dir', fresh[('initial', mode)]],
                 'NUMERIC_REPORT_REVIEW.json', 'initial/review/NUMERIC_REPORT_REVIEW.json'),
                ('native-independent', 'compute/review-response/check_response.py',
                 ['--repo-root', self.repo, '--results-dir', fresh[('response', mode)]],
                 'RESPONSE_REVIEW.json', 'compute/review-response/RESPONSE_REVIEW.json'),
                ('coefficient-independent', 'structure/response-review/check_response_coefficients.py',
                 ['--repo-root', self.repo, '--capsule-root', self.root,
                  '--normal-results', fresh[('response', 'normal')], '--optimized-results', fresh[('response', 'optimized')]],
                 'COEFFICIENT_SAMPLE_CHECKS.json', 'structure/response-review/COEFFICIENT_SAMPLE_CHECKS.json'),
                ('scope-independent', 'structure/response-review/check_scope.py',
                 ['--capsule-root', self.root, '--normal-results', fresh[('response', 'normal')],
                  '--optimized-results', fresh[('response', 'optimized')]],
                 'SCOPE_CHECKS.json', 'structure/response-review/SCOPE_CHECKS.json'),
                ('TT-independent', 'tensor/review/check_tt_total_report.py',
                 ['--repo-root', self.repo, '--report', self.output/'tensor'/('TT_TOTAL_'+mode+'.json'),
                  '--evaluator', self.root/'tensor/tt_response.py'],
                 'TT_TOTAL_REVIEW.json', 'tensor/review/TT_TOTAL_REVIEW_'+mode.upper()+'.json'),
                ('TT-tail-geometry', 'tensor/review/check_tt_tail_geometry.py',
                 ['--repo-root', self.repo], 'TT_TAIL_GEOMETRY_REVIEW.json',
                 'tensor/review/TT_TAIL_GEOMETRY_REVIEW_'+mode.upper()+'.json'),
            ]
            for identity, script, arguments, filename, reference in jobs:
                produced = folder/filename
                self.run(identity+'-'+mode, script, arguments+['--output', produced], optimized)
                self.golden(produced, reference, identity+'-'+mode)
        guard = reviews/'TT_EVALUATOR_GUARD_REVIEW.json'
        self.run('TT-io-guards', 'tensor/review/check_tt_io_guards.py',
                 ['--repo-root', self.repo, '--report', self.output/'tensor/TT_TOTAL_normal.json',
                  '--evaluator', self.root/'tensor/tt_response.py', '--output', guard])
        for mode, _ in MODES:
            self.golden(guard, 'tensor/review/TT_EVALUATOR_GUARD_REVIEW_'+mode.upper()+'.json', 'TT-io-guards')
        data = json.loads(guard.read_text())
        need(data['status'] == 'PASS_ALL_EVALUATOR_IO_GUARDS' and data['rejected_case_count'] == 20,
             'TT_GUARD_BOTH_MODES')
        self.nested.append({'wrapper': 'tensor/review/check_tt_io_guards.py', 'program': 'tensor/tt_response.py',
                            'negative_modes': ['normal', 'optimized'], 'negative_runs': 20,
                            'controls': data['rejected_cases'], 'source_oracle_executed': False})

    def program_index(self):
        direct = {}
        for run in self.runs:
            direct.setdefault(run['script'], []).append(run['id'])
        nested = {}
        for row in self.nested:
            nested.setdefault(row['program'], []).append(row)
        ancillary = {
            'initial/review/check_output.py': 'Ancillary saved support review; fresh primary byte equality and independent numeric review cover its scientific checks. Its legacy default writes into source, so it is not executed.',
            'initial/review/check_nodal_norm.py': 'Ancillary point/nodal comparison, outside the final scientific target ledger; archived result authenticated. Its legacy default writes into source, so it is not executed.',
            'initial/review/check_mode_parity.py': 'Superseded here by mandatory fresh normal/optimized primary byte comparisons; source-layout-specific optimized archive is not required.',
            'compute/review-response/check_mode_parity.py': 'Superseded here by mandatory fresh normal/optimized response byte comparisons.',
            'tensor/review/check_tt_portability.py': 'Historical pre-portability comparison needs external historical reference files. Its accepted provenance report is authenticated; current portable canonical outputs are freshly reproduced here.',
        }
        index = []
        for path in sorted(self.root.rglob('*.py')):
            relative = path.relative_to(self.root).as_posix()
            row = {'path': relative, 'sha256': sha(path)}
            if relative == 'reproduce.py':
                row['coverage'] = 'Executing capsule entry point; complete manifest/upstream authentication, archive checks and seven source-corruption controls.'
            elif relative in direct or relative in nested:
                row['direct_run_ids'] = sorted(direct.get(relative, []))
                row['nested_invocations'] = nested.get(relative, [])
                row['coverage'] = 'Executed directly and/or by the authenticated indexed control wrappers.'
            elif relative == 'compute/angular_ridge.py':
                row['coverage'] = 'Library imported by fresh native response producers and independent full23 kernel/sample reviews; not a standalone CLI.'
            elif relative in ancillary:
                row['coverage'] = ancillary[relative]
            else:
                raise ValueError('UNINDEXED_CAPSULE_PROGRAM:'+relative)
            index.append(row)
        return index


def output_index(output):
    return [{'path': path.relative_to(output).as_posix(), 'bytes': path.stat().st_size, 'sha256': sha(path)}
            for path in sorted(output.rglob('*')) if path.is_file()
            and path.name not in ('REPRODUCTION_RECEIPT.json', 'REPRODUCTION_FAILURE.json')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--archive-sha256')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    repo = args.repo_root.resolve()
    output = args.output_dir.resolve()
    protected = [root]+[(repo/'research'/name).resolve() for name, _ in UPSTREAM.values()]
    need(not args.output_dir.exists() and not args.output_dir.is_symlink()
         and all(output != path and path not in output.parents for path in protected),
         'FRESH_EXTERNAL_OUTPUT_DIRECTORY')
    need(bool(args.archive) == bool(args.archive_sha256), 'ARCHIVE_AND_EXTERNAL_PIN_REQUIRED_TOGETHER')
    before = authenticate(root, args.manifest_sha256)  # No project import/execution precedes this.
    upstream_before = authenticate_upstream(repo)  # All eight complete roots, including precision35.
    archive = archive_check(root, args.archive, args.archive_sha256) if args.archive else None
    output.mkdir(parents=True)
    replay = Replay(root, repo, output)
    try:
        controls = corruption_controls(root, args.manifest_sha256, output)
        fresh = replay.produce()
        replay.wrappers()
        replay.reviews(fresh)
        index = replay.program_index()
        after = authenticate(root, args.manifest_sha256)
        upstream_after = authenticate_upstream(repo)
        need(before == after and upstream_before == upstream_after, 'SOURCE_CHANGED_DURING_REPLAY')
        receipt = {
            'schema': 'HDBLAST_FULL23_ANGULAR_REPRODUCTION_V1', 'status': 'PASS',
            'python_version': sys.version, 'complete_new_capsule_before': before,
            'complete_new_capsule_after': after, 'complete_eight_upstream_roots_before': upstream_before,
            'complete_eight_upstream_roots_after': upstream_after, 'archive': archive,
            'source_corruption_controls': controls,
            'full_new_bounded_initial_native_response_and_TT_numerics': True,
            'scientific_normal_optimized_golden_matches_mandatory': True,
            'golden_matches': sorted(replay.matches, key=lambda row: (row['actor'], row['golden'])),
            'direct_runs': sorted(replay.runs, key=lambda row: row['id']),
            'nested_control_programs': replay.nested, 'complete_program_index': index,
            'output_file_index': output_index(output), 'source_unchanged_after_replay': True,
            'log_scope': 'Logs and execution receipts may contain environment paths or timing. Scientific JSONs and canonical independent reviews are compared byte for byte. Wrapper execution receipts are checked structurally, not byte-compared.',
            'replay_scope': 'Only new bounded producers, independent checks and declared negative controls were executed. Frozen upstream functions are imported through their authenticated closures; no old full campaign producer was rerun.',
            'limits': ['Conditional registered-gauge formal second response and coefficient bounds.',
                       'No C2 pointwise Einstein residual, nonlinear family, physical origin or novelty claim.',
                       'Historical portability/elapsed-time provenance is authenticated rather than re-created.'],
        }
        with (output/'REPRODUCTION_RECEIPT.json').open('x') as handle:
            handle.write(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
        print(json.dumps({'status': 'PASS', 'direct_runs': len(replay.runs),
                          'scientific_golden_matches': len(replay.matches), 'source_corruption_controls': len(controls),
                          'indexed_programs': len(index)}), flush=True)
    except BaseException as error:
        failure = {'status': 'FAIL', 'error': str(error), 'completed_direct_runs': replay.runs,
                   'completed_golden_matches': replay.matches}
        try:
            failure['capsule_unchanged_after_failure'] = authenticate(root, args.manifest_sha256) == before
            failure['upstream_unchanged_after_failure'] = authenticate_upstream(repo) == upstream_before
        except Exception as auth_error:
            failure['post_failure_authentication_error'] = str(auth_error)
        (output/'REPRODUCTION_FAILURE.json').write_text(json.dumps(failure, indent=2, sort_keys=True)+'\n')
        raise


if __name__ == '__main__':
    main()
