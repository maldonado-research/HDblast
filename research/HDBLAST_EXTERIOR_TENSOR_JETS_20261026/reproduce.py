#!/usr/bin/env python3
"""Authenticate and replay the new bounded exterior tensor two-jet capsule.

Run the packaged copy at the capsule root.  An external manifest SHA-256 is
mandatory.  Scientific JSONs must reproduce byte for byte in normal and
optimized isolated standard-library Python.  No frozen old campaign is run.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import zipfile

sys.dont_write_bytecode = True
MANIFEST_SCHEMA = 'HDBLAST_EXTERIOR_TENSOR_JETS_MANIFEST_V1'
CAPSULE_NAME = 'HDBLAST_EXTERIOR_TENSOR_JETS_20261026'
UPSTREAM = {
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
MODES = (('normal', False), ('optimized', True))
PROVENANCE_PROGRAMS = {
    'review-inputs/generate_exact_fixture.py': 'SymPy provenance generator. Its exact gradient fixture is authenticated and checked by the standard-library oracle review; this generator is never executed during replay.',
    'review-inputs/dependency_review.py': 'SymPy provenance derivation of the complete source dependency and principal symbol. Its accepted normal/optimized reports are authenticated; the independent standard-library dependency/theory audits and oracle fixture checks are freshly executed.',
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def digest_string(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def member_path(value):
    need(isinstance(value, str) and value and '\\' not in value and '\x00' not in value, 'MANIFEST_MEMBER_PATH')
    path = PurePosixPath(value)
    need(value == path.as_posix() and not path.is_absolute() and value not in ('.','..')
         and '..' not in path.parts and '.' not in path.parts, 'MANIFEST_MEMBER_PATH')
    return path


def authenticate(root, pin, new_capsule=False):
    """Complete file closure; no project module is imported or executed."""
    root = Path(root)
    need(digest_string(pin), 'MANIFEST_PIN_FORMAT')
    need(root.is_dir() and not root.is_symlink(), 'ROOT_SYMLINK_OR_TYPE')
    manifest = root/'MANIFEST.json'
    need(manifest.is_file() and not manifest.is_symlink() and sha(manifest) == pin, 'MANIFEST_PIN')
    data = json.loads(manifest.read_text())
    need(isinstance(data, dict) and isinstance(data.get('files'), list), 'MANIFEST_SCHEMA')
    if new_capsule:
        need(data.get('schema') == MANIFEST_SCHEMA, 'MANIFEST_SCHEMA')
    seen, total = {'MANIFEST.json'}, manifest.stat().st_size
    for row in data['files']:
        need(isinstance(row,dict), 'MANIFEST_ENTRY')
        relative = row.get('path')
        member_path(relative)
        need(relative not in seen, 'MANIFEST_MEMBER_PATH')
        seen.add(relative)
        need(type(row.get('bytes')) is int and row['bytes'] >= 0 and digest_string(row.get('sha256')), 'MANIFEST_ENTRY_FIELDS')
        path = root/relative
        need(path.is_file() and not path.is_symlink(), 'MANIFEST_MEMBER_TYPE')
        need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'MANIFEST_MEMBER_BYTES')
        total += row['bytes']
    actual = {'MANIFEST.json'}
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'TREE_SYMLINK')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            need(path.is_dir(), 'TREE_SPECIAL_MEMBER')
    need(actual == seen, 'MANIFEST_COMPLETE_CLOSURE')
    return {'manifest_sha256':pin, 'schema':data.get('schema','legacy_manifest_without_schema'),
            'files_including_manifest':len(seen), 'declared_bytes_including_manifest':total}


def authenticate_upstream(repo):
    return {key:authenticate(repo/'research'/name,pin) for key,(name,pin) in UPSTREAM.items()}


def guard_output(root, repo, output):
    output = Path(output)
    need(not output.is_symlink(), 'OUTPUT_DIRECTORY_SYMLINK')
    need(not output.exists(), 'FRESH_EXTERNAL_OUTPUT_DIRECTORY')
    target = output.resolve()
    protected = [root.resolve()]+[(repo/'research'/name).resolve() for name,_ in UPSTREAM.values()]
    need(all(target != path and path not in target.parents for path in protected), 'OUTPUT_INSIDE_PROTECTED_CAPSULE')
    return target


def archive_check(root, archive, pin):
    need(archive.is_file() and not archive.is_symlink() and digest_string(pin) and sha(archive) == pin, 'ARCHIVE_PIN')
    with zipfile.ZipFile(archive) as packed:
        need(packed.testzip() is None, 'ARCHIVE_CRC')
        names = packed.namelist()
        prefix = root.name+'/'
        expected = {prefix+path.relative_to(root).as_posix() for path in root.rglob('*') if path.is_file()}
        need(len(names) == len(set(names)) and set(names) == expected, 'ARCHIVE_EXACT_MEMBERS')
        for entry in packed.infolist():
            member_path(entry.filename)
            need(not entry.is_dir() and not stat.S_ISLNK(entry.external_attr >> 16), 'ARCHIVE_MEMBER_TYPE')
            need(packed.read(entry.filename) == (root/entry.filename[len(prefix):]).read_bytes(), 'ARCHIVE_MEMBER_BYTES')
    return {'sha256':pin, 'members':len(names), 'crc_pass':True, 'all_source_bytes_equal':True}


def corruption_controls(root, pin, output):
    """Seven reason-specific rejections on a disposable external source copy."""
    with tempfile.TemporaryDirectory(prefix='capsule-corruption-',dir=output) as temporary:
        copy = Path(temporary)/'copy'
        shutil.copytree(root,copy)
        manifest = copy/'MANIFEST.json'
        original = manifest.read_bytes()
        rows = json.loads(original)['files']
        need(bool(rows), 'CORRUPTION_CONTROL_NEEDS_MEMBER')
        target = copy/rows[0]['path']
        saved = target.read_bytes()
        passed = []

        def reject(name, marker, supplied_pin=pin):
            try:
                authenticate(copy,supplied_pin,new_capsule=True)
            except ValueError as error:
                need(str(error) == marker, 'WRONG_CORRUPTION_REASON:'+name+':'+str(error))
                passed.append({'control':name,'rejected':True,'reason':marker})
            else:
                raise ValueError('CORRUPTION_ACCEPTED:'+name)

        reject('wrong_external_manifest_pin','MANIFEST_PIN','0'*64)
        target.write_bytes(saved+b'!')
        reject('changed_member_bytes','MANIFEST_MEMBER_BYTES')
        target.write_bytes(saved)
        extra = copy/'_UNDECLARED_REPRODUCTION_CONTROL_'
        need(not extra.exists(), 'CORRUPTION_EXTRA_COLLISION')
        extra.write_bytes(b'undeclared')
        reject('undeclared_member','MANIFEST_COMPLETE_CLOSURE')
        extra.unlink()
        target.unlink()
        target.symlink_to(root/rows[0]['path'])
        reject('member_symlink','MANIFEST_MEMBER_TYPE')
        target.unlink()
        target.write_bytes(saved)
        for name,change in (
            ('duplicate_manifest_member',lambda data:data['files'].append(data['files'][0].copy())),
            ('parent_traversal',lambda data:data['files'][0].update(path='../outside')),
            ('absolute_member_path',lambda data:data['files'][0].update(path='/outside')),
        ):
            data = json.loads(original)
            change(data)
            manifest.write_text(json.dumps(data))
            reject(name,'MANIFEST_MEMBER_PATH',sha(manifest))
            manifest.write_bytes(original)
        need(authenticate(copy,pin,new_capsule=True) == authenticate(root,pin,new_capsule=True), 'RESTORED_CORRUPTION_COPY')
        return passed


class Replay:
    def __init__(self,root,repo,output):
        self.root,self.repo,self.output = root,repo,output
        self.runs,self.matches,self.controls = [],[],[]
        self.lock = threading.Lock()
        self.input = root/'inputs/actual_third_jets.py'
        self.evaluator = root/'compute/tt_twojet.py'
        self.input_sha,self.evaluator_sha = sha(self.input),sha(self.evaluator)
        (output/'logs').mkdir()

    def run(self,identity,script,arguments,optimized=False,reject_reason=None,absent=None):
        path = self.root/script
        need(path.is_file(), 'MISSING_REPLAY_PROGRAM:'+script)
        command = [sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(path)]+list(map(str,arguments))
        log = self.output/'logs'/(identity+'.log')
        with self.lock:
            print(json.dumps({'phase':'start','run':identity,'script':script,'optimized':optimized}),flush=True)
        with log.open('x') as handle:
            result = subprocess.run(command,stdout=handle,stderr=subprocess.STDOUT,text=True)
        content = log.read_text()
        if reject_reason is None:
            need(result.returncode == 0, 'REPLAY_FAILED:'+identity+':'+content[-6000:])
        else:
            need(result.returncode != 0 and reject_reason in content, 'WRONG_NEGATIVE_CONTROL:'+identity+':'+content[-6000:])
            if absent is not None:
                need(not absent.exists() and not absent.is_symlink(), 'NEGATIVE_CONTROL_WROTE_OUTPUT:'+identity)
        record = {'id':identity,'script':script,'optimized':optimized,'exit_code':result.returncode,
                  'log':log.relative_to(self.output).as_posix(),'log_sha256':sha(log)}
        if reject_reason is not None:
            record.update(expected_rejection=reject_reason,reason_verified=True,output_absent=True)
        with self.lock:
            self.runs.append(record)
            if reject_reason is not None:
                self.controls.append(record)
            print(json.dumps({'phase':'pass','run':identity,'expected_rejection':reject_reason}),flush=True)
        return record

    def golden(self,produced,reference,actor):
        source = self.root/reference
        need(produced.is_file() and source.is_file(), 'GOLDEN_MEMBER_MISSING:'+reference)
        need(produced.read_bytes() == source.read_bytes(), 'GOLDEN_BYTES:'+reference+':'+actor)
        row = {'actor':actor,'output':produced.relative_to(self.output).as_posix(),'golden':reference,
               'sha256':sha(produced),'bytes':produced.stat().st_size,'byte_identical':True}
        self.matches.append(row)
        # Native CI logs retain every comparison even when remote artifact
        # publishing is unavailable.  These are exact hashes, not display floats.
        print(json.dumps({'phase':'golden','actor':actor,'golden':reference,'sha256':row['sha256'],
                          'bytes':row['bytes'],'byte_identical':True}),flush=True)

    def pair(self,identity,script,arguments,references):
        outputs = {mode:self.output/(identity+'-'+mode+'.json') for mode,_ in MODES}
        with ThreadPoolExecutor(max_workers=2) as executor:
            jobs = [executor.submit(self.run,identity+'-'+mode,script,
                                    arguments(mode,outputs[mode]),optimized) for mode,optimized in MODES]
            for job in jobs:
                job.result()
        for mode,_ in MODES:
            self.golden(outputs[mode],references[mode],identity+'-'+mode)
        need(outputs['normal'].read_bytes() == outputs['optimized'].read_bytes(), 'MODE_PARITY:'+identity)
        return outputs

    def produce(self):
        self.pair('source-Fr','inputs/actual_third_jets.py',
            lambda mode,out:['--repo-root',self.repo,'--output',out],
            {mode:'inputs/actual_source_radial_final_'+mode+'.json' for mode,_ in MODES})
        self.pair('source-dependency','inputs/c4_need_audit.py',
            lambda mode,out:['--repo-root',self.repo,'--output',out],
            {mode:'inputs/c4_need_audit_'+mode+'.json' for mode,_ in MODES})
        self.pair('source-output-guards','inputs/check_output_guards.py',
            lambda mode,out:['--repo-root',self.repo,'--output',out],
            {mode:'inputs/output_guards_'+mode+'.json' for mode,_ in MODES})
        numeric = self.pair('actual-twojet','compute/tt_twojet.py',
            lambda mode,out:['--repo-root',self.repo,'--input-module',self.input,
                             '--input-module-sha256',self.input_sha,'--output',out],
            {mode:'compute/TT_TWOJET_'+mode.upper()+'.json' for mode,_ in MODES})
        self.pair('polynomial-controls','compute/check_polynomial_controls.py',
            lambda mode,out:['--repo-root',self.repo,'--evaluator',self.evaluator,
                             '--evaluator-sha256',self.evaluator_sha,'--output',out],
            {'normal':'compute/POLYNOMIAL_CONTROLS.json','optimized':'compute/POLYNOMIAL_CONTROLS_OPTIMIZED.json'})
        self.pair('theory-controls','theory/check_theory_controls.py',
            lambda mode,out:['--repo-root',self.repo,'--output',out],
            {mode:'theory/THEORY_CONTROLS_'+mode.upper()+'.json' for mode,_ in MODES})
        return numeric

    def reviews(self,numeric):
        self.pair('independent-source','review-inputs/oracle_review.py',
            lambda mode,out:['--repo-root',self.repo,'--output',out],
            {mode:'review-inputs/ORACLE_FINAL_'+mode.upper()+'.json' for mode,_ in MODES})
        self.pair('independent-twojet','review-numeric/check_tt_twojet.py',
            lambda mode,out:['--repo-root',self.repo,'--evaluator',self.evaluator,
                '--evaluator-sha256',self.evaluator_sha,'--input-module',self.input,
                '--input-module-sha256',self.input_sha,'--report',numeric[mode],
                '--report-sha256',sha(numeric[mode]),'--output',out],
            {mode:'review-numeric/REVIEW_'+mode.upper()+'.json' for mode,_ in MODES})
        for mode,optimized in MODES:
            for mutation,reason in (
                ('omit_yy_conversion','INDEPENDENT_PHYSICAL_CHAIN:yy:0'),
                ('wrong_negative_power_derivative','INDEPENDENT_FULL_SOURCE_LEIBNIZ:0'),
            ):
                out = self.output/('rejected-source-'+mutation+'-'+mode+'.json')
                self.run('source-mutant-'+mutation+'-'+mode,'review-inputs/oracle_review.py',
                    ['--repo-root',self.repo,'--inject-error',mutation,'--output',out],optimized,
                    reject_reason=reason,absent=out)
        # This control ledger has no paths or elapsed time and is itself a
        # canonical golden, so regenerate it exactly from the fresh failures.
        cases = []
        for row in self.controls:
            if row['script'] != 'review-inputs/oracle_review.py':
                continue
            mode = 'optimized' if row['optimized'] else 'normal'
            mutation = row['id'][len('source-mutant-'):-(len(mode)+1)]
            cases.append({'exit_code':row['exit_code'],
                          'flags':['-I','-S','-B']+(['-O'] if row['optimized'] else []),
                          'injection':mutation,'mode':mode,'output_created':False,
                          'reason':'ValueError: '+row['expected_rejection']})
        ledger = {'all_checks_active_under_optimization':True,'cases':cases,
                  'reviewer_source_sha256':sha(self.root/'review-inputs/oracle_review.py'),
                  'status':'PASS_INDEPENDENT_TT_FR_WRONG_FORMULA_CONTROLS'}
        destination = self.output/'source-math-control-ledger.json'
        destination.write_text(json.dumps(ledger,indent=2,sort_keys=True)+'\n')
        self.golden(destination,'review-inputs/CONTROL_LEDGER.json','source-math-control-ledger')

    def provenance(self):
        # Symbolic derivation is provenance, not an undeclared runtime SymPy
        # dependency.  Its two independently retained modes must still agree.
        left = self.root/'review-inputs/DEPENDENCY_NORMAL.json'
        right = self.root/'review-inputs/DEPENDENCY_OPTIMIZED.json'
        need(left.read_bytes() == right.read_bytes(), 'SYMBOLIC_PROVENANCE_MODE_PARITY')
        return {'symbolic_dependency_reports_byte_identical':True,'dependency_report_sha256':sha(left),
                'exact_fixture_sha256':sha(self.root/'review-inputs/EXACT_FIXTURE.json'),
                'runtime_sympy_imported':False,'provenance_generators_executed':False,
                'fresh_independent_checks':['source-dependency','theory-controls','independent-source']}

    def program_index(self):
        direct = {}
        for run in self.runs:
            direct.setdefault(run['script'],[]).append(run['id'])
        rows = []
        for path in sorted(self.root.rglob('*.py')):
            relative = path.relative_to(self.root).as_posix()
            row = {'path':relative,'sha256':sha(path)}
            if relative == 'reproduce.py':
                row['coverage'] = 'Executing standard-library entry point: complete capsule/nine-upstream authentication, optional ZIP checks, seven reason-specific source-corruption controls, isolated replay and exact golden comparisons.'
            elif relative in direct:
                row.update(direct_run_ids=sorted(direct[relative]),coverage='Executed in normal and optimized isolated standard-library Python; declared negative controls are separately recorded.')
                if relative == 'inputs/actual_third_jets.py':
                    row['nested_module_imports_by'] = ['compute/tt_twojet.py','review-inputs/oracle_review.py','review-numeric/check_tt_twojet.py']
                elif relative == 'compute/tt_twojet.py':
                    row['nested_module_imports_by'] = ['compute/check_polynomial_controls.py','review-numeric/check_tt_twojet.py:isolated_synthetic_driver_guards_only']
            elif relative in PROVENANCE_PROGRAMS:
                row.update(coverage=PROVENANCE_PROGRAMS[relative],runtime_execution=False)
            else:
                raise ValueError('UNINDEXED_CAPSULE_PROGRAM:'+relative)
            rows.append(row)
        need(set(PROVENANCE_PROGRAMS).issubset({row['path'] for row in rows}), 'MISSING_PROVENANCE_PROGRAM')
        return rows


def output_index(output):
    return [{'path':path.relative_to(output).as_posix(),'bytes':path.stat().st_size,'sha256':sha(path)}
            for path in sorted(output.rglob('*')) if path.is_file()
            and path.name not in ('REPRODUCTION_RECEIPT.json','REPRODUCTION_FAILURE.json')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256',required=True)
    parser.add_argument('--repo-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--archive-sha256')
    args = parser.parse_args()
    root,repo = Path(__file__).resolve().parent,args.repo_root.resolve()
    output = guard_output(root,repo,args.output_dir)
    need(bool(args.archive) == bool(args.archive_sha256), 'ARCHIVE_AND_EXTERNAL_PIN_REQUIRED_TOGETHER')
    before = authenticate(root,args.manifest_sha256,new_capsule=True)
    upstream_before = authenticate_upstream(repo)
    archive = archive_check(root,args.archive,args.archive_sha256) if args.archive else None
    output.mkdir(parents=True)
    replay = Replay(root,repo,output)
    try:
        controls = corruption_controls(root,args.manifest_sha256,output)
        fresh = replay.produce()
        replay.reviews(fresh)
        provenance = replay.provenance()
        index = replay.program_index()
        after,upstream_after = authenticate(root,args.manifest_sha256,new_capsule=True),authenticate_upstream(repo)
        need(before == after and upstream_before == upstream_after, 'SOURCE_CHANGED_DURING_REPLAY')
        positive_runs = sum('expected_rejection' not in row for row in replay.runs)
        negative_runs = len(replay.controls)
        final_pins = {'capsule_manifest_sha256':args.manifest_sha256,
                      'input_module_sha256':replay.input_sha,'evaluator_sha256':replay.evaluator_sha,
                      'final_numeric_report_sha256':sha(fresh['normal']),
                      'independent_numeric_reviewer_sha256':sha(root/'review-numeric/check_tt_twojet.py')}
        receipt = {'schema':'HDBLAST_EXTERIOR_TENSOR_JETS_REPRODUCTION_V1','status':'PASS',
            'python_version':sys.version,'complete_new_capsule_before':before,'complete_new_capsule_after':after,
            'complete_nine_upstream_roots_before':upstream_before,'complete_nine_upstream_roots_after':upstream_after,
            'all_complete_authentication_before_any_project_execution':True,'archive':archive,
            'source_corruption_controls':controls,'direct_runs':sorted(replay.runs,key=lambda row:row['id']),
            'positive_direct_runs':positive_runs,'expected_nonzero_direct_runs':negative_runs,'final_pins':final_pins,
            'scientific_golden_matches_mandatory':True,'golden_matches':sorted(replay.matches,key=lambda row:(row['actor'],row['golden'])),
            'independent_math_negative_controls':sorted(replay.controls,key=lambda row:row['id']),
            'complete_program_index':index,'symbolic_provenance':provenance,'output_file_index':output_index(output),
            'source_unchanged_after_replay':True,
            'log_scope':'Scientific JSONs and canonical independent reviews are byte-compared. Environment paths, Python version, stdout logs and execution receipt metadata are indexed separately and are not golden scientific reports.',
            'replay_scope':'Only new bounded producers, exact controls and independent checks run. Frozen upstream functions are imported through authenticated closures; old full campaign CLIs and SymPy provenance generators are never rerun.',
            'limits':['Conditional registered-gauge ordinary coefficient two-jet with all actual source errors and complete scalar feedback retained.',
                      'The declared second-derivative sign diagnostic may remain inconclusive; finite enclosures are the primary result.',
                      'No curvature, gauge-invariant observable, nonlinear solution family, origin or novelty claim.']}
        with (output/'REPRODUCTION_RECEIPT.json').open('x') as handle:
            json.dump(receipt,handle,indent=2,sort_keys=True)
            handle.write('\n')
        print(json.dumps({'status':'PASS','direct_runs':len(replay.runs),
                          'positive_direct_runs':positive_runs,'expected_nonzero_direct_runs':negative_runs,
                          'scientific_golden_matches':len(replay.matches),'final_pins':final_pins,
                          'source_corruption_controls':len(controls),'independent_math_negative_controls':len(replay.controls),
                          'indexed_programs':len(index),'complete_upstream_closures':len(UPSTREAM)}),flush=True)
    except BaseException as error:
        failure = {'status':'FAIL','error':str(error),'completed_direct_runs':replay.runs,'completed_golden_matches':replay.matches}
        try:
            failure['capsule_unchanged_after_failure'] = authenticate(root,args.manifest_sha256,new_capsule=True) == before
            failure['upstream_unchanged_after_failure'] = authenticate_upstream(repo) == upstream_before
        except Exception as auth_error:
            failure['post_failure_authentication_error'] = str(auth_error)
        (output/'REPRODUCTION_FAILURE.json').write_text(json.dumps(failure,indent=2,sort_keys=True)+'\n')
        raise


if __name__ == '__main__':
    main()
