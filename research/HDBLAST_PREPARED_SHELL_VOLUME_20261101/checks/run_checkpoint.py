#!/usr/bin/env python3
"""Authenticate the closed prepared-shell-volume capsule and replay only its new work."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,os
from pathlib import Path,PurePosixPath
import subprocess,sys
HERE=Path(__file__).resolve().parent
CAP=HERE.parent
class Failure(Exception):pass
def need(ok,witness):
    if not ok:raise Failure(witness)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)
def relative(x):
    need(isinstance(x,str) and bool(x),'SAFE_RELATIVE_PATH');p=PurePosixPath(x)
    need(not p.is_absolute() and '..' not in p.parts and p.as_posix()==x and x!='.','SAFE_RELATIVE_PATH');return Path(*p.parts)
def authenticate(root,pin=None,old=False):
    prefix='SOURCE' if old else 'ROOT'
    need(root.is_dir() and not root.is_symlink(),prefix+'_CAPSULE_ROOT')
    manifest=root/'MANIFEST.json';need(regular(manifest),prefix+'_MANIFEST_REGULAR');actual_pin=sha(manifest)
    need(pin is None or actual_pin==pin,prefix+'_MANIFEST_PIN')
    rows=read(manifest)['files']
    if isinstance(rows,dict):rows=[{'path':k,**v} for k,v in rows.items()]
    need(isinstance(rows,list),prefix+'_MANIFEST_ROWS');seen={'MANIFEST.json'}
    for row in rows:
        rel=relative(row['path']);need(str(rel) not in seen,prefix+'_DUPLICATE_MEMBER');seen.add(str(rel));p=root/rel
        need(regular(p) and type(row['bytes']) is int and row['bytes']>=0 and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],prefix+'_MEMBER_BYTES')
    actual=set()
    for p in root.rglob('*'):
        need(not p.is_symlink(),prefix+'_SYMLINK')
        if p.is_file():actual.add(str(p.relative_to(root)))
        else:need(p.is_dir(),prefix+'_SPECIAL_FILE')
    need(actual==seen,prefix+'_CLOSED_INVENTORY')
    return {'manifest_sha256':actual_pin,'members_including_manifest':len(seen)}
def sources(repo,plan):
    result={name:authenticate(repo/'research'/relative(name),pin,True) for name,pin in sorted(plan['source_capsules'].items())}
    for rel,row in sorted(plan['selected_source_files'].items()):
        p=repo/relative(rel);need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SELECTED_OLD_SOURCE')
    for rel,row in sorted(plan['critical_new_files'].items()):
        p=CAP/relative(rel);need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SELECTED_NEW_SOURCE')
    return result
def dig(value,path):
    for key in path.split('.'):
        value=value[int(key)] if isinstance(value,list) else value[key]
    return value
def match_fields(document,fields,witness):
    for key,value in fields.items():need(dig(document,key)==value,witness)
def control_contract(document,contract):
    match_fields(document,contract.get('fields',{}),'CONTROL_REQUIRED_FIELD')
    rows=dig(document,contract['rows'])
    if 'has_field'in contract:rows=[row for row in rows if contract['has_field']in row]
    if 'filter' in contract:rows=[row for row in rows if all(dig(row,k)==v for k,v in contract['filter'].items())]
    need(len(rows)==contract['count'],'CONTROL_COMPLETE_COUNT')
    cases={}
    for row in rows:
        case=dig(row,contract['case']);mode=contract['mode_constant'] if 'mode_constant' in contract else dig(row,contract['mode']);need(case in contract['expected_by_case'],'CONTROL_CASE_ROSTER')
        match_fields(row,contract.get('row_fields',{}),'CONTROL_ROW_REQUIRED_FIELD')
        match_fields(row,contract.get('row_fields_by_mode',{}).get(str(mode),{}),'CONTROL_MODE_REQUIRED_FIELD')
        match_fields(row,contract.get('row_fields_by_case',{}).get(case,{}),'CONTROL_CASE_REQUIRED_FIELD')
        for left,right in contract.get('row_equal_fields',[]):need(dig(row,left)==dig(row,right),'CONTROL_ARTIFACT_CONTRACT')
        expected=contract['expected_by_case'][case]
        need(dig(row,contract['observed_witness'])==expected,'EXACT_FINAL_CONTROL_WITNESS')
        if 'returncode' in contract:need(type(dig(row,contract['returncode'])) is int and dig(row,contract['returncode'])==contract.get('expected_returncode',2),'CONTROL_EXACT_RETURN_CODE')
        cases.setdefault(case,[]).append(mode)
    need(set(cases)==set(contract['expected_by_case']) and all(len(v)==len(contract['modes']) and set(v)==set(contract['modes']) for v in cases.values()),'CONTROL_BOTH_MODES_AND_CASES')
    positives={}
    for path,count in contract.get('positive_counts',{}).items():
        v=dig(document,path);need((len(v) if isinstance(v,(list,dict)) else v)==count,'CONTROL_POSITIVE_COUNT');positives[path]=count
    return {'actual_rejected_subprocesses':len(rows),'distinct_rejection_cases':len(cases),'separately_counted_positive_checks':positives}
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path)
    ap.add_argument('--manifest-sha256');ap.add_argument('--jobs',type=int,choices=(1,2,3,4,5,6),default=3);ap.add_argument('--authenticate-only',action='store_true');a=ap.parse_args()
    repo=a.repo_root.absolute();out=a.output_dir.absolute();need(all(not p.is_symlink() for p in [repo,*repo.parents,out,*out.parents]),'PATH_SYMLINK')
    repo=repo.resolve(strict=True);out=out.resolve();need(not out.exists(),'FRESH_OUTPUT_DIRECTORY')
    need(all(base!=out and base not in out.parents for base in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    own=authenticate(CAP,a.manifest_sha256);plan=read(HERE/'REPLAY_PLAN.json')
    need(plan['schema']=='HDBLAST_PREPARED_SHELL_VOLUME_ROOT_REPLAY_V1','REPLAY_PLAN_SCHEMA');before=sources(repo,plan)
    if a.authenticate_only:
        need(authenticate(CAP,own['manifest_sha256'])==own and sources(repo,plan)==before,'AUTHENTICATION_POSTCHECK')
        print(json.dumps({'status':'PASS','witness':'CLOSED_SOURCE_AUTHENTICATION','old_capsules':len(before)},sort_keys=True));return
    # No third-party or project scientific module has been imported above.
    need(sys.version.split()[0]==plan['python_version'],'PINNED_PYTHON_VERSION')
    env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
    libraries=sorted(plan['requirements']);need(all(x.isidentifier() for x in libraries),'REQUIREMENTS_MODULE_NAMES')
    command=[sys.executable,'-I','-B','-c','import importlib,json;print(json.dumps({k:importlib.import_module(k).__version__ for k in '+repr(libraries)+'},sort_keys=True))']
    version=subprocess.run(command,capture_output=True,text=True,env=env,timeout=30)
    need(version.returncode==0 and not version.stderr,'PINNED_REQUIREMENTS_MISSING');versions=json.loads(version.stdout);need(versions==plan['requirements'],'PINNED_REQUIREMENTS_VERSION')
    out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir();bindings={'CAP':str(CAP),'REPO':str(repo),'OUT':str(out)}
    resolve=lambda x:x.format_map(bindings)
    def execute(group):
        reports={};modes={};controls={};inprocess={};oracle_calls=0;oracle_kinds={};positive_subprocesses=0
        for step in group['steps']:
            flags=['-I','-B'] if step.get('third_party') else ['-I','-S','-B']
            command=[sys.executable,*flags,*(['-O'] if step.get('optimized') else []),*[resolve(x) for x in step['argv']]]
            if step.get('create_output_parent',False):Path(resolve(step['output'])).parent.mkdir(parents=True,exist_ok=True)
            proc=subprocess.run(command,capture_output=True,env=env,cwd=repo,timeout=step.get('timeout_seconds',600))
            log=out/'logs'/(group['name']+'-'+step['name']);log.with_suffix('.stdout').write_bytes(proc.stdout);log.with_suffix('.stderr').write_bytes(proc.stderr)
            need(proc.returncode==0,'REPLAY_STEP_FAILED:'+group['name']+'/'+step['name'])
            if step.get('final_stdout'):
                lines=proc.stdout.decode().splitlines();need(len(lines)==1 and not proc.stderr,'STEP_EXACT_FINAL_STDOUT')
                final=json.loads(lines[-1]);match_fields(final,step['final_stdout'],'STEP_FINAL_WITNESS')
            kind=step['kind']
            if kind in ['canonical','controls']:
                # One scientific process may produce several canonical artifacts.
                # Verify every artifact without rerunning its actual initial oracle.
                for artifact in step.get('artifacts',[step]):
                    produced=Path(resolve(artifact['output']));need(regular(produced),'PRODUCED_REPORT_REGULAR');document=read(produced)
                    if artifact.get('reference'):
                        expected=CAP/relative(artifact['reference']);need(regular(expected) and sha(expected)==artifact['expected_sha256'] and produced.read_bytes()==expected.read_bytes(),'CANONICAL_BYTE_IDENTITY:'+group['name']+'/'+step['name'])
                    match_fields(document,artifact.get('fields',{}),'CANONICAL_REQUIRED_FIELD')
                    for path,count in artifact.get('inprocess_controls',{}).items():
                        value=dig(document,path);need((len(value) if isinstance(value,(list,dict)) else value)==count,'INPROCESS_CONTROL_ROSTER');inprocess[step['name']+':'+artifact.get('canonical_key','controls')+':'+path]=count
                    if kind=='canonical':
                        key=artifact['canonical_key'];pin=sha(produced);need(key not in reports or reports[key]==pin,'NORMAL_OPTIMIZED_BYTE_IDENTITY');reports[key]=pin;modes.setdefault(key,[]).append(bool(step.get('optimized')))
                    if 'contract' in artifact:
                        controls[step['name']]=control_contract(document,artifact['contract'])
            elif kind=='success':need(bool(step.get('final_stdout')),'SUCCESS_REQUIRES_EXACT_FINAL_WITNESS')
            else:raise Failure('UNKNOWN_STEP_KIND')
            if 'oracle_counter' in step:
                counter=step['oracle_counter'];actual=dig(read(Path(resolve(counter['report']))),counter['path'])
                need(type(actual)is int and actual>=0 and actual==step['new_actual_oracle_calls'],'OBSERVED_ORACLE_CALL_COUNT')
                oracle_calls+=actual
                category=counter['kind'];oracle_kinds[category]=oracle_kinds.get(category,0)+actual
            else:
                need(step.get('new_actual_oracle_calls',0)==0,'NONZERO_ORACLE_COUNT_NEEDS_OBSERVED_COUNTER')
            if 'positive_counter'in step:
                counter=step['positive_counter'];actual=dig(read(Path(resolve(counter['report']))),counter['path'])
                need(type(actual)is int and actual>=0 and actual==step['positive_subprocesses'],'OBSERVED_POSITIVE_SUBPROCESS_COUNT')
                positive_subprocesses+=actual
            else:need(step.get('positive_subprocesses',0)==0,'POSITIVE_COUNT_NEEDS_OBSERVED_COUNTER')
            print(json.dumps({'stage':'completed','group':group['name'],'step':step['name']},sort_keys=True),flush=True)
        need(all(len(v)==2 and set(v)=={False,True} for v in modes.values()),'EVERY_CANONICAL_REPORT_BOTH_MODES')
        return {'canonical':reports,'strict_subprocess_controls':controls,'separately_counted_inprocess_checks':inprocess,'new_actual_oracle_calls_executed':oracle_calls,'oracle_call_counts_by_kind':oracle_kinds,'positive_control_subprocesses_executed':positive_subprocesses}
    groups={g['name']:g for g in plan['groups']};need(len(groups)==len(plan['groups']),'UNIQUE_GROUP_NAMES');results={}
    while len(results)<len(groups):
        ready=[g for name,g in groups.items() if name not in results and set(g.get('depends_on',[]))<=set(results)]
        need(bool(ready),'GROUP_DEPENDENCY_CYCLE')
        with ThreadPoolExecutor(max_workers=a.jobs) as pool:
            fresh=list(pool.map(execute,ready))
        for group,result in zip(ready,fresh):results[group['name']]=result
    pairs=sum(len(v['canonical']) for v in results.values());rejections=sum(c['actual_rejected_subprocesses'] for v in results.values() for c in v['strict_subprocess_controls'].values())
    oracle=sum(v['new_actual_oracle_calls_executed'] for v in results.values());positives=sum(v['positive_control_subprocesses_executed'] for v in results.values())
    kinds={}
    for group in results.values():
        for kind,count in group['oracle_call_counts_by_kind'].items():kinds[kind]=kinds.get(kind,0)+count
    need(pairs==plan['expected_new_canonical_report_pairs'],'COMPLETE_CANONICAL_REPORT_ROSTER')
    need(rejections==plan['expected_rejected_subprocesses'],'COMPLETE_REJECTED_SUBPROCESS_ROSTER')
    need(oracle==plan['expected_new_actual_oracle_calls_executed'],'COMPLETE_NEW_ORACLE_CALL_ACCOUNTING')
    need(positives==plan['expected_positive_control_subprocesses'],'COMPLETE_POSITIVE_SUBPROCESS_ACCOUNTING')
    need(kinds==plan['expected_oracle_call_counts_by_kind'],'COMPLETE_TYPED_ORACLE_CALL_ACCOUNTING')
    need(authenticate(CAP,own['manifest_sha256'])==own and sources(repo,plan)==before,'ALL_SOURCES_UNCHANGED_AFTER_REPLAY')
    result={'schema':'HDBLAST_PREPARED_SHELL_VOLUME_ROOT_REPLAY_RECEIPT_V1','status':'PASS_COMPLETE_CONDITIONAL_PREPARED_SHELL_VOLUME_CHECKPOINT',
        'capsule_manifest_sha256':own['manifest_sha256'],'canonical_report_pairs_reproduced':pairs,'actual_rejected_subprocesses':rejections,
        'positive_control_subprocesses_executed':positives,'new_actual_oracle_calls_executed_including_positive_controls':oracle,
        'oracle_call_counts_by_kind':kinds,
        'primary_numeric_actual_oracle_calls_per_report_mode':plan['primary_numeric_actual_oracle_calls_per_report_mode'],
        'all_sources_authenticated_before_scientific_imports_and_after':True,'inherited_closed_capsules_authenticated':before,
        'inherited_checkpoint_proof_replay_executed':False,'inherited_source_authentication_is_proof_replay':False,
        'requirements':versions,'python_version':plan['python_version'],'groups':results,'scope':plan['scope']}
    receipt=out/'VERIFICATION.json';receipt.write_bytes(canonical(result))
    print(json.dumps({'status':'PASS','witness':'COMPLETE_NEW_PREPARED_SHELL_VOLUME_REPLAY','canonical_report_pairs':pairs,'actual_rejected_subprocesses':rejections,'verification_sha256':sha(receipt)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
