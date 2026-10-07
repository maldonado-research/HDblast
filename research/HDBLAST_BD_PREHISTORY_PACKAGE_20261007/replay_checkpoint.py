#!/usr/bin/env python3
"""Portable additive replay; source authorization remains in frozen custodian.

Install this helper beside, never inside, the frozen BD prehistory checkpoint.
All expectations originate in externally supplied pin hashes. Verification
alone evaluates no source; a replay calls the already frozen bounded entry.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import importlib.util
import json
import os
import re
import stat
import subprocess
import sys
import time

sys.dont_write_bytecode=True
LIMIT=20*1024*1024

class ReplayStepFailure(RuntimeError):
    def __init__(self,label,code,log):
        super().__init__('Replay step failed '+label+'; inspect '+str(log))
        self.exit_code=code if code>=0 else 128-code

def require(value,label):
    if not value:
        raise ValueError(label)

def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()

def sha_file(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def load(raw):
    def unique(items):
        result={}
        for key,value in items:
            require(key not in result,'Duplicate JSON key')
            result[key]=value
        return result
    def invalid(value):
        raise ValueError('Nonfinite JSON constant '+value)
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)

def actual_directory(path):
    path=Path(os.path.abspath(path))
    require(path.is_dir() and path.resolve()==path,'Real directory required')
    require(all(not p.is_symlink() for p in (path,*path.parents)),'Directory symlink ancestor forbidden')
    return path

def safe_file(root,name):
    require(type(name) is str and name and '\\' not in name,'Canonical relative filename required')
    parsed=PurePosixPath(name)
    require(not parsed.is_absolute() and parsed.as_posix()==name and
            all(p not in ('','.','..') for p in name.split('/')),'Unsafe relative filename')
    path=root
    for part in parsed.parts:
        path=path/part
        require(not path.is_symlink(),'Payload symlink forbidden')
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode),'Regular payload file required '+name)
    return path

def pinned_json(path,expected):
    require(type(expected) is str and re.fullmatch('[0-9a-f]{64}',expected),'External SHA256 expectation required')
    path=Path(os.path.abspath(path))
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<=LIMIT and
            all(not p.is_symlink() for p in path.parents),'Bounded real JSON input required')
    raw=path.read_bytes()
    require(sha_bytes(raw)==expected,'Externally pinned JSON bytes differ')
    return load(raw),raw

def verify_payload(root,manifest):
    require(type(manifest) is dict and set(manifest)=={'schema_version','files'} and
            type(manifest['schema_version']) is int and manifest['schema_version']==1,'Exact payload manifest schema required')
    files=manifest['files']
    require(type(files) is dict and 'FULL_REGISTRATION.json' in files,'Complete checkpoint payload inventory required')
    for name,item in files.items():
        require(type(item) is dict and set(item)=={'bytes','sha256'} and
                type(item['bytes']) is int and item['bytes']>=0 and
                type(item['sha256']) is str and re.fullmatch('[0-9a-f]{64}',item['sha256']),
                'Invalid checkpoint payload identity')
        path=safe_file(root,name)
        require(path.stat().st_size==item['bytes'] and sha_file(path)==item['sha256'],
                'Checkpoint payload size/SHA256 differs '+name)
    actual=set()
    for path in root.rglob('*'):
        require(not path.is_symlink(),'Checkpoint symlink forbidden')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
    require(actual==set(files),'Payload inventory omits or adds checkpoint leaves')
    return len(files)

def prepare(root,pins_path,pins_sha,manifest_path,manifest_sha,go_path):
    root=actual_directory(root)
    pins,pins_raw=pinned_json(pins_path,pins_sha)
    required={'schema_version','payload_manifest_sha256','registration_sha256','public_go_sha256',
              'reference_data','reference_continuous_certificate'}
    require(type(pins) is dict and set(pins)==required and type(pins['schema_version']) is int and
            pins['schema_version']==1,'Exact replay pin schema required')
    for key in ('payload_manifest_sha256','registration_sha256','public_go_sha256'):
        require(type(pins[key]) is str and re.fullmatch('[0-9a-f]{64}',pins[key]),'Invalid replay pin '+key)
    require(pins['payload_manifest_sha256']==manifest_sha,'External and replay manifest pins differ')
    manifest,manifest_raw=pinned_json(manifest_path,manifest_sha)
    count=verify_payload(root,manifest)
    require(sha_file(safe_file(root,'FULL_REGISTRATION.json'))==pins['registration_sha256'],
            'Pinned complete registration differs')
    require(type(pins['reference_data']) is dict and set(pins['reference_data'])=={'primary','independent'},
            'Both original scientific outputs must be pinned')
    for relative in [*pins['reference_data'].values(),pins['reference_continuous_certificate']]:
        require(relative in manifest['files'],'Scientific reference omitted from complete payload')
        safe_file(root,relative)
    # Only after every payload leaf has been authenticated may frozen guard
    # code be loaded. Its runtime version and complete public freeze checks
    # occur before numerical route imports in every bounded subprocess too.
    execution=root/'execution'
    sys.path.insert(0,str(execution))
    spec=importlib.util.spec_from_file_location('portable_frozen_registration_guard',safe_file(root,'execution/registration_guard.py'))
    guard=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(guard)
    verified=guard.authenticate(root,go_path,pins['public_go_sha256'],pins['registration_sha256'])
    runtime=guard.runtime_identity()
    saved_review=verify_saved_science(root,pins,verified)
    return {'root':root,'pins':pins,'manifest':manifest,'payload_files':count,
            'guard':guard,'verified':verified,'runtime':runtime,
            'saved_scientific_review':saved_review,
            'pins_sha256':sha_bytes(pins_raw),'manifest_sha256':sha_bytes(manifest_raw)}

def verified_module(root,relative,label):
    # Caller has already authenticated every complete payload and registered
    # source byte before using this import helper.
    spec=importlib.util.spec_from_file_location(label,safe_file(root,relative))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def verify_saved_science(root,pins,verified):
    validate=verified_module(root,'execution/validate_outputs.py','portable_frozen_output_validator')
    post=verified_module(root,'analysis/certify_continuous_model.py','portable_frozen_continuous_postprocessor')
    collected={}
    for route in ('primary','independent'):
        data_path=safe_file(root,pins['reference_data'][route])
        parent=data_path.parent
        budget_path=safe_file(root,(parent/'BUDGET.json').relative_to(root).as_posix())
        entry_path=safe_file(root,(parent/'ENTRY_RECEIPT.json').relative_to(root).as_posix())
        execution=load(safe_file(root,(parent/'EXECUTION.json').relative_to(root).as_posix()).read_bytes())
        entry=load(entry_path.read_bytes())
        require(execution.get('status')=='PASS_BOUNDED_REGISTERED_ENTRY' and execution.get('exit_code')==0 and
                execution.get('fabricated_only') is False and execution.get('registration_sha256')==pins['registration_sha256'],
                'Original authoritative execution receipt differs')
        require(entry.get('status')=='PASS_AUTHENTICATED_OFFLINE_REGISTERED_ENTRY' and
                entry.get('route')==route and entry.get('fabricated_only') is False and
                entry.get('registration_sha256')==pins['registration_sha256'] and
                entry.get('receipt_sha256')==pins['public_go_sha256'] and
                entry.get('freeze_commit')==verified['receipt']['freeze_commit'],'Original registered entry differs')
        for name,path in [('DATA.json',data_path),('BUDGET.json',budget_path)]:
            require(entry.get('files',{}).get(name)=={'bytes':path.stat().st_size,'sha256':sha_file(path)},
                    'Original entry does not bind scientific bytes')
        callbacks=entry.get('source_callbacks',{})
        require(callbacks.get('physical_source_callbacks')==22 and callbacks.get('archive_arrays_decoded')==0,
                'Original callback inventory differs')
        events=callbacks.get('events')
        journal=safe_file(root,(parent/'SOURCE_ATTEMPTS.jsonl').relative_to(root).as_posix())
        require(type(events) is list and len(events)==22 and
                events==[load(line) for line in journal.read_bytes().splitlines()],
                'Original source journal does not match entry events')
        collected[route]=(load(data_path.read_bytes()),load(budget_path.read_bytes()),entry_path,budget_path)
    intersections=validate.compare_outputs(collected['primary'][0],collected['primary'][1],
                                          collected['independent'][0],collected['independent'][1],False,True)
    certificate=load(safe_file(root,pins['reference_continuous_certificate']).read_bytes())
    provenance=certificate.get('provenance',{})
    require(provenance.get('independent_budget_sha256')==sha_file(collected['independent'][3]) and
            provenance.get('independent_entry_receipt_sha256')==sha_file(collected['independent'][2]) and
            provenance.get('registration_sha256')==pins['registration_sha256'] and
            provenance.get('freeze_commit')==verified['receipt']['freeze_commit'] and
            provenance.get('postprocessor_sha256')==sha_file(root/'analysis/certify_continuous_model.py'),
            'Saved continuous certificate provenance differs')
    derived=post.build_certificate(collected['independent'][1])
    core=dict(certificate);core.pop('provenance',None)
    require(derived==core,'Saved continuous certificate mathematical core differs')
    return {'status':'PASS_SAVED_36_INTERSECTIONS_AND_CONTINUOUS_MODEL_CORE',
            'overlap_rectangles':intersections['overlap_rectangles'],'component_intersections':intersections['component_intersections'],
            'continuous_model_status':derived['status'],'new_physical_source_callbacks':0,'saved_array_decodes':0}

def python_command(python,optimized,script,arguments):
    return [str(python),'-B']+(['-O'] if optimized else [])+[str(script),*map(str,arguments)]

def run_step(label,cmd,out,steps,timeout=930):
    started=time.monotonic()
    log=out/(label+'.log')
    print(json.dumps({'step':label,'status':'RUNNING'}),flush=True)
    with log.open('xb') as handle:
        try:
            process=subprocess.run(cmd,stdout=handle,stderr=subprocess.STDOUT,timeout=timeout,
                                   env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','PYTHONHASHSEED':'0',
                                        'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1',
                                        'MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1'})
            code=process.returncode
        except subprocess.TimeoutExpired:
            code=124
    receipt={'step':label,'command':cmd,'exit_code':code,'wall_seconds':time.monotonic()-started,
             'log_sha256':sha_file(log)}
    steps.append(receipt)
    if code:
        raise ReplayStepFailure(label,code,log)
    print(json.dumps({'step':label,'status':'PASS'}),flush=True)

def compare_bytes(left,right,label):
    require(left.read_bytes()==right.read_bytes(),'Exact replay bytes differ '+label)
    return sha_file(left)

def fresh_output_directory(out,protected):
    out=Path(os.path.abspath(out))
    require(not out.exists() and all(out!=p and p not in out.parents for p in protected) and
            all(not p.is_symlink() for p in (out,*out.parents)),
            'Fresh output directory outside checkpoint and delivery package required')
    return out

def scientific_budget_core(budget,route):
    require(type(budget) is dict and route in ('primary','independent'),'Known route budget required')
    core=dict(budget)
    if route=='primary':
        resource={key:core.pop(key) for key in ('wall_seconds','peak_rss_kib')}
    else:
        require(type(core.get('resource')) is dict and set(core['resource'])=={'wall_seconds','peak_rss_kib'},
                'Only registered independent resource measurements may vary')
        resource=core['resource'];core['resource']={}
    require(type(resource['wall_seconds']) in (int,float) and resource['wall_seconds']>=0 and
            type(resource['peak_rss_kib']) is int and resource['peak_rss_kib']>=0,
            'Valid recorded resource measurements required')
    return core

def compare_budget_cores(left,right,route):
    left_core=scientific_budget_core(load(left.read_bytes()),route)
    right_core=scientific_budget_core(load(right.read_bytes()),route)
    canonical=lambda value:(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    left_raw=canonical(left_core);right_raw=canonical(right_core)
    require(left_raw==right_raw,'Exact canonical JSON scientific budget core differs '+route)
    return sha_bytes(left_raw)

def replay(prepared,go_path,python,out):
    root=prepared['root'];pins=prepared['pins']
    require(os.getuid()!=0,'Replay must run as existing nonroot user')
    require(sys.platform=='linux','Frozen worker requires Linux/libseccomp.so.2')
    python=Path(os.path.abspath(python))
    require(python.is_file() and os.access(python,os.X_OK),'Explicit executable Python required')
    out=fresh_output_directory(out,(root,Path(__file__).resolve().parent))
    out.mkdir(parents=True)
    steps=[]
    summary={'status':'RUNNING','scope':'REGISTERED_BD_PREHISTORY_REPLAY_ONLY',
             'payload_files':prepared['payload_files'],'replay_pins_sha256':prepared['pins_sha256'],
             'payload_manifest_sha256':prepared['manifest_sha256'],
             'registration_sha256':pins['registration_sha256'],'public_go_sha256':pins['public_go_sha256'],
             'helper_sha256':sha_file(Path(__file__)),'runtime':prepared['runtime'],
             'saved_scientific_review':prepared['saved_scientific_review'],
             'original_binary80_state_error':'NOT_ENCLOSED','full_twelve_case_certificate':'UNRESOLVED',
             'metric_calibration':'FAIL','higher_dimensional_Big_Bang_origin':'NOT_ESTABLISHED','steps':steps,
             'exact_scientific_budget_core_sha256':{}}
    try:
        for mode,optimized in [('normal',False),('optimized',True)]:
            for route in ('primary','independent'):
                arguments=['--root',root,'--receipt',Path(go_path).absolute(),
                    '--receipt-sha256',pins['public_go_sha256'],'--registration-sha256',pins['registration_sha256'],
                    '--python',python,'--route',route,'--output-directory',out/(route+'-'+mode)]
                if optimized: arguments.append('--optimized')
                run_step(route+'-'+mode,python_command(python,optimized,root/'execution/execute_bounded.py',arguments),out,steps)
                execution=load((out/(route+'-'+mode)/'EXECUTION.json').read_bytes())
                require(execution.get('status')=='PASS_BOUNDED_REGISTERED_ENTRY' and execution.get('exit_code')==0 and
                        execution.get('fabricated_only') is False and execution.get('optimized')==optimized,
                        'Authoritative bounded execution failed '+route+' '+mode)
                compare_bytes(out/(route+'-'+mode)/'DATA.json',safe_file(root,pins['reference_data'][route]),
                              'original '+route+' '+mode)
                reference_budget=safe_file(root,str(PurePosixPath(pins['reference_data'][route]).parent/'BUDGET.json'))
                core_sha=compare_budget_cores(out/(route+'-'+mode)/'BUDGET.json',reference_budget,route)
                summary['exact_scientific_budget_core_sha256'][route+'-'+mode]=core_sha
            run_step('intersections-'+mode,python_command(python,optimized,root/'execution/validate_outputs.py',[
                '--primary-data',out/('primary-'+mode)/'DATA.json','--primary-budget',out/('primary-'+mode)/'BUDGET.json',
                '--independent-data',out/('independent-'+mode)/'DATA.json','--independent-budget',out/('independent-'+mode)/'BUDGET.json',
                '--output',out/('INTERSECTIONS_'+mode+'.json')]),out,steps)
            intersections=load((out/('INTERSECTIONS_'+mode+'.json')).read_bytes())
            require(intersections.get('status')=='PASS_36_INDEPENDENT_BD_PREHISTORY_RECTANGLE_INTERSECTIONS' and
                    intersections.get('overlap_rectangles')==36 and intersections.get('component_intersections')==72,
                    'Complete independent intersection check failed')
            proof_tasks=[
                ('representation',root/'theory/verify_incoming_representation.py',[], 'evidence/preparation/PRIMARY_EXACT_'+mode.upper()+'.json'),
                ('independent-cap',root/'independent/verify_independent_ir_and_cap.py',['--output'], 'evidence/actual/INDEPENDENT_CAP_CURRENT_'+mode.upper()+'.json'),
                ('continuous-selftest',root/'analysis/certify_continuous_model.py',['--self-test','--output'], 'evidence/preparation/CONTINUOUS_MODEL_SYNTHETIC_'+mode.upper()+'.json'),
                ('primary-fabricated',root/'primary/test_fabricated.py',[], 'evidence/preparation/PRIMARY_FABRICATED_'+mode.upper()+'.json'),
                ('outer-controls',root/'execution/test_controls.py',['--contract',root/'REGISTRATION_CONTRACT.json','--output'], 'evidence/preparation/OUTER_CONTROLS_'+mode.upper()+'.json')]
            for label,script,arguments,reference in proof_tasks:
                dest=out/(label+'-'+mode+'.json')
                run_step(label+'-'+mode,python_command(python,optimized,script,[*arguments,dest]),out,steps)
                compare_bytes(dest,safe_file(root,reference),label+' '+mode)
            independent=out/('independent-'+mode)
            dest=out/('CONTINUOUS_MODEL_'+mode+'.json')
            run_step('continuous-model-'+mode,python_command(python,optimized,root/'analysis/certify_continuous_model.py',[
                '--budget',independent/'BUDGET.json','--budget-sha256',sha_file(independent/'BUDGET.json'),
                '--entry-receipt',independent/'ENTRY_RECEIPT.json','--entry-receipt-sha256',sha_file(independent/'ENTRY_RECEIPT.json'),
                '--output',dest]),out,steps)
            certificate=load(dest.read_bytes())
            reference=load(safe_file(root,pins['reference_continuous_certificate']).read_bytes())
            certificate.pop('provenance',None);reference.pop('provenance',None)
            require(certificate==reference,'Continuous mathematical certificate differs; execution provenance is retained separately')
        summary['exact_original_and_normal_optimized_data_sha256']={route:compare_bytes(
            out/(route+'-normal')/'DATA.json',out/(route+'-optimized')/'DATA.json',route+' normal/O')
            for route in ('primary','independent')}
        verify_payload(root,prepared['manifest'])
        prepared['guard'].authenticate(root,go_path,pins['public_go_sha256'],pins['registration_sha256'])
        summary['status']='PASS_COMPLETE_PORTABLE_BD_PREHISTORY_REPLAY'
        summary['source_callbacks_in_four_registered_replays']=88
        summary['saved_array_decodes']=0
    except BaseException as error:
        summary['status']='REPLAY_FAILED'
        summary['failure_type']=type(error).__name__
        summary['failure_message']=str(error)
        if isinstance(error,ReplayStepFailure):
            summary['underlying_exit_code']=error.exit_code
        raise
    finally:
        (out/'REPLAY_RECEIPT.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    return summary

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sibling=Path(__file__).resolve().parent
    parser.add_argument('--checkpoint',default=str(sibling.parent/'HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY'))
    parser.add_argument('--pins',default=str(sibling/'REPLAY_PINS.json'))
    parser.add_argument('--pins-sha256',required=True)
    parser.add_argument('--payload-manifest',default=str(sibling/'PAYLOAD_MANIFEST.json'))
    parser.add_argument('--payload-manifest-sha256',required=True)
    parser.add_argument('--public-go',default=str(sibling/'PUBLIC_GO.json'))
    parser.add_argument('--python',default=sys.executable)
    parser.add_argument('--output-directory')
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args()
    prepared=prepare(args.checkpoint,args.pins,args.pins_sha256,args.payload_manifest,args.payload_manifest_sha256,args.public_go)
    if args.verify_only:
        require(args.output_directory is None,'Verify-only creates no replay directory')
        print(json.dumps({'status':'PASS_COMPLETE_PAYLOAD_AND_PUBLIC_FREEZE_AUTHENTICATION',
                          'payload_files':prepared['payload_files'],'saved_scientific_review':prepared['saved_scientific_review'],
                          'new_physical_source_callbacks':0,'saved_array_decodes':0}))
        return
    require(args.output_directory is not None,'Fresh external output directory required for replay')
    result=replay(prepared,args.public_go,args.python,args.output_directory)
    print(json.dumps({key:result[key] for key in ('status','exact_original_and_normal_optimized_data_sha256',
                      'source_callbacks_in_four_registered_replays','saved_array_decodes')},sort_keys=True))

if __name__=='__main__':
    try:
        main()
    except ReplayStepFailure as error:
        print(str(error),file=sys.stderr)
        sys.exit(error.exit_code)
