#!/usr/bin/env python3
"""Synthetic portable replay guards; never execute a physical child."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

STUBS=('derive_metric_contacts.py','verify_primary_preflight.py','metric_contact_algebra.py','verify_metric_ward_trace.py','verify_metric_response.py','verify_primary_inventory_agreement.py','verify_cse_specialization.py','verify_prerequisites.py','verify_metric_tail_algebra.py','test_metric_validator_guards.py','verify_raw_metric_baselines.py','verify_inherited_sources.py','certify_metric_derivative_bounds.py','metric_primary.py','validate_metric.py')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):path.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def fixture(root,driver):
    (root/'code').mkdir(parents=True);(root/'independent').mkdir();(root/'postrun').mkdir();(root/'theory').mkdir()
    shutil.copyfile(driver,root/'code/replay_metric.py')
    for name in STUBS:(root/'code'/name).write_text('raise RuntimeError("Synthetic guard forbids every child command")\n')
    for name in ('summarize_metric.py','plot_metric.py'):(root/'postrun'/name).write_text('raise RuntimeError("Synthetic guard forbids every child command")\n')
    (root/'independent/forced_metric.py').write_text('raise RuntimeError("Synthetic guard forbids every physical child")\n')
    write(root/'theory/FINITE_K_CONTACT_INVENTORY.json',{})
    experiment={'runtime':{'python_major_minor':[3,12],'physical_producers_optimized':False,'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'},'gates':{'synthetic':1}}
    write(root/'EXPERIMENT.json',experiment)
    write(root/'independent/MANIFEST.json',{'files':[{'path':'forced_metric.py','sha256':sha(root/'independent/forced_metric.py')}]})
    registration={'schema_version':1,'files':{p.relative_to(root).as_posix():sha(p) for p in sorted(root.rglob('*')) if p.is_file()},'independent_manifest_sha256':sha(root/'independent/MANIFEST.json'),'frozen_configuration':experiment,'frozen_gates':experiment['gates']}
    write(root/'FULL_REGISTRATION.json',registration)
    write(root/'FREEZE_RECEIPT.json',{'public_freeze_commit':'a'*40,'registration_sha256':sha(root/'FULL_REGISTRATION.json'),'independent_manifest_sha256':sha(root/'independent/MANIFEST.json')})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Fresh synthetic guard evidence required')
    args.output.mkdir(parents=True)
    driver=Path(__file__).with_name('replay_metric.py')
    reports=[]
    def case(name,mutate=None,flags=(),expected=1,output_selector=None):
        root=args.output/(name+'_source');fixture(root,driver)
        if mutate:mutate(root)
        destination=args.output/(name+'_result') if output_selector is None else output_selector(root)
        command=[sys.executable,str(root/'code/replay_metric.py'),'--plan-only','--output',str(destination),*flags]
        result=subprocess.run(command,capture_output=True,text=True,timeout=30)
        (args.output/(name+'.stdout.log')).write_text(result.stdout);(args.output/(name+'.stderr.log')).write_text(result.stderr)
        if (result.returncode==0)!=(expected==0):raise RuntimeError('Guard exit mismatch '+name+': '+result.stderr)
        state=json.loads((destination/'VALIDATION.json').read_text()) if (destination/'VALIDATION.json').exists() else None
        if state:
            if state['attempted_commands']!=0 or state['physical_producer_commands_attempted']!=0:raise RuntimeError('Guard launched a child '+name)
            if expected==0:
                plan=json.loads((destination/'COMMAND_PLAN.json').read_text())
                if len(plan)!=32 or sum(x['physical_producer'] for x in plan)!=2:raise RuntimeError('Unexpected command plan')
                if any(x['command'][1]=='-O' for x in plan if x['physical_producer']):raise RuntimeError('Optimized producer')
        reports.append({'name':name,'passed':True,'exit_code':result.returncode,'physical_children_executed':0,'attempted_commands':0})
    case('verified_plan',expected=0)
    case('changed_registration',lambda r:(r/'FULL_REGISTRATION.json').write_text('{}\n'))
    case('changed_frozen_source',lambda r:(r/'code/metric_primary.py').write_text('changed\n'))
    case('missing_freeze_receipt',lambda r:(r/'FREEZE_RECEIPT.json').unlink())
    case('wrong_external_commit',flags=('--freeze-commit','b'*40))
    case('changed_independent_manifest',lambda r:(r/'independent/MANIFEST.json').write_text('{}\n'))
    case('missing_registered_command',lambda r:(r/'code/verify_metric_response.py').unlink())
    case('symlink_registered_source',lambda r:((r/'code/metric_primary.py').unlink(),(r/'code/metric_primary.py').symlink_to(r/'code/validate_metric.py')))
    case('output_inside_checkpoint',output_selector=lambda r:r/'forbidden_output')
    def existing(root):
        path=args.output/'existing_result';path.mkdir();return path
    case('existing_output',output_selector=existing)
    def escape(root):
        reg=json.loads((root/'FULL_REGISTRATION.json').read_text());reg['files']['../escape.py']='0'*64;write(root/'FULL_REGISTRATION.json',reg)
        rec=json.loads((root/'FREEZE_RECEIPT.json').read_text());rec['registration_sha256']=sha(root/'FULL_REGISTRATION.json');write(root/'FREEZE_RECEIPT.json',rec)
    case('registered_path_escape',escape)
    def manifest_without_receipt(root):
        files={p.relative_to(root).as_posix():sha(p) for p in root.rglob('*') if p.is_file() and p.name!='FREEZE_RECEIPT.json'}
        write(root/'MANIFEST.json',{'files':files,'excluded':['MANIFEST.json','FREEZE_RECEIPT.json']})
    case('package_excludes_freeze_receipt',manifest_without_receipt)
    def valid_manifest(root):
        write(root/'MANIFEST.json',{'files':{p.relative_to(root).as_posix():sha(p) for p in root.rglob('*') if p.is_file()},'excluded':['MANIFEST.json']})
    case('verified_package_plan',valid_manifest,expected=0)
    def extra(root):valid_manifest(root);(root/'unregistered.txt').write_text('unexpected payload\n')
    case('extra_package_payload',extra)
    receipt={'status':'PASS','scope':'Synthetic source trees, authenticated plan-only commands and rejection guards. No physical child or source evaluation.','check_count':len(reports),'checks':reports,'python_optimization':sys.flags.optimize,'driver_sha256':sha(driver),'verifier_sha256':sha(__file__),'physical_children_executed':0}
    write(args.output/'GUARDS.json',receipt)
    print(json.dumps({'status':'PASS','check_count':len(reports),'physical_children_executed':0}))

if __name__=='__main__':main()
