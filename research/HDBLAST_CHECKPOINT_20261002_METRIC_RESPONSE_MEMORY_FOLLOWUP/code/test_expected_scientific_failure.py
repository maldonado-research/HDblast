#!/usr/bin/env python3
"""Test the expected-failure classifier with fabricated metadata/NPZ only."""
from __future__ import annotations
import argparse
import copy
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import sys
import subprocess
import zipfile

import numpy as np

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('expected_science_failure',HERE/'reproduce_metric_expected_scientific_failure.py')
helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
BINDINGS={'public_freeze_commit':'a'*40,'registration_sha256':'b'*64,'independent_manifest_sha256':'c'*64}

def write(path,value): helper.write(path,value)

def fabricated_npz(path,run,schema,missing=None,wrong_dtype=None):
    """Stream synthetic zero bytes; never import or call a physical producer."""
    coarse=run['setting']=='coarse';nodes=8192 if coarse else 16384;points=577 if coarse else 1153
    zeros=b'\0'*(64*1024)
    with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for member in schema['members']:
            name=member['name']
            if name==missing: continue
            shape=tuple(nodes if x==16384 else points if x==1153 else x for x in member['shape'])
            dtype=np.dtype('<f8' if name==wrong_dtype else member['dtype'])
            header=io.BytesIO();np.lib.format.write_array_header_1_0(header,{'descr':dtype.str,'fortran_order':False,'shape':shape})
            size=math.prod(shape)*dtype.itemsize
            with archive.open(name+'.npy','w',force_zip64=True) as stream:
                stream.write(header.getvalue())
                while size:
                    count=min(size,len(zeros));stream.write(zeros[:count]);size-=count

def make_fixture(root,schema,frozen_experiment):
    root.mkdir();fresh=root/'fresh';(fresh/'primary').mkdir(parents=True);(fresh/'independent').mkdir()
    experiment=copy.deepcopy(frozen_experiment);experiment['FABRICATED']=True
    config=experiment['independent_configuration'];gates=experiment['gates']['independent']
    helper.require(experiment['primary']['producer_wall_seconds']==900,'Actual frozen primary budget key/value required')
    primary={'status':'completed','rows':[{'source':s,'eta':e} for s in helper.SOURCES for e in helper.OBSERVATIONS],
             'elapsed_seconds':1.0,'provenance':{k:BINDINGS[k] for k in ('public_freeze_commit','registration_sha256')},'FABRICATED':True}
    write(fresh/'primary/results.json',primary)
    names=['fabricated_prerequisite_'+str(n) for n in range(26)]+['primary','independent_modes']+['validation','validation_optimized','summary','figures']
    plan=[{'check':name,'physical_producer':name in ('primary','independent_modes'),'arguments':['/FABRICATED/'+name],
           'command':[sys.executable,'/FABRICATED/'+name],'timeout_seconds':1000,'cwd':'/FABRICATED/checkpoint'} for name in names]
    commands=[{**c,'status':'PASS','exit_code':0,'timed_out':False} for c in plan[:28]]
    commands[-1].update(status='FAIL',exit_code=1)
    ledger={'status':'FAIL','source_verification_before':'PASS','source_verification_after':'PASS','planned_commands':32,
            'attempted_commands':28,'successful_commands':27,'commands':commands,
            'physical_producer_commands_attempted':2,'physical_producer_commands_completed':1,
            'failure':{'exception':"RuntimeError('independent_modes failed; preserved stdout/stderr/exit and fresh outputs')"},'FABRICATED':True}
    write(root/'VALIDATION.json',ledger);write(root/'COMMAND_PLAN.json',plan);write(root/'INPUT_MANIFEST.json',{**BINDINGS,'FABRICATED':True})
    runs=[];combined=[]
    for source,setting in helper.RUNS:
        coarse=setting=='coarse';steps=576 if coarse else 1152;nodes=8192 if coarse else 16384
        rows=[]
        for eta in helper.OBSERVATIONS:
            points=[]
            for cutoff in helper.CUTOFFS:
                witness=(source,setting,eta,cutoff)==('positive_B','fine',-1.5,256)
                ledger_value=1e-5 if witness else 0.0
                points.append({'K':cutoff,'values':{k:0.0 for k in helper.QUANTITIES},
                               'ward':{'ledger':ledger_value,'direct_endpoint':0.0,'endpoint_difference':abs(ledger_value)}})
            rows.append({'source':source,'eta':eta,'finite_k':points})
        run={'source':source,'setting':setting,'time_steps':steps,'momentum_nodes':nodes,
             'direct_stress_history_points':steps+1,'dt':1/(128 if coarse else 256),'momentum_panel_width':0.5 if coarse else 0.25,
             'wronskian_max_scaled':0.0,'physical_wronskian_max_scaled':0.0,'rows':rows}
        path=fresh/'independent'/f'metric_modes_{source}_{setting}.npz';fabricated_npz(path,run,schema)
        run['archive']={'path':path.name,'sha256':helper.sha(path)};runs.append(run)
        if setting=='fine': combined.extend(copy.deepcopy(rows))
    failures=['positive_B/-1.5/256 Ward endpoint','positive_B/-1.5/256 Ward refinement']
    modes={'status':'failed','active_run':None,'freeze_commit':BINDINGS['public_freeze_commit'],
           'provenance':{'registration_sha256':BINDINGS['registration_sha256'],'manifest_sha256':BINDINGS['independent_manifest_sha256']},
           'configuration':config,'resources':{'elapsed_seconds':2.0,'peak_rss_kib':100000},'runs':runs,'rows':combined,
           'gate_failures':failures,'exception':{'type':'RuntimeError','message':'Registered internal gates failed: '+', '.join(failures)},
           'maximum_refinement_differences':{k:0.0 for k in helper.QUANTITIES},
           'maximum_ward_endpoint_difference':1e-5,'maximum_ward_refinement_difference':1e-5,'FABRICATED':True}
    write(fresh/'independent/results.json',modes)
    ledger['fresh_output_sha256']={p.relative_to(fresh).as_posix():helper.sha(p) for p in sorted(fresh.rglob('*')) if p.is_file()}
    write(root/'VALIDATION.json',ledger)
    return experiment

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--schema',type=Path,required=True);parser.add_argument('--experiment',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=args.output.resolve();helper.require(not output.exists(),'Fresh synthetic output required');output.mkdir(parents=True)
    schema=helper.read(args.schema);fixture=output/'FABRICATED_REPLAY';experiment=make_fixture(fixture,schema,helper.read(args.experiment))
    passed=[]
    def classify(exit_code=1,timed_out=False): return helper.qualify_saved_failure(fixture,exit_code,BINDINGS,experiment,schema,timed_out)
    result=classify();helper.require(result['control_status']=='EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED' and result['underlying_scientific_status']=='FAIL' and result['registered_metric_calibration_passed'] is False,'Fabricated positive classification failed');passed.append('fabricated_complete_known_scientific_failure')
    def reject(name,action):
        try: action()
        except (RuntimeError,KeyError,FileNotFoundError,zipfile.BadZipFile,ValueError): passed.append(name)
        else: raise RuntimeError('Mutation accepted: '+name)
    def mutate_json(name,path,mutator):
        original=path.read_bytes();value=helper.read(path);mutator(value);write(path,value)
        receipt=fixture/'VALIDATION.json';receipt_original=receipt.read_bytes()
        if path.is_relative_to(fixture/'fresh'):
            d=helper.read(receipt);d['fresh_output_sha256'][path.relative_to(fixture/'fresh').as_posix()]=helper.sha(path);write(receipt,d)
        try: reject(name,classify)
        finally:
            path.write_bytes(original)
            if path.is_relative_to(fixture/'fresh'): receipt.write_bytes(receipt_original)
    validation=fixture/'VALIDATION.json';modes=fixture/'fresh/independent/results.json';primary=fixture/'fresh/primary/results.json'
    reject('unexpected_success_exit',lambda:classify(0));reject('signal_exit',lambda:classify(-9));reject('timeout',lambda:classify(1,True))
    mutate_json('replay_unexpected_pass',validation,lambda d:d.update(status='PASS'))
    mutate_json('source_changed_after',validation,lambda d:d.update(source_verification_after='FAIL'))
    mutate_json('wrong_command_prefix',validation,lambda d:d['commands'][0].update(check='wrong'))
    mutate_json('wrong_command_arguments',validation,lambda d:d['commands'][0].update(arguments=['/wrong']))
    mutate_json('earlier_prerequisite_failure',validation,lambda d:d['commands'][0].update(status='FAIL',exit_code=1))
    mutate_json('independent_signal_exit',validation,lambda d:d['commands'][-1].update(exit_code=-9))
    mutate_json('independent_timeout',validation,lambda d:d['commands'][-1].update(timed_out=True))
    mutate_json('wrong_failed_command',validation,lambda d:d['commands'][-1].update(check='validation'))
    mutate_json('primary_incomplete',primary,lambda d:d['rows'].pop())
    mutate_json('primary_wall_budget',primary,lambda d:d.update(elapsed_seconds=901))
    value=experiment['primary'].pop('producer_wall_seconds')
    try: reject('missing_actual_primary_producer_wall_seconds_key',classify)
    finally: experiment['primary']['producer_wall_seconds']=value
    mutate_json('primary_provenance_changed',primary,lambda d:d['provenance'].update(registration_sha256='d'*64))
    mutate_json('missing_fourth_run',modes,lambda d:d['runs'].pop())
    mutate_json('duplicate_run_identity',modes,lambda d:d['runs'][-1].update(source='positive_B'))
    mutate_json('wrong_run_grid',modes,lambda d:d['runs'][0].update(momentum_nodes=1))
    mutate_json('missing_run_observation',modes,lambda d:d['runs'][0]['rows'].pop())
    mutate_json('rss_budget_exceeded',modes,lambda d:d['resources'].update(peak_rss_kib=262145))
    mutate_json('wall_budget_exceeded',modes,lambda d:d['resources'].update(elapsed_seconds=901))
    mutate_json('nonfinite_resources',modes,lambda d:d['resources'].update(peak_rss_kib=-1))
    mutate_json('other_exception',modes,lambda d:d['exception'].update(message='Frozen 256-MiB peak-memory budget exceeded'))
    mutate_json('missing_known_gate_failure',modes,lambda d:d['gate_failures'].pop(0))
    mutate_json('fabricated_maximum',modes,lambda d:d.update(maximum_ward_endpoint_difference=1.0))
    mutate_json('wronskian_failure',modes,lambda d:d['runs'][0].update(wronskian_max_scaled=1.0))
    mutate_json('saved_archive_hash_changed',modes,lambda d:d['runs'][0]['archive'].update(sha256='e'*64))
    mutate_json('combined_values_disagree',modes,lambda d:d['rows'][0]['finite_k'][0]['values'].update(q=1.0))
    mutate_json('combined_Ward_disagrees',modes,lambda d:d['rows'][-1]['finite_k'][-1]['ward'].update(ledger=1.0))
    mutate_json('fresh_output_hash_changed',validation,lambda d:d['fresh_output_sha256'].update({'primary/results.json':'0'*64}))
    (fixture/'fresh/report').mkdir()
    try: reject('unexpected_post_failure_report',classify)
    finally: (fixture/'fresh/report').rmdir()
    archive=fixture/'fresh/independent/metric_modes_positive_B_coarse.npz';original=archive.read_bytes();saved_modes=modes.read_bytes();saved_ledger=validation.read_bytes();run=helper.read(modes)['runs'][0]
    def rebind_outputs():
        d=helper.read(validation);d['fresh_output_sha256']={p.relative_to(fixture/'fresh').as_posix():helper.sha(p) for p in sorted((fixture/'fresh').rglob('*')) if p.is_file()};write(validation,d)
    try:
        fabricated_npz(archive,run,schema,missing='sub_deltaW2_0');d=helper.read(modes);d['runs'][0]['archive']['sha256']=helper.sha(archive);write(modes,d)
        rebind_outputs()
        reject('missing_full_contact_archive_member',classify)
        fabricated_npz(archive,run,schema,wrong_dtype='k');d['runs'][0]['archive']['sha256']=helper.sha(archive);write(modes,d)
        rebind_outputs()
        reject('wrong_archive_dtype',classify)
    finally: archive.write_bytes(original);modes.write_bytes(saved_modes);validation.write_bytes(saved_ledger)
    archive.rename(archive.with_suffix('.hidden'))
    try: reject('missing_archive',classify)
    finally: archive.with_suffix('.hidden').rename(archive)
    auth=output/helper.CHECKPOINT;(auth/'code').mkdir(parents=True);(auth/'independent').mkdir()
    (auth/'code/metric_primary.py').write_text('# FABRICATED placeholder; never executed\n')
    (auth/'independent/forced_metric.py').write_text('# FABRICATED placeholder; never executed\n')
    write(auth/'independent/PRESERVED_ARCHIVE_SCHEMA.json',schema)
    write(auth/'independent/MANIFEST.json',{'FABRICATED':True,'files':[]})
    driver_source='''# FABRICATED pure authentication stub, no producers/subprocesses.\nimport json\ndef verify_inputs(root,override=None):\n receipt=json.loads((root/'FREEZE_RECEIPT.json').read_text())\n registration=json.loads((root/'FULL_REGISTRATION.json').read_text())\n return dict(receipt,frozen_files=registration['files'])\n'''
    (auth/'code/replay_metric.py').write_text(driver_source)
    files={p.relative_to(auth).as_posix():helper.sha(p) for p in sorted(auth.rglob('*')) if p.is_file()}
    write(auth/'FULL_REGISTRATION.json',{'files':files,'FABRICATED':True})
    auth_bindings={'public_freeze_commit':'f'*40,'registration_sha256':helper.sha(auth/'FULL_REGISTRATION.json'),
                   'independent_manifest_sha256':helper.sha(auth/'independent/MANIFEST.json')}
    write(auth/'FREEZE_RECEIPT.json',auth_bindings)
    driver,inputs=helper.load_authenticated_driver(auth,auth_bindings)
    helper.require(all(inputs[k]==v for k,v in auth_bindings.items()),'Fabricated pure authentication failed');passed.append('explicit_public_bindings_authenticate_before_pure_import')
    for key,length in [('public_freeze_commit',40),('registration_sha256',64),('independent_manifest_sha256',64)]:
        wrong=dict(auth_bindings);wrong[key]='0'*length
        reject('explicit_binding_mismatch_'+key,lambda wrong=wrong:helper.load_authenticated_driver(auth,wrong))
    receipt=(auth/'FREEZE_RECEIPT.json').read_bytes()
    write(auth/'FREEZE_RECEIPT.json',{**auth_bindings,'public_freeze_commit':'1'*40})
    try: reject('receipt_public_freeze_mismatch',lambda:helper.load_authenticated_driver(auth,auth_bindings))
    finally: (auth/'FREEZE_RECEIPT.json').write_bytes(receipt)
    for name in ('code/replay_metric.py','independent/PRESERVED_ARCHIVE_SCHEMA.json'):
        path=auth/name;original=path.read_bytes();path.write_bytes(original+b'\n# MUTATED\n')
        try: reject('registered_bytes_changed_'+Path(name).name,lambda:helper.load_authenticated_driver(auth,auth_bindings))
        finally: path.write_bytes(original)
    path=auth/'code/metric_primary.py';path.rename(path.with_suffix('.hidden'));path.symlink_to(path.with_suffix('.hidden'))
    try: reject('registered_symlink',lambda:helper.load_authenticated_driver(auth,auth_bindings))
    finally: path.unlink();path.with_suffix('.hidden').rename(path)
    prior=Path('/workspace/hdblast-research-work/metric-followup-complete-replay-001')
    if prior.is_dir():
        prior_receipt=helper.read(prior/'checkpoint/FREEZE_RECEIPT.json')
        reject('saved_prior_memory_failure_is_not_expected_scientific_failure',lambda:helper.qualify_saved_failure(prior,1,prior_receipt,helper.read(prior/'checkpoint/EXPERIMENT.json'),schema))
    real_popen=helper.subprocess.Popen;real_killpg=helper.os.killpg;events=[]
    class FabricatedProcess:
        pid=424242;returncode=-9
        def __init__(self,*a,**kw):
            helper.require(kw['start_new_session'] is True,'Replay must start a new process group');self.calls=0;events.append('new_session')
        def __enter__(self): return self
        def __exit__(self,*a): return False
        def wait(self,timeout=None):
            self.calls+=1
            if self.calls==1:
                helper.require(timeout==2400,'Replay outer deadline differs');raise subprocess.TimeoutExpired('FABRICATED_NO_SUBPROCESS',timeout)
            events.append('wait_after_group_kill');return -9
    try:
        helper.subprocess.Popen=FabricatedProcess
        helper.os.killpg=lambda pid,sig:events.append(('kill_group',pid,sig))
        try: helper.run_with_process_group(['FABRICATED_NO_SUBPROCESS'],output,{},None,None)
        except subprocess.TimeoutExpired: pass
        else: raise RuntimeError('Fabricated timeout was swallowed')
        helper.require(events==['new_session',('kill_group',424242,helper.signal.SIGKILL),'wait_after_group_kill'],'Timeout must kill and reap complete process group');passed.append('fabricated_timeout_kills_and_reaps_process_group')
    finally: helper.subprocess.Popen=real_popen;helper.os.killpg=real_killpg
    name='code/reproduce_metric_expected_scientific_failure.py';pin=helper.sha(HERE/'reproduce_metric_expected_scientific_failure.py')
    helper.verify_registered_helper({'frozen_files':{name:pin}});passed.append('running_helper_bound_to_registration')
    reject('running_helper_omitted_from_registration',lambda:helper.verify_registered_helper({'frozen_files':{}}))
    reject('running_helper_hash_differs_from_registration',lambda:helper.verify_registered_helper({'frozen_files':{name:'0'*64}}))
    report={'status':'PASS','python_optimization':sys.flags.optimize,'checks':len(passed),'passed_checks':passed,
            'physical_evaluations':0,'replay_subprocesses':0,'fixtures':'Explicitly FABRICATED metadata and streamed zero-valued NPZ arrays only',
            'helper_sha256':helper.sha(HERE/'reproduce_metric_expected_scientific_failure.py'),
            'schema_sha256':helper.sha(args.schema),'production_control_executed':False}
    report['actual_frozen_experiment_sha256']=helper.sha(args.experiment)
    report['actual_primary_wall_key']='producer_wall_seconds'
    write(output/'CHECKS.json',report);print(json.dumps(report,indent=2))

if __name__=='__main__': main()
