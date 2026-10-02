#!/usr/bin/env python3
"""Replay the frozen memory followup and qualify an expected scientific FAIL.

This control passes only when the unchanged replay reaches completed primary
and all four independent runs, preserves all complete archives within budget,
and fails its registered internal science gates including the declared Ward
witness. Control success never changes the experiment's scientific FAIL.
Public freeze/registration/manifest bindings must be supplied explicitly.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback
import zipfile

CHECKPOINT = 'HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP'
VERSIONS = {'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}
SOURCES = ('positive_B','signed_uB')
OBSERVATIONS = (-5.5,-4.5,-4.0,-3.5,-2.5,-1.5)
CUTOFFS = (64,128,256)
QUANTITIES = ('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
RUNS = tuple((s,t) for s in SOURCES for t in ('coarse','fine'))

def require(ok,message):
    if not ok: raise RuntimeError(message)
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''): digest.update(block)
    return digest.hexdigest()
def read(path): return json.loads(Path(path).read_text())
def write(path,value): Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def stamp(): return datetime.now(timezone.utc).isoformat()
def finite(value): return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)
def near(first,second): return finite(first) and finite(second) and math.isclose(first,second,rel_tol=1e-12,abs_tol=1e-15)

def run_with_process_group(command,root,env,stdout,stderr):
    """Terminate the replay and its producer descendants before finalizing."""
    with subprocess.Popen(command,cwd=root,env=env,stdout=stdout,stderr=stderr,start_new_session=True) as process:
        try:
            process.wait(timeout=2400)
        except BaseException:
            try: os.killpg(process.pid,signal.SIGKILL)
            except ProcessLookupError: pass
            process.wait()
            raise
        return process.returncode

def safe_file(root,name):
    relative=Path(name)
    require(not relative.is_absolute() and '..' not in relative.parts,'Unsafe input path: '+str(name))
    target=root/relative
    require(target.is_file() and target.resolve().is_relative_to(root.resolve()),'Missing/escaping input: '+str(name))
    require(not any((root/Path(*relative.parts[:i])).is_symlink() for i in range(1,len(relative.parts)+1)),
            'Symlink input: '+str(name))
    return target

def load_authenticated_driver(root,bindings):
    """Authenticate the driver's bytes before importing its pure verifier."""
    require(root.name==CHECKPOINT,'Memory-followup checkpoint identity required')
    for key,length in (('public_freeze_commit',40),('registration_sha256',64),('independent_manifest_sha256',64)):
        value=bindings[key]
        require(len(value)==length and all(c in '0123456789abcdef' for c in value),'Invalid explicit public binding: '+key)
    receipt=read(safe_file(root,'FREEZE_RECEIPT.json'))
    require(all(receipt.get(k)==v for k,v in bindings.items()),'Freeze receipt differs from explicit public bindings')
    require(sha(safe_file(root,'FULL_REGISTRATION.json'))==bindings['registration_sha256'],'Registration differs from explicit public hash')
    registration=read(root/'FULL_REGISTRATION.json')
    files=registration['files'];require(isinstance(files,dict),'Malformed registration file table')
    for name,pin in files.items(): require(sha(safe_file(root,name))==pin,'Frozen input changed: '+name)
    for name in ('code/replay_metric.py','code/metric_primary.py','independent/forced_metric.py','independent/PRESERVED_ARCHIVE_SCHEMA.json'):
        require(name in files,'Required frozen executable/schema missing: '+name)
    sys.dont_write_bytecode=True
    spec=importlib.util.spec_from_file_location('memory_metric_frozen_replay',root/'code/replay_metric.py')
    driver=importlib.util.module_from_spec(spec);spec.loader.exec_module(driver)
    inputs=driver.verify_inputs(root,bindings['public_freeze_commit'])
    require(all(inputs[k]==v for k,v in bindings.items()),'Frozen replay verifier differs from explicit bindings')
    return driver,inputs

def verify_registered_helper(inputs):
    name='code/reproduce_metric_expected_scientific_failure.py'
    require(name in inputs['frozen_files'] and sha(__file__)==inputs['frozen_files'][name],
            'Running control helper is not the frozen registered helper')

def verify_archive(path,run,schema):
    """Inspect every NPY header and stream CRCs; perform no mode evaluation."""
    import numpy as np
    require(schema['member_count']==631 and len(schema['members'])==631,'Frozen complete archive schema differs')
    coarse=run['setting']=='coarse';nodes=8192 if coarse else 16384;points=577 if coarse else 1153
    expected={m['name']:m for m in schema['members']}
    require(len(expected)==631,'Duplicate frozen archive schema member')
    with zipfile.ZipFile(path) as archive:
        names=archive.namelist()
        require(len(names)==631 and len(set(names))==631 and set(names)=={k+'.npy' for k in expected},'Incomplete/extra archive members: '+path.name)
        require(archive.testzip() is None,'Archive CRC failure: '+path.name)
        for name,member in expected.items():
            shape=tuple(nodes if x==16384 else points if x==1153 else x for x in member['shape'])
            with archive.open(name+'.npy') as stream:
                version=np.lib.format.read_magic(stream)
                require(version in ((1,0),(2,0)),'Unsupported NPY version')
                actual_shape,fortran,dtype=(np.lib.format.read_array_header_1_0(stream) if version==(1,0) else np.lib.format.read_array_header_2_0(stream))
                require(actual_shape==shape and fortran==member['fortran_order'] and dtype.str==member['dtype'] and not dtype.hasobject,
                        'Archive header differs: '+path.name+'/'+name)
                require(archive.getinfo(name+'.npy').file_size==stream.tell()+math.prod(shape)*dtype.itemsize,'Truncated/extra NPY payload: '+name)
    return {'path':path.name,'sha256':sha(path),'members':631,'header_contract':'PASS','streamed_crc':'PASS'}

def qualify_saved_failure(directory,exit_code,bindings,experiment,schema,timed_out=False,expected_plan=None):
    """Inspect saved outputs only. The helper does not loosen scientific gates."""
    directory=Path(directory);fresh=directory/'fresh'
    require(isinstance(exit_code,int) and exit_code>0 and not timed_out,'Expected ordinary nonzero replay exit, not success/signal/timeout')
    ledger=read(directory/'VALIDATION.json');plan=read(directory/'COMMAND_PLAN.json')
    if expected_plan is not None: require(plan==expected_plan,'Saved command plan differs from authenticated frozen driver')
    require(ledger['status']=='FAIL' and ledger['source_verification_before']==ledger['source_verification_after']=='PASS','Underlying replay/source status differs')
    require(len(plan)==ledger['planned_commands']==32,'Exactly 32 frozen planned commands required')
    commands=ledger['commands']
    require(len(commands)==ledger['attempted_commands']==28 and ledger['successful_commands']==27,'Failure must follow 27 completed commands')
    require([c['check'] for c in commands]==[c['check'] for c in plan[:28]],'Executed command prefix differs from frozen plan')
    require(all(all(c.get(k)==p.get(k) for k in ('check','arguments','physical_producer','timeout_seconds','command')) for c,p in zip(commands,plan)),
            'Executed command arguments differ from saved frozen plan')
    require(all(c['status']=='PASS' and c['exit_code']==0 and c['timed_out'] is False for c in commands[:-1]),'A prerequisite/primary command failed')
    last=commands[-1]
    require(last['check']=='independent_modes' and last['status']=='FAIL' and isinstance(last['exit_code'],int) and last['exit_code']>0 and last['timed_out'] is False,'Unexpected failing command or exit')
    require([c['check'] for c in commands if c['physical_producer']]==['primary','independent_modes'],'Physical command coverage differs')
    require(ledger['physical_producer_commands_attempted']==2 and ledger['physical_producer_commands_completed']==1,'Physical completion accounting differs')
    require(ledger['failure']['exception']=="RuntimeError('independent_modes failed; preserved stdout/stderr/exit and fresh outputs')",'Replay failure class/message differs')
    for path in ('validation','validation_optimized','report','figures'):
        require(not (fresh/path).exists(),'Post-failure calibration output unexpectedly exists: '+path)
    require(read(directory/'INPUT_MANIFEST.json')['registration_sha256']==bindings['registration_sha256'],'Replay input binding differs')
    for label in ('public_freeze_commit','registration_sha256','independent_manifest_sha256'):
        require(read(directory/'INPUT_MANIFEST.json')[label]==bindings[label],'Replay input manifest differs: '+label)
    output_pins=ledger['fresh_output_sha256']
    actual={p.relative_to(fresh).as_posix() for p in fresh.rglob('*') if p.is_file()}
    require(actual==set(output_pins),'Saved fresh output inventory differs from replay ledger')
    for name,pin in output_pins.items(): require(sha(safe_file(fresh,name))==pin,'Saved fresh output changed: '+name)
    primary=read(fresh/'primary/results.json');modes=read(fresh/'independent/results.json')
    expected_rows=[(s,e) for s in SOURCES for e in OBSERVATIONS]
    require(primary['status']=='completed' and [(r['source'],r['eta']) for r in primary['rows']]==expected_rows,'Primary must complete all 12 rows')
    require(experiment['primary']['producer_wall_seconds']==900,'Frozen primary wall contract differs')
    require(finite(primary['elapsed_seconds']) and 0<=primary['elapsed_seconds']<=experiment['primary']['producer_wall_seconds'],'Primary wall budget failed')
    for key in ('public_freeze_commit','registration_sha256'):
        require(primary['provenance'][key]==bindings[key],'Primary provenance differs: '+key)
    require(modes['status']=='failed' and modes['active_run'] is None,'Independent producer must fail only after full runs')
    require(modes['freeze_commit']==bindings['public_freeze_commit'] and modes['provenance']['registration_sha256']==bindings['registration_sha256'] and modes['provenance']['manifest_sha256']==bindings['independent_manifest_sha256'],'Independent provenance differs')
    config=experiment['independent_configuration'];gates=experiment['gates']['independent']
    require(modes['configuration']==config,'Independent frozen configuration differs')
    require(config['wall_budget_seconds']==900 and config['rss_budget_kib']==262144 and gates['ward_endpoint']==2e-6,'Frozen resource/Ward contract differs')
    resources=modes['resources']
    require(finite(resources['elapsed_seconds']) and 0<=resources['elapsed_seconds']<=900 and finite(resources['peak_rss_kib']) and 0<resources['peak_rss_kib']<=262144,'Independent resource budget failed')
    require([(r['source'],r['setting']) for r in modes['runs']]==list(RUNS),'All four ordered independent runs required')
    require([(r['source'],r['eta']) for r in modes['rows']]==expected_rows,'Combined independent row coverage differs')
    archives=[];lookup={};maxima={k:0.0 for k in QUANTITIES};ward_max=ward_refine=0.0;failures=[]
    for run in modes['runs']:
        source,setting=run['source'],run['setting'];lookup[(source,setting)]=run
        require([(r['source'],r['eta']) for r in run['rows']]==[(source,e) for e in OBSERVATIONS],'Incomplete independent run observations')
        coarse=setting=='coarse';steps=576 if coarse else 1152;nodes=8192 if coarse else 16384
        require(run['time_steps']==steps and run['momentum_nodes']==nodes and run['direct_stress_history_points']==steps+1 and run['dt']==1/(128 if coarse else 256) and run['momentum_panel_width']==(0.5 if coarse else 0.25),'Independent run grid differs')
        require(all(finite(run[k]) and 0<=run[k]<=gates['wronskian'] for k in ('wronskian_max_scaled','physical_wronskian_max_scaled')),'Wronskian failure is not the expected control')
        name=f'metric_modes_{source}_{setting}.npz'
        require(run['archive']['path']==name,'Unexpected archive path')
        path=safe_file(fresh/'independent',name)
        require(sha(path)==run['archive']['sha256'],'Saved archive hash differs')
        archives.append(verify_archive(path,run,schema))
    combined={(r['source'],r['eta']):r for r in modes['rows']}
    for source in SOURCES:
        for coarse,fine in zip(lookup[(source,'coarse')]['rows'],lookup[(source,'fine')]['rows']):
            require([p['K'] for p in coarse['finite_k']]==[p['K'] for p in fine['finite_k']]==list(CUTOFFS),'Incomplete cutoff grid')
            for c,f in zip(coarse['finite_k'],fine['finite_k']):
                view=next(p for p in combined[(source,fine['eta'])]['finite_k'] if p['K']==f['K'])
                require(view['values']==f['values'],'Combined independent values differ from fine run')
                require(all(view['ward'][k]==f['ward'][k] for k in ('ledger','direct_endpoint','endpoint_difference')),'Combined Ward values differ from fine run')
                for key in QUANTITIES:
                    require(finite(c['values'][key]) and finite(f['values'][key]),'Nonfinite saved observable')
                    value=abs(f['values'][key]-c['values'][key]);maxima[key]=max(maxima[key],value)
                    if key in gates['refinement'] and value>gates['refinement'][key]: failures.append(f"{source}/{fine['eta']}/{f['K']}/{key} refinement")
                for point in (c,f):
                    ward=point['ward'];require(near(ward['endpoint_difference'],abs(ward['ledger']-ward['direct_endpoint'])) and near(ward['direct_endpoint'],point['values']['rho']),'Inconsistent saved Ward endpoint')
                wd=abs(f['ward']['ledger']-c['ward']['ledger']);ward_max=max(ward_max,f['ward']['endpoint_difference']);ward_refine=max(ward_refine,wd)
                if f['ward']['endpoint_difference']>gates['ward_endpoint']: failures.append(f"{source}/{fine['eta']}/{f['K']} Ward endpoint")
                if wd>gates['ward_refinement']: failures.append(f"{source}/{fine['eta']}/{f['K']} Ward refinement")
    require(failures and modes['gate_failures']==failures,'Internal failure list differs from saved data and frozen gates')
    require(modes['exception']['type']=='RuntimeError' and modes['exception']['message']=='Registered internal gates failed: '+', '.join(failures),'Unexpected producer failure class/message')
    require(all(near(modes['maximum_refinement_differences'][k],v) for k,v in maxima.items()) and near(modes['maximum_ward_endpoint_difference'],ward_max) and near(modes['maximum_ward_refinement_difference'],ward_refine),'Declared maxima differ from saved runs')
    witness=next(p for p in lookup[('positive_B','fine')]['rows'][-1]['finite_k'] if p['K']==256)
    require(witness['ward']['endpoint_difference']>2e-6 and 'positive_B/-1.5/256 Ward endpoint' in failures,'Known fine Ward failure was not reproduced')
    return {'control_status':'EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED','underlying_scientific_status':'FAIL',
            'registered_metric_calibration_passed':False,'primary_completed_rows':12,'independent_completed_runs':4,
            'planned_commands':32,'attempted_commands':28,'successful_commands':27,'resources':resources,
            'gate_failures':failures,'maximum_refinement_differences':maxima,'maximum_ward_endpoint_difference':ward_max,
            'maximum_ward_refinement_difference':ward_refine,'known_Ward_witness':{'source':'positive_B','setting':'fine','eta':-1.5,'K':256,'difference':witness['ward']['endpoint_difference'],'frozen_gate':2e-6},
            'complete_archives':archives,'validation_run_status':'NOT_REACHED','cross_route_calibration_status':'NOT_ESTABLISHED',
            'interpretation':'Control PASS reproduces a registered scientific FAIL; no source, gate, grid, model, or numerical value is changed.'}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--public-freeze-commit',required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--independent-manifest-sha256',required=True)
    args=parser.parse_args();root=args.checkpoint.resolve();output=args.output.resolve()
    bindings={k:getattr(args,k) for k in ('public_freeze_commit','registration_sha256','independent_manifest_sha256')}
    require(not output.exists() and not output.is_relative_to(root),'Fresh external control output required')
    output.mkdir(parents=True)
    state={'status':'RUNNING','started_utc':stamp(),'helper_sha256':sha(__file__),'explicit_public_bindings':bindings,
           'underlying_scientific_status':'NOT_YET_TESTED','registered_metric_calibration_passed':False,
           'source_verification_before':'NOT_STARTED','source_verification_after':'NOT_STARTED'}
    driver=inputs=None
    try:
        require(sys.version_info[:2]==(3,12) and sys.flags.optimize==0,'Normal Python3.12 required')
        for name,pin in VERSIONS.items(): require(importlib.metadata.version(name)==pin,'Pinned dependency differs: '+name)
        driver,inputs=load_authenticated_driver(root,bindings);state['source_verification_before']='PASS'
        verify_registered_helper(inputs)
        write(output/'INPUT_MANIFEST.json',inputs)
        replay=output/'replay'
        require(len(driver.command_plan(root,replay/'fresh',inputs))==32,'Frozen replay plan must contain exactly 32 commands')
        expected_plan=[dict(c,command=[sys.executable,*c['arguments']],cwd=str(replay/'checkpoint'))
                       for c in driver.command_plan(replay/'checkpoint',replay/'fresh',inputs)]
        command=[sys.executable,str(root/'code/replay_metric.py'),'--output',str(replay),'--freeze-commit',bindings['public_freeze_commit']]
        state['command']=command
        env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',MPLBACKEND='Agg')
        record={'timed_out':False};started=time.monotonic()
        with (output/'replay.stdout.log').open('xb') as stdout,(output/'replay.stderr.log').open('xb') as stderr:
            try:
                record['exit_code']=run_with_process_group(command,root,env,stdout,stderr)
            except subprocess.TimeoutExpired:
                record.update(exit_code=None,timed_out=True);raise
            finally:
                record['elapsed_seconds']=time.monotonic()-started;write(output/'replay.exit.json',record)
        require(driver.verify_inputs(root,bindings['public_freeze_commit'])==inputs,'Frozen sources changed during control')
        require(driver.verify_inputs(replay/'checkpoint',bindings['public_freeze_commit'])==inputs,'Replay copied sources changed')
        state['source_verification_after']='PASS'
        state['underlying_replay_status']=read(replay/'VALIDATION.json')['status']
        state['underlying_scientific_status']='FAIL_OR_INCOMPLETE_UNQUALIFIED'
        if (replay/'fresh/independent/results.json').is_file():
            state['underlying_independent_status']=read(replay/'fresh/independent/results.json')['status']
        state.update(qualify_saved_failure(replay,record['exit_code'],bindings,read(root/'EXPERIMENT.json'),read(root/'independent/PRESERVED_ARCHIVE_SCHEMA.json'),record['timed_out'],expected_plan))
        state['status']='PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL'
    except BaseException as exc:
        state.update(status='FAIL_EXPECTED_SCIENTIFIC_FAILURE_CONTROL',failure={'exception':repr(exc),'traceback':traceback.format_exc()});raise
    finally:
        if inputs is not None and driver is not None:
            try:
                require(driver.verify_inputs(root,bindings['public_freeze_commit'])==inputs,'Frozen sources changed during control');state['source_verification_after']='PASS'
            except BaseException as exc:
                state.update(status='FAIL_EXPECTED_SCIENTIFIC_FAILURE_CONTROL',source_verification_after='FAIL',source_failure=repr(exc))
        state['completed_utc']=stamp()
        state['output_sha256']={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='EXPECTED_SCIENTIFIC_FAILURE.json'}
        write(output/'EXPECTED_SCIENTIFIC_FAILURE.json',state)
    require(state['status']=='PASS_EXPECTED_SCIENTIFIC_FAILURE_CONTROL','Expected scientific failure control did not pass')
    print(json.dumps({k:state[k] for k in ('status','control_status','underlying_scientific_status','registered_metric_calibration_passed')}))

if __name__=='__main__': main()
