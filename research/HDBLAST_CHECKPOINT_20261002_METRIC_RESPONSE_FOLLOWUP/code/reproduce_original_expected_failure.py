#!/usr/bin/env python3
"""Reproduce the original frozen primary failure, without reclassifying its science.

A control PASS means the original IntegrationWarning failure was reproduced.
It does not mean the original experiment passed. The amended followup replay
is a distinct experiment with its own prospective registration and receipt.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

ORIGINAL='HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE'
EXPECTED={'public_freeze_commit':'57668b8fadd75df8738565e0bbd1eb852c1ebae8',
          'registration_sha256':'f8d6bbd17b540b46c5e4382ab21fd74e15f0b283956991a250f49994553a476d',
          'independent_manifest_sha256':'4760f2a7e864dbe81ad8911797499dcd9fcc9e37a45421a607a101d35e6ce507'}
PRODUCER_SHA='732345e68c99ed48f14715102aadc1c7ce2e50c30808466631630cfc71f901ba'
DRIVER_SHA='bb4980e49e30a3167220eebc546d21e6e65aa9c464e821bbb2437f4486766cd3'
VERSIONS={'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}


def require(ok,message):
    if not ok:raise RuntimeError(message)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
def stamp():return datetime.now(timezone.utc).isoformat()


def qualify_saved_failure(directory,exit_code,timed_out=False):
    """Pure inspection of saved evidence; never invokes a producer."""
    directory=Path(directory)
    require(exit_code is not None and exit_code!=0 and not timed_out,'Expected ordinary nonzero primary exit, not success or timeout')
    require(not (directory/'results.json').exists(),'Original primary unexpectedly completed')
    failure=read(directory/'failure.json');partial=read(directory/'partial_results.json');started=read(directory/'started.json')
    require(failure['exception'].startswith('IntegrationWarning(') and 'roundoff error' in failure['exception'],'Original failure class/message differs')
    require('log_history' in failure['traceback'] and "weight='alg-logb'" in failure['traceback'],'Original failure was not the frozen weighted-log route')
    require(partial['active']=={'source':'positive_B','eta':-4.5},'Original failure source/time differs')
    require([(row['source'],row['eta']) for row in partial['rows']]==[('positive_B',-5.5)],'Unexpected original completed row coverage')
    for name in ('registration_sha256','public_freeze_commit'):require(started[name]==EXPECTED[name],'Original started provenance differs: '+name)
    require(started['producer_sha256']==PRODUCER_SHA and started['python'].startswith('3.12.'),'Original producer/runtime pin differs')
    require(started['dependencies']=={k:VERSIONS[k] for k in ('numpy','scipy','mpmath','sympy')},'Original dependency pins differ')
    return {'control_status':'EXPECTED_FAILURE_CONFIRMED','original_scientific_status':'ORIGINAL_EXPERIMENT_FAILED',
            'failure_class':'IntegrationWarning','failure_kind':'weighted-log quadrature roundoff',
            'active_source':'positive_B','active_eta':-4.5,'completed_original_rows':1,
            'failure_sha256':sha(directory/'failure.json'),'partial_results_sha256':sha(directory/'partial_results.json'),
            'started_sha256':sha(directory/'started.json')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-checkpoint',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();root=args.original_checkpoint.resolve();output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(root),'Fresh external control output required')
    output.mkdir(parents=True)
    state={'status':'RUNNING','started_utc':stamp(),'helper_sha256':sha(__file__),
           'original_scientific_status':'ORIGINAL_EXPERIMENT_FAILED','followup_scientific_status':'NOT_TESTED_BY_THIS_CONTROL',
           'expected_original_pins':EXPECTED,'reproduced_physical_producer':'original primary only; no original independent-mode rerun',
           'source_verification_before':'NOT_STARTED','source_verification_after':'NOT_STARTED'}
    inputs=None;verifier=None
    try:
        require(sys.version_info[:2]==(3,12) and sys.flags.optimize==0,'Normal Python3.12 required')
        for name,pin in VERSIONS.items():require(importlib.metadata.version(name)==pin,'Pinned original dependency differs: '+name)
        require(root.name==ORIGINAL,'Original checkpoint identity required')
        require(sha(root/'code/metric_primary.py')==PRODUCER_SHA and sha(root/'code/replay_metric.py')==DRIVER_SHA,'Original frozen executable differs')
        sys.dont_write_bytecode=True
        spec=importlib.util.spec_from_file_location('original_metric_input_authenticator',root/'code/replay_metric.py')
        verifier=importlib.util.module_from_spec(spec);spec.loader.exec_module(verifier)
        inputs=verifier.verify_inputs(root)
        require(all(inputs[key]==value for key,value in EXPECTED.items()),'Original receipt/registration differs from expected failed experiment')
        write(output/'ORIGINAL_INPUT_MANIFEST.json',inputs);state['source_verification_before']='PASS'
        fresh=output/'primary'
        command=[sys.executable,str(root/'code/metric_primary.py'),'--output',str(fresh),
                 '--registration-sha256',EXPECTED['registration_sha256'],'--public-freeze-commit',EXPECTED['public_freeze_commit']]
        state['command']=command
        env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',MPLBACKEND='Agg')
        started=time.monotonic();record={'expected':'nonzero IntegrationWarning roundoff at positive_B eta-4.5','timed_out':False}
        with (output/'original-primary.stdout.log').open('xb') as stdout,(output/'original-primary.stderr.log').open('xb') as stderr:
            try:
                result=subprocess.run(command,cwd=root,env=env,stdout=stdout,stderr=stderr,timeout=1000)
                record['exit_code']=result.returncode
            except subprocess.TimeoutExpired:
                record.update(exit_code=None,timed_out=True);raise
            finally:
                record['elapsed_seconds']=time.monotonic()-started;write(output/'original-primary.exit.json',record)
        require(verifier.verify_inputs(root)==inputs,'Original sources changed during reproduction')
        state['source_verification_after']='PASS'
        state.update(qualify_saved_failure(fresh,record['exit_code'],record['timed_out']))
        state['status']='PASS_EXPECTED_FAILURE_CONTROL'
    except BaseException as exc:
        state.update(status='FAIL_EXPECTED_FAILURE_CONTROL',failure={'exception':repr(exc),'traceback':traceback.format_exc()});raise
    finally:
        if inputs is not None and verifier is not None:
            try:
                require(verifier.verify_inputs(root)==inputs,'Original sources changed during reproduction');state['source_verification_after']='PASS'
            except BaseException as exc:state.update(status='FAIL_EXPECTED_FAILURE_CONTROL',source_verification_after='FAIL',source_failure=repr(exc))
        state['completed_utc']=stamp();state['output_sha256']={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file() and p.name!='EXPECTED_FAILURE_CONTROL.json'}
        write(output/'EXPECTED_FAILURE_CONTROL.json',state)
    require(state['status']=='PASS_EXPECTED_FAILURE_CONTROL','Expected original failure control did not pass')
    print(json.dumps({k:state[k] for k in ('status','control_status','original_scientific_status','followup_scientific_status')}))

if __name__=='__main__':main()
