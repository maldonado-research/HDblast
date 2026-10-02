#!/usr/bin/env python3
"""Self-contained replay of the publicly frozen matched-stress checkpoint.

Install pinned requirements with Python 3.12. Put this driver in checkpoint/code
and support inputs in checkpoint/replay_inputs. No repository path is needed.
Each physical producer runs once, normally; both validators and six analytic/
format checks run normally and optimized. No results or gates are altered.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

FREEZE='4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d'
REGISTRATION='4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893'
INDEPENDENT='4933c5fe2bd9d3fc2ea154371ce0fdf9d8066e3546a672975d1e007a9151cf42'
SUPPORT={
 'derive_local_trace.py':'25a9703f720fe068e97b1fd6731808751d66d90960fd40b0607e79a48276e0f0',
 'prior_variance_tail_bounds.py':'642fe91d3c241d38db0c874599edb49217a5ef05f266dff3a830141fb9e27b35'}
SOURCE_MAP={
 'research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/COMMON_ACTION_SMOOTH_FRW.md':'inherited/COMMON_ACTION_SMOOTH_FRW.md',
 'research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/code/derive_local_trace.py':'replay_inputs/derive_local_trace.py',
 'research/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/theory/PROSPECTIVE_MATCHED_STRESS_RESPONSE.md':'inherited/MATCHED_STRESS_RESPONSE_PROTOCOL.md',
 'research/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/theory/verify_matched_stress_proposal.py':'inherited/verify_prior_stress_algebra.py',
 'research/HDBLAST_CHECKPOINT_20261002_CAUSAL_RESPONSE/theory/tail_bounds.py':'replay_inputs/prior_variance_tail_bounds.py'}
EXPECTED_COUNTS={'core_comparisons':180,'extras':108,'direct_stress_closed':72,'finite_K_trace':36,
 'raw_reconstructed_observation_K_points':72,'raw_Ward_endpoints':72,'Ward_refinement_points':36,
 'core_refinement_comparisons':180,'continuum_comparisons':288,'controls':20}


def require(condition,message):
    if not condition:raise RuntimeError(message)

def read(path):return json.loads(Path(path).read_text())
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def stamp():return datetime.now(timezone.utc).isoformat()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')

def safe_path(root,name):
    path=(root/name).resolve()
    require(not Path(name).is_absolute() and path.is_relative_to(root.resolve()),'Input path escape: '+name)
    require(path.is_file() and not (root/name).is_symlink(),'Missing/nonregular/symlink input: '+name)
    return path

def verify_inputs(package):
    require(digest(package/'FULL_REGISTRATION.json')==REGISTRATION,'Full prospective registration changed')
    frozen=read(package/'FULL_REGISTRATION.json')['files']
    for name,pin in frozen.items():require(digest(safe_path(package,name))==pin,'Frozen file changed: '+name)
    require(digest(package/'independent/MANIFEST.json')==INDEPENDENT,'Independent manifest changed')
    manifest=read(package/'independent/MANIFEST.json')
    for entry in manifest['files']:
        path=(package/'independent'/entry['path']).resolve()
        require(path.is_relative_to(package.resolve()) and path.is_file(),'Independent input escape/missing')
        require(digest(path)==entry['sha256'],'Independent input changed: '+entry['path'])
    for item in read(package/'PREREGISTRATION.json')['inherited_files'].values():
        require(digest(safe_path(package,item['local_copy']))==item['sha256'],'Inherited input changed')
    for name,pin in SUPPORT.items():require(digest(safe_path(package,'replay_inputs/'+name))==pin,'Standalone replay support changed: '+name)
    pins=read(package/'theory/SOURCE_PINS.json')['files']
    require({p['repository_relative_path'] for p in pins}==set(SOURCE_MAP),'Tail algebra source-pin membership changed')
    for item in pins:require(digest(safe_path(package,SOURCE_MAP[item['repository_relative_path']]))==item['sha256'],'Tail source-pin input changed')
    payload=None
    if (package/'MANIFEST.json').exists():
        payload=read(package/'MANIFEST.json');excluded=set(payload['excluded'])
        actual={p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.relative_to(package).as_posix() not in excluded}
        require(actual==set(payload['files']),'Packaged payload membership changed')
        for name,pin in payload['files'].items():require(digest(safe_path(package,name))==pin,'Payload changed: '+name)
    return {'registration_sha256':REGISTRATION,'independent_manifest_sha256':INDEPENDENT,
      'frozen_file_count':len(frozen),'frozen_files':frozen,'support_files':SUPPORT,
      'package_manifest_sha256':digest(package/'MANIFEST.json') if payload else None}

def build_source_mirror(work,mirror):
    for item in read(work/'theory/SOURCE_PINS.json')['files']:
        path=mirror/item['repository_relative_path'];path.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(work/SOURCE_MAP[item['repository_relative_path']],path)
        require(digest(path)==item['sha256'],'Source mirror digest mismatch')

def command_plan(work,fresh,mirror):
    commands=[]
    def pair(name,args):
        commands.append((name,args));commands.append((name+'_optimized',['-O',*args]))
    for optimize in (False,True):
        suffix='_optimized' if optimize else '';flag=['-O'] if optimize else []
        commands.extend([
          ('primary_algebra'+suffix,[*flag,'code/verify_primary_algebra.py','--output',str(fresh/('PRIMARY_ALGEBRA'+suffix+'.json'))]),
          ('source_jet_bridge'+suffix,[*flag,'code/verify_source_jet_bridge.py','--tail-module',str(work/'theory/pulse-norms/pulse_derivative_budgets.py'),'--output',str(fresh/('SOURCE_JET_BRIDGE'+suffix+'.json'))]),
          ('stress_tail_algebra'+suffix,[*flag,'theory/verify_stress_tail_algebra.py','--repo-root',str(mirror),'--output',str(fresh/('STRESS_TAIL_ALGEBRA'+suffix+'.json'))]),
          ('pulse_exact_proof'+suffix,[*flag,'theory/pulse-norms/verify_pulse_derivative_budgets.py','--output',str(fresh/('PULSE_EXACT_PROOF'+suffix+'.json'))]),
          ('finite_cutoff_ward_algebra'+suffix,[*flag,'theory/independent_review/derive_finite_cutoff_ward.py']),
          ('validator_guards'+suffix,[*flag,'code/test_validator_guards.py'])])
    commands.extend([
      ('primary',['code/stress_primary.py','--output',str(fresh/'primary'),'--registration-sha256',REGISTRATION,'--public-freeze-commit',FREEZE]),
      ('independent_modes',['independent/forced_stress.py','--output-dir',str(fresh/'independent'),'--freeze-commit',FREEZE,'--manifest-sha256',INDEPENDENT])])
    validator=['code/validate_stress.py','--primary',str(fresh/'primary/results.json'),'--modes',str(fresh/'independent/results.json'),'--registration-sha256',REGISTRATION,'--public-freeze-commit',FREEZE]
    commands.extend([('validation',[*validator,'--output',str(fresh/'validation')]),('validation_optimized',['-O',*validator,'--output',str(fresh/'validation_optimized')])])
    require(len(commands)==16,'Unexpected registered command plan')
    return commands

def compare_checks(normal,optimized):
    require(normal['status']==optimized['status']=='PASS','Validator status failed')
    require(normal['counts']==optimized['counts']==EXPECTED_COUNTS,'Validator coverage changed')
    n={k:v for k,v in normal.items() if k!='python_optimization'}
    o={k:v for k,v in optimized.items() if k!='python_optimization'}
    require(n==o,'Normal/optimized validation evidence differs')
    require(normal['python_optimization']==0 and optimized['python_optimization']==1,'Optimization modes changed')
    controls=normal['controls']
    require(len(controls)==20 and sum(c.get('rejected',False) for c in controls)==5,'Control classes/counts changed')
    require(all(c.get('nonzero_algebraic_sensitivity_witness',False) for c in controls if 'maximum_absolute_residual' in c),'Control sensitivity missing')

def final_checks(fresh,logs):
    totals={}
    for suffix in ('','_optimized'):
        algebra=read(fresh/('PRIMARY_ALGEBRA'+suffix+'.json'));bridge=read(fresh/('SOURCE_JET_BRIDGE'+suffix+'.json'))
        tail=read(fresh/('STRESS_TAIL_ALGEBRA'+suffix+'.json'));pulse=read(fresh/('PULSE_EXACT_PROOF'+suffix+'.json'))
        require(algebra['passed'] and algebra['check_count']==30,'Primary algebra coverage')
        require(bridge['passed'] and bridge['check_count']==12 and all(c['passed'] for c in bridge['checks']),'Independent source-jet bridge coverage')
        require(tail['status']=='PASS' and tail['identity_count']==25 and tail['mutation_count']==6,'Stress tail algebra coverage')
        require(pulse['status']=='PASS' and pulse['exact_identity_count']==18 and pulse['critical_root_certificate_count']==6 and pulse['detected_mutation_count']==7,'Pulse proof coverage')
        guard=(logs/('validator_guards'+suffix+'.log')).read_text()
        require('PASS 21 synthetic algebra/format/static checks' in guard,'Guard check coverage')
        totals[suffix or 'normal']={'primary_algebra':algebra['check_count'],'source_jet_bridge':bridge['check_count'],'tail_identities':tail['identity_count'],'tail_mutations':tail['mutation_count'],'pulse_identities':pulse['exact_identity_count'],'pulse_root_certificates':pulse['critical_root_certificate_count'],'pulse_mutations':pulse['detected_mutation_count'],'synthetic_guards':21}
    require(len(read(fresh/'primary/results.json')['rows'])==12 and len(read(fresh/'independent/results.json')['rows'])==12 and len(read(fresh/'independent/results.json')['runs'])==4,'Producer coverage changed')
    compare_checks(read(fresh/'validation/CHECKS.json'),read(fresh/'validation_optimized/CHECKS.json'))
    return {'response_points':12,'independent_runs':4,'validator':EXPECTED_COUNTS,'analytic_checks_by_mode':totals}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--plan-only',action='store_true',help='Authenticate and stage the exact command plan without executing any child command')
    args=parser.parse_args();require(not sys.flags.optimize,'Run replay driver normally; it invokes -O checks itself')
    package=Path(__file__).resolve().parents[1];output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(package),'Choose a fresh external replay output directory')
    # A fresh receipt retains even dependency/input failures, before any child runs.
    inputs=None;output.mkdir(parents=True,exist_ok=False)
    work=output/'checkpoint';fresh=output/'fresh';logs=output/'logs';mirror=output/'source_pin_mirror'
    state={'status':'RUNNING','started_utc':stamp(),'public_freeze_commit':FREEZE,'registration_sha256':REGISTRATION,'independent_manifest_sha256':INDEPENDENT,'driver_sha256':digest(__file__),'commands':[],'scope':'Fixed geometry, linear homogeneous minimal stress/current only','coupled_geometry_executed':False,'stability_established':False,'heating_established':False,'certified_total_numerical_error':False}
    write(output/'VALIDATION.json',state)
    try:
        require(sys.version_info[:2]==(3,12),'Python 3.12 required')
        for name,version in {'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0'}.items():require(importlib.metadata.version(name)==version,'Pinned dependency mismatch: '+name)
        inputs=verify_inputs(package);write(output/'INPUT_MANIFEST.json',inputs)
        state['source_verification_before']='PASS'
        shutil.copytree(package,work,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        require(verify_inputs(work)==inputs,'Copied checkpoint input verification differs')
        fresh.mkdir();logs.mkdir();mirror.mkdir();build_source_mirror(work,mirror)
        env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',MPLBACKEND='Agg')
        commands=command_plan(work,fresh,mirror);state['planned_commands']=len(commands)
        write(output/'COMMAND_PLAN.json',[{'check':name,'command':[sys.executable,*command],'cwd':str(work),'timeout_seconds':1000} for name,command in commands])
        if args.plan_only:
            state.update(status='VERIFIED_PLAN_ONLY',executed_commands=0,physical_producer_commands_executed=0)
        else:
            for name,command in commands:
                receipt={'check':name,'command':[sys.executable,*command],'started_utc':stamp(),'status':'RUNNING','expected_exit_code':0,'timeout_seconds':1000}
                state['commands'].append(receipt);write(output/'EXECUTION.json',state['commands']);start=time.monotonic()
                with (logs/(name+'.log')).open('x') as stream:
                    try:
                        result=subprocess.run(receipt['command'],cwd=work,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=1000)
                        receipt.update(exit_code=result.returncode,timed_out=False,status='PASS' if result.returncode==0 else 'FAIL')
                    except subprocess.TimeoutExpired:
                        stream.write('\nReplay wrapper timed out at 1000 seconds; producer budget is unchanged at 900 seconds.\n')
                        receipt.update(exit_code=None,timed_out=True,status='FAIL')
                    except BaseException as exc:
                        receipt.update(exit_code=None,timed_out=False,status='FAIL',exception=repr(exc));raise
                    finally:
                        receipt.update(elapsed_seconds=time.monotonic()-start,completed_utc=stamp());write(output/'EXECUTION.json',state['commands'])
                require(receipt['status']=='PASS',name+' failed; retained log '+str(logs/(name+'.log')))
                require(verify_inputs(work)==inputs and verify_inputs(package)==inputs,'Frozen source changed during command '+name)
                print(name+' PASS',flush=True)
            state.update(status='PASS',counts=final_checks(fresh,logs),executed_commands=len(state['commands']),physical_producer_commands_executed=2)
    except BaseException as exc:
        state.update(status='FAIL',failure={'exception':repr(exc),'traceback':traceback.format_exc()})
        raise
    finally:
        try:
            if inputs is None:
                state['source_verification_after']='NOT_STARTED; dependency or initial input validation failed before execution'
            else:
                require(verify_inputs(package)==inputs,'Original checkpoint changed during replay')
                if work.exists():require(verify_inputs(work)==inputs,'Copied checkpoint changed during replay')
                state['source_verification_after']='PASS'
        except BaseException as exc:
            state.update(status='FAIL',source_verification_after='FAIL',source_verification_failure={'exception':repr(exc),'traceback':traceback.format_exc()})
        state.update(completed_utc=stamp(),attempted_commands=len(state['commands']),successful_commands=sum(c['status']=='PASS' for c in state['commands']),physical_producer_commands_attempted=sum(c['check'] in ('primary','independent_modes') for c in state['commands']),physical_producer_commands_completed=sum(c['check'] in ('primary','independent_modes') and c['status']=='PASS' for c in state['commands']))
        if fresh.exists():state['fresh_output_sha256']={p.relative_to(fresh).as_posix():digest(p) for p in sorted(fresh.rglob('*')) if p.is_file()}
        write(output/'VALIDATION.json',state)
    require(state['status'] in ('PASS','VERIFIED_PLAN_ONLY'),'Replay post-run source verification failed')
    print(json.dumps({'status':state['status'],'planned_commands':state['planned_commands'],'attempted_commands':state['attempted_commands'],'successful_commands':state['successful_commands']},indent=2))

if __name__=='__main__':main()
