#!/usr/bin/env python3
"""Portable replay; authenticate all inputs before launching any child.

Run normally with --output NEW_EXTERNAL. --plan-only authenticates and writes
the exact command plan while executing zero subprocesses or physical calls.
"""
from __future__ import annotations
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

VERSIONS={'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}
QUANTITIES=('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
COUNTS={'core_comparisons':180,'extras':144,'raw_reconstructed_observation_K_points':72,'raw_Ward_endpoints':72,'Ward_refinement_points':36,'core_refinement_comparisons':180,'current_refinement_comparisons':36,'primary_memory_contact_reductions':432,'primary_differentiated_Ward':48,'continuum_tail_comparisons':324,'controls':15}


def require(ok,message):
    if not ok:raise RuntimeError(message)
def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,data):Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')
def stamp():return datetime.now(timezone.utc).isoformat()


def safe_path(root,name):
    target=(root/name).resolve()
    require(not Path(name).is_absolute() and target.is_relative_to(root.resolve()),'Input path escapes checkpoint: '+name)
    require(target.is_file() and not (root/name).is_symlink(),'Missing or symlink input: '+name)
    return target


def verify_inputs(root,override=None):
    root=root.resolve()
    receipt=read(safe_path(root,'FREEZE_RECEIPT.json'))
    freeze=receipt['public_freeze_commit']
    require(len(freeze)==40 and all(c in '0123456789abcdef' for c in freeze),'Invalid public freeze commit')
    require(override is None or override==freeze,'Explicit freeze override differs from receipt')
    require(sha(safe_path(root,'FULL_REGISTRATION.json'))==receipt['registration_sha256'],'Registration differs from freeze receipt')
    registration=read(root/'FULL_REGISTRATION.json')
    require(isinstance(registration['files'],dict),'Malformed frozen files')
    for name,pin in registration['files'].items():require(sha(safe_path(root,name))==pin,'Frozen input changed: '+name)
    require(sha(safe_path(root,'independent/MANIFEST.json'))==receipt['independent_manifest_sha256']==registration['independent_manifest_sha256'],'Independent manifest pin differs')
    manifest=read(root/'independent/MANIFEST.json')
    entries=manifest['files']
    entries=[{'path':k,'sha256':v} for k,v in entries.items()] if isinstance(entries,dict) else entries
    for item in entries:require(sha(safe_path(root,'independent/'+item['path']))==item['sha256'],'Independent input changed: '+item['path'])
    experiment=read(safe_path(root,'EXPERIMENT.json'))
    require(registration['frozen_configuration']==experiment and registration['frozen_gates']==experiment['gates'],'Frozen experiment or gates differ')
    require(experiment['runtime']['python_major_minor']==[3,12] and experiment['runtime']['physical_producers_optimized'] is False,'Frozen runtime differs')
    for name,version in VERSIONS.items():require(experiment['runtime'][name]==version,'Frozen dependency differs: '+name)
    payload_pin=None
    if (root/'MANIFEST.json').exists():
        payload=read(root/'MANIFEST.json')
        require(isinstance(payload['files'],dict) and 'FREEZE_RECEIPT.json' in payload['files'],'Final package manifest must protect freeze receipt')
        excluded=set(payload['excluded'])
        require('FREEZE_RECEIPT.json' not in excluded,'Freeze receipt cannot be excluded from package')
        actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.relative_to(root).as_posix() not in excluded}
        require(actual==set(payload['files']),'Package payload membership changed')
        for name,pin in payload['files'].items():require(sha(safe_path(root,name))==pin,'Package payload changed: '+name)
        payload_pin=sha(root/'MANIFEST.json')
    return {'public_freeze_commit':freeze,'registration_sha256':receipt['registration_sha256'],'independent_manifest_sha256':receipt['independent_manifest_sha256'],'freeze_receipt_sha256':sha(root/'FREEZE_RECEIPT.json'),'frozen_files':registration['files'],'package_manifest_sha256':payload_pin}


def command_plan(work,fresh,inputs):
    files=inputs['frozen_files']
    def locate(name):
        hits=[k for k in files if Path(k).name==name and 'reference_inputs' not in Path(k).parts and 'development' not in Path(k).parts]
        require(len(hits)==1,'Registered command filename must resolve uniquely: '+name)
        return str(work/hits[0])
    commands=[]
    def add(name,script,args=(),optimized=False,physical=False):
        commands.append({'check':name,'arguments':(['-O'] if optimized else [])+[locate(script),*map(str,args)],'physical_producer':physical,'timeout_seconds':1000})
    inventory=locate('FINITE_K_CONTACT_INVENTORY.json')
    for optimized in (False,True):
        suffix='_optimized' if optimized else ''
        algebra=fresh/('primary_algebra'+suffix)
        add('primary_algebra'+suffix,'derive_metric_contacts.py',['--output',algebra],optimized)
        add('primary_preflight'+suffix,'verify_primary_preflight.py',['--output',fresh/('PRIMARY_PREFLIGHT'+suffix+'.json'),'--algebra',algebra/'METRIC_CONTACT_ALGEBRA.json','--independent-inventory',inventory,'--experiment',work/'EXPERIMENT.json'],optimized)
        for check,script in [('primary_contact_action','metric_contact_algebra.py'),('primary_ward_trace','verify_metric_ward_trace.py'),('independent_generic_wkb','verify_metric_response.py'),('independent_inventory_agreement','verify_primary_inventory_agreement.py'),('generic_cse_agreement','verify_cse_specialization.py'),('actual_root_algebra','verify_prerequisites.py'),('metric_tail_algebra','verify_metric_tail_algebra.py')]:
            add(check+suffix,script,['--output',fresh/'proofs'/(check+suffix+'.json')],optimized)
        add('validator_guards'+suffix,'test_metric_validator_guards.py',optimized=optimized)
        add('raw_baseline_algebra'+suffix,'verify_raw_metric_baselines.py',optimized=optimized)
        add('inherited_source_pins'+suffix,'verify_inherited_sources.py',['--output',fresh/'proofs'/('INHERITED_SOURCES'+suffix+'.json')],optimized)
    add('independent_preflight','forced_metric.py',['--preflight'])
    add('metric_derivative_certificate','certify_metric_derivative_bounds.py',['--output',fresh/'proofs/DERIVATIVE_CERTIFICATE.json'])
    reg=inputs['registration_sha256'];freeze=inputs['public_freeze_commit'];manifest=inputs['independent_manifest_sha256']
    add('primary','metric_primary.py',['--output',fresh/'primary','--registration-sha256',reg,'--public-freeze-commit',freeze],physical=True)
    add('independent_modes','forced_metric.py',['--output-dir',fresh/'independent','--registration',work/'FULL_REGISTRATION.json','--registration-sha256',reg,'--manifest-sha256',manifest,'--freeze-commit',freeze],physical=True)
    validation=['--primary',fresh/'primary/results.json','--modes',fresh/'independent/results.json','--registration-sha256',reg,'--public-freeze-commit',freeze]
    for optimized in (False,True):
        suffix='_optimized' if optimized else ''
        add('validation'+suffix,'validate_metric.py',[*validation,'--output',fresh/('validation'+suffix)],optimized)
    add('summary','summarize_metric.py',['--checkpoint',work,'--primary',fresh/'primary/results.json','--modes',fresh/'independent/results.json','--checks',fresh/'validation/CHECKS.json','--checks-optimized',fresh/'validation_optimized/CHECKS.json','--output',fresh/'report'])
    add('figures','plot_metric.py',['--summary',fresh/'report/SUMMARY.json','--output',fresh/'figures'])
    optional=[k for k in files if Path(k).name=='audit_saved_metric.py' and 'development' not in Path(k).parts]
    if optional:
        require(len(optional)==1,'Archive audit source ambiguous')
        for optimized in (False,True):
            suffix='_optimized' if optimized else ''
            add('saved_mode_audit'+suffix,'audit_saved_metric.py',['--checkpoint',work,'--modes',fresh/'independent/results.json','--output',fresh/('SAVED_MODE_AUDIT'+suffix+'.json')],optimized)
    require(sum(c['physical_producer'] for c in commands)==2,'Exactly two physical producer commands required')
    require(all(not c['arguments'][0]=='-O' for c in commands if c['physical_producer']),'Physical producers must be unoptimized')
    return commands


def final_checks(work,fresh,logs):
    analytic=[]
    for path in sorted((fresh/'proofs').glob('*.json')):
        report=read(path)
        require(report.get('passed') is True or report.get('status')=='PASS','Analytic evidence failed: '+path.name)
        analytic.append(path.name)
    require(len(analytic)==17,'Expected analytic evidence coverage differs')
    for suffix in ('','_optimized'):
        algebra=read(fresh/('primary_algebra'+suffix)/'METRIC_CONTACT_ALGEBRA.json')
        preflight=read(fresh/('PRIMARY_PREFLIGHT'+suffix+'.json'))
        require(algebra['passed'] and len(algebra['checks'])==20,'Primary exact algebra coverage differs')
        require(preflight['passed'] and preflight['check_count']==45,'Primary exact preflight coverage differs')
        guards=read(logs/('validator_guards'+suffix+'.stdout.log'))
        baseline=read(logs/('raw_baseline_algebra'+suffix+'.stdout.log'))
        require(guards['status']=='PASS' and guards['checks']>=66 and guards['physical_evaluations']==0,'Validator synthetic guards missing')
        require(baseline['status']=='PASS' and len(baseline['identities'])==3 and len(baseline['mutations'])==6 and baseline['physical_evaluations']==0,'Stable raw baseline proof coverage differs')
    independent_preflight=read(logs/'independent_preflight.stdout.log')
    require(independent_preflight['status']=='pure_input_preflight_passed' and independent_preflight['physical_evaluations']==0 and len(independent_preflight['scheduled_runs'])==4,'Independent input preflight coverage differs')
    primary=read(fresh/'primary/results.json');modes=read(fresh/'independent/results.json')
    require(primary['status']=='completed' and modes['status']=='passed_internal_gates','Producer completion failed')
    require(len(primary['rows'])==len(modes['rows'])==12 and len(modes['runs'])==4,'Producer coverage differs')
    normal=read(fresh/'validation/CHECKS.json');optimized=read(fresh/'validation_optimized/CHECKS.json')
    require(normal['status']==optimized['status']=='PASS','Validator failed')
    require(normal['counts']==optimized['counts']==COUNTS,'Validator coverage differs')
    require(normal['python_optimization']==0 and optimized['python_optimization']==1,'Validator optimization modes differ')
    require({k:v for k,v in normal.items() if k!='python_optimization'}=={k:v for k,v in optimized.items() if k!='python_optimization'},'Normal/optimized validation evidence differs')
    summary=read(fresh/'report/SUMMARY.json')
    require(summary['status']=='registered_metric_calibration_passed' and len(summary['registered_rows'])==12,'Summary classification or coverage differs')
    require(summary['postrun_only'] and summary['new_physical_evaluations']==0,'Summary must only consume saved data')
    require(summary['normalizations']==primary['normalization'] and summary['epsilon']==primary['epsilon'],'Summary normalization differs')
    require(summary['continuum_tail_comparisons']==normal['continuum_tail_comparisons'] and summary['controls']==normal['controls'],'Summary continuum/control evidence differs')
    mode_rows={(r['source'],r['eta']):r for r in modes['rows']}
    presented={(r['source'],r['eta']):r for r in summary['registered_rows']}
    for row in primary['rows']:
        view=presented[(row['source'],row['eta'])]
        require(view['continuum_values']==row['continuum']['values'] and view['source_jet_over_epsilon']==row['source_jet_over_epsilon'] and view['canonical_forcing_jet_over_epsilon']==row['canonical_forcing_jet_over_epsilon'],'Summary source/continuum differs')
        direct=mode_rows[(row['source'],row['eta'])]
        for point,first,second in zip(view['finite_k'],row['finite_k'],direct['finite_k']):
            require(point['K']==first['K']==second['K'] and point['primary_values']==first['values'] and point['direct_mode_values']==second['values'],'Summary finiteK values differ')
            require(point['quadrature_error_estimates']==first['quadrature_error_estimates'] and point['mode_refinement_differences']==second['refinement_differences'] and point['ward']==second['ward'],'Summary uncertainty/ledger differs')
    certificate=read(fresh/'proofs/DERIVATIVE_CERTIFICATE.json');frozen_certificate=read(work/'theory/DERIVATIVE_CERTIFICATE.json')
    require(certificate==frozen_certificate,'Regenerated derivative certificate differs from frozen input')
    expected={n+'.'+e for n in ('metric_response','metric_cutoff','metric_checks') for e in ('png','svg','pdf')}
    provenance=read(fresh/'figures/FIGURE_PROVENANCE.json')
    require(set(provenance['figures_sha256'])==expected,'Exactly nine figure artifacts required')
    require({p.name for p in (fresh/'figures').iterdir() if p.suffix in ('.png','.svg','.pdf')}==expected,'Figure directory artifact coverage differs')
    require(provenance['summary_sha256']==sha(fresh/'report/SUMMARY.json') and provenance['renderer_sha256']==sha(work/'postrun/plot_metric.py'),'Figure provenance differs')
    for name,pin in provenance['figures_sha256'].items():require(sha(fresh/'figures'/name)==pin,'Figure output changed')
    require((fresh/'report/NUMERICAL_RESULTS.md').is_file(),'Numerical narrative missing')
    return {'response_points':12,'finite_K_points':36,'independent_runs':4,'validator':COUNTS,'analytic_output_files':len(analytic),'figure_artifacts':9}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--freeze-commit')
    parser.add_argument('--plan-only',action='store_true')
    args=parser.parse_args()
    require(sys.flags.optimize==0,'Replay driver must run normally')
    package=Path(__file__).resolve().parents[1];output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(package),'Choose fresh output outside checkpoint')
    output.mkdir(parents=True,exist_ok=False)
    work=output/'checkpoint';fresh=output/'fresh';logs=output/'logs'
    state={'status':'RUNNING','started_utc':stamp(),'driver_sha256':sha(__file__),'commands':[],'physical_experiments':2,'scope':'Homogeneous linear metric calibration at H1,r2,xi0 with fixed physical scalar and incoming BD state','coupled_evolution':False,'heating_or_stability_established':False,'certified_total_numerical_error':False}
    inputs=None
    write(output/'VALIDATION.json',state)
    try:
        require(sys.version_info[:2]==(3,12),'Python3.12 required')
        for name,version in VERSIONS.items():require(importlib.metadata.version(name)==version,'Pinned dependency mismatch: '+name)
        inputs=verify_inputs(package,args.freeze_commit);write(output/'INPUT_MANIFEST.json',inputs)
        state['source_verification_before']='PASS'
        shutil.copytree(package,work,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        require(verify_inputs(work,args.freeze_commit)==inputs,'Copied source authentication differs')
        fresh.mkdir();logs.mkdir();(fresh/'proofs').mkdir()
        commands=command_plan(work,fresh,inputs)
        write(output/'COMMAND_PLAN.json',[dict(c,command=[sys.executable,*c['arguments']],cwd=str(work)) for c in commands])
        state['planned_commands']=len(commands)
        if args.plan_only:
            state.update(status='VERIFIED_PLAN_ONLY',executed_commands=0,physical_producer_commands_executed=0)
        else:
            env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONOPTIMIZE='0',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',BLIS_NUM_THREADS='1',MPLBACKEND='Agg')
            for command in commands:
                name=command['check'];started=time.monotonic()
                record=dict(command,command=[sys.executable,*command['arguments']],started_utc=stamp(),status='RUNNING',expected_exit_code=0)
                state['commands'].append(record);write(output/'EXECUTION.json',state['commands'])
                with (logs/(name+'.stdout.log')).open('xb') as stdout,(logs/(name+'.stderr.log')).open('xb') as stderr:
                    try:
                        result=subprocess.run(record['command'],cwd=work,env=env,stdout=stdout,stderr=stderr,timeout=1000)
                        record.update(exit_code=result.returncode,timed_out=False,status='PASS' if result.returncode==0 else 'FAIL')
                    except subprocess.TimeoutExpired:
                        record.update(exit_code=None,timed_out=True,status='FAIL')
                    except BaseException as exc:
                        record.update(exit_code=None,timed_out=False,status='FAIL',exception=repr(exc));raise
                    finally:
                        record.update(elapsed_seconds=time.monotonic()-started,completed_utc=stamp());write(output/'EXECUTION.json',state['commands'])
                        (logs/(name+'.exit.json')).write_text(json.dumps({k:record.get(k) for k in ('status','exit_code','timed_out','elapsed_seconds')},indent=2)+'\n')
                require(verify_inputs(package,args.freeze_commit)==inputs and verify_inputs(work,args.freeze_commit)==inputs,'Source changed during '+name)
                require(record['status']=='PASS',name+' failed; preserved stdout/stderr/exit and fresh outputs')
                print(name+' PASS',flush=True)
            state.update(status='PASS',counts=final_checks(work,fresh,logs),executed_commands=len(commands),physical_producer_commands_executed=2)
    except BaseException as exc:
        state.update(status='FAIL',failure={'exception':repr(exc),'traceback':traceback.format_exc()});raise
    finally:
        try:
            if inputs is not None:
                require(verify_inputs(package,args.freeze_commit)==inputs,'Original source changed during replay')
                if work.exists():require(verify_inputs(work,args.freeze_commit)==inputs,'Copied source changed during replay')
                state['source_verification_after']='PASS'
            else:state['source_verification_after']='NOT_STARTED'
        except BaseException as exc:
            state.update(status='FAIL',source_verification_after='FAIL',source_verification_failure=repr(exc))
        state.update(completed_utc=stamp(),attempted_commands=len(state['commands']),successful_commands=sum(c['status']=='PASS' for c in state['commands']),physical_producer_commands_attempted=sum(c['physical_producer'] for c in state['commands']),physical_producer_commands_completed=sum(c['physical_producer'] and c['status']=='PASS' for c in state['commands']))
        if fresh.exists():state['fresh_output_sha256']={p.relative_to(fresh).as_posix():sha(p) for p in sorted(fresh.rglob('*')) if p.is_file()}
        write(output/'VALIDATION.json',state)
    require(state['status'] in ('PASS','VERIFIED_PLAN_ONLY'),'Post-run source verification failed')
    print(json.dumps({k:state[k] for k in ('status','planned_commands','attempted_commands','successful_commands')},indent=2))


if __name__=='__main__':main()
