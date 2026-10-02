#!/usr/bin/env python3
"""Self-contained replay of the publicly frozen matched-stress checkpoint.

Install pinned requirements with Python 3.12. Put this driver in checkpoint/code
and support inputs in checkpoint/replay_inputs. No repository path is needed.
Each physical producer runs once, normally; validators, six analytic/format
checks and the separate saved-data audit/future analytic proof run in normal
and optimized modes. Saved-data reports/figures add no physical evaluations.
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
POSTRUN={'postrun/summarize_stress.py':'412508e4506d5454e906806271c0563279ccc70f48e1d075c6396337bf9c6857','postrun/plot_stress.py':'bdfdf9ad6861f607c363f48cde3663f5a7ceda648a4a063a345d3308eee0673a'}
AUDIT_SOURCE='94f1d3f1d6197c5dd9b1865889daa69a07c90d44d48de399c90b0fd85a0a78d9'
FUTURE_PROOF_SOURCE='d40c08f12c3d9c7f287e7b6940bd4f7bf5137857c314f1d039a1facc35c1b780'
FUTURE_SOURCE_PINS='f7168f15e63dc60d40173bf0a1f2d4010c35823ecfbf560f1a6c9f68112b6beb'
AUDIT_COUNTS={'mode_snapshots':24,'reconstructed_quantity_sums':576,'registered_ward_endpoints':72,'ward_refinements':36,'panel_polynomial_moments':98304}
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
    for name,pin in POSTRUN.items():require(digest(safe_path(package,name))==pin,'Postrun presentation source changed: '+name)
    for name,pin in SUPPORT.items():require(digest(safe_path(package,'replay_inputs/'+name))==pin,'Standalone replay support changed: '+name)
    pins=read(package/'theory/SOURCE_PINS.json')['files']
    require({p['repository_relative_path'] for p in pins}==set(SOURCE_MAP),'Tail algebra source-pin membership changed')
    for item in pins:require(digest(safe_path(package,SOURCE_MAP[item['repository_relative_path']]))==item['sha256'],'Tail source-pin input changed')
    require(digest(safe_path(package,'review/archive_audit/audit_saved_stress.py'))==AUDIT_SOURCE,'Postrun independent audit source changed')
    require(digest(safe_path(package,'theory/second_order/verify_second_order_energy_algebra.py'))==FUTURE_PROOF_SOURCE,'Future analytic proof source changed')
    require(digest(safe_path(package,'theory/second_order/SOURCE_PINS.json'))==FUTURE_SOURCE_PINS,'Future analytic source pins changed')
    future_pins=read(package/'theory/second_order/SOURCE_PINS.json')
    require(future_pins['status']=='INHERITED_INPUTS_FOR_UNREGISTERED_ANALYTIC_FOLLOWUP' and future_pins['preceding_stress_freeze']==FREEZE and future_pins['no_new_physical_evaluation'] and len(future_pins['sources'])==9,'Future analytic input classification/count')
    for item in future_pins['sources']:require(digest(safe_path(package,'theory/second_order/'+item['local_copy']))==item['sha256'],'Future analytic reference input changed')
    payload=None
    if (package/'MANIFEST.json').exists():
        payload=read(package/'MANIFEST.json');excluded=set(payload['excluded'])
        actual={p.relative_to(package).as_posix() for p in package.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.relative_to(package).as_posix() not in excluded}
        require(actual==set(payload['files']),'Packaged payload membership changed')
        for name,pin in payload['files'].items():require(digest(safe_path(package,name))==pin,'Payload changed: '+name)
    return {'registration_sha256':REGISTRATION,'independent_manifest_sha256':INDEPENDENT,
      'frozen_file_count':len(frozen),'frozen_files':frozen,'support_files':SUPPORT,'postrun_sources':POSTRUN,'postrun_saved_mode_audit_source':AUDIT_SOURCE,'future_analytic_proof_source':FUTURE_PROOF_SOURCE,'future_analytic_source_pins':FUTURE_SOURCE_PINS,
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
    commands.extend([('summary',['postrun/summarize_stress.py','--checkpoint',str(work),'--primary',str(fresh/'primary/results.json'),'--modes',str(fresh/'independent/results.json'),'--checks',str(fresh/'validation/CHECKS.json'),'--checks-optimized',str(fresh/'validation_optimized/CHECKS.json'),'--output',str(fresh/'report')]),('figures',['postrun/plot_stress.py','--summary',str(fresh/'report/SUMMARY.json'),'--output',str(fresh/'figures')])])
    for optimized in (False,True):
        suffix='_optimized' if optimized else '';flag=['-O'] if optimized else []
        commands.extend([('postrun_stress_audit'+suffix,[*flag,'review/archive_audit/audit_saved_stress.py','--checkpoint',str(work),'--modes',str(fresh/'independent/results.json'),'--output',str(fresh/('POSTRUN_STRESS_AUDIT'+suffix+'.json'))]),('second_order_analytic_proof'+suffix,[*flag,'theory/second_order/verify_second_order_energy_algebra.py','--output',str(fresh/('SECOND_ORDER_ANALYTIC_PROOF'+suffix+'.json'))])])
    require(len(commands)==22,'Unexpected registered, postrun and separate analytic command plan')
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

def presentation_checks(work,fresh):
    summary=read(fresh/'report/SUMMARY.json');primary=read(fresh/'primary/results.json');modes=read(fresh/'independent/results.json')
    checks=read(fresh/'validation/CHECKS.json')
    require(summary['status']=='registered_matched_stress_calibration_passed' and summary['postrun_only'] and summary['new_physical_source_response_or_mode_evaluations']==0,'Postrun summary status/scope')
    require(summary['public_freeze_commit']==FREEZE and summary['registration_sha256']==REGISTRATION,'Summary prospective pins')
    require(summary['frozen_source_sha256']==read(work/'FULL_REGISTRATION.json')['files'],'Summary frozen-source pins')
    require(summary['normalizations']==primary['normalization'] and summary['epsilon']==primary['epsilon'],'Summary units/amplitude')
    require(summary['metrics']['counts']==EXPECTED_COUNTS and summary['controls']==checks['controls'] and summary['Ward_endpoints']==checks['raw_direct_history_Ward_endpoints'] and summary['trace_comparisons']==checks['finite_K_trace_comparisons'] and summary['core_cross_route_comparisons']==checks['core_comparisons'],'Summary counts and evidence')
    paths={'FULL_REGISTRATION.json':work/'FULL_REGISTRATION.json','outputs/primary/results.json':fresh/'primary/results.json','outputs/independent/results.json':fresh/'independent/results.json','outputs/validation/CHECKS.json':fresh/'validation/CHECKS.json','outputs/validation_optimized/CHECKS.json':fresh/'validation_optimized/CHECKS.json','postrun/summarize_stress.py':work/'postrun/summarize_stress.py'}
    for run in modes['runs']:paths['outputs/independent/'+run['archive']['path']]=fresh/'independent'/run['archive']['path']
    for name in ('started.json','partial_results.json','EXECUTION.json'):paths['outputs/primary/'+name]=fresh/'primary'/name
    require(summary['raw_and_presentation_sha256']=={name:digest(path) for name,path in paths.items()},'Summary raw/presentation hashes')
    indexed={(r['source'],r['eta']):r for r in summary['registered_rows']}
    require(len(indexed)==len(summary['registered_rows'])==12,'Summary registered point coverage')
    mode_rows={(r['source'],r['eta']):r for r in modes['rows']}
    for row in primary['rows']:
        presented=indexed[(row['source'],row['eta'])];direct=mode_rows[(row['source'],row['eta'])]
        require(presented['a']==row['a'] and presented['source_jet_over_epsilon']==row['source_jet_over_epsilon'] and presented['continuum_values']==row['continuum']['values'] and presented['continuum_components']==row['continuum']['components'] and presented['continuum_quadrature_error_estimates']==row['continuum']['quadrature_error_estimates'],'Summary continuum/source values')
        a=row['a'];epsilon=primary['epsilon'];v=row['continuum']['values']
        physical={'delta_Q':epsilon*v['q']/a**2,'delta_rho':epsilon*v['rho']/a**4,'delta_p':epsilon*v['p']/a**4,'delta_j':epsilon*v['current']/a**2}
        require(presented['physical_continuum_response']==physical,'Summary physical normalization conversion')
        require(len(presented['finite_k'])==3,'Summary cutoff count')
        for point,finite,mode in zip(presented['finite_k'],row['finite_k'],direct['finite_k']):
            require(point['K']==finite['K']==mode['K'] and point['primary_values']==finite['values'] and point['direct_mode_values']==mode['values'],'Summary finite-band raw values')
            require(point['quadrature_error_estimates']==finite['quadrature_error_estimates'] and point['mode_refinement_differences']==mode['refinement_differences'] and point['ward']==mode['ward'],'Summary uncertainty and Ward values')
            names=('q','q_prime','q_second','rho','p','Q0','anomaly','current')
            require(point['analytic_omitted_band_bounds']=={n:finite['analytic_tail_bounds']['Q0_difference' if n=='Q0' else n] for n in names},'Summary actual UV bounds')
            for label,first,second in (('primary_to_continuum_differences',finite['values'],v),('direct_to_continuum_differences',mode['values'],v),('direct_to_primary_differences',mode['values'],finite['values'])):
                require(point[label]=={n:abs(first[n]-second[n]) for n in names},'Summary actual differences')
    metrics=summary['metrics']
    for metric,entries in (('maximum_core_cross_route_difference',checks['core_comparisons']),('maximum_extra_cross_route_difference',checks['baseline_anomaly_current_comparisons']),('maximum_direct_stress_closed_difference',checks['direct_stress_closed_comparisons'])):
        expected={name:max(x['gap'] for x in entries if x['quantity']==name) for name in {x['quantity'] for x in entries}}
        require(metrics[metric]==expected,'Summary actual cross-route maxima')
    names=('q','q_prime','q_second','rho','p','Q0','anomaly','current');final=[r['finite_k'][-1] for r in summary['registered_rows']]
    for metric,label in (('maximum_K256_direct_to_continuum_difference','direct_to_continuum_differences'),('maximum_K256_primary_to_continuum_difference','primary_to_continuum_differences'),('maximum_K256_analytic_tail_bound','analytic_omitted_band_bounds')):
        require(metrics[metric]=={n:max(point[label][n] for point in final) for n in names},'Summary K256 maxima')
    quad={n:max([r['continuum']['quadrature_error_estimates'][n] for r in primary['rows']]+[p['quadrature_error_estimates'][n] for r in primary['rows'] for p in r['finite_k']]) for n in names}
    require(metrics['maximum_primary_quadrature_estimate']==quad and metrics['maximum_mode_refinement_difference']==modes['maximum_refinement_differences'],'Summary numerical uncertainty maxima')
    values={'maximum_fine_Ward_endpoint_difference':modes['maximum_ward_endpoint_difference'],'maximum_Ward_refinement_difference':modes['maximum_ward_refinement_difference'],'maximum_finite_K_trace_difference':max(x['gap'] for x in checks['finite_K_trace_comparisons']),'maximum_linear_Wronskian_over_epsilon':max(r['wronskian_max_scaled'] for r in modes['runs']),'maximum_physical_Wronskian_over_epsilon':max(r['physical_wronskian_max_scaled'] for r in modes['runs']),'primary_elapsed_seconds':primary['elapsed_seconds'],'independent_elapsed_seconds':modes['resources']['elapsed_seconds'],'independent_peak_rss_kib':modes['resources']['peak_rss_kib'],'algebraic_or_saved_operator_sensitivity_controls':sum(c.get('nonzero_algebraic_sensitivity_witness',False) for c in checks['controls']),'sensitivity_controls_exceeding_operational_cross_gate':sum(c.get('operational_gate_exceeded_somewhere',False) for c in checks['controls']),'synthetic_guard_mutations_rejected':sum(c.get('rejected',False) for c in checks['controls'])}
    require(all(metrics[name]==value for name,value in values.items()),'Summary recorded residual/resource/control metrics')
    figure_directory=fresh/'figures';provenance=read(figure_directory/'FIGURE_PROVENANCE.json')
    expected={name+'.'+extension for name in ('stress_response','stress_uv_bounds','stress_checks') for extension in ('svg','png','pdf')}
    actual={p.name for p in figure_directory.iterdir() if p.is_file() and p.suffix in ('.svg','.png','.pdf')}
    require(actual==expected and set(provenance['figures_sha256'])==expected,'Exactly nine figure artifacts required')
    require(provenance['summary_sha256']==digest(fresh/'report/SUMMARY.json') and provenance['renderer_sha256']==digest(work/'postrun/plot_stress.py') and provenance['matplotlib']=='3.10.1' and provenance['postrun_only'] and provenance['new_physical_evaluations']==0,'Figure source/scope/version provenance')
    for name,pin in provenance['figures_sha256'].items():require(digest(figure_directory/name)==pin,'Figure output hash mismatch: '+name)
    require((fresh/'report/NUMERICAL_RESULTS.md').is_file(),'Numerical report missing')
    return {'summary_registered_rows':12,'summary_finite_band_points':36,'summary_raw_presentation_hashes':len(paths),'figure_files':len(expected),'matplotlib':'3.10.1'}

def additional_postrun_checks(work,fresh):
    modes=read(fresh/'independent/results.json');runs={(r['source'],r['setting']):r for r in modes['runs']}
    reports=[]
    for suffix,opt in (('',0),('_optimized',1)):
        report=read(fresh/('POSTRUN_STRESS_AUDIT'+suffix+'.json'));reports.append(report)
        require(report['status']=='passed' and report['classification']=='post_run_saved_data_audit' and report['counts']==AUDIT_COUNTS,'Second saved-data audit status/counts')
        runtime=report['runtime'];require(runtime['python_optimization']==opt and runtime['python'].startswith('3.12.') and runtime['numpy']=='2.2.6' and runtime['longdouble_nmant']>=63,'Second audit runtime/optimization')
        require(report['provenance']=={'freeze_commit':FREEZE,'registration_sha256':REGISTRATION,'mode_results_sha256':digest(fresh/'independent/results.json'),'audit_source_sha256':AUDIT_SOURCE},'Second audit source/input provenance')
        require(len(report['runs'])==4 and {(r['source'],r['setting']) for r in report['runs']}==set(runs),'Second audit run coverage')
        for audit in report['runs']:
            original=runs[(audit['source'],audit['setting'])]
            require(audit['archive_sha256']==original['archive']['sha256']==digest(fresh/'independent'/original['archive']['path']),'Second audit raw archive pin')
            require(len(audit['snapshots'])==6 and len(audit['registered_ward_endpoints'])==18,'Second audit snapshot/endpoint coverage')
            for snapshot,row in zip(audit['snapshots'],original['rows']):
                require(snapshot['eta']==row['eta'] and len(snapshot['finite_k'])==3,'Second audit saved observation grid')
                for point,raw in zip(snapshot['finite_k'],row['finite_k']):
                    require(point['K']==raw['K'] and set(point['values'])==set(raw['values']),'Second audit cutoff/quantity grid')
                    for name,value in point['values'].items():require(abs(value-raw['values'][name])<=2e-11+2e-13*max(1,abs(raw['values'][name])),'Second audit actual reconstructed value '+name)
        require(len(report['ward_refinements'])==36,'Second audit refinement point coverage')
    normal={k:v for k,v in reports[0].items() if k!='runtime'};optimized={k:v for k,v in reports[1].items() if k!='runtime'}
    require(normal==optimized,'Second audit normal/optimized numerical evidence differs')
    require({k:v for k,v in reports[0]['runtime'].items() if k!='python_optimization'}=={k:v for k,v in reports[1]['runtime'].items() if k!='python_optimization'},'Second audit normal/optimized runtime differs')
    proofs=[]
    for suffix in ('','_optimized'):
        proof=read(fresh/('SECOND_ORDER_ANALYTIC_PROOF'+suffix+'.json'));proofs.append(proof)
        require(proof['status']=='PASS' and proof['proposal_status']=='UNREGISTERED_AND_NUMERICALLY_UNEXECUTED','Future analytic proof incorrectly classified')
        require(proof['checks_count']==len(proof['checks'])==len(set(proof['checks']))==21 and proof['mutations_rejected_count']==len(proof['mutations_rejected'])==len(set(proof['mutations_rejected']))==8 and proof['source_pins_checked']==9,'Future analytic identity/mutation/source count')
        require(proof['verifier_sha256']==FUTURE_PROOF_SOURCE and proof['sympy']=='1.14.0' and proof['python'].startswith('3.12.') and proof['checkout_sources_checked'] is False,'Future self-contained analytic proof provenance/runtime')
    require(proofs[0]==proofs[1],'Future exact proof normal/optimized evidence differs')
    return {'saved_mode_audit':{'classification':'post_run_saved_data_audit','counts':AUDIT_COUNTS,'normal_and_optimized':'PASS'},'future_second_order_proposal':{'classification':'UNREGISTERED_AND_NUMERICALLY_UNEXECUTED','physical_evaluations':0,'exact_identities':21,'algebraic_mutations_rejected':8,'local_source_pins':9,'normal_and_optimized':'PASS','optimization_modes_recorded_in':'EXECUTION.json commands; the exact proof JSON does not include an optimization flag'}}

def final_checks(work,fresh,logs):
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
    return {'response_points':12,'independent_runs':4,'validator':EXPECTED_COUNTS,'analytic_checks_by_mode':totals,'presentation':presentation_checks(work,fresh),'additional_postrun':additional_postrun_checks(work,fresh)}

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
        for name,version in {'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}.items():require(importlib.metadata.version(name)==version,'Pinned dependency mismatch: '+name)
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
            state.update(status='PASS',counts=final_checks(work,fresh,logs),executed_commands=len(state['commands']),physical_producer_commands_executed=2)
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
