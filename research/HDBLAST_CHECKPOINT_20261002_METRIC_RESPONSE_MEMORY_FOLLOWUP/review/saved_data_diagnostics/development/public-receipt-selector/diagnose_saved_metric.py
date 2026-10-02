#!/usr/bin/env python3
"""Saved metric archives only: honest partial/full diagnostics of failed science.

No physical producer import, source-grid experiment, ODE evolution or response
quadrature. The frozen raw auditor reconstructs saved observations and Simpson
Ward ledgers. Scientific gates are reported without converting FAIL to PASS.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import traceback
import numpy as np

CORE=('q','q_prime','q_second','rho','p')
EXTRAS=('Q0','rho0','p0','current')
NAMES=CORE+EXTRAS

def require(ok,message):
    if not ok:raise RuntimeError(message)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def read(path):
    def unique(items):
        output={}
        for key,value in items:
            require(key not in output,'Duplicate JSON key '+key);output[key]=value
        return output
    return json.loads(Path(path).read_text(),object_pairs_hook=unique,
                      parse_constant=lambda value:(_ for _ in ()).throw(RuntimeError('Nonfinite JSON '+value)))

def write(path,value):Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')

def hexadecimal(value,n):return isinstance(value,str) and len(value)==n and all(c in '0123456789abcdef' for c in value)

def relative(root,name):
    p=Path(name);resolved=(root/p).resolve()
    require(not p.is_absolute() and '..' not in p.parts and resolved.is_relative_to(root),'Input path escape')
    return resolved

def verify_sources(root,pin,freeze,public_path):
    require(hexadecimal(pin,64) and hexadecimal(freeze,40),'Explicit immutable registration/public freeze pins required')
    regpath=root/'FULL_REGISTRATION.json';require(sha(regpath)==pin,'Registration digest')
    registration=read(regpath);require(registration.get('schema_version')==1,'Registration schema')
    files=registration['files'];require(isinstance(files,dict) and files,'Registered source inventory')
    for name,digest in files.items():require(hexadecimal(digest,64) and sha(relative(root,name))==digest,'Frozen input mismatch '+name)
    public=read(public_path)
    require(public.get('status')=='PASS_PUBLIC_TREE_VERIFIED_BEFORE_PHYSICAL_EXECUTION' and public.get('public_freeze_commit')==freeze and public.get('remote_branch_commit')==freeze and public.get('public_repository')=='maldonado-research/HDblast','Faithful actual prospective public receipt required')
    require(public.get('registration_sha256')==pin and public.get('independent_manifest_sha256')==registration['independent_manifest_sha256'],'Public receipt pins')
    require(type(public.get('physical_evaluations_so_far')) is int and public['physical_evaluations_so_far']==0,'Public verification preceded this experiment')
    checks=public['checks'];entries=[r for r in checks if r.get('path','').endswith('/FULL_REGISTRATION.json')]
    require(len(entries)==1 and entries[0]['sha256']==pin,'Public registration entry')
    prefix=entries[0]['path'][:-len('FULL_REGISTRATION.json')]
    expected={prefix+name:(relative(root,name),digest) for name,digest in files.items()}
    expected[prefix+'FULL_REGISTRATION.json']=(regpath,pin)
    require(type(public['files_verified']) is int and public['files_verified']==len(checks)==len(expected),'Public receipt complete coverage')
    seen=set()
    for record in checks:
        name=record['path'];require(name in expected and name not in seen,'Duplicate/unexpected public receipt path');seen.add(name)
        path,digest=expected[name];body=path.read_bytes()
        object_digest=hashlib.sha1(b'blob '+str(len(body)).encode()+b'\0'+body).hexdigest()
        require(record['sha256']==digest and record['git_blob_sha']==object_digest and type(record['bytes']) is int and record['bytes']==len(body),'Public receipt SHA/bytes mismatch')
    require(seen==set(expected),'Missing public receipt path')
    verified=datetime.fromisoformat(public['utc']);require(verified.tzinfo is not None and verified.utcoffset().total_seconds()==0,'UTC public receipt')
    return registration,verified

def load_validator(root):
    sys.path.insert(0,str(root/'code'))
    spec=importlib.util.spec_from_file_location('saved_diagnostic_frozen_validator',root/'code/validate_metric.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module,sys.modules['raw_metric_audit']

def verify_independent(modes,registration,experiment,root,pin,freeze):
    require(modes['schema_version']==1 and modes['route']=='independent_metric_forced_modes_direct_minimal_stress','Independent schema/route')
    require(modes['status'] in ('failed','failed_internal_gates','passed_internal_gates'),'No running producer snapshot')
    require(modes['freeze_commit']==freeze and modes['provenance']['registration_sha256']==pin,'Independent prospective pins')
    require(modes['epsilon']==1e-4 and modes['mass_law_b']==1 and modes['configuration']==registration['frozen_configuration']['independent_configuration'],'Independent source/configuration')
    require(modes['configuration']['state']=='incoming_BD' and modes['configuration']['phi']=='fixed' and modes['configuration']['wronskian_projection'] is False,'State/fixed source/no projection')
    manifest_path=root/'independent/MANIFEST.json';manifest=read(manifest_path)
    require(sha(manifest_path)==registration['independent_manifest_sha256']==modes['provenance']['manifest_sha256'] and modes['provenance']['files']==manifest['files'],'Independent manifest pins')
    seen=set()
    for record in manifest['files']:
        require(record['path'] not in seen,'Duplicate independent source path');seen.add(record['path'])
        require(sha(relative(manifest_path.parent,record['path']))==record['sha256'],'Independent source bytes')
    runtime=modes['runtime'];require(runtime['python'].startswith('3.12.') and runtime['numpy']=='2.2.6' and runtime['python_optimization']==0,'Recorded physical runtime')
    require(np.__version__=='2.2.6' and np.finfo(np.longdouble).nmant>=63,'Diagnostic extended precision/runtime')
    require(runtime['real_dtype']==np.dtype(np.longdouble).name and runtime['complex_dtype']==np.dtype(np.clongdouble).name and runtime['longdouble_eps']==str(np.finfo(np.longdouble).eps),'Recorded dtype precision')
    names={('positive_B','coarse'),('positive_B','fine'),('signed_uB','coarse'),('signed_uB','fine')}
    actual={(r['source'],r['setting']) for r in modes['runs']}
    require(len(actual)==len(modes['runs']) and actual<=names and actual,'Distinct registered saved runs')
    return names,actual

def comparison(source,eta,K,name,left,right,tolerance,label):
    require(all(math.isfinite(float(x)) for x in (left,right,tolerance)) and tolerance>0,'Nonfinite comparison')
    gap=abs(float(left)-float(right))
    return {'source':source,'eta':eta,'K':K,'quantity':name,'left':float(left),'right':float(right),'absolute_gap':gap,'frozen_tolerance':tolerance,'gate_ratio':gap/tolerance,'gate_result':'MEETS_GATE' if gap<=tolerance else 'EXCEEDS_GATE','comparison':label}

def reconstructed_run(raw,run,directory,experiment):
    """Capture pure reconstruction output without changing frozen auditor bytes."""
    path=raw.archive_path(directory,run['archive'])
    with np.load(path,allow_pickle=False) as arc:weights=arc['momentum_weights'].copy()
    original=raw.reconstructed;points={};wronskians=[]
    def capture(eta,k,u,w,hjet,epsilon):
        returned=original(eta,k,u,w,hjet,epsilon)
        data,archive,W,physical,mutations,sub=returned
        wronskians.append({'eta':eta,'linear_scaled':W,'physical_scaled':physical})
        for K in experiment['cutoffs']:
            count=int(np.count_nonzero(k<np.longdouble(K)))
            points[(eta,K)]={name:float(np.sum(data[name][:count]*weights[:count],dtype=np.longdouble)) for name in NAMES}
        return returned
    raw.reconstructed=capture
    try:report=raw.audit_run(run,directory,experiment)
    finally:raw.reconstructed=original
    require(len(points)==18 and len(wronskians)==6,'Complete saved run reconstruction')
    report['saved_observation_wronskians']=wronskians
    report['reconstructed_values']=[dict(eta=eta,K=K,values=value) for (eta,K),value in points.items()]
    return report,points

def plot_report(report,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups=[('Cross route',report['cross_core']+report['cross_baseline_current']),('Coarse/fine',report['core_refinement']+report['current_refinement']),('Fine Ward',report['fine_Ward_gate_comparisons'])]
    fig,axes=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for ax,(label,items) in zip(axes,groups):
        ratios=[max(item['gate_ratio'],1e-16) for item in items]
        ax.scatter(range(len(ratios)),ratios,s=12,color=['#c62828' if r>1 else '#455a64' for r in ratios]);ax.axhline(1,color='#c62828',linestyle='--')
        ax.set_yscale('log');ax.set_xlabel('Saved comparison index');ax.set_ylabel('Gap / unchanged gate');ax.set_title(label)
    fig.suptitle('Underlying science FAIL — saved-data diagnostic ('+report['coverage']['classification']+')')
    fig.savefig(output/'FAILED_SAVED_DATA_GATE_RATIOS.png',dpi=180);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--primary',type=Path,required=True);p.add_argument('--modes',type=Path,required=True)
    p.add_argument('--registration-sha256',required=True);p.add_argument('--public-freeze-commit',required=True);p.add_argument('--public-freeze-evidence',type=Path,required=True)
    p.add_argument('--failure-evidence',type=Path,required=True);p.add_argument('--science-status',choices=['FAIL'],default='FAIL')
    p.add_argument('--output',type=Path,required=True);p.add_argument('--plots',action='store_true');args=p.parse_args()
    root=args.checkpoint.resolve();output=args.output.resolve()
    require(not output.exists() and not output.is_relative_to(root),'Fresh diagnostic output outside checkpoint required');output.mkdir(parents=True)
    report={'diagnostic_status':'STARTED','underlying_science_status':'FAIL','full_metric_prerequisite':'BLOCKED','scientific_validator_PASS':False,'physical_evolution_or_response_quadratures':0,'python_optimization':sys.flags.optimize,'created_at_utc':datetime.now(timezone.utc).isoformat()}
    write(output/'STARTED.json',report)
    try:
        inputs={str(path.resolve()):sha(path) for path in [args.primary,args.modes,args.public_freeze_evidence,args.failure_evidence,root/'FULL_REGISTRATION.json']}
        registration,verified=verify_sources(root,args.registration_sha256,args.public_freeze_commit,args.public_freeze_evidence)
        validator,raw=load_validator(root);experiment=read(root/'EXPERIMENT.json')
        require(registration['frozen_configuration']==experiment and registration['frozen_gates']==experiment['gates'],'Frozen configuration/gates')
        validator.check_experiment(experiment)
        primary,modes,failure=read(args.primary),read(args.modes),read(args.failure_evidence)
        require(failure.get('status') in ('FAIL','failed','failed_internal_gates') and (failure.get('exception') or (type(failure.get('exit_code')) is int and failure['exit_code']!=0) or failure.get('gate_failures')),'Actual retained FAIL evidence required')
        validator.finite_tree(primary,'primary');validator.finite_tree(modes,'modes')
        validator.primary_provenance(primary,args.registration_sha256,args.public_freeze_commit,experiment)
        require(datetime.fromisoformat(primary['provenance']['started_at_utc'])>verified,'Primary execution followed public verification')
        expected,present=verify_independent(modes,registration,experiment,root,args.registration_sha256,args.public_freeze_commit)
        pr=validator.row_index(primary,experiment,'primary')
        audits=[];errors=[];values={};archives={}
        for run in modes['runs']:
            ident=(run['source'],run['setting']);path=raw.archive_path(args.modes.parent,run['archive']);archives[str(path)]=sha(path)
            try:
                audit,points=reconstructed_run(raw,run,args.modes.parent,experiment);audits.append(audit);values[ident]=points
            except Exception as exc:errors.append({'source':ident[0],'setting':ident[1],'archive_sha256':sha(path),'audit_error':repr(exc),'traceback':traceback.format_exc(),'reconstructed_values_accepted':False})
            write(output/'PROGRESS.json',dict(report,audited_runs=len(audits),raw_audit_errors=errors))
        g=experiment['gates'];cross=[];extras=[];refine=[];current=[];baseline_refine=[];ward=[];fineward=[];wardref=[]
        for audit in audits:
            for point in audit['ward_endpoints']:
                entry=dict(point,source=audit['source'],setting=audit['setting']);ward.append(entry)
                if audit['setting']=='fine':fineward.append(comparison(audit['source'],point['eta'],point['K'],'Ward_endpoint',point['residual'],0.0,g['ward_endpoint'],'reconstructed saved-history fine endpoint'))
        for source in experiment['source']['ids']:
            fine=values.get((source,'fine'));coarse=values.get((source,'coarse'))
            if fine:
                for (eta,K),point in fine.items():
                    primary_point=next(x for x in pr[(source,eta)]['finite_k'] if x['K']==K)['values']
                    for name in CORE:cross.append(comparison(source,eta,K,name,primary_point[name],point[name],g['cross_'+name],'primary vs reconstructed fine saved modes'))
                    for name in EXTRAS:extras.append(comparison(source,eta,K,name,primary_point[name],point[name],g['cross_'+name],'primary vs reconstructed fine saved modes'))
            if coarse and fine:
                for (eta,K),point in fine.items():
                    for name in CORE:refine.append(comparison(source,eta,K,name,point[name],coarse[(eta,K)][name],g['refinement_'+name],'reconstructed fine vs coarse'))
                    current.append(comparison(source,eta,K,'current',point['current'],coarse[(eta,K)]['current'],g['independent']['refinement']['current'],'reconstructed fine vs coarse'))
                    for name in ('Q0','rho0','p0'):baseline_refine.append({'source':source,'eta':eta,'K':K,'quantity':name,'absolute_difference':abs(point[name]-coarse[(eta,K)][name]),'gate':'none registered for baseline refinement'})
                wc={(x['eta'],x['K']):x for x in ward if x['source']==source and x['setting']=='coarse'}
                wf={(x['eta'],x['K']):x for x in ward if x['source']==source and x['setting']=='fine'}
                for (eta,K),point in wf.items():wardref.append(comparison(source,eta,K,'Ward_ledger',point['ledger'],wc[(eta,K)]['ledger'],g['ward_refinement'],'reconstructed saved-history fine vs coarse Simpson ledger'))
        grouped={'cross_core':cross,'cross_baseline_current':extras,'core_refinement':refine,'current_refinement':current,'fine_Ward_gate_comparisons':fineward,'Ward_refinement':wardref}
        stats={name:{'count':len(items),'exceeds_gate':sum(x['gate_result']=='EXCEEDS_GATE' for x in items),'maximum_gate_ratio':max((x['gate_ratio'] for x in items),default=None),'maximum_absolute_gap':max((x['absolute_gap'] for x in items),default=None)} for name,items in grouped.items()}
        resources=modes.get('resources',{})
        report.update(diagnostic_status='COMPLETED_SAVED_DATA_DIAGNOSTIC',**grouped,raw_archive_audits=audits,raw_audit_errors=errors,reconstructed_Ward_endpoints=ward,baseline_refinement_descriptive=baseline_refine,gate_statistics=stats,
          coverage={'classification':'complete four-run saved coverage' if present==expected and not errors else 'partial saved coverage','present_runs':[list(x) for x in sorted(present)],'missing_runs':[list(x) for x in sorted(expected-present)],'audited_runs':len(audits),'raw_reconstructed_points':sum(x['reconstructed_observation_points'] for x in audits),'raw_Ward_endpoints':len(ward),'expected_complete_raw_points':72,'expected_complete_raw_Ward_endpoints':72,'expected_complete_cross_core':180,'expected_complete_cross_baseline_current':144,'expected_complete_core_refinement':180,'expected_complete_current_refinement':36,'expected_complete_Ward_refinement':36},
          recorded_producer_status={'primary':primary['status'],'independent':modes['status'],'independent_exception':modes.get('exception'),'independent_declared_gate_failures':modes.get('gate_failures',[]),'independent_resources':resources,'independent_elapsed_within_registered_budget':resources.get('elapsed_seconds',math.inf)<=experiment['independent']['producer_wall_seconds'],'independent_peak_rss_within_registered_budget':resources.get('peak_rss_kib',math.inf)<=experiment['independent']['peak_rss_kib']},
          provenance={'diagnostic_source_sha256':sha(__file__),'registration_sha256':args.registration_sha256,'public_freeze_commit':args.public_freeze_commit,'registered_files_verified':len(registration['files']),'raw_auditor_sha256':sha(root/'code/raw_metric_audit.py'),'frozen_validator_sha256':sha(root/'code/validate_metric.py'),'input_sha256':inputs,'archive_sha256':archives},
          scope='Read-only reconstruction from saved complex observation modes and direct saved histories; no physical producer import/evolution/response quadrature. The source-pinned raw auditor enforces its original structural/operator/subtraction/Wronskian checks. Its returned Ward endpoints and unchanged cross/refinement gates are diagnostic classifications only. Scientific FAIL, original failures and full-metric BLOCKED remain unchanged. Global all-time Wronskians remain producer declarations; raw observation Wronskians are reconstructed. No continuum tail/precision certificate, complete original-validator PASS, stability or heating claim.')
        for name,digest in {**inputs,**archives}.items():require(sha(name)==digest,'Saved input changed during diagnosis')
        write(output/'FAILED_SAVED_DATA_DIAGNOSTIC.json',report)
        if args.plots:plot_report(report,output)
        lines=['# Saved metric diagnostic: underlying science FAIL','',report['scope'],'',f'Coverage: {report["coverage"]["classification"]}; {len(audits)} audited runs, {report["coverage"]["raw_reconstructed_points"]} raw points and {len(ward)} Ward endpoints.','', '| Comparison | Count | Exceeds unchanged gate | Maximum gate ratio |','|---|---:|---:|---:|']
        for name,stat in stats.items():lines.append(f'| {name} | {stat["count"]} | {stat["exceeds_gate"]} | {stat["maximum_gate_ratio"]} |')
        lines+=['','Individual comparison rows and all raw audit failures are retained in the JSON. Coarse Ward endpoints and baseline refinement differences have no extra registered gate. A successful diagnostic process does not imply scientific validation.']
        (output/'FAILED_SAVED_DATA_DIAGNOSTIC.md').write_text('\n'.join(lines)+'\n')
        print(json.dumps({'diagnostic_status':report['diagnostic_status'],'underlying_science_status':'FAIL','coverage':report['coverage'],'gate_statistics':stats},indent=2))
        return 0 if not errors else 2
    except Exception as exc:
        report.update(diagnostic_status='BLOCKED_SAVED_DATA_DIAGNOSTIC',exception=repr(exc),traceback=traceback.format_exc());write(output/'DIAGNOSTIC_FAILURE.json',report);raise

if __name__=='__main__':sys.exit(main())
