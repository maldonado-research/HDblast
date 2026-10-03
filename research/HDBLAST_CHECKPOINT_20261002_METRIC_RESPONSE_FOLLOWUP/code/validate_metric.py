#!/usr/bin/env python3
"""Frozen metric-response evidence validator. No physical producer imports.

The primary saved memory coefficients are independently reduced. All saved
complex observation modes, operator/contact arrays, weighted finite bands and
direct-history Ward ledgers are reconstructed separately. Pure action/contact
algebra is source-pinned; it is not an independent numerical producer.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import traceback
import numpy as np
from raw_metric_audit import NAMES,LD,audit_run,jets,forcing,near_array,near_json,require

ROOT=Path(__file__).resolve().parents[1]
CORE=('q','q_prime','q_second','rho','p')
EXTRAS=('Q0','rho0','p0','current')
CONTACT_MODULE=None
TAIL_MODULE=None
EXPECTED_GATES={'cross_q':2e-9,'cross_q_prime':2e-9,'cross_q_second':2e-8,
 'cross_rho':2e-7,'cross_p':2e-7,'cross_Q0':2e-11,'cross_rho0':2e-11,
 'cross_p0':2e-11,'cross_current':2e-9,'refinement_q':2e-10,
 'refinement_q_prime':2e-10,'refinement_q_second':1e-8,'refinement_rho':1e-8,
 'refinement_p':1e-8,'ward_endpoint':2e-6,'ward_refinement':1e-6,
 'wronskian_over_epsilon':1e-10,'numerical_roundoff_allowance':1e-9,
 'primary_normalized_quad_estimate_max':1e-9}

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path,data): Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')

def close(actual,expected,tolerance,label):
    require(math.isfinite(actual) and math.isfinite(expected) and math.isfinite(tolerance) and tolerance>=0,label+': nonfinite/tolerance')
    gap=abs(actual-expected);require(gap<=tolerance,label+': mismatch; residual='+str(gap)+' tolerance='+str(tolerance))
    return gap

def hexadecimal(value,length,label):
    require(isinstance(value,str) and len(value)==length and all(c in '0123456789abcdef' for c in value),label+': full hexadecimal pin')

def finite_tree(value,label):
    if isinstance(value,dict):
        for key,item in value.items(): finite_tree(item,label+'/'+str(key))
    elif isinstance(value,list):
        for i,item in enumerate(value): finite_tree(item,label+'/'+str(i))
    elif isinstance(value,(int,float)) and not isinstance(value,bool): require(math.isfinite(value),label+': nonfinite')

def check_registration(pin):
    hexadecimal(pin,64,'Registration SHA');path=ROOT/'FULL_REGISTRATION.json'
    require(sha(path)==pin,'Registration hash mismatch');registration=json.loads(path.read_text())
    for name,digest in registration['files'].items():
        relative=Path(name);target=(ROOT/relative).resolve()
        require(not relative.is_absolute() and '..' not in relative.parts and target.is_relative_to(ROOT),'Registration path escape')
        hexadecimal(digest,64,'Registered file digest');require(sha(target)==digest,'Frozen file changed: '+name)
    needed={'EXPERIMENT.json','code/metric_primary.py','code/source_jet.py','code/metric_contact_coefficients.py','code/validate_metric.py','code/raw_metric_audit.py','code/test_metric_validator_guards.py','code/verify_raw_metric_baselines.py','independent/forced_metric.py','independent/MANIFEST.json','theory/metric_wkb.py','theory/metric_tail_bounds.py','theory/DERIVATIVE_CERTIFICATE.json','theory/METRIC_TAIL_INPUTS.json'}
    require(needed<=set(registration['files']),'Required registered producer/validator/pure-theory input absent')
    experiment=json.loads((ROOT/'EXPERIMENT.json').read_text())
    require(registration['frozen_configuration']==experiment and registration['frozen_gates']==experiment['gates'],'Registration configuration/gates')
    require(registration['independent_manifest_sha256']==sha(ROOT/'independent/MANIFEST.json'),'Registration independent manifest pin')
    return registration,experiment

def check_experiment(experiment):
    geometry=experiment['geometry'];require(geometry['H']==1 and geometry['r']==2 and geometry['xi']==0 and geometry['geometry_response'] is True,'Registered metric geometry')
    require(experiment['initial_eta']==-6 and experiment['final_eta']==-1.5,'Registered domain')
    require(experiment['source']['ids']==['positive_B','signed_uB'] and experiment['source']['support']==[-5,-3] and experiment['source']['epsilon']==1e-4 and experiment['source']['mass_law_translation_b']==1,'Registered source/current law')
    require(experiment['observation_eta']==[-5.5,-4.5,-4,-3.5,-2.5,-1.5] and experiment['cutoffs']==[64,128,256],'Registered observations/cutoffs')
    for name,value in EXPECTED_GATES.items(): require(experiment['gates'].get(name)==value,'Prospective gate changed: '+name)
    require(experiment['primary']['producer_wall_seconds']==900 and experiment['independent']['producer_wall_seconds']==900 and experiment['independent']['peak_rss_kib']==256*1024,'Registered resource budgets')
    require(experiment['independent']['settings']=={'coarse':{'dt':'1/128','momentum_panel_width':'1/2'},'fine':{'dt':'1/256','momentum_panel_width':'1/4'}},'Registered independent resolution settings')
    require(experiment['independent']['momentum_order']==16 and experiment['independent']['time_order']==8,'Registered Gauss rules')

def row_index(document,experiment,label):
    rows=document['rows'];require(isinstance(rows,list),label+': rows list')
    expected={(s,t) for s in experiment['source']['ids'] for t in experiment['observation_eta']}
    indexed={(r['source'],r['eta']):r for r in rows}
    require(len(indexed)==len(rows) and set(indexed)==expected,label+': source/time grid or duplicate')
    for row in indexed.values():
        require([p['K'] for p in row['finite_k']]==experiment['cutoffs'],label+': cutoff grid')
        close(row['a'],-1/row['eta'],1e-15,label+': scale factor')
        declared=np.array([row['source_jet_over_epsilon'][f'h{n}'] for n in range(6)],dtype=LD)
        expected_jet=jets(np.array([row['eta']],dtype=LD),row['source'])[:,0]
        near_json(declared,expected_jet,label+': source derivatives')
        name='canonical_forcing_jet_over_epsilon' if label=='primary' else 'forcing_jet_over_epsilon'
        near_json([row[name][f'g{n}'] for n in range(4)],forcing(LD(str(row['eta'])),expected_jet),label+': canonical forcing derivatives')
    return indexed

def contact_coefficients(L,v,hjet):
    global CONTACT_MODULE
    if CONTACT_MODULE is None:
        spec=importlib.util.spec_from_file_location('registered_metric_contacts',ROOT/'code/metric_contact_coefficients.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);CONTACT_MODULE=module
    return CONTACT_MODULE.local_coefficients(L,v,hjet)

def metric_tail(source,eta,K):
    global TAIL_MODULE
    if TAIL_MODULE is None:
        spec=importlib.util.spec_from_file_location('registered_metric_tail',ROOT/'theory/metric_tail_bounds.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);TAIL_MODULE=module
    return TAIL_MODULE.metric_tail_bounds(source,eta,K)

def reduced_block(eta,hjet,gjet,memory,cutoff=None):
    """Saved primary history integrals plus action-pinned local algebra only."""
    L=a=-1/eta;v=1.0 if cutoff is None else cutoff/math.hypot(cutoff,math.sqrt(2)*a)
    g,g1,g2=gjet[:3];pref=1/(8*math.pi**2)
    if cutoff is None:
        q,qp,qpp=-pref*memory[0],-pref*(memory[1]+L*g),-pref*(memory[2]+2*L*g1+L*L*g)
    else:
        A=math.asinh(cutoff/(math.sqrt(2)*a))-v;Ap=-L*v**3;App=L*L*v**3*(2-3*v*v)
        q,qp,qpp=pref*(g*A-memory[0]),pref*(g1*A+g*Ap-memory[1]),pref*(g2*A+2*g1*Ap+g*App-memory[2])
    c=contact_coefficients(L,v,hjet);Q=c['Q0']
    # Independently written compact baseline variance; avoids defining it by c.
    close(Q,(v*v/(2*(1+v))-v**3/48-v**5/16)/(2*math.pi**2),1e-14,'Baseline variance primitive')
    j5,j7,j9=v**3/6,(v**3/3-v**5/5)/4,(v**3/3-2*v**5/5+v**7/7)/8
    cr=(3*L*L*g-L*g1)*v**3/(96*math.pi**2)
    cp=(70*L*L*g*j9-(30*L*L*g+10*L*g1)*j7+(g2-L*g1-9*L*L*g)*j5)/(48*math.pi**2)
    scalarR=(3*L*L*q-L*qp)/2+a*a*Q*g/2+cr
    scalarP=(qpp-3*L*qp-3*L*L*q)/6-a*a*Q*g/6+cp
    values={'q':q+c['q_local'],'q_prime':qp+c['q_local_prime'],'q_second':qpp+c['q_local_second'],
            'rho':scalarR+c['rho_local'],'p':scalarP+c['p_local'],'Q0':Q,'rho0':c['rho0'],'p0':c['p0'],'current':q+c['q_local']}
    Qp=-L*v*(1-v*v)*(v*(2+v)/(2*(1+v)**2)-v*v/16-5*v**4/16)/(2*math.pi**2)
    crp=((6*L**3*g+2*L*L*g1-L*g2)*v**3-3*L*(3*L*L*g-L*g1)*v**3*(1-v*v))/(96*math.pi**2)
    Rp=3*L**3*q+L*L*qp-L*qpp/2+a*a*(L*Q*g+Qp*g/2+Q*g1/2)+crp+c['rho_local_prime']
    components={'scalar_memory_q':q,'scalar_memory_q_prime':qp,'scalar_memory_q_second':qpp,'scalar_mass_density':scalarR,'scalar_mass_pressure':scalarP,'full_metric_local_variance':c['q_local'],'full_metric_local_density':c['rho_local'],'full_metric_local_pressure':c['p_local']}
    return values,components,dict(c,v=v),Rp

def primary_provenance(primary,registration_pin,freeze,experiment):
    require(primary['status']=='completed' and primary['schema_version']==1 and primary['route']=='scalar_canonical_memory_and_independently_integrated_rational_full_metric_contacts','Primary route/status')
    p=primary['provenance'];require(p['registration_sha256']==registration_pin and p['public_freeze_commit']==freeze,'Primary prospective pins')
    for field,path in (('producer_sha256','code/metric_primary.py'),('source_jet_sha256','code/source_jet.py'),('contact_coefficients_sha256','code/metric_contact_coefficients.py'),('experiment_sha256','EXPERIMENT.json')):
        require(p[field]==sha(ROOT/path),'Primary input digest '+field)
    require(p['dependencies']=={'numpy':'2.2.6','scipy':'1.15.3','mpmath':'1.3.0','sympy':'1.14.0'} and p['python'].startswith('3.12.'),'Primary runtime dependencies')
    require(primary['normalization']==experiment['normalization'] and primary['epsilon']==1e-4 and primary['mass_law_b']==1,'Primary units/amplitude/current law')
    require(0<=primary['elapsed_seconds']<=experiment['primary']['producer_wall_seconds'],'Primary wall budget')

def independent_provenance(modes,registration_pin,freeze,experiment,registration):
    require(modes['status']=='passed_internal_gates' and modes['schema_version']==1 and modes['route']=='independent_metric_forced_modes_direct_minimal_stress' and modes['gate_failures']==[],'Independent route/status')
    require(modes['freeze_commit']==freeze and modes['epsilon']==1e-4 and modes['mass_law_b']==1,'Independent prospective pin/amplitude/current')
    require(modes['provenance']['registration_sha256']==registration_pin and modes['provenance']['manifest_sha256']==registration['independent_manifest_sha256'],'Independent prospective registration/manifest')
    manifest_path=ROOT/'independent/MANIFEST.json';manifest=json.loads(manifest_path.read_text())
    require(sha(manifest_path)==registration['independent_manifest_sha256'] and modes['provenance']['files']==manifest['files'],'Independent input manifest pin')
    seen=set()
    for entry in manifest['files']:
        relative=Path(entry['path']);path=(manifest_path.parent/relative).resolve()
        require(not relative.is_absolute() and '..' not in relative.parts and path.is_relative_to(manifest_path.parent) and entry['path'] not in seen,'Independent manifest path/duplicate')
        seen.add(entry['path']);require(sha(path)==entry['sha256'],'Independent source digest')
    require({'forced_metric.py','metric_wkb.py','stable_baselines.py'}<=seen,'Independent executable/subtraction/stable-baseline manifest missing')
    require(modes['configuration']==registration['frozen_configuration']['independent_configuration'],'Independent configuration')
    require(modes['configuration']['wronskian_projection'] is False and modes['configuration']['state']=='incoming_BD' and modes['configuration']['phi']=='fixed','Independent state/no projection')
    runtime=modes['runtime'];require(runtime['python'].startswith('3.12.') and runtime['numpy']=='2.2.6' and runtime['python_optimization']==0,'Independent pinned runtime')
    require(runtime['real_dtype']==np.dtype(np.longdouble).name and runtime['complex_dtype']==np.dtype(np.clongdouble).name and runtime['longdouble_nmant']>=63 and runtime['longdouble_eps']==str(np.finfo(LD).eps),'Independent precision')
    expected={'q':'a0^2 deltaQ/epsilon','q_prime':"(a0^2 deltaQ)'/epsilon",'q_second':"(a0^2 deltaQ)''/epsilon",'rho':'a0^4 delta_rho/epsilon','p':'a0^4 delta_p/epsilon','Q0':'physical Q0,K','rho0':'physical rho0,K','p0':'physical p0,K','current':'a0^2 delta_j/epsilon=q, fixed phi, b=1, x_phi=2'}
    require(modes['normalizations']==expected,'Independent normalization declarations')
    resources=modes['resources'];require(0<=resources['elapsed_seconds']<=900 and 0<=resources['peak_rss_kib']<=256*1024,'Independent resources')

def mutation(records,name,residual,gate,kind):
    require(math.isfinite(residual) and residual>=0,'Nonfinite control diagnostic')
    old=records.setdefault(name,{'name':name,'kind':kind,'maximum_absolute_residual':0.0,'operational_tolerance':gate,'operational_gate_exceeded_somewhere':False})
    old['maximum_absolute_residual']=max(old['maximum_absolute_residual'],residual)
    old['operational_gate_exceeded_somewhere'] |= residual>gate

def guard_control(records,name,call):
    try:call()
    except RuntimeError: records[name]={'name':name,'kind':'synthetic metadata/raw-format mutation; no altered-state physical run','rejected':True}
    else: raise RuntimeError('Guard mutation survived: '+name)

def check(primary,modes,experiment,directory):
    check_experiment(experiment);finite_tree(primary,'primary');finite_tree(modes,'independent')
    g=experiment['gates'];pr=row_index(primary,experiment,'primary');mr=row_index(modes,experiment,'independent')
    runs={(r['source'],r['setting']):r for r in modes['runs']}
    require(len(runs)==len(modes['runs'])==4 and set(runs)=={(s,t) for s in experiment['source']['ids'] for t in ('coarse','fine')},'Four distinct registered independent runs')
    raw=[];ward=[];controls={};cross=[];extras=[];refinement=[];current_refinement=[];primary_reductions=[];primary_ward=[];continuum=[];tails=[];ward_refinement=[]
    for (source,setting),run in runs.items():
        for name in ('wronskian_max_scaled','physical_wronskian_max_scaled'):require(0<=run[name]<=g['wronskian_over_epsilon'],'Declared global Wronskian')
        audit=audit_run(run,directory,experiment);raw.append(audit)
        for entry in audit['ward_endpoints']:
            if setting=='fine':require(entry['residual']<=g['ward_endpoint'],'Raw fine Ward endpoint')
            ward.append(dict(entry,source=source,setting=setting))
        for name,residual in audit['raw_mutations'].items():
            gate=g['ward_endpoint'] if 'Ward' in name else g['cross_q'] if 'variance' in name else g['cross_p'] if 'pressure' in name else g['cross_rho']
            mutation(controls,name,residual,gate,'same saved modes/history, alternate operator/subtraction/prefactor; no additional physical run')
    by_run={(s,t):{r['eta']:r for r in run['rows']} for (s,t),run in runs.items()}
    max_refinement={n:0.0 for n in NAMES};max_ward=max_ward_ref=0.0
    for (source,eta),row in pr.items():
        mrow=mr[(source,eta)];hjet=[row['source_jet_over_epsilon'][f'h{n}'] for n in range(6)];gjet=[row['canonical_forcing_jet_over_epsilon'][f'g{n}'] for n in range(4)];a=L=-1/eta
        for K,block in [(None,row['continuum'])]+[(b['K'],b) for b in row['finite_k']]:
            values=block['values'];require(set(values)==set(block['quadrature_error_estimates'])==set(NAMES),'Primary quantity/error keys')
            for name,error in block['quadrature_error_estimates'].items(): require(0<=error<=g['primary_normalized_quad_estimate_max'],'Primary quadrature estimate '+name)
            history=block['log_history_integrals_over_epsilon'] if K is None else block['finite_history_integrals_over_epsilon']
            require([m['forcing_derivative_order'] for m in history]==([1,2,3] if K is None else [0,1,2]),'Primary memory derivative order')
            require(all(math.isfinite(m['quadrature_error_estimate']) and m['quadrature_error_estimate']>=0 for m in history),'Raw memory quadrature estimate')
            if K is not None:
                left,right=max(0.,eta+3),max(0.,eta+5);panels=0 if eta<=-5 else max(1,math.ceil(2*K*(right-left)/math.pi))
                require(all(m['time_panels']==panels for m in history),'Primary finite-memory panel count')
            expected,components,local,Rp=reduced_block(eta,hjet,gjet,[m['value'] for m in history],K)
            e0,e1,e2=[m['quadrature_error_estimate']/(8*math.pi**2) for m in history]
            error_map={'q':e0,'q_prime':e1,'q_second':e2,'rho':(3*L*L*e0+L*e1)/2,'p':(e2+3*L*e1+3*L*L*e0)/6,'Q0':0.0,'rho0':0.0,'p0':0.0,'current':e0}
            for name,error in error_map.items():close(block['quadrature_error_estimates'][name],error,0,'Primary propagated error estimate '+name)
            for name in NAMES:
                gap=close(values[name],expected[name],g['numerical_roundoff_allowance'],'Saved primary memory reduction '+name)
                primary_reductions.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap})
            require(set(block['components'])==set(components) and set(block['local_coefficients'])==set(local),'Primary component/local keys')
            for name,value in components.items():close(block['components'][name],value,g['numerical_roundoff_allowance'],'Primary component '+name)
            for name,value in local.items():close(block['local_coefficients'][name],value,g['numerical_roundoff_allowance'],'Primary local coefficient '+name)
            close(block['density_derivative']['scaled_density_prime_over_epsilon'],Rp,g['numerical_roundoff_allowance'],'Differentiated scaled primary density')
            require(0<=block['density_derivative']['quadrature_error_estimate']<=g['primary_normalized_quad_estimate_max'],'Primary density derivative error estimate')
            close(block['density_derivative']['quadrature_error_estimate'],3*L**3*e0+L*L*e1+L*e2/2,0,'Primary propagated density-derivative error estimate')
            close(block['trace_from_direct_stresses'],-values['rho']+3*values['p'],g['numerical_roundoff_allowance'],'Primary trace diagnostic')
            ledger=L*(values['rho']-3*values['p'])-3*hjet[1]*a**4*(values['rho0']+values['p0'])
            gap=close(Rp,ledger,g['numerical_roundoff_allowance'],'Primary differentiated metric Ward identity')
            primary_ward.append({'source':source,'eta':eta,'K':K,'gap':gap})
            if eta<=-5:require(all(values[n]==0 for n in (*CORE,'current')),'Primary initial causality')
            close(values['current'],values['q'],0,'Fixed-source metric scalar current')
            if K is None:continue
            point=next(p for p in mrow['finite_k'] if p['K']==K);mv=point['values'];require(set(mv)==set(NAMES),'Independent quantity keys')
            tail=metric_tail(source,eta,K);require(tail['interval_dps']==60,'Registered metric tail interval precision');finite_tree(tail,'metric tail');tails.append(tail)
            for name in CORE:
                gap=close(values[name],mv[name],g['cross_'+name],'Cross-route metric '+name)
                cross.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'tolerance':g['cross_'+name]})
            for name in EXTRAS:
                gap=close(values[name],mv[name],g['cross_'+name],'Cross-route baseline/current '+name)
                extras.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'tolerance':g['cross_'+name]})
            coarse=next(p for p in by_run[(source,'coarse')][eta]['finite_k'] if p['K']==K);fine=next(p for p in by_run[(source,'fine')][eta]['finite_k'] if p['K']==K)
            for name in NAMES:
                close(mv[name],fine['values'][name],0,'Fine run/summary '+name);close(point['coarse_values'][name],coarse['values'][name],0,'Coarse run/summary '+name)
                difference=abs(fine['values'][name]-coarse['values'][name]);max_refinement[name]=max(max_refinement[name],difference)
                close(point['refinement_differences'][name],difference,0,'Declared refinement '+name);close(point['estimated_numerical_errors'][name],2*difference,0,'Declared numerical estimate '+name)
                if name in CORE:
                    require(difference<=g['refinement_'+name],'Core refinement '+name);refinement.append({'source':source,'eta':eta,'K':K,'quantity':name,'difference':difference,'tolerance':g['refinement_'+name]})
                if name=='current':
                    require(difference<=2e-10,'Current refinement');current_refinement.append({'source':source,'eta':eta,'K':K,'difference':difference,'tolerance':2e-10})
                target=row['continuum']['values'][name];gap=abs(target-mv[name])
                analytic=tail[name];quad=row['continuum']['quadrature_error_estimates'][name];numerical=point['estimated_numerical_errors'][name];cross_allowance=g['cross_'+name];arithmetic=g['numerical_roundoff_allowance']
                allowance=analytic+quad+numerical+cross_allowance+arithmetic
                require(analytic>=0 and gap<=allowance,'Metric omitted-band comparison '+name)
                continuum.append({'source':source,'eta':eta,'K':K,'quantity':name,'absolute_gap':gap,'primary_finite_K_gap':abs(target-values[name]),'analytic_tail_upper_bound':analytic,'primary_continuum_quadrature_estimate':quad,'independent_refinement_estimate':numerical,'cross_route_allowance':cross_allowance,'arithmetic_allowance':arithmetic,'total_comparison_allowance':allowance,'status':'omitted UV band enclosed separately; loose bound, empirical finite errors; no certified total error or precise continuum sign'})
            wf, wc=fine['ward'],coarse['ward'];difference=abs(wf['ledger']-wc['ledger'])
            close(point['ward']['refinement_difference'],difference,0,'Ward coarse/fine summary');close(point['ward']['coarse_ledger'],wc['ledger'],0,'Ward coarse ledger summary');close(point['ward']['coarse_endpoint_difference'],wc['endpoint_difference'],0,'Ward coarse residual summary')
            require(wf['endpoint_difference']<=g['ward_endpoint'] and difference<=g['ward_refinement'],'Fine Ward/refinement gate')
            max_ward=max(max_ward,wf['endpoint_difference']);max_ward_ref=max(max_ward_ref,difference)
            ward_refinement.append({'source':source,'eta':eta,'K':K,'difference':difference,'tolerance':g['ward_refinement']})
    for name in NAMES:close(modes['maximum_refinement_differences'][name],max_refinement[name],0,'Maximum refinement summary '+name)
    close(modes['maximum_ward_endpoint_difference'],max_ward,0,'Maximum Ward summary');close(modes['maximum_ward_refinement_difference'],max_ward_ref,0,'Maximum Ward refinement summary')
    for name,call in [('nonfinite_summary',lambda:finite_tree({'x':float('nan')},'synthetic')),('insufficient_raw_precision',lambda:require(np.dtype('float64')==np.dtype(LD),'synthetic dtype')),('changed_incoming_state',lambda:require('alternative'=='incoming_BD','synthetic state')),('Wronskian_projection_declared',lambda:require(True is False,'synthetic projection')),('moving_reference_scale',lambda:require(3==experiment['geometry']['r'],'synthetic reference'))]:guard_control(controls,name,call)
    require(len(controls)==15,'Ten archived operator/history and five synthetic controls required')
    require(all(c.get('maximum_absolute_residual',1)>0 for c in controls.values()),'Saved-mode control has no sensitivity witness')
    require(len(cross)==180 and len(extras)==144 and len(refinement)==180 and len(ward_refinement)==36 and sum(r['reconstructed_observation_points'] for r in raw)==72 and len(ward)==72,'Registered comparison coverage')
    return {'core_comparisons':cross,'baseline_current_comparisons':extras,'primary_memory_contact_reductions':primary_reductions,'primary_differentiated_Ward_checks':primary_ward,'raw_archive_audits':raw,'raw_direct_history_Ward_endpoints':ward,'Ward_refinement_points':ward_refinement,'core_refinement_comparisons':refinement,'current_refinement_comparisons':current_refinement,'continuum_tail_comparisons':continuum,'metric_tail_envelopes':tails,'controls':list(controls.values()),'counts':{'core_comparisons':len(cross),'extras':len(extras),'primary_memory_contact_reductions':len(primary_reductions),'primary_differentiated_Ward':len(primary_ward),'raw_reconstructed_observation_K_points':sum(r['reconstructed_observation_points'] for r in raw),'raw_Ward_endpoints':len(ward),'Ward_refinement_points':len(ward_refinement),'core_refinement_comparisons':len(refinement),'current_refinement_comparisons':len(current_refinement),'continuum_tail_comparisons':len(continuum),'controls':len(controls)}}

def validate_output_path(path):
    output=Path(path).resolve();require(not output.exists() and not output.is_relative_to(ROOT),'Fresh validator output must be outside frozen checkpoint');return output

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary',type=Path,required=True);parser.add_argument('--modes',type=Path,required=True)
    parser.add_argument('--registration-sha256',required=True);parser.add_argument('--public-freeze-commit',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=validate_output_path(args.output);output.mkdir(parents=True)
    try:
        require(np.__version__=='2.2.6' and sys.version_info[:2]==(3,12) and np.finfo(LD).nmant>=63,'Validator pinned runtime/precision')
        hexadecimal(args.public_freeze_commit,40,'Public freeze commit');registration,experiment=check_registration(args.registration_sha256)
        primary=json.loads(args.primary.read_text());modes=json.loads(args.modes.read_text());finite_tree(primary,'primary');finite_tree(modes,'independent')
        primary_provenance(primary,args.registration_sha256,args.public_freeze_commit,experiment);independent_provenance(modes,args.registration_sha256,args.public_freeze_commit,experiment,registration)
        result=check(primary,modes,experiment,args.modes.resolve().parent)
        result.update(status='PASS',schema_version=1,python_optimization=sys.flags.optimize,registration_sha256=args.registration_sha256,public_freeze_commit=args.public_freeze_commit,input_sha256={'primary':sha(args.primary),'modes':sha(args.modes)},scope='Linear homogeneous metric-response calibration at H1,x=r2,xi0 on prescribed geometry and fixed incoming BD state; fixed source phi. No actual shifted-root propagator, coupled evolution, stability, particle yield or heating.',uncertainty='Quadrature, resolution refinement and floating-point allowances are empirical, separate numerical diagnostics. The separately frozen directed module encloses only omitted UV bands, with deliberately loose derivative caps; no total-error certificate, precise continuum pressure or sign is inferred. Observation Wronskians and all direct-history Ward endpoints are reconstructed; global all-time maxima remain producer declarations. Pure contact algebra is shared source-pinned action evidence, not a third independent physical implementation. Stable baseline rationalizations are independently proved against generic WKB; literal bare/subtraction differences retain leading-term ulp diagnostics.')
        write(output/'CHECKS.json',result);print(json.dumps({'status':result['status'],'counts':result['counts']}))
    except BaseException as exc:
        write(output/'FAILURE.json',{'status':'FAIL','exception':repr(exc),'traceback':traceback.format_exc(),'python_optimization':sys.flags.optimize});raise

if __name__=='__main__':main()
