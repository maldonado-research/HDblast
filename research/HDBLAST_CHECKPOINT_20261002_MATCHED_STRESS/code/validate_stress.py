#!/usr/bin/env python3
"""Registered matched-stress evidence validator, with explicit guards under -O.

No numerical producer is imported. Run only after a publicly verified freeze.
Numerical integration/refinement uncertainties remain estimates. Raw-mode
observations and direct-history Ward endpoints are separately reconstructed.
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
from raw_stress_audit import NAMES, audit_run, jets, near_array

ROOT=Path(__file__).resolve().parents[1]
CORE=('q','q_prime','q_second','rho','p')

def require(condition,message):
    if not condition:raise RuntimeError(message)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def close(actual,expected,tolerance,label):
    require(math.isfinite(actual) and math.isfinite(expected),label+': nonfinite')
    difference=abs(actual-expected)
    require(difference<=tolerance,label+': mismatch; residual='+str(difference)+' tolerance='+str(tolerance))
    return difference

def write(path,data):Path(path).write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')

def check_registration(pin):
    path=ROOT/'FULL_REGISTRATION.json'
    require(sha(path)==pin,'Registration hash mismatch')
    reg=json.loads(path.read_text())
    for name,digest in reg['files'].items():
        target=(ROOT/name).resolve()
        require(target.is_relative_to(ROOT),'Registration path escape')
        require(sha(target)==digest,'Frozen input changed: '+name)
    needed={'EXPERIMENT.json','code/stress_primary.py','code/source_jet.py','code/validate_stress.py','code/raw_stress_audit.py','independent/forced_stress.py','independent/MANIFEST.json','theory/stress_tail_bounds.py'}
    require(needed<=set(reg['files']),'Required registered validator/producer/theory input absent')
    return reg

def row_index(document,experiment,label):
    required={(s,t) for s in experiment['source']['ids'] for t in experiment['observation_eta']}
    indexed={(r['source'],r['eta']):r for r in document['rows']}
    require(len(indexed)==len(document['rows']) and set(indexed)==required,label+': source/time grid')
    for row in indexed.values():
        require([p['K'] for p in row['finite_k']]==experiment['cutoffs'],label+': cutoff grid')
        close(row['a'],-1/row['eta'],1e-15,label+': scale factor')
        declared=np.array([row['source_jet_over_epsilon'][f'f{n}'] for n in range(6)],dtype=np.longdouble)
        expected=jets(np.array([row['eta']],dtype=np.longdouble),row['source'])[:,0]
        near_array(declared,expected,label+': registered source derivatives')
    return indexed

def local(eta,K=None):
    a=L=-1/eta
    v=1.0 if K is None else K/math.hypot(K,math.sqrt(2)*a)
    Q=1/(12*math.pi**2) if K is None else (v*v/(2*(1+v))-v**3/48-v**5/16)/(2*math.pi**2)
    return {'Q0':Q,'J5':v**3/6,'J7':(v**3/3-v**5/5)/4,'J9':(v**3/3-2*v**5/5+v**7/7)/8,'v':v}

def closed(eta,jet,values,K=None):
    """Independent fixed-r integrated stress/contact expressions."""
    L=a=-1/eta;f,f1,f2=jet[:3];q,q1,q2=(values[n] for n in CORE[:3])
    c=local(eta,K);Q=c['Q0'];j5,j7,j9=(c[n] for n in ('J5','J7','J9'))
    rho_variance=(3*L*L*q-L*q1)/2
    rho_mass=a*a*Q*f/2
    rho_local=(3*L*L*f-L*f1)*c['v']**3/(96*math.pi**2)
    p_variance=(q2-3*L*q1-3*L*L*q)/6
    p_mass=-a*a*Q*f/6
    p_local=(70*L*L*f*j9-(30*L*L*f+10*L*f1)*j7+(f2-L*f1-9*L*L*f)*j5)/(48*math.pi**2)
    anomaly=((f2-12*L*L*f)*j5-(30*L*L*f+10*L*f1)*j7+70*L*L*f*j9)/(16*math.pi**2)
    return {'rho':rho_variance+rho_mass+rho_local,'p':p_variance+p_mass+p_local,
            'Q0':Q,'anomaly':anomaly,'current':q+Q*f/4,
            'rho_variance':rho_variance,'rho_mass':rho_mass,'rho_local':rho_local,
            'p_variance':p_variance,'p_mass':p_mass,'p_local':p_local,'current_contact':Q*f/4}

def primary_provenance(primary,registration_pin,freeze,experiment):
    require(primary['status']=='completed' and primary['route']=='logarithmic_variance_derivatives_and_direct_closed_stress','Primary route/status')
    p=primary['provenance']
    require(p['registration_sha256']==registration_pin and p['public_freeze_commit']==freeze,'Primary prospective pins')
    for field,name in (('producer_sha256','code/stress_primary.py'),('source_jet_sha256','code/source_jet.py'),('experiment_sha256','EXPERIMENT.json')):
        require(p[field]==sha(ROOT/name),'Primary evidence input pin '+field)
    require(p['dependencies']=={'numpy':'2.2.6','scipy':'1.15.3','mpmath':'1.3.0','sympy':'1.14.0'},'Primary dependencies')
    require(p['python'].startswith('3.12.'),'Primary Python version')
    require(primary['normalization']==experiment['normalization'] and primary['epsilon']==experiment['source']['epsilon'] and primary['mass_law_b']==1,'Primary units/amplitude/current law')
    require(0<=primary['elapsed_seconds']<=experiment['primary']['producer_wall_seconds'],'Primary wall budget')

def independent_provenance(modes,experiment,freeze):
    require(modes['status']=='passed_internal_gates' and modes['route']=='independent_forced_modes_direct_minimal_stress' and modes['gate_failures']==[],'Independent route/status')
    require(modes['freeze_commit']==freeze and modes['epsilon']==experiment['source']['epsilon'] and modes['mass_law_b']==1,'Independent public pin/amplitude')
    runtime=modes['runtime'];require(runtime['numpy']=='2.2.6' and runtime['python'].startswith('3.12.'),'Independent dependencies')
    require(runtime['real_dtype']==np.dtype(np.longdouble).name and runtime['complex_dtype']==np.dtype(np.clongdouble).name and runtime['longdouble_nmant']>=63,'Independent dtype')
    require(runtime['longdouble_eps']==str(np.finfo(np.longdouble).eps),'Independent arithmetic precision')
    normal={'q':'a^2 deltaQ/epsilon','q_prime':"(a^2 deltaQ)'/epsilon",'q_second':"(a^2 deltaQ)''/epsilon",'rho':'a^4 delta_rho/epsilon','p':'a^4 delta_p/epsilon','Q0':'physical Q0,K','anomaly':'a^4 deltaA_K/epsilon','current':'a^2 delta_j/epsilon'}
    require(modes['normalizations']==normal,'Independent normalization declarations')
    manifest_path=ROOT/'independent/MANIFEST.json';manifest=json.loads(manifest_path.read_text())
    require(modes['provenance']['manifest_sha256']==sha(manifest_path) and modes['provenance']['files']==manifest['files'],'Independent input manifest pin')
    seen=set()
    for entry in manifest['files']:
        target=(manifest_path.parent/entry['path']).resolve()
        require(target.is_relative_to(ROOT) and entry['path'] not in seen,'Independent manifest path/duplicate')
        seen.add(entry['path']);require(sha(target)==entry['sha256'],'Independent registered file hash')
    require('forced_stress.py' in seen and '../EXPERIMENT.json' in seen,'Independent registration source/spec missing')
    resources=modes['resources']
    require(0<=resources['elapsed_seconds']<=experiment['independent']['producer_wall_seconds'],'Independent wall budget')
    require(0<=resources['peak_rss_kib']<=experiment['independent']['peak_rss_kib'],'Independent memory budget')

def mutation(records,name,residual,gate,kind):
    """A sensitivity witness, NOT a claim every physical cross-gate rejects it."""
    require(math.isfinite(residual) and residual>=0,'Nonfinite control diagnostic')
    old=records.setdefault(name,{'name':name,'kind':kind,'maximum_absolute_residual':0.0,'operational_tolerance':gate,'operational_gate_exceeded_somewhere':False})
    old['maximum_absolute_residual']=max(old['maximum_absolute_residual'],residual)
    old['operational_gate_exceeded_somewhere'] |= residual>gate

def guard_control(records,name,call):
    try:call()
    except RuntimeError:
        records[name]={'name':name,'kind':'synthetic metadata/raw-format mutation; no altered-state physical run','rejected':True}
    else:raise RuntimeError('Guard mutation survived: '+name)

def check(primary,modes,experiment,mode_directory):
    require(experiment['geometry']=={'H':1.0,'r':2.0,'xi':0.0,'a':'-1/eta','incoming_state':'fixed BD','geometry_response':False},'Registered fixed geometry/state')
    require(experiment['source']['epsilon']==1e-4 and experiment['source']['mass_law_translation_b']==1,'Registered amplitude/mass law')
    g=experiment['gates'];p_rows=row_index(primary,experiment,'primary');m_rows=row_index(modes,experiment,'independent')
    runs={(r['source'],r['setting']):r for r in modes['runs']}
    require(len(runs)==len(modes['runs'])==4 and set(runs)=={(s,t) for s in experiment['source']['ids'] for t in ('coarse','fine')},'Four registered independent runs required')
    raw=[];controls={};ward=[];cross=[];extras=[];continuum=[];trace=[];closed_checks=[];refinement=[]
    for (source,setting),run in runs.items():
        for name in ('wronskian_max_scaled','physical_wronskian_max_scaled'):
            require(0<=run[name]<=g['wronskian_over_epsilon'],'Declared global Wronskian')
        audit=audit_run(run,mode_directory,experiment);raw.append(audit)
        for entry in audit['ward_endpoints']:
            if setting=='fine':require(entry['residual']<=g['ward_endpoint'],'Raw fine direct-history Ward endpoint')
            ward.append(dict(entry,source=source,setting=setting))
        for name,residual in audit['raw_mutations'].items():
            gate=g['cross_rho' if name.endswith(('density','operator')) and 'pressure' not in name else 'cross_p']
            mutation(controls,name,residual,gate,'same saved modes, explicit alternate operator/subtraction; no additional numerical run')
    fine_runs={(s,t):{r['eta']:r for r in v['rows']} for (s,t),v in runs.items()}
    maxima={n:0.0 for n in NAMES};maxward=maxwardref=0.0
    # Trusted publicly frozen analytic tail module, not a numerical producer.
    spec=importlib.util.spec_from_file_location('registered_stress_tail',ROOT/'theory/stress_tail_bounds.py')
    tail_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(tail_module)
    for key,row in p_rows.items():
        source,eta=key;mr=m_rows[key];jet=[row['source_jet_over_epsilon'][f'f{n}'] for n in range(6)];a=L=-1/eta
        blocks=[(None,row['continuum'])]+[(x['K'],x) for x in row['finite_k']]
        for K,block in blocks:
            v=block['values'];expected=closed(eta,jet,v,K)
            for name,error in block['quadrature_error_estimates'].items():
                require(name in NAMES and math.isfinite(error) and 0<=error<=g['primary_normalized_quad_estimate_max'],'Primary quadrature estimate '+name)
            require(set(v)==set(block['quadrature_error_estimates'])==set(NAMES),'Primary quantity keys')
            for name in ('rho','p','Q0','anomaly','current'):
                close(v[name],expected[name],g['numerical_roundoff_allowance'],'Primary independent closed formula '+name)
            component_names={'rho_variance_response':'rho_variance','rho_reduced_baseline_mass_contact':'rho_mass','rho_local_subtraction_contact':'rho_local','pressure_variance_response':'p_variance','pressure_reduced_baseline_mass_contact':'p_mass','pressure_direct_subtraction_contact':'p_local','current_mass_law_contact':'current_contact'}
            require(set(block['components'])==set(component_names),'Primary component keys')
            for component,key in component_names.items():close(block['components'][component],expected[key],g['numerical_roundoff_allowance'],'Primary component '+component)
            local_values=local(eta,K)
            for name in ('Q0','J5','J7','J9','v'):close(block['local_coefficients'][name],local_values[name],1e-14,'Primary local coefficient '+name)
            vcoef=local_values['v'];vprime=0.0 if K is None else -L*vcoef*(1-vcoef*vcoef)
            dQdv=(vcoef*(2+vcoef)/(2*(1+vcoef)**2)-vcoef*vcoef/16-5*vcoef**4/16)/(2*math.pi**2)
            Qprime=dQdv*vprime
            base=(3*L*L*jet[0]-L*jet[1])/(96*math.pi**2)
            baseprime=(6*L**3*jet[0]+2*L*L*jet[1]-L*jet[2])/(96*math.pi**2)
            Cprime=baseprime*vcoef**3+base*3*vcoef*vcoef*vprime
            physical_rhoprime=(-3*L**3*v['q']+3*L*L*v['q_prime']-L*v['q_second']/2+a*a*((Qprime/2-L*expected['Q0'])*jet[0]+expected['Q0']*jet[1]/2)+Cprime-4*L*expected['rho_local'])
            close(block['density_derivative']['a4_delta_rho_prime_over_epsilon'],physical_rhoprime,g['numerical_roundoff_allowance'],'Independent differentiated density')
            derivative_error=block['density_derivative']['quadrature_error_estimate']
            require(math.isfinite(derivative_error) and 0<=derivative_error<=g['primary_normalized_quad_estimate_max'],'Density derivative numerical estimate')
            near_array(np.array(block['variance_jet_over_epsilon']),np.array([v['q']/a**2,(v['q_prime']-2*L*v['q'])/a**2,(v['q_second']-4*L*v['q_prime']+2*L*L*v['q'])/a**2]),'Physical variance derivatives')
            if eta<=experiment['source']['support'][0]:
                require(all(v[n]==0 for n in NAMES if n!='Q0'),'Primary retarded causality')
            memories=block['log_history_integrals_over_epsilon'] if K is None else block['finite_history_integrals_over_epsilon']
            require([m['source_derivative_order'] for m in memories]==([1,2,3] if K is None else [0,1,2]),'Primary memory derivative registration')
            require(all(math.isfinite(m['value']) and 0<=m['quadrature_error_estimate']<=g['primary_normalized_quad_estimate_max'] for m in memories),'Primary raw memory estimates')
            memory=[m['value'] for m in memories];pref=1/(8*math.pi**2)
            if K is None:
                reduced=[-pref*memory[0],-pref*(memory[1]+L*jet[0]),-pref*(memory[2]+2*L*jet[1]+L*L*jet[0])]
            else:
                lc=block['local_coefficients'];A=math.asinh(K/(math.sqrt(2)*a))-vcoef
                Ap=-L*vcoef**3;App=L*L*vcoef**3*(2-3*vcoef*vcoef)
                for name,value in (('A',A),('A_prime',Ap),('A_second',App)):close(lc[name],value,1e-14,'Finite variance coefficient '+name)
                reduced=[pref*(jet[0]*A-memory[0]),pref*(jet[1]*A+jet[0]*Ap-memory[1]),pref*(jet[2]*A+2*jet[1]*Ap+jet[0]*App-memory[2])]
                right=max(0.0,eta-experiment['source']['support'][0]);left=max(0.0,eta-experiment['source']['support'][1])
                panels=0 if eta<=experiment['source']['support'][0] else max(1,math.ceil(2*K*(right-left)/math.pi))
                require(all(m['time_panels']==panels for m in memories),'Finite memory panels changed')
            near_array(np.array([v[n] for n in CORE[:3]]),np.array(reduced),'Primary memory/summary variance derivatives')
            if K is None:continue
            mp=next(x for x in mr['finite_k'] if x['K']==K);mv=mp['values']
            for name in CORE:
                gap=close(v[name],mv[name],g['cross_'+name],'Physical cross-route '+name)
                cross.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'tolerance':g['cross_'+name]})
            for name in ('Q0','anomaly','current'):
                gap=close(v[name],mv[name],g['cross_'+name],'Baseline/anomaly/current cross-route '+name)
                extras.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'tolerance':g['cross_'+name]})
            expected_direct=closed(eta,jet,mv,K)
            for name in ('rho','p'):
                gap=close(mv[name],expected_direct[name],g['cross_'+name],'Independent direct stress vs closed '+name)
                closed_checks.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'tolerance':g['cross_'+name]})
            exact_trace=mv['q_second']/2-L*mv['q_prime']-3*L*L*mv['q']-a*a*mv['Q0']*jet[0]+mv['anomaly']
            gap=close(-mv['rho']+3*mv['p'],exact_trace,g['finite_K_trace'],'Finite-band matched trace')
            trace.append({'source':source,'eta':eta,'K':K,'gap':gap,'tolerance':g['finite_K_trace']})
            coarse=next(x for x in fine_runs[(source,'coarse')][eta]['finite_k'] if x['K']==K)
            fine=next(x for x in fine_runs[(source,'fine')][eta]['finite_k'] if x['K']==K)
            for name in NAMES:
                close(mv[name],fine['values'][name],0,'Fine raw/summary value')
                close(mp['coarse_values'][name],coarse['values'][name],0,'Coarse raw/summary value')
                d=abs(fine['values'][name]-coarse['values'][name]);maxima[name]=max(maxima[name],d)
                close(mp['refinement_differences'][name],d,0,'Refinement summary')
                close(mp['estimated_numerical_errors'][name],2*d,0,'Recorded numerical estimate')
                if name in CORE:
                    require(d<=g['refinement_'+name],'Registered refinement '+name)
                    refinement.append({'source':source,'eta':eta,'K':K,'quantity':name,'difference':d,'tolerance':g['refinement_'+name]})
            W=mp['ward'];wr=abs(fine['ward']['ledger']-coarse['ward']['ledger']);maxwardref=max(maxwardref,wr);maxward=max(maxward,fine['ward']['endpoint_difference'])
            close(W['refinement_difference'],wr,0,'Ward refinement summary');require(wr<=g['ward_refinement'],'Ward refinement')
            for name,exp in (('ledger',fine['ward']['ledger']),('direct_endpoint',fine['ward']['direct_endpoint']),('endpoint_difference',fine['ward']['endpoint_difference']),('coarse_ledger',coarse['ward']['ledger']),('coarse_endpoint_difference',coarse['ward']['endpoint_difference'])):close(W[name],exp,0,'Ward merged summary '+name)
            bound=tail_module.stress_tail_bounds(source,eta,K)
            require(block['analytic_tail_bounds']==bound,'Actual directed-interval stress/derivative tails changed')
            for name in NAMES:
                bn='Q0_difference' if name=='Q0' else name
                allowance=block['quadrature_error_estimates'][name]+row['continuum']['quadrature_error_estimates'][name]+g['numerical_roundoff_allowance']+mp['estimated_numerical_errors'][name]
                gap=close(row['continuum']['values'][name],mv[name],bound[bn]+allowance,'Removed-cutoff enclosure plus numerical allowances '+name)
                continuum.append({'source':source,'eta':eta,'K':K,'quantity':name,'gap':gap,'analytic_omitted_band_bound':bound[bn],'numerical_allowance_estimate':allowance,'total_allowance':bound[bn]+allowance})
            # Specific exact-contact mutations are diagnostic, separate from cross gates.
            for name,residual,gate in (
                ('omit_density_local_contact',abs(expected['rho_local']),g['cross_rho']),
                ('omit_pressure_local_contact',abs(expected['p_local']),g['cross_p']),
                ('wrong_density_contact_sign',2*abs(expected['rho_local']),g['cross_rho']),
                ('wrong_pressure_contact_sign',2*abs(expected['p_local']),g['cross_p']),
                ('omit_current_quadratic_mass_law_contact',abs(expected['current_contact']),g['cross_current']),
                ('wrong_current_contact_sign',2*abs(expected['current_contact']),g['cross_current']),
                ('continuum_Q0_at_finite_K',abs(local(eta)['Q0']-expected['Q0']),g['cross_Q0']),
                ('continuum_anomaly_at_finite_K',abs(closed(eta,jet,v)['anomaly']-expected['anomaly']),g['cross_anomaly']),
                ('wrong_variance_derivative_sign',abs(L*v['q_prime']),g['cross_rho']),
                ('missing_stress_conformal_factor',abs(expected['rho']*(a**4-1)),g['cross_rho'])):
                mutation(controls,name,residual,gate,'algebraic formula/contact diagnostic on registered values; no alternate physical run')
    require(len(cross)==180 and len(ward)==72 and len(trace)==36 and len(closed_checks)==72,'Audit count mismatch')
    for name,value in maxima.items():close(modes['maximum_refinement_differences'][name],value,0,'Global refinement maximum')
    close(modes['maximum_ward_endpoint_difference'],maxward,0,'Global Ward endpoint maximum');close(modes['maximum_ward_refinement_difference'],maxwardref,0,'Global Ward refinement maximum')
    for control in controls.values():
        require(control['maximum_absolute_residual']>0,'Unwitnessed algebraic/raw-operator mutation '+control['name'])
        control['nonzero_algebraic_sensitivity_witness']=True
    guard_control(controls,'moving_reference_scale_declaration',lambda:require(3==experiment['geometry']['r'],'Synthetic moving reference scale'))
    guard_control(controls,'nonfinite_raw_value',lambda:near_array(np.array([float('nan')]),np.array([0.]),'Synthetic nonfinite'))
    guard_control(controls,'insufficient_raw_precision',lambda:require(np.dtype('float64')==np.dtype(np.longdouble),'Synthetic dtype'))
    guard_control(controls,'altered_initial_data',lambda:require(np.count_nonzero(np.array([0.,1.]))==0,'Synthetic initial state'))
    guard_control(controls,'altered_raw_Wronskian',lambda:require(abs(2*.1-.0/2)<=g['wronskian_over_epsilon'],'Synthetic Wronskian'))
    require(len(controls)==20,'Exactly twenty prospective controls required')
    return {'core_comparisons':cross,'baseline_anomaly_current_comparisons':extras,'direct_stress_closed_comparisons':closed_checks,'finite_K_trace_comparisons':trace,'raw_archive_audits':raw,'raw_direct_history_Ward_endpoints':ward,'refinement_comparisons':refinement,'continuum_comparisons':continuum,'controls':list(controls.values()),'counts':{'core_comparisons':len(cross),'extras':len(extras),'direct_stress_closed':len(closed_checks),'finite_K_trace':len(trace),'raw_reconstructed_observation_K_points':sum(r['reconstructed_observation_points'] for r in raw),'raw_Ward_endpoints':len(ward),'Ward_refinement_points':36,'core_refinement_comparisons':len(refinement),'continuum_comparisons':len(continuum),'controls':len(controls)}}

def validate_output_path(path):
    output=Path(path).resolve()
    require(not output.exists() and not output.is_relative_to(ROOT),'Fresh validator output must be outside frozen checkpoint')
    return output

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary',type=Path,required=True);parser.add_argument('--modes',type=Path,required=True)
    parser.add_argument('--registration-sha256',required=True);parser.add_argument('--public-freeze-commit',required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=validate_output_path(args.output)
    output.mkdir(parents=True)
    try:
        reg=check_registration(args.registration_sha256);experiment=json.loads((ROOT/'EXPERIMENT.json').read_text())
        require(len(args.public_freeze_commit)==40 and all(c in '0123456789abcdef' for c in args.public_freeze_commit),'Explicit full public freeze SHA required')
        primary=json.loads(args.primary.read_text());modes=json.loads(args.modes.read_text())
        primary_provenance(primary,args.registration_sha256,args.public_freeze_commit,experiment);independent_provenance(modes,experiment,args.public_freeze_commit)
        result=check(primary,modes,experiment,args.modes.resolve().parent)
        result.update(status='PASS',schema_version=1,python_optimization=sys.flags.optimize,registration_sha256=args.registration_sha256,public_freeze_commit=args.public_freeze_commit,input_sha256={'primary':sha(args.primary),'modes':sha(args.modes)},scope='Fixed-geometry linear homogeneous minimal stress/current calibration only. No metric kernels, shell dynamics, stability, particle yield or heating.',uncertainty='Analytic tails bound omitted bands only. Quadrature, refinement and arithmetic allowances are separately recorded numerical estimates. Global all-time Wronskian maxima are producer declarations; archived observation modes and complete direct-stress Ward histories are reconstructed.')
        write(output/'CHECKS.json',result)
        print(json.dumps({'status':result['status'],'counts':result['counts']}))
    except BaseException as exc:
        write(output/'FAILURE.json',{'status':'FAIL','exception':repr(exc),'traceback':traceback.format_exc(),'python_optimization':sys.flags.optimize})
        raise

if __name__=='__main__':main()
