#!/usr/bin/env python3
"""Independent residual/source and empirical refinement checks; no assert gates."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import traceback
import mpmath as mp
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DESITTER='research/HDBLAST_CHECKPOINT_20261002_DESITTER'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_module(path):
    spec=importlib.util.spec_from_file_location('validation_action',path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value)
    return value

def write_json(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')

def richardson(fun,at):
    h=abs(at)*mp.mpf('0.0001')
    full=(fun(at+h)-fun(at-h))/(2*h)
    half=(fun(at+h/2)-fun(at-h/2))/h
    return (4*half-full)/3

def validate(repo,runs_path,output):
    if output.exists():
        raise RuntimeError('Refusing to overwrite validator output')
    data=json.loads(runs_path.read_text())
    protocol=json.loads((ROOT/'PREREGISTRATION.json').read_text())
    experiment=protocol['experiment'];tol=experiment['acceptance']
    fixed=json.loads((ROOT/'FIXED_MODEL.json').read_text())
    mp.mp.dps=80
    ref=mp.mpf(fixed['r_ref']);eta0=mp.mpf(fixed['eta_ref']);z0=mp.mpf(fixed['z_ref'])
    delta=mp.mpf(fixed['delta']);c=mp.mpf(fixed['c'])
    action=load_module(repo/DESITTER/'code/desitter_sources.py')
    gates=[];diagnostics=[]
    def gate(name,condition,**details):
        gates.append({'name':name,'pass':bool(condition),**details})
    def close(name,observed,expected,absolute,relative=0.):
        error=abs(mp.mpf(observed)-mp.mpf(expected))
        threshold=mp.mpf(absolute)+mp.mpf(relative)*abs(mp.mpf(expected))
        gate(name,error<=threshold,absolute_error=float(error),threshold=float(threshold))
    gate('registration_hash',data['registration_sha256']==sha(ROOT/'PREREGISTRATION.json'))
    gate('producer_hash',data['producer_sha256']==protocol['local_sha256']['code/stationary_shell.py'])
    for name,digest in protocol['local_sha256'].items():
        gate('frozen_local/'+name,sha(ROOT/name)==digest)
    for name,digest in protocol['source_sha256'].items():
        gate('pinned_input/'+name,sha(repo/name)==digest)
    gate('execution_complete',data.get('execution_status')=='COMPLETE')
    gate('no_run_failures',not data.get('failures'))
    expected=[]
    for resolution in experiment['resolutions']:
        grid=experiment['gamma_values'] if resolution['coverage']=='all' else [0,experiment['gamma_values'][-1]]
        expected.extend((resolution['name'],b,g) for b in experiment['b_values'] for g in grid)
    observed=[(r['resolution'],r['b'],r['gamma']) for r in data['runs']]
    gate('exact_registered_coverage',len(observed)==len(expected) and set(observed)==set(expected),
         expected_count=len(expected),observed_count=len(observed))
    index={(r['resolution'],r['b'],r['gamma']):r for r in data['runs']}
    for row in data['runs']:
        label=f"{row['resolution']}/b{row['b']}/gamma{row['gamma']}"
        try:
            e=mp.mpf(row['eta']);w=mp.mpf(row['w']);R=mp.mpf(row['R']);z=1/R**2;b=row['b'];gamma=mp.mpf(row['gamma'])
            factor=1+b*(e-eta0)/2
            x=ref*factor*factor;xphi=b*ref*factor;xphiphi=b*b*ref/2
            src=action.action_sources(x,z,ref,dps=80,N=256)
            W=mp.mpf(1)/3+e**2+e**3/3;p=2*e+e**2
            U=p*p/2-2*W*W/3
            sigma=2*W+delta*(1+c*(1+e))
            sigmap=2*p+delta*c
            q=mp.sqrt(z+w*w/12-U/6)
            density=gamma*src['rho'];current=gamma*xphi*src['Q']/2
            C=z+w*w/12-U/6-(sigma+density)**2/36
            B=w+sigmap/2+current/2
            E1=q-(sigma+density)/6
            gate(label+'/metric_residual',abs(C)<=tol['C_absolute'],absolute_residual=float(abs(C)),threshold=tol['C_absolute'])
            gate(label+'/scalar_residual',abs(B)<=tol['B_absolute'],absolute_residual=float(abs(B)),threshold=tol['B_absolute'])
            close(label+'/stable_C',row['C'],C,3e-18)
            close(label+'/original_metric_residual',row['residual'][0],E1,3e-17)
            close(label+'/H2_radius',row['H2'],z,1e-18)
            gate(label+'/positive_tension',sigma+density>0)
            gate(label+'/negative_cone_displacement',row['eh']<0)
            gate(label+'/declared_domain',mp.mpf('.75')<=z/z0<=mp.mpf('1.25') and abs(e-eta0)<=mp.mpf('.1')*abs(eta0) and 0<x/z<mp.mpf('4.5') and factor>0)
            determinant=float(np.linalg.det(np.asarray(row['jacobian'])))
            gate(label+'/nonsingular_numerical_jacobian',math.isfinite(determinant) and abs(determinant)>1e-12,
                 determinant=determinant,condition_number=row['jacobian_condition'])
            gate(label+'/constraint_jacobian_identity',max(map(abs,row['jacobian_constraint_discrepancy']))<=tol['constraint_jacobian_absolute'],
                 maximum_absolute_error=max(map(abs,row['jacobian_constraint_discrepancy'])))
            for key in ['W','rho','Q']:
                close(label+'/source_refinement/'+key,row['source'][key],src[key],tol['source_refinement_absolute'],tol['source_refinement_relative'])
                gate(label+'/source_tail/'+key,src['tail_bound'][key]<=tol['source_refinement_absolute'],
                     bound=float(src['tail_bound'][key]),threshold=tol['source_refinement_absolute'])
            for key,bound in row['source']['derivative_tail_bound'].items():
                gate(label+'/feedback_derivative_tail/'+key,bound<=tol['source_derivative_absolute'],
                     bound=bound,threshold=tol['source_derivative_absolute'])
            close(label+'/mass_law',row['source']['x'],x,1e-19)
            close(label+'/paired_current',row['source']['j'],xphi*src['Q']/2,1e-22,1e-10)
            h=x-ref
            anomaly=(h*h/2-2*h*z+mp.mpf(29)*z*z/15)/(16*mp.pi**2)
            close(label+'/independent_trace',-4*src['rho']+x*src['Q'],anomaly,1e-35,1e-25)
            close(label+'/digamma_Q',action.digamma_Q(x,z,ref),src['Q'],1e-30,1e-20)
            energy=max(math.sqrt(float(z)),math.sqrt(float(x)),abs(float(q)),
                       row['hierarchy']['max_sampled_sqrt_abs_kT'],row['hierarchy']['max_sampled_sqrt_abs_kN'])
            close(label+'/hierarchy_formula',row['hierarchy']['energy_over_M5'],energy*row['gamma']**(1/3),1e-13)
            gate(label+'/hierarchy_classification',row['hierarchy']['tenfold_screen_pass']==(energy*row['gamma']**(1/3)<=.1))
            if row['resolution']=='refined' and row['gamma'] in (0,experiment['gamma_values'][-1]):
                S=lambda xx,zz:action.action_sources(xx,zz,ref,dps=80,N=256)
                deriv={'rho_x':richardson(lambda v:S(v,z)['rho'],x),
                       'rho_z':richardson(lambda v:S(x,v)['rho'],z),
                       'Q_x':richardson(lambda v:S(v,z)['Q'],x),
                       'Q_z':richardson(lambda v:S(x,v)['Q'],z)}
                for name,independent in deriv.items():
                    close(label+'/action_derivative/'+name,row['source'][name],independent,
                          tol['source_derivative_absolute'],tol['source_derivative_relative'])
                close(label+'/action_pairing',deriv['rho_x'],src['Q']/2-z*deriv['Q_z']/4,1e-20,1e-8)
                close(label+'/current_eta_chain_rule',row['source']['j_eta'],(xphiphi*src['Q']+xphi*xphi*deriv['Q_x'])/2,1e-20,1e-8)
                close(label+'/current_z_chain_rule',row['source']['j_z'],xphi*deriv['Q_z']/2,1e-20,1e-8)
        except Exception:
            gate(label+'/validation_exception',False,traceback=traceback.format_exc())
    for row in data['runs']:
        if row['resolution']=='refined':continue
        reference=index.get(('refined',row['b'],row['gamma']))
        if reference is None:continue
        label=f"refinement/{row['resolution']}/b{row['b']}/gamma{row['gamma']}"
        close(label+'/H2',row['H2'],reference['H2'],tol['cross_resolution_H2_absolute'])
        close(label+'/eta',row['eta'],reference['eta'],tol['cross_resolution_eta_absolute'])
        for j in range(2):
            close(label+'/coordinate'+str(j),row['coordinates'][j],reference['coordinates'][j],tol['cross_resolution_coordinate_absolute'])
    shifts=[]
    for b in experiment['b_values']:
        control=index.get(('refined',b,0))
        if control is None:continue
        for gamma in experiment['gamma_values'][1:]:
            row=index.get(('refined',b,gamma))
            if row is None:continue
            record={'b':b,'gamma':gamma,'observables':{},'hierarchy':row['hierarchy']}
            for field,prediction_index,floor in [('eta',0,1e-17),('H2',1,1e-18)]:
                delta_value=row[field]-control[field]
                candidates=[]
                for resolution in experiment['resolutions']:
                    endpoint=index.get((resolution['name'],b,gamma));zero=index.get((resolution['name'],b,0))
                    if endpoint is not None and zero is not None:
                        candidates.append(endpoint[field]-zero[field])
                spread=max(candidates)-min(candidates)
                threshold=tol['resolution_spread_multiple']*spread+floor
                linear=row['linear_prediction_delta_observables'][prediction_index]
                nonlinear=delta_value-linear
                nonlinear_threshold=max(threshold,tol['nonlinear_fraction_of_shift_min']*abs(delta_value))
                record['observables'][field]={'delta_from_own_solved_control':delta_value,
                    'empirical_spread':spread,'resolution_count':len(candidates),
                    'resolution_threshold':threshold,'shift_resolved':abs(delta_value)>threshold,
                    'archived_linear_prediction':linear,'nonlinear_departure':nonlinear,
                    'nonlinear_fraction':abs(nonlinear/delta_value) if delta_value else None,
                    'nonlinear_detection_threshold':nonlinear_threshold,
                    'nonlinear_resolved':abs(nonlinear)>=nonlinear_threshold}
            shifts.append(record)
    mutation_names=[r['mutation'] for r in data['negative_controls']]
    gate('negative_control_coverage',len(mutation_names)==len(experiment['negative_controls']) and set(mutation_names)==set(experiment['negative_controls']))
    for row in data['negative_controls']:
        e=mp.mpf(row['eta']);w=mp.mpf(row['w']);z=1/mp.mpf(row['R'])**2
        b=row['b'];gamma=mp.mpf(row['gamma']);f=1+b*(e-eta0)/2
        x=ref*f*f;xphi=b*ref*f
        src=action.action_sources(x,z,ref,dps=80,N=256)
        W=mp.mpf(1)/3+e*e+e**3/3;p=2*e+e*e
        U=p*p/2-2*W*W/3;sigma=2*W+delta*(1+c*(1+e))
        C=z+w*w/12-U/6-(sigma+gamma*src['rho'])**2/36
        B=w+p+delta*c/2+gamma*xphi*src['Q']/4
        factor=tol['negative_control_residual_factor']
        close('negative_control_reconstruction/'+row['mutation']+'/C',row['correct_equation_residual'][0],C,3e-18)
        close('negative_control_reconstruction/'+row['mutation']+'/B',row['correct_equation_residual'][1],B,3e-17)
        gate('negative_control/'+row['mutation'],abs(C)>factor*tol['C_absolute'] or abs(B)>factor*tol['B_absolute'],
             correct_equation_residual=[float(C),float(B)],thresholds=[factor*tol['C_absolute'],factor*tol['B_absolute']])
        gate('wrong_model_root/'+row['mutation'],abs(row['C'])<=tol['C_absolute'] and abs(row['B'])<=tol['B_absolute'])
    passed=all(g['pass'] for g in gates)
    result={'status':'PASS' if passed else 'FAIL','checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256':sha(runs_path),'registration_sha256':sha(ROOT/'PREREGISTRATION.json'),
        'validator_sha256':sha(__file__),'gates':gates,'shift_diagnostics':shifts,
        'run_count':len(data['runs']),'failures':data.get('failures',[]),
        'scope':'Empirical numerical stationary closure in a declared dimensionless model family; hierarchy screen and resolved nonlinear feedback are separately classified.',
        'high_gamma_nonlinear_resolved':any(any(v['nonlinear_resolved'] for v in s['observables'].values()) for s in shifts if s['gamma']==experiment['gamma_values'][-1]),
        'hierarchy_screen_failures':[{'b':s['b'],'gamma':s['gamma']} for s in shifts if not s['hierarchy']['tenfold_screen_pass']]}
    output.parent.mkdir(parents=True,exist_ok=True)
    write_json(output,result)
    print(result['status'],len(gates),'gates;',sum(not g['pass'] for g in gates),'failed;',len(shifts),'shift rows')
    return 0 if passed else 1

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo',type=Path,default=Path('/workspace/HDblast'))
    ap.add_argument('--runs',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    return validate(args.repo,args.runs,args.output)

if __name__=='__main__':
    sys.exit(main())
