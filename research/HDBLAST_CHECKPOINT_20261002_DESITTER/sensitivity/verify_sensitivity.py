#!/usr/bin/env python3
"""Symbolic derivation, refinement report and registered negative controls."""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import sympy as s
ROOT=Path(__file__).resolve().parent

def symbolic():
    q,h,w,f,g,z,v,up=s.symbols('q h w f g z v up')
    S=v+w*f/3
    Sprime=2*h*z-2*w*g/3+(up-4*q*w)*f/3+w*g/3
    identity=s.expand(s.expand(Sprime+4*q*S+2*h*z).subs(up*f,w*g-12*q*v-12*h*z))
    assert identity==0
    e,ph=s.symbols('e ph');W=1-ph+ph**3/3
    pot=s.expand((s.diff(W,ph)**2/2-s.Rational(2,3)*W**2).subs(ph,1+e))
    expected=-s.Rational(2,27)+s.Rational(14,9)*e**2+s.Rational(50,27)*e**3-e**4/6-s.Rational(4,9)*e**5-s.Rational(2,27)*e**6
    assert s.expand(pot-expected)==0
    y,u,u1,u2=s.symbols('y u u1 u2');a,b,c,d=s.symbols('a b c d')
    r=y+a*y**3+b*y**5;eta=c*y**2+d*y**4;w0=s.diff(eta,y)
    grav=s.series(s.diff(r,y,2)+r*(w0*w0/4+(u+u1*eta+u2*eta**2/2)/6),y,0,4).removeO()
    scalar=s.series(s.diff(eta,y,2)-(u1+u2*eta)+4*s.diff(r,y)/r*w0,y,0,3).removeO()
    coeff={a:-u/36,b:u*u/4320-u1*u1/750,c:u1/10,d:u1*(u2/280+u/630)}
    assert s.simplify(grav.subs(coeff))==0 and s.simplify(scalar.subs(coeff))==0
    return {'constraint_to_weighted_variation_identity':True,'shifted_potential_coefficients':True,'regular_cone_series_orders_checked':['R y^3,y^5','phi y^2,y^4']}

def main():
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,default=ROOT/'RUNS.json');p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    data=json.loads(a.runs.read_text());runs=[r for r in data['runs'] if 'J' in r];cs=[r for r in data['runs'] if r['method']=='complex_step']
    ref=next(r for r in runs if r['method']=='DOP853' and r['max_step']==.05)
    gates={};stats={}
    for key in ['phi_ell','H2_ell','E1_ell_integral','E1_ell_direct','J','coordinate_response','observable_response']:
        arr=np.array([r[key] for r in runs]);delta=np.max(abs(arr-np.array(ref[key])),axis=0);span=np.ptp(arr,axis=0)
        stats[key]={'full_span':span.tolist(),'max_difference_from_reference':delta.tolist(),'conservative_empirical_error_scale':(10*np.maximum(delta,span)).tolist()}
    # Independent CS entries are included in the componentwise acceptance.
    csdiff={k:max(abs(r[k]-ref[k if k!='E1_ell' else 'E1_ell_integral']) for r in cs) for k in ('phi_ell','H2_ell','E1_ell')}
    csdiff['E2_ell']=max(abs(r['E2_ell']-ref['J'][1][0]) for r in cs)
    gates['phi_ell_absolute']=max(stats['phi_ell']['full_span'],csdiff['phi_ell'])<2e-10
    gates['H2_ell_absolute']=max(stats['H2_ell']['full_span'],csdiff['H2_ell'])<2e-14
    gates['ordinary_J_absolute']=np.max(np.array(stats['J']['full_span']))<2e-10
    gates['tiny_J00_identity']=max(abs(r['integral_identity_error']) for r in runs)<max(1e-19,1e-5*abs(ref['E1_ell_integral']))
    gates['tiny_J00_independent_methods']=max(stats['E1_ell_integral']['full_span'],csdiff['E1_ell'])<max(1e-19,1e-5*abs(ref['E1_ell_integral']))
    gates['tiny_J00_positive_all']=all(r['E1_ell_integral']>0 for r in runs) and all(r['E1_ell']>0 for r in cs)
    obs_sp=np.array(stats['observable_response']['full_span']);gates['phi_response_absolute']=bool(np.max(obs_sp[0])<1e-7);gates['H2_response_absolute']=bool(np.max(obs_sp[1])<1e-10)
    gates['constraint_response_absolute']=max(np.max(abs(np.array(r['constraint_response_error']))) for r in runs)<1e-14
    # Compute y0 envelope only among finest same-step primary runs.
    starts=[r for r in runs if r['method']=='DOP853' and r['rtol']==3e-14 and r['max_step']==.1]
    cone={k:np.ptp(np.array([r[k] for r in starts]),axis=0).tolist() for k in ('phi_ell','H2_ell','E1_ell_integral','observable_response')}
    J=np.array(ref['J']);Y=np.array(ref['output_J']);T=np.diag([1/6,-1/2]);coord=np.array(ref['coordinate_response']);obs=np.array(ref['observable_response'])
    wrong_current=np.linalg.solve(J,np.diag([1/6,1/2]));scalar_sign_error=float(abs((J@wrong_current-T)[1,1]));assert scalar_sign_error>.9
    missing_moving=np.outer(Y[:,0],coord[0,:]);missing_error=abs(missing_moving-obs);assert missing_error[0,0]>.1 and missing_error[1,0]>.03
    sigma1=4*ref['eta']+2*ref['eta']**2+.001*.5975949350280132
    missing_scalar_junction=ref['E1_ell_direct']+sigma1*ref['phi_ell']/6
    omission_error=abs(missing_scalar_junction-ref['E1_ell_integral']);assert omission_error>1e-9
    gates['negative_controls_detected']=True;gates['no_run_failures']=len(data['failures'])==0
    # Newton correction estimates quantify how close fixed archived coordinates are to an exact vacuum.
    correction=-np.linalg.solve(J,np.array(ref['residual']))
    source_absolute=np.abs(obs);source_relative=source_absolute[1]/ref['H2']
    out={'status':'PASS' if all(gates.values()) else 'FAIL','scope':'classical local differential response; empirical floating-point evidence, no interval error certification','symbolic':symbolic(),'gates':{k:bool(v) for k,v in gates.items()},'reference':ref,'refinement_statistics':stats,'complex_step_max_difference':csdiff,'finest_start_spread':cone,'baseline_newton_correction_estimate_ell_y':correction.tolist(),'negative_controls':{'drop_scalar_junction_metric_column':{'wrong':missing_scalar_junction,'identity_discrepancy':omission_error},'reverse_current_sign':{'target_residual_discrepancy':scalar_sign_error},'omit_moving_endpoint':{'wrong_observable_response':missing_moving.tolist(),'absolute_discrepancy':missing_error.tolist()}},'small_response_coefficients':{'absolute_phi':source_absolute[0].tolist(),'relative_H2':source_relative.tolist()},'run_count':len(data['runs']),'failures':data['failures'],'inputs':{'runs_sha256':hashlib.sha256(a.runs.read_bytes()).hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}}
    a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':out['status'],'gates':out['gates'],'finest_start_spread':cone,'complex_step_max_difference':csdiff,'small_response_coefficients':out['small_response_coefficients'],'baseline_correction':correction.tolist()},indent=2));assert out['status']=='PASS'

if __name__=='__main__':main()
