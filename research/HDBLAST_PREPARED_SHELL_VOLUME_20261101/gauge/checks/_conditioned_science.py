#!/usr/bin/env python3
"""Exact initial pullback, monopole wave, canonical constraint, and adapter tests."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def need(ok,message):
    if not ok:raise RuntimeError(message)


def load(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def run(mutation=None):
    original=load('mean_initial');conditioned=load('conditioned_mean_initial')
    rows=[]
    for n in (1,3,7):
        rho=Q(129+n,1);a=Q(4+n,97);p=-Q(2+n,10000);Vp=Q(3+n,100000)
        ay=-1/rho**2-p*p/3;py=Vp-4*a*p
        xi=Q(5+n,101);xiy=-Q(7+n,103);xiyy=Q(11+n,107)
        f=-Q(13+n,109);fy=Q(17+n,113)
        u=Q(19+n,127);uy=-Q(23+n,131);ut=Q(29+n,137);qm=Q(31+n,139)
        A,Ay,B,By=original.radial_metric(f,fy,xi,xiy,xiyy,a,ay)
        args=dict(rho=rho,a=a,u=u,ut=ut,uy=uy,qm=qm)
        old=original.initial_means(**args,A=A,Ay=Ay,B=B,By=By)
        new=conditioned.conditioned_initial_means(**args,f=f,fy=fy,xi=xi,xiy=xiy,p=p,py=py)
        ell=original.ELL
        if mutation=="omit_scalar_compensation":new['U0']=new['U0'][:4]+(Q(0),)
        # Independent direct Lie derivative initial spacetime values.
        X=-xi;Xt=Q(0);b=Q(0);bt=a*xi
        if mutation=="omit_time_jacobian":bt=Q(0)
        delta_yy=-2*xiy;delta_tt=-2*rho*rho*(a*X+bt)
        delta_angular=2*a*X;delta_phi=p*X
        need(delta_tt==0,'same initial lapse-square value')
        need(tuple(new['U0'][i]-old['U0'][i] for i in range(5))==
            (delta_yy,0,delta_tt/rho**2,delta_angular,delta_phi),'pullback displacement')
        # Use general covariant radial-vector wave at t=0, H=0.
        Xtt=rho*rho*(-xiyy-4*a*xiy)-2*a*rho*rho*bt-rho*rho*(4*ay+8*a*a)*X
        if mutation=="omit_vector_potential":Xtt+=rho*rho*(4*ay+8*a*a)*X
        btt=Q(0)
        delta_gty_t=Xtt-rho*rho*(ay*xi+a*xiy)
        need(new['V0'][1]-old['V0'][1]==ell*delta_gty_t/rho,'mixed velocity vector-wave join')
        need(all(new['V0'][i]==old['V0'][i] for i in (0,2,3,4)),'other velocities unchanged')
        need(btt==0 and Xt*p==0,'scalar and timelike velocity')
        need(new['U0_r'][4]==ell*rho*(py*xi+p*xiy),'conditioned scalar slope')
        # Linearized Hamiltonian of a pure spatial gauge pair vanishes.
        PgLie=p*p*xiy+(4*a*p*p+2*p*py)*xi
        PphiLie=-p*(py*xi+p*xiy)-Vp*p*xi
        need(PgLie+PphiLie==0,'metric-scalar Hamiltonian cancellation')
        need(PgLie!=0,'control: removing metric Lie without scalar fails')
        # Omitting bt=a xi fails the lapse condition; omitting wave curvature fails the join.
        need(-2*rho*rho*a*X!=0,'control: missing time Jacobian')
        badXtt=rho*rho*(-xiyy-4*a*xiy)-2*a*rho*rho*bt
        need(badXtt!=Xtt,'control: missing vector Ricci potential')
        rows.append({'label':n,'pure_second_metric_scalar_pullback':True,
            'same_lapse_shift_values':True,'harmonic_wave_velocity_join':True,
            'other_velocities_and_scalar_slope':True,'linear_Hamiltonian_gauge_identity':True,
            'three_nonzero_mutation_witnesses':True})
    ell=Q(1,1000000);R=Q(131);a=Q(14,125);ay=Q(1,10000);H=Q(3,1000000)
    Az=((3*ell*H,2*ell*a*R*R),(2*ell*a,3*ell*H))
    Bz=((3*ell*a*R,Q(0)),(Q(0),5*ell*a*R))
    Dz=((ell*ell*R*R*(4*ay+8*a*a),6*ell*ell*a*R*R*H),(Q(0),3*ell*ell*(1+H*H)))
    rowcaps=[3*sum(Az[i][j]+Bz[i][j]+Q(3,2)*Dz[i][j] for j in range(2)) for i in range(2)]
    need(max(rowcaps)<Q(3,250),'fresh monopole contraction <0.012')
    return {'status':'PASS','exact_checks':rows,'in_process_control_count':9,'real_subprocess_control_count':0,
        'monopole_wave_row_contraction_caps':[str(v) for v in rowcaps],
        'monopole_wave_contraction_strictly_below':'3/250',
        'scope':'Exact finite algebra and coefficient majorants only; no evolved numerical response.',
        'adapter_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ('mean_initial.py','conditioned_mean_initial.py')}}

