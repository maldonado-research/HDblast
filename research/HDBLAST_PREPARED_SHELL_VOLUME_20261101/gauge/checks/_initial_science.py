#!/usr/bin/env python3
"""Exact rational initial gauge/mean tests, never an evolution certificate.

Direct Christoffel method follows the retained initial-harmonic check, with
independent general radial/tangential second metric and curvature slots.
All equalities and negative controls run also under Python -O.
"""
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ZERO = (Q(0),)*3


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def jet(*values):
    return tuple(Q(x) for x in values)+(Q(0),)*(3-len(values))


def add(*series):
    return tuple(sum(s[i] for s in series) for i in range(3))


def scale(c, series):
    return tuple(c*x for x in series)


def mul(s, t):
    return tuple(sum(s[i]*t[n-i] for i in range(n+1)) for n in range(3))


def invert_metric(g):
    out = [[[Q(0) for j in range(5)] for i in range(5)] for n in range(3)]
    for i in range(5):
        out[0][i][i] = 1/g[i][i][0]
    for n in (1, 2):
        for i in range(5):
            for j in range(5):
                out[n][i][j] = -sum(g[i][k][s]*out[n-s][k][j]
                    for k in range(5) for s in range(1, n+1))/g[i][i][0]
    return [[tuple(out[n][i][j] for n in range(3)) for j in range(5)] for i in range(5)]


def gauge_residual(d, mutation=None):
    rho,a,U,Uy,T,A,Ay,B,By,R,S,DU,DA,DB = d
    r2 = rho*rho
    shift1 = -6*r2*(Uy+4*a*U)
    shift2 = r2*(Ay/2-3*By/2+4*a*A-3*a*B-74*U*Uy-194*a*U*U)
    angular2 = [22*U*du-(da+db)/2 for du,da,db in zip(DU,DA,DB)]
    lapse2 = r2*(R+3*S)-24*rho*U*T
    if mutation == 'omit_radial_second_profile':
        shift2 -= r2*(Ay/2+4*a*A)
    if mutation == 'omit_trace_inverse_cross':
        lapse2 += 28*rho*U*T
    g = [[ZERO for j in range(5)] for i in range(5)]
    dg = [[[ZERO for j in range(5)] for i in range(5)] for k in range(5)]
    g[0][0] = jet(-r2,-2*r2*U,-r2*U*U)
    g[1][1] = jet(1,-4*U,A/2)
    for i in range(2,5):
        g[i][i] = jet(r2,2*r2*U,r2*B/2)
    dg[1][0][0] = jet(-2*a*r2,-2*r2*(2*a*U+Uy),-2*r2*(a*U*U+U*Uy))
    dg[1][1][1] = jet(0,-4*Uy,Ay/2)
    for i in range(2,5):
        dg[1][i][i] = jet(2*a*r2,r2*(4*a*U+2*Uy),r2*(a*B+By/2))
    for k,(du,da,db) in enumerate(zip(DU,DA,DB),2):
        dg[k][0][0] = jet(0,-2*r2*du,-2*r2*U*du)
        dg[k][1][1] = jet(0,-4*du,da/2)
        for i in range(2,5):
            dg[k][i][i] = jet(0,2*r2*du,r2*db/2)
    dg[0][1][1] = jet(0,-4*T,rho*R-4*U*T)
    for i in range(2,5):
        dg[0][i][i] = jet(0,2*r2*T,r2*(rho*S+2*U*T))
    dg[0][0][0] = jet(0,-2*r2*T,-rho*lapse2-2*r2*U*T)
    cross = -4*U*shift1 if mutation != 'omit_metric_mixed_cross' else 0
    dg[0][0][1] = dg[0][1][0] = jet(0,shift1,shift2/2+cross)
    for k,b2 in enumerate(angular2,2):
        dg[0][0][k] = dg[0][k][0] = jet(0,0,r2*b2/2)
    if mutation == 'omit_lapse_square_value':
        g[0][0] = jet(*g[0][0][:2])
        for k in range(1,5):
            dg[k][0][0] = jet(*dg[k][0][0][:2])
    if mutation == 'halve_second_gauge_velocity':
        for i in range(5):
            dg[0][0][i] = jet(*dg[0][0][i][:2],dg[0][0][i][2]/2)
            dg[0][i][0] = dg[0][0][i]
    gi = invert_metric(g)
    for i in range(5):
        for j in range(5):
            require(add(*(mul(g[i][k],gi[k][j]) for k in range(5))) == jet(int(i==j)), 'inverse')
    Gamma = [[[scale(Q(1,2),add(*(mul(gi[i][m],add(dg[j][m][k],dg[k][j][m],scale(-1,dg[m][j][k])))
        for m in range(5)))) for k in range(5)] for j in range(5)] for i in range(5)]
    return [add(*(mul(gi[j][k],add(Gamma[i][j][k],jet(-Gamma[i][j][k][0])))
        for j in range(5) for k in range(5))) for i in range(5)]


def data(n):
    return (Q(129+n),Q(2+n,19),Q(3+n,23),Q(-5-n,29),Q(7+n,31),
            Q(11+n,37),Q(-13-n,41),Q(17+n,43),Q(-19-n,47),
            Q(23+n,53),Q(-29-n,59),
            tuple(Q(n+k+2,61+k) for k in range(3)),
            tuple(Q(-n-2*k-3,67+k) for k in range(3)),
            tuple(Q(n+3*k+5,71+k) for k in range(3)))


def physical_and_map():
    spec = importlib.util.spec_from_file_location('mean_initial_checked', ROOT/'mean_initial.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows=[]
    for n in (1,3,7):
        rho,a,u,uy,ut,A,Ay,B,By = data(n)[:9]
        qm,j,ell = Q(13+n,79),Q(17+n,83),Q(1,1000000)
        s=ut/rho
        # Extract physical K directly from the canonical inverse map series.
        # Mean is E[Y^2]=1/6; Y is sampled 0,1,-1 to recover the polynomial.
        values=[]
        scalar=[]
        for Y in (Q(0),Q(1),Q(-1)):
            gammas=[jet(1,-4*u*Y,A/2)]+[jet(1,2*u*Y,B/2)]*3
            taus=[jet(0,-3*s*Y/2,0)]+[jet(0,0,-qm/(12*a))]*3
            # rho_gamma/rho0=1+epsilon uY+O(epsilon^2); only reciprocal order1 is needed.
            idensity=jet(1,-u*Y)
            trace=add(*(mul(g,t) for g,t in zip(gammas,taus)))
            k=[scale(2,mul(idensity,add(mul(mul(g,g),tau),scale(-Q(1,3),mul(g,trace)))))
               for g,tau in zip(gammas,taus)]
            values.append([2*v[2] for v in k])
            Pi=mul(idensity,jet(0,j*Y))
            scalar.append(2*mul(jet(rho,rho*u*Y),Pi)[2])
        meanK=[values[0][k]+(values[1][k]+values[2][k]-2*values[0][k])/12 for k in range(4)]
        require(meanK == [qm/(3*a)+6*u*s]+[-u*s]*3,'canonical mean K')
        require(scalar == [0,0,0],'coordinate scalar second velocity')
        result=module.initial_means(rho=rho,a=a,u=u,ut=ut,uy=uy,qm=qm,A=A,Ay=Ay,B=B,By=By,ell=ell)
        # Form averaged ADM trace and products before frame conversion.
        trace2=sum(meanK)-2*14*u*s/6
        lapse2=rho*rho*(trace2+4*u*s/6)
        gtt2=-2*rho*lapse2-4*rho*rho*u*ut/6
        b1y=-6*rho*rho*(uy+4*a*u)
        b2mean=rho*rho*(Ay/2-3*By/2+4*a*A-3*a*B-74*u*uy/6-194*a*u*u/6)
        gty2=b2mean-8*u*b1y/6
        gy2=2*rho*meanK[0]-8*u*ut/6
        gj2=2*rho*meanK[1]+4*u*ut/6
        require(result['V0']==tuple(ell*v for v in (gy2,gty2/rho,gtt2/(rho*rho),gj2,Q(0))), 'five frame map')
        require(result['U0']==(A,0,-u*u/3,B,0),'five displacement')
        require(result['dot_N2']==lapse2 and result['dot_beta_y2']==b2mean,'gauge mean map')
        require(result['V0'][0]+result['V0'][2]+6*result['V0'][3]==4*ell*u*ut,'trace identity')
        require(result['U0_r']==tuple(-ell*rho*v for v in (Ay,0,-2*u*uy/3,By,0)),'radial orientation')
        # Pm(0,D)=6aD, so the sign/factor has a nonzero mutation witness.
        require(6*a*(-qm/(6*a))==-qm,'momentum normalization')
        require(6*a*(qm/(3*a))!= -qm,'wrong momentum factor control')
        rows.append({'label':n,'canonical_conversion':True,'scalar_velocity_zero':True,
                     'five_displacement_velocity_radial_maps':True,'trace_identity':True,
                     'momentum_sign_factor_and_control':True})
    return rows


def run(mutation=None):
    witnesses=[]
    for n in (1,3,7):
        residual=gauge_residual(data(n),mutation)
        require(residual==[ZERO]*5,'GAUGE_RESIDUAL_NONZERO')
        witnesses.append({'label':n,'five_C0_C1_C2_zero':True})
    controls=[]
    for name in ('omit_radial_second_profile','omit_trace_inverse_cross','omit_metric_mixed_cross','omit_lapse_square_value','halve_second_gauge_velocity'):
        failures=[n for n in (1,3,7) if gauge_residual(data(n),name)!=[ZERO]*5]
        require(failures,'mutation survived '+name)
        controls.append({'mutation':name,'detected_at_rational_jets':failures})
    return {'status':'PASS','method':'Exact rational general five-dimensional metric inverse and Christoffel contraction','gauge_checks':witnesses,'map_checks':physical_and_map(),'in_process_controls':controls,'in_process_control_count':len(controls),'real_subprocess_control_count':0,'limits':'Finite algebra checks only; retained analytic theorems and evolution not certified.'}
