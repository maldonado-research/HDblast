#!/usr/bin/env python3
"""Prospective phase-aware active Ward diagnostic; import evaluates no sources."""
from __future__ import annotations
import argparse
import ast
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import resource
import sys
import time
import mpmath
import numpy as np

LD, CD = np.longdouble, np.clongdouble
SOURCES=('positive_B','signed_uB')
SETTINGS=('coarse','fine')
CUTOFFS=(64,128,256)
QUANTITIES=('q','q_prime','q_second','rho','p','Q0','rho0','p0','current')
MEMBERS=('k','momentum_weights','u_1','w_1','u_2','w_2','u_3','w_3',
         'observation_eta','history_eta','history_values','history_ledger_integrand',
         'history_baseline_contact','history_source_jet','history_forcing_jet','history_geometry',
         'history_cutoffs','history_quantity_names','history_geometry_names')
ORDERS=(24,32)
LEVELS=((80,266),(100,333))
EPS_RATIO=(3777893186295716171,37778931862957161709568)
PI_RATIO=(14488038916154245685,4611686018427387904)
AUTHORIZED=False


def require(ok,message):
    if not ok: raise ValueError(message)


def budget(started):
    require(time.monotonic()-started<=900,'Entire primary route exceeds900 seconds')
    require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<=262144,'Primary peakRSS exceeds262144KiB')


def represented(ctx,x):
    require(isinstance(x,LD) and np.isfinite(x),'Exact native binary80 input required')
    n,d=x.as_integer_ratio()
    require(d>0 and d&(d-1)==0,'Nonbinary represented denominator')
    return ctx.mpf(n)/ctx.mpf(d)


def complex_represented(ctx,z):
    require(isinstance(z,CD),'Exact native complex160 input required')
    return ctx.mpc(represented(ctx,z.real),represented(ctx,z.imag))


def round_shift(n,p):
    if p<=0:return n<<(-p)
    q=n>>p
    r=n-(q<<p)
    return q+int(r>(1<<(p-1)) or (r==(1<<(p-1)) and q&1))


def quantize(x,p):
    sign,man,exp,_=x._mpf_
    return round_shift(-man if sign else man,-exp-p)


def gauss_ld(n):
    """Newton Legendre roots/weights generated natively in binary80."""
    pi=np.arccos(LD(-1))
    x=np.cos(pi*(np.arange(n,dtype=LD)+LD('0.75'))/(LD(n)+LD('0.5')))
    for iteration in range(20):
        p0=np.ones(n,dtype=LD);p1=x.copy()
        for j in range(2,n+1):
            p0,p1=p1,((2*j-1)*x*p1-(j-1)*p0)/j
        dp=n*(x*p1-p0)/(x*x-1)
        increment=p1/dp
        x-=increment
        if np.max(np.abs(increment))<=4*np.finfo(LD).eps:break
    else:raise ValueError('Binary80 Legendre nodes failed to converge')
    # Recompute derivative at the returned roots.
    p0=np.ones(n,dtype=LD);p1=x.copy()
    for j in range(2,n+1):p0,p1=p1,((2*j-1)*x*p1-(j-1)*p0)/j
    dp=n*(x*p1-p0)/(x*x-1)
    weights=2/((1-x*x)*dp*dp)
    order=np.argsort(x)
    return x[order],weights[order]


def bump_polynomials():
    poly=[1];result=[]
    for n in range(6):
        result.append(tuple(poly));new=[0]*(len(poly)+4)
        for j,c in enumerate(poly):
            if j:
                for shift,fac in ((0,1),(2,-2),(4,1)):new[j-1+shift]+=j*c*fac
            new[j+1]+=(4*n-2)*c;new[j+3]-=4*n*c
        while len(new)>1 and new[-1]==0:new.pop()
        poly=new
    return tuple(result)


POLYNOMIALS=bump_polynomials()


def source_jet(times,source):
    require(AUTHORIZED,'Physical source evaluation requires published registration')
    require(source in SOURCES,'Unexpected source')
    times=np.asarray(times,dtype=LD);u=times+4;inside=np.abs(u)<1
    z=u[inside];den=1-z*z;b=np.exp(1-1/den)
    out=np.zeros((6,)+times.shape,dtype=LD)
    for j,poly in enumerate(POLYNOMIALS):
        value=np.zeros_like(z)
        for c in reversed(poly):value=value*z+c
        out[j][inside]=b*value/den**(2*j)
    if source=='signed_uB':
        original=out.copy()
        for j in range(6):out[j]=u*original[j]+(j*original[j-1] if j else 0)
    return out


def forcing_jet(times,hjet):
    """Pure canonical g^(0..3), valid for arbitrary consistent source jets."""
    L=-LD(1)/np.asarray(times,dtype=LD)
    return np.stack([sum(math.comb(n,j)*(4*math.factorial(j+1)*L**(j+2)*hjet[n-j]
                        -2*math.factorial(j)*L**(j+1)*hjet[n-j+1])
                        for j in range(n+1))-hjet[n+2] for n in range(4)])


class IntegerAccumulator:
    """Exact integer products/sums of represented modes and MP-seeded weights."""
    def __init__(self,k,weight,dps,p):
        self.ctx=mpmath.mp.clone();self.ctx.dps=dps;self.p=p
        ctx=self.ctx;eps=ctx.mpf(EPS_RATIO[0])/EPS_RATIO[1];pi=ctx.mpf(PI_RATIO[0])/PI_RATIO[1]
        self.coefficients=[];self.M=[ctx.zero]*3
        for ki,wi in zip(k,weight):
            km,wm=represented(ctx,ki),represented(ctx,wi)
            b=0 if ki<64 else (1 if ki<128 else 2)
            coefC=wm*km/(4*pi*pi*eps)
            coefE=wm*km**3/(2*pi*pi*eps)
            coefJ=wm*km**2/(4*pi*pi*eps)
            self.coefficients.append((b,quantize(coefC,p),quantize(coefE,p),quantize(coefJ,p),quantize(wm*km*km/(2*pi*pi),p)))
            self.M[b]+=wm*km/(4*pi*pi)
        self.M=[sum(self.M[:j+1]) for j in range(3)]
        self.quantization_accumulation_error_bound=ctx.zero

    def moments(self,u,w):
        bins=[[0,0,0,0] for _ in range(3)]
        max_abs=self.ctx.zero
        for ui,wi,coeffs in zip(u,w,self.coefficients):
            b,cc,ce,cj,cm=coeffs
            un,ud=ui.real.as_integer_ratio();wn,wd=wi.real.as_integer_ratio();jn,jd=wi.imag.as_integer_ratio()
            row=bins[b]
            row[0]+=round_shift(cc*un,ud.bit_length()-1)
            row[1]+=round_shift(cc*wn,wd.bit_length()-1)
            row[2]+=round_shift(ce*un,ud.bit_length()-1)
            row[3]+=round_shift(cj*jn,jd.bit_length()-1)
        ctx=self.ctx;den=1<<self.p;sumints=[0,0,0,0];out=[]
        for row in bins:
            sumints=[a+b for a,b in zip(sumints,row)]
            out.append([ctx.mpf(v)/den for v in sumints])
        # Coefficient quantization + product rounding only; phase/source arithmetic excluded.
        native_magnitude=max(np.max(np.abs(u.real)),np.max(np.abs(w.real)),np.max(np.abs(w.imag)))
        bound=ctx.mpf(len(u))*(represented(ctx,LD(native_magnitude))+1)/(2*den)
        self.quantization_accumulation_error_bound=max(self.quantization_accumulation_error_bound,bound)
        return out


def reduce_discrete(acc,arrays):
    bins=[[0]*len(arrays) for _ in range(3)]
    for i,coeffs in enumerate(acc.coefficients):
        b,cc,ce,cj,cm=coeffs
        for j,array in enumerate(arrays):
            n,d=array[i].as_integer_ratio();bins[b][j]+=round_shift(cm*n,d.bit_length()-1)
    cumulative=[0]*len(arrays);out=[]
    for row in bins:
        cumulative=[x+y for x,y in zip(cumulative,row)]
        out.append([acc.ctx.mpf(v)/(1<<acc.p) for v in cumulative])
    return out


def pure_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,'Pure reference module unavailable')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def discrete_action_contacts(k,eta,h,wkb,baseline):
    L=-LD(1)/eta;M=2*L*L
    mass=[math.factorial(n+1)*M*L**n for n in range(5)]
    dm=[2*sum(math.comb(n,j)*mass[j]*h[n-j] for j in range(n+1)) for n in range(5)]
    omega=[np.sqrt(k*k+M)];delta=[dm[0]/(2*omega[0])]
    for n in range(1,5):
        omega.append((mass[n]-sum(math.comb(n,j)*omega[j]*omega[n-j] for j in range(1,n)))/(2*omega[0]))
        delta.append((dm[n]-sum(math.comb(n,j)*(delta[j]*omega[n-j]+omega[j]*delta[n-j])
                      for j in range(1,n))-2*delta[0]*omega[n])/(2*omega[0]))
    C=(2*L*L,4*L**3,12*L**4)
    dC=(h[2]+2*L*h[1],h[3]+2*L*h[2]+2*L*L*h[1],
        h[4]+2*L*h[3]+4*L*L*h[2]+4*L**3*h[1])
    sub=wkb.evaluate_wkb(k*k,L,M,omega,delta,C,dC,h[1],2*M*h[0])
    _,R0,P0=baseline.baseline_differences(k,omega[0],L)
    CR=(L*h[1]+2*L*L*h[0])/(2*k)-sub['deltaR']-4*h[0]*R0
    CP=(L*h[1]-2*L*L*h[0])/(2*k)-sub['deltaP']-4*h[0]*P0
    return (CR,CP,R0/L**4,P0/L**4)


def signed_flow_and_triangle(acc,k,weights,pred_u,pred_w,observed_u,observed_w,eta):
    ctx=acc.ctx;eps=ctx.mpf(EPS_RATIO[0])/EPS_RATIO[1];pi=ctx.mpf(PI_RATIO[0])/PI_RATIO[1]
    L=-1/represented(ctx,eta);bins=[[ctx.zero,ctx.zero] for _ in CUTOFFS]
    for i,kn in enumerate(k):
        km=represented(ctx,kn);wm=represented(ctx,weights[i]);mu=wm*km*km/(2*pi*pi)
        du=complex_represented(ctx,observed_u[i])-complex_represented(ctx,pred_u[i])
        dw=complex_represented(ctx,observed_w[i])-complex_represented(ctx,pred_w[i])
        b=0 if kn<64 else (1 if kn<128 else 2)
        bins[b][0]+=mu*((2*km*km+3*L*L)*du.real-L*dw.real-km*dw.imag)/(2*km*eps)
        bins[b][1]+=mu*((2*km*km+3*L*L)*abs(du)+(L+km)*abs(dw))/(2*km*eps)
    return [[sum(bins[b][j] for b in range(ci+1)) for j in (0,1)] for ci in range(3)]


def phase_plan(k,dt,order):
    nodes,weights=gauss_ld(order)
    offset=(nodes+1)*dt/2;local_weight=weights*dt/2
    ik=CD(2j)*k
    rotation=np.exp(ik*dt);drift=np.expm1(ik*dt)/ik
    exponent=ik[:,None]*(dt-offset)[None,:]
    force_w=np.exp(exponent)*local_weight[None,:]
    force_u=np.expm1(exponent)/ik[:,None]*local_weight[None,:]
    return {'offset':offset,'weights':local_weight,'rotation':rotation,'drift':drift,
            'force_w':force_w,'force_u':force_u,'order':order}


def phase_step(u,w,plan,g):
    require(g.dtype==np.dtype(LD),'Forcing must be native binary80')
    eps=LD('0.0001')
    return (u+plan['drift']*w-eps*np.einsum('ij,j->i',plan['force_u'],g,optimize=False),
            plan['rotation']*w-eps*np.einsum('ij,j->i',plan['force_w'],g,optimize=False))


class RationalContacts:
    """Execute frozen action contacts with every literal division in MP."""
    def __init__(self,path,ctx):
        class ToMP(ast.NodeTransformer):
            def visit_Constant(self,node):
                if type(node.value) is int:
                    return ast.copy_location(ast.Call(func=ast.Attribute(value=ast.Name(id='ctx',ctx=ast.Load()),
                         attr='mpf',ctx=ast.Load()),args=[ast.Constant(str(node.value))],keywords=[]),node)
                return node
            def visit_Attribute(self,node):
                if isinstance(node.value,ast.Name) and node.value.id=='math' and node.attr=='pi':
                    return ast.copy_location(ast.Name(id='pi',ctx=ast.Load()),node)
                return self.generic_visit(node)
        tree=ast.parse(Path(path).read_text());tree=ast.fix_missing_locations(ToMP().visit(tree))
        space={'ctx':ctx,'pi':ctx.mpf(PI_RATIO[0])/PI_RATIO[1]}
        exec(compile(tree,str(path),'exec'),space)
        self.evaluate=space['local_coefficients'];self.ctx=ctx

    def full(self,eta,K,h,g,M):
        ctx=self.ctx;L=-1/eta;pi=ctx.mpf(PI_RATIO[0])/PI_RATIO[1]
        v=K/ctx.sqrt(K*K+2*L*L);contact=self.evaluate(L,v,h)
        A=ctx.asinh(K/(ctx.sqrt(2)*L))-v
        Ap=-L*v**3;App=L*L*v**3*(2-3*v*v);pref=1/(8*pi*pi)
        qsub=pref*g[0]*A
        qsubp=pref*(g[1]*A+g[0]*Ap)
        qsubpp=pref*(g[2]*A+2*g[1]*Ap+g[0]*App)
        J5=v**3/6;J7=(v**3/3-v**5/5)/4;J9=(v**3/3-2*v**5/5+v**7/7)/8
        density=(3*L*L*g[0]-L*g[1])*v**3/(96*pi*pi)
        pressure=(70*L*L*g[0]*J9-(30*L*L*g[0]+10*L*g[1])*J7
                  +(g[2]-L*g[1]-9*L*L*g[0])*J5)/(48*pi*pi)
        CR=(3*L*L*qsub-L*qsubp)/2+L*L*contact['Q0']*g[0]/2+density+contact['rho_local']
        CP=(qsubpp-3*L*qsubp-3*L*L*qsub)/6-L*L*contact['Q0']*g[0]/6+pressure+contact['p_local']-(K*K/(8*pi*pi))*g[0]/3
        baseline=3*h[1]*L**4*(contact['rho0']+contact['p0'])
        CF=L*(CR-3*CP)-baseline
        return {'R':CR,'P':CP,'F':CF,'baseline_contact':baseline,'Lg':L*g[0],
                'rho0':contact['rho0'],'p0':contact['p0']}


class NativeContext:
    mpf=staticmethod(LD)
    sqrt=staticmethod(lambda value:np.sqrt(np.asarray(value,dtype=LD)))
    asinh=staticmethod(lambda value:np.arcsinh(np.asarray(value,dtype=LD)))


def direct_stresses(ctx,moments,time,contact):
    C,B,E,J=moments;L=-1/time
    R=E+3*L*L*C-L*B-J+contact['R']
    P=E/3-L*L*C-L*B-J+contact['P']
    return {'R':R,'P':P,'F':L*(R-3*P)-contact['baseline_contact'],
            'G':3*L*L*C-L*B,'canonical_constant':E-J}


def serialize(ctx,x):
    if isinstance(x,list):return [serialize(ctx,v) for v in x]
    if isinstance(x,dict):return {k:serialize(ctx,v) for k,v in x.items()}
    if type(x) in (str,int,bool) or x is None:return x
    require(ctx.isfinite(x),'Nonfinite scientific scalar')
    text=ctx.nstr(x,n=ctx.dps+12,strip_zeros=False)
    sign,man,exp,_=x._mpf_;num=-man if sign else man
    exact=Fraction(num<<exp,1) if exp>=0 else Fraction(num,1<<(-exp))
    error=abs(exact-Fraction(Decimal(text)))
    require(error<=Fraction('1e-12'),'Scientific serialization exceeds1e-12')
    ctx.serialization_max_error=max(getattr(ctx,'serialization_max_error',Fraction(0)),error)
    return text


def simpson_to(history,end,dt):
    require(end%2==0,'Global Simpson endpoint must be even')
    return dt/3*(history[0]+history[end]+4*np.sum(history[1:end:2],axis=0,dtype=LD)
                 +2*np.sum(history[2:end:2],axis=0,dtype=LD))


def augment_saved(a,setting,ctx,row,ci,times,measured_end):
    dt=LD(1)/(128 if setting=='coarse' else 256)
    indices=np.flatnonzero((a['history_eta']>=LD('-4.5'))&(a['history_eta']<=LD('-3.5')))
    require(len(indices)==len(times) and np.array_equal(a['history_eta'][indices],times),
            'Inherited history does not exactly match fixed active interval')
    ia,ib=int(indices[0]),int(indices[-1]);rr,pp=QUANTITIES.index('rho'),QUANTITIES.index('p')
    hist=a['history_ledger_integrand']
    S=simpson_to(hist,ib,dt)-simpson_to(hist,ia,dt) # both prefixes/subtraction remain binary80
    reset=simpson_to(hist[ia:ib+1],ib-ia,dt)
    stored_R=[represented(ctx,v) for v in a['history_values'][indices,ci,rr]]
    stored_P=[represented(ctx,v) for v in a['history_values'][indices,ci,pp]]
    stored_F=[represented(ctx,v) for v in hist[indices,ci]]
    deltaR=represented(ctx,a['history_values'][ib,ci,rr])-represented(ctx,a['history_values'][ia,ci,rr])
    row.update({'stored_R':stored_R,'stored_P':stored_P,'stored_F':stored_F,
                'stored_baseline_contact':[represented(ctx,v) for v in a['history_baseline_contact'][indices,ci]],
                'DeltaR':deltaR,'S_ab':represented(ctx,S[ci]),'S_direct_reset':represented(ctx,reset[ci])})
    row['D_S']=deltaR-row['S_ab'];row['D_cont']=deltaR-row['I_ab'];row['E_Q']=row['I_ab']-row['S_ab']
    row['S_global_minus_reset']=row['S_ab']-row['S_direct_reset']
    row['decomposition_error']=row['D_S']-row['D_cont']-row['E_Q']
    row['contact_endpoint_defect']=(stored_R[-1]-measured_end['R'])-(stored_R[0]-row['R'][0])
    row['flow_operator_closure']=row['D_cont']-row['E_flow']
    row['flow_contact_decomposition_error']=row['D_cont']-row['E_flow']-row['contact_endpoint_defect']-row['predicted_Ward_residual']
    for name,stored in (('R',stored_R),('P',stored_P),('F',stored_F)):
        row[name+'_profile_max_error']=max(abs(x-y) for x,y in zip(row[name],stored))
    row['baseline_contact_profile_max_error']=max(abs(x-represented(ctx,y)) for x,y in
              zip(row['baseline_contact'],a['history_baseline_contact'][indices,ci]))
    if setting=='fine':
        double=simpson_to(hist[::2],ib//2,2*dt)-simpson_to(hist[::2],ia//2,2*dt)
        row['S_ab_double_global']=represented(ctx,double[ci]);row['S_native_minus_double']=row['S_ab']-row['S_ab_double_global']
    # Midpoint is fixed prospectively; it measures inherited trajectory discrepancy,
    # independent of the final endpoint and cannot be selected after seeing results.
    mid=len(times)//2
    row['midpoint_contact_plus_flow_density_gap']=stored_R[mid]-row['R'][mid]
    return row


def run_case(a,source,setting,contacts_path,started,force_function=None,discrete_modules=None):
    """All source callbacks guarded externally; this function also accepts fabricated jets."""
    dt=LD(1)/(128 if setting=='coarse' else 256)
    expected_count=129 if setting=='coarse' else 257
    times=LD('-4.5')+np.arange(expected_count,dtype=LD)*dt
    k=a['k'];weight=a['momentum_weights'];results=[]
    discrete=[]
    if discrete_modules is not None:
        wkb=pure_module(discrete_modules[0],'primary_registered_pointwise_wkb')
        baseline=pure_module(discrete_modules[1],'primary_registered_baseline')
        for ti in (0,len(times)//2,len(times)-1):
            hn=(source_jet(np.asarray([times[ti]],dtype=LD),source) if force_function is None
                else force_function(np.asarray([times[ti]],dtype=LD)))[:,0]
            discrete.append(discrete_action_contacts(k,times[ti],hn,wkb,baseline))
    for order in ORDERS:
        budget(started)
        plan=phase_plan(k,dt,order)
        accum=[IntegerAccumulator(k,weight,dps,p) for dps,p in LEVELS]
        contacts=[RationalContacts(contacts_path,ac.ctx) for ac in accum]
        native_contacts=RationalContacts(contacts_path,NativeContext())
        profiles=[[{'R':[],'P':[],'F':[],'G':[],'canonical_constant':[],'contact_R':[],
                    'contact_P':[],'contact_F':[],'baseline_contact':[],'rho0':[],'p0':[]} for _ in CUTOFFS] for _ in LEVELS]
        sums=[[[ac.ctx.zero,ac.ctx.zero] for _ in CUTOFFS] for ac in accum]
        u=a['u_1'].copy();w=a['w_1'].copy()
        native_h_history=[];native_g_history=[];mid_u=mid_w=None
        for step,eta in enumerate(times):
            if step:
                previous=eta-dt;local_times=previous+plan['offset']
                local_h=(source_jet(local_times,source) if force_function is None else force_function(local_times))
                local_g=forcing_jet(local_times,local_h)
                u,w=phase_step(u,w,plan,local_g[0])
                # Direct ledger contacts are evaluated in native binary80 at
                # actual analytic source nodes. MP80/100 covers exact-ratio
                # accumulation of these represented integrands, not their evaluation.
                local_contact=[native_contacts.full(local_times,LD(K),list(local_h),list(local_g),LD(0))
                               for K in CUTOFFS]
                for li,ac in enumerate(accum):
                    ctx=ac.ctx
                    for j,nt in enumerate(local_times):
                        mw=represented(ctx,LD(plan['weights'][j]))
                        for ci,K in enumerate(CUTOFFS):
                            sums[li][ci][0]+=mw*represented(ctx,LD(local_contact[ci]['Lg'][j]))
                            sums[li][ci][1]+=mw*represented(ctx,LD(local_contact[ci]['F'][j]))
            ht_native=(source_jet(np.asarray([eta],dtype=LD),source) if force_function is None
                       else force_function(np.asarray([eta],dtype=LD)))
            gt_native=forcing_jet(np.asarray([eta],dtype=LD),ht_native)
            native_h_history.append(ht_native[:,0].copy());native_g_history.append(gt_native[:,0].copy())
            if step==len(times)//2:mid_u,mid_w=u.copy(),w.copy()
            for li,(ac,ct) in enumerate(zip(accum,contacts)):
                ctx=ac.ctx;moments=ac.moments(u,w)
                ht=[represented(ctx,LD(v)) for v in ht_native[:,0]]
                gt=[represented(ctx,LD(v)) for v in gt_native[:,0]];mt=represented(ctx,LD(eta))
                for ci,K in enumerate(CUTOFFS):
                    c=ct.full(mt,ctx.mpf(K),ht,gt,ac.M[ci]);st=direct_stresses(ctx,moments[ci],mt,c)
                    for name,value in st.items():profiles[li][ci][name].append(value)
                    for target,name in (('contact_R','R'),('contact_P','P'),('contact_F','F'),('baseline_contact','baseline_contact'),('rho0','rho0'),('p0','p0')):
                        profiles[li][ci][target].append(c[name])
            if step%16==0:budget(started)
        levels={}
        for li,ac in enumerate(accum):
            ctx=ac.ctx;final_observed=ac.moments(a['u_3'],a['w_3']);rows=[]
            midpoint_observed=ac.moments(a['u_2'],a['w_2']) if 'u_2' in a else None
            discrete_sums=[reduce_discrete(ac,d) for d in discrete]
            end_flow=signed_flow_and_triangle(ac,k,weight,u,w,a['u_3'],a['w_3'],times[-1])
            mid_flow=signed_flow_and_triangle(ac,k,weight,mid_u,mid_w,a['u_2'],a['w_2'],times[len(times)//2]) if 'u_2' in a else None
            for ci,K in enumerate(CUTOFFS):
                row=profiles[li][ci]
                Ia=row['G'][-1]-row['G'][0]-ac.M[ci]*sums[li][ci][0]+sums[li][ci][1]
                # Independently expand the same operator on measured final modes;
                # no measured final density contributes to Ia.
                measured_end=direct_stresses(ctx,final_observed[ci],represented(ctx,times[-1]),
                    {'R':row['contact_R'][-1],'P':row['contact_P'][-1],
                     'baseline_contact':row['baseline_contact'][-1]})
                row.update({'K':K,'times':[represented(ctx,t) for t in times],'I_ab':Ia,'source_work_integral':sums[li][ci][0],
                            'full_contact_ledger_integral':sums[li][ci][1],
                            'E_flow':measured_end['R']-row['R'][-1],
                            'predicted_delta_R':row['R'][-1]-row['R'][0],
                            'predicted_Ward_residual':row['R'][-1]-row['R'][0]-Ia,
                            'canonical_constant_drift':row['canonical_constant'][-1]-row['canonical_constant'][0],
                            'discrete_source_moment':ac.M[ci],
                            'analytic_source_moment':ctx.mpf(K*K)/(8*(ctx.mpf(PI_RATIO[0])/PI_RATIO[1])**2),
                            'mixed_momentum_work':(ctx.mpf(K*K)/(8*(ctx.mpf(PI_RATIO[0])/PI_RATIO[1])**2)-ac.M[ci])*sums[li][ci][0],
                            'moment_quantization_only_bound':ac.quantization_accumulation_error_bound})
                row['E_flow_signed_projection']=end_flow[ci][0];row['E_flow_triangle_bound']=end_flow[ci][1]
                row['flow_exact_projection_gap']=end_flow[ci][0]-row['E_flow']
                row['primitive_contact_identity_gap']=Ia-(row['G'][-1]-row['G'][0]+row['contact_R'][-1]-row['contact_R'][0]+row['mixed_momentum_work'])
                if 'history_eta' in a:augment_saved(a,setting,ctx,row,ci,times,measured_end)
                if discrete_sums:
                    row['discrete_contact_R_anchors']=[d[ci][0] for d in discrete_sums]
                    row['discrete_contact_P_anchors']=[d[ci][1] for d in discrete_sums]
                    row['discrete_rho0_anchors']=[d[ci][2] for d in discrete_sums]
                    row['discrete_p0_anchors']=[d[ci][3] for d in discrete_sums]
                    dc=(discrete_sums[-1][ci][0]-row['contact_R'][-1])-(discrete_sums[0][ci][0]-row['contact_R'][0])
                    row['contact_momentum_endpoint_difference']=dc
                    row['E_momentum']=dc-row['mixed_momentum_work']
                    if 'history_eta' in a:
                        row['E_operator']=row['contact_endpoint_defect']-dc
                        row['source_primitive_residual']=row['predicted_Ward_residual']-row['canonical_constant_drift']+row['mixed_momentum_work']
                        row['E_reconstruction']=row['canonical_constant_drift']+row['source_primitive_residual']
                        row['complete_flow_momentum_operator_decomposition_error']=row['D_cont']-row['E_flow']-row['E_momentum']-row['E_operator']-row['E_reconstruction']
                        mid=len(times)//2
                        old_mid=direct_stresses(ctx,midpoint_observed[ci],represented(ctx,times[mid]),
                           {'R':row['contact_R'][mid],'P':row['contact_P'][mid],'baseline_contact':row['baseline_contact'][mid]})
                        row['midpoint_E_flow']=old_mid['R']-row['R'][mid]
                        row['midpoint_E_momentum']=discrete_sums[1][ci][0]-row['contact_R'][mid]
                        row['midpoint_E_operator']=row['stored_R'][mid]-old_mid['R']-row['midpoint_E_momentum']
                        row['midpoint_decomposition_error']=row['stored_R'][mid]-row['R'][mid]-row['midpoint_E_flow']-row['midpoint_E_momentum']-row['midpoint_E_operator']
                        row['midpoint_flow_triangle_bound']=mid_flow[ci][1]
                        row['midpoint_flow_exact_projection_gap']=mid_flow[ci][0]-row['midpoint_E_flow']
                        knot=(0,mid,len(times)-1)
                        old_R=(row['R'][0],old_mid['R'],measured_end['R'])
                        old_P=(row['P'][0],old_mid['P'],measured_end['P'])
                        row['discrete_operator_R_residual_anchors']=[row['stored_R'][ti]-old_R[j]-(discrete_sums[j][ci][0]-row['contact_R'][ti]) for j,ti in enumerate(knot)]
                        row['discrete_operator_P_residual_anchors']=[row['stored_P'][ti]-old_P[j]-(discrete_sums[j][ci][1]-row['contact_P'][ti]) for j,ti in enumerate(knot)]
                        row['E_operator_initial']=row['discrete_operator_R_residual_anchors'][0]
                        row['E_operator_end']=row['discrete_operator_R_residual_anchors'][-1]
                        row['analytic_rho0_anchors']=[row['rho0'][ti] for ti in knot]
                        row['analytic_p0_anchors']=[row['p0'][ti] for ti in knot]
                        row['discrete_contact_R_minus_analytic_anchors']=[discrete_sums[j][ci][0]-row['contact_R'][ti] for j,ti in enumerate(knot)]
                        row['discrete_contact_P_minus_analytic_anchors']=[discrete_sums[j][ci][1]-row['contact_P'][ti] for j,ti in enumerate(knot)]
                        indices=np.flatnonzero((a['history_eta']>=LD('-4.5'))&(a['history_eta']<=LD('-3.5')))
                        row['discrete_baseline_contact_residual_anchors']=[represented(ctx,a['history_baseline_contact'][indices[ti],ci])-3*represented(ctx,native_h_history[ti][1])*(-1/represented(ctx,times[ti]))**4*(discrete_sums[j][ci][2]+discrete_sums[j][ci][3]) for j,ti in enumerate(knot)]
                        row['source_jet_profile_max_error']=max(abs(represented(ctx,x)-represented(ctx,y)) for xs,ys in zip(native_h_history,a['history_source_jet'][indices]) for x,y in zip(xs,ys))
                        row['forcing_jet_profile_max_error']=max(abs(represented(ctx,x)-represented(ctx,y)) for xs,ys in zip(native_g_history,a['history_forcing_jet'][indices]) for x,y in zip(xs,ys))

                        row['I_A']=row['I_ab'];row['D_cont_A']=row['D_cont'];row['E_Q_A']=row['E_Q']
                        row['I_ab']=row['I_A']+row['E_momentum']
                        row['D_cont']=row['DeltaR']-row['I_ab'];row['E_Q']=row['I_ab']-row['S_ab']
                        row['matched_decomposition_error']=row['D_cont']-row['E_flow']-row['E_operator']-row['E_reconstruction']
                        row['matched_integral_correction']=row['E_momentum']
                        row['raw_analytic_decomposition_error']=row['decomposition_error']
                        row['decomposition_error']=row['D_S']-row['D_cont']-row['E_Q']
                rows.append(serialize(ctx,row))
            levels[str(ctx.dps)]=rows
        results.append({'quadrature_order':order,'levels':levels,
                        'serialization_max_error_by_dps':{str(ac.ctx.dps):str(getattr(ac.ctx,'serialization_max_error',Fraction(0))) for ac in accum}})
        del plan,accum,contacts
    return {'source':source,'setting':setting,'controls':results}


def synthetic_jets(times):
    """Fabricated polynomial h, unrelated to the registered compact source."""
    z=np.asarray(times,dtype=LD)+4
    return np.stack((LD('0.2')+LD('0.07')*z+LD('0.03')*z*z,
                     LD('0.07')+LD('0.06')*z,np.full(z.shape,LD('0.06')),
                     np.zeros_like(z),np.zeros_like(z),np.zeros_like(z)))


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(65536),b''):h.update(block)
    return h.hexdigest()


def authenticated_helper(root,registration_sha256,freeze_commit):
    require(isinstance(registration_sha256,str) and re.fullmatch(r'[0-9a-f]{64}',registration_sha256),
            'Explicit registration SHA256 required')
    require(isinstance(freeze_commit,str) and re.fullmatch(r'[0-9a-f]{40}',freeze_commit),
            'Explicit public freeze commit required')
    root=Path(root).absolute()
    for path in (root,*root.parents):require(not path.is_symlink(),'Checkpoint ancestry symlink rejected')
    registration=root/'FULL_REGISTRATION.json';helper=root/'code/active_integrity.py'
    entry=root/'code/diagnostic_primary.py'
    for path in (registration,helper,entry,root/'code'):require(not path.is_symlink(),'Bootstrap source symlink rejected')
    require(registration.is_file() and sha(registration)==registration_sha256,'Bootstrap registration hash mismatch')
    table=json.loads(registration.read_text())['files']
    for relative,path in (('code/active_integrity.py',helper),('code/diagnostic_primary.py',entry)):
        require(path.is_file() and sha(path)==table.get(relative),'Bootstrap registered source mismatch '+relative)
    require(Path(__file__).absolute()==entry,'Execute the registered primary entry point')
    module=pure_module(helper,'active_integrity')
    return module.verify_frozen(root,registration_sha256,freeze_commit)


def validate_arrays(a,setting):
    n=8192 if setting=='coarse' else 16384;count=577 if setting=='coarse' else 1153
    require(set(a)==set(MEMBERS),'Exactly19 declared capsule members required')
    for key in ('k','momentum_weights'):
        require(a[key].shape==(n,) and a[key].dtype==np.dtype(LD) and np.all(np.isfinite(a[key]))
                and np.all(a[key]>0),'Invalid node/weight schema')
    require(np.all(np.diff(a['k'])>0) and a['k'][-1]<256,'Invalid momentum ordering/range')
    for key in ('u_1','w_1','u_2','w_2','u_3','w_3'):
        require(a[key].shape==(n,) and a[key].dtype==np.dtype(CD) and np.all(np.isfinite(a[key])),
                'Invalid unconstrained mode precision/schema '+key)
    require(tuple(a['history_quantity_names'])==QUANTITIES,'History quantity labels differ')
    require(tuple(a['history_cutoffs'])==CUTOFFS,'Cutoff labels differ')
    require(tuple(a['observation_eta'])==(-5.5,-4.5,-4,-3.5,-2.5,-1.5),'Observation labels differ')
    require(a['history_eta'].shape==(count,) and a['history_eta'].dtype==np.dtype(LD), 'History precision/schema differs')
    dt=LD(1)/(128 if setting=='coarse' else 256)
    require(np.array_equal(a['history_eta'],LD(-6)+np.arange(count,dtype=LD)*dt),'History grid differs')
    shapes={'history_values':(count,3,9),'history_ledger_integrand':(count,3),
            'history_baseline_contact':(count,3),'history_source_jet':(count,6),
            'history_forcing_jet':(count,4),'history_geometry':(count,10)}
    for key,shape in shapes.items():
        require(a[key].shape==shape and a[key].dtype==np.dtype(LD) and np.all(np.isfinite(a[key])),
                'History payload precision/schema differs '+key)
    prefixes=(2048,4096,8192) if setting=='coarse' else (4096,8192,16384)
    require(tuple(int(np.count_nonzero(a['k']<cut)) for cut in CUTOFFS)==prefixes,'Momentum prefixes differ')


def recursive_gap(a,b):
    require(type(a)==type(b),'Scientific control schema differs')
    if isinstance(a,list):
        require(len(a)==len(b),'Scientific control length differs')
        return max((recursive_gap(x,y) for x,y in zip(a,b)),default=Fraction(0))
    if isinstance(a,dict):
        require(a.keys()==b.keys(),'Scientific control keys differ')
        return max((recursive_gap(a[k],b[k]) for k in a if k!='K'),default=Fraction(0))
    return abs(Fraction(Decimal(a))-Fraction(Decimal(b)))


def analyze_records(records):
    failures=[];fatal=[];attribution=[];precision_gap=Fraction(0);order_gap=Fraction(0)
    for case in records:
        controls=case['controls'];require(tuple(c['quadrature_order'] for c in controls)==ORDERS,'Fixed controls differ')
        for control in controls:
            gap=recursive_gap(control['levels']['80'],control['levels']['100']);precision_gap=max(precision_gap,gap)
            if gap>Fraction('1e-12'):fatal.append({'source':case['source'],'setting':case['setting'],'gate':'precision_gap','exact_gap':str(gap)})
            for dps in ('80','100'):
                for row in control['levels'][dps]:
                    for key in ('decomposition_error','flow_contact_decomposition_error','complete_flow_momentum_operator_decomposition_error','midpoint_decomposition_error','matched_decomposition_error','flow_exact_projection_gap','midpoint_flow_exact_projection_gap'):
                        if abs(Fraction(Decimal(row[key])))>Fraction('1e-12'):
                            fatal.append({'source':case['source'],'setting':case['setting'],'K':row['K'],
                                          'order':control['quadrature_order'],'dps':dps,'gate':key,'value':row[key]})
                    for key in ('R_profile_max_error','P_profile_max_error','F_profile_max_error','D_cont','E_flow','E_operator','E_momentum','E_reconstruction','midpoint_E_flow','midpoint_E_operator','midpoint_E_momentum','baseline_contact_profile_max_error','canonical_constant_drift','source_primitive_residual','primitive_contact_identity_gap','source_jet_profile_max_error','forcing_jet_profile_max_error'):
                        if abs(Fraction(Decimal(row[key])))>Fraction('2e-7'):
                            failures.append({'source':case['source'],'setting':case['setting'],'K':row['K'],
                                             'order':control['quadrature_order'],'dps':dps,'gate':key,'value':row[key]})
                    for key in ('discrete_operator_R_residual_anchors','discrete_operator_P_residual_anchors','discrete_baseline_contact_residual_anchors','discrete_contact_R_minus_analytic_anchors','discrete_contact_P_minus_analytic_anchors'):
                        if max(abs(Fraction(Decimal(x))) for x in row[key])>Fraction('2e-7'):
                            failures.append({'source':case['source'],'setting':case['setting'],'K':row['K'],
                                             'order':control['quadrature_order'],'dps':dps,'gate':key,'values':row[key]})
        for dps in ('80','100'):
            gap=recursive_gap(controls[0]['levels'][dps],controls[1]['levels'][dps]);order_gap=max(order_gap,gap)
            if gap>Fraction('2e-7'):failures.append({'source':case['source'],'setting':case['setting'],'dps':dps,'gate':'quadrature_control_gap','exact_gap':str(gap)})
    if not failures and not fatal:
        for case in records:
            for row in case['controls'][-1]['levels']['100']:
                ds=abs(Fraction(Decimal(row['D_S'])))
                if ds>Fraction('2e-6') and abs(Fraction(Decimal(row['D_cont'])))<=ds/10 and abs(Fraction(Decimal(row['E_Q'])))>=9*ds/10:
                    attribution.append({'source':case['source'],'setting':case['setting'],'K':row['K']})
    return {'classification':'CONSISTENCY_FAILURE' if failures or fatal else ('LEDGER_ERROR_DEMONSTRATED' if attribution else 'NO_GATE_SCALE_ATTRIBUTION'),
            'failures':failures,'fatal_failures':fatal,'attribution_cases':attribution,
            'precision_gap_max_exact_rational':str(precision_gap),'quadrature_control_gap_max_exact_rational':str(order_gap)}


def main():
    global AUTHORIZED
    parser=argparse.ArgumentParser();parser.add_argument('--synthetic',action='store_true')
    parser.add_argument('--contacts',type=Path);parser.add_argument('--wkb',type=Path);parser.add_argument('--baseline',type=Path)
    parser.add_argument('--output',type=Path);parser.add_argument('--checkpoint-root',type=Path)
    parser.add_argument('--registration-sha256');parser.add_argument('--freeze-commit');parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args();started=time.monotonic();records=[];executed_source_sha256=sha(__file__)
    require(np.finfo(LD).nmant==63 and np.dtype(LD).itemsize==16,'Native x87 binary80 platform required')
    require(LD('0.0001').as_integer_ratio()==EPS_RATIO and np.arccos(LD(-1)).as_integer_ratio()==PI_RATIO,'Represented constants differ')
    if args.synthetic:
        require(args.checkpoint_root is None and args.output_dir is None,'Synthetic path cannot load registered inputs')
        require(args.contacts is not None and args.wkb is not None and args.baseline is not None and args.output is not None,'Synthetic paths required')
        require(not args.output.exists(),'Output must be fresh')
        for source in SOURCES:
            for setting in SETTINGS:
                n=8192 if setting=='coarse' else 16384
                k=(np.arange(n,dtype=LD)+LD('0.5'))*LD(256)/n;weight=np.full(n,LD(256)/n,dtype=LD)
                a={'k':k,'momentum_weights':weight,'u_1':np.zeros(n,dtype=CD),'w_1':np.zeros(n,dtype=CD),
                   'u_3':np.zeros(n,dtype=CD),'w_3':np.zeros(n,dtype=CD)}
                count=577 if setting=='coarse' else 1153
                history_eta=LD(-6)+np.arange(count,dtype=LD)/(128 if setting=='coarse' else 256)
                a.update(u_2=np.zeros(n,dtype=CD),w_2=np.zeros(n,dtype=CD),history_eta=history_eta,
                   history_values=np.zeros((count,3,9),dtype=LD),history_ledger_integrand=np.zeros((count,3),dtype=LD),
                   history_baseline_contact=np.zeros((count,3),dtype=LD),history_source_jet=synthetic_jets(history_eta).T,
                   history_forcing_jet=forcing_jet(history_eta,synthetic_jets(history_eta)).T)
                records.append(run_case(a,source,setting,args.contacts,started,synthetic_jets,(args.wkb,args.baseline)))
                print(json.dumps({'event':'synthetic_case_completed','source':source,'setting':setting,'seconds':time.monotonic()-started}),flush=True)
        result={'status':'SYNTHETIC_FULL_SHAPE_ONLY','physical_arrays_loaded':False,'physical_source_evaluated':False,'records':records,'fabricated_analysis':analyze_records(records)}
    else:
        require(sys.flags.optimize==0 and sys.version_info[:3]==(3,12,14) and np.__version__=='2.2.6' and mpmath.__version__=='1.3.0','Pinned normal physical runtime required')
        require(args.checkpoint_root is not None and args.output_dir is not None and not args.output_dir.exists(),'Registered path and fresh output directory required')
        require(args.contacts is None and args.wkb is None and args.baseline is None and args.output is None,'Physical paths are fixed by checkpoint')
        auth=authenticated_helper(args.checkpoint_root,args.registration_sha256,args.freeze_commit)
        root=args.checkpoint_root;args.contacts=root/'reference_inputs/metric_contact_coefficients.py'
        args.wkb=root/'reference_inputs/metric_wkb.py';args.baseline=root/'reference_inputs/stable_baselines.py'
        manifest=json.loads((root/'inputs/INPUT_MANIFEST.json').read_text())
        require(manifest['scientific_calibration_status']=='FAIL_UNCHANGED','Earlier scientific FAIL must remain')
        AUTHORIZED=True;args.output_dir.mkdir();args.output=args.output_dir/'diagnostic.json'
        for source in SOURCES:
            for setting in SETTINGS:
                path=root/'inputs'/('metric_modes_'+source+'_'+setting+'.npz')
                with np.load(path,allow_pickle=False) as z:a={key:z[key] for key in MEMBERS}
                validate_arrays(a,setting)
                records.append(run_case(a,source,setting,args.contacts,started,None,(args.wkb,args.baseline)))
                with (args.output_dir/'progress.jsonl').open('a') as handle:
                    handle.write(json.dumps({'event':'case_completed','source':source,'setting':setting,'seconds':time.monotonic()-started})+'\n')
                del a
        result={'schema_version':1,'route':'primary_phase_green_primitive_integer_reductions',
                'status':'COMPLETED_SAVED_DATA_DIAGNOSTIC','records':records,'old_metric_status':'FAIL',
                'fixed_cases':12,'interval':['-4.5','-3.5'],'registration_sha256':args.registration_sha256,
                'freeze_commit':args.freeze_commit,'input_manifest_sha256':sha(root/'inputs/INPUT_MANIFEST.json'),
                **analyze_records(records)}
        authenticated_helper(root,args.registration_sha256,args.freeze_commit)
        AUTHORIZED=False
    budget(started)
    require(sha(__file__)==executed_source_sha256,'Primary source changed during execution')
    result['producer_sha256']=executed_source_sha256
    result['source_unchanged_during_execution']=True
    result['resources']={'seconds':time.monotonic()-started,'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         'scope':'entire primary route all12cases both24/32 controls both80/100 accumulation/contact levels','limit_seconds':900,'limit_peak_rss_kib':262144}
    result['arithmetic']={'phase_source':'native binary80/complex160','quadrature_contact_integrands':'native binary80 direct action contacts','decimal_precisions':[80,100],
       'fractional_bits':[266,333],'MP_scope':'profile/anchor rational contacts, exact-ratio direct native-contact quadrature accumulation and MP-seeded exact-ratio mode accumulation only',
       'quad_error':'fixed24/32 comparison empirical; no rigorous quadrature or total arithmetic certificate',
       'moment_quantization_bound_scope':'each of C/B/E/J accumulation only; excludes final R/P/I, MP seed error, phase/source arithmetic',
       'gap_universe':'every non-K scientific scalar/profile under controls.levels'}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    budget(started)
    print(json.dumps({'status':result['status'],'output':str(args.output),'resources':result['resources']},sort_keys=True))
    return 2 if result.get('fatal_failures') else 0


if __name__=='__main__':sys.exit(main())
