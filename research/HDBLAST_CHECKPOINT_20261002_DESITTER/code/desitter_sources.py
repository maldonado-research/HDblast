#!/usr/bin/env python3
"""Convergent, bounded four-point S4 determinant and common-action sources."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback
from functools import lru_cache
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
SOURCE_HEAD = '0205cc651bfb614c32e39dfe833d93d957264229'
INPUTS = [
 'research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/static_bridge/STATIC_DESITTER_SOURCE_BRIDGE.md',
 'research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/ACTION_AND_JUNCTION_SPECIFICATION.md',
 'research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/COMMON_ACTION_SMOOTH_FRW.md',
]
LOCAL = ['REGISTRATION.md', 'code/desitter_sources.py', 'code/validate_sources.py']

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def utc():
    return datetime.now(timezone.utc).isoformat()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, sort_keys=True)+'\n')

def freeze(repo):
    target = ROOT/'PREREGISTRATION.json'
    if target.exists() or any((ROOT/'outputs').iterdir()):
        raise RuntimeError('Refusing to overwrite registration or register after outputs exist')
    head = subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'], text=True).strip()
    if head != SOURCE_HEAD:
        raise RuntimeError(f'Expected source HEAD {SOURCE_HEAD}; got {head}')
    doc = {'registered_at_utc':utc(), 'source_commit':head,
           'grid':{'x_over_r':['1','2'],'H2_over_r':['0.5','1'],'r':'1'},
           'resolutions':[{'dps':50,'N':128},{'dps':80,'N':256}],
           'local_sha256':{f:sha(ROOT/f) for f in LOCAL},
           'source_sha256':{f:sha(repo/f) for f in INPUTS},
           'claim':'Four mathematical points only; no physical shell parameters or BVP'}
    write_json(target,doc)
    return doc

def verify_inputs(repo):
    p=ROOT/'PREREGISTRATION.json'
    doc=json.loads(p.read_text())
    for name,digest in doc['local_sha256'].items():
        if sha(ROOT/name) != digest:
            raise RuntimeError(f'Prospective local hash mismatch: {name}')
    for name,digest in doc['source_sha256'].items():
        if sha(repo/name) != digest:
            raise RuntimeError(f'Inherited source hash mismatch: {name}')
    return doc

@lru_cache(maxsize=4)
def coefficients(dps, N):
    with mp.workdps(dps):
        a=mp.mpf(3)/2
        z1=mp.zeta(-1,a)
        c0=2*mp.diff(lambda s:mp.zeta(s,a),-3)-mp.diff(lambda s:mp.zeta(s,a),-1)/2
        c1=-z1-mp.digamma(a)/4
        c2=mp.mpf(1)/4-mp.digamma(a)/2-mp.zeta(3,a)/8
        ak=tuple(mp.zeta(2*k-3,a)-mp.zeta(2*k-1,a)/4 for k in range(3,N+1))
        return c0,c1,c2,ak

def action_sources(x,y,r=1,*,dps=80,N=256):
    """No use of anomaly/trace relation in W, metric rho, or mass Q."""
    with mp.workdps(dps):
        x,y,r=map(mp.mpf,(x,y,r))
        if min(x,y,r)<=0:
            raise ValueError('Positive physical mass squared, curvature, and reference required')
        a=mp.mpf(3)/2
        q=x/y-mp.mpf(9)/4
        t=abs(q)/(a*a)
        if t>=1:
            raise ValueError('Outside registered convergent series domain')
        c0,c1,c2,ak=coefficients(dps,N)
        B0=mp.zeta(-3,a)-mp.zeta(-1,a)/4+q/8+q*q/4
        B0q=mp.mpf(1)/8+q/2
        B1=c0+c1*q+c2*q*q
        B1q=c1+2*c2*q
        for k,A in enumerate(ak,3):
            B1+=(-q)**k*A/k
            B1q+=(-1)**k*q**(k-1)*A
        P=B1-B0*mp.log(y)
        Pq=B1q-B0q*mp.log(y)
        A=x*x/2-2*x*y+mp.mpf(29)*y*y/15
        C=1/(16*mp.pi**2)
        W=-C*y*y*P+C*(r*x-r*r/4-2*r*y-A*mp.log(r))/2
        Wx=-C*y*Pq+C*(r-(x-2*y)*mp.log(r))/2
        Wy=-C*y*(2*P-(q+mp.mpf(9)/4)*Pq-B0)+C*(-2*r-(-2*x+mp.mpf(58)*y/15)*mp.log(r))/2
        rho=W-y*Wy/2
        Q=2*Wx
        common=a**3*(1+a/(2*N-2))*t**(N+1)/(1-t)
        eB=common/(N+1)
        eD=common/abs(q) if q else mp.mpf(0)
        eW=C*y*y*eB
        eQ=2*C*y*eD
        eWy=C*y*(2*eB+abs(q+mp.mpf(9)/4)*eD)
        erho=eW+y*eWy/2
        return {'W':W,'rho':rho,'p':-rho,'Q':Q,'Wx':Wx,'Wy':Wy,
                'q':q,'tail_bound':{'W':eW,'rho':erho,'Q':eQ}}

def digamma_Q(x,y,r=1):
    x,y,r=map(mp.mpf,(x,y,r));u=x/y;v=r/y
    nu=mp.sqrt(mp.mpf(9)/4-u)
    psi=mp.digamma(mp.mpf(3)/2+nu)+mp.digamma(mp.mpf(3)/2-nu)
    return mp.re(y*((u-2)*(psi-mp.log(v))-u+v+mp.mpf(4)/3)/(16*mp.pi**2))

def richardson(fun,at,h=mp.mpf('0.0001')):
    coarse=(fun(at+h)-fun(at-h))/(2*h)
    fine=(fun(at+h/2)-fun(at-h/2))/h
    return (4*fine-coarse)/3,abs(fine-coarse)

def encoded(obj):
    if isinstance(obj,dict): return {k:encoded(v) for k,v in obj.items()}
    if isinstance(obj,list): return [encoded(v) for v in obj]
    if isinstance(obj,(mp.mpf,mp.mpc)): return mp.nstr(obj,45)
    return obj

def benchmark(repo,out):
    protocol=verify_inputs(repo)
    out.mkdir(parents=True,exist_ok=False)
    write_json(out/'started.json',{'started_at_utc':utc(),'registration_sha256':sha(ROOT/'PREREGISTRATION.json'),
               'source_commit':SOURCE_HEAD,'observed_checkout_head':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
               'python':platform.python_version(),'mpmath':mp.__version__,'producer_sha256':sha(__file__)})
    gates=[];points=[]
    def gate(label,residual,scale=0,atol='1e-11',rtol='1e-9'):
        threshold=mp.mpf(atol)+mp.mpf(rtol)*abs(scale)
        ok=abs(residual)<=threshold
        gates.append({'name':label,'pass':bool(ok),'absolute_residual':abs(residual),'threshold':threshold})
        return ok
    try:
        mp.mp.dps=80
        for xs in protocol['grid']['x_over_r']:
            for ys in protocol['grid']['H2_over_r']:
                x,y=mp.mpf(xs),mp.mpf(ys);label=f'x{xs}_y{ys}'
                low=action_sources(x,y,dps=50,N=128)
                high=action_sources(x,y,dps=80,N=256)
                for resolution,values in [('primary',low),('refined',high)]:
                    for key in ['W','rho','Q']:
                        gate(f'{label}/{resolution}/{key}_tail',values['tail_bound'][key],atol='1e-14',rtol='0')
                for key in ['W','rho','Q']:
                    gate(f'{label}/{key}_precision_refinement',high[key]-low[key],high[key])
                dx,dxe=richardson(lambda z:action_sources(z,y)['W'],x)
                dy,dye=richardson(lambda z:action_sources(x,z)['W'],y)
                fdQ=2*dx;fdrho=high['W']-y*dy/2
                gate(f'{label}/mass_variation',fdQ-high['Q'],high['Q'],'1e-10','1e-8')
                gate(f'{label}/metric_variation',fdrho-high['rho'],high['rho'],'1e-10','1e-8')
                rhox,_=richardson(lambda z:action_sources(z,y)['rho'],x)
                Qy,_=richardson(lambda z:action_sources(x,z)['Q'],y)
                pairing=rhox-high['Q']/2+y*Qy/4
                gate(f'{label}/common_action_pairing',pairing,rhox,'1e-10','1e-8')
                independentQ=digamma_Q(x,y)
                gate(f'{label}/independent_digamma_Q',independentQ-high['Q'],high['Q'])
                h=x-1;c=h*h/2-2*h*y+mp.mpf(29)*y*y/15
                anomaly=c/(16*mp.pi**2)
                trace=-4*high['rho']+x*high['Q']-anomaly
                gate(f'{label}/independent_trace_check',trace,max(abs(4*high['rho']),abs(x*high['Q']),abs(anomaly)))
                points.append({'x_over_r':xs,'H2_over_r':ys,'r':'1',
                               'primary':low,'refined':high,
                               'independent_action_differences':{'Q':fdQ,'rho':fdrho,'Wx_step_change':dxe,'Wy_step_change':dye},
                               'digamma_Q':independentQ,'pairing_residual':pairing,'trace_residual':trace})
                print(label,'W=',mp.nstr(high['W'],16),'rho=',mp.nstr(high['rho'],16),'Q=',mp.nstr(high['Q'],16),flush=True)
        # Deliberate mutations must be rejected by the same source checks.
        p=points[0];h=p['refined'];x=mp.mpf(p['x_over_r']);y=mp.mpf(p['H2_over_r'])
        negatives=[]
        for name,wrong,correct in [('omit_metric_variation',h['W'],h['rho']),('flip_scalar_current',-h['Q'],h['Q'])]:
            residual=abs(wrong-correct);threshold=mp.mpf('1e-10')+mp.mpf('1e-8')*abs(correct)
            negatives.append({'name':name,'detected':bool(residual>threshold),'residual':residual,'threshold':threshold})
        # Local Euclidean W_L=V−12yF with V=x²/2,F=x.
        local_rhox=x-6*y;local_Q=2*x-24*y;local_Qy=-24
        gate('local_V_F_action_pairing',local_rhox-local_Q/2+y*local_Qy/4,atol='1e-30',rtol='0')
        omit_F_Q=2*x
        defect=local_rhox-omit_F_Q/2 # wrongly Qy=0
        negatives.append({'name':'omit_F_scalar_current','detected':bool(abs(defect)>mp.mpf('1e-10')),'residual':abs(defect),'threshold':mp.mpf('1e-10')})
        accepted=all(g['pass'] for g in gates) and all(g['detected'] for g in negatives)
        result={'completed_at_utc':utc(),'registration_sha256':sha(ROOT/'PREREGISTRATION.json'),
                'status':'PASS' if accepted else 'FAIL','points':points,'gates':gates,'negative_controls':negatives,
                'error_claim':'Explicit analytic zeta-series truncation bounds; arithmetic/derivative convergence checks are not interval certification.',
                'scope':'Only four declared mathematical points; no physical shell parameters, BVP, heating or dynamics.'}
        write_json(out/'results.json',encoded(result))
        print('terminal_status',result['status'],f'({len(gates)} gates; {len(negatives)} negative controls)',flush=True)
        return 0 if accepted else 1
    except Exception as exc:
        write_json(out/'failure.json',{'failed_at_utc':utc(),'error':str(exc),'traceback':traceback.format_exc(),'gates_so_far':encoded(gates),'completed_points':encoded(points)})
        raise

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,default=Path('/workspace/HDblast'))
    ap.add_argument('--freeze',action='store_true');ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    if args.freeze:
        print(json.dumps(freeze(args.repo),indent=2));return 0
    if not args.output:ap.error('--output required; every run requires a new directory')
    return benchmark(args.repo,args.output)

if __name__=='__main__':sys.exit(main())
