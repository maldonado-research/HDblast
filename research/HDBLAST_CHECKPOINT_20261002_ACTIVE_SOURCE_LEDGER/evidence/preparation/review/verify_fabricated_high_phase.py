#!/usr/bin/env python3
"""Independent high-phase kernel checks using fabricated inputs, never B jets."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import mpmath as mp
import numpy as np


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def represented(ctx,value):
    n,d=value.as_integer_ratio()
    return ctx.mpf(n)/ctx.mpf(d)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--primary-source',type=Path,required=True)
    p.add_argument('--independent-kernel',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise RuntimeError('Fresh fabricated audit receipt required')
    primary=load(args.primary_source,'review_primary_fabricated_only')
    independent=load(args.independent_kernel,'review_independent_fabricated_only')
    LD,CD=np.longdouble,np.clongdouble
    # Explicitly artificial values, not any retained study modes or nodes.
    k=np.array([LD('0.125'),LD('37.5'),LD(12345)/LD(49)],dtype=LD)
    dt=LD(1)/128
    a=LD('-8.5')
    epsilon=LD('0.0001')
    force=LD('0.00073')
    u0=np.array([CD('0.17')+CD('0.23j')]*3,dtype=CD)
    w0=np.array([CD('0.03')+CD('0.08j')]*3,dtype=CD)
    rows=[]
    for dps in (80,100):
        ctx=mp.mp.clone();ctx.dps=dps
        em=represented(ctx,epsilon);fm=represented(ctx,force)
        am=represented(ctx,a);dm=represented(ctx,dt)
        def exact(ki,t):
            km=represented(ctx,ki);omega=2*km
            ua=ctx.mpc(represented(ctx,u0[0].real),represented(ctx,u0[0].imag))
            wa=ctx.mpc(represented(ctx,w0[0].real),represented(ctx,w0[0].imag))
            E=ctx.exp(1j*omega*t)
            phi=ctx.expm1(1j*omega*t)/(1j*omega)
            return ua+phi*wa-em*fm*(phi-t)/(1j*omega),E*wa-em*fm*phi
        def G(ki,t):
            km=represented(ctx,ki);u,w=exact(ki,t);L=-1/(am+t)
            return (3*L*L*u.real-L*w.real)/(2*km*em)
        for order in (24,32):
            plan=primary.phase_plan(k,dt,order)
            ue,we=primary.phase_step(u0,w0,plan,np.full(order,force,dtype=LD))
            gaps=[]
            for j,ki in enumerate(k):
                ut,wt=exact(ki,dm)
                un=ctx.mpc(represented(ctx,ue[j].real),represented(ctx,ue[j].imag))
                wn=ctx.mpc(represented(ctx,we[j].real),represented(ctx,we[j].imag))
                gaps.append(max(abs(un-ut),abs(wn-wt)))
            largest=max(gaps)
            if largest>ctx.mpf('1e-15'):raise RuntimeError('Primary fabricated exact Duhamel mismatch')
            rows.append({'route':'primary','order':order,'decimal_context':dps,
                         'state_gap':ctx.nstr(largest,35),'state_gate':'1e-15'})
        for source_order,ledger_order in ((16,16),(16,24),(24,24)):
            kernel=independent.DenseKernel(k,dt,source_order,ledger_order)
            ue,we=kernel.step(u0,w0,np.full(source_order,force,dtype=LD),epsilon)
            ud,wd=kernel.dense(u0,w0,np.full((ledger_order,source_order),force,dtype=LD),epsilon)
            state=[];ledger=[]
            for j,ki in enumerate(k):
                km=represented(ctx,ki)
                ut,wt=exact(ki,dm)
                un=ctx.mpc(represented(ctx,ue[j].real),represented(ctx,ue[j].imag))
                wn=ctx.mpc(represented(ctx,we[j].real),represented(ctx,we[j].imag))
                state.append(max(abs(un-ut),abs(wn-wt)))
                integral=ctx.zero
                for q,offset in enumerate(kernel.t):
                    tm=represented(ctx,offset);L=-1/(am+tm)
                    U=ctx.mpc(represented(ctx,ud[j,q].real),represented(ctx,ud[j,q].imag))/em
                    W=ctx.mpc(represented(ctx,wd[j,q].real),represented(ctx,wd[j,q].imag))/em
                    # R and P independently use the direct bare operators.
                    R=((2*km*km+3*L*L)*U.real-km*W.imag-L*W.real)/(2*km)
                    P=((2*km*km/3-L*L)*U.real-km*W.imag-L*W.real)/(2*km)
                    integral+=represented(ctx,kernel.tw[q])*L*(R-3*P)
                source_work=fm*ctx.log(abs(am)/abs(am+dm))/(2*km)
                exact_integral=G(ki,dm)-G(ki,ctx.zero)-source_work
                ledger.append(abs(integral-exact_integral))
            sg,lg=max(state),max(ledger)
            if sg>ctx.mpf('1e-15'):raise RuntimeError('Independent fabricated Duhamel mismatch')
            if lg>ctx.mpf('1e-9'):raise RuntimeError('Independent fabricated direct-F phase integral mismatch')
            rows.append({'route':'independent','source_order':source_order,'ledger_order':ledger_order,'decimal_context':dps,
                         'state_gap':ctx.nstr(sg,35),'state_gate':'1e-15',
                         'direct_F_gap':ctx.nstr(lg,35),'direct_F_gate':'1e-9'})
    result={'status':'PASS_FABRICATED_HIGH_PHASE_REVIEW','physical_evaluations':0,
            'physical_arrays_loaded':0,'study_source_invocations':0,'rows':rows,
            'primary_sha256':hashlib.sha256(args.primary_source.read_bytes()).hexdigest(),
            'independent_kernel_sha256':hashlib.sha256(args.independent_kernel.read_bytes()).hexdigest(),
            'review_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python_optimization':sys.flags.optimize,'precision_scope':'Native LD kernels versus independently computed MP exact constant-force Duhamel truth; fabricated study-independent state and momentum nodes.'}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('status','physical_evaluations','study_source_invocations')}))


if __name__=='__main__':main()
