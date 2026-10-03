#!/usr/bin/env python3
"""Independent fabricated high-phase truth; never uses study source or arrays."""
import importlib.util
import json
from pathlib import Path
import sys
import mpmath
import numpy as np

path=Path(__file__).with_name('diagnostic_primary.py')
spec=importlib.util.spec_from_file_location('tested_primary',path)
P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)
checks=[]
for dps in (80,100):
    ctx=mpmath.mp.clone();ctx.dps=dps
    for order in P.ORDERS:
        k=np.asarray([P.LD('0.13'),P.LD('250.25'),P.LD('255.75')],dtype=P.LD)
        dt=P.LD(1)/128;plan=P.phase_plan(k,dt,order)
        u=np.asarray([P.CD('0.02')+P.CD('0.013')*1j]*3,dtype=P.CD)
        w=np.asarray([P.CD('0.07')-P.CD('0.011')*1j]*3,dtype=P.CD)
        z=plan['offset'];g=P.LD('0.43')+P.LD('-0.37')*z+P.LD('0.17')*z*z
        up,wp=P.phase_step(u,w,plan,g)
        flipped=P.phase_step(u,w,plan,-g)
        for i,kn in enumerate(k):
            km=P.represented(ctx,kn);tm=P.represented(ctx,dt);eps=ctx.mpf(P.EPS_RATIO[0])/P.EPS_RATIO[1]
            um=P.complex_represented(ctx,u[i]);wm=P.complex_represented(ctx,w[i]);lam=2*ctx.j*km
            a,b,c=[P.represented(ctx,P.LD(x)) for x in ('0.43','-0.37','0.17')]
            forcing=lambda s:a+b*s+c*s*s
            Ew=ctx.quad(lambda s:ctx.exp(lam*(tm-s))*forcing(s),[0,tm])
            Eu=ctx.quad(lambda s:ctx.expm1(lam*(tm-s))/lam*forcing(s),[0,tm])
            truthw=ctx.exp(lam*tm)*wm-eps*Ew;truthu=um+ctx.expm1(lam*tm)/lam*wm-eps*Eu
            error=max(abs(P.complex_represented(ctx,wp[i])-truthw),abs(P.complex_represented(ctx,up[i])-truthu))
            P.require(error<ctx.mpf('1e-17'),'Fabricated high-phase flow truth mismatch')
            wrong=max(abs(P.complex_represented(ctx,flipped[1][i])-truthw),abs(P.complex_represented(ctx,flipped[0][i])-truthu))
            P.require(wrong>ctx.mpf('1e-8'),'Changed forcing sign did not reject')
            checks.append({'dps':dps,'order':order,'k':str(kn),'truth_error':ctx.nstr(error,dps),
                           'changed_sign_error':ctx.nstr(wrong,dps),'nonzero_unprojected_cR':True})
# The study source is guarded even when the module has been imported.
try:P.source_jet(np.asarray([P.LD('-4.1')]),'positive_B')
except ValueError:checks.append({'physical_source_guard':'PASS'})
else:raise ValueError('Physical source bypass accepted')
result={'status':'PASS','check_count':len(checks),'physical_source_or_arrays_used':False,
        'python_optimization':sys.flags.optimize,'producer_sha256':P.sha(path),'checks':checks}
output=Path(__file__).with_name('SYNTHETIC_TRUTH_'+('OPTIMIZED' if sys.flags.optimize else 'NORMAL')+'_'+P.sha(path)[:12]+'.json')
P.require(not output.exists(),'Synthetic truth output already exists')
output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
