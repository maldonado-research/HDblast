#!/usr/bin/env python3
"""Read-only causal sampling diagnostics and exact paired-profile identities.

No PDE evolution or shooting solve. Requires NumPy and SymPy.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import sympy as S

base=Path(__file__).resolve().parent/'inputs'
rp,vp,sp,U0,dU,dU1,U10,E,k,II,eps,b,dv=S.symbols('rho0 v0 s0 U0 dU dU1 U10 exp_a delta_s I epsilon bump delta_v',real=True)
# E represents exp(a)>0. The reference has I0=0 and the positive radial branch.
v02=1+rp**2*(sp**2/12-U0/6)
full2=1+rp**2*E**2*((sp+k)**2/12-(U0+dU)/6)-II/(rp**2*E**2)
DX=rp**2*((E**2-1)*(sp**2/12-U0/6)+E**2*((2*sp*k+k*k)/12-dU/6))-II/(rp**2*E**2)
checks=[]
def check(name,expr,zero=True):
    out=S.factor(S.expand(expr));ok=(out==0) if zero else (out!=0)
    checks.append(dict(name=name,pass_check=bool(ok),negative_control=not zero))
    if not ok:raise RuntimeError(name+': '+str(out))
check('difference of squared radial slopes',full2-v02-DX)
check('rationalized slope difference',(vp+dv)**2-vp**2-dv*(2*vp+dv))
check('scalar difference equation',rp*E*(sp+k)-rp*sp-rp*((E-1)*sp+E*k))
full_s_z=rp*E*(U10+dU1)-4*(vp+dv)*(sp+k)+eps*rp*E*b
ref_s_z=rp*U10-4*vp*sp
stable_s_z=rp*((E-1)*U10+E*dU1)-4*(vp*k+dv*(sp+k))+eps*rp*E*b
check('proper scalar slope difference equation',full_s_z-ref_s_z-stable_s_z)
check('proper coordinate difference equation',rp*E-rp-rp*(E-1))
check('mass integral conformal derivative',eps*(rp*E)**4*(sp+k)*b/6*(rp*E)-eps*rp**5*E**5*(sp+k)*b/6)
check('reject missing conformal rho in source',full_s_z-ref_s_z-(stable_s_z-eps*rp*E*b+eps*b),False)
check('reject opposite mass-integral sign',full2-v02-(DX+2*II/(rp**2*E**2)),False)
check('reject omitting radial slope cross term',full_s_z-ref_s_z-(stable_s_z+4*dv*k),False)

profile_file=base/'BALANCED_CONSTRAINT_SEED_PROFILES.npz'
meta_file=base/'BALANCED_CONSTRAINT_SEED_RESULTS.json'
meta=json.loads(meta_file.read_text())
rows=[]
with np.load(profile_file,allow_pickle=False) as f:
    ref=f['0.0'];r_yb=float(ref[0,-1]);r_phib=float(ref[3,-1]-1)
    for name in f.files:
        y,r,v,p,s,z,I=f[name]
        monotone=bool(np.all(np.diff(z)>0) and np.all(np.diff(p)>0) and np.all(r>0))
        if not monotone:raise RuntimeError('Unexpected profile monotonicity')
        edges={str(a):float(np.interp(a,p,z)) for a in [.2,.5,.55,.85]}
        rows.append(dict(epsilon=float(name),monotone_positive_profile=monotone,
            conformal_support_edges=edges,earliest_compact_influence_time=-edges['0.85'],
            conformal_table_start=float(z[0]),brane_y=float(y[-1]),brane_phi=float(p[-1]-1),
            brane_y_difference=float(y[-1]-r_yb),brane_phi_difference=float(p[-1]-1-r_phib),
            leading_coordinate_shift_contribution=float(s[-1]*(y[-1]-r_yb)),
            interpolation_scope='Linear interpolation of archived proper-distance samples for causal-scale diagnostics only; do not use this interpolation to initialize the PDE.'))
out={'status':'PASS','positive_symbolic_checks':sum(not c['negative_control'] for c in checks),
     'negative_controls':sum(c['negative_control'] for c in checks),'symbolic_checks':checks,
     'profiles':rows,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [profile_file,meta_file]},
     'scope':'Exact candidate difference equations and archived-profile causal diagnostics; not an implemented alternative integrator or an evolved solution.'}
Path(__file__).with_name('SAMPLING_DIFFERENCE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','positive_symbolic_checks','negative_controls','profiles']},indent=2))
