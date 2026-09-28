#!/usr/bin/env python3
"""Independent initialization checks; no time evolution. NumPy/SciPy required."""
from pathlib import Path
import hashlib,importlib.util,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'evolution/evolve_balanced.py'
spec=importlib.util.spec_from_file_location('audited_wrapper',path)
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)
rows=[]
for eps,h in [(0.,2e-4),(.01,2e-4),(.01,1e-4)]:
    s,state,meta=E.setup(eps,h,6.)
    # Recompute dense seed independently of the sampled deviation field.
    sample,seedmeta=E.profile(eps)
    r,v,p,sy,zz,I=sample(s.z)
    D=I/r**2
    k=seedmeta['compensation']
    bump=np.array([E.B.local_bump(x,.35,.15)-k*E.B.local_bump(x,.70,.15) for x in p])
    exact_acc=np.array([-2*D,4*D,eps*r*r*bump])
    acc=s.rhs(state)[3:]
    # No taper is active in the local comparison domain.
    mask=s.z>-.15
    err=np.abs(acc-exact_acc)
    near=p>.85
    diag,HH,MM,CC=E.diagnostics(s,state,0.)
    row=dict(epsilon=eps,hmin=h,nodes=s.n,
        shell_scalar_displacement=meta['initial_shell_scalar_displacement'],
        max_acceleration_error_near_wall=np.max(err[:,mask],axis=1).tolist(),
        max_exact_acceleration_near_wall=np.max(np.abs(exact_acc[:,mask]),axis=1).tolist(),
        shell_acceleration=acc[:,-1].tolist(),analytic_shell_acceleration=exact_acc[:,-1].tolist(),
        initial_normalized_H=diag['core_H'],initial_normalized_M=diag['core_M'],
        analytic_boundary_slope_minus_target=meta['analytic_boundary_slope_minus_target'],
        phi_gap_to_support=float(p[-1]-.85),
        nearest_support_edge_time=float(-np.interp(.85,p,s.z)),
        reference_transport_formulation=bool(meta['reference']['transport_formulation']),
        seed_transport_formulation=bool(seedmeta['transport_formulation']))
    if eps==0:
        row['zero_seed_exact_array']=bool(np.count_nonzero(state)==0)
        row['zero_seed_rhs_max']=float(np.max(np.abs(s.rhs(state))))
        if not row['zero_seed_exact_array'] or row['zero_seed_rhs_max']!=0:raise RuntimeError('zero-seed control failed')
    if max(abs(x) for x in row['analytic_boundary_slope_minus_target'])>1e-9:raise RuntimeError('Boundary slope failure')
    rows.append(row);print(json.dumps(row),flush=True)
coarse,fine=rows[1:]
ratio=[c/f if f else None for c,f in zip(coarse['max_acceleration_error_near_wall'],fine['max_acceleration_error_near_wall'])]
if any(r is None or r<4 for r in ratio):raise RuntimeError('Initial acceleration errors did not converge adequately')
result=dict(status='PASS_INITIALIZATION_CONTROLS_NO_EVOLUTION',rows=rows,
            acceleration_error_reduction_factor=ratio,
            wrapper_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            scope='Floating independent dense-profile/initial-acceleration controls only; does not prove time-evolution stability.')
Path(__file__).with_name('WRAPPER_INITIALIZATION_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
