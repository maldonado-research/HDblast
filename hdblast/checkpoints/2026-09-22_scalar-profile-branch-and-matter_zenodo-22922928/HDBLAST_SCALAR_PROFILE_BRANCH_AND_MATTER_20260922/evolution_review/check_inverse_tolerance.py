#!/usr/bin/env python3
"""Initial-data-only inverse-map tolerance experiment via function override.

No producer source is modified and no PDE time step is taken. All six setups
share the same two dense radial solutions; only brentq's inverse tolerance varies.
"""
from pathlib import Path
import hashlib,importlib.util,json,time
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'evolution/evolve_balanced.py'
spec=importlib.util.spec_from_file_location('inverse_tolerance_audit_wrapper',path)
E=importlib.util.module_from_spec(spec);spec.loader.exec_module(E)
cache={};inverse_checks=[];xtol_current=2e-14;rtol_inverse=4*np.finfo(float).eps
def override_profile(epsilon,rtol=8e-14):
    key=(epsilon,rtol)
    if key not in cache:cache[key]=E.B.run(epsilon,rtol=rtol,transport=True)
    sol,meta=cache[key];zb=sol.y[4,-1]
    def sample(z):
        targets=np.asarray(z)+zb
        if np.min(targets)<sol.y[4,0]:raise ValueError('Grid exceeds dense table')
        ix=np.clip(np.searchsorted(sol.y[4],targets)-1,0,len(sol.t)-2)
        ys=np.array([brentq(lambda y:sol.sol(y)[4]-target,sol.t[i],sol.t[i+1],xtol=xtol_current,rtol=rtol_inverse)
                     for target,i in zip(targets,ix)])
        vals=sol.sol(ys)
        inverse_checks.append(dict(epsilon=epsilon,xtol=xtol_current,max_z_equation_residual=float(np.max(abs(vals[4]-targets)))))
        return vals
    return sample,meta
E.profile=override_profile
rows=[];states={};start=time.monotonic()
for h in [2e-4,1e-4,5e-5]:
    for xtol in [2e-14,1e-300]:
        xtol_current=xtol
        s,v,meta=E.setup(.01,h,3.,stretch=2.)
        diag,H,M,C=E.diagnostics(s,v,0.)
        mask=(s.z>-.8*s.L)&(np.arange(s.n)>5);ii=np.flatnonzero(mask)
        scale=1+6*s.hc*s.hc+s.phz*s.phz;j=ii[np.argmax(abs(H[mask])/scale[mask])]
        row=dict(hmin=h,xtol=xtol,rtol_inverse=rtol_inverse,nodes=s.n,
                 global_H_background_norm=diag['hamiltonian_max'],core_H_background_norm=diag['core_H'],
                 weighted_C_max=diag['weighted_characteristic_max'],peak_z=float(s.z[j]),
                 peak_raw_H=float(H[j]),core_raw_H=diag['core_raw_H'],
                 analytic_boundary_slope_minus_target=meta['analytic_boundary_slope_minus_target'],
                 seed_mass_defect=meta['seed']['mass_defect_b'],
                 initial_scalar_acceleration_at_shell=float(s.rhs(v)[5,-1]),
                 max_inverse_z_residual=max(r['max_z_equation_residual'] for r in inverse_checks[-2:]))
        if xtol==2e-14:states[h]=(v.copy(),s.rho.copy(),s.phi.copy())
        else:
            old,oldrho,oldphi=states[h]
            row['difference_from_default_max_by_component']=np.max(abs(v-old),axis=1).tolist()
            row['reference_rho_change_max']=float(np.max(abs(s.rho-oldrho)))
            row['reference_phi_change_max']=float(np.max(abs(s.phi-oldphi)))
        rows.append(row);print(json.dumps(row),flush=True)
out=dict(status='COMPLETED_INITIALIZATION_ONLY_TOLERANCE_EXPERIMENT',epsilon=.01,L=3.,stretch=2.,
         dense_ode_rtol=8e-14,distinct_dense_solutions=len(cache),rows=rows,
         runtime_seconds=time.monotonic()-start,wrapper_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
         scope='Function override in this process only; no source edits, no PDE evolution; tolerance change cannot remove all dense-profile differentiation error.')
Path(__file__).with_name('INVERSE_TOLERANCE_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
