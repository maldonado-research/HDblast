#!/usr/bin/env python3
"""Matched-time state comparison; replays only ten already specified RK4 steps.

The coarse saved state is at .4992. Its ten steps of dt=.00008 reconstruct the
unrecorded .5 state using the original RHS. Fine states are read at .5 directly.
No background shooting, new seed or full evolution is executed.
"""
from pathlib import Path
import hashlib,importlib.util,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('frozen_comparison_solver',ROOT/'frozen/registered_solver.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)

def restore(file):
    meta=json.loads(file.with_suffix('.json').read_text())
    with np.load(file,allow_pickle=False) as f:
        z=f['z'].copy();bg=f['background'].copy()
        target_step=round(.5/meta['dt']);name='state_'+str(target_step)
        if name in f:state=f[name].copy();start=.5;steps=0
        elif abs(meta['records'][-1]['time']-.5)<1e-12:state=f['final'].copy();start=.5;steps=0
        else:state=f['state_6240'].copy();start=6240*meta['dt'];steps=round((.5-start)/meta['dt'])
    if steps not in [0,10]:raise RuntimeError('Unapproved step count '+str(steps))
    s=R.Solver.__new__(R.Solver);s.z=z;s.n=len(z);s.hmin=meta['hmin'];s.L=meta['L'];s.stretch=meta.get('stretch',.05)
    s.zp=s.hmin*(1-z/s.stretch);s.zpp=-s.hmin*s.zp/s.stretch
    s.rho,s.phi,s.hc,s.phz=bg;s.rho2=s.rho**2;s.rb=s.rho[-1];s.pb=s.phi[-1]
    s.tdet=.001;s.ko=0.;s.pot=R.derivatives(s.phi);s.s0=2*R.W(s.pb)+s.tdet*(1+R.C*s.pb)
    s.s10=2*(s.pb*s.pb-1)+s.tdet*R.C;s.s20=4*s.pb
    # Reconstruct the original analytic map identically, rather than its
    # algebraically equivalent zp=h*(1-z/stretch), for the bounded replay.
    zr,zpr,zppr=R.grid(s.hmin,s.L,s.stretch)
    if np.max(abs(zr-z))>1e-13:raise RuntimeError('Stored map does not match source')
    s.zp=zpr;s.zpp=zppr
    if steps:
        s.D1,s.D2,s.g1,s.g2,s.ghost=R.matrices(s.zp,s.zpp)
        dt=meta['dt']
        for i in range(steps):
            k1=s.rhs(state);k2=s.rhs(state+.5*dt*k1);k3=s.rhs(state+.5*dt*k2);k4=s.rhs(state+dt*k3)
            state+=dt*(k1+2*k2+2*k3+k4)/6
        diag,H,M,C=s.diagnostics(state,.5)
    else:
        record=next(r for r in meta['records'] if abs(r['time']-.5)<1e-12)
        diag={k:record[k] for k in ['momentum_max','hamiltonian_max','boundary_velocity_slope_error','weighted_characteristic_max']}
    total=state.copy();total[0]+=np.log(s.rho);total[1]+=np.log(s.rho);total[2]+=s.phi;total[3]+=1
    info=dict(source=file.name,source_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),hmin=s.hmin,nodes=s.n,
              start_time=start,additional_RK4_steps=steps,dt=meta['dt'],comparison_time=.5,diagnostics=diag)
    return z,bg,state,total,info

files=[p for p in (ROOT/'runs').glob('balanced_epsp01_wide_*.npz') if p.with_suffix('.json').exists()]
files.sort(key=lambda p:json.loads(p.with_suffix('.json').read_text())['hmin'],reverse=True)
states=[restore(p) for p in files]
if len(states)<2:raise RuntimeError('Need at least two complete wide runs')
zc=states[0][0];aligned=[];metadata=[]
for z,bg,state,total,info in states:
    ratio=round((len(z)-1)/(len(zc)-1));inds=np.arange(len(zc))*ratio
    if len(z)==len(zc):inds=np.arange(len(zc))
    if np.max(abs(z[inds]-zc))>1e-13:raise RuntimeError('Grids are not nested')
    aligned.append(total[:,inds]);info['coarse_node_stride']=ratio;metadata.append(info)
if metadata[0]['additional_RK4_steps']:
    np.savez_compressed(Path(__file__).with_name('MATCHED_TIME_COARSE_STATE.npz'),z=states[0][0],background=states[0][1],state=states[0][2],time=.5)
pairwise=[]
for i in range(len(states)-1):
    delta=aligned[i]-aligned[i+1]
    pairwise.append(dict(coarser_h=metadata[i]['hmin'],finer_h=metadata[i+1]['hmin'],
        Linf_by_component=np.max(abs(delta),axis=1).tolist(),
        L2_dz_by_component=np.sqrt(np.trapezoid(delta*delta,zc,axis=1)).tolist()))
ratio=None
if len(pairwise)==2:
    ratio=(np.array(pairwise[0]['L2_dz_by_component'])/np.array(pairwise[1]['L2_dz_by_component'])).tolist()
out=dict(status='MATCHED_TIME_FLOATING_REFINEMENT_DIAGNOSTIC',time=.5,runs=metadata,pairwise=pairwise,
         L2_reduction_factor_three_grids=ratio,
         components=['A-t','B','phi','A_t','B_t','phi_t'],
         norm_scope='Per-component coordinate L2(dz) and Linf on identical coarse nodes; not a physical energy norm or error bound.',
         replay_scope='Only the coarse ten-step RK4 replay is newly integrated; all finer states are directly archived at .5.',
         frozen_solver_sha256=hashlib.sha256((ROOT/'frozen/registered_solver.py').read_bytes()).hexdigest())
Path(__file__).with_name('MATCHED_TIME_REFINEMENT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
