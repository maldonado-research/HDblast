#!/usr/bin/env python3
"""Post-processing: EFT-validity diagnostics in 5D units, case-A boundary approach, merged results JSON."""
import numpy as np, json
from blast_common import *
from blast_homogeneous import evolve, first
R = dict(homogeneous=json.load(open("BLAST_HOMOGENEOUS_RESULTS.json")), particles=json.load(open("BLAST_PARTICLE_RESULTS.json")),
         reheating=json.load(open("BLAST_REHEATING_RESULTS.json")), table_checks=json.load(open("BLAST_TABLE_CHECKS.json")))
for key in list(R['particles'].keys()):
    if isinstance(R['particles'][key], dict):
        for kk in ('N_k', 'kappa'): R['particles'][key].pop(kk, None)     # full spectra stay in BLAST_PARTICLE_RESULTS.json
X = {}
d = np.load("traj_throat_k0.npz"); f0 = 2*Ip; v0 = (1)/f0**2
for tt in (1e-3, 1e-6):
    st = np.sqrt(tt); pj = np.abs(d['phidotJ_over_sqrt_t'])*st; HJ = np.abs(d['hJ'])*np.sqrt(f0*tt*v0/3); sp = np.abs(d['sinhpsi_over_sqrt_t'])*st
    e = {}
    for nm, arr, thr in (("|dphi_b/dtau_J|>0.1 (5D units)", pj, 0.1), ("|dphi_b/dtau_J|>0.3", pj, 0.3), ("|H_J|>0.1", HJ, 0.1), ("naive |sinh psi|>1 (leaf rapidity)", sp, 1.0)):
        i = first(arr > thr); e[nm] = None if i is None else dict(tau=float(d['tau'][i]), phi_b=float(d['phi'][i]), V_E_over_V0=float(d['u'][i]))
    ic = first(d['Th'] >= 0.5479058809936374)
    e['at phi_b=-1 crossing'] = dict(dphi_dtauJ=float(pj[ic]), H_J=float(HJ[ic]), H_J_times_l_minus=float(HJ[ic]*9/5), one_plus_phi_where_naive_sinhpsi_is_1=float(pj[ic]/2))
    X["throat side, t=%g" % tt] = e
R['EFT_validity_throat_side'] = X
# case A boundary approach
T, C, Th_end = build_tables(); b = 1 + c_star
LA = Landscape(lambda p: (1-p)*(1+b*p), lambda p: -(1+b*p) + b*(1-p), T, C, Th_end)
sel = (LA.phi > 1-1e-7) & (LA.phi < 1-1e-30) & (LA.v > 0)
dTh = LA.Th[sel]-LA.Th_plus_end; R['reheating']['A']['dlnV_E_dln_DeltaTheta_asymptotic_fit'] = float(np.polyfit(np.log(dTh[::50]), np.log(LA.v[sel][::50]), 1)[0])
mu2 = R['reheating']['A']['mu2_hilltop']; sA = (-3+np.sqrt(9-4*mu2))/2; dq = np.sqrt(1e-3*LA.v0/3)/(2*np.pi)
o = evolve(LA, -dq, -sA*dq, 0, 60.0, dtau=5e-4); sp = np.abs(o['sinhpsi_over_sqrt_t'])
A = {}
for tt in (1e-3, 1e-6):
    for thr in (0.3, 1.0):
        i = first(sp*np.sqrt(tt) > thr); A["t=%g |sinh psi|>%.1f" % (tt, thr)] = dict(tau=float(o['tau'][i]), phi_b_minus_1=float(o['phi'][i]-1), y_b=float(o['y'][i]), DeltaTheta=float(o['Th'][i]-LA.Th_plus_end), w=float(((o['K']-o['u'])/(o['K']+o['u']))[i]))
R['reheating']['A']['boundary_approach'] = A
print(json.dumps(X, indent=1)); print(json.dumps(R['reheating']['A'], indent=1))
json.dump(R, open("BLAST_RESULTS.json", "w"), indent=1)
