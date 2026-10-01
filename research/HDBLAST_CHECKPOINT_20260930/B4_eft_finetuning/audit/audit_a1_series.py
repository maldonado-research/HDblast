#!/usr/bin/env python3
"""Audit: independent reads of the A1 5D time series (read-only) to check B4 claims about A1:
late vac value, phi_s, Weyl a^4 drop between phi_b=0.95 and plateau, e-folds with |Omega_vac|<=0.1,
recollapse-time closed form, and the recipe-on-5D r.  Output: audit_a1_series.json"""
import json, math, sys
from pathlib import Path
import numpy as np
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum')
sys.path.insert(0, str(A1))
HERE = Path(__file__).resolve().parent
out = {}
for tag in ('main_dstar_Y1_dc1e-2_dzf5e-4', 'main_dstar_Y1_dc1e-2_dzf1e-3', 'main_dstar_Y0.3_dc1e-2_dzf5e-4'):
    d = np.load(A1/'runs'/'main'/(tag+'_timeseries.npz')); cols = [str(c) for c in d['cols']]; R = d['rec']
    g = {c: R[:, i] for i, c in enumerate(cols)}
    H = g['H_over_H0']; lna = g['ln_a']; tau = g['H0tau']
    Wa4 = g['Wy']*np.exp(4*lna); Ra4 = g['R']*np.exp(4*lna)
    Omv = g['vac']/H**2; Omr = g['rad']/H**2
    i95 = int(np.flatnonzero(g['phi_b'] >= 0.95)[0])
    # A1 plateau via A1's own function
    import a1_analyze as A
    dd = dict(ln_a=lna, H_over_H0=H, Wa4=Wa4, Ra4=Ra4)
    p, j0 = A.plateau_index(dd)
    # e-folds with |Omega_vac|<=0.1, H>0 (whole run, longest contiguous)
    m = (np.abs(Omv) <= 0.1) & (H > 0) & (g['phi_b'] > 1.0) & (np.abs(g['v_over_H0']) < 0.01)
    best = 0.0; st = None
    for i in range(len(m)):
        if m[i] and st is None: st = i
        if (not m[i] or i == len(m)-1) and st is not None:
            en = i if m[i] else i-1; best = max(best, float(lna[en]-lna[st])); st = None
    # also counting from first record with v small (scalar stopped)
    istop = int(np.flatnonzero((np.abs(g['v_over_H0']) < 0.01) & (g['phi_b'] > 1.0))[0])
    out[tag] = dict(n=len(tau), tau_end=float(tau[-1]), vac_end=float(g['vac'][-1]), phi_end=float(g['phi_b'][-1]),
                    vac_at_stop=float(g['vac'][istop]), tau_stop=float(tau[istop]), H_stop=float(H[istop]),
                    Omega_vac_stop_actual=float(Omv[istop]),
                    Omega_vac_at_stop_using_lam=float(-3.3992e-4/H[istop]**2),
                    tau_phi095=float(tau[i95]), Wa4_phi095=float(Wa4[i95]),
                    plateau=None if p is None else dict(tau=float(tau[p]), r=float(g['Wy'][p]/g['rad'][p]), Omega_vac=float(Omv[p]),
                                                        Omega_r=float(Omr[p]), Wa4=float(Wa4[p]), lna_since_stop=float(lna[p]-lna[istop])),
                    Wa4_ratio_phi095_over_plateau=None if p is None else float(Wa4[i95]/Wa4[p]),
                    efolds_absOmega_vac_le_0p1_after_stop=best,
                    vac_series_tail=[float(x) for x in g['vac'][-5:]], H_end=float(H[-1]))
json.dump(out, open(HERE/'audit_a1_series.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
