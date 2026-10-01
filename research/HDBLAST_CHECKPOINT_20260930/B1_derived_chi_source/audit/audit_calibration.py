#!/usr/bin/env python3
"""Audit: compare K1/K2 timeseries with the A1 archives directly (no producer analysis code). Output audit/AUDIT_CALIBRATION.json"""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
A1 = '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main'
def load(p):
    z = np.load(p); return {str(c): z['rec'][:, i] for i, c in enumerate(z['cols'])}
out = {}
for b1, a1 in [('K1_chioff_Y0_dc1e-2_dzf1e-3', 'main_dstar_Y0_dc1e-2_dzf1e-3'), ('K2_toy_Y1_dc1e-2_dzf1e-3', 'main_dstar_Y1_dc1e-2_dzf1e-3')]:
    x = load(os.path.join(ROOT, 'runs/cal', b1 + '_timeseries.npz')); y = load(os.path.join(A1, a1 + '_timeseries.npz'))
    n = min(len(x['T']), len(y['T']))
    d = {k: float(np.nanmax(np.abs(x[k][:n] - y[k][:n]))) for k in ['T', 'phi_b', 'H_over_H0', 'Wy', 'ln_a'] if k in x and k in y}
    rx = x['Wy'][:n]/np.where(x['rad'][:n] > 0, x['rad'][:n], np.nan); ry = y['Wy'][:n]/np.where(y['rad'][:n] > 0, y['rad'][:n], np.nan)
    d['r_maxreldiff'] = float(np.nanmax(np.abs(rx-ry)/np.abs(ry))) if np.any(np.isfinite(ry)) else None
    out[b1] = dict(n_B1=len(x['T']), n_A1=len(y['T']), maxabs=d, end_tau=(float(x['H0tau'][-1]), float(y['H0tau'][-1])))
a1sum = json.load(open(os.path.join(A1, 'main_dstar_Y1_dc1e-2_dzf1e-3_summary.json')))
json.dump(out, open(os.path.join(HERE, 'AUDIT_CALIBRATION.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
