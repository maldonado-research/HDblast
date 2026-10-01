#!/usr/bin/env python3
"""Audit: compare audit reruns with the producer's saved grid run. Output audit/AUDIT_RERUNS.json"""
import json, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
def load(p):
    z = np.load(p); return {str(c): z['rec'][:, i] for i, c in enumerate(z['cols'])}
orig = load(os.path.join(ROOT, 'runs/main/main_ps0.5_G100_y1_bcut_dzf1e-3_timeseries.npz'))
fine = load(os.path.join(ROOT, 'runs/main/main_ps0.5_G100_y1_bcut_dzf5e-4_timeseries.npz'))
out = {}
for tag in ['rerun_ps0.5_G100_y1_bcut_dzf1e-3', 'rerun_ps0.5_G100_y1_bcut_dzf7.5e-4', 'rerun_ps0.5_G100_y1_bcut_wrongsign_dzf1e-3']:
    x = load(os.path.join(HERE, 'runs', tag + '_timeseries.npz'))
    def at_end(t):
        k = len(t['T'])-1; H = t['H_over_H0'][k]
        return dict(H0tau=float(t['H0tau'][k]), r=float(t['Wy'][k]/t['rad'][k]), Omega_r=float(t['rad'][k]/H**2),
                    Wa4=float(t['Wy'][k]*np.exp(4*t['ln_a'][k])), Ra4=float(t['R'][k]*np.exp(4*t['ln_a'][k])))
    o = dict(end=at_end(x), orig_dzf1e3_end=at_end(orig), orig_dzf5e4_end=at_end(fine))
    n = min(len(x['T']), len(orig['T']))
    if 'dzf1e-3' in tag:
        o['max_abs_diff_vs_orig'] = {k: float(np.nanmax(np.abs(x[k][:n]-orig[k][:n]))) for k in ['H0tau', 'phi_b', 'H_over_H0', 'Wy', 'R']}
    # ledger residual: d(X+R)/dtau + 3H(X+R+P) - J v  (all in run units; use T derivative)
    out[tag] = o
json.dump(out, open(os.path.join(HERE, 'AUDIT_RERUNS.json'), 'w'), indent=1); print(json.dumps(out, indent=1))
