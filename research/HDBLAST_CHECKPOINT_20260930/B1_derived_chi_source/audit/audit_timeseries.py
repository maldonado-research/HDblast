#!/usr/bin/env python3
"""Independent recomputation (audit) of descriptive B1 numbers from the saved timeseries.
Does not import the producer's analysis code. Output: audit/AUDIT_TIMESERIES.json"""
import json, glob, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)

def load(npz):
    z = np.load(npz); cols = [str(c) for c in z['cols']]; rec = z['rec']
    return {c: rec[:, i] for i, c in enumerate(cols)}

def rel_end(t, thr=0.05):
    bad = np.flatnonzero((np.nan_to_num(t['Hmax_near']) > thr) | (np.nan_to_num(t['Mmax_near']) > thr))
    return int(bad[0]) - 1 if len(bad) else len(t['T']) - 1

def longest_plateau(t, iend, tol=0.05):
    """brute force: max Delta ln a of any window ending at p with Wa4,Ra4 spread <= tol*|value at p| and H>0"""
    lna = t['ln_a'][:iend+1]; H = t['H_over_H0'][:iend+1]
    Wa4 = (t['Wy']*np.exp(4*t['ln_a']))[:iend+1]; Ra4 = (t['R']*np.exp(4*t['ln_a']))[:iend+1]
    best = 0.0
    for p in range(len(lna)):
        for j in range(p, -1, -1):
            if H[j] <= 0: break
            w = Wa4[j:p+1]; rr = Ra4[j:p+1]
            if (w.max()-w.min()) > tol*abs(Wa4[p]) or Ra4[p] <= 0 or (rr.max()-rr.min()) > tol*Ra4[p]: break
            best = max(best, lna[p]-lna[j])
    return best

out = {}
for npz in sorted(glob.glob(os.path.join(ROOT, 'runs', '*', '*_timeseries.npz'))):
    tag = os.path.basename(npz).replace('_timeseries.npz', '')
    if tag.startswith('pre'): continue
    t = load(npz); i = rel_end(t); n = len(t['T']) - 1
    H = t['H_over_H0']
    def vals(k):
        rad = t['rad'][k]; return dict(H0tau=float(t['H0tau'][k]), r=float(t['Wy'][k]/rad) if rad > 0 else None,
            Omega_r=float(rad/H[k]**2), Omega_chi=float(t['chi'][k]/H[k]**2), Omega_vac=float(t['vac'][k]/H[k]**2),
            Wa4=float(t['Wy'][k]*np.exp(4*t['ln_a'][k])), Ra4=float(t['R'][k]*np.exp(4*t['ln_a'][k])), ln_a=float(t['ln_a'][k]))
    # post-decay window: chi energy < 1e-3 of R-energy proxy (rad)
    chi = t['chi'][:i+1]; rad = t['rad'][:i+1]
    Rr = t['R'][:i+1]; started = np.cumsum(Rr > 0) > 0
    dec = np.flatnonzero(started & (Rr > 0) & (chi < 1e-2*np.maximum(rad, 1e-300)))
    post = None
    if len(dec):
        k0 = dec[0]; Wa4 = (t['Wy']*np.exp(4*t['ln_a']))[k0:i+1]; Ra4 = (t['R']*np.exp(4*t['ln_a']))[k0:i+1]
        post = dict(H0tau_start=float(t['H0tau'][k0]), dlna=float(t['ln_a'][i]-t['ln_a'][k0]),
                    Wa4_spread=float((Wa4.max()-Wa4.min())/abs(Wa4[-1])), Ra4_spread=float((Ra4.max()-Ra4.min())/Ra4[-1]))
    closure = t['vac']+t['rad']+t['chi']+t['kin']+t['fric']+t['Wy']-H**2
    out[tag] = dict(n_rec=n+1, reliable_end_idx=i, at_reliable_end=vals(i), at_run_end=vals(n),
                    max_near_constraint_reliable=float(np.nanmax(np.maximum(t['Hmax_near'][:i+1], t['Mmax_near'][:i+1]))),
                    longest_plateau_dlna=float(longest_plateau(t, i)), post_decay=post,
                    max_Rj=float(np.nanmax(np.abs(t['J'][:i+1])/np.abs(t['sigma1'][:i+1]))),
                    closure_max_abs=float(np.nanmax(np.abs(closure[:i+1]))))
json.dump(out, open(os.path.join(HERE, 'AUDIT_TIMESERIES.json'), 'w'), indent=1)
for k, v in out.items():
    a = v['at_reliable_end']; b = v['at_run_end']; p = v['post_decay']
    print(f"{k:48s} rel_end={a['H0tau']:.2f} r={a['r'] if a['r'] is None else round(a['r'],4)} Om_r={a['Omega_r']:.3g} | run_end={b['H0tau']:.2f} r={b['r'] if b['r'] is None else round(b['r'],4)} Om_r={b['Omega_r']:.3g} | plat={v['longest_plateau_dlna']:.3f} con={v['max_near_constraint_reliable']:.2g} Rj={v['max_Rj']:.3g} post={None if p is None else (round(p['dlna'],3), round(p['Wa4_spread'],4), round(p['Ra4_spread'],4))}")
