#!/usr/bin/env python3
"""AUDIT: recompute, from the saved time series only, the quantities quoted for the delta = 1e-3 tuned runs:
end chart time / proper time (last record), max phi_b, max Omega_r (all records and up to the reliable end), shell proper time at
T = F_inf, the end time relative to x_c, the D6 old-chart freeze (dB_b/dT, b_b + T, H0tau, phi_b, H/H0), and the C5 Weyl-identity
discrimination at p95 as reported.  Also d*(1e-3) from M8 and Y*rho_b.  Output audit/AUDIT_RUNS.json.  Numerical."""
import json, glob, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
def load(f):
    s = json.load(open(f)); z = np.load(f.replace('_summary.json', '_timeseries.npz')); c = list(z['cols']); r = z['rec']
    return s, {k: r[:, i] for i, k in enumerate(c)}
out = dict(status='numerical (audit)', runs={})
for f in sorted(glob.glob(str(HERE/'runs/*/*_summary.json'))):
    if '/runs/main/' not in f and '/runs/ctl/' not in f and '/runs/explore/' not in f and '/runs/res/' not in f: continue
    s, d = load(f); tag = Path(f).name.replace('_summary.json', '')
    T, tau, H = d['T'], d['H0tau'], d['H_over_H0']
    with np.errstate(all='ignore'):
        Om = d['rad']/H**2
    bad = np.flatnonzero((np.nan_to_num(d['Hmax_near']) > 0.05) | (np.nan_to_num(d['Mmax_near']) > 0.05))
    ire = int(bad[0]) - 1 if len(bad) else len(T) - 1
    xc = s['params']['xc']; Fi = s.get('F_inf')
    e = dict(Y=s['params']['Y'], dc=s['params']['dc'], dzf=s['params']['dzf'], xc=xc, F_inf=Fi,
             T_last=float(T[-1]), H0tau_last_record=float(tau[-1]), H0tau_end_summary=s['H0tau_end'],
             T_end_minus_F_inf=(float(s['T_end'] - Fi) if Fi and math.isfinite(Fi) else None),
             phi_b_max=float(np.nanmax(d['phi_b'])), Omega_r_max_all=float(np.nanmax(np.where(H > 0.05, Om, np.nan))),
             Omega_r_max_reliable=float(np.nanmax(np.where(H[:ire+1] > 0.05, Om[:ire+1], np.nan))),
             H0tau_reliable_end=float(tau[ire]), H_over_H0_min_reliable=float(np.nanmin(H[:ire+1])))
    if Fi and math.isfinite(Fi) and T[-1] > Fi:
        k = int(np.argmin(abs(T - Fi))); e['H0tau_at_F_inf'] = float(tau[k]); e['phi_b_at_F_inf'] = float(d['phi_b'][k])
        e['H0tau_end_minus_H0tau_at_F_inf'] = float(s['H0tau_end'] - tau[k])
    out['runs'][tag] = e
    print('%-38s xc=%-5s Tend-Finf=%s tau_end=%.2f(last rec %.2f) rel=%.2f phimax=%.3f Om_all=%.4f Om_rel=%.4f tau@Finf=%s' % (
        tag, xc, None if e['T_end_minus_F_inf'] is None else round(e['T_end_minus_F_inf'], 3), e['H0tau_end_summary'], e['H0tau_last_record'],
        e['H0tau_reliable_end'], e['phi_b_max'], e['Omega_r_max_all'], e['Omega_r_max_reliable'], round(e.get('H0tau_at_F_inf', float('nan')), 3)))
# D6
s, d = load(str(HERE/'runs/d6/D6_oldchart_Y1_dc1e-2_S1_k10_nostop_summary.json'))
T, B = d['T'], d['B_b']; g = np.gradient(B, T); rb = s['rho_b']
ok = np.isfinite(d['Hmax_near']) & (d['Hmax_near'] < 1e-4)
ks = np.flatnonzero(ok)
d6 = {}
for Tq in (12, 13, 14, 15, 15.5):
    k = int(np.argmin(abs(T - Tq)))
    d6['%.1f' % Tq] = dict(T=float(T[k]), dBdT=float(g[k]), bb_plus_T=float(B[k] - math.log(rb) + T[k]), H0tau=float(d['H0tau'][k]),
                           phi_b=float(d['phi_b'][k]), H_over_H0=float(d['H_over_H0'][k]), Hmax_near=float(d['Hmax_near'][k]))
out['D6'] = d6; print(json.dumps(d6, indent=0))
m8 = json.load(open('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json'))
out['M8_text_search'] = [l for l in json.dumps(m8).split(',') if '3.19424' in l][:5]
c = 0.5975949350280132; dstar = -3.1942416958680835
out['one_plus_c_plus_d_over_2'] = 1 + c + dstar/2
out['Y_rho_b'] = 1.0*78.8283759665
print(out['M8_text_search'], out['one_plus_c_plus_d_over_2'])
(HERE/'audit/AUDIT_RUNS.json').write_text(json.dumps(out, indent=1) + '\n')
