#!/usr/bin/env python3
"""B2 diagnosis D7: where and how the delta = 1e-3 runs end, compared with the A1 delta = 0.1 runs that reached a plateau.
For every run: F_inf (= x_c, the chart's light-cone time), T_end, H0tau_end, stop reason; phi_b, H/H0, lapse at T = F_inf;
the growth rate of ln(lapse) per unit chart time over [F_inf, T_end] (median of the finite-difference slope);
the chart time spent beyond F_inf; the maximum phi_b reached.  Also the old-chart continuation run D6 (no stop): the asymptotic
slope d B_b/dT, the frozen shell proper time H0tau_freeze, phi_b and H/H0 at the freeze, and the asymptotic b_b + T (the x_c the
criterion approximates).  Output diag/D7_CHART_END.json.  Numerical (reads saved time series only)."""
import json, glob, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
def load(f):
    s = json.load(open(f)); z = np.load(f.replace('_summary.json', '_timeseries.npz')); cols = list(z['cols']); r = z['rec']
    return s, {c: r[:, i] for i, c in enumerate(cols)}
out = dict(status='numerical', purpose=__doc__, runs={}, old_chart_continuation={})
files = sorted(glob.glob(str(HERE/'runs/main/*_summary.json')) + glob.glob(str(HERE/'runs/ctl/*_summary.json')) + glob.glob(str(HERE/'runs/explore/*_summary.json')) + glob.glob(str(HERE/'runs/res/*_summary.json')))
files += [str(A1/('main_dstar_Y1_dc1e-%d_dzf%s_summary.json' % (k, z))) for k in (2, 4) for z in ('1e-3', '5e-4')]
for f in files:
    s, d = load(f); tag = ('A1_delta0.1/' if 'A1_tuned' in f else '') + Path(f).name.replace('_summary.json', '')
    Fi = s.get('F_inf'); T = d['T']; lap = d['lapse']
    e = dict(delta=s['params']['delta'], Y=s['params']['Y'], dc=s['params']['dc'], dzf=s['params']['dzf'], xc=s['params']['xc'], F_inf=Fi,
             stop_reason=s['stop_reason'], T_end=s['T_end'], H0tau_end=s['H0tau_end'], phi_b_max=float(np.nanmax(d['phi_b'])),
             H_over_H0_end=float(d['H_over_H0'][-1]), runtime_s=s['runtime_s'])
    if Fi is not None and math.isfinite(Fi) and T[-1] > Fi:
        k = int(np.argmin(np.abs(T - Fi)))
        e.update(at_F_inf=dict(phi_b=float(d['phi_b'][k]), H_over_H0=float(d['H_over_H0'][k]), lapse=float(lap[k]), H0tau=float(d['H0tau'][k])),
                 T_beyond_F_inf=float(T[-1] - Fi))
        m = (T >= Fi) & np.isfinite(lap) & (lap > 0)
        if m.sum() > 5:
            sl = np.gradient(np.log(lap[m]), T[m]); e['dlnlapse_dT_beyond_F_inf_median'] = float(np.median(sl))
    out['runs'][tag] = e
for f in sorted(glob.glob(str(HERE/'runs/d6/*_summary.json'))):
    s, d = load(f); T = d['T']; B = d['B_b']; g = np.gradient(B, T); rb = s['rho_b']
    good = np.isfinite(d['Hmax_near']) & (d['Hmax_near'] < 1e-4)
    Tgood = T[good].max()
    k = int(np.argmin(np.abs(T - Tgood)))
    out['old_chart_continuation'][Path(f).name.replace('_summary.json', '')] = dict(
        last_T_with_near_shell_residual_lt_1em4=float(Tgood), dBdT_there=float(g[k]), H0tau_there=float(d['H0tau'][k]),
        phi_b_there=float(d['phi_b'][k]), H_over_H0_there=float(d['H_over_H0'][k]), lapse_there=float(d['lapse'][k]),
        b_b_plus_T_there=float(B[k] - math.log(rb) + T[k]), stop_reason=s['stop_reason'], H0tau_end=s['H0tau_end'])
(HERE/'diag/D7_CHART_END.json').write_text(json.dumps(out, indent=1) + '\n')
for k, v in out['runs'].items():
    print('%-44s xc=%-5s Tend=%.2f beyond=%s tau_end=%.2f phimax=%.3f dlnlapse/dT=%s %s' % (k, v['xc'], v['T_end'], round(v.get('T_beyond_F_inf', float('nan')), 2), v['H0tau_end'], v['phi_b_max'], round(v.get('dlnlapse_dT_beyond_F_inf_median', float('nan')), 2), v['stop_reason'][:40]))
print(json.dumps(out['old_chart_continuation'], indent=1))
