#!/usr/bin/env python3
"""A1 diagnosis of the 'chart freeze' in the archived pilot runs (read-only).  Output: D0_FREEZE_DIAGNOSIS.json

In the conformal chart ds^2 = e^{2B}(-dt^2+dz^2)+..., shell fixed at z=0, the shell lapse is e^{B_b(t)}.
If the shell crosses the future light cone of the static solution's vertex (w = -U_K V_K = 0, a Rindler-type
horizon), the regular null coordinate is U_K ~ -e^{-u}, so e^{B_b} must decay like e^{-t}: d B_b/dt -> -1 exactly,
and the shell proper time saturates at a finite tau_max = tau(t) + e^{B_b(t)}/1 (in model units).
This script measures d b_b/dt (b_b = B_b - ln rho_b) at the end of each archived pilot run and the implied tau_max.
Numerical (centred differences of the saved 0.01-t samples)."""
import json, glob, hashlib
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PIL = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/runs')
out = dict(status='numerical (archived pilot time series, read-only)', runs={})
for f in sorted(PIL.glob('*_timeseries.npz')):
    tag = f.name.replace('_timeseries.npz', '')
    s = json.load(open(str(f).replace('_timeseries.npz', '_summary.json')))
    r = np.load(f)['rec']; t, bb, tau = r[:, 0], r[:, 3], r[:, 8]
    k = np.gradient(bb, t)
    # H0 tau = int e^{b_b} dt ; remaining proper time if e^{b_b} ~ e^{-(t - t_end)} continues:  e^{b_b(t_end)}/|k_end|
    kend = float(np.mean(k[-20:]))
    rem = float(np.exp(bb[-1])/abs(kend)) if kend < -0.05 else None
    out['runs'][tag] = dict(d_quad=s.get('d_quad', 0.0), Y=s['Y'], dc=s['dc'], t_end=float(t[-1]), H0tau_end=float(tau[-1]),
                            dbdt_end=kend, dbdt_min=float(k.min()), b_b_end=float(bb[-1]),
                            H0tau_max_extrapolated=(float(tau[-1]) + rem) if rem is not None else None,
                            H_over_H0_end=float(r[-1, 2]),
                            t_first_dbdt_below_m0p5=float(t[np.argmax(k < -0.5)]) if (k < -0.5).any() else None,
                            asymptotic_bb_plus_t=float(bb[-1] + t[-1]))
    print(tag, out['runs'][tag])
vals = [v['dbdt_end'] for v in out['runs'].values() if v['t_end'] - (v['t_first_dbdt_below_m0p5'] or 1e9) > 4]
out['summary'] = dict(runs_followed_4_or_more_time_units_past_onset=len(vals),
                      dbdt_end_values=vals, max_abs_deviation_from_minus_1=float(max(abs(np.array(vals) + 1))) if vals else None,
                      interpretation='d B_b/dt -> -1 (to the precision listed) is the signature of the shell reaching the future light cone '
                                     'of the static vertex (Rindler-type horizon) at finite proper time; the conformal chart with the shell '
                                     'at z=0 cannot cover the shell beyond tau_max. Continuing requires a chart regular across that '
                                     'light cone (null coordinate ~ U_K), not a relabelling of the old slices.')
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'D0_FREEZE_DIAGNOSIS.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out['summary'], indent=1))
