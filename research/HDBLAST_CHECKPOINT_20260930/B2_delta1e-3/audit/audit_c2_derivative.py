#!/usr/bin/env python3
"""AUDIT (independent of the producer's estimators): growth rate of the C2 runs from the log-slope of d(dev)/dtau,
which removes any constant (neutral-mode) offset B in dev = A exp(lam tau) + B exactly, with no fit of B.
Also: registered-estimator window rates recomputed independently (np.polyfit of ln|dev| vs H0tau) to check the producer's numbers.
Reads saved runs/cal/*_timeseries.npz only.  Output audit/AUDIT_C2.json.  Numerical."""
import json, glob
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
TARGET = 1.65719
out = dict(status='numerical (audit)', target=TARGET, predicted_gap_ratio=float(np.exp(-TARGET/2)), runs={})
for f in sorted(glob.glob(str(HERE/'runs/cal/*_summary.json'))):
    s = json.load(open(f)); tag = Path(f).name.replace('_summary.json', '')
    z = np.load(f.replace('_summary.json', '_timeseries.npz')); cols = list(z['cols']); rec = z['rec']
    tau = rec[:, cols.index('H0tau')]; dev = rec[:, cols.index('phi_b')] - s['phi_b_static']
    k = np.concatenate(([True], np.diff(tau) > 1e-9)); tau, dev = tau[k], dev[k]
    dev = dev*np.sign(dev[-1])
    reg = {}
    for a in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0):
        m = (tau >= a) & (tau <= a+1) & (np.abs(dev) < 1e-3) & (dev > 0)
        if m.sum() > 10 and tau[m][-1] > a+0.95: reg['%.1f' % a] = float(np.polyfit(tau[m], np.log(dev[m]), 1)[0])
    regset = [reg[x] for x in ('2.0', '2.5', '3.0', '3.5', '4.0') if x in reg]
    # derivative estimator: centered differences on the (non-uniform) record grid
    ddev = np.gradient(dev, tau)
    der = {}
    for a in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0):
        m = (tau >= a) & (tau <= a+1) & (ddev > 0) & (np.abs(dev) < 1e-3)
        if m.sum() > 10 and tau[m][-1] > a+0.95: der['%.1f' % a] = float(np.polyfit(tau[m], np.log(ddev[m]), 1)[0])
    out['runs'][tag] = dict(dc=s['params']['dc'], dzf=s['params']['dzf'], grid=s['params'].get('grid_kind'),
        registered_windows=reg, registered_rate=float(np.mean(regset)) if len(regset) == 5 else None,
        registered_spread=float(np.ptp(regset)) if len(regset) == 5 else None,
        derivative_windows=der,
        derivative_rate_2to5=float(np.mean([der[x] for x in ('2.0','2.5','3.0','3.5','4.0') if x in der])) if der else None,
        derivative_spread_2to5=float(np.ptp([der[x] for x in ('2.0','2.5','3.0','3.5','4.0') if x in der])) if der else None)
    r = out['runs'][tag]; print(tag, r['registered_rate'], r['registered_spread'], r['derivative_rate_2to5'], r['derivative_spread_2to5'])
    print('   deriv windows', {k: round(v, 5) for k, v in der.items()})
(HERE/'audit/AUDIT_C2.json').write_text(json.dumps(out, indent=1) + '\n')
