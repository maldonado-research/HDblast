#!/usr/bin/env python3
"""B2 diagnosis D4 (post-registration, supplementary to C2): model the C2 growth curve as an unstable mode plus a
neutral (constant) offset plus the leading nonlinear term,
    dev(tau) = A exp(lam tau) + B + C exp(2 lam tau),     dev = phi_b - phi_b(static, tension c),
fitted by nonlinear least squares in log space on tau in [t0, t1] (records with |dev| < 2e-3).  Fits on several windows test
stability.  The pure-exponential window slopes (registered estimator) are biased low at early tau by B (gap ~ exp(-lam tau),
ratio exp(-lam/2) = 0.437 per 0.5 in tau) and at late tau by C.  Output diag/D4_C2_MODE_FIT.json.  Numerical."""
import json, glob, sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
HERE = Path(__file__).resolve().parent.parent
TARGET = 1.65719
out = dict(status='numerical', purpose=__doc__, target=TARGET, runs={})
for f in sorted(glob.glob(str(HERE/'runs/cal/*_summary.json')) + glob.glob(str(HERE/'runs/diag/*_summary.json'))):
    s = json.load(open(f)); tag = Path(f).name.replace('_summary.json', '')
    z = np.load(f.replace('_summary.json', '_timeseries.npz')); cols = list(z['cols']); rec = z['rec']
    tau = rec[:, cols.index('H0tau')]; dev = rec[:, cols.index('phi_b')] - s['phi_b_static']
    keep = np.concatenate(([True], np.diff(tau) > 1e-9)); tau, dev = tau[keep], dev[keep]
    sg = np.sign(dev[-1]); dev = sg*dev
    res = {}
    for (t0, t1, model) in [(1.0, 5.0, 'AB'), (1.5, 5.5, 'AB'), (1.0, 7.0, 'ABC'), (1.5, 7.5, 'ABC'), (2.0, 7.0, 'ABC'), (3.0, 7.5, 'ABC'), (1.0, 6.0, 'ABC')]:
        m = (tau >= t0) & (tau <= t1) & (dev > 0) & (dev < 2e-3)
        if m.sum() < 50 or tau[m][-1] < t1 - 0.3: continue
        t, y = tau[m], dev[m]
        def f_(p):
            A, lam, B = np.exp(p[0]), p[1], p[2]*1e-8
            C = (p[3]*1e-3 if model == 'ABC' else 0.0)
            return np.log(np.abs(A*np.exp(lam*t) + B + C*(A*np.exp(lam*t))**2)) - np.log(y)
        p0 = [np.log(y[-1]) - 1.657*t[-1], 1.657, 0.0] + ([0.0] if model == 'ABC' else [])
        r = least_squares(f_, p0, x_scale='jac', xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=20000)
        res['%s %.1f-%.1f' % (model, t0, t1)] = dict(lam=float(r.x[1]), diff=float(r.x[1] - TARGET), B_over_dev0=float(r.x[2]*1e-8/max(dev[0], 1e-300)),
                                                     C_rel=(float(r.x[3]*1e-3) if model == 'ABC' else None), rms_log_resid=float(np.sqrt(np.mean(r.fun**2))), n=int(m.sum()))
    lams = [v['lam'] for v in res.values()]
    out['runs'][tag] = dict(dc=s['params']['dc'], grid=s['params'].get('grid_kind'), dzf=s['params']['dzf'], table_dx=s['params'].get('table_dx'),
                           points_per_wall=s.get('points_per_wall'), T_end=s['T_end'], dev0=float(dev[0]), fits=res,
                           lam_range=[min(lams), max(lams)] if lams else None)
    print(tag, out['runs'][tag]['lam_range'], {k: round(v['lam'], 6) for k, v in res.items()})
(HERE/'diag/D4_C2_MODE_FIT.json').write_text(json.dumps(out, indent=1) + '\n')
