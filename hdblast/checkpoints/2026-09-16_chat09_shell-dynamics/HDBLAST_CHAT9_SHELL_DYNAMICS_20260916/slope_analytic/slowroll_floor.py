#!/usr/bin/env python3
"""By-product: with sigma_t = 2W + t(1 + c_star phi + d phi^2/2) at FINITE t the O(t) (four-derivative) shift of the equilibrium
turns the tuned point d = d0 into an avoided saddle-node: a hilltop (phi_b > 0) and a minimum (phi_b < 0) coexist and the
hilltop tachyon mass cannot be tuned to zero.  Scan of d with the resummed predictor + one full 5D verification."""
import sys, json, math, numpy as np
import slope_formula as sf, resummed_predictor as rp, validate_5d as v5
r0 = sf.slope(0.0, sf.linear_detuning); c = r0['Vf'].c
tab = rp.Table(-0.03, 0.05); out = {}
for t in (1e-3, 2.5e-4):
    best = None; rows = []
    for d in np.linspace(1.100, 1.125, 101):
        Vf = lambda p, d=d: (1 + c*p + d*p*p/2, c + d*p, d)
        p = rp.predict(t, Vf, tab, 6e-3*math.sqrt(t/1e-3))
        if p['phi_b'] > 0 and p['mu2'] < 0:
            rows.append((float(d), p['phi_b'], p['mu2']))
            if best is None or p['mu2'] > best[2]: best = rows[-1]
    print("t=%g: least tachyonic hilltop: d=%.5f phi_b=%.5e mu2=%.6f   (mu2/sqrt(t) = %.4f)" % (t, *best, best[2]/math.sqrt(t)))
    out[str(t)] = dict(best=best, rows=rows)
t = 1e-3; d = out[str(t)]['best'][0]
Vf = lambda q: (1 + c*q + d*q*q/2, c + d*q, d)
sol = v5.solve(t, Vf, v5.guess(t, out[str(t)]['best'][1], r0['x1'], 0.0, v5.throat_kappa()))
f = lambda m: v5.shoot(m, sol['phi_h'], sol['y_b'], t, Vf, n=16000)
mu5 = v5.root(f, -0.15, 0.0)
print("5D check at t=1e-3, d=%.5f: phi_b=%.6e mu2=%.9f resid=%s" % (d, sol['phi_b'], mu5, sol['resid']))
out['check_5D'] = dict(t=t, d=d, phi_b=sol['phi_b'], mu2=float(mu5))
json.dump(out, open("SLOWROLL_FLOOR.json", "w"), indent=1)
