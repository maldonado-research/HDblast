#!/usr/bin/env python3
"""Referee test: is the tachyon a property of dS shells in this model, or of the registered LINEAR detuning?
EFT predicts mu^2 = 3(d-d0)/Z_E at phi_b=0 for sigma_t = 2W + t(1 + c phi + d phi^2/2). For d > d0 = 1.1135 the static
shell should be linearly STABLE.  We solve the full 5D problem (slope_analytic/validate_5d.py machinery) for a few d at t=1e-3
and scan mu^2 in [-60, 2.2] for bound states.  numpy only; float, not certified."""
import sys, json, os
import numpy as np
HERE = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/"
sys.path.insert(0, HERE + "slope_analytic"); sys.path.insert(0, HERE + "background")
import validate_5d as v5, slope_formula as sf
t = 1e-3; out = []
kap = v5.throat_kappa()
for d in (0.0, 1.0, 1.3, 1.4, 1.6):
    r = sf.slope(0.0, lambda p, I, I1: sf.linear_detuning(p, I, I1, d=d)); V = r['Vf']
    sol = v5.solve(t, V, v5.guess(t, 0.0, r['x1'], r['beta'], kap))
    f = lambda m: v5.shoot(m, sol['phi_h'], sol['y_b'], t, V, n=6000)
    g = np.linspace(-60, 2.2, 312); val = f(g)
    # normalise sign changes; refine each
    idx = np.where(val[:-1]*val[1:] < 0)[0]; roots = []
    for j in idx:
        roots.append(float(v5.root(f, g[j], g[j+1])))
    row = dict(d=d, mu0_EFT=r['mu0'], slope_EFT=r['slope'], EFT_pred=r['mu0'] + r['slope']*t, phi_b=sol['phi_b'], rho_b=sol['rho_b'],
               resid=sol['resid'], roots_5D=roots)
    out.append(row); print(json.dumps(row), flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "REFEREE_STABLE_SHELL_TEST.json"), "w"), indent=1)
