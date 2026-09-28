#!/usr/bin/env python3
"""Check of the curvature expansion of the BACKGROUND against the full 5D cone-to-shell solution (t = 1e-3, c = c_star):
   rho'/rho = W/3 + X I + X^2 z2 + O(X^3),   phi' = -W_phi + X p1 + X^2 p2 + O(X^3),   X = 1/rho^2,  z2 = -l_-^3/8 = -0.729."""
import sys, json, math
import numpy as np
HERE = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/"
sys.path.insert(0, HERE + "background")
import hdblast_background as bg, slope_formula as sf
sol = json.load(open(HERE + "background/REGISTERED_SHELL_FLOAT_SOLUTION.json"))
_, tr = bg.integrate(sol["phi_h"], sol["v_b"], dv=2e-4, record=True)
idx = [int(len(tr)*f) for f in (0.55, 0.65, 0.75, 0.85, 0.92, 0.97)] + [len(tr) - 1]
phis = [float(tr[i, 2]) for i in idx]; tab = sf.wall(phis); rows = []
print("   y        phi          X=1/rho^2    dH=H-W/3      dH-XI        dH-XI-X^2z2   (last)/X^3 |  F=phi'+W_phi   F-Xp1      F-Xp1-X^2p2  (last)/X^3")
for i, p in zip(idx, phis):
    v, rho, phi, s = tr[i]; X = 1/rho**2
    Hh = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))/rho
    L = sf.local(p, *tab[p], lambda q: (1.0, 0.0, 0.0))
    dH = Hh - bg.W(phi)/3; F = s + bg.W1(phi)
    r = (dH - X*L['I'] - X*X*L['z2'], F - X*L['p1'] - X*X*L['p2'])
    print("%8.4f %+.6e %.6e %.6e %+.3e %+.3e %+.3f | %+.6e %+.3e %+.3e %+.3f" % (v, phi, X, dH, dH - X*L['I'], r[0], r[0]/X**3, F, F - X*L['p1'], r[1], r[1]/X**3))
    rows.append(dict(y=v, phi=phi, X=X, dH=dH, dH_res1=dH - X*L['I'], dH_res2=r[0], F=F, F_res1=F - X*L['p1'], F_res2=r[1]))
json.dump(rows, open("BACKGROUND_EXPANSION_CHECK.json", "w"), indent=1)
