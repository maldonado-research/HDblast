#!/usr/bin/env python3
"""Stability window of the quadratic-detuning family sigma_t = 2W + t(1 + c* phi + d phi^2/2) at t = 1e-3:
where does the modulus stop being tachyonic (mu2 = 0) and where does it merge into the continuum (mu2 = 9/4)?"""
import sys, math, json
import numpy as np
sys.path.insert(0, "../../background")
import hdblast_background as bg
from t_scan_orchestrator import shoot_vec, roots_from_grid
from stable_shell_recheck import solve
Ip = bg.I_plus(); c = 2/Ip - 4/3; t = 1e-3; ZE = c/2 + 3*c*c/8; d0 = 1.1134965631668923
rows = []
for d in [1.05, 1.10, 1.11, 1.12, 1.13, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.42, 1.43, 1.44, 1.45, 1.46, 1.50]:
    sol = solve(t, c, d, (math.log10(8.7855e-7), 8.2282))
    f = lambda m: shoot_vec(m, sol["phi_h"], sol["y_b"], t, c, d2=d, n=6000)
    grid = np.concatenate([np.linspace(-3, 2.0, 60), np.linspace(2.0, 2.2499, 40)[1:]])
    roots = roots_from_grid(grid, f(grid), f)
    rows.append(dict(d=d, phi_b=sol["phi_b"], mu2_5D=[float(r) for r in roots], mu2_EFT_t0=3*(d - d0)/ZE))
    print("d=%.2f  phi_b=%+.3e  5D mu2=%s  EFT(t->0)=%+.5f" % (d, sol["phi_b"], [round(float(r), 6) for r in roots], 3*(d - d0)/ZE), flush=True)
json.dump(dict(t=t, c=c, d0_t0=d0, Z_E=ZE, d_merge_EFT_t0=d0 + 0.75*ZE, rows=rows), open("STABILITY_WINDOW_SCAN.json", "w"), indent=1)
