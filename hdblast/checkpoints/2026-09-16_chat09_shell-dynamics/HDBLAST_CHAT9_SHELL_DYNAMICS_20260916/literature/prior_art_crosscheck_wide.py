#!/usr/bin/env python3
"""Robustness controls for prior_art_crosscheck.py at the registered point t = 1e-3 (run: python3 prior_art_crosscheck_wide.py).
 1. wide scan of the Frolov-Kofman boundary mismatch for mu^2 in [-40, 2.2] (below the continuum threshold 9/4);
 2. step-halving and start-point (v0) control of the tachyonic eigenvalue.
Floating-point diagnostic only."""
import math, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prior_art_crosscheck as pc
bg = pc.bg
I = bg.I_plus(); c = 2/I - 4/3; t = 1e-3
res = {}
for dv in (4e-4, 2e-4):
    sol = bg.solve_shell(t, c, guess=(math.log10(8.7855e-7), 8.6), dv=dv)
    r = pc.eigen(sol['phi_h'], sol['v_b'], t, c, -9.0, -6.0, dv)
    print("dv = %g : mu^2 = %.8f" % (dv, r)); res["dv_%g" % dv] = r
sol = bg.solve_shell(t, c, guess=(math.log10(8.7855e-7), 8.6), dv=4e-4)
# v0 control: re-implement shoot start by monkeypatching default v0
for v0 in (2e-2, 5e-3):
    f = lambda m: pc.shoot(m, sol['phi_h'], sol['v_b'], t, c, 4e-4, v0)[0]
    lo, hi = -9.0, -6.0; flo = f(lo)
    for _ in range(45):
        mid = 0.5*(lo + hi); fm = f(mid)
        if flo*fm <= 0: hi = mid
        else: lo, flo = mid, fm
    print("v0 = %g (perturbation start; background series still started at same v0): mu^2 = %.8f" % (v0, 0.5*(lo + hi)))
    res["v0_%g" % v0] = 0.5*(lo + hi)
grid = [-40 + 0.4*i for i in range(0, 106)]
grid = [g for g in grid if g < 2.24]
vals = [pc.shoot(m, sol['phi_h'], sol['v_b'], t, c, 4e-4)[0] for m in grid]
sc = [(a, b) for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]) if fa*fb < 0]
print("sign changes of the shell mismatch for mu^2 in [-40, 2.2]:", sc)
res["sign_change_brackets"] = sc
res["mismatch_samples"] = list(zip(grid[::5], vals[::5]))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "prior_art_crosscheck_wide_output.json"), "w"), indent=1)
