#!/usr/bin/env python3
"""Orchestrator re-run of the referee's finding: with quadratic detuning sigma_t = 2W + t(1 + c phi + d phi^2/2)
the dS shell has NO tachyon for d >= ~1.3.  Independent shell solve (own Newton with the d-term) + own shooting code."""
import sys, math, json
import numpy as np
sys.path.insert(0, "../../background")
import hdblast_background as bg
from t_scan_orchestrator import shoot_vec, roots_from_grid

def shell_res(x, t, c, d, dv=4e-4):
    phi_h = -1 + 10**x[0]
    rho, phi, s = bg.integrate(phi_h, x[1], dv)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
    sig = 2*bg.W(phi) + t*(1 + c*phi + d*phi*phi/2); dsig = 2*bg.W1(phi) + t*(c + d*phi)
    return np.array([rp/rho - sig/6, s + dsig/2]), (rho, phi)

def solve(t, c, d, guess):
    x = np.array(guess, float)
    for it in range(40):
        F, _ = shell_res(x, t, c, d); J = np.zeros((2, 2))
        for j in range(2):
            xp = x.copy(); xp[j] += 1e-7; xm = x.copy(); xm[j] -= 1e-7
            J[:, j] = (shell_res(xp, t, c, d)[0] - shell_res(xm, t, c, d)[0])/2e-7
        dx = np.linalg.solve(J, -F); lam = 1.0
        while lam > 1e-4 and np.linalg.norm(shell_res(x + lam*dx, t, c, d)[0]) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < 1e-13: break
    F, (rho, phi) = shell_res(x, t, c, d)
    return dict(phi_h=-1 + 10**x[0], y_b=x[1], rho_b=rho, phi_b=phi, resid=[float(F[0]), float(F[1])])

Ip = bg.I_plus(); c = 2/Ip - 4/3; t = 1e-3
ZE = c/2 + 3*c*c/8; d0 = 1.1134965631668923
out = []
for d in [0.0, 1.0, 1.3, 1.4, 1.6, 2.0]:
    sol = solve(t, c, d, (math.log10(8.7855e-7), 8.2282))
    f = lambda m: shoot_vec(m, sol["phi_h"], sol["y_b"], t, c, d2=d, n=6000)
    grid = np.linspace(-60, 2.2, 200); roots = roots_from_grid(grid, f(grid), f)
    eft = 3*(d - d0)/ZE
    print("d=%.2f  phi_b=%+.3e  resid=(%.1e,%.1e)  5D bound states mu2 = %s   EFT(t->0) 3(d-d0)/Z_E = %+.5f" % (
        d, sol["phi_b"], sol["resid"][0], sol["resid"][1], [round(float(r), 7) for r in roots], eft), flush=True)
    out.append(dict(d=d, phi_b=sol["phi_b"], mu2_5D=[float(r) for r in roots], mu2_EFT_t0=eft))
json.dump(out, open("STABLE_SHELL_RECHECK.json", "w"), indent=1)
