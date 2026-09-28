#!/usr/bin/env python3
"""5D check of the SECOND static-shell branch predicted by the 4D effective theory (transcritical structure at d0)."""
import sys, math, json
import numpy as np
PREV = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916"
sys.path.insert(0, PREV + "/background"); sys.path.insert(0, PREV + "/stability/orchestrator_derivation")
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
    for it in range(60):
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
E = {r["d"]: r["stationary"] for r in json.load(open("EFT_BRANCHES.json"))["rows"]}
out = []
for d in [0.9, 1.0, 1.05, 1.3, 1.6, 2.0]:
    for q in E[d]:
        pb = q["phi_b"]
        guess = (math.log10(8.7855e-7), 8.2282 + math.atanh(pb))
        try:
            sol = solve(t, c, d, guess)
        except Exception as e:
            print("d=%.2f target phi_b=%+.4f : solve failed (%s)" % (d, pb, e)); continue
        f = lambda m: shoot_vec(m, sol["phi_h"], sol["y_b"], t, c, d2=d, n=6000)
        grid = np.concatenate([np.linspace(-12, 2.0, 100), np.linspace(2.0, 2.2499, 30)[1:]])
        roots = roots_from_grid(grid, f(grid), f)
        print("d=%.2f  EFT: phi_b=%+.5f mu2=%+.4f | 5D: phi_b=%+.5f h=%.6e resid=(%.0e,%.0e) mu2=%s" % (
            d, pb, q["mu2"], sol["phi_b"], 1/sol["rho_b"]**2, sol["resid"][0], sol["resid"][1], [round(float(r), 5) for r in roots]), flush=True)
        out.append(dict(d=d, eft=q, phi_b_5D=sol["phi_b"], h_5D=1/sol["rho_b"]**2, mu2_5D=[float(r) for r in roots]))
json.dump(out, open("BRANCHES_5D_CHECK.json", "w"), indent=1)
