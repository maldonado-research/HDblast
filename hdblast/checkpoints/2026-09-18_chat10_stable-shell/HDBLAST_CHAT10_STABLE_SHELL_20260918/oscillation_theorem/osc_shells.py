#!/usr/bin/env python3
"""FLOAT (non-rigorous) shell solutions for sigma_t = 2W + t(1 + c phi + d phi^2/2), t = 1e-3, c = 0.5975949350280,
for several d.  Output shells.json with B = phi''/phi' + sigma_t''/2, beta = rho^4 phi'^2 B, and a check that
phi' > 0 and phi < phi_0 = 0.4088 (root of U_phi) on (0, y_b] (hypotheses H2 of OSCILLATION_THEOREM.md).
RUNTIME: about 1 minute per shell (7-8 min for the 8 shells; result cached in shells.json, run once).
Shell Newton copied from Chat 9 stability/orchestrator_derivation/stable_shell_recheck.py (guidance only)."""
import sys, math, json
import numpy as np
BG = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/background"
sys.path.insert(0, BG)
import hdblast_background as bg
T = 1e-3; C = 0.5975949350280

def shell_res(x, t, c, d, dv=4e-4):
    phi_h = -1 + 10**x[0]
    rho, phi, s = bg.integrate(phi_h, x[1], dv)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
    sig = 2*bg.W(phi) + t*(1 + c*phi + d*phi*phi/2); dsig = 2*bg.W1(phi) + t*(c + d*phi)
    return np.array([rp/rho - sig/6, s + dsig/2]), (rho, phi, s, rp)

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
    F, (rho, phi, s, rp) = shell_res(x, t, c, d)
    H = rp/rho; spp = bg.U1(phi) - 4*H*s
    sig2 = 2*bg.W2(phi) + t*d
    B = spp/s + sig2/2
    _, traj = bg.integrate(-1 + 10**x[0], x[1], 4e-4, record=True)
    return dict(d=d, phi_h=-1 + 10**x[0], y_b=x[1], rho_b=rho, phi_b=phi, s_b=s, H_b=H, B=B, beta=rho**4*s*s*B,
                resid=[float(F[0]), float(F[1])], min_phi_prime_on_grid=float(traj[:, 3].min()),
                max_phi_on_grid=float(traj[:, 2].max()), min_phi_on_grid=float(traj[:, 2].min()))

if __name__ == "__main__":
    out = []
    for d in [0.0, 0.4, 0.7, 1.0, 1.3, 1.4, 1.6, 2.0]:
        r = solve(T, C, d, (math.log10(8.7855e-7), 8.2282)); out.append(r)
        print("d=%.2f phi_h+1=%.6e y_b=%.6f rho_b=%.4f phi_b=%+.4e phi'_b=%.4e B=%+.4e beta=%+.4e resid=%s min phi'=%.2e max phi=%+.2e"
              % (d, r["phi_h"] + 1, r["y_b"], r["rho_b"], r["phi_b"], r["s_b"], r["B"], r["beta"], r["resid"],
                 r["min_phi_prime_on_grid"], r["max_phi_on_grid"]), flush=True)
    json.dump(out, open("shells.json", "w"), indent=1)
