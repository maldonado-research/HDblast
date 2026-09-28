#!/usr/bin/env python3
"""Chat 12 - floating-point spot check (system python3 + numpy) of the tensor sector on the actual backgrounds:
   h'' + 4H h' + m2 h/rho^2 = 0,  regular branch h ~ y^(-3/2 + sqrt(9/4 - m2)),  shell condition h'(y_b) = 0.
   Reports D(m2) = y_b h'(y_b)/h(y_b) on a grid of m2 in [-20, 2.2499] and min over the bulk of V2 - 9/4 = rho^2(9 phi'^2/16 - U/8) / rho^2."""
import sys, math, json
import numpy as np
P146 = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916"
sys.path.insert(0, P146 + "/background"); import hdblast_background as bg
Ip = bg.I_plus(); c_star = 2/Ip - 4/3
def shell(t, c, d, guess):
    def res(x):
        rho, phi, s = bg.integrate(-1 + 10**x[0], x[1], 4e-4)
        rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
        return np.array([rp/rho - (2*bg.W(phi) + t*(1 + c*phi + d*phi*phi/2))/6, s + (2*bg.W1(phi) + t*(c + d*phi))/2])
    x = np.array(guess, float)
    for _ in range(40):
        F = res(x); J = np.zeros((2, 2))
        for j in range(2):
            xp = x.copy(); xp[j] += 1e-7; xm = x.copy(); xm[j] -= 1e-7; J[:, j] = (res(xp) - res(xm))/2e-7
        dx = np.linalg.solve(J, -F); x = x + dx
        if np.linalg.norm(dx) < 1e-13: break
    return -1 + 10**x[0], x[1]
def D_of_m2(m2, phi_h, y_b, y0=1e-3, n=32000):   # n=8000 gives a spurious RK4 'node' for m2 -> 9/4 (first step ~ y0); converged for n >= 16000
    m2 = np.asarray(m2, float); aexp = -1.5 + np.sqrt(2.25 - m2)
    b = bg.series(phi_h, y0); M = m2.size
    Y = [np.full(M, b[0]), np.full(M, b[1]), np.full(M, b[2]), np.ones(M), aexp/y0]
    def rhs(Y):
        rho, phi, s, h, dh = Y
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho
        return [rp, s, bg.U1(phi) - 4*H*s, dh, -4*H*dh - m2*h/rho**2]
    add = lambda Y, K, a: [u + a*k_ for u, k_ in zip(Y, K)]
    dv = (y_b - y0)/n; nodes = np.zeros(M, int)
    for _ in range(n):
        k1 = rhs(Y); k2 = rhs(add(Y, k1, dv/2)); k3 = rhs(add(Y, k2, dv/2)); k4 = rhs(add(Y, k3, dv))
        Ynew = [u + dv/6*(a + 2*b_ + 2*c_ + d_) for u, a, b_, c_, d_ in zip(Y, k1, k2, k3, k4)]
        nodes += (Ynew[3]*Y[3] < 0); Y = Ynew
    return y_b*Y[4]/Y[3], nodes
out = {}
for name, (t, c, d) in {"S_8/5": (1e-3, 5975949350280/1e13, 1.6), "registered (d=0)": (1e-3, c_star, 0.0)}.items():
    phi_h, y_b = shell(t, c, d, (math.log10(8.7855e-7), 8.2282))
    grid = np.concatenate([np.linspace(-20, -0.05, 60), [0.0], np.linspace(0.02, 2.2499, 120)])
    Dv, nodes = D_of_m2(grid, phi_h, y_b)
    sign_changes = [(float(grid[i]), float(grid[i+1])) for i in range(len(grid)-1) if Dv[i]*Dv[i+1] < 0]
    i0 = int(np.argmin(np.abs(grid)))
    yend, rec = bg.integrate(phi_h, y_b, 4e-4, record=True)
    gap = 9*rec[:, 3]**2/16 - bg.U(rec[:, 2])/8
    pos = grid > 0.01; neg = grid < -0.01
    out[name] = {"phi_h": phi_h, "y_b": y_b, "D_at_m2_0": float(Dv[i0]),
                 "D_min_for_0<m2<9/4": float(Dv[pos].min()), "D_max_for_0<m2<9/4": float(Dv[pos].max()),
                 "D_range_for_m2<0": [float(Dv[neg].min()), float(Dv[neg].max())], "sign_changes": sign_changes, "max_nodes": int(nodes.max()),
                 "min_over_bulk_of_(V2-9/4)/rho^2": float(gap.min()), "max_phi_in_bulk": float(rec[:, 2].max())}
    print(name, json.dumps(out[name]), flush=True)
json.dump(out, open("TENSOR_FLOAT_SPOTCHECK.json", "w"), indent=1)
