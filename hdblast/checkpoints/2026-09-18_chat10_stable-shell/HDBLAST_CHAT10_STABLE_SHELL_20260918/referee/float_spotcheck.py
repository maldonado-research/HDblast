#!/usr/bin/env python3
"""Referee float spot-check (GUIDANCE ONLY, not a proof).  Independent wrapper: own RK4 loop (n=20000, y0=2e-3, differs from
Task B / Chat 9 settings), tracks sign changes of psi, scans m(mu^2) = chi' + 2 phi' psi + sigma_t'' chi/2 on -200 <= mu^2 < 9/4."""
import sys, json, numpy as np
sys.path.insert(0, "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/background")
import hdblast_background as bg
T = 1e-3
def scan(phi_h, y_b, d, mu2, y0=2e-3, n=20000):
    mu2 = np.asarray(mu2, float); M = mu2.size
    b = bg.series(phi_h, y0)
    Y = np.zeros((5, M)); Y[0], Y[1], Y[2] = b[0], b[1], b[2]; Y[3] = 1.0; Y[4] = (0.5 + np.sqrt(2.25 - mu2))/y0
    def rhs(Y):
        rho, phi, s, psi, dpsi = Y
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho
        spp = bg.U1(phi) - 4*H*s; g = spp/s
        return np.array([rp, s, spp, dpsi, -2*(H - g)*dpsi - (-(4/3)*s*s - 4*H*g + (2 + mu2)/rho**2)*psi])
    h = (y_b - y0)/n; Z = np.zeros(M, int); minphip = 1e9
    for _ in range(n):
        k1 = rhs(Y); k2 = rhs(Y + h/2*k1); k3 = rhs(Y + h/2*k2); k4 = rhs(Y + h*k3)
        Yn = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
        Z += (np.sign(Yn[3]) != np.sign(Y[3])); Y = Yn
        Y[3:] /= np.abs(Y[3]) + 1e-300   # rescale (linear eq.) to avoid overflow for very negative mu^2
        minphip = min(minphip, Y[2].min())
    rho, phi, s, psi, dpsi = Y
    rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho
    G = bg.U1(phi)/s; X = -(dpsi + 2*H*psi); R = X/psi
    B = G - 4*H + 2*phi + T*d/2
    bmis = (mu2 + 4)/rho**2 + B*R          # m = (3 psi/s) b
    return dict(Z=Z, b=bmis, B=float(B[0]), rho_b=float(rho[0]), phi_b=float(phi[0]), minphip=float(minphip),
                junc=[float(H[0] - (2*bg.W(phi[0]) + T*(1 + C*phi[0] + d*phi[0]**2/2))/6)])
C = 0.5975949350280
shells = {s["d"]: s for s in json.load(open("../oscillation_theorem/shells.json"))}
cases = [("S_8/5 (Task A certified centre)", 1.6, -0.99999912156299811104, 8.22819172487514216),
         ("d=1.6 (Task B float shell)", 1.6, shells[1.6]["phi_h"], shells[1.6]["y_b"]),
         ("d=1.3 (Task B float shell)", 1.3, shells[1.3]["phi_h"], shells[1.3]["y_b"])]
grid = np.concatenate([np.linspace(-200, -10, 96), np.linspace(-10, 2.2499, 400)[1:]])
out = {}
for name, d, ph, yb in cases:
    r = scan(ph, yb, d, grid)
    sc = np.where(r["b"][:-1]*r["b"][1:] < 0)[0]
    roots = []
    for i in sc:
        a, c_ = grid[i], grid[i+1]
        for _ in range(3):
            g = np.linspace(a, c_, 17); v = scan(ph, yb, d, g)["b"]; j = np.where(v[:-1]*v[1:] < 0)[0][0]; a, c_ = g[j], g[j+1]
        roots.append(0.5*(a + c_))
    out[name] = dict(B=r["B"], rho_b=r["rho_b"], phi_b=r["phi_b"], junction_G1_resid=r["junc"][0], min_phi_prime=r["minphip"], max_psi_zeros=int(r["Z"].max()),
                     n_sign_changes_of_b=len(sc), roots_mu2=roots, b_at_mu2_0=float(np.interp(0, grid, r["b"])), b_at_2p2499=float(r["b"][-1]),
                     b_max_on_grid=float(r["b"].max()), b_min_on_grid=float(r["b"].min()))
    print(name, json.dumps(out[name]))
json.dump(out, open("float_spotcheck_output.json", "w"), indent=1)
