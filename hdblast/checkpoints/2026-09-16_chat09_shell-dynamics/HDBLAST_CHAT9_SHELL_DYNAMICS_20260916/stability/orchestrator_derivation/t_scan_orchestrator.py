#!/usr/bin/env python3
"""Vectorised (over mu2) shooting; scan of detuning t at c = c_star, to test the closed-form HJ limit."""
import sys, math, json
import numpy as np
sys.path.insert(0, "../../background")
import hdblast_background as bg

def shoot_vec(mu2, phi_h, y_b, t, c, d2=0.0, y0=1e-3, n=12000):
    mu2 = np.asarray(mu2, float); M = mu2.size
    alpha = 0.5 + np.sqrt(2.25 - mu2)
    b = bg.series(phi_h, y0)
    Y = np.zeros((5, M)); Y[0] = b[0]; Y[1] = b[1]; Y[2] = b[2]; Y[3] = 1.0; Y[4] = alpha/y0
    def rhs(Y):
        rho, phi, s, psi, dpsi = Y
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
        spp = bg.U1(phi) - 4*Hh*s; g = spp/s
        ddpsi = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
        return np.array([rp, s, spp, dpsi, ddpsi])
    dv = (y_b - y0)/n
    for _ in range(n):
        k1 = rhs(Y); k2 = rhs(Y + dv/2*k1); k3 = rhs(Y + dv/2*k2); k4 = rhs(Y + dv*k3)
        Y = Y + dv/6*(k1 + 2*k2 + 2*k3 + k4)
    rho, phi, s, psi, dpsi = Y
    rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
    spp = bg.U1(phi) - 4*Hh*s; g = spp/s
    ddpsi = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
    dHh = -1/rho**2 - s*s/3
    chi = -3*(dpsi + 2*Hh*psi)/s
    dchi = -3*((ddpsi + 2*dHh*psi + 2*Hh*dpsi)/s - spp*(dpsi + 2*Hh*psi)/s**2)
    sig2 = 2*bg.W2(phi) + t*d2
    Bc = dchi + 2*s*psi + 0.5*sig2*chi
    return Bc/(np.abs(dchi) + np.abs(2*s*psi) + 1e-300)

def roots_from_grid(grid, vals, f):
    out = []
    for i in range(len(grid)-1):
        if vals[i]*vals[i+1] < 0:
            a, b_ = grid[i], grid[i+1]
            for _ in range(6):                       # vectorised bracketing refinement: 33 points per sweep
                g = np.linspace(a, b_, 33); v = f(g)
                j = np.where(v[:-1]*v[1:] < 0)[0][0]; a, b_ = g[j], g[j+1]
            out.append(0.5*(a + b_))
    return out

if __name__ == "__main__":
    Ip = bg.I_plus(); c = 2/Ip - 4/3
    pred = -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))
    print("I_plus=%.13f c_star=%.13f  HJ closed-form t->0 limit mu2 = %.9f" % (Ip, c, pred), flush=True)
    out = []
    for t in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
        guess = (math.log10(0.2207*t**1.8), 8.2282 + 0.9*math.log(1e-3/t))
        sol = bg.solve_shell(t, c, guess, dv=4e-4)
        res = {}
        for n in (6000, 12000):
            f = lambda m: shoot_vec(m, sol["phi_h"], sol["v_b"], t, c, n=n)
            grid = np.linspace(-14, 2.2, 82); vals = f(grid)
            res[n] = roots_from_grid(grid, vals, f)
        print("t=%.0e h=%.6e phi_b=%+.3e y_b=%.5f resid=(%.1e,%.1e) roots(n=6000)=%s roots(n=12000)=%s  (mu2-pred)/t=%s" % (
            t, sol["h"], sol["phi_b"], sol["v_b"], sol["residual"][0], sol["residual"][1], res[6000], res[12000],
            [(r - pred)/t for r in res[12000]]), flush=True)
        out.append(dict(t=t, h=sol["h"], phi_b=sol["phi_b"], y_b=sol["v_b"], mu2_n6000=res[6000], mu2_n12000=res[12000]))
    json.dump(dict(I_plus=Ip, c_star=c, hj_closed_form_limit=pred, scan=out), open("T_SCAN_ORCHESTRATOR.json", "w"), indent=1)
