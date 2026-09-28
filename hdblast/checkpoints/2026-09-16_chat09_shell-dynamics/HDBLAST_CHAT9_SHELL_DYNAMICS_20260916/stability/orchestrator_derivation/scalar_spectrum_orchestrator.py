#!/usr/bin/env python3
"""Orchestrator's own longitudinal-gauge scalar-sector shooting code (independent of the workflow agents).
Master equation (y = proper distance from the cone, unit-dS harmonics Box Y = mu2 Y):
  psi'' + 2(rho'/rho - phi''/phi') psi' + [ -(4/3) phi'^2 - 4 (rho'/rho)(phi''/phi') + (2+mu2)/rho^2 ] psi = 0
  chi = -3 (psi' + 2 (rho'/rho) psi)/phi'                      (momentum constraint, xi = -2 psi)
Shell BC (no bending for mu2 != -4):  chi' + 2 phi' psi + (sigma_t''/2) chi = 0.
Cone: psi ~ y^alpha, alpha = 1/2 + sqrt(9/4 - mu2).
"""
import sys, math, json
import numpy as np
sys.path.insert(0, "../../background")
import hdblast_background as bg

def full_rhs(Y, mu2):
    rho, phi, s, psi, dpsi = Y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
    Hh = rp/rho
    spp = bg.U1(phi) - 4*Hh*s          # phi''
    g = spp/s
    ddpsi = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
    return np.array([rp, s, spp, dpsi, ddpsi])

def shoot(mu2, phi_h, y_b, t, c, d2=0.0, y0=1e-3, n=40000):
    alpha = 0.5 + math.sqrt(2.25 - mu2)
    b = bg.series(phi_h, y0)
    Y = np.array([b[0], b[1], b[2], y0**alpha, alpha*y0**(alpha-1)])
    # normalise to avoid overflow
    Y[3:] /= y0**alpha
    dv = (y_b - y0)/n
    for _ in range(n):
        k1 = full_rhs(Y, mu2); k2 = full_rhs(Y + dv/2*k1, mu2); k3 = full_rhs(Y + dv/2*k2, mu2); k4 = full_rhs(Y + dv*k3, mu2)
        Y = Y + dv/6*(k1 + 2*k2 + 2*k3 + k4)
    rho, phi, s, psi, dpsi = Y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
    spp = bg.U1(phi) - 4*Hh*s; g = spp/s
    ddpsi = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
    dHh = -1/rho**2 - s*s/3
    chi = -3*(dpsi + 2*Hh*psi)/s
    dchi = -3*((ddpsi + 2*dHh*psi + 2*Hh*dpsi)/s - spp*(dpsi + 2*Hh*psi)/s**2)
    sig2 = 2*bg.W2(phi) + t*d2       # sigma_t'' (d2 = optional quadratic detuning coefficient)
    Bc = dchi + 2*s*psi + 0.5*sig2*chi
    scale = abs(dchi) + abs(2*s*psi) + abs(0.5*sig2*chi) + 1e-300
    return Bc/scale, (psi, dpsi, chi, dchi)

if __name__ == "__main__":
    sol = json.load(open("../../background/REGISTERED_SHELL_FLOAT_SOLUTION.json"))
    t, c, phi_h, y_b = sol["t"], sol["c"], sol["phi_h"], sol["v_b"]
    grid = np.concatenate([np.linspace(-12, 2.2, 72)])
    vals = []
    for m in grid:
        B, _ = shoot(m, phi_h, y_b, t, c, n=8000)
        vals.append(B); print("mu2=%8.4f  B=%+.6e" % (m, B))
    roots = []
    for i in range(len(grid)-1):
        if vals[i]*vals[i+1] < 0:
            a, b_ = grid[i], grid[i+1]; fa = vals[i]
            for _ in range(60):
                mid = 0.5*(a+b_); fm, _ = shoot(mid, phi_h, y_b, t, c, n=8000)
                if fa*fm <= 0: b_ = mid
                else: a, fa = mid, fm
            roots.append(0.5*(a+b_))
    print("sign-change roots (n=8000):", roots)
    # refine with more steps
    for r in roots:
        for n in (16000, 32000):
            a, b_ = r-1e-3, r+1e-3
            fa, _ = shoot(a, phi_h, y_b, t, c, n=n)
            for _ in range(50):
                mid = 0.5*(a+b_); fm, _ = shoot(mid, phi_h, y_b, t, c, n=n)
                if fa*fm <= 0: b_ = mid
                else: a, fa = mid, fm
            print("  refined root n=%d : mu2 = %.10f" % (n, 0.5*(a+b_)))
