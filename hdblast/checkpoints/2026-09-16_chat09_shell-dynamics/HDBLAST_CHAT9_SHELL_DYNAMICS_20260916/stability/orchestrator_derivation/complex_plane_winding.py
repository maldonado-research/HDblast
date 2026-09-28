#!/usr/bin/env python3
"""Argument-principle count of zeros of the shell boundary mismatch B(mu2) in the complex mu2 plane
(checks that the eigenvalue-dependent boundary condition produces no complex-conjugate unstable pairs)."""
import sys, json, math
import numpy as np
sys.path.insert(0, "../../background")
import hdblast_background as bg

def B_complex(mu2, phi_h, y_b, t, y0=1e-3, n=8000):
    mu2 = np.asarray(mu2, complex); M = mu2.size
    alpha = 0.5 + np.sqrt(2.25 - mu2)             # principal branch: decaying (normalizable) cone behaviour, Re alpha > 1/2
    b = bg.series(phi_h, y0)
    rho, phi, s = b; psi = np.ones(M, complex); dpsi = alpha/y0
    Y = [np.full(M, rho, float), np.full(M, phi, float), np.full(M, s, float), psi, dpsi]
    def rhs(Y):
        rho, phi, s, psi, dpsi = Y
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
        spp = bg.U1(phi) - 4*Hh*s; g = spp/s
        dd = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
        return [rp, s, spp, dpsi, dd]
    add = lambda Y, K, a: [y + a*k for y, k in zip(Y, K)]
    dv = (y_b - y0)/n
    for _ in range(n):
        k1 = rhs(Y); k2 = rhs(add(Y, k1, dv/2)); k3 = rhs(add(Y, k2, dv/2)); k4 = rhs(add(Y, k3, dv))
        Y = [y + dv/6*(a + 2*b_ + 2*c_ + d_) for y, a, b_, c_, d_ in zip(Y, k1, k2, k3, k4)]
    rho, phi, s, psi, dpsi = Y
    rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
    spp = bg.U1(phi) - 4*Hh*s; g = spp/s
    dd = -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi
    dHh = -1/rho**2 - s*s/3
    chi = -3*(dpsi + 2*Hh*psi)/s
    dchi = -3*((dd + 2*dHh*psi + 2*Hh*dpsi)/s - spp*(dpsi + 2*Hh*psi)/s**2)
    Bc = dchi + 2*s*psi + 0.5*(2*bg.W2(phi))*chi
    return Bc/psi                                   # divide by psi(y_b): removes the growth of the solution, zeros unchanged unless psi_b=0

sol = json.load(open("../../background/REGISTERED_SHELL_FLOAT_SOLUTION.json"))
def winding(path):
    v = B_complex(path, sol["phi_h"], sol["v_b"], sol["t"])
    ph = np.unwrap(np.angle(v)); return (ph[-1] - ph[0])/(2*math.pi), v
res = {}
for (x0, x1, yy) in [(-60.0, 2.0, 40.0), (-400.0, 2.0, 300.0)]:
    m = 400
    bottom = np.linspace(x0, x1, m) - 1j*yy; right = x1 + 1j*np.linspace(-yy, yy, m)
    top = np.linspace(x1, x0, m) + 1j*yy; left = x0 + 1j*np.linspace(yy, -yy, m)
    path = np.concatenate([bottom, right, top, left, bottom[:1]])
    w, v = winding(path)
    print("rectangle Re in [%g,%g], |Im|<=%g : winding number of B/psi_b = %.6f  (min |B/psi_b| on contour = %.3e)" % (x0, x1, yy, w, np.abs(v).min()))
    res["rect_%g_%g_%g" % (x0, x1, yy)] = w
json.dump(res, open("COMPLEX_PLANE_WINDING.json", "w"), indent=1)
