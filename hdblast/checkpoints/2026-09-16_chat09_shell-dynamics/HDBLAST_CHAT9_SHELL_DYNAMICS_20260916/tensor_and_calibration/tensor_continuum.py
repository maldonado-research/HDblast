#!/usr/bin/env python3
"""HDBLAST Chat 9 - KK continuum (m^2 = (9/4 + k^2) H^2) weight on the registered shell.
Integrates h_ss + (4 y rho'/rho - 1) h_s + m^2 (y/rho)^2 h = 0 in s = ln y inward from the shell (h=1, h_s=0) to y0 = 1e-4,
background taken from the forward RK4 solution by cubic Hermite interpolation.  u = rho^{3/2} h -> A cos(k z + delta) at the cone.
w(k) = (2/pi) rho_b^2 I_kept / A^2 = |u_k(shell)|^2 / |u_0(shell)|^2 per unit dk  (continuum normalised to delta(k-k'))."""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tensor_spectrum as ts
bg = ts.bg

def hermite(x, xg, f, fp):
    i = np.clip(np.searchsorted(xg, x) - 1, 0, len(xg) - 2)
    hh = xg[i+1] - xg[i]; tt = (x - xg[i])/hh
    h00 = 2*tt**3 - 3*tt**2 + 1; h10 = tt**3 - 2*tt**2 + tt; h01 = -2*tt**3 + 3*tt**2; h11 = tt**3 - tt**2
    return h00*f[i] + h10*hh*fp[i] + h01*f[i+1] + h11*hh*fp[i+1]

def run(k, ds, y0=1e-4, B=None):
    B = B or ts.background(1e-4, v0=y0)
    y, rho, phi, s, rp = B["y"], B["rho"], B["phi"], B["s"], B["rp"]
    spp = bg.U1(phi) - 4*rp/rho*s; rpp = rho*(B["F"] - s*s/3)
    n = int(math.ceil(math.log(y[-1]/y0)/ds)); sg = np.linspace(math.log(y[-1]), math.log(y0), 2*n + 1)   # half-step grid
    yg = np.exp(sg); yg[0] = y[-1]; yg[-1] = y0
    r = hermite(yg, y, rho, rp); rpg = hermite(yg, y, rp, rpp)
    P = 4*yg*rpg/r - 1; Q = (yg/r)**2
    k = np.asarray(k, float); m2 = 2.25 + k*k
    h = np.ones_like(k); g = np.zeros_like(k); d = sg[2] - sg[0]
    f = lambda i, h, g: (g, -P[i]*g - m2*Q[i]*h)
    for j in range(n):
        i = 2*j
        k1 = f(i, h, g); k2 = f(i+1, h + d/2*k1[0], g + d/2*k1[1]); k3 = f(i+1, h + d/2*k2[0], g + d/2*k2[1]); k4 = f(i+2, h + d*k3[0], g + d*k3[1])
        h, g = h + d/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]), g + d/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    r0, rp0 = r[-1], rpg[-1]
    u = r0**1.5*h; uz = r0**1.5*(1.5*rp0*h + (r0/y0)*g)       # u_z = rho d/dy (rho^{3/2} h), h' = h_s/y
    return u*u + uz*uz/(k*k), B

if __name__ == "__main__":
    B = ts.background(1e-4, v0=1e-4); rho_b = B["rho"][-1]; H = 1/rho_b
    Ik = (ts.simpson(B["rho"]**2, B["y"]) + B["y"][0]**3/3)/rho_b**2; Ip = ts.sol["I_plus"]
    print("I_kept (v0=1e-4 grid) = %.10f" % Ik)
    k = np.array([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 44, 60, 79, 100, 150, 200, 300, 500])
    A2, _ = run(k, 1e-4, B=B); A2b, _ = run(k, 2e-4, B=B)
    w = (2/np.pi)*rho_b**2*Ik/A2; wb = (2/np.pi)*rho_b**2*Ik/A2b
    print("     k        m/H      m L0       w(k)          w/(k H^2)=l_eff^2[L0^2]   w/k^2       step-halving rel.diff")
    for kk, a, b in zip(k, w, wb):
        m = math.sqrt(2.25 + kk*kk); print("  %7.2f  %8.3f  %8.4f  %.6e   %.5f    %.5e   %.1e" % (kk, m, m*H, a, a/(kk*H*H), a/kk**2, abs(b/a - 1)))
    kf = np.exp(np.linspace(math.log(0.02), math.log(600), 500)); wf = (2/np.pi)*rho_b**2*Ik/run(kf, 2e-4, B=B)[0]
    dl = np.diff(np.log(wf))/np.diff(np.log(kf))
    print("fine log grid k in [0.02,600]: w monotone increasing: %s ; d ln w/d ln k range [%.3f, %.3f] (2 at threshold, 1 = RS-like, 0 = flat 5D)" % (bool((np.diff(wf) > 0).all()), dl.min(), dl.max()))
    print("reference: l_-^2 = 3.24 ; 1/(2 k_-^3 I_+) = %.5f ; large-k 5D limit: w -> (2/pi) I_kept H = %.5e" % (1/(2*(5/9)**3*Ip), (2/np.pi)*Ik*H))
    # KK-sum contribution to a static-like brane-brane exchange at separation r >> 1/H ignored; integrated light weight:
    print("integrated weight int_0^K w dk:", [(K, float(np.trapz(wf[kf <= K], kf[kf <= K]))) for K in (1, 10, 79)])
    json.dump(dict(k=k.tolist(), w=w.tolist(), leff2=(w/(k*H*H)).tolist(), kfine=kf.tolist(), wfine=wf.tolist()), open(os.path.join(HERE, "TENSOR_CONTINUUM.json"), "w"))
