#!/usr/bin/env python3
"""HDBLAST Chat 9 - transverse-traceless (graviton) sector on the registered cone-to-shell background.

ds^2 = dy^2 + rho(y)^2 (gamma_mn + h_mn),  h_mn TT w.r.t. unit dS_4, (Box_gamma - 2) h_mn = m^2 h_mn  (m^2 in units of H^2):
    h'' + 4 (rho'/rho) h' + (m^2/rho^2) h = 0 ,   h'(y_b) = 0   (Z2 shell, no anisotropic stress).
Schroedinger form: dz = dy/rho, u = rho^{3/2} h:  -u_zz + V1 u = m^2 u,
    V1 = (9/4) rho'^2 + (3/2) rho rho'' = 9/4 + rho^2 [(15/4) F - phi'^2/2],   F = phi'^2/12 - U/6,
    = Q^+ Q with Q = d_z - (3/2) rho';  partner V2 = (9/4) rho'^2 - (3/2) rho rho'' = 9/4 + rho^2 [(9/16) phi'^2 - U/8].
numpy + stdlib only; floating point (RK4), not interval arithmetic.
"""
import json, math, sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "background"))
import hdblast_background as bg

sol = json.load(open(os.path.join(HERE, "..", "background", "REGISTERED_SHELL_FLOAT_SOLUTION.json")))
t, c = sol["t"], sol["c"]

def background(dv, v0=1e-2, phi_h=None, v_b=None):
    v_b = sol["v_b"] if v_b is None else v_b
    n = int(math.ceil((v_b - v0)/dv)); n += n % 2                      # even number of intervals (RK4 on the stored grid uses double steps)
    (rho, phi, s), rec = bg.integrate(sol["phi_h"] if phi_h is None else phi_h, v_b, dv=(v_b - v0)/(n - 0.5), v0=v0, record=True)
    assert (len(rec) - 1) == n
    y, rho, phi, s = rec.T
    F = s*s/12 - bg.U(phi)/6
    rp = np.sqrt(1 + rho*rho*F)
    return dict(y=y, rho=rho, phi=phi, s=s, F=F, rp=rp)

def trap(f, x): return float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(x)))
def simpson(f, x):
    n = len(x) - 1
    if n % 2: return simpson(f[:-1], x[:-1]) + 0.5*(f[-1] + f[-2])*(x[-1] - x[-2])
    hh = (x[-1] - x[0])/n
    return float(hh/3*(f[0] + f[-1] + 4*np.sum(f[1:-1:2]) + 2*np.sum(f[2:-1:2])))

def shoot_from_cone(B, m2):
    """vectorised in m2 (< 9/4): normalisable Frobenius branch h ~ y^a, a = -3/2 + sqrt(9/4 - m2). Returns rho_b h'/h at shell and node count of h."""
    y, rho, rp = B["y"], B["rho"], B["rp"]
    m2 = np.asarray(m2, float); a = -1.5 + np.sqrt(2.25 - m2)
    u0 = bg.U(B["phi"][0]); v0 = y[0]
    c2 = u0*(4*a - m2)/(18*((a + 2)*(a + 5) + m2))
    h = v0**a*(1 + c2*v0**2); hp = v0**(a - 1)*(a + (a + 2)*c2*v0**2)
    # normalise
    hp = hp/h; h = np.ones_like(h); nodes = np.zeros(m2.shape, int)
    P = 4*rp/rho; Qc = 1/rho**2
    def f(i, h, hp): return hp, -P[i]*hp - m2*Qc[i]*h
    n = len(y)
    for i in range(0, n - 2, 2):
        d = y[i+2] - y[i]
        k1 = f(i, h, hp); k2 = f(i+1, h + d/2*k1[0], hp + d/2*k1[1]); k3 = f(i+1, h + d/2*k2[0], hp + d/2*k2[1]); k4 = f(i+2, h + d*k3[0], hp + d*k3[1])
        hn = h + d/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); hp = hp + d/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        nodes += (hn*h < 0); h = hn
        if i % 2000 == 0:
            sc = np.abs(h) + 1e-300; big = sc > 1e50; h = np.where(big, h/sc, h); hp = np.where(big, hp/sc, hp)
    return B["rho"][-1]*hp/h, nodes

def continuum(B, k):
    """m2 = 9/4 + k^2; integrate inward from the shell (h=1,h'=0); amplitude A of u = rho^{3/2} h ~ A cos(k z + delta) at the cone.
    Returns w(k) = (2/pi) rho_b^2 I_kept / A^2 = |u_k(shell)|^2/|u_0(shell)|^2 per unit dk (delta(k-k') normalised continuum)."""
    y, rho, rp = B["y"], B["rho"], B["rp"]
    k = np.asarray(k, float); m2 = 2.25 + k*k
    h = np.ones_like(k); hp = np.zeros_like(k)
    P = 4*rp/rho; Qc = 1/rho**2
    def f(i, h, hp): return hp, -P[i]*hp - m2*Qc[i]*h
    n = len(y); start = n - 1 if (n - 1) % 2 == 0 else n - 2   # need even number of grid intervals
    assert start == n - 1
    for i in range(n - 1, 1, -2):
        d = y[i-2] - y[i]
        k1 = f(i, h, hp); k2 = f(i-1, h + d/2*k1[0], hp + d/2*k1[1]); k3 = f(i-1, h + d/2*k2[0], hp + d/2*k2[1]); k4 = f(i-2, h + d*k3[0], hp + d*k3[1])
        h, hp = h + d/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]), hp + d/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    r0, rp0 = rho[0], rp[0]
    u = r0**1.5*h; uz = r0**1.5*(1.5*rp0*h + r0*hp)
    A2 = u*u + uz*uz/(k*k)
    return A2

if __name__ == "__main__":
    out = {}
    for dv in (2e-4, 1e-4):
        B = background(dv)
        rho_b = B["rho"][-1]; y = B["y"]
        # junction residuals on this grid
        sig = 2*bg.W(B["phi"][-1]) + t*(1 + c*B["phi"][-1]); dsig = 2*bg.W1(B["phi"][-1]) + t*c
        res = (B["rp"][-1]/rho_b - sig/6, B["s"][-1] + dsig/2)
        # I_kept: add the analytic piece from 0 to v0 (rho ~ y): v0^3/3
        Ik = (simpson(B["rho"]**2, y) + y[0]**3/3)/rho_b**2
        print("dv=%g: junction residuals %.2e %.2e ; rho_b=%.8f H=1/rho_b=%.8f ; I_kept=%.10f ; M4^2=2 I_kept=%.10f" % (dv, res[0], res[1], rho_b, 1/rho_b, Ik, 2*Ik))
    Ip = sol["I_plus"]; H = 1/rho_b
    out.update(I_kept=Ik, M4sq=2*Ik, two_I_plus=2*Ip, ratio_Ip_over_Ikept=Ip/Ik, H=H, rho_b=rho_b, eightpiG_zero_mode=1/(2*Ik))
    print("2 I_plus = %.10f ; I_plus/I_kept = F_shell^2 = %.8f ; 8piG_N = 1/(2 I_kept) = %.6f (Theorem 4 quotes 0.48362)" % (2*Ip, Ip/Ik, 1/(2*Ik)))
    # sum rule Theorem 4
    Z = B["rp"]/B["rho"] - bg.W(B["phi"])/3; Fd = B["s"] + bg.W1(B["phi"])
    quad = simpson((B["rho"]/rho_b)**4*(2*Z*Z - Fd*Fd/6), y)
    dlam = t*(1 + c*B["phi"][-1])
    print("sum rule: dlam/6 = %.10e ; h I_kept = %.10e ; quadratic = %.4e ; residual = %.2e" % (dlam/6, Ik/rho_b**2, quad, dlam/6 - Ik/rho_b**2 - quad))
    out.update(sumrule=dict(lhs=dlam/6, hI=Ik/rho_b**2, quad=quad, residual=dlam/6 - Ik/rho_b**2 - quad))
    # potentials
    V1 = 2.25 + B["rho"]**2*(3.75*B["F"] - B["s"]**2/2)
    V2 = 2.25 + B["rho"]**2*((9/16)*B["s"]**2 - bg.U(B["phi"])/8)
    print("V1: value at cone end %.8f, min %.6f, max (at shell) %.4f ; delta-well strength 3 rho'_b = %.5f" % (V1[0], V1.min(), V1[-1], 3*B["rp"][-1]))
    print("partner V2 - 9/4: min over background = %.6e (at y=%.4f) ; value at shell %.4f   => V2 >= 9/4 everywhere: %s" % ((V2 - 2.25).min(), y[np.argmin(V2)], V2[-1] - 2.25, bool((V2 >= 2.25).all())))
    out.update(V1_min=float(V1.min()), V1_shell=float(V1[-1]), V2minus94_min=float((V2 - 2.25).min()), delta_strength=float(3*B["rp"][-1]))
    # bound-state search
    m2 = np.concatenate([np.linspace(-60, -0.5, 120), np.linspace(-0.495, 2.2499, 550)])
    D, nodes = shoot_from_cone(B, m2)
    sgn = np.sign(D); zc = np.where(sgn[:-1]*sgn[1:] < 0)[0]
    print("shooting D(m2) = rho_b h'/h at shell for m2 in [-60, 2.2499], %d samples: sign changes at" % len(m2), [(float(m2[i]), float(m2[i+1])) for i in zc], "; max nodes of h:", int(nodes.max()))
    i0 = np.argmin(np.abs(m2)); print("   D at m2 = %.4f : %.3e ;  D(-0.5)=%.4f  D(0.5)=%.4f  D(1)=%.4f D(2)=%.4f D(2.2499)=%.4f" % (m2[i0], D[i0], *[float(D[np.argmin(np.abs(m2 - v))]) for v in (-0.5, 0.5, 1, 2, 2.2499)]))
    Dz, _ = shoot_from_cone(B, np.array([0.0, 1e-6, -1e-6]))
    print("   D(0) = %.3e, dD/dm2 at 0 = %.6f  (exact: D = -m2 * int rho^2 h dy/(rho_b^3 h_b) => -I_kept/rho_b = %.6f)" % (Dz[0], (Dz[1] - Dz[2])/2e-6, -Ik/rho_b))
    out.update(bound_state_scan=dict(m2_range=[-60, 2.2499], sign_changes=[(float(m2[i]), float(m2[i+1])) for i in zc], monotone=bool((np.diff(D) < 0).all())))
    print("   D monotone decreasing on the scan:", bool((np.diff(D) < 0).all()))
    # continuum weights
    k = np.array([0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30, 44, 60, 80, 100, 150, 200])
    A2 = continuum(B, k); w = (2/np.pi)*rho_b**2*Ik/A2
    kf = np.linspace(0.05, 120, 800); wf = (2/np.pi)*rho_b**2*Ik/continuum(B, kf)
    print("\ncontinuum brane weight w(k) = |u_k(b)|^2/|u_0(b)|^2 per dk  (m = H sqrt(9/4+k^2)); RS-like expectation w ~ l_eff^2 H^2 k")
    print("     k        m/H        w(k)          w/(k H^2) [= l_eff^2 in L0^2]    w/k^2 (small k)")
    for kk, ww in zip(k, w): print("  %7.2f  %8.3f   %.6e   %.5f    %.5e" % (kk, math.sqrt(2.25 + kk*kk), ww, ww/(kk*H*H), ww/kk**2))
    print("reference: l_-^2 = (9/5)^2 = 3.24 ; thick-wall estimate 1/(2 k_-^3 I_+) = %.5f ; wall scale k ~ 1/(H L0) = %.1f" % (1/(2*(5/9)**3*Ip), rho_b))
    dl = np.diff(np.log(wf))/np.diff(np.log(kf)); print("max |d ln w/d ln k| on fine grid (resonance check) = %.3f ; w monotone increasing: %s" % (np.abs(dl).max(), bool((np.diff(wf) > 0).all())))
    out.update(continuum=dict(k=k.tolist(), w=w.tolist(), leff2=(w/(k*H*H)).tolist(), monotone=bool((np.diff(wf) > 0).all())))
    json.dump(out, open(os.path.join(HERE, "TENSOR_SPECTRUM.json"), "w"), indent=1)
