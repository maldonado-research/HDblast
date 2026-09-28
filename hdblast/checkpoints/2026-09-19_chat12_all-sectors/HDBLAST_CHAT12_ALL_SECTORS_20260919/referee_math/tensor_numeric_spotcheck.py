#!/usr/bin/env python3
"""Referee (math) - FLOATING-POINT spot check (numpy only; NOT a certificate) of the tensor sector on
  (A) S_8/5  (t=1e-3, c=0.597594935028, d=8/5; certified phi_h, y_b from Chat 10), and
  (B) the Chat 9 registered shell (d=0, c=c_star, solved with the Chat 9 float Newton solver).
 1. background in the conformal coordinate z (d/dz = rho d/dy), junction residuals, range of phi, sign of U, min(V2-9/4), min(V1-9/4), B.
 2. Riccati shooting for h''+4Hh'+m2 h/rho^2=0:  R = rho h'/h,  R_z = -3 rho' R - m2 - R^2,  R(-inf) = s_+ = (-3+sqrt(9-4 m2))/2 (normalisable branch);
    a Neumann bound state needs R(z_b) = 0.  Scan m2 in (-5, 9/4).
 3. finite-volume spectrum of -(rho^4 h')' = m2 rho^2 h on a geometric y-grid [y_min, y_b], Dirichlet at y_min, Neumann h'=0 at y_b.\n    (A first attempt on a uniform z-grid with dz=3e-3 was under-resolved: the wall next to the shell has z-thickness ~1/rho_b ~ 0.013.)"""
import json, math, sys, os
import numpy as np
sys.path.insert(0, "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/background")
import hdblast_background as hb
W, W1, W2, U, U1 = hb.W, hb.W1, hb.W2, hb.U, hb.U1

def background_z(phi_h, y_b, y0=1e-4, n=24000):
    """integrate (y, rho, phi, s) in z from y0 to y_b; the last step is adjusted to land on y_b."""
    st = np.array([y0, *hb.series(phi_h, y0)])
    def f(v):
        yy, rho, phi, s = v
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
        return np.array([rho, rho*rp, rho*s, rho*(U1(phi) - 4*(rp/rho)*s)])
    # estimate total z-length by a coarse pass
    def run(dz, record):
        v = st.copy(); z = math.log(y0); out = [(z, *v)]
        while v[0] < y_b:
            k1 = f(v)
            def rk(step):
                k2 = f(v + step/2*k1); k3 = f(v + step/2*k2); k4 = f(v + step*k3)
                return v + step/6*(k1 + 2*k2 + 2*k3 + k4)
            vn = rk(dz)
            if vn[0] < y_b:
                v = vn; z += dz
                if record: out.append((z, *v))
                continue
            step = (y_b - v[0])/k1[0]
            for _ in range(60):                                   # Newton on the last step length so that y lands on y_b
                vn = rk(step); err = vn[0] - y_b
                if abs(err) < 1e-13: break
                step -= err/vn[1]
            v = vn; z += step
            if record: out.append((z, *v))
            break
        return z, np.array(out) if record else None
    zb, _ = run(2e-3, False)
    dz = (zb - math.log(y0))/n
    zb, arr = run(dz, True)
    return arr

def analyse(tag, phi_h, y_b, t, c, d):
    arr = background_z(phi_h, y_b)
    z, yy, rho, phi, s = arr.T
    rp = np.sqrt(1 + rho**2*(s**2/12 - U(phi)/6))
    sig = 2*W(phi[-1]) + t*(1 + c*phi[-1] + d*phi[-1]**2/2); dsig = 2*W1(phi[-1]) + t*(c + d*phi[-1]); d2sig = 2*W2(phi[-1]) + t*d
    res_j = [rp[-1]/rho[-1] - sig/6, s[-1] + dsig/2]
    phipp = U1(phi[-1]) - 4*rp[-1]/rho[-1]*s[-1]
    B = phipp/s[-1] + d2sig/2
    V2m = rho**2*(9*s**2/16 - U(phi)/8); V1m = rho**2*(-5*U(phi)/8 - 3*s**2/16)
    out = dict(y_end=float(yy[-1]), rho_b=float(rho[-1]), phi_b=float(phi[-1]), junction_residuals=[float(x) for x in res_j], B=float(B),
               phi_min=float(phi.min()), phi_max=float(phi.max()), phi_nondecreasing=bool(np.all(np.diff(phi) >= 0)), phi_strictly_increasing_for_y_gt_0p05=bool(np.all(np.diff(phi[yy > 0.05]) > 0)), Umax_along_bulk=float(U(phi).max()),
               min_V2_minus_9over4_over_rho2=float((V2m/rho**2).min()), min_V1_minus_9over4_over_rho2=float((V1m/rho**2).min()),
               rho_prime_min=float(rp.min()))
    # Riccati shooting on the stored grid (RK4 with midpoint values by linear interpolation of rho')
    def shoot(m2):
        R = (-3 + math.sqrt(9 - 4*m2))/2; nodes = 0
        for i in range(len(z) - 1):
            dz = z[i+1] - z[i]; a0, a1 = rp[i], rp[i+1]; am = 0.5*(a0 + a1)
            g = lambda R_, a_: -3*a_*R_ - m2 - R_*R_
            k1 = g(R, a0); k2 = g(R + dz/2*k1, am); k3 = g(R + dz/2*k2, am); k4 = g(R + dz*k3, a1)
            R += dz/6*(k1 + 2*k2 + 2*k3 + k4)
            if not np.isfinite(R) or R < -1e6: nodes += 1; return None, nodes
        return R, nodes
    scan = {}
    grid = [-5.0, -2.0, -0.5, -1e-3, 1e-6, 1e-3, 0.1, 0.5, 1.0, 1.5, 2.0, 2.2, 2.24, 2.2499]
    for m2 in grid:
        R, nodes = shoot(m2); scan["%g" % m2] = None if R is None else float(R)
    out["riccati_R_at_shell_vs_m2"] = scan
    vals = [v for k, v in scan.items() if v is not None]
    out["riccati_no_zero_for_m2_in_(0,9/4)"] = bool(all(scan["%g" % m2] is not None and scan["%g" % m2] < 0 for m2 in grid if m2 > 0))
    # FD spectrum: finite-volume discretisation of -(rho^4 h')' = m2 rho^2 h on a geometric y-grid (resolves the cone logarithmically and the wall
    # near the shell with dy ~ 0.01-0.03); natural (Neumann) condition h' = 0 at y_b, Dirichlet h = 0 at y_min.
    spec = {}
    for ymin, N in ((1e-4, 3000), (1e-3, 3000), (1e-4, 1500)):
        yg = ymin*(yy[-1]/ymin)**(np.arange(N + 1)/N); rg = np.interp(yg, yy, rho)
        ym = 0.5*(yg[1:] + yg[:-1]); pm = np.interp(ym, yy, rho)**4; dl = np.diff(yg)
        Mw = np.zeros(N + 1); Mw[1:-1] = rg[1:-1]**2*0.5*(dl[1:] + dl[:-1]); Mw[-1] = rg[-1]**2*0.5*dl[-1]
        n_ = N                                                    # unknowns h_1..h_N (h_0 = 0)
        K = np.zeros((n_, n_))
        for i in range(1, N + 1):
            K[i-1, i-1] = pm[i-1]/dl[i-1] + (pm[i]/dl[i] if i < N else 0.0)
            if i < N: K[i-1, i] = K[i, i-1] = -pm[i]/dl[i]
        sc = 1/np.sqrt(Mw[1:]); S = K*sc[:, None]*sc[None, :]
        ev = np.linalg.eigvalsh(S)[:4]
        spec["y_min=%g,N=%d" % (ymin, N)] = dict(lowest_eigenvalues=[float(e) for e in ev], dy_at_shell=float(dl[-1]))
    out["fd_spectrum"] = spec
    print(tag, json.dumps(out, indent=1)); return out

if __name__ == "__main__":
    R = {}
    R["S_8/5"] = analyse("S_8/5", -1 + 8.78437001888956254016e-7, 8.22819172487514216035, 1e-3, 5975949350280/1e13, 1.6)
    Ip = hb.I_plus(); c_star = 2/Ip - 4/3
    sol = hb.solve_shell(1e-3, c_star, guess=(math.log10(8.7855e-7), 8.6))
    R["registered_Chat9"] = analyse("registered", sol["phi_h"], sol["v_b"], 1e-3, c_star, 0.0)
    R["registered_Chat9"]["solver"] = {k: sol[k] for k in ("phi_h", "v_b", "phi_b", "residual")}
    json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "TENSOR_NUMERIC_SPOTCHECK.json"), "w"), indent=1)
