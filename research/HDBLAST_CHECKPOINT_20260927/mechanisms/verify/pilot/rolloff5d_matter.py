#!/usr/bin/env python3
"""PILOT: Chat14 v2 shell roll-off (copied numerics) + a shell radiation fluid fed by a dissipative scalar junction.

Derived from rolloff5d_v2.py (Chat14, sha256 7dc0cd58...) with ONLY these changes (kappa_5 = 1, model units):
  * shell radiation density R = kappa_5^2 rho, p = R/3, and scalar source kappa_5^2 j = Y v, v = dphi_b/dtau
    (a friction-type coupling: the matter action is not specified; this is a phenomenological closure that
    respects the exact ledger  dR/dtau + 4 H R = j v = Y v^2  and  d(sigma+R)/dtau + 4HR = (sigma' + Y v) v);
  * the three junctions of the 22 Sept matter extension:
        nA = (sigma + R)/6,   nB = (sigma - 3R)/6,   n.phi = -(sigma' + Y v)/2 ;
    A and B now get different Neumann data (Chat14 used one gA for both);
  * dR/dt = -4 (1 + a_t) R + Y phi_t^2 e^{-B}  (t = conformal time, e^B = rho_b e^b at the shell);
  * neumann_t (used only by the weak Kreiss-Oliger boundary ghosts) omits the Y d(phi_t)/dt and dR/dt pieces
    (documented approximation; KO strength is 'effectively off' per the Chat14 review).
Records extra columns: R, Weyl scalar Wy/H0^2 (H^2 identity with matter), energy-ledger residual.
Y = 0 must reproduce the vacuum run exactly (control).  Floating point; a pilot, not a converged result."""
import json, math, os, sys, time
import numpy as np
import rolloff5d_v1 as R_

def build_grid(L, dz_fine, dz_coarse, z_fine=0.06, width=0.02):
    sig = lambda z: 0.5*(1 + np.tanh((z + z_fine)/width))
    dsig = lambda z: 0.5*(1 - np.tanh((z + z_fine)/width)**2)/width
    s = lambda z: dz_coarse + (dz_fine - dz_coarse)*sig(z)
    sp = lambda z: (dz_fine - dz_coarse)*dsig(z)
    zs = [0.0]; z = 0.0
    while z > -L:
        k1 = -s(z); k2 = -s(z + 0.5*k1); k3 = -s(z + 0.5*k2); k4 = -s(z + k3)
        z = z + (k1 + 2*k2 + 2*k3 + k4)/6; zs.append(z)
    zs = np.array(zs[::-1]); zp = s(zs); zpp = sp(zs)*zp
    return zs, zp, zpp

def run(Y=0.0, dc=1e-4, dz_fine=1e-3, dz_coarse=8e-3, z_fine=0.4, L=12.0, t_final=12.0, t_det=0.1, tag="run", out_dir=".",
        ko_eps=0.05, stop_phi=(-1.6, 0.995), closure=2, taper=0.15, width=0.02, log=print, c_override=None, d_quad=0.0):
    Ip = 1.0357712571566784; c = 2/Ip - 4/3 if c_override is None else c_override
    ten = R_.Tension(t_det, c, d_quad); ten_seed = R_.Tension(t_det, c + dc, d_quad)
    g0 = (math.log10(0.2207*t_det**1.8), 8.2282 + 0.9*math.log(1e-3/t_det))
    ph_h, y_b, r0 = R_.solve_shell(ten, guess=g0); ph_h2, y_b2, r2 = R_.solve_shell(ten_seed, guess=(math.log10(ph_h + 1), y_b))
    z, zp, zpp = build_grid(L, dz_fine, dz_coarse, z_fine, width); N = len(z) - 1
    S = R_.static_on_grid(ph_h, y_b, z, hz=min(dz_fine/4, 2.5e-5)); S2 = R_.static_on_grid(ph_h2, y_b2, z, hz=min(dz_fine/4, 2.5e-5))
    rho, Hc, phs, phz = S["rho"], S["Hc"], S["phi"], S["phiz"]; rho2 = rho*rho; Hb, rb, pb = Hc[-1], rho[-1], phs[-1]
    a = S2["lnrho"] - S["lnrho"]; b = a.copy(); f = S2["phi"] - S["phi"]
    if taper > 0:
        w = np.ones_like(z); zt = -L*(1 - taper); m = z < zt
        w[m] = 0.5*(1 - np.cos(np.pi*(z[m] + L)/(L*taper))); a *= w; b *= w; f *= w
    pa = np.zeros_like(a); pbv = np.zeros_like(a); pf = np.zeros_like(a); Rm = np.zeros(1)
    h = 1.0
    def ghosts(F, gR):
        gx = gR*zp[-1]
        if closure == 2:
            return F[-2] + 2*h*gx, F[-3] + 4*h*gx
        e1 = 4*h*gx - 2*F[-3] + 6*F[-2] - (10/3)*F[-1] + F[-4]/3
        e2 = 5*e1 - 10*F[-1] + 10*F[-2] - 5*F[-3] + F[-4]
        return e1, e2
    def ext(F, gR):
        e1, e2 = ghosts(F, gR); return np.concatenate(([0.0, 0.0], F, [e1, e2]))
    def D1x(E): return (E[:-4] - 8*E[1:-3] + 8*E[3:-1] - E[4:])/(12*h)
    def D2x(E): return (-E[:-4] + 16*E[1:-3] - 30*E[2:-2] + 16*E[3:-1] - E[4:])/(12*h*h)
    def dz1(F, gR):
        E = ext(F, gR); return D1x(E)/zp
    def dz2(F, gR):
        E = ext(F, gR); fx = D1x(E); fxx = D2x(E); return fxx/zp**2 - fx*zpp/zp**3
    def ko(F, gR):
        E = ext(F, gR); E = np.concatenate(([0.0], E, [2*E[-1] - E[-2]]))
        return ko_eps*(E[:-6] - 6*E[1:-5] + 15*E[2:-4] - 20*E[3:-3] + 15*E[4:-2] - 6*E[5:-1] + E[6:])/(64*h)
    def neumann(bb, ff, pff, Rv):
        eB = rb*math.exp(bb); p = pb + ff
        gA = ((ten.s(p) + Rv)*eB - ten.s(pb)*rb)/6
        gB = ((ten.s(p) - 3*Rv)*eB - ten.s(pb)*rb)/6
        gF = -((ten.s1(p))*eB + Y*pff - ten.s1(pb)*rb)/2   # s1 includes the quadratic term t*d*p
        return gA, gB, gF
    def neumann_t(bb, ff, pbb, pff):
        eB = rb*math.exp(bb); p = pb + ff; s2 = 2*R_.W2(p) + ten.t*ten.d
        return (ten.s1(p)*pff + ten.s(p)*pbb)*eB/6, -(s2*pff + ten.s1(p)*pbb)*eB/2
    def rhs(state):
        a, b, f, pa, pbv, pf, Rm = state
        gA, gB, gF = neumann(b[-1], f[-1], pf[-1], Rm[0])
        az, bz, fz = dz1(a, gA), dz1(b, gB), dz1(f, gF); azz, bzz, fzz = dz2(a, gA), dz2(b, gB), dz2(f, gF)
        dU, dU1, U0, U10 = R_.dU_exact(phs, f); e2b = np.expm1(2*b)
        SU = rho2*(e2b*(U0 + dU) + dU); SU1 = rho2*(e2b*(U10 + dU1) + dU1)
        kin = 2*pa + pa*pa; grad = 2*Hc*az + az*az
        ra = azz - 3*kin + 3*grad + (2/3)*SU
        rbb = bzz + 3*kin - 3*grad - 0.5*pf*pf + 0.5*(2*phz*fz + fz*fz) - (1/3)*SU
        rf = fzz - 3*(1 + pa)*pf + 3*((Hc + az)*(phz + fz) - Hc*phz) - SU1
        gAt, gFt = neumann_t(b[-1], f[-1], pbv[-1], pf[-1])
        eB = rb*math.exp(b[-1])
        dR = np.array([-4*(1 + pa[-1])*Rm[0] + Y*pf[-1]**2/eB])
        out = [pa + ko(a, gA), pbv + ko(b, gB), pf + ko(f, gF), ra + ko(pa, gAt), rbb + ko(pbv, gAt), rf + ko(pf, gFt)]
        for o in out: o[0:2] = 0.0
        return out + [dR]
    dt = 0.5*zp.min(); nsteps = int(round(t_final/dt)); state = [a, b, f, pa, pbv, pf, Rm]
    rec = []; t = 0.0; tau = 0.0; t0 = time.time(); every = max(1, int(round(0.01/dt))); stop_reason = "t_final"
    Wm = 0.0     # cumulative matter work  int Y v^2 dtau  (in H0 units at record)
    for n in range(nsteps + 1):
        if n % every == 0:
            a_, b_, f_, pa_, pbv_, pf_, R_v = state
            eB = rb*math.exp(b_[-1]); HJ = (1 + pa_[-1])/eB; phib = pb + f_[-1]; v = pf_[-1]/eB; Rv = R_v[0]
            sig = ten.s(phib); s1 = ten.s1(phib)
            Wy = HJ**2 - (sig + Rv)**2/36 - v*v/12 + (s1 + Y*v)**2/48 - R_.U(phib)/6
            rec.append((t, phib, HJ*rb, b_[-1], pa_[-1], v*rb, Rv, Wy*rb*rb, tau, float(np.min(phs + f_)), float(np.max(phs + f_))))
            if n % (every*100) == 0:
                log("t=%6.2f H0tau=%7.4f phi_b=%+.5f H/H0=%+.5f R=%.4e Wy/H0^2=%+.4e [%.0fs]" % (t, tau, phib, HJ*rb, Rv, Wy*rb*rb, time.time() - t0))
            if not np.isfinite(phib) or not np.isfinite(HJ): stop_reason = "non-finite"; break
            if phib < stop_phi[0] or phib > stop_phi[1]: stop_reason = "phi_b left [%g,%g]" % stop_phi; break
        if n == nsteps: break
        tau += dt*0.5*math.exp(state[1][-1])
        k1 = rhs(state); k2 = rhs([s_ + dt/2*k for s_, k in zip(state, k1)]); k3 = rhs([s_ + dt/2*k for s_, k in zip(state, k2)]); k4 = rhs([s_ + dt*k for s_, k in zip(state, k3)])
        state = [s_ + dt/6*(q1 + 2*q2 + 2*q3 + q4) for s_, q1, q2, q3, q4 in zip(state, k1, k2, k3, k4)]; t += dt
        tau += dt*0.5*math.exp(state[1][-1])
    rec = np.array(rec)
    cols = ['t', 'phi_b', 'H_over_H0', 'b_b', 'a_t_b', 'v_over_H0', 'R', 'Wy_over_H0sq', 'H0_tau', 'phi_min', 'phi_max']
    np.savez_compressed(os.path.join(out_dir, tag + "_timeseries.npz"), rec=rec, cols=np.array(cols))
    # radiation ledger check along the record: dR/dtau + 4 H R - Y v^2  (tau, H in model units)
    s_ = rec[:, 8]; ok = np.concatenate(([True], np.diff(s_) > 1e-12)); rr = rec[ok]
    vv = rr[:, 5]/rb
    dRdtau = np.gradient(rr[:, 6], rr[:, 8]*rb); Hh = rr[:, 2]/rb
    ledger = dRdtau + 4*Hh*rr[:, 6] - Y*vv*vv
    ledger_scale = np.abs(dRdtau) + np.abs(4*Hh*rr[:, 6]) + np.abs(Y*vv*vv) + 1e-300
    summ = dict(tag=tag, Y=Y, c=c, d_quad=d_quad, z_fine=z_fine, dc=dc, t_det=t_det, dz_fine=float(zp[-1]), dz_coarse=float(zp[0]), n_points=N + 1, closure=closure, L=L,
                rho_b=float(rb), phi_b_static=float(pb), t_end=float(rec[-1, 0]), H0_tau_end=float(rec[-1, 8]), stop_reason=stop_reason,
                phi_b_end=float(rec[-1, 1]), H_over_H0_end=float(rec[-1, 2]), R_end=float(rec[-1, 6]), R_max=float(rec[:, 6].max()),
                Wy_over_H0sq_end=float(rec[-1, 7]),
                radiation_ledger_median_rel=float(np.median(np.abs(ledger[3:-3])/ledger_scale[3:-3])) if Y > 0 else 0.0,
                runtime_s=round(time.time() - t0, 1))
    json.dump(summ, open(os.path.join(out_dir, tag + "_summary.json"), "w"), indent=1); log(json.dumps(summ)); return summ

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--Y", type=float, default=0.0); ap.add_argument("--dc", type=float, default=1e-4)
    ap.add_argument("--dzf", type=float, default=1e-3); ap.add_argument("--dzc", type=float, default=8e-3)
    ap.add_argument("--tdet", type=float, default=0.1); ap.add_argument("--L", type=float, default=12.0); ap.add_argument("--tf", type=float, default=12.0)
    ap.add_argument("--closure", type=int, default=2); ap.add_argument("--tag", default="run"); ap.add_argument("--phimax", type=float, default=0.995); ap.add_argument("--zfine", type=float, default=0.4); ap.add_argument("--c", type=float, default=None); ap.add_argument("--d", type=float, default=0.0)
    A_ = ap.parse_args()
    run(Y=A_.Y, dc=A_.dc, dz_fine=A_.dzf, dz_coarse=A_.dzc, z_fine=A_.zfine, L=A_.L, t_final=A_.tf, t_det=A_.tdet, tag=A_.tag, closure=A_.closure,
        stop_phi=(-1.6, A_.phimax), log=lambda s: print(s, flush=True), c_override=A_.c, d_quad=A_.d)
