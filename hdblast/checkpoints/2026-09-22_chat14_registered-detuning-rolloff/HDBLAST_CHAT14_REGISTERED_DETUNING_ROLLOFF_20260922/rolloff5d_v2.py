#!/usr/bin/env python3
"""HDBLAST Chat 14 - nonlinear 1+1 shell evolution, version 2: STRETCHED GRID so the registered detuning t = 1e-3 can be run.

Same equations, gauge, deviation formulation and junction treatment as Chat 13's rolloff5d.py (imported for the model and the
static background).  New: (1) a smooth coordinate map z = z(xi), uniform xi, with dz/dxi = s(z) = dz_fine near the shell and dz_coarse
far from it, so the wall (width ~1/rho_b in z) is resolved by ~100-300 points while the cone region stays cheap;
(2) the seed deviation is tapered smoothly to zero over the outer 15% of the domain (Chat 13 clamped it, which launched a weak
inward front); (3) optional 4th-order one-sided Neumann closure at the shell (--closure 4) instead of the symmetric (even-extension) one;
(4) the time-derivative fields get the consistent Neumann slopes d/dt(g) at the shell (Chat 13 used 0: an O(h) boundary error, and it made the
    shell-layer constraint monitor meaningless - the momentum constraint is identically satisfied at the shell by the junction data).
Derivatives in z:  f_z = f_xi / z',   f_zz = f_xi xi / z'^2 - f_xi z'' / z'^3   (z' = s(z), z'' = s'(z) s(z), both exact on the nodes).
NOTES AFTER REVIEW (Chat 14 referee): (a) the Kreiss-Oliger term is applied with h = 1 in the index coordinate, so its damping rate is
eps per unit t (about 2,500 times weaker than Chat 13's eps/dz) - effectively no dissipation; all runs were stable and free of grid noise
regardless, but a future version should scale it as eps/dxi_phys.  (b) closure 2 makes the second derivative at the shell first-order
accurate; closure 4 makes the first derivative fourth-order and the second third-order.  (c) the far region behind the receding wall
must be resolved on the +1 side: use --dzc 2e-3 (default 8e-3 under-resolves it after H0 tau ~ 6.5 and pulls H_J down by ~5%)."""
import json, math, os, sys, time
import numpy as np
import rolloff5d_v1 as R

def build_grid(L, dz_fine, dz_coarse, z_fine=0.06, width=0.02):
    """nodes z_i <= 0 from the ODE dz/dxi = s(z), s = dz_coarse + (dz_fine - dz_coarse) * sig(z), integrated from the shell outward"""
    sig = lambda z: 0.5*(1 + np.tanh((z + z_fine)/width))            # 1 near the shell (z > -z_fine), 0 far away
    dsig = lambda z: 0.5*(1 - np.tanh((z + z_fine)/width)**2)/width
    s = lambda z: dz_coarse + (dz_fine - dz_coarse)*sig(z)
    sp = lambda z: (dz_fine - dz_coarse)*dsig(z)
    zs = [0.0]; z = 0.0
    while z > -L:                                                       # RK4 in xi with unit step, outward (decreasing z)
        k1 = -s(z); k2 = -s(z + 0.5*k1); k3 = -s(z + 0.5*k2); k4 = -s(z + k3)
        z = z + (k1 + 2*k2 + 2*k3 + k4)/6; zs.append(z)
    zs = np.array(zs[::-1]); zp = s(zs); zpp = sp(zs)*zp
    return zs, zp, zpp

def run(d=0.0, dc=1e-6, dz_fine=1.5e-4, dz_coarse=8e-3, z_fine=0.06, L=8.0, t_final=12.0, t_det=1e-3, tag="run", out_dir=".", ko_eps=0.05,
        stop_phi=(-1.6, 0.995), closure=2, taper=0.15, log=print, probe_seed=False, width=0.02):
    Ip = 1.0357712571566784; c = 2/Ip - 4/3
    ten = R.Tension(t_det, c, d); ten_seed = R.Tension(t_det, c + dc, d)
    g0 = (math.log10(0.2207*t_det**1.8), 8.2282 + 0.9*math.log(1e-3/t_det))
    ph_h, y_b, r0 = R.solve_shell(ten, guess=g0); ph_h2, y_b2, r2 = R.solve_shell(ten_seed, guess=(math.log10(ph_h + 1), y_b))
    z, zp, zpp = build_grid(L, dz_fine, dz_coarse, z_fine, width); N = len(z) - 1
    S = R.static_on_grid(ph_h, y_b, z, hz=min(dz_fine/4, 2.5e-5)); S2 = R.static_on_grid(ph_h2, y_b2, z, hz=min(dz_fine/4, 2.5e-5))
    rho, Hc, phs, phz = S["rho"], S["Hc"], S["phi"], S["phiz"]; rho2 = rho*rho; Hb, rb, pb = Hc[-1], rho[-1], phs[-1]
    log("grid: %d points, dz from %.2e (shell) to %.2e (far), wall width in z = %.4f -> %.0f points across it" % (N + 1, zp[-1], zp[0], 1/rb, (1/rb)/zp[-1]))
    log("static shell: phi_h+1=%.6e y_b=%.8f rho_b=%.6f phi_b=%+.4e junction residuals=(%.1e,%.1e)" % (ph_h + 1, y_b, rb, pb, r0[0], r0[1]))
    log("junction check on grid: Hc_b - sigma rho_b/6 = %.2e ; phiz_b + sigma' rho_b/2 = %.2e" % (Hb - ten.s(pb)*rb/6, phz[-1] + ten.s1(pb)*rb/2))
    a = S2["lnrho"] - S["lnrho"]; b = a.copy(); f = S2["phi"] - S["phi"]
    if taper > 0:                                                       # smooth taper of the seed over the outer part of the domain
        w = np.ones_like(z); zt = -L*(1 - taper); m = z < zt
        w[m] = 0.5*(1 - np.cos(np.pi*(z[m] + L)/(L*taper))); a *= w; b *= w; f *= w
    pa = np.zeros_like(a); pbv = np.zeros_like(a); pf = np.zeros_like(a)
    log("seed: dc=%.1e  max|a|=%.2e max|f|=%.2e  f_b=%+.3e" % (dc, np.abs(a).max(), np.abs(f).max(), f[-1]))
    h = 1.0                                                             # uniform xi spacing
    def ghosts(F, gR):
        """ghost values at xi = N+1, N+2 for prescribed z-derivative gR at the shell (converted to xi-derivative)"""
        gx = gR*zp[-1]
        if closure == 2:
            return F[-2] + 2*h*gx, F[-3] + 4*h*gx
        e1 = 4*h*gx - 2*F[-3] + 6*F[-2] - (10/3)*F[-1] + F[-4]/3       # 4th-order: 5-point D1 = gx with quartic extrapolation for the 2nd ghost
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
    def neumann(bb, ff):
        eB = rb*math.exp(bb); p = pb + ff
        return (ten.s(p)*eB - ten.s(pb)*rb)/6, -(ten.s1(p)*eB - ten.s1(pb)*rb)/2
    def neumann_t(bb, ff, pbb, pff):
        """time derivatives of the Neumann data (consistent slopes for the time-derivative fields at the shell)"""
        eB = rb*math.exp(bb); p = pb + ff; s2 = 2*R.W2(p) + ten.t*ten.d
        return (ten.s1(p)*pff + ten.s(p)*pbb)*eB/6, -(s2*pff + ten.s1(p)*pbb)*eB/2
    def rhs(state):
        a, b, f, pa, pbv, pf = state
        gA, gF = neumann(b[-1], f[-1])
        az, bz, fz = dz1(a, gA), dz1(b, gA), dz1(f, gF); azz, bzz, fzz = dz2(a, gA), dz2(b, gA), dz2(f, gF)
        dU, dU1, U0, U10 = R.dU_exact(phs, f); e2b = np.expm1(2*b)
        SU = rho2*(e2b*(U0 + dU) + dU); SU1 = rho2*(e2b*(U10 + dU1) + dU1)
        kin = 2*pa + pa*pa; grad = 2*Hc*az + az*az
        ra = azz - 3*kin + 3*grad + (2/3)*SU
        rbb = bzz + 3*kin - 3*grad - 0.5*pf*pf + 0.5*(2*phz*fz + fz*fz) - (1/3)*SU
        rf = fzz - 3*(1 + pa)*pf + 3*((Hc + az)*(phz + fz) - Hc*phz) - SU1
        gAt, gFt = neumann_t(b[-1], f[-1], pbv[-1], pf[-1])
        out = [pa + ko(a, gA), pbv + ko(b, gA), pf + ko(f, gF), ra + ko(pa, gAt), rbb + ko(pbv, gAt), rf + ko(pf, gFt)]
        for o in out: o[0:2] = 0.0
        return out
    def constraint(state):
        a, b, f, pa, pbv, pf = state
        gA, gF = neumann(b[-1], f[-1])
        gAt, gFt = neumann_t(b[-1], f[-1], pbv[-1], pf[-1])
        az, bz, fz = dz1(a, gA), dz1(b, gA), dz1(f, gF); paz = dz1(pa, gAt)
        M = -3*paz - 3*(1 + pa)*(az - bz) + 3*(Hc + az)*pbv - pf*(phz + fz)
        scale = 3*np.abs(paz) + 3*np.abs(Hc + az)*(np.abs(1 + pa) + np.abs(pbv)) + 3*np.abs(1 + pa)*np.abs(Hc + bz) + np.abs(pf)*np.abs(phz + fz) + 1e-300
        r = np.abs(M)/scale; i = int(np.argmax(r[10:-6])) + 10; return float(r[i]), float(np.max(r[-6:])), float(z[i])
    if probe_seed:
        a_, b_, f_, pa_, pbv_, pf_ = a, b, f, pa, pbv, pf
        gA, gF = neumann(b_[-1], f_[-1]); gAt, gFt = neumann_t(b_[-1], f_[-1], pbv_[-1], pf_[-1])
        az, bz, fz = dz1(a_, gA), dz1(b_, gA), dz1(f_, gF); paz = dz1(pa_, gAt)
        M = -3*paz - 3*(1 + pa_)*(az - bz) + 3*(Hc + az)*pbv_ - pf_*(phz + fz)
        scale = 3*np.abs(paz) + 3*np.abs(Hc + az)*(np.abs(1 + pa_) + np.abs(pbv_)) + 3*np.abs(1 + pa_)*np.abs(Hc + bz) + np.abs(pf_)*np.abs(phz + fz) + 1e-300
        r0 = rhs([a_, b_, f_, pa_, pbv_, pf_])
        return dict(z=z, zp=zp, a=a_, f=f_, Mrel=np.abs(M)/scale, Mabs=np.abs(M), rho=rho, phis=phs, acc_a=r0[3], acc_b=r0[4], acc_f=r0[5])
    dt = 0.5*zp.min(); nsteps = int(round(t_final/dt)); state = [a, b, f, pa, pbv, pf]
    rec = []; t = 0.0; tau = 0.0; t0 = time.time(); every = max(1, int(round(0.01/dt))); stop_reason = "t_final"
    for n in range(nsteps + 1):
        if n % every == 0:
            a, b, f, pa, pbv, pf = state
            eB = rb*math.exp(b[-1]); HJ = (1 + pa[-1])/eB; phib = pb + f[-1]
            Mb, Ms, zM = constraint(state) if n % (every*10) == 0 else (float("nan"), float("nan"), float("nan"))
            rec.append((t, phib, HJ*rb, b[-1], pa[-1], pf[-1]/eB*rb, Mb, float(np.max(np.abs(f))), float(np.max(np.abs(b))), tau, float(np.min(phs + f)), float(np.max(phs + f)), Ms))
            if n % (every*50) == 0:
                log("t=%6.2f H0tau=%7.4f phi_b=%+.6e H_J/H0=%+.6f b_b=%+.3e phi_min=%+.4f M_bulk=%.1e@z=%.3f M_shell=%.1e  [%.0fs]" % (t, tau, phib, HJ*rb, b[-1], np.min(phs + f), Mb, zM, Ms, time.time() - t0))
            if not np.isfinite(phib) or not np.isfinite(HJ): stop_reason = "non-finite"; break
            if phib < stop_phi[0] or phib > stop_phi[1]: stop_reason = "phi_b left [%g,%g]" % stop_phi; break
        if n == nsteps: break
        tau += dt*0.5*math.exp(state[1][-1])
        k1 = rhs(state); k2 = rhs([s_ + dt/2*k for s_, k in zip(state, k1)]); k3 = rhs([s_ + dt/2*k for s_, k in zip(state, k2)]); k4 = rhs([s_ + dt*k for s_, k in zip(state, k3)])
        state = [s_ + dt/6*(q1 + 2*q2 + 2*q3 + q4) for s_, q1, q2, q3, q4 in zip(state, k1, k2, k3, k4)]; t += dt
        tau += dt*0.5*math.exp(state[1][-1])
    rec = np.array(rec)
    np.savez_compressed(os.path.join(out_dir, tag + "_timeseries.npz"), rec=rec, z=z, phis=phs, rho=rho, f_end=state[2], b_end=state[1], a_end=state[0])
    tt, dphi = rec[:, 0], np.abs(rec[:, 1] - pb); seed0 = abs(S2["phi"][-1] - S["phi"][-1]) if dc != 0 else 1e-13
    m = (dphi > 30*seed0) & (dphi < 2e-3) & (tt > 1.0); fit = None
    if m.sum() > 20:
        sl = np.polyfit(tt[m], np.log(dphi[m]), 1); fit = dict(growth_rate=float(sl[0]), window=[float(tt[m][0]), float(tt[m][-1])], n=int(m.sum()))
    summ = dict(tag=tag, version=2, d=d, dc=dc, dz_fine=float(zp[-1]), dz_coarse=float(zp[0]), n_points=N + 1, points_across_wall=float((1/rb)/zp[-1]), closure=closure, taper=taper, L=L, t_det=t_det,
                t_end=float(rec[-1, 0]), H0_tau_end=float(rec[-1, 9]), stop_reason=stop_reason, phi_b_static=float(pb), rho_b=float(rb), phi_b_end=float(rec[-1, 1]), HJ_over_H0_end=float(rec[-1, 2]),
                growth_fit=fit, max_bulk_constraint=float(np.nanmax(rec[:, 6])), max_shell_constraint=float(np.nanmax(rec[:, 12])), runtime_s=round(time.time() - t0, 1))
    json.dump(summ, open(os.path.join(out_dir, tag + "_summary.json"), "w"), indent=1); log(json.dumps(summ)); return summ

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=float, default=0.0); ap.add_argument("--dc", type=float, default=1e-6); ap.add_argument("--dzf", type=float, default=1.5e-4); ap.add_argument("--dzc", type=float, default=8e-3)
    ap.add_argument("--zfine", type=float, default=0.06); ap.add_argument("--tdet", type=float, default=1e-3); ap.add_argument("--L", type=float, default=8.0); ap.add_argument("--tf", type=float, default=12.0)
    ap.add_argument("--tag", default="run"); ap.add_argument("--ko", type=float, default=0.05); ap.add_argument("--closure", type=int, default=2); ap.add_argument("--taper", type=float, default=0.15)
    ap.add_argument("--width", type=float, default=0.02); ap.add_argument("--phimin", type=float, default=-1.6); ap.add_argument("--phimax", type=float, default=0.995)
    A_ = ap.parse_args()
    run(d=A_.d, dc=A_.dc, dz_fine=A_.dzf, dz_coarse=A_.dzc, z_fine=A_.zfine, L=A_.L, t_final=A_.tf, t_det=A_.tdet, tag=A_.tag, ko_eps=A_.ko, closure=A_.closure, taper=A_.taper,
        stop_phi=(A_.phimin, A_.phimax), width=A_.width, log=lambda s: print(s, flush=True))
