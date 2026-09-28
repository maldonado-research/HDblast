#!/usr/bin/env python3
"""HDBLAST Chat 13 - fully NON-LINEAR 1+1 evolution of the Z2 shell in the 5D Einstein-scalar model (numpy only).

Chart (2D conformal gauge, flat 3-slices; the shell is kept at z = 0, bulk on z < 0):
    ds^2 = e^{2B(t,z)} (-dt^2 + dz^2) + e^{2A(t,z)} dx_3^2 ,    phi(t,z)
    A_tt = A_zz - 3A_t^2 + 3A_z^2 + (2/3) e^{2B} U
    B_tt = B_zz + 3A_t^2 - 3A_z^2 - phi_t^2/2 + phi_z^2/2 - (1/3) e^{2B} U
    phi_tt = phi_zz - 3A_t phi_t + 3A_z phi_z - e^{2B} U'
    shell (Israel + scalar junction, pure tension):  A_z = B_z = (sigma_t(phi)/6) e^{B},   phi_z = -(sigma_t'(phi)/2) e^{B}   at z = 0
    constraints (monitored, not imposed):  M = -3A_tz - 3A_tA_z + 3A_tB_z + 3A_zB_t - phi_t phi_z,
                                           Hc = -2Ue^{2B} + 6A_t^2 + 6A_tB_t - 12A_z^2 + 6A_zB_z - 6A_zz - phi_t^2 - phi_z^2
(all derived and checked symbolically in derive_evolution_equations.py).
Static registered solution:  B = ln rho(z), A = ln rho(z) + t, phi = phi0(z), dz = dy/rho.  One unit of t is one Hubble time of the static shell.

Numerical method.  The code evolves the DEVIATION (a, b, f) = (A - A_static, B - B_static, phi - phi_static) with the static derivatives
supplied analytically from the background ODEs, so the static solution is static to round-off and truncation error acts only on the deviation.
The system is nevertheless fully non-linear.  4th-order centred differences, RK4, Kreiss-Oliger dissipation, Neumann data at the shell through
ghost points (even extension plus the linear odd part), Dirichlet (deviation = 0) at the far end z = -L, which is causally isolated for t < L.
Seed: the exact static solution for a slightly different tension slope c' = c + dc (bulk constraints hold exactly; the junction is off by O(dc)).
Floating point; not a certificate."""
import json, math, os, sys, time
import numpy as np

# ---------------- model
def W(p):  return 1 - p + p**3/3
def W1(p): return p*p - 1
def W2(p): return 2*p
def U(p):  return 0.5*W1(p)**2 - (2/3)*W(p)**2
_Uc = np.array([-1/6, 4/3, -5/3, -4/9, 17/18, 0.0, -2/27])                 # U = sum c_k phi^k
def _poly_derivs(p, order):
    """derivatives U^(j)(p), j = 0..order, of the registered sextic (exact polynomial arithmetic)"""
    out = []; c = _Uc.copy()
    for j in range(order + 1):
        out.append(sum(ck*p**k for k, ck in enumerate(c)))
        c = np.array([k*c[k] for k in range(1, len(c))]) if len(c) > 1 else np.array([0.0])
    return out
def dU_exact(ps, f):
    """U(ps+f) - U(ps) and U'(ps+f) - U'(ps) without cancellation (finite Taylor sums of the sextic)"""
    D = _poly_derivs(ps, 6)
    dU = sum(D[j]*f**j/math.factorial(j) for j in range(1, 7))
    dU1 = sum(D[j+1]*f**j/math.factorial(j) for j in range(1, 6))
    return dU, dU1, D[0], D[1]

class Tension:
    def __init__(self, t, c, d): self.t, self.c, self.d = t, c, d
    def s(self, p):  return 2*W(p) + self.t*(1 + self.c*p + self.d*p*p/2)
    def s1(self, p): return 2*W1(p) + self.t*(self.c + self.d*p)

# ---------------- static background in the conformal coordinate
def cone_series(phi_h, y):
    u, u1 = U(phi_h), (W1(phi_h)*W2(phi_h) - (4/3)*W(phi_h)*W1(phi_h))
    u2 = W2(phi_h)**2 + 2*W1(phi_h) - (4/3)*(W1(phi_h)**2 + W(phi_h)*W2(phi_h))
    k4 = u1*(u2/280 + u/630)
    rho = y - u*y**3/36; phi = phi_h + u1*y*y/10 + k4*y**4; s = u1*y/5 + 4*k4*y**3
    return rho, phi, s
def U1f(p): return W1(p)*W2(p) - (4/3)*W(p)*W1(p)
def bg_rhs_y(Y):
    rho, phi, s = Y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    return np.array([rp, s, U1f(phi) - 4*(rp/rho)*s])
def integrate_y(phi_h, y_end, y0=1e-2, dv=2e-4):
    Y = np.array(cone_series(phi_h, y0)); n = max(1, int(math.ceil((y_end - y0)/dv))); h = (y_end - y0)/n
    for _ in range(n):
        k1 = bg_rhs_y(Y); k2 = bg_rhs_y(Y + h/2*k1); k3 = bg_rhs_y(Y + h/2*k2); k4 = bg_rhs_y(Y + h*k3)
        Y = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return Y
def solve_shell(ten, guess=(math.log10(8.7855e-7), 8.2282)):
    def res(x):
        rho, phi, s = integrate_y(-1 + 10**x[0], x[1])
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
        return np.array([rp/rho - ten.s(phi)/6, s + ten.s1(phi)/2])
    x = np.array(guess, float)
    for _ in range(60):
        F = res(x); J = np.zeros((2, 2))
        for j in range(2):
            xp = x.copy(); xp[j] += 1e-7; xm = x.copy(); xm[j] -= 1e-7; J[:, j] = (res(xp) - res(xm))/2e-7
        dx = np.linalg.solve(J, -F); lam = 1.0
        while lam > 1e-4 and np.linalg.norm(res(x + lam*dx)) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < 1e-13: break
    return -1 + 10**x[0], x[1], res(x)

def static_on_grid(phi_h, y_b, zgrid, hz=None):
    """background (rho, rho_y, phi, phi_y) at the conformal grid points z <= 0 (shell at z = 0).  Integrates the ODE in z outward from the cone."""
    dzg = abs(zgrid[1] - zgrid[0]); hz = hz or dzg/4
    y0 = 1e-3
    def rhs(Yv):
        y, rho, phi, s = Yv
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
        return np.array([rho, rho*rp, rho*s, rho*(U1f(phi) - 4*(rp/rho)*s)])
    r0, p0, s0 = cone_series(phi_h, y0); Yv = np.array([y0, r0, p0, s0]); zr = 0.0
    tab = [(zr, *Yv)]
    while Yv[0] < y_b:
        h = min(hz, 2e-3/Yv[1])                          # keep the proper-distance step below 2e-3 near the shell (rho ~ 79)
        k1 = rhs(Yv); k2 = rhs(Yv + h/2*k1); k3 = rhs(Yv + h/2*k2); k4 = rhs(Yv + h*k3)
        Ynew = Yv + h/6*(k1 + 2*k2 + 2*k3 + k4)
        if Ynew[0] > y_b:                                   # land exactly on the shell: shrink the last step (secant on y)
            lo, hi = 0.0, h
            for _ in range(60):
                mid = 0.5*(lo + hi)
                k1 = rhs(Yv); k2 = rhs(Yv + mid/2*k1); k3 = rhs(Yv + mid/2*k2); k4 = rhs(Yv + mid*k3)
                Ym = Yv + mid/6*(k1 + 2*k2 + 2*k3 + k4)
                if Ym[0] > y_b: hi = mid
                else: lo = mid
            h = 0.5*(lo + hi); k1 = rhs(Yv); k2 = rhs(Yv + h/2*k1); k3 = rhs(Yv + h/2*k2); k4 = rhs(Yv + h*k3)
            Ynew = Yv + h/6*(k1 + 2*k2 + 2*k3 + k4); zr += h; Yv = Ynew; tab.append((zr, *Yv)); break
        zr += h; Yv = Ynew; tab.append((zr, *Yv))
    tab = np.array(tab); zt = tab[:, 0] - tab[-1, 0]                                   # shell at z = 0
    yt, rt, pt, st = tab[:, 1], tab[:, 2], tab[:, 3], tab[:, 4]
    rpt = np.sqrt(1 + rt*rt*(st*st/12 - U(pt)/6))
    out = {}
    # Hermite interpolation in z on the ODE table; series region below the table
    zmin_tab = zt[0]; u = U(phi_h)
    inside = zgrid >= zmin_tab
    idx = np.clip(np.searchsorted(zt, zgrid[inside]) - 1, 0, len(zt) - 2)
    h = zt[idx+1] - zt[idx]; tt = (zgrid[inside] - zt[idx])/h
    h00 = 2*tt**3 - 3*tt**2 + 1; h10 = tt**3 - 2*tt**2 + tt; h01 = -2*tt**3 + 3*tt**2; h11 = tt**3 - tt**2
    def herm(f, fz): return h00*f[idx] + h10*h*fz[idx] + h01*f[idx+1] + h11*h*fz[idx+1]
    lnr = np.empty_like(zgrid); ph = np.empty_like(zgrid); yy = np.empty_like(zgrid)
    lnr[inside] = herm(np.log(rt), rpt); ph[inside] = herm(pt, rt*st); yy[inside] = herm(yt, rt)
    zo = zgrid[~inside]
    if zo.size:
        c0 = zmin_tab - (math.log(y0) + u*y0*y0/72)                                   # z = c0 + ln y + u y^2/72  (from dz = dy/rho)
        yv = np.exp(zo - c0)
        for _ in range(3): yv = np.exp(zo - c0 - u*yv*yv/72)
        r_, p_, s_ = cone_series(phi_h, yv); lnr[~inside] = np.log(r_); ph[~inside] = p_; yy[~inside] = yv
    rho = np.exp(lnr)
    # phi_y from the ODE table / series
    sy = np.empty_like(zgrid); sy[inside] = herm(st, rt*(U1f(pt) - 4*(rpt/rt)*st))
    if zo.size: sy[~inside] = cone_series(phi_h, yy[~inside])[2]
    rp = np.sqrt(1 + rho*rho*(sy*sy/12 - U(ph)/6))
    return dict(rho=rho, Hc=rp, phi=ph, phiz=rho*sy, lnrho=lnr)

# ---------------- finite differences (4th order) with ghost points on the right (shell) and zero deviation on the left
def d1(f, h, gR):
    """first derivative; gR = Neumann value f'(0) at the last point (shell)"""
    F = np.concatenate(([0.0, 0.0], f, [f[-2] + 2*h*gR, f[-3] + 4*h*gR]))      # left: deviation vanishes; right: even + linear-odd extension
    return (F[:-4] - 8*F[1:-3] + 8*F[3:-1] - F[4:])/(12*h)
def d2(f, h, gR):
    F = np.concatenate(([0.0, 0.0], f, [f[-2] + 2*h*gR, f[-3] + 4*h*gR]))
    return (-F[:-4] + 16*F[1:-3] - 30*F[2:-2] + 16*F[3:-1] - F[4:])/(12*h*h)
def ko(f, h, gR, eps):
    # Kreiss-Oliger dissipation  +eps * delta^6 f / (64 h): on the grid-scale mode delta^6 f = -64 f, i.e. damping at rate eps/h.
    # (A first version had the opposite sign and AMPLIFIED that mode at rate eps/h = 25 for dz = 2e-3; found with formulation_lab.py.)
    F = np.concatenate(([0.0, 0.0, 0.0], f, [f[-2] + 2*h*gR, f[-3] + 4*h*gR, f[-4] + 6*h*gR]))
    return +eps*(F[:-6] - 6*F[1:-5] + 15*F[2:-4] - 20*F[3:-3] + 15*F[4:-2] - 6*F[5:-1] + F[6:])/(64*h)

def run(d=0.0, dc=1e-6, dz=4e-4, L=8.0, t_final=12.0, t_det=1e-3, tag="run", out_dir=".", ko_eps=0.05, stop_phi=(-1.6, 0.995), log=print, pslope=False):
    """REFEREE COPY.  Changes relative to ../rolloff5d.py (all diagnostic, plus one optional switch):
       * constraint monitors: momentum constraint on the FULL grid (one-sided 4th-order pa_z at the last two points), reported
         separately for the interior (code's slice 10:-6) and for the shell layer (last 6 points); Hamiltonian constraint likewise;
       * a_b (deviation of A at the shell) is recorded so the brane scale factor can be reconstructed;
       * pslope=True: the ghost slope of the time-derivative fields (pa, pbv, pf) in the KO operator is the time derivative of the
         junction data, d/dt gA and d/dt gF, instead of zero (the original uses 0, which is inconsistent with A_tz = d/dt[(sigma/6)e^B]).
       Everything else is byte-identical to the original."""
    Ip = 1.0357712571566784; c = 2/Ip - 4/3
    ten = Tension(t_det, c, d); ten_seed = Tension(t_det, c + dc, d)
    g0 = (math.log10(0.2207*t_det**1.8), 8.2282 + 0.9*math.log(1e-3/t_det))
    ph_h, y_b, r0 = solve_shell(ten, guess=g0); ph_h2, y_b2, r2 = solve_shell(ten_seed, guess=(math.log10(ph_h + 1), y_b))
    N = int(round(L/dz)); z = -L + dz*np.arange(N + 1)
    S = static_on_grid(ph_h, y_b, z); S2 = static_on_grid(ph_h2, y_b2, z)
    rho, Hc, phs, phz = S["rho"], S["Hc"], S["phi"], S["phiz"]
    rho2 = rho*rho; Hb, rb, pb = Hc[-1], rho[-1], phs[-1]
    log("static shell: phi_h+1=%.6e y_b=%.8f rho_b=%.6f phi_b=%+.4e junction residuals=(%.1e,%.1e)" % (ph_h + 1, y_b, rb, pb, r0[0], r0[1]))
    log("junction check on grid: Hc_b - sigma rho_b/6 = %.2e ; phiz_b + sigma' rho_b/2 = %.2e" % (Hb - ten.s(pb)*rb/6, phz[-1] + ten.s1(pb)*rb/2))
    # deviation fields: a, b, f and their time derivatives
    a = S2["lnrho"] - S["lnrho"]; b = a.copy(); f = S2["phi"] - S["phi"]
    a[0:2] = b[0:2] = f[0:2] = 0.0
    pa = np.zeros_like(a); pbv = np.zeros_like(a); pf = np.zeros_like(a)
    log("seed: dc=%.1e  max|a|=%.2e max|f|=%.2e  f_b=%+.3e" % (dc, np.abs(a).max(), np.abs(f).max(), f[-1]))
    def neumann(bb, ff):
        eB = rb*math.exp(bb); p = pb + ff
        gA = (ten.s(p)*eB - ten.s(pb)*rb)/6; gF = -(ten.s1(p)*eB - ten.s1(pb)*rb)/2
        return gA, gF
    def neumann_t(bb, ff, pbb, pff):
        """time derivative of the junction data: slopes of pa (= of pbv) and of pf at the shell"""
        eB = rb*math.exp(bb); p = pb + ff
        gAt = (ten.s1(p)*pff + ten.s(p)*pbb)*eB/6
        gFt = -((2*W2(p) + ten.d)*pff + ten.s1(p)*pbb)*eB/2
        return gAt, gFt
    def rhs(state):
        a, b, f, pa, pbv, pf = state
        gA, gF = neumann(b[-1], f[-1])
        gAt, gFt = neumann_t(b[-1], f[-1], pbv[-1], pf[-1]) if pslope else (0.0, 0.0)
        az, bz, fz = d1(a, dz, gA), d1(b, dz, gA), d1(f, dz, gF)
        azz, bzz, fzz = d2(a, dz, gA), d2(b, dz, gA), d2(f, dz, gF)
        dU, dU1, U0, U10 = dU_exact(phs, f)
        e2b = np.expm1(2*b)
        SU = rho2*(e2b*(U0 + dU) + dU)                       # e^{2B}U(phi) - rho^2 U(phi_s)
        SU1 = rho2*(e2b*(U10 + dU1) + dU1)
        kin = 2*pa + pa*pa; grad = 2*Hc*az + az*az
        ra = azz - 3*kin + 3*grad + (2/3)*SU
        rbb = bzz + 3*kin - 3*grad - 0.5*pf*pf + 0.5*(2*phz*fz + fz*fz) - (1/3)*SU
        rf = fzz - 3*(1 + pa)*pf + 3*((Hc + az)*(phz + fz) - Hc*phz) - SU1
        out = [pa + ko(a, dz, gA, ko_eps), pbv + ko(b, dz, gA, ko_eps), pf + ko(f, dz, gF, ko_eps),
               ra + ko(pa, dz, gAt, ko_eps), rbb + ko(pbv, dz, gAt, ko_eps), rf + ko(pf, dz, gFt, ko_eps)]
        for o in out: o[0:2] = 0.0
        return out
    def d1_onesided_last(f, h):
        """4th-order one-sided first derivative at the last two points (no ghost assumption)"""
        g1 = (25*f[-1] - 48*f[-2] + 36*f[-3] - 16*f[-4] + 3*f[-5])/(12*h)
        g2 = (3*f[-1] + 10*f[-2] - 18*f[-3] + 6*f[-4] - f[-5])/(12*h)
        return g1, g2
    def constraints(state):
        """returns (M_interior, M_shell_layer, H_interior, H_shell_layer), all normalised by the sum of |terms|"""
        a, b, f, pa, pbv, pf = state
        gA, gF = neumann(b[-1], f[-1])
        az, bz, fz = d1(a, dz, gA), d1(b, dz, gA), d1(f, dz, gF); paz = d1(pa, dz, 0.0)
        paz[-1], paz[-2] = d1_onesided_last(pa, dz)                                    # referee: proper one-sided pa_z at the shell
        M = -3*paz - 3*(1 + pa)*(az - bz) + 3*(Hc + az)*pbv - pf*(phz + fz)
        scale = 3*np.abs(paz) + 3*np.abs(Hc + az)*(np.abs(1 + pa) + np.abs(pbv)) + 3*np.abs(1 + pa)*np.abs(Hc + bz) + np.abs(pf)*np.abs(phz + fz) + 1e-300
        azz = d2(a, dz, gA); dU, dU1, U0, U10 = dU_exact(phs, f); e2b = np.expm1(2*b)
        SU = rho2*(e2b*(U0 + dU) + dU)
        # Hc = -2Ue^{2B} + 6A_t^2 + 6A_tB_t - 12A_z^2 + 6A_zB_z - 6A_zz - phi_t^2 - phi_z^2, static part subtracted analytically:
        # static: -2 rho^2 U_s + 6 - 12 Hc^2 + 6 Hc^2 - 6 A_szz - phz^2 = 0 with A_szz = Hc^2 - 1 - phz^2/3 (background ODE)
        At_ = 1 + pa; Az_ = Hc + az; Bz_ = Hc + bz
        Hm = -2*SU + 6*(At_**2 - 1) + 6*At_*pbv - 12*(Az_**2 - Hc**2) + 6*(Az_*Bz_ - Hc**2) - 6*azz - pf**2 - (2*phz*fz + fz**2)
        Aszz = Hc**2 - 1 - phz**2/3
        hs = 2*np.abs(rho2*(e2b + 1)*(U0 + dU)) + 6*At_**2 + 6*np.abs(At_*pbv) + 12*Az_**2 + 6*np.abs(Az_*Bz_) + 6*np.abs(azz + Aszz) + pf**2 + (phz + fz)**2 + 1e-300
        sl = slice(10, -6); se = slice(-6, None)
        return (float(np.max(np.abs(M[sl])/scale[sl])), float(np.max(np.abs(M[se])/scale[se])),
                float(np.max(np.abs(Hm[sl])/hs[sl])), float(np.max(np.abs(Hm[se])/hs[se])))
    dt = 0.5*dz; nsteps = int(round(t_final/dt)); state = [a, b, f, pa, pbv, pf]
    rec = []; t = 0.0; tau = 0.0; t0 = time.time(); every = max(1, int(round(0.01/dt))); stop_reason = "t_final"
    snaps = {}; snap_times = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    for n in range(nsteps + 1):
        if n % every == 0:
            a, b, f, pa, pbv, pf = state
            eB = rb*math.exp(b[-1]); HJ = (1 + pa[-1])/eB                       # brane-frame Hubble rate  (dA/dt) e^{-B}
            phib = pb + f[-1]; dphib_dtau = pf[-1]/eB
            Mx, Ms, Hx, Hs = constraints(state) if n % (every*10) == 0 else (float("nan"),)*4
            rec.append((t, phib, HJ*rb, b[-1], pa[-1], dphib_dtau*rb, Mx, float(np.max(np.abs(f))), float(np.max(np.abs(b))), tau, float(np.min(phs + f)), float(np.max(phs + f)),
                        Ms, Hx, Hs, a[-1]))
            if n % (every*50) == 0:
                log("t=%6.2f H0tau=%7.4f phi_b=%+.6e H_J/H0=%+.6f b_b=%+.3e a_b=%+.3e phi_min=%+.4f M=%.1e Mshell=%.1e H=%.1e Hshell=%.1e  [%.0fs]" % (t, tau, phib, HJ*rb, b[-1], a[-1], np.min(phs + f), Mx, Ms, Hx, Hs, time.time() - t0))
            if not np.isfinite(phib) or not np.isfinite(HJ): stop_reason = "non-finite"; break
            if phib < stop_phi[0] or phib > stop_phi[1]: stop_reason = "phi_b left [%g,%g]" % stop_phi; break
            for ts in snap_times:
                if abs(t - ts) < dt/2 and ts not in snaps: snaps[ts] = (f[::8].copy(), b[::8].copy(), a[::8].copy())
        if n == nsteps: break
        k1 = rhs(state); k2 = rhs([s_ + dt/2*k for s_, k in zip(state, k1)]); k3 = rhs([s_ + dt/2*k for s_, k in zip(state, k2)]); k4 = rhs([s_ + dt*k for s_, k in zip(state, k3)])
        tau += dt*0.5*(math.exp(state[1][-1]))
        state = [s_ + dt/6*(q1 + 2*q2 + 2*q3 + q4) for s_, q1, q2, q3, q4 in zip(state, k1, k2, k3, k4)]; t += dt
        tau += dt*0.5*(math.exp(state[1][-1]))
        if n % 20000 == 0 and n > 0:
            np.savez_compressed(os.path.join(out_dir, tag + "_checkpoint.npz"), t=t, z=z, state=np.array(state), rec=np.array(rec))
    rec = np.array(rec)
    np.savez_compressed(os.path.join(out_dir, tag + "_timeseries.npz"), rec=rec, z=z[::8], phis=phs[::8], rho=rho[::8],
                        **{"snap_f_%g" % k: v[0] for k, v in snaps.items()}, **{"snap_b_%g" % k: v[1] for k, v in snaps.items()})
    # growth-rate fit on the linear stage
    tt, dphi = rec[:, 0], np.abs(rec[:, 1] - pb); m = (dphi > 30*abs(f[-1] if False else (S2["phi"][-1] - S["phi"][-1]))) & (dphi < 2e-3) & (tt > 1.0)
    fit = None
    if m.sum() > 20:
        sfit = np.polyfit(tt[m], np.log(dphi[m]), 1); fit = dict(growth_rate=float(sfit[0]), window=[float(tt[m][0]), float(tt[m][-1])], n=int(m.sum()))
    summ = dict(H0_tau_end=float(rec[-1, 9]), tag=tag, d=d, dc=dc, dz=dz, L=L, t_det=t_det, t_end=float(rec[-1, 0]), stop_reason=stop_reason, phi_b_static=float(pb), rho_b=float(rb),
                phi_b_end=float(rec[-1, 1]), HJ_over_H0_end=float(rec[-1, 2]), growth_fit=fit, max_momentum_constraint=float(np.nanmax(rec[:, 6])), max_M_shell=float(np.nanmax(rec[:, 12])), max_H_interior=float(np.nanmax(rec[:, 13])), max_H_shell=float(np.nanmax(rec[:, 14])), pslope=pslope,
                runtime_s=round(time.time() - t0, 1))
    json.dump(summ, open(os.path.join(out_dir, tag + "_summary.json"), "w"), indent=1); log(json.dumps(summ))
    return summ

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--d", type=float, default=0.0); ap.add_argument("--dc", type=float, default=1e-6); ap.add_argument("--dz", type=float, default=4e-4)
    ap.add_argument("--tdet", type=float, default=1e-3); ap.add_argument("--L", type=float, default=8.0); ap.add_argument("--tf", type=float, default=12.0); ap.add_argument("--tag", default="run"); ap.add_argument("--ko", type=float, default=0.05); ap.add_argument("--phimin", type=float, default=-1.6); ap.add_argument("--phimax", type=float, default=0.995); ap.add_argument("--pslope", type=int, default=0)
    A_ = ap.parse_args()
    run(d=A_.d, dc=A_.dc, dz=A_.dz, L=A_.L, t_final=A_.tf, t_det=A_.tdet, tag=A_.tag, ko_eps=A_.ko, stop_phi=(A_.phimin, A_.phimax), log=lambda s: print(s, flush=True), pslope=bool(A_.pslope))
