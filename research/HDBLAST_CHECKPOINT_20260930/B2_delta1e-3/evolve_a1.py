#!/usr/bin/env python3
"""A1 solver: coupled shell roll-off (5D Einstein-scalar + shell radiation) in a chart that continues through the
old 'chart freeze'.  Floating point; numerical results only.

Derived from research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/rolloff5d_matter.py (read-only) with these changes:
  * NO phi_b > 0.995 stop (the tuned endpoint is phi_b ~ 1.036).
  * Chart.  The pilot's conformal chart (shell fixed at z = 0) loses the shell at finite proper time: the lapse decays as
    e^{-t} (D0_FREEZE_DIAGNOSIS.json: d B_b/dt -> -1.000), i.e. the shell reaches the future light cone of the static
    vertex.  Here the same conformal-gauge equations are solved in the chart U = F(u), V = F(v) (one common F on both null
    coordinates, as in PROPER_CLOCK_PROOF.md, folder 152), F(x) = ln(1 + e^{-xc}) - ln(e^{-x} + e^{-xc}).  F' ~ e^{-(x-xc)}
    at late times cancels the e^{-t} decay, so the shell lapse stays O(1) (a near-proper clock: dtau/dT -> const), and the
    chart is regular across the vertex light cone.  xc = inf reproduces the pilot's chart exactly.
  * Initial data.  The pilot's seed (static shell with c + dc, zero velocity deviation) is used without taper: in the new
    chart it is the exact static solution evaluated on the T = 0 slice (static_w.reference), including the part of the slice
    beyond the vertex light cone (w < 0).  The far boundary Z = -L carries the exact static data, which is exact while
    V = T - L < 0 (the shell is unaffected for T < 2L).
  * Deviation form.  For T < T_switch the solver evolves d = X - S(T) (S = exact seed solution in the new chart) with analytic
    derivatives of S and finite differences of d, so the static seed is preserved to round-off (as in Chat 13/14).  At
    T_switch (before S becomes singular at V = F_inf) the state is converted to full fields X and evolved with finite
    differences of X; the far nodes keep exact data.
  * Shell junctions of the 22 Sept matter extension (as in the pilot): n.A = (sigma+R)/6, n.B = (sigma-3R)/6,
    n.phi = -(sigma' + Y v)/2; radiation dR/dtau + 4 H R = Y v^2.  The time derivatives of the Neumann data (used for the
    velocity ghosts: momentum constraint and damping) now include B_T, R_T and the Y phi_TT source term (missing in the pilot).
  * Operators: 4th-order centred differences on the mapped grid with a degree-5 Hermite ghost at the shell (lab.py order 4).
  * Remedy A (lab.py): discrete Hamiltonian-constraint projection of the initial warp (A and B shifted equally) on a
    window next to the shell.  Remedy B (lab.py): outgoing-characteristic damping -kappa (H + 2M)/(6 (A_T + A_Z)) added to
    B_TT, switched off smoothly for Z > -damp_off[0] (C^2 step to full kappa at Z < -damp_off[1]).
Units: kappa_5 = 1, model length; H0 = 1/rho_b of the unperturbed static shell (tension c, d).

B2 copy (30 Sept 2026, HDBLAST_CHECKPOINT_20260930/B2_delta1e-3): changes relative to A1 evolve_a1.py, all numerical-method only:
  * --table_dx: spacing of the static-reference tables (A1: 2.5e-4; B2 default 2.5e-5).  At delta = 1e-3 the A1 tables violate
    the shell junction by 1.8e-10 (a spurious tension shift dc_eff = 3.5e-7, 35x the growth-calibration seed; diag/D1_TABLE_KICK.json).
  * identity chart (xc = inf): the reference is exactly static apart from A -> A + T, so it is evaluated once and cached
    (checked against a fresh evaluation at the first call; speed only).
  * --damp_off z1 z2 (Remedy B switch zone), --zwidth (fine/coarse transition width), optional middle grid level --dzm/--zmid.
  * records the growth-rate helpers unchanged; output goes to this folder.
  * (1 Oct, dated note 3) --smax exposes the bounded chart's saturation constant s_max (default 20 as in A1; exploratory runs only).
"""
import json, math, os, sys, time, hashlib
from pathlib import Path
import numpy as np
from scipy import sparse
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import static_w as S

# ---------------------------------------------------------------- grid and operators (after lab.py / rolloff5d_v2.py)
def build_grid(L, dz_fine, dz_coarse, z_fine=0.4, width=0.02, dz_mid=None, z_mid=None, width_mid=None):
    sig = lambda z: 0.5*(1 + np.tanh((z + z_fine)/width))
    dsig = lambda z: 0.5*(1 - np.tanh((z + z_fine)/width)**2)/width
    if dz_mid is None:
        s = lambda z: dz_coarse + (dz_fine - dz_coarse)*sig(z)
        sp_ = lambda z: (dz_fine - dz_coarse)*dsig(z)
    else:                       # three levels: dz_fine (z > -z_fine), dz_mid (-z_mid < z < -z_fine), dz_coarse beyond
        wm = width if width_mid is None else width_mid
        sig2 = lambda z: 0.5*(1 + np.tanh((z + z_mid)/wm))
        dsig2 = lambda z: 0.5*(1 - np.tanh((z + z_mid)/wm)**2)/wm
        s = lambda z: dz_coarse + (dz_mid - dz_coarse)*sig2(z) + (dz_fine - dz_mid)*sig(z)
        sp_ = lambda z: (dz_mid - dz_coarse)*dsig2(z) + (dz_fine - dz_mid)*dsig(z)
    zs = [0.0]; z = 0.0
    while z > -L:
        k1 = -s(z); k2 = -s(z + 0.5*k1); k3 = -s(z + 0.5*k2); k4 = -s(z + k3)
        z = z + (k1 + 2*k2 + 2*k3 + k4)/6; zs.append(z)
    zs = np.array(zs[::-1]); zp = s(zs); zpp = sp_(zs)*zp
    return zs, zp, zpp

def build_grid_rho(L, dz_fine, dz_coarse, rho_fn, dz_cap=None, z_cap=3.0, w_cap=0.5):
    """B2 grid: spacing tied to the static background's conformal factor, s(z) = 1/(rho(z)/(rho_b dz_fine) + 1/dz_coarse),
    i.e. a fixed number of points per local wall/curvature width (which scales as 1/rho(z) in the conformal coordinate) and a
    smooth cap dz_coarse far from the shell.  rho_fn(z) -> (rho, h = rho_z/rho) of the static background (z <= 0)."""
    rb = float(rho_fn(np.array([0.0]))[0][0])
    # optional (B2 dated note): the cap is dz_cap for z > -z_cap and rises smoothly to dz_coarse beyond (width w_cap)
    def cap(z):
        if dz_cap is None: return np.full_like(z, dz_coarse), np.zeros_like(z)
        t = np.tanh((z + z_cap)/w_cap)
        return dz_cap + (dz_coarse - dz_cap)*0.5*(1 - t), -(dz_coarse - dz_cap)*0.5*(1 - t*t)/w_cap
    def s(z):
        z = np.atleast_1d(np.asarray(z, float)); r, h = rho_fn(z); q = r/(rb*dz_fine); cp, _ = cap(z)
        return 1.0/(q + 1.0/cp)
    def sp_(z):
        z = np.atleast_1d(np.asarray(z, float)); r, h = rho_fn(z); q = r/(rb*dz_fine); cp, cp1 = cap(z)
        return -(q*h - cp1/cp**2)/(q + 1.0/cp)**2
    zs = [0.0]; z = 0.0
    while z > -L:
        k1 = -s(z)[0]; k2 = -s(z + 0.5*k1)[0]; k3 = -s(z + 0.5*k2)[0]; k4 = -s(z + k3)[0]
        z = z + (k1 + 2*k2 + 2*k3 + k4)/6; zs.append(z)
    zs = np.array(zs[::-1]); zp = s(zs); zpp = sp_(zs)*zp
    return zs, zp, zpp

D1W = {4: [1/12, -8/12, 0, 8/12, -1/12], 6: [-1/60, 9/60, -45/60, 0, 45/60, -9/60, 1/60]}
D2W = {4: [-1/12, 16/12, -30/12, 16/12, -1/12], 6: [2/180, -27/180, 270/180, -490/180, 270/180, -27/180, 2/180]}
def hermite_ghost_weights(order):
    deg = order + 1; nv = order + 1
    V = np.zeros((deg + 1, deg + 1)); xx = -np.arange(nv, dtype=float)
    V[:nv] = xx[:, None]**np.arange(deg + 1)[None, :]; V[nv, 1] = 1.0
    r = order//2
    X = np.array([(k + 1.0)**np.arange(deg + 1) for k in range(r)])
    return X @ np.linalg.inv(V)
def build_operators(zp, zpp, order=4):
    """identical construction to evolution_constraints/lab.py build_operators (checked in t2_solver_controls.py)"""
    n = len(zp); r = order//2; gh = hermite_ghost_weights(order); nv = order + 1
    rows = list(range(r, n + r)); cols = list(range(n)); vals = [1.0]*n
    gv = np.zeros(n + 2*r)
    for k in range(r):
        for j in range(nv):
            rows.append(n + r + k); cols.append(n - 1 - j); vals.append(gh[k, j])
        gv[n + r + k] = gh[k, nv]*zp[-1]
    E_ = sparse.coo_matrix((vals, (rows, cols)), shape=(n + 2*r, n)).tocsr()
    D1x = sum(c*E_[k:k + n] for k, c in enumerate(D1W[order]) if c != 0)
    D2x = sum(c*E_[k:k + n] for k, c in enumerate(D2W[order]) if c != 0)
    g1 = sum(c*gv[k:k + n] for k, c in enumerate(D1W[order]))
    g2 = sum(c*gv[k:k + n] for k, c in enumerate(D2W[order]))
    D1 = sparse.diags(1/zp) @ D1x
    D2 = sparse.diags(1/zp**2) @ D2x - sparse.diags(zpp/zp**3) @ D1x
    return D1.tocsr(), D2.tocsr(), g1/zp, g2/zp**2 - g1*zpp/zp**3

# ---------------------------------------------------------------- potential differences without cancellation
_Uc = np.array([-1/6, 4/3, -5/3, -4/9, 17/18, 0.0, -2/27])      # U = sum c_k phi^k (registered sextic)
def _poly_derivs(p, order):
    out = []; c = _Uc.copy()
    for j in range(order + 1):
        out.append(sum(ck*p**k for k, ck in enumerate(c)))
        c = np.array([k*c[k] for k in range(1, len(c))]) if len(c) > 1 else np.array([0.0])
    return out
def dU_exact(ps, f):
    D = _poly_derivs(ps, 6)
    dU = sum(D[j]*f**j/math.factorial(j) for j in range(1, 7))
    dU1 = sum(D[j + 1]*f**j/math.factorial(j) for j in range(1, 6))
    return dU, dU1, D[0], D[1]

class Run:
    def __init__(self, delta=0.1, c=S.C_REG, d=0.0, dc=1e-2, Y=0.0, dzf=1e-3, dzc=4e-3, zfine=0.4, width=0.02, L=16.0,
                 xc=math.inf, T_switch=None, kappa=0.0, damp_off=(0.02, 0.05), project=None, cfl=0.5, order=4,
                 shells=None, log=print, chart_kind='bounded', table_dx=2.5e-5, dzm=None, zmid=None, width_mid=None, grid_kind='tanh', dzcap=None, zcap=3.0, s_max=20.0):
        self.p = dict(delta=delta, c=c, d=d, dc=dc, Y=Y, dzf=dzf, dzc=dzc, zfine=zfine, width=width, L=L, xc=xc, chart_kind=chart_kind,
                      T_switch=T_switch, kappa=kappa, damp_off=list(damp_off), project=project, cfl=cfl, order=order,
                      table_dx=table_dx, dzm=dzm, zmid=zmid, width_mid=width_mid, grid_kind=grid_kind, dzcap=dzcap, zcap=zcap, s_max=s_max)
        self.log = log; self.Y = Y
        self.ten = S.Tension(delta, c, d)
        t0 = time.time()
        cache_f = HERE/'shell_cache'/('shells_%r_%r_%r_%r_%r.pkl' % (delta, c, d, dc, table_dx))
        if shells is None and cache_f.exists():          # B2 (speed only): static shells are cached on disk after the first solve
            import pickle
            shells = pickle.load(open(cache_f, 'rb')); log('static shells loaded from %s' % cache_f.name)
        if shells is None:
            import rolloff5d_v1 as R_
            g0 = (math.log10(0.2207*delta**1.8), 8.2282 + 0.9*math.log(1e-3/delta))
            ph_h, y_b, _ = R_.solve_shell(R_.Tension(delta, c, d), guess=g0)
            ph_h2, y_b2, _ = R_.solve_shell(R_.Tension(delta, c + dc, d), guess=(math.log10(ph_h + 1), y_b))
            self.bg = S.StaticShell(S.Tension(delta, c, d), phh_guess=ph_h, table_dx=table_dx)
            self.seed = S.StaticShell(S.Tension(delta, c + dc, d), phh_guess=ph_h2, table_dx=table_dx) if dc != 0 else self.bg
            try:
                import pickle; cache_f.parent.mkdir(exist_ok=True)
                pickle.dump((self.bg, self.seed), open(str(cache_f) + '.tmp', 'wb')); os.replace(str(cache_f) + '.tmp', cache_f)
            except Exception as err: log('shell cache not written: %s' % err)
        else:
            self.bg, self.seed = shells
        self.rb = self.bg.rhob; self.H0 = 1/self.rb
        log('static shells: rho_b=%.10f phi_b=%.6e (bg); seed phi_b=%.6e rho_b=%.10f  [%.0fs]' % (self.bg.rhob, self.bg.phb, self.seed.phb, self.seed.rhob, time.time() - t0))
        self.chart = S.Chart(xc, chart_kind, s_max=s_max)
        Fi = self.chart.F_inf
        self.T_switch = (math.inf if not math.isfinite(Fi) else Fi - 1.0) if T_switch is None else T_switch
        if math.isfinite(Fi) and self.T_switch > Fi - 0.3: raise ValueError('T_switch too close to F_inf')
        if grid_kind == 'rho':
            def rho_fn(zz, sh=self.bg):
                r, h, _, _ = sh.profile_z(np.minimum(zz, 0.0)); return r, h
            z, zp, zpp = build_grid_rho(L, dzf, dzc, rho_fn, dz_cap=dzcap, z_cap=zcap)
        else:
            z, zp, zpp = build_grid(L, dzf, dzc, zfine, width, dzm, zmid, width_mid)
        self.z, self.zp, self.zpp, self.n = z, zp, zpp, len(z)
        self.D1, self.D2, self.g1, self.g2 = build_operators(zp, zpp, order)
        self.nfix = 2
        x = np.clip((-z - damp_off[0])/(damp_off[1] - damp_off[0]), 0, 1)
        self.kz = kappa*x**3*(10 - 15*x + 6*x*x)
        self.kappa = kappa
        self.inner = slice(6, self.n - 6)
        self.mode = 'dev' if self.T_switch > 0 else 'full'
        log('grid (%s): %d points, dz %.2e (shell) .. %.2e (far), L=%.2f, %.1f points per wall width 1/rho_b; chart xc=%s F_inf=%s; T_switch=%s; kappa=%g' % (
            grid_kind, self.n, zp[-1], zp[0], L, (1/self.rb)/zp[-1], xc, Fi, self.T_switch, kappa))

    # ------------------------------------------------------------ reference and state helpers
    def ref(self, T):
        if self.chart.kind == 'identity':
            if getattr(self, '_rf0', None) is None:
                self._rf0 = S.reference(self.seed, self.chart, np.zeros(self.n), self.z)
                chk = S.reference(self.seed, self.chart, np.full(self.n, 1.234), self.z)
                err = max(float(np.max(np.abs(chk[k] - (self._rf0[k] + (1.234 if k == 'A' else 0.0)))/(1.0 + np.abs(chk[k])))) for k in chk)
                self.ref_cache_check = err
                self.log('identity-chart reference cache: max |fresh - cached|/(1+|value|) at T=1.234: %.2e (round-off)' % err)
                if err > 1e-10: raise RuntimeError('reference is not static in the identity chart')
            r = dict(self._rf0); r['A'] = self._rf0['A'] + T
            return r
        # B2 (speed only): RK4 evaluates the reference at T, T+h/2 (twice) and T+h, and T+h is the next step's T; keep the last 4
        if not hasattr(self, '_rcache'): self._rcache = {}
        key = float(T)
        r = self._rcache.get(key)
        if r is None:
            r = S.reference(self.seed, self.chart, np.full(self.n, T), self.z)
            if len(self._rcache) >= 4: self._rcache.pop(next(iter(self._rcache)))
            self._rcache[key] = r
        return r
    def shell_B(self, T, Yv):
        if self.mode == 'dev':
            return float(Yv[1, -1] + self.ref_nodes(T, np.array([self.n - 1]))['B'][0])
        return float(Yv[1, -1])
    def ref_nodes(self, T, idx):
        if self.chart.kind == 'identity':
            r0 = self.ref(0.0)
            r = {k: v[idx] for k, v in r0.items()}; r['A'] = r['A'] + T
            return r
        return S.reference(self.seed, self.chart, np.full(len(idx), T), self.z[idx])

    def neumann(self, Bb, phib, phiTb, R):
        eB = math.exp(Bb); t = self.ten
        return (eB*(t.s(phib) + R)/6, eB*(t.s(phib) - 3*R)/6, -(eB*t.s1(phib) + self.Y*phiTb)/2)

    def neumann_t(self, Bb, BTb, phib, phiTb, phiTTb, R, RT):
        eB = math.exp(Bb); t = self.ten; s0, s1, s2 = t.s(phib), t.s1(phib), t.s2(phib)
        gAt = eB*(BTb*(s0 + R) + s1*phiTb + RT)/6
        gBt = eB*(BTb*(s0 - 3*R) + s1*phiTb - 3*RT)/6
        gFt = -(eB*(BTb*s1 + s2*phiTb) + self.Y*phiTTb)/2
        return gAt, gBt, gFt

    def fields(self, T, Yv, R, need_constraints=False, ref=None):
        """full fields and derivatives on the grid.  Yv: (6, n) state (deviation or full, per self.mode)."""
        if self.mode == 'dev':
            rf = self.ref(T) if ref is None else ref
            A = rf['A'] + Yv[0]; B = rf['B'] + Yv[1]; ph = rf['phi'] + Yv[2]
            AT = rf['A_T'] + Yv[3]; BT = rf['B_T'] + Yv[4]; phT = rf['phi_T'] + Yv[5]
        else:
            rf = None
            A, B, ph, AT, BT, phT = Yv
        gA, gB, gF = self.neumann(B[-1], ph[-1], phT[-1], R)
        if self.mode == 'dev':
            gd = np.array([gA - rf['A_Z'][-1], gB - rf['B_Z'][-1], gF - rf['phi_Z'][-1]])
            d1 = (self.D1 @ Yv[:3].T).T + self.g1[None, :]*gd[:, None]
            d2 = (self.D2 @ Yv[:3].T).T + self.g2[None, :]*gd[:, None]
            AZ, BZ, phZ = rf['A_Z'] + d1[0], rf['B_Z'] + d1[1], rf['phi_Z'] + d1[2]
            AZZ, BZZ, phZZ = rf['A_ZZ'] + d2[0], rf['B_ZZ'] + d2[1], rf['phi_ZZ'] + d2[2]
        else:
            g = np.array([gA, gB, gF])
            d1 = (self.D1 @ Yv[:3].T).T + self.g1[None, :]*g[:, None]
            d2 = (self.D2 @ Yv[:3].T).T + self.g2[None, :]*g[:, None]
            AZ, BZ, phZ = d1; AZZ, BZZ, phZZ = d2
        return dict(rf=rf, A=A, B=B, phi=ph, A_T=AT, B_T=BT, phi_T=phT, A_Z=AZ, B_Z=BZ, phi_Z=phZ, A_ZZ=AZZ, B_ZZ=BZZ, phi_ZZ=phZZ,
                    gA=gA, gB=gB, gF=gF)

    def accel(self, F):
        """second time derivatives from the conformal-gauge equations (full fields F)"""
        e2B = np.exp(2*F['B']); Uv = S.U(F['phi']); U1v = S.U1(F['phi'])
        Att = F['A_ZZ'] - 3*F['A_T']**2 + 3*F['A_Z']**2 + (2/3)*e2B*Uv
        Btt = F['B_ZZ'] + 3*F['A_T']**2 - 3*F['A_Z']**2 - 0.5*F['phi_T']**2 + 0.5*F['phi_Z']**2 - e2B*Uv/3
        ftt = F['phi_ZZ'] - 3*F['A_T']*F['phi_T'] + 3*F['A_Z']*F['phi_Z'] - e2B*U1v
        return Att, Btt, ftt

    def accel_dev(self, F):
        """d_TT = N(S + d) - N(S) with cancellation-free differences (S analytic)"""
        rf = F['rf']
        a_T, a_Z = F['A_T'] - rf['A_T'], F['A_Z'] - rf['A_Z']
        f_T, f_Z = F['phi_T'] - rf['phi_T'], F['phi_Z'] - rf['phi_Z']
        b = F['B'] - rf['B']; f = F['phi'] - rf['phi']
        U0 = S.U(rf['phi']); U10 = S.U1(rf['phi'])
        dU = S.U(F['phi']) - U0; dU1 = S.U1(F['phi']) - U10          # exactly 0 when f = 0
        e2r = np.exp(2*rf['B']); e2b = np.expm1(2*b)
        SU = e2r*(e2b*(U0 + dU) + dU); SU1 = e2r*(e2b*(U10 + dU1) + dU1)
        kin = a_T*(2*rf['A_T'] + a_T); grad = a_Z*(2*rf['A_Z'] + a_Z)
        ra = (F['A_ZZ'] - rf['A_ZZ']) - 3*kin + 3*grad + (2/3)*SU
        rb = (F['B_ZZ'] - rf['B_ZZ']) + 3*kin - 3*grad - 0.5*f_T*(2*rf['phi_T'] + f_T) + 0.5*f_Z*(2*rf['phi_Z'] + f_Z) - SU/3
        rfv = (F['phi_ZZ'] - rf['phi_ZZ']) - 3*(F['A_T']*F['phi_T'] - rf['A_T']*rf['phi_T']) + 3*(F['A_Z']*F['phi_Z'] - rf['A_Z']*rf['phi_Z']) - SU1
        return ra, rb, rfv

    def constraints(self, F, velocity_state, R, RT, phiTTb):
        """Hamiltonian and momentum constraints (full), with the velocity ghost from the complete Neumann_t."""
        gAt, gBt, gFt = self.neumann_t(F['B'][-1], F['B_T'][-1], F['phi'][-1], F['phi_T'][-1], phiTTb, R, RT)
        if self.mode == 'dev':
            rf = F['rf']
            AtZ = rf['A_TZ'] + self.D1 @ velocity_state + self.g1*(gAt - rf['A_TZ'][-1])
        else:
            AtZ = self.D1 @ velocity_state + self.g1*gAt
        e2B = np.exp(2*F['B']); Uv = S.U(F['phi'])
        M = -3*AtZ - 3*F['A_T']*F['A_Z'] + 3*F['A_T']*F['B_Z'] + 3*F['A_Z']*F['B_T'] - F['phi_T']*F['phi_Z']
        H = (-2*Uv*e2B + 6*F['A_T']**2 + 6*F['A_T']*F['B_T'] - 12*F['A_Z']**2 + 6*F['A_Z']*F['B_Z'] - 6*F['A_ZZ']
             - F['phi_T']**2 - F['phi_Z']**2)
        sM = 3*np.abs(AtZ) + 3*np.abs(F['A_T']*F['A_Z']) + 3*np.abs(F['A_T']*F['B_Z']) + 3*np.abs(F['A_Z']*F['B_T']) + np.abs(F['phi_T']*F['phi_Z']) + 1e-300
        sH = (2*np.abs(Uv)*e2B + 6*F['A_T']**2 + 6*np.abs(F['A_T']*F['B_T']) + 12*F['A_Z']**2 + 6*np.abs(F['A_Z']*F['B_Z'])
              + 6*np.abs(F['A_ZZ']) + F['phi_T']**2 + F['phi_Z']**2 + 1e-300)
        return H, M, sH, sM

    def rhs(self, T, Yv, R):
        F = self.fields(T, Yv, R)
        if self.mode == 'dev':
            ra, rb, rfv = self.accel_dev(F)
        else:
            ra, rb, rfv = self.accel(F)
        eBb = math.exp(F['B'][-1])
        RT = -4*F['A_T'][-1]*R + self.Y*F['phi_T'][-1]**2/eBb
        if self.kappa:
            phiTTb = rfv[-1] + (F['rf']['phi_TT'][-1] if self.mode == 'dev' else 0.0)
            H, M, _, _ = self.constraints(F, Yv[3], R, RT, phiTTb)
            den = F['A_T'] + F['A_Z']
            den = np.where(np.abs(den) < 0.1, np.sign(den)*0.1 + (den == 0)*0.1, den)
            rb = rb - self.kz*(H + 2*M)/(6*den)
        out = np.empty_like(Yv)
        out[:3] = Yv[3:]; out[3] = ra; out[4] = rb; out[5] = rfv
        if self.mode == 'dev':
            out[:, :self.nfix] = 0.0
        else:
            idx = np.arange(self.nfix)
            rn = self.ref_nodes(T, idx)
            out[0, idx] = rn['A_T']; out[1, idx] = rn['B_T']; out[2, idx] = rn['phi_T']
            out[3, idx] = rn['A_TT']; out[4, idx] = rn['B_TT']; out[5, idx] = rn['phi_TT']
        return out, RT, F

    # ------------------------------------------------------------ Remedy A
    def project(self, Yv, R, z_left, tol=1e-12, maxit=10):
        """Newton solve of the discrete H_i = 0 (window z >= z_left, excluding frozen nodes) for a warp shift delta added to
        both A and B (velocities fixed).  Jacobian by 5-colouring (the operator has bandwidth 2 plus the shell ghosts)."""
        idx = np.flatnonzero(self.z >= z_left); idx = idx[idx >= self.nfix + 2]
        v = Yv.copy(); logs = []
        def Hwin(vv):
            F = self.fields(0.0, vv, R)
            eBb = math.exp(F['B'][-1]); RT = -4*F['A_T'][-1]*R + self.Y*F['phi_T'][-1]**2/eBb
            H, M, sH, sM = self.constraints(F, vv[3], R, RT, 0.0)
            return H, M, sH, sM
        for it in range(maxit):
            H, M, sH, sM = Hwin(v)
            res = H[idx]
            logs.append(dict(iteration=it, max_rel_H_window=float(np.max(np.abs(res)/sH[idx])), max_rel_M_window=float(np.max(np.abs(M[idx])/sM[idx]))))
            if np.max(np.abs(res)/sH[idx]) < tol: break
            eps = 1e-7
            pos = {j: k for k, j in enumerate(idx)}
            rows_, cols_, vals_ = [], [], []
            for col in range(5):
                cols = idx[(idx % 5) == col]
                vp = v.copy(); vp[0, cols] += eps; vp[1, cols] += eps
                Hp = Hwin(vp)[0][idx]
                dH = (Hp - res)/eps
                for j in cols:
                    near = set(range(j - 2, j + 3))
                    # the Hermite ghosts use the last 5 nodes, so rows n-2, n-1 depend on columns n-5..n-1;
                    # within any 5 consecutive columns the colours are distinct, so the colouring stays exact
                    if j >= self.n - 5: near |= {self.n - 2, self.n - 1}
                    for i in near:
                        if i in pos: rows_.append(pos[i]); cols_.append(pos[j]); vals_.append(dH[pos[i]])
            J = sparse.csc_matrix((vals_, (rows_, cols_)), shape=(len(idx), len(idx)))
            from scipy.sparse.linalg import spsolve
            dd = spsolve(J, -res)
            v[0, idx] += dd; v[1, idx] += dd
        H, M, sH, sM = Hwin(v)
        logs.append(dict(final=True, max_rel_H_window=float(np.max(np.abs(H[idx])/sH[idx])), max_rel_M_window=float(np.max(np.abs(M[idx])/sM[idx])),
                         max_abs_shift=float(np.max(np.abs(v[0] - Yv[0])))))
        return v, logs

    # ------------------------------------------------------------ evolution
    def evolve(self, T_final=12.0, rec_dT=0.01, tag='run', out_dir='.', snap_dT=1.0, wall_limit=1500.0,
               stop_H_over_H0=(-20.0, 20.0), constraint_every=10, stop_dbdt=None, ckpt_path=None, resume=False, ckpt_every=240.0):
        # B2: checkpoint/restart (each process stays below the per-run wall limit; a run may be continued with --resume)
        import pickle
        t_start = time.time(); elapsed0 = 0.0; resumed_from = None
        Yv = np.zeros((6, self.n)); R = 0.0; T = 0.0
        if self.mode == 'full':      # full fields from the reference at T = 0
            rf = self.ref(0.0); Yv = np.array([rf['A'], rf['B'], rf['phi'], rf['A_T'], rf['B_T'], rf['phi_T']])
        proj_log = None
        dt = self.p['cfl']*self.zp.min(); nsteps = int(math.ceil(T_final/dt))
        every = max(1, int(round(rec_dT/dt)))
        rows = []; snaps = {}; tau = 0.0; stop = 'T_final'
        A0 = None; ln_a = 0.0
        nsnap = 0; n = 0
        if resume and ckpt_path and os.path.exists(ckpt_path):
            ck = pickle.load(open(ckpt_path, 'rb'))
            if ck['n_points'] != self.n or abs(ck['dt'] - dt) > 1e-15: raise RuntimeError('checkpoint does not match this grid')
            Yv, R, T, tau, rows, snaps, A0, nsnap, n, self.mode, proj_log, elapsed0 = (ck[k] for k in
                ['Yv', 'R', 'T', 'tau', 'rows', 'snaps', 'A0', 'nsnap', 'n', 'mode', 'proj_log', 'elapsed'])
            resumed_from = dict(T=T, H0tau=tau, elapsed_s=elapsed0, n_resumes=ck.get('n_resumes', 0) + 1)
            self.n_resumes = resumed_from['n_resumes']
            self.log('resumed from %s at T=%.4f H0tau=%.4f (mode %s, %d records)' % (ckpt_path, T, tau, self.mode, len(rows)))
        elif self.p['project'] is not None:
            Yv, proj_log = self.project(Yv, R, -abs(self.p['project']))
            self.log('Remedy A projection: %s' % proj_log[-1])
        last_ckpt = time.time()
        def save_ckpt():
            if not ckpt_path: return
            tmp = ckpt_path + '.tmp'
            pickle.dump(dict(Yv=Yv, R=R, T=T, tau=tau, rows=rows, snaps=snaps, A0=A0, nsnap=nsnap, n=n, mode=self.mode, proj_log=proj_log,
                             elapsed=elapsed0 + time.time() - t_start, n_points=self.n, dt=dt, n_resumes=getattr(self, 'n_resumes', 0),
                             params=self.p), open(tmp, 'wb'))
            os.replace(tmp, ckpt_path)
        def record(T, Yv, R, tau, want_c):
            F = self.fields(T, Yv, R)
            Bb, Ab, phb, ATb, phTb = F['B'][-1], F['A'][-1], F['phi'][-1], F['A_T'][-1], F['phi_T'][-1]
            eB = math.exp(Bb); H = ATb/eB; v = phTb/eB; t = self.ten
            sg, s1 = t.s(phb), t.s1(phb)
            Wy = H*H - (sg + R)**2/36 - v*v/12 + (s1 + self.Y*v)**2/48 - S.U(phb)/6
            vac = sg**2/36 - s1**2/48 + S.U(phb)/6
            rad = sg*R/18 + R*R/36
            kin = v*v/12; fric = -(2*s1*self.Y*v + self.Y**2*v*v)/48
            row = dict(T=T, H0tau=tau, phi_b=phb, H_over_H0=H*self.rb, lapse=eB/self.rb, A_b=Ab, v_over_H0=v*self.rb, R=R,
                       Wy=Wy*self.rb**2, vac=vac*self.rb**2, rad=rad*self.rb**2, kin=kin*self.rb**2, fric=fric*self.rb**2,
                       phi_min=float(F['phi'].min()), phi_max=float(F['phi'].max()), B_b=Bb)
            if want_c:
                eBb = eB; RT = -4*ATb*R + self.Y*phTb**2/eBb
                if self.mode == 'dev': ra, rb_, rfv = self.accel_dev(F); phiTTb = rfv[-1] + F['rf']['phi_TT'][-1]
                else: ra, rb_, rfv = self.accel(F); phiTTb = rfv[-1]
                H_, M_, sH, sM = self.constraints(F, Yv[3], R, RT, phiTTb)
                hr = np.abs(H_)/sH; mr = np.abs(M_)/sM; ii = self.inner
                near = self.z > -1.0
                iH = int(np.argmax(hr[ii])) + 6
                zi = self.z[ii]; nm = near[ii]
                row.update(Hnear_z=float(zi[nm][int(np.argmax(hr[ii][nm]))]), Mnear_z=float(zi[nm][int(np.argmax(mr[ii][nm]))]))
                row.update(Hmax=float(hr[ii].max()), Mmax=float(mr[ii].max()), Hmax_z=float(self.z[iH]),
                           Hmax_near=float(hr[ii][near[ii]].max()), Mmax_near=float(mr[ii][near[ii]].max()),
                           Hshell=float(hr[-6:].max()), Mshell=float(mr[-6:].max()))
            return row, F
        last_log = time.time()
        while True:
            if n % every == 0 and time.time() - last_ckpt > ckpt_every:
                save_ckpt(); last_ckpt = time.time()
            if n % every == 0:
                want_c = (n % (every*constraint_every) == 0)
                row, F = record(T, Yv, R, tau, want_c)
                if A0 is None: A0 = row['A_b']
                row['ln_a'] = row['A_b'] - A0
                row['mode'] = self.mode
                rows.append(row)
                if T >= nsnap*snap_dT - 1e-12:
                    snaps['%.3f' % T] = dict(A=F['A'][::4].copy(), B=F['B'][::4].copy(), phi=F['phi'][::4].copy(), A_T=F['A_T'][::4].copy())
                    nsnap += 1
                if time.time() - last_log > 30 or n == 0:
                    self.log('T=%7.3f H0tau=%7.4f phi_b=%+.5f H/H0=%+.5f lapse=%.4f R=%.3e Wy/H0^2=%+.4e rad=%.4e mode=%s Hmax=%s [%.0fs]' % (
                        T, tau, row['phi_b'], row['H_over_H0'], row['lapse'], R, row['Wy'], row['rad'], self.mode,
                        ('%.1e' % row['Hmax']) if 'Hmax' in row else '-', time.time() - t_start))
                    last_log = time.time()
                if not all(np.isfinite([row['phi_b'], row['H_over_H0'], R])) or not np.all(np.isfinite(Yv)):
                    stop = 'non-finite'; break
                if row['H_over_H0'] < stop_H_over_H0[0] or row['H_over_H0'] > stop_H_over_H0[1]:
                    stop = 'H/H0 left [%g,%g]' % stop_H_over_H0; break
                if row['lapse'] < 1e-4 or row['lapse'] > 1e4:
                    stop = 'shell lapse out of [1e-4, 1e4] (chart failure)'; break
                if time.time() - t_start > wall_limit:
                    stop = 'wall-clock limit %.0fs' % wall_limit
                    rows.pop()                       # the record at this T is re-made on resume
                    save_ckpt(); break
                if stop_dbdt is not None and len(rows) > 5:
                    dbdt = (rows[-1]['B_b'] - rows[-3]['B_b'])/(rows[-1]['T'] - rows[-3]['T'])
                    if dbdt < stop_dbdt:
                        self.xc_suggest = round(rows[-2]['B_b'] - math.log(self.rb) + rows[-2]['T'], 1)
                        stop = 'd b_b/dT < %g: suggested xc = %.1f' % (stop_dbdt, self.xc_suggest); break
            if T >= T_final - 1e-12: break
            # switch deviation -> full fields
            if self.mode == 'dev' and T + dt > self.T_switch:
                rf = self.ref(T)
                Yv = np.array([rf['A'] + Yv[0], rf['B'] + Yv[1], rf['phi'] + Yv[2], rf['A_T'] + Yv[3], rf['B_T'] + Yv[4], rf['phi_T'] + Yv[5]])
                self.mode = 'full'
                self.log('switched to full-field mode at T=%.4f' % T)
            h = min(dt, T_final - T)
            Bb0 = self.shell_B(T, Yv)
            try:
                k1, r1, F1 = self.rhs(T, Yv, R)
                k2, r2, _ = self.rhs(T + h/2, Yv + h/2*k1, R + h/2*r1)
                k3, r3, _ = self.rhs(T + h/2, Yv + h/2*k2, R + h/2*r2)
                k4, r4, _ = self.rhs(T + h, Yv + h*k3, R + h*r3)
            except (OverflowError, FloatingPointError, ValueError) as err:
                stop = 'step failed at T=%.4f (%s: %s)' % (T, type(err).__name__, err); break
            Yv = Yv + h/6*(k1 + 2*k2 + 2*k3 + k4); R = R + h/6*(r1 + 2*r2 + 2*r3 + r4)
            T = T + h
            if self.mode == 'full':
                rn = self.ref_nodes(T, np.arange(self.nfix))
                for k, nm in enumerate(['A', 'B', 'phi', 'A_T', 'B_T', 'phi_T']): Yv[k, :self.nfix] = rn[nm]
            Bb1 = self.shell_B(T, Yv)
            if not (abs(Bb1) < 40 and np.isfinite(Bb1)):
                stop = 'shell lapse diverged or vanished (|B_b| > 40) at T=%.4f' % T; n += 1; break
            tau += h*0.5*(math.exp(Bb0) + math.exp(Bb1))/self.rb      # trapezoid in e^{B_b}; O(h^2), h ~ 5e-4
            n += 1
        # ------------------------------------------------------------ output
        cols = sorted(set().union(*[r.keys() for r in rows]) - {'mode'})
        arr = np.array([[r.get(c, np.nan) for c in cols] for r in rows], dtype=float)
        modes = np.array([r['mode'] for r in rows])
        os.makedirs(out_dir, exist_ok=True)
        np.savez_compressed(os.path.join(out_dir, tag + '_timeseries.npz'), rec=arr, cols=np.array(cols), modes=modes, z=self.z[::4],
                            **{'snap_%s_%s' % (k, nm): v for k, d_ in snaps.items() for nm, v in d_.items()})
        last = rows[-1]
        summ = dict(tag=tag, params=self.p, xc_suggest=getattr(self, 'xc_suggest', None), n_points=self.n, dt=dt, stop_reason=stop, T_end=T, H0tau_end=tau,
                    runtime_s=round(elapsed0 + time.time() - t_start, 1), resumed_from=resumed_from, points_per_wall=float((1/self.rb)/self.zp[-1]),
                    rho_b=self.rb, phi_b_static=self.bg.phb, seed_phi_b=self.seed.phb, seed_rho_b=self.seed.rhob,
                    T_switch=self.T_switch, F_inf=self.chart.F_inf,
                    end={k: v for k, v in last.items() if k != 'mode'}, projection_log=proj_log,
                    code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    static_w_sha256=hashlib.sha256((HERE/'static_w.py').read_bytes()).hexdigest())
        json.dump(summ, open(os.path.join(out_dir, tag + '_summary.json'), 'w'), indent=1, default=float)
        self.log(json.dumps(dict(stop=stop, T_end=T, H0tau_end=tau, runtime=summ['runtime_s'])))
        return summ, rows

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--delta', type=float, default=0.1); ap.add_argument('--d', type=float, default=0.0)
    ap.add_argument('--dstar', action='store_true', help='use d = d*(delta) from M8_QUADRATIC_TENSION_TUNING.json (MODEL CHANGE)')
    ap.add_argument('--dc', type=float, default=1e-2); ap.add_argument('--Y', type=float, default=0.0)
    ap.add_argument('--dzf', type=float, default=1e-3); ap.add_argument('--dzc', type=float, default=4e-3)
    ap.add_argument('--zfine', type=float, default=0.4); ap.add_argument('--L', type=float, default=16.0)
    ap.add_argument('--xc', type=float, default=math.inf); ap.add_argument('--Tswitch', type=float, default=None)
    ap.add_argument('--kappa', type=float, default=0.0); ap.add_argument('--project', type=float, default=None)
    ap.add_argument('--Tf', type=float, default=12.0); ap.add_argument('--cfl', type=float, default=0.5)
    ap.add_argument('--tag', default='run'); ap.add_argument('--out', default='runs'); ap.add_argument('--wall', type=float, default=1500.0)
    ap.add_argument('--chart', default='bounded', choices=['bounded', 'softplus', 'asinh'])
    ap.add_argument('--prerun', action='store_true', help='old chart, stop when d b_b/dT < -0.9 and report the suggested xc')
    ap.add_argument('--xc_from', default=None, help='take xc from the xc_suggest of this pre-run summary json')
    ap.add_argument('--table_dx', type=float, default=2.5e-5); ap.add_argument('--damp_off', type=float, nargs=2, default=[0.02, 0.05])
    ap.add_argument('--zwidth', type=float, default=0.02); ap.add_argument('--dzm', type=float, default=None)
    ap.add_argument('--zmid', type=float, default=None); ap.add_argument('--zwidth_mid', type=float, default=None)
    ap.add_argument('--grid', default='tanh', choices=['tanh', 'rho']); ap.add_argument('--resume', action='store_true')
    ap.add_argument('--dzcap', type=float, default=None); ap.add_argument('--zcap', type=float, default=3.0)
    ap.add_argument('--smax', type=float, default=20.0, help='bounded chart saturation s_max (B2 dated note 3; default 20 = A1)'); ap.add_argument('--no_ckpt', action='store_true'); ap.add_argument('--stop_dbdt', type=float, default=-0.9)
    a = ap.parse_args()
    d = a.d
    if a.dstar:
        M8 = json.load(open('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json'))
        d = M8['d_star'][{0.1: '0.1', 0.01: '0.01', 0.001: '0.001'}[a.delta]]
    sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
    xc = a.xc
    if a.xc_from:
        xs = json.load(open(a.xc_from))['xc_suggest']
        xc = math.inf if xs is None else xs
    if a.prerun: xc = math.inf
    run = Run(delta=a.delta, d=d, dc=a.dc, Y=a.Y, dzf=a.dzf, dzc=a.dzc, zfine=a.zfine, width=a.zwidth, L=a.L, xc=xc, T_switch=a.Tswitch,
              kappa=a.kappa, damp_off=tuple(a.damp_off), project=a.project, cfl=a.cfl, log=lambda s: print(s, flush=True), chart_kind=a.chart,
              table_dx=a.table_dx, dzm=a.dzm, zmid=a.zmid, width_mid=a.zwidth_mid, grid_kind=a.grid, dzcap=a.dzcap, zcap=a.zcap, s_max=a.smax)
    od = os.path.join(str(HERE), a.out)
    os.makedirs(od, exist_ok=True)
    run.evolve(T_final=a.Tf, tag=a.tag, out_dir=od, wall_limit=a.wall, stop_dbdt=(a.stop_dbdt if a.prerun else None),
               ckpt_path=(None if a.no_ckpt else os.path.join(od, a.tag + '_ckpt.pkl')), resume=a.resume)
