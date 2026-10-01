#!/usr/bin/env python3
"""Static O(4,1)-symmetric shell of the (possibly tuned) HDBLAST model, written in the invariant w = -U_K V_K.

Model (kappa_5^2 = 1): W = 1 - phi + phi^3/3, U = W'^2/2 - (2/3) W^2,
    sigma(phi) = 2W + delta (1 + c phi + d phi^2/2)      (d = 0: registered model; d != 0: MODEL CHANGE, A1)
Static solution in the conformal chart:  B = ln rho(z), A = ln rho(z) + t, phi = phi(z), shell at z = 0 (bulk z < 0).
With x = ln|w|, w = e^{2(z - c0)} > 0 outside the vertex light cone and w < 0 inside it (Milne region), the fields are
    rho^2 = w e^{l(w)},  z = c0 + (1/2) ln w,
and (l, lam = w l', phi, psi = w phi') solve (exact derivation: S0_SYMBOLIC_CHECKS.json, item 3)
    l_x = lam,  lam_x = -3 lam - 1.5 lam^2 - w e^l U(phi)/3,  phi_x = psi,  psi_x = w e^l U'(phi)/4 - 1.5 psi - 1.5 lam psi,
with the regular vertex behaviour l = l1 w + ..., phi = phi_h + p1 w + ... and the first-order constraint
    2 lam + lam^2 = psi^2/3 - w e^l U(phi)/6.
Junctions (pure tension): e^{-B} A_z = sigma/6 and e^{-B} phi_z = -sigma'/2, i.e. (1 + lam)/rho = sigma/6, 2 psi/rho = -sigma'/2.

Chart map (same F on both null coordinates, as in PROPER_CLOCK_PROOF.md):  U = F(u), V = F(v),
    F(x) = ln(1 + e^{-xc}) - ln(e^{-x} + e^{-xc})    (F(0) = 0, F' = 1/(1 + e^{x - xc}) > 0, F -> x for x << xc).
Kruskal-type labels of the static solution: U_K = -e^{-(u + c0)} = e^{-c0}(e^{-xc} - a e^{-U}),  V_K = e^{v - c0} = e^{-c0}/(a e^{-V} - e^{-xc}),
a = 1 + e^{-xc}.  U_K is analytic across U_K = 0 (the future light cone of the vertex), so this chart continues through
the surface where the old conformal chart (xc = infinity) freezes.  Reference fields at (T, Z) (U = T - Z, V = T + Z):
    A = l(w)/2 + ln V_K,  B = l(w)/2 + ln(U_K' V_K')/2,  phi = Phi(w),  w = -U_K V_K.
Floating point (DOP853, rtol 1e-13); tables with cubic Hermite interpolation; series for |w| < w_series.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

C_REG = 2/1.0357712571566784 - 4/3

def W(p): return 1 - p + p*p*p/3
def W1(p): return p*p - 1
# U = W'^2/2 - (2/3) W^2 = -1/6 + 4/3 p - 5/3 p^2 - 4/9 p^3 + 17/18 p^4 + 0 p^5 - 2/27 p^6   (Horner forms; checked in t2)
def U(p): return ((((((-2/27)*p + 0.0)*p + 17/18)*p - 4/9)*p - 5/3)*p + 4/3)*p - 1/6
def U1(p): return (((((-12/27)*p + 0.0)*p + 68/18)*p - 12/9)*p - 10/3)*p + 4/3
def U_direct(p): return 0.5*(p*p - 1)**2 - (2/3)*W(p)**2
def U1_direct(p): return 2*p*(p*p - 1) - (4/3)*W(p)*(p*p - 1)
def U2(p):
    return 2*(p*p - 1) + 4*p*p - (4/3)*((p*p - 1)**2 + W(p)*2*p)
def U3(p):
    # d/dp of U2: 12 p - (4/3)(4p(p^2-1) + 2p(p^2-1) + 2W)  = 12p - (4/3)(6p(p^2-1) + 2W)
    return 12*p - (4/3)*(6*p*(p*p - 1) + 2*W(p))

class Tension:
    def __init__(self, delta, c, d=0.0): self.t, self.c, self.d = delta, c, d
    def s(self, p): return 2*W(p) + self.t*(1 + self.c*p + self.d*p*p/2)
    def s1(self, p): return 2*W1(p) + self.t*(self.c + self.d*p)
    def s2(self, p): return 4*p + self.t*self.d

def rhs_x(x, y, sgn):
    l, lam, ph, psi = y
    w = sgn*math.exp(x); el = math.exp(l)
    return [lam, -3*lam - 1.5*lam*lam - w*el*U(ph)/3, psi, 0.25*w*el*U1(ph) - 1.5*psi - 1.5*lam*psi]

def rhs_x_vec(w, l, lam, ph, psi):
    el = np.exp(l)
    return (lam, -3*lam - 1.5*lam*lam - w*el*U(ph)/3, psi, 0.25*w*el*U1(ph) - 1.5*psi - 1.5*lam*psi)

def series_coeffs(phh):
    u0, u1, u2, u3 = U(phh), U1(phh), U2(phh), U3(phh)
    l1 = -u0/12; p1 = u1/10
    l2 = (-3*l1*l1 - (2/3)*(l1*u0 + u1*p1))/20
    p2 = (l1*u1 + u2*p1 - 6*l1*p1)/28
    C2 = u1*p2 + 0.5*u2*p1*p1 + l1*u1*p1 + (l2 + 0.5*l1*l1)*u0
    l3 = (-6*l1*l2 - C2/3)/18
    D2 = u2*p2 + 0.5*u3*p1*p1 + l1*u2*p1 + (l2 + 0.5*l1*l1)*u1
    p3 = (0.25*D2 - 3*l1*p2 - 3*l2*p1)/13.5
    return (l1, l2, l3), (p1, p2, p3)

class StaticShell:
    """Solve the static shell for tension ten; build tables of l, lam, phi, psi versus x = ln|w| on both sides of w = 0."""
    def __init__(self, ten, phh_guess, rtol=1e-13, atol=1e-16, w_start=1e-8, bracket=0.3, table_dx=2.5e-4):
        self.ten, self.rtol, self.atol, self.w_start = ten, rtol, atol, w_start
        self.table_dx = table_dx
        eta0 = math.log(phh_guess + 1)
        f = lambda eta: self._shoot(-1 + math.exp(eta))[0]
        lo, hi = eta0 - bracket, eta0 + bracket
        flo, fhi = f(lo), f(hi)
        k = 0
        while flo*fhi > 0 and k < 20:
            lo -= bracket; hi += bracket; flo, fhi = f(lo), f(hi); k += 1
        eta = brentq(f, lo, hi, xtol=1e-15, rtol=1e-15, maxiter=200)
        self.phh = -1 + math.exp(eta)
        r2, sol, xb = self._shoot(self.phh, dense=True)
        self.xb = xb; self.wb = math.exp(xb)
        yb = sol.sol(xb); self.lb, self.lamb, self.phb, self.psib = yb
        self.rhob = math.sqrt(self.wb*math.exp(self.lb))
        self.c0 = -0.5*xb                      # z = c0 + ln(w)/2 and z_b = 0
        self.junction_residuals = ((1 + self.lamb)/self.rhob - ten.s(self.phb)/6, 2*self.psib/self.rhob + ten.s1(self.phb)/2)
        self.constraint_at_shell = 2*self.lamb + self.lamb**2 - (self.psib**2/3 - self.wb*math.exp(self.lb)*U(self.phb)/6)
        self.lc, self.pc = series_coeffs(self.phh)
        self._tables(dx=self.table_dx)

    def _y0(self, phh, sgn):
        (l1, l2, l3), (p1, p2, p3) = series_coeffs(phh)
        w = sgn*self.w_start
        return [l1*w + l2*w*w + l3*w**3, l1*w + 2*l2*w*w + 3*l3*w**3, phh + p1*w + p2*w*w + p3*w**3, p1*w + 2*p2*w*w + 3*p3*w**3]

    def _shoot(self, phh, dense=False):
        ten = self.ten
        def ev(x, y, sgn):
            l, lam, ph, psi = y; rho = math.sqrt(math.exp(x + l))
            return (1 + lam) - ten.s(ph)*rho/6
        ev.terminal = True; ev.direction = -1
        sol = solve_ivp(rhs_x, (math.log(self.w_start), 12.0), self._y0(phh, 1), args=(1,), method='DOP853',
                        rtol=self.rtol, atol=self.atol, events=ev, dense_output=dense)
        if len(sol.t_events[0]) == 0:
            return (1e3, sol, None) if not dense else (1e3, sol, None)
        xb = sol.t_events[0][0]; l, lam, ph, psi = sol.y_events[0][0]
        rho = math.sqrt(math.exp(xb + l))
        return 2*psi/rho + ten.s1(ph)/2, sol, xb

    def _tables(self, dx=2.5e-4, x_int_max=None):
        x0 = math.log(self.w_start)
        self.x_lo = x0
        # exterior: from the vertex to slightly beyond the shell (w up to 1.05 w_b)
        xe = np.arange(x0, self.xb + 0.05 + dx, dx)
        se = solve_ivp(rhs_x, (x0, xe[-1]), self._y0(self.phh, 1), args=(1,), method='DOP853', rtol=self.rtol, atol=self.atol, t_eval=xe)
        # interior (Milne region, w < 0): up to |w| = 25 w_b (the interior recollapses; its crunch is at |w| ~ 6e3 w_b)
        xi_max = self.xb + math.log(25.0) if x_int_max is None else x_int_max
        xi = np.arange(x0, xi_max + dx, dx)
        si = solve_ivp(rhs_x, (x0, xi[-1]), self._y0(self.phh, -1), args=(-1,), method='DOP853', rtol=self.rtol, atol=self.atol, t_eval=xi)
        self.tab = {}
        for sgn, s in [(1, se), (-1, si)]:
            x = s.t; Y = s.y; w = sgn*np.exp(x)
            d = np.array(rhs_x_vec(w, *Y))
            self.tab[sgn] = dict(x=x, Y=Y, dY=d, dx=dx, x0=x[0], n=len(x))
        self.x_hi = {1: se.t[-1], -1: si.t[-1]}

    # ---- evaluation of l(w), l'(w), l''(w), Phi, Phi', Phi'' (arrays)
    def fields_w(self, w):
        """returns l, lhat=l'(w), lhat_w=l''(w), phi, phat=Phi'(w), phat_w=Phi''(w) for an array w"""
        w = np.asarray(w, float)
        out = [np.empty_like(w) for _ in range(6)]
        aw = np.abs(w)
        ser = aw < self.w_start*1.0000001
        (l1, l2, l3), (p1, p2, p3) = self.lc, self.pc
        if ser.any():
            ws = w[ser]
            out[0][ser] = l1*ws + l2*ws**2 + l3*ws**3; out[1][ser] = l1 + 2*l2*ws + 3*l3*ws**2; out[2][ser] = 2*l2 + 6*l3*ws
            out[3][ser] = self.phh + p1*ws + p2*ws**2 + p3*ws**3; out[4][ser] = p1 + 2*p2*ws + 3*p3*ws**2; out[5][ser] = 2*p2 + 6*p3*ws
        for sgn in (1, -1):
            m = (~ser) & ((w > 0) if sgn == 1 else (w < 0))
            if not m.any(): continue
            tb = self.tab[sgn]; x = np.log(aw[m])
            if np.any(x > tb['x'][-1] + 1e-12):
                raise ValueError('w outside the static table (sgn=%d, max x=%.4f > %.4f)' % (sgn, x.max(), tb['x'][-1]))
            i = np.clip(((x - tb['x0'])/tb['dx']).astype(int), 0, tb['n'] - 2)
            h = tb['x'][i + 1] - tb['x'][i]; s = (x - tb['x'][i])/h
            h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s; h01 = -2*s**3 + 3*s**2; h11 = s**3 - s**2
            Y = [h00*tb['Y'][k][i] + h10*h*tb['dY'][k][i] + h01*tb['Y'][k][i + 1] + h11*h*tb['dY'][k][i + 1] for k in range(4)]
            wm = w[m]
            l, lam, ph, psi = Y
            dl, dlam, dph, dpsi = rhs_x_vec(wm, l, lam, ph, psi)
            out[0][m] = l; out[1][m] = lam/wm; out[2][m] = (dlam - lam)/wm**2
            out[3][m] = ph; out[4][m] = psi/wm; out[5][m] = (dpsi - psi)/wm**2
        return out

    def profile_z(self, z):
        """exterior static profile at conformal z <= 0 (for controls): rho, h = rho_z/rho, phi, phi_z"""
        w = np.exp(2*(np.asarray(z) - self.c0))
        l, lh, lhw, ph, pht, phw = self.fields_w(w)
        rho = np.sqrt(w*np.exp(l))
        return rho, 1 + w*lh, ph, 2*w*pht


class Chart:
    """Common null-label map U = F(u), V = F(v).
    kind='asinh' (default since the 28 Sept dated note in REGISTRATION.md):
        F(x) = C - arcsinh(e^{-x}/(2b)),  C = arcsinh(1/(2b)),  b = e^{-xc};  F(0) = 0,  F' = 1/sqrt(1 + 4 b^2 e^{2x}).
        In terms of U_K = -e^{-(u + c0)}:  U = C - arcsinh(-U_K e^{c0}/(2b)), an entire function of U_K: the chart is regular
        across the vertex light cone (U_K = 0, U = C) and for all U_K > 0.
    kind='bounded' (default since the second 28 Sept dated note): s = e^{-u} = h(sigma), sigma = 2 b sinh(C - U),
        C = arcsinh(1/(2b)), h(sigma) = s_max sigma/(s_max + eps ln(1 + e^{-sigma/eps})), s_max = 20, eps = b.
        h = sigma for sigma >> eps (so the map equals 'asinh' until the light cone), h is smooth and monotone across it,
        and h -> -s_max + s_max^2/|sigma| beyond: U -> infinity at U_K = s_max e^{-c0} with U_K' ~ e^{-U}, so every static
        (V < 0) grid point keeps |w| <= ~s_max w_b (away from the crunch of the light-cone interior at |w| ~ 6e3 w_b) and
        the far region stays smooth (B ~ -U/2 there).
    kind='softplus' (first version): F(x) = ln(1 + b) - ln(e^{-x} + b); singular at U_K = b e^{-c0} > 0 (found by the
        registered-c calibration run, which reached it at H0 tau = 6.59).
    xc = inf gives the identity map (the pilot's chart)."""
    def __init__(self, xc, kind='bounded', s_max=20.0, eps=None):
        self.xc, self.kind, self.s_max = xc, kind, s_max
        self.b = 0.0 if not math.isfinite(xc) else math.exp(-xc)
        self.eps = self.b if eps is None else eps
        if self.b == 0.0:
            self.kind = 'identity'; self.F_inf = math.inf; self.C = math.inf
        elif kind == 'bounded':
            self.C = math.asinh(1/(2*self.b)); self.F_inf = self.C
        elif kind == 'softplus':
            self.a = 1 + self.b; self.F_inf = math.log(self.a) - math.log(self.b); self.C = self.F_inf
        else:
            self.C = math.asinh(1/(2*self.b)); self.F_inf = self.C
    def h_of_sigma(self, sig):
        """h(sigma) = s_max sigma / D,  D = s_max + eps ln(1 + e^{-sigma/eps}):  h = sigma (to ~e^{-sigma/eps}) for sigma >> eps,
        h -> -s_max + s_max^2/|sigma| for sigma -> -infinity; smooth and monotone.  Returns h and its first three derivatives."""
        sm, eps = self.s_max, self.eps
        x = -sig/eps
        L = np.where(x > 0, x + np.log1p(np.exp(-np.abs(x))), np.log1p(np.exp(-np.abs(x))))
        pp = np.where(x > 0, 1/(1 + np.exp(-np.abs(x))), np.exp(-np.abs(x))/(1 + np.exp(-np.abs(x))))
        D = sm + eps*L; D1 = -pp; D2 = pp*(1 - pp)/eps; D3 = -pp*(1 - pp)*(1 - 2*pp)/eps**2
        h0 = sm*sig/D
        h1 = sm*(D - sig*D1)/D**2
        h2 = sm*(-sig*D2/D**2 - 2*(D - sig*D1)*D1/D**3)
        h3 = sm*(-D2/D**2 - sig*D3/D**2 + 4*sig*D1*D2/D**3 - 2*(D - sig*D1)*D2/D**3 + 6*(D - sig*D1)*D1**2/D**4)
        return h0, h1, h2, h3
    def s_of_U(self, U):
        """bounded map: s = e^{-u} = h(sigma(U)), sigma = 2 b sinh(C - U); returns s and its first three U-derivatives"""
        b, C = self.b, self.C
        x = C - U; sig = 2*b*np.sinh(x); sig1 = -2*b*np.cosh(x); sig2 = sig; sig3 = sig1
        h0, h1, h2, h3 = self.h_of_sigma(sig)
        s1 = h1*sig1
        s2 = h2*sig1**2 + h1*sig2
        s3 = h3*sig1**3 + 3*h2*sig1*sig2 + h1*sig3
        return h0, s1, s2, s3
    def F(self, x):
        x = np.asarray(x, float)
        if self.kind == 'identity': return x
        if self.kind == 'bounded':
            hv = np.exp(-x); sig = hv.copy()                          # invert h by Newton (h' > 0)
            for _ in range(60):
                h0, h1, _, _ = self.h_of_sigma(sig); sig = sig - (h0 - hv)/h1
            return self.C - np.arcsinh(sig/(2*self.b))
        if self.kind == 'softplus': return math.log(self.a) - np.log(np.exp(-x) + self.b)
        return self.C - np.arcsinh(np.exp(-x)/(2*self.b))
    def Fp(self, x):
        x = np.asarray(x, float)
        if self.kind == 'identity': return np.ones_like(x)
        if self.kind == 'bounded':
            U = self.F(x); s0, s1, _, _ = self.s_of_U(U)
            return -s0/s1                      # dU/du = -(ds/du)/(ds/dU) ... with ds/du = -s
        if self.kind == 'softplus': return 1/(1 + self.b*np.exp(x))
        return 1/np.sqrt(1 + 4*self.b**2*np.exp(2*x))

def _lnsinh(x):
    """ln sinh(x) for x > 0, stable for large x"""
    return x + np.log1p(-np.exp(-2*x)) - math.log(2)

def maps(chart, c0, Uc, Vc):
    """U_K, U_K', U_K'' and (ln U_K')' , (ln U_K')'' ; ln V_K, (ln V_K)', (ln V_K)'', ln V_K', (ln V_K')', (ln V_K')''"""
    ec = math.exp(-c0)
    if chart.kind == 'identity':
        eU = np.exp(-Uc)
        UK = -ec*eU; UK1 = ec*eU; UK2 = -UK1
        lU1p = -np.ones_like(Uc); lU1pp = np.zeros_like(Uc)
        lnVK = Vc - c0; r1 = np.ones_like(Vc); r2 = np.zeros_like(Vc)
        lnVK1 = lnVK; q1 = np.ones_like(Vc); q2 = np.zeros_like(Vc)
    elif chart.kind == 'bounded':
        s0, s1, s2, s3 = chart.s_of_U(Uc)
        UK = -ec*s0; UK1 = -ec*s1; UK2 = -ec*s2
        lU1p = s2/s1; lU1pp = s3/s1 - (s2/s1)**2
        v0, v1, v2, v3 = chart.s_of_U(Vc)
        if np.any(v0 <= 0): raise ValueError('reference evaluated at V >= F_inf')
        lnVK = -c0 - np.log(v0)
        r1 = -v1/v0; r2 = -v2/v0 + (v1/v0)**2
        lnVK1 = lnVK + np.log(r1)
        q1 = v2/v1 - 2*v1/v0; q2 = (v3/v1 - (v2/v1)**2) - 2*(v2/v0 - (v1/v0)**2)
    elif chart.kind == 'softplus':
        a, b = chart.a, chart.b
        eU = np.exp(-Uc)
        UK = ec*(b - a*eU); UK1 = ec*a*eU; UK2 = -UK1
        lU1p = -np.ones_like(Uc); lU1pp = np.zeros_like(Uc)
        eV = np.exp(-Vc); D = a*eV - b
        if np.any(D <= 0): raise ValueError('reference evaluated at V >= F_inf')
        lnVK = -c0 - np.log(D); r1 = a*eV/D; r2 = a*b*eV/D**2
        lnVK1 = lnVK + np.log(r1); q1 = -1 + 2*a*eV/D; q2 = 2*a*b*eV/D**2
    else:
        b, C = chart.b, chart.C
        xu = C - Uc
        UK = -2*b*ec*np.sinh(xu); UK1 = 2*b*ec*np.cosh(xu); UK2 = UK
        lU1p = -np.tanh(xu); lU1pp = 1/np.cosh(xu)**2
        xv = C - Vc
        if np.any(xv <= 0): raise ValueError('reference evaluated at V >= F_inf')
        lnVK = -c0 - math.log(2*b) - _lnsinh(xv)
        cth = 1/np.tanh(xv); csch2 = 1/np.sinh(xv)**2
        r1 = cth; r2 = csch2
        lnVK1 = lnVK + np.log(cth)
        s2x = np.sinh(2*xv)
        q1 = cth + 2/s2x; q2 = csch2 + 4*np.cosh(2*xv)/s2x**2
    return UK, UK1, UK2, lU1p, lU1pp, lnVK, r1, r2, lnVK1, q1, q2

def reference(shell, chart, T, Z, want2=True):
    """Static solution `shell` evaluated at chart points (T, Z) (arrays).  Returns dict of A, B, phi and their first
    (T, Z) derivatives and, if want2, second derivatives (ZZ, TZ, TT).  Valid for V = T + Z < F_inf (V_K finite)."""
    Uc = T - Z; Vc = T + Z
    UK, UK1, UK2, lU1p, lU1pp, lnVK, r1, r2, lnVK1, q1, q2 = maps(chart, shell.c0, Uc, Vc)
    VK = np.exp(lnVK); VK1 = VK*r1; VK2 = VK*(r2 + r1*r1)
    w = -UK*VK
    l, lh, lhw, ph, pht, phw = shell.fields_w(w)
    wU = -UK1*VK; wV = -UK*VK1
    lU = lh*wU; lV = lh*wV
    pU = pht*wU; pV = pht*wV
    A = 0.5*l + lnVK
    Bv = 0.5*l + 0.5*np.log(UK1) + 0.5*lnVK1
    AU = 0.5*lU; AV = 0.5*lV + r1
    BU = 0.5*lU + 0.5*lU1p; BV = 0.5*lV + 0.5*q1
    out = dict(A=A, B=Bv, phi=ph, w=w,
               A_T=AU + AV, A_Z=AV - AU, B_T=BU + BV, B_Z=BV - BU, phi_T=pU + pV, phi_Z=pV - pU)
    if want2:
        wUU = -UK2*VK; wVV = -UK*VK2; wUV = -UK1*VK1
        lUU = lhw*wU*wU + lh*wUU; lVV = lhw*wV*wV + lh*wVV; lUV = lhw*wU*wV + lh*wUV
        pUU = phw*wU*wU + pht*wUU; pVV = phw*wV*wV + pht*wVV; pUV = phw*wU*wV + pht*wUV
        AUU = 0.5*lUU; AVV = 0.5*lVV + r2; AUV = 0.5*lUV
        BUU = 0.5*lUU + 0.5*lU1pp; BVV = 0.5*lVV + 0.5*q2; BUV = 0.5*lUV
        for nm, (fUU, fVV, fUV) in dict(A=(AUU, AVV, AUV), B=(BUU, BVV, BUV), phi=(pUU, pVV, pUV)).items():
            out[nm + '_ZZ'] = fUU - 2*fUV + fVV
            out[nm + '_TZ'] = fVV - fUU
            out[nm + '_TT'] = fUU + 2*fUV + fVV
    return out
