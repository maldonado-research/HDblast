"""Numerical core: backgrounds and the gauge-invariant scalar/tensor boundary problems.

Everything is written in the shifted scalar e = phi - v (v = +1 or -1 is the vacuum at the
regular cone), so that cone displacements as small as 1e-28 are represented exactly.

Background (kappa_5^2 = 1):   rho' = sqrt(1 + rho^2 (e'^2/12 - U/6)),   e'' = U_phi - 4 (rho'/rho) e'
Shell (Z2, bulk y<y_b):       rho'/rho = sigma/6 ,  e' = -sigma'/2 ,  sigma = 2W + delta (1 + c phi)

Scalar sector (derived and checked in derive_linearized.py, DERIVATION_RESULTS.json):
    Z' = -(2H + g) Z - X/3,     X' = g X + (3 lam/rho^2 - 2 e'^2) Z,     g = phi''/phi',  lam = mu2 + 4
    shell:  M(mu2) = (g + sigma''/2) X + 3 lam Z/rho^2 = 0
    X = gauge-invariant scalar perturbation, Psi = phi' Z = gauge-invariant curvature perturbation.
    Regular cone data:  X ~ y^s,  Z ~ -y^(s+1)/(3(s+4)),  s = -3/2 + sqrt(9/4 - mu2).
Tensor sector:  h'' + 4H h' + mu2 h/rho^2 = 0,   h'(y_b) = 0,  h ~ y^s at the cone.
"""
import math
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root

C_REG = 2/1.0357712571566784 - 4/3

def _poly_coeffs(v):
    e, f = sp.symbols('e f')
    W = 1 - f + f**3/3
    U = sp.expand(((sp.diff(W, f))**2/2 - sp.Rational(2, 3)*W**2).subs(f, v + e))
    P = sp.Poly(U, e)
    c = [sp.Rational(P.coeff_monomial(e**k)) for k in range(7)]
    return c

class Potential:
    """U(v+e) and derivatives from exact rational Taylor coefficients (Horner)."""
    def __init__(self, v):
        self.v = v
        self.c = _poly_coeffs(v)
        cf = [float(x) for x in self.c]
        self.c0 = np.array(cf)
        self.c1 = np.array([k*cf[k] for k in range(1, 7)])
        self.c2 = np.array([k*(k-1)*cf[k] for k in range(2, 7)])
    @staticmethod
    def _h(c, e):
        r = 0*e + c[-1]
        for a in c[-2::-1]:
            r = r*e + a
        return r
    def U(self, e): return self._h(self.c0, e)
    def U1(self, e): return self._h(self.c1, e)
    def U2(self, e): return self._h(self.c2, e)
    def U3(self, e): return self._h(np.array([k*(k-1)*(k-2)*float(self.c[k]) for k in range(3, 7)]), e)

def sigma_fns(v, delta, c):
    """sigma, sigma', sigma'' as functions of e (phi = v+e)."""
    def W(e): f = v + e; return 1 - f + f**3/3
    s0 = lambda e: 2*W(e) + delta*(1 + c*(v + e))
    s1 = lambda e: 2*(2*v*e + e*e) + delta*c          # 2(phi^2-1) + delta c
    s2 = lambda e: 4*(v + e)
    return s0, s1, s2

class Background:
    def __init__(self, v, delta, c=C_REG, e_h=None, y_b=None, y0=1e-4, rtol=1e-12, polish=True):
        self.v, self.delta, self.c, self.y0, self.rtol = v, delta, c, y0, rtol
        self.pot = Potential(v)
        self.s0, self.s1, self.s2 = sigma_fns(v, delta, c)
        self.e_h, self.y_b = e_h, y_b
        self.polish_info = None
        if polish:
            self.polish()
    # regular cone series (checked in the source package's independent review)
    def cone_state(self, e_h, y0):
        P = self.pot; u, up, upp = P.U(e_h), P.U1(e_h), P.U2(e_h)
        aa = -u/36; bb = u*u/4320 - up*up/750; cc = up/10; dd = up*(upp/280 + u/630)
        return np.array([y0 + aa*y0**3 + bb*y0**5, e_h + cc*y0*y0 + dd*y0**4, 2*cc*y0 + 4*dd*y0**3])
    def rhs_bg(self, y, st):
        rho, e, ep = st
        rad = 1 + rho*rho*(ep*ep/12 - self.pot.U(e)/6)
        rp = math.sqrt(rad)
        return [rp, ep, self.pot.U1(e) - 4*rp*ep/rho]
    def integrate(self, e_h, y_b, dense=False, y0=None, rtol=None):
        y0 = self.y0 if y0 is None else y0; rtol = self.rtol if rtol is None else rtol
        st = self.cone_state(e_h, y0)
        sc = abs(st[2]) + 1e-300
        sol = solve_ivp(self.rhs_bg, (y0, y_b), st, method='DOP853', rtol=rtol,
                        atol=[1e-14, 1e-3*abs(e_h)*1e-12 + 1e-300, sc*1e-12], max_step=0.1, dense_output=dense)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol
    def junction_residual(self, e_h, y_b):
        rho, e, ep = self.integrate(e_h, y_b).y[:, -1]
        rp = math.sqrt(1 + rho*rho*(ep*ep/12 - self.pot.U(e)/6))
        return np.array([rp/rho - self.s0(e)/6, ep + self.s1(e)/2])
    def polish(self):
        sgn = math.copysign(1.0, self.e_h)
        f = lambda x: self.junction_residual(sgn*10**x[0], x[1])/max(self.delta, 1e-3)
        before = self.junction_residual(self.e_h, self.y_b)
        fit = root(f, [math.log10(abs(self.e_h)), self.y_b], tol=1e-13, options={'eps': 1e-9})
        eh, yb = sgn*10**fit.x[0], fit.x[1]
        after = self.junction_residual(eh, yb)
        if np.max(np.abs(after)) <= np.max(np.abs(before)):
            self.e_h, self.y_b = eh, yb
        else:
            after = before
        self.polish_info = dict(residual_before=before.tolist(), residual_after=after.tolist(),
                                e_h=self.e_h, y_b=self.y_b)
    def shell_values(self):
        rho, e, ep = self.integrate(self.e_h, self.y_b).y[:, -1]
        P = self.pot
        rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6))
        H = rp/rho
        g = P.U1(e)/ep - 4*H
        return dict(rho_b=rho, phi_b=self.v + e, eta_b=e, phi_y_b=ep, H2_brane=1/rho**2,
                    g_b=g, sigma2_b=self.s2(e), B=g + self.s2(e)/2,
                    metric_junction=rp/rho - self.s0(e)/6, scalar_junction=ep + self.s1(e)/2)

    # ------------------------------------------------------------------ scalar sector
    def _pert_rhs(self, mu2):
        lam = mu2 + 4.0; P = self.pot
        def f(y, st):
            rho, e, ep = st[0].real, st[1].real, st[2].real
            X, Z = st[3], st[4]
            rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6)); H = rp/rho
            g = P.U1(e)/ep - 4*H
            return [rp, ep, P.U1(e) - 4*rp*ep/rho,
                    g*X + (3*lam/(rho*rho) - 2*ep*ep)*Z,
                    -(2*H + g)*Z - X/3.0]
        return f
    def scalar_M(self, mu2, y0=None, rtol=None, full=False):
        """Mismatch M(mu2) for the regular solution normalised X ~ y^s at the cone."""
        y0 = self.y0 if y0 is None else y0; rtol = self.rtol if rtol is None else rtol
        cplx = np.iscomplexobj(mu2) or isinstance(mu2, complex)
        s = -1.5 + np.sqrt(2.25 - mu2 + 0j) if cplx else -1.5 + math.sqrt(2.25 - mu2)
        bg = self.cone_state(self.e_h, y0)
        X0 = y0**s; Z0 = -y0**(s + 1)/(3*(s + 4))
        st = np.array([bg[0], bg[1], bg[2], X0, Z0], dtype=complex if cplx else float)
        sol = solve_ivp(self._pert_rhs(mu2), (y0, self.y_b), st, method='DOP853', rtol=rtol,
                        atol=[1e-14, 1e-300, 1e-300, 1e-300, 1e-300], max_step=0.1)
        if not sol.success:
            raise RuntimeError(sol.message)
        rho, e, ep, X, Z = sol.y[:, -1]
        rho, e, ep = rho.real, e.real, ep.real
        P = self.pot
        rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6)); H = rp/rho
        g = P.U1(e)/ep - 4*H; B = g + self.s2(e)/2; lam = mu2 + 4
        M = B*X + 3*lam*Z/rho**2
        Mn = M/(abs(B*X) + abs(3*lam*Z/rho**2))
        if full:
            return dict(M=M, Mhat=Mn, X_b=X, Z_b=Z, B=B, s=s, nfev=sol.nfev)
        return Mn

    def scalar_XZ_renorm(self, mu2, nseg=60, y0=None, rtol=None):
        """Real mu2 only: integrate the (X, Z) system in segments, rescaling (X, Z) by a positive factor at
        each segment end (allowed: the problem is linear and only the ratio Z_b/X_b enters the sign of M).
        Avoids overflow of y^s for very negative mu2.  Returns (X_b, Z_b) up to a positive factor."""
        y0 = self.y0 if y0 is None else y0; rtol = self.rtol if rtol is None else rtol
        s = -1.5 + math.sqrt(2.25 - mu2)
        bg = self.cone_state(self.e_h, y0)
        st = np.array([bg[0], bg[1], bg[2], 1.0, -y0/(3*(s + 4))])
        edges = np.geomspace(y0, self.y_b, nseg + 1); edges[-1] = self.y_b
        f = self._pert_rhs(mu2)
        for a, b in zip(edges[:-1], edges[1:]):
            sol = solve_ivp(f, (a, b), st, method='DOP853', rtol=rtol,
                            atol=[1e-14, 1e-300, 1e-300, 1e-300, 1e-300], max_step=0.1)
            if not sol.success:
                raise RuntimeError(sol.message)
            st = sol.y[:, -1].copy(); sc = abs(st[3]) + abs(st[4]); st[3:] /= sc
        return st

    # ------------------------------------------------------------------ tensor sector
    def tensor_mismatch(self, mu2, y0=None, rtol=None):
        """h'(y_b)/ (|h'| + |h|/rho_b) for the regular solution h ~ y^s."""
        y0 = self.y0 if y0 is None else y0; rtol = self.rtol if rtol is None else rtol
        s = -1.5 + math.sqrt(2.25 - mu2)
        bg = self.cone_state(self.e_h, y0)
        P = self.pot
        def f(y, st):
            rho, e, ep, h, hp = st
            rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6))
            return [rp, ep, P.U1(e) - 4*rp*ep/rho, hp, -4*rp/rho*hp - mu2*h/rho**2]
        st = [bg[0], bg[1], bg[2], y0**s, s*y0**(s - 1)]
        sol = solve_ivp(f, (y0, self.y_b), st, method='DOP853', rtol=rtol,
                        atol=[1e-14, 1e-300, 1e-300, 1e-300, 1e-300], max_step=0.1)
        rho, e, ep, h, hp = sol.y[:, -1]
        return hp/(abs(hp) + abs(h)/rho)

    # ------------------------------------------------------------------ cross-check: raw longitudinal gauge
    def longitudinal_M(self, mu2, y0=None, rtol=None):
        """Integrate the UNREDUCED longitudinal-gauge equations (psi, chi, chi') taken directly from the
        sympy output (yt constraint + scalar equation, N=-2psi) and evaluate the scalar junction
        chi' + sigma''/2 chi + 2 phi' psi at the shell.  Also returns the Gauss-constraint residual."""
        y0 = self.y0 if y0 is None else y0; rtol = self.rtol if rtol is None else rtol
        s = -1.5 + math.sqrt(2.25 - mu2); lam = mu2 + 4
        bg = self.cone_state(self.e_h, y0); P = self.pot
        rho0, e0, ep0 = bg
        rp0 = math.sqrt(1 + rho0**2*(ep0**2/12 - P.U(e0)/6)); H0 = rp0/rho0; g0 = P.U1(e0)/ep0 - 4*H0
        X0 = y0**s; Z0 = -y0**(s + 1)/(3*(s + 4))
        Xp0 = g0*X0 + (3*lam/rho0**2 - 2*ep0**2)*Z0
        def f(y, st):
            rho, e, ep, psi, chi, chip = st
            rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6)); H = rp/rho
            psip = -2*H*psi - ep*chi/3
            chipp = -4*H*chip - (mu2/rho**2 - P.U2(e))*chi - 4*P.U1(e)*psi - 6*ep*psip
            return [rp, ep, P.U1(e) - 4*rp*ep/rho, psip, chip, chipp]
        st = [rho0, e0, ep0, ep0*Z0, X0, Xp0]
        sol = solve_ivp(f, (y0, self.y_b), st, method='DOP853', rtol=rtol,
                        atol=[1e-14, 1e-300, 1e-300, 1e-300, 1e-300, 1e-300], max_step=0.1)
        rho, e, ep, psi, chi, chip = sol.y[:, -1]
        rp = math.sqrt(1 + rho*rho*(ep*ep/12 - P.U(e)/6)); H = rp/rho
        M = chip + self.s2(e)/2*chi + 2*ep*psi
        Mn = M/(abs(chip) + abs(self.s2(e)/2*chi) + abs(2*ep*psi))
        # Gauss constraint (from sympy): phi' chi' = U1 chi - 4 H phi' chi + (3 lam/rho^2 - 2 phi'^2) psi
        gauss = ep*chip - (P.U1(e)*chi - 4*H*ep*chi + (3*lam/rho**2 - 2*ep*ep)*psi)
        gscale = abs(ep*chip) + abs(P.U1(e)*chi) + abs(4*H*ep*chi) + abs(3*lam/rho**2*psi)
        return Mn, gauss/gscale

def growth_exponent(mu2):
    """Mode ~ e^{p H tau} with p^2 + 3p + mu2 = 0 (largest root)."""
    return -1.5 + math.sqrt(2.25 - mu2)
