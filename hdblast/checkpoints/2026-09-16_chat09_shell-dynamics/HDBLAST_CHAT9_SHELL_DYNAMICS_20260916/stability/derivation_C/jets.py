#!/usr/bin/env python3
"""HDBLAST Chat 9 / derivation C  -  truncated multivariate Taylor ("jet") arithmetic and
brute-force differential geometry.  numpy + stdlib only.

Purpose: machine-precision, formula-free verification of every hand-derived equation in
DERIVATION_C.md.  A jet is the array of Taylor coefficients (NOT derivatives: coefficient of
prod dx_v^e_v, i.e. derivative / prod e_v!) of a function of `nvars` variables around a base point,
truncated at total degree `order`.  Complex coefficients are allowed, so exact first-order
linearisation is obtained by the complex-step trick  d/d(eps) F = Im F(i*1e-30) / 1e-30
(no subtractive cancellation, no truncation error).

Conventions of the geometry routines (MTW):
    Gamma^a_bc = 1/2 g^ad (d_b g_dc + d_c g_db - d_d g_bc)
    Riem^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb + Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb
    Ric_bd     = Riem^a_bad          (unit dS_4 has Ric = +3 gamma, unit sphere S^n has Ric = (n-1) g)
Field equations of the registered model (kappa_5 = 1):
    E_AB  := R_AB - d_A phi d_B phi - (2/3) U(phi) g_AB = 0 ,      E_phi := Box phi - U_phi = 0 .
"""
import itertools
import math
import numpy as np


class JetSpace:
    def __init__(self, nvars, order):
        self.nv, self.N = nvars, order
        monos = [e for e in itertools.product(range(order + 1), repeat=nvars) if sum(e) <= order]
        monos.sort(key=lambda e: (sum(e), tuple(-x for x in e)))
        self.monos = monos
        self.n = n = len(monos)
        self.idx = {e: i for i, e in enumerate(monos)}
        T = np.zeros((n, n, n))
        for i, ei in enumerate(monos):
            for j, ej in enumerate(monos):
                ek = tuple(a + b for a, b in zip(ei, ej))
                if sum(ek) <= order:
                    T[i, j, self.idx[ek]] = 1.0
        self.T = T
        self.D = []
        for v in range(nvars):
            Dm = np.zeros((n, n))
            for i, ei in enumerate(monos):
                if ei[v] > 0:
                    ek = list(ei); ek[v] -= 1
                    Dm[self.idx[tuple(ek)], i] = ei[v]
            self.D.append(Dm)

    # ---- constructors -------------------------------------------------------------------
    def const(self, c):
        a = np.zeros(self.n, dtype=complex); a[0] = c; return a

    def var(self, v, x0=0.0):
        a = np.zeros(self.n, dtype=complex); a[0] = x0
        e = [0] * self.nv; e[v] = 1; a[self.idx[tuple(e)]] = 1.0; return a

    def from_derivs(self, derivs):
        """derivs: dict {exponent tuple: partial derivative value}; converts to Taylor coefficients."""
        a = np.zeros(self.n, dtype=complex)
        for e, val in derivs.items():
            fac = 1.0
            for k in e: fac *= math.factorial(k)
            a[self.idx[tuple(e)]] = val / fac
        return a

    def coef(self, a, e):
        return a[..., self.idx[tuple(e)]]

    def deriv_value(self, a, e):
        fac = 1.0
        for k in e: fac *= math.factorial(k)
        return a[..., self.idx[tuple(e)]] * fac

    # ---- algebra ------------------------------------------------------------------------
    def mul(self, a, b):
        return np.einsum('...i,...j,ijk->...k', a, b, self.T)

    def d(self, a, v):
        return a @ self.D[v].T

    def _series(self, a, coefs):
        """sum_k coefs[k] * (a - a[0])^k for a scalar jet a (1-D array)."""
        delta = a.copy(); delta[0] = 0.0
        out = self.const(coefs[0]); pw = self.const(1.0)
        for k in range(1, self.N + 1):
            pw = self.mul(pw, delta)
            out = out + coefs[k] * pw
        return out

    def exp(self, a):
        e0 = np.exp(a[0]); return self._series(a, [e0 / math.factorial(k) for k in range(self.N + 1)])

    def recip(self, a):
        return self._series(a, [(-1.0) ** k / a[0] ** (k + 1) for k in range(self.N + 1)])

    def sqrt(self, a):
        r = np.sqrt(a[0]); c = [r]; binom = 1.0
        for k in range(1, self.N + 1):
            binom *= (0.5 - (k - 1)) / k
            c.append(r * binom / a[0] ** k)
        return self._series(a, c)

    def div(self, a, b):
        return self.mul(a, self.recip(b))

    def compose(self, taylor, a, base):
        """f(a) for a univariate f given by Taylor coefficients `taylor` about the point `base`;
        a[0]-base may be a (tiny, complex-step) non-zero number, so ALL supplied powers are used."""
        delta = a.copy(); delta[0] = a[0] - base
        out = self.const(taylor[0]); pw = self.const(1.0)
        for k in range(1, len(taylor)):
            pw = self.mul(pw, delta)
            out = out + taylor[k] * pw
        return out

    def poly(self, a, pcoefs):
        """ordinary polynomial sum_k pcoefs[k] a^k (exact in jet arithmetic)."""
        out = self.const(pcoefs[0]); pw = self.const(1.0)
        for k in range(1, len(pcoefs)):
            pw = self.mul(pw, a); out = out + pcoefs[k] * pw
        return out

    # ---- jet matrices -------------------------------------------------------------------
    def matmul(self, A, B):
        return np.einsum('aci,cbj,ijk->abk', A, B, self.T)

    def matinv(self, G):
        Dn = G.shape[0]
        X = np.zeros_like(G); X[:, :, 0] = np.linalg.inv(G[:, :, 0])
        two = np.zeros_like(G); two[:, :, 0] = 2 * np.eye(Dn)
        for _ in range(3):                       # exact to jet degree 2^3-1 >= order
            X = self.matmul(X, two - self.matmul(G, X))
        return X


# ---------------------------------------------------------------------------------------------
# registered model
# ---------------------------------------------------------------------------------------------
W_POLY = [1.0, -1.0, 0.0, 1.0 / 3.0]                     # W = 1 - p + p^3/3


def W(p):  return 1 - p + p ** 3 / 3
def W1(p): return p * p - 1
def W2(p): return 2 * p
def W3(p): return 2.0
def U(p):  return 0.5 * W1(p) ** 2 - (2 / 3) * W(p) ** 2
def U1(p): return W1(p) * W2(p) - (4 / 3) * W(p) * W1(p)
def U2(p): return W2(p) ** 2 + W1(p) * W3(p) - (4 / 3) * (W1(p) ** 2 + W(p) * W2(p))
def U3(p): return 3 * W2(p) * W3(p) - (4 / 3) * (3 * W1(p) * W2(p) + W(p) * W3(p))


def jet_W(JS, p):  return JS.poly(p, W_POLY)
def jet_W1(JS, p): return JS.poly(p, [-1.0, 0.0, 1.0])
def jet_W2(JS, p): return JS.poly(p, [0.0, 2.0])
def jet_U(JS, p):
    w, w1 = jet_W(JS, p), jet_W1(JS, p)
    return 0.5 * JS.mul(w1, w1) - (2 / 3) * JS.mul(w, w)
def jet_U1(JS, p):
    w, w1, w2 = jet_W(JS, p), jet_W1(JS, p), jet_W2(JS, p)
    return JS.mul(w1, w2) - (4 / 3) * JS.mul(w, w1)


def background_taylor(rho, phi, dphi, nder=3, sign=+1):
    """Taylor data of an exact background through a point, from the ODEs
         rho'^2 = 1 + rho^2 (phi'^2/12 - U/6),  rho'' = -rho (phi'^2/4 + U/6),  phi'' = U_phi - 4 (rho'/rho) phi'.
    Returns derivative lists [rho, rho', rho'', rho'''], [phi, phi', phi'', phi''']."""
    rp = sign * math.sqrt(1 + rho * rho * (dphi ** 2 / 12 - U(phi) / 6))
    rpp = -rho * (dphi ** 2 / 4 + U(phi) / 6)
    ppp = U1(phi) - 4 * (rp / rho) * dphi
    rppp = -rp * (dphi ** 2 / 4 + U(phi) / 6) - rho * (dphi * ppp / 2 + U1(phi) * dphi / 6)
    App = rpp / rho - (rp / rho) ** 2
    pppp = U2(phi) * dphi - 4 * (App * dphi + (rp / rho) * ppp)
    return [rho, rp, rpp, rppp][:nder + 1], [phi, dphi, ppp, pppp][:nder + 1]


def taylor_from_derivs(d):
    return [d[k] / math.factorial(k) for k in range(len(d))]


# ---------------------------------------------------------------------------------------------
# brute-force geometry at the base point
# ---------------------------------------------------------------------------------------------
def geometry(JS, g, varmap):
    """g: (D,D,n) jet matrix of a metric; varmap[c] = jet variable index of coordinate c or None.
    Returns dict with ginv (jets), Gamma (jets, valid to degree order-1), Gamma0, Ricci (values)."""
    Dn = g.shape[0]
    ginv = JS.matinv(g)
    dg = np.zeros((Dn, Dn, Dn, JS.n), dtype=complex)               # dg[a,b,c] = d_c g_ab
    for c in range(Dn):
        if varmap[c] is not None:
            dg[:, :, c, :] = JS.d(g, varmap[c])
    Gl = 0.5 * (dg.transpose(0, 2, 1, 3) + dg - dg.transpose(2, 0, 1, 3))   # Gamma_{a,bc}
    Gam = np.einsum('adi,dbcj,ijk->abck', ginv, Gl, JS.T)
    dGam = np.zeros((Dn, Dn, Dn, Dn), dtype=complex)                # dGam[a,b,c,e] = d_e Gamma^a_bc
    for e in range(Dn):
        if varmap[e] is not None:
            dGam[:, :, :, e] = JS.d(Gam, varmap[e])[..., 0]
    G0 = Gam[..., 0]
    Ric = (np.einsum('adba->bd', dGam) - np.einsum('aabd->bd', dGam)
           + np.einsum('aae,edb->bd', G0, G0) - np.einsum('ade,eab->bd', G0, G0))
    return dict(ginv=ginv, Gam=Gam, Gam0=G0, Ric=Ric, g0=g[..., 0], ginv0=ginv[..., 0])


def einstein_tensor(JS, g, varmap):
    geo = geometry(JS, g, varmap)
    R = np.einsum('ab,ab->', geo['ginv0'], geo['Ric'])
    return geo['Ric'] - 0.5 * R * geo['g0'], geo


def field_residuals(JS, g, phi, varmap):
    """E_AB = R_AB - d_A phi d_B phi - (2/3) U g_AB   and   E_phi = Box phi - U_phi   at the base point."""
    Dn = g.shape[0]
    geo = geometry(JS, g, varmap)
    dphi = np.zeros(Dn, dtype=complex); ddphi = np.zeros((Dn, Dn), dtype=complex)
    for a in range(Dn):
        if varmap[a] is None: continue
        da = JS.d(phi, varmap[a]); dphi[a] = da[0]
        for b in range(Dn):
            if varmap[b] is None: continue
            ddphi[a, b] = JS.d(da, varmap[b])[0]
    p0 = phi[0]
    E = geo['Ric'] - np.outer(dphi, dphi) - (2 / 3) * U(p0) * geo['g0']
    box = np.einsum('ab,ab->', geo['ginv0'], ddphi - np.einsum('cab,c->ab', geo['Gam0'], dphi))
    return E, box - U1(p0), geo


def junction_residuals(JS, g, phi, varmap, sigma, dsigma, ynorm=4):
    """Z2 shell on the coordinate surface x^ynorm = const, bulk kept on the side of DEcreasing x^ynorm,
    unit normal n_A = delta_A^y / sqrt(g^yy) pointing towards increasing x^ynorm (out of the kept bulk).
        J_mn  := K_mn - sigma(phi)/6 * g_mn   (tangential m,n),   K_mn = -Gamma^y_mn / sqrt(g^yy)
        J_phi := n^A d_A phi + sigma'(phi)/2
    Both vanish on the registered background.  sigma, dsigma: callables on complex numbers."""
    Dn = g.shape[0]
    geo = geometry(JS, g, varmap)
    gyy = geo['ginv0'][ynorm, ynorm]
    tang = [a for a in range(Dn) if a != ynorm]
    K = -geo['Gam0'][ynorm][np.ix_(tang, tang)] / np.sqrt(gyy)
    h = geo['g0'][np.ix_(tang, tang)]
    dphi = np.zeros(Dn, dtype=complex)
    for a in range(Dn):
        if varmap[a] is not None: dphi[a] = JS.d(phi, varmap[a])[0]
    nphi = np.einsum('a,a->', geo['ginv0'][ynorm], dphi) / np.sqrt(gyy)
    p0 = phi[0]
    return K - sigma(p0) / 6 * h, nphi + dsigma(p0) / 2, K, h


if __name__ == "__main__":
    # self-test 1: unit round S^3 in hyperspherical coordinates has Ric = 2 g
    JS = JetSpace(2, 2)
    ch, th = JS.var(0, 0.7), JS.var(1, 1.1)
    def jsin(a):
        s0, c0 = np.sin(a[0]), np.cos(a[0]); return JS._series(a, [s0, c0, -s0 / 2])
    g = np.zeros((3, 3, JS.n), dtype=complex)
    g[0, 0] = JS.const(1.0); s1 = jsin(ch); g[1, 1] = JS.mul(s1, s1); s2 = jsin(th)
    g[2, 2] = JS.mul(g[1, 1], JS.mul(s2, s2))
    geo = geometry(JS, g, [0, 1, None])
    print("S^3 test  max|Ric - 2 g| =", np.max(np.abs(geo['Ric'] - 2 * geo['g0'])))
    # self-test 2: unit dS_4 in flat slicing has Ric = 3 gamma
    JS = JetSpace(1, 2); tau = JS.var(0, 0.3); e2 = JS.exp(2 * tau)
    g = np.zeros((4, 4, JS.n), dtype=complex); g[0, 0] = JS.const(-1.0)
    for i in (1, 2, 3): g[i, i] = e2
    geo = geometry(JS, g, [0, None, None, None])
    print("dS_4 test max|Ric - 3 g| =", np.max(np.abs(geo['Ric'] - 3 * geo['g0'])))
