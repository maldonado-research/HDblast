#!/usr/bin/env python3
"""Shared SymPy helper: first-order perturbation of R_AB - d_A phi d_B phi - (2/3) U g_AB on the HDBLAST background
   ds^2 = dy^2 + rho(y)^2 gamma,  gamma = unit dS_4 in flat slicing  (-dtau^2 + e^{2 tau} dx^2),  phi = phi0(y).
General (non-diagonal) metric perturbation dg (5x5, first order).  Everything is truncated at first order in eps."""
import sympy as sp
y, tau, x1, x2, x3, eps, k = sp.symbols('y tau x1 x2 x3 epsilon k')
X = [y, tau, x1, x2, x3]
rho, ph0 = sp.Function('rho')(y), sp.Function('phi0')(y)
U0, U1, U2 = sp.symbols('U0 U1 U2')
rp, pp = sp.symbols('rho1 phi1')
gam = sp.diag(-1, sp.exp(2*tau), sp.exp(2*tau), sp.exp(2*tau))
def background_metric():
    g = sp.zeros(5, 5); g[0, 0] = 1
    for m in range(4): g[m+1, m+1] = rho**2*gam[m, m]
    return g
def lin(e):
    e = sp.expand(e); return sp.expand(e.coeff(eps, 0) + eps*e.coeff(eps, 1))
def first_order_equations(dg, dphi=0):
    """returns dict {(A,B): first-order coefficient of R_AB - dphi dphi - (2/3) U g_AB}  (dg, dphi are the O(eps) parts)"""
    g0 = background_metric(); g0i = g0.inv()
    g = g0 + eps*dg
    gi = (g0i - eps*g0i*dg*g0i).applyfunc(sp.expand)
    n = 5
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = sum(gi[a, d]*(sp.diff(g[d, c], X[b]) + sp.diff(g[d, b], X[c]) - sp.diff(g[b, c], X[d])) for d in range(n))
                v = lin(sp.Rational(1, 2)*s); Gam[a][b][c] = Gam[a][c][b] = v
    out = {}
    phi = ph0 + eps*dphi; dph = [sp.diff(phi, v) for v in X]
    for b in range(n):
        for c in range(b, n):
            r = 0
            for a in range(n):
                r += sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
                for d in range(n):
                    r += Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a]
            e = lin(r - dph[b]*dph[c] - sp.Rational(2, 3)*(U0 + eps*U1*dphi)*g[b, c])
            out[(b, c)] = e.coeff(eps, 1)
    return out
def bg_reduce(e):
    r2 = (sp.diff(rho, y)**2/rho - 1/rho - rho*sp.diff(ph0, y)**2/3)
    p2 = (U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))
    for _ in range(2): e = e.subs({sp.diff(rho, y, 2): r2, sp.diff(ph0, y, 2): p2})
    e = e.subs({sp.diff(rho, y): rp, sp.diff(ph0, y): pp})
    return e.subs(U0, 6*(1 - rp**2)/rho**2 + pp**2/2)
def chr4(l, m, n_):
    xs = X[1:]; gi = gam.inv()
    return sum(sp.Rational(1, 2)*gi[l, s]*(sp.diff(gam[s, m], xs[n_]) + sp.diff(gam[s, n_], xs[m]) - sp.diff(gam[m, n_], xs[s])) for s in range(4))
