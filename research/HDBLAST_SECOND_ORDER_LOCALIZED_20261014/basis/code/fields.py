"""Constraint maps and the inherited first-order reconstruction on T=0."""
import sympy as sp
from engine import (X, O, F, U, J3, register, d, inv_metric_jet, christoffel,
                    ricci, reduce_bg, S2)

O0 = O[0]

def gamma0():
    return sp.diag(*([O0**2] * 4))

def euler(expr):
    """D = x.grad (the slice restriction of E = X^A d_A at T=0)."""
    return sum(X[i] * d(expr, i) for i in range(4))

def harmonic_quadratic(A):
    """x^T A x for a symmetric traceless 4x4 sympy matrix A."""
    xv = sp.Matrix(X)
    return sp.expand((xv.T * A * xv)[0, 0])

def reconstruction(q, p):
    """Inherited ADM reconstruction (eq. (22) of the shell-Noether note).

    Returns dict with h (4x4 Matrix), dphi, dpi (4x4 Matrix, densitized,
    upper indices), dvarpi (densitized).
    """
    xv = sp.Matrix(X)
    h = 2 * q * (S2 * sp.eye(4) - 3 * xv * xv.T)
    dphi = -sp.Rational(3, 2) * (euler(q) + 4 * q) / (O0**2 * F[1])
    dpi = -sp.Rational(3, 2) * p / O0 * (xv * xv.T)
    dvarpi = -sp.Rational(3, 2) * O0 * (euler(p) + 5 * p) / F[1]
    return dict(h=sp.expand(h), dphi=sp.expand(dphi), dpi=sp.expand(dpi),
                dvarpi=sp.expand(dvarpi))

def phi0_grad():
    return [2 * X[i] * F[1] for i in range(4)]

def constraints(h, dphi, dpi, dvarpi, hh=None, psi2=None, pi2=None, varpi2=None,
                expand=True):
    """eps-jets of the Hamiltonian (non-densitized) and momentum (densitized)
    constraints along x0 + eps v + eps^2/2 w.

    Hamiltonian:
      H = R/2 + (2/det g)(pi^2/3 - |pi|^2) - varpi^2/(2 det g)
          - |D phi|^2/2 - U(phi)
    Momentum (densitized pi, varpi):
      M_i = -2 (d_j pi^j_i - Gamma^k_ji pi^j_k) + varpi d_i phi
    """
    g0 = gamma0()
    hh = hh if hh is not None else sp.zeros(4)
    g = [[J3(g0[i, j], h[i, j], hh[i, j]) for j in range(4)] for i in range(4)]
    # inverse metric jet with second-order part including hh
    gi0 = g0.inv()
    first = -gi0 * h * gi0
    from engine import J3 as _J3
    second = (sp.zeros(4) if _J3.FIRST_ONLY else
              2 * gi0 * h * gi0 * h * gi0 - gi0 * hh * gi0)
    gi = [[J3(gi0[i, j], first[i, j], second[i, j]) for j in range(4)]
          for i in range(4)]
    if expand:
        ex = lambda e: sp.expand(e)
        g = [[c.map(ex) for c in row] for row in g]
        gi = [[c.map(ex) for c in row] for row in gi]
    Gam = christoffel(g, gi)
    if expand:
        Gam = [[[c.map(sp.expand) for c in row] for row in mat] for mat in Gam]
    Ric = ricci(Gam)
    R = J3(0)
    for i in range(4):
        for j in range(4):
            R = R + gi[i][j] * Ric[i][j]
    phi = [J3(0, d(dphi, i), d(psi2, i) if psi2 is not None else 0) for i in range(4)]
    gp = phi0_grad()
    phi = [J3(gp[i]) + phi[i] for i in range(4)]
    grad2 = J3(0)
    for i in range(4):
        for j in range(4):
            grad2 = grad2 + gi[i][j] * phi[i] * phi[j]
    psi2v = psi2 if psi2 is not None else 0
    Upot = J3(U[0], U[1] * dphi, U[1] * psi2v + U[2] * dphi**2)
    pi2v = pi2 if pi2 is not None else sp.zeros(4)
    varpi2v = varpi2 if varpi2 is not None else 0
    pi = [[J3(0, dpi[i, j], pi2v[i, j]) for j in range(4)] for i in range(4)]
    varpi = J3(0, dvarpi, varpi2v)
    # kinetic terms only need the background metric (pi is O(eps))
    detg0 = O0**8
    trpi = J3(0)
    sqpi = J3(0)
    for i in range(4):
        for j in range(4):
            trpi = trpi + g[i][j] * pi[i][j]
            for k in range(4):
                for l in range(4):
                    sqpi = sqpi + g[i][k] * g[j][l] * pi[i][j] * pi[k][l]
    Hm = R * sp.Rational(1, 2) + (trpi * trpi * sp.Rational(1, 3) - sqpi) * (2 / detg0) \
        - varpi * varpi * (1 / (2 * detg0)) - grad2 * sp.Rational(1, 2) - Upot
    # momentum constraint
    pimix = [[J3(0) for _ in range(4)] for _ in range(4)]  # pi^j_i
    for jj in range(4):
        for ii in range(4):
            acc = J3(0)
            for k in range(4):
                acc = acc + pi[jj][k] * g[k][ii]
            pimix[jj][ii] = acc
    M = []
    for i in range(4):
        acc = J3(0)
        for jj in range(4):
            acc = acc + pimix[jj][i].d(jj)
            for k in range(4):
                acc = acc - Gam[k][jj][i] * pimix[jj][k]
        M.append(acc * (-2) + varpi * phi[i])
    return Hm, M
