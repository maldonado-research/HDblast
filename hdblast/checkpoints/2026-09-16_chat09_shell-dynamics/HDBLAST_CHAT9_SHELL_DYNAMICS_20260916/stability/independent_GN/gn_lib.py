#!/usr/bin/env python3
"""Independent Gaussian-normal-gauge scalar perturbation solver for the HDBLAST registered shell.
numpy + stdlib only.  See GN_DERIVATION.md for the equations.

Perturbed metric: ds^2 = dy^2 + rho^2[(1+2 psi Y) gamma_mn + 2 E nabla_m nabla_n Y] , phi = phi0 + chi Y,
Box_gamma Y = mu2 Y.   Variables integrated: (chi, chi', Q=E', psi).
  psi'  = Q - phi' chi/3                                  (momentum constraint)
  Q'    = -4 H Q - 2 psi/rho^2                            (trace-free evolution eq.)
  chi'' = -4 H chi' - (4 psi' + mu2 Q) phi' - (mu2/rho^2 - U'') chi   (scalar eq.)
Hamiltonian constraint (monitored):  3H(4 psi' + mu2 Q) + 3(4+mu2) psi/rho^2 - phi' chi' + U' chi = 0
Trace evolution eq. (monitored):     psi'' + 8 H psi' + mu2 H Q + (6+mu2) psi/rho^2 + (2/3) U' chi = 0
"""
import math, sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "background"))
from hdblast_background import W, W1, W2, U, U1, U2, series, solve_shell, I_plus

def bg_rhs(y):
    rho, phi, s = y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    return np.array([rp, s, U1(phi) - 4*(rp/rho)*s])

def make_nodes(y0, y_b, n, y_sw=0.5):
    """Main nodes: geometric from y0 to y_sw (step/y <= ~ (y_b-y_sw)/(n*y_sw)), then uniform to y_b."""
    d = (y_b - y_sw)/n; r = d/y_sw
    ng = int(math.ceil(math.log(y_sw/y0)/math.log1p(r)))
    geo = y0*np.exp(np.linspace(0, math.log(y_sw/y0), ng+1))
    return np.concatenate([geo[:-1], np.linspace(y_sw, y_b, n+1)])

def background(phi_h, y_b, n, y0=2e-3):
    """RK4 background on main nodes + midpoints (2*len-1 nodes) from y0 to y_b. Returns rows (y, rho, phi, s)."""
    main = make_nodes(y0, y_b, n); ys = np.empty(2*len(main)-1); ys[0::2] = main; ys[1::2] = 0.5*(main[1:] + main[:-1])
    out = np.zeros((len(ys), 4)); st = series(phi_h, y0); out[0] = (y0, *st)
    for i in range(len(ys)-1):
        d = ys[i+1] - ys[i]
        k1 = bg_rhs(st); k2 = bg_rhs(st + d/2*k1); k3 = bg_rhs(st + d/2*k2); k4 = bg_rhs(st + d*k3)
        st = st + d/6*(k1 + 2*k2 + 2*k3 + k4); out[i+1] = (ys[i+1], *st)
    return out

def bg_coeffs(bg):
    y, rho, phi, s = bg.T
    H = np.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))/rho
    return dict(y=y, rho=rho, phi=phi, s=s, H=H, ir2=1/rho**2, U1=U1(phi), U2=U2(phi), spp=U1(phi) - 4*H*s)

def pert_rhs(i, X, mu2, C):
    chi, dchi, Q, psi = X
    H, ir2, s = C['H'][i], C['ir2'][i], C['s'][i]
    P = Q - s*chi/3
    return np.array([dchi, -4*H*dchi - (4*P + mu2*Q)*s - (mu2*ir2 - C['U2'][i])*chi, -4*H*Q - 2*psi*ir2, P])

def cone_start(mu2, C, phi_h):
    """Normalisable Frobenius branch chi ~ y^s+, s+ = -3/2 + sqrt(9/4 - mu2), with O(y^2) correction."""
    mu2 = np.asarray(mu2, float); y0 = C['y'][0]
    sp = -1.5 + np.sqrt(2.25 - mu2)
    k2 = -U(phi_h)/6; a = U1(phi_h)/5; u2 = U2(phi_h)
    c2 = (u2 + mu2*k2/3 - 4*sp*k2/3)/((sp+2)*(sp+5) + mu2)
    chi = y0**sp*(1 + c2*y0**2); dchi = y0**(sp-1)*(sp + (sp+2)*c2*y0**2)
    Q = 2*a/(3*(sp+3)*(sp+4))*y0**(sp+1)
    psi = -(a/3)*(sp+5)/((sp+3)*(sp+4))*y0**(sp+2)
    return np.array([chi, dchi, Q, psi]), sp

def integrate_pert(mu2, C, phi_h, record_every=0):
    X, sp = cone_start(mu2, C, phi_h)
    m = len(C['y']) - 1; n = m//2
    rec = []
    for j in range(n):
        i = 2*j; d = C['y'][i+2] - C['y'][i]
        k1 = pert_rhs(i, X, mu2, C); k2 = pert_rhs(i+1, X + d/2*k1, mu2, C)
        k3 = pert_rhs(i+1, X + d/2*k2, mu2, C); k4 = pert_rhs(i+2, X + d*k3, mu2, C)
        X = X + d/6*(k1 + 2*k2 + 2*k3 + k4)
        if record_every and (j+1) % record_every == 0: rec.append((i+2, X.copy()))
    return X, rec

def residuals(i, X, mu2, C):
    """Normalised residuals of Hamiltonian constraint, trace evolution eq. and yy equation at node i."""
    chi, dchi, Q, psi = X
    H, ir2, s, u1, spp = C['H'][i], C['ir2'][i], C['s'][i], C['U1'][i], C['spp'][i]
    dX = pert_rhs(i, X, mu2, C); ddchi, dQ, P = dX[1], dX[2], dX[3]
    dP = dQ - spp*chi/3 - s*dchi/3
    Hp = -ir2 - s*s/3
    hc_terms = [3*H*(4*P + mu2*Q), 3*(4+mu2)*psi*ir2, -s*dchi, u1*chi]
    ev_terms = [dP, 8*H*P, mu2*H*Q, (6+mu2)*psi*ir2, (2/3)*u1*chi]
    yy_terms = [-(4*dP + mu2*dQ), -2*H*(4*P + mu2*Q), -2*s*dchi, -(2/3)*u1*chi]
    f = lambda T: sum(T)/(sum(abs(t) for t in T) + 1e-300)
    return f(hc_terms), f(ev_terms), f(yy_terms)

def shell_data(sol):
    t, c, pb = sol['t'], sol['c'], sol['phi_b']
    return dict(sig=2*W(pb) + t*(1 + c*pb), dsig=2*W1(pb) + t*c, ddsig=2*W2(pb))

def mismatch(X, C, sol):
    """Gauge-invariant mismatch M and auxiliary shell-gauge quantities (eps fixed by J1: E_s'=0)."""
    chi, dchi, Q, psi = X
    rho, s, spp, H = C['rho'][-1], C['s'][-1], C['spp'][-1], C['H'][-1]
    sd = shell_data(sol)
    eps = rho*rho*Q
    chi_s = chi + s*eps
    M = dchi + spp*eps + 0.5*sd['ddsig']*chi_s
    return M, chi_s, eps

def metric_mismatch(X, C, sol):
    """Alternative: fix eps by the scalar junction (J2), return E_s', psi_s at the shell."""
    chi, dchi, Q, psi = X
    rho, s, spp, H = C['rho'][-1], C['s'][-1], C['spp'][-1], C['H'][-1]
    sd = shell_data(sol)
    eps = -(dchi + 0.5*sd['ddsig']*chi)/(spp + 0.5*sd['ddsig']*s)
    return Q - eps/rho**2, psi + H*eps, chi + s*eps, eps

def polish_shell(t, c, alpha, y_b, n=4000, y0=2e-3, itmax=12):
    """Newton polish of (ln alpha, y_b) with THIS file's integrator so that junctions hold to ~1e-13."""
    def F(x):
        bg = background(-1 + math.exp(x[0]), x[1], n, y0); y, rho, phi, s = bg[-1]
        H = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))/rho
        return np.array([H - (2*W(phi) + t*(1 + c*phi))/6, s + (2*W1(phi) + t*c)/2])
    x = np.array([math.log(alpha), y_b])
    for it in range(itmax):
        f = F(x); J = np.zeros((2, 2))
        for j in range(2):
            e = 1e-6; xp = x.copy(); xp[j] += e; xm = x.copy(); xm[j] -= e
            J[:, j] = (F(xp) - F(xm))/(2*e)
        dx = np.linalg.solve(J, -f); x = x + dx
        if np.linalg.norm(dx) < 1e-13: break
    bg = background(-1 + math.exp(x[0]), x[1], n, y0); y, rho, phi, s = bg[-1]
    sol = dict(t=t, c=c, alpha=math.exp(x[0]), phi_h=-1 + math.exp(x[0]), y_b=x[1], rho_b=rho, phi_b=phi, s_b=s,
               h=1/rho**2, residual=[float(v) for v in F(x)])
    return sol, bg
