#!/usr/bin/env python3
"""HDBLAST Chat 9 - global registered background from the first cone (horizon) to the de Sitter shell.

Model (kappa_5 = 1):  W = 1 - phi + phi^3/3,  U = W_phi^2/2 - (2/3) W^2.
Chart: 5D metric  ds^2 = dvr^2 + rho(vr)^2 dS_4(unit Hubble),  vr = proper distance from the first cone.
    rho' = sqrt(1 + rho^2 (s^2/12 - U/6)),   phi' = s,   s' = U_phi - 4 (rho'/rho) s.
Regular cone data: rho ~ vr, phi = phi_h + U_phi(phi_h) vr^2/10 + O(vr^4).
Shell (Z2 brane) at vr_b with tension sigma_t(phi) = 2W + t (1 + c phi):
    rho'/rho = sigma_t/6,      s = -sigma_t'(phi)/2        (both evaluated at the shell)
and the brane Hubble^2 is h = 1/rho_b^2 (M462 normalisation q = sqrt(h) rho = 1 at the shell).
Only numpy + stdlib. Floating point diagnostic, not an interval certificate.
"""
import json, math, sys
import numpy as np

def W(p):   return 1 - p + p**3/3
def W1(p):  return p*p - 1
def W2(p):  return 2*p
def U(p):   return 0.5*W1(p)**2 - (2/3)*W(p)**2
def U1(p):  return W1(p)*W2(p) - (4/3)*W(p)*W1(p)
def U2(p):  return W2(p)**2 + W1(p)*2 - (4/3)*(W1(p)**2 + W(p)*W2(p))

def rhs(y):
    rho, phi, s = y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    return np.array([rp, s, U1(phi) - 4*(rp/rho)*s])

def series(phi_h, v0):
    u, u1, u2 = U(phi_h), U1(phi_h), U2(phi_h)
    # M465 jets: b = v + u v^3/36..., with rho_FRW analytic continuation: here spacelike side => signs flip: rho = v - u v^3/36
    rho = v0 - u*v0**3/36
    phi = phi_h + u1*v0**2/10 + u1*(u2/280 + u/630)*v0**4   # spacelike continuation of M465 series (tau^2 -> -v^2)
    s   = u1*v0/5 + 4*u1*(u2/280 + u/630)*v0**3
    return np.array([rho, phi, s])

def integrate(phi_h, v_end, dv=2e-4, v0=1e-2, record=False):
    y = series(phi_h, v0); v = v0
    n = max(1, int(math.ceil((v_end - v0)/dv))); dv = (v_end - v0)/n
    out = [(v, *y)] if record else None
    for _ in range(n):
        k1 = rhs(y); k2 = rhs(y + dv/2*k1); k3 = rhs(y + dv/2*k2); k4 = rhs(y + dv*k3)
        y = y + dv/6*(k1 + 2*k2 + 2*k3 + k4); v += dv
        if record: out.append((v, *y))
    return (y, np.array(out)) if record else y

def shell_residuals(phi_h, v_b, t, c, dv=2e-4):
    rho, phi, s = integrate(phi_h, v_b, dv)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    sig  = 2*W(phi) + t*(1 + c*phi)
    dsig = 2*W1(phi) + t*c
    return np.array([rp/rho - sig/6, s + dsig/2]), (rho, phi, s, rp)

def solve_shell(t, c, guess, dv=2e-4, tol=1e-13, itmax=40):
    """Newton in x = (log10(phi_h + 1), v_b)."""
    x = np.array(guess, float)
    for it in range(itmax):
        f = lambda xx: shell_residuals(-1 + 10**xx[0], xx[1], t, c, dv)[0]
        F = f(x); J = np.zeros((2, 2))
        for j, e in enumerate([1e-7, 1e-7]):
            xp = x.copy(); xp[j] += e; xm = x.copy(); xm[j] -= e
            J[:, j] = (f(xp) - f(xm))/(2*e)
        dx = np.linalg.solve(J, -F)
        # damp
        lam = 1.0
        while lam > 1e-4 and np.linalg.norm(f(x + lam*dx)) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < tol: break
    res, (rho, phi, s, rp) = shell_residuals(-1 + 10**x[0], x[1], t, c, dv)
    return dict(t=t, c=c, phi_h=-1 + 10**x[0], alpha=10**x[0], v_b=x[1], rho_b=rho, h=1/rho**2,
                phi_b=phi, s_b=s, rhop_b=rp, residual=[float(res[0]), float(res[1])], iters=it + 1)

def I_plus():
    A = lambda y: -(1/3)*(y + (2/3)*np.log(np.cosh(y)) - 1/(6*np.cosh(y)**2) + 1/6)
    x, w = np.polynomial.legendre.leggauss(100); tot = 0.0
    ed = np.linspace(0, 80, 801)
    for i in range(800):
        m = (ed[i] + ed[i+1])/2; hh = (ed[i+1] - ed[i])/2
        tot += hh*np.sum(w*np.exp(2*A(m + hh*x)))
    return tot

if __name__ == "__main__":
    Ip = I_plus(); c_star = 2/Ip - 4/3; t = 1e-3
    print("I_plus =", repr(Ip), " c_star =", repr(c_star), " h_lin =", 1/(6*Ip))
    sol = solve_shell(t, c_star, guess=(math.log10(8.7855e-7), 8.6))
    sol["delta_M462"] = sol["alpha"]/t**1.8
    sol["eta_M462"] = (sol["h"]/(t/(6*Ip)) - 1)/t
    sol["I_plus"] = Ip; sol["c_star"] = c_star
    print(json.dumps(sol, indent=1))
    print("registered certificate centre: eta = 0.12006349132995513, delta = 0.2206821504276063")
    # step-halving control
    sol2 = solve_shell(t, c_star, guess=(math.log10(sol['alpha']), sol['v_b']), dv=1e-4)
    print("step-halving: d(delta) = %.3e, d(eta) = %.3e" % (sol2['alpha']/t**1.8 - sol['delta_M462'],
          (sol2['h']/(t/(6*Ip)) - 1)/t - sol['eta_M462']))
    json.dump(sol, open("REGISTERED_SHELL_FLOAT_SOLUTION.json", "w"), indent=1)
