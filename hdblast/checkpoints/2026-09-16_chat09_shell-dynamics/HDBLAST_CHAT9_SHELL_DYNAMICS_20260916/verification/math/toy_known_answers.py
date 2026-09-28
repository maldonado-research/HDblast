#!/usr/bin/env python3
"""VERIFICATION (math agent) - known-answer tests.

A. Exact residual-gauge solution at mu2 = -4 fed through the ORCHESTRATOR's own full_rhs (master equation M_y).
   For harmonics with nabla_mu nabla_nu Y = -gamma Y (mu2 = -4) the gauge vector eps^y = -rho^2 e'(y) Y, eps_par = e(y) Y
   preserves the longitudinal form; demanding xi = -2 psi gives  rho^2 e'' + 4 rho rho' e' + 2 e = 0  and
        psi_g = -(rho'/rho) rho^2 e' - e,   chi_g = -phi' rho^2 e'.
   psi_g must solve (M_y) at mu2 = -4 and chi_g must equal the Codazzi expression.  (Any coefficient error in M_y that is
   not proportional to (mu2+4) would show up here.)  Background with phi' = O(1) (phi_h = -0.5) so every term matters.
B. Sturm-Liouville / Rayleigh identity for the registered shell.  With P = rho^2 psi, p = 1/(rho^2 phi'^2), w = 1/(rho^4 phi'^2),
   lam = mu2 + 4 the system is   -(p P')' + (2/(3 rho^2)) P = lam w P,   shell:  p P' = (lam/beta) P,
   beta = rho^4 phi' (phi'' + sigma'' phi'/2)  at the shell.  Hence
        int [p P'^2 + 2P^2/(3rho^2)] dy  =  lam * ( int w P^2 dy + P_b^2/beta )          (R)
   LHS > 0: a mode with mu2 < -4 exists only if beta < 0 and then has negative "w-norm" N, but lam*N > 0 for every mode.
C. Toy with exact answer: test scalar of mass M on dS-sliced AdS5 (rho = sinh y) with a Robin brane condition
   chi' + L chi = 0 at y_b (Langlois-Sasaki hep-th/0302069 type problem).  Exact regular solution
        chi = sinh^b(y) 2F1(a1, a2; b + 5/2; -sinh^2 y),  b = -3/2 + sqrt(9/4 - mu2),  a_{1,2} = (b + 2 +- sqrt(4 + M^2))/2
   (verified below against the ODE by finite differences); roots of the exact mismatch are compared with roots found by
   the orchestrator's shooting logic (leading power start at y0 = 1e-3, fixed-step RK4, normalised mismatch, bisection).
"""
import numpy as np, math, sys, os, json, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "background"))
import hdblast_background as bg
spec = importlib.util.spec_from_file_location("orch", os.path.join(ROOT, "stability", "orchestrator_derivation", "scalar_spectrum_orchestrator.py"))
orch = importlib.util.module_from_spec(spec); spec.loader.exec_module(orch)

def rk4(f, Y, a, b, n):
    h = (b - a)/n
    for _ in range(n):
        k1 = f(Y); k2 = f(Y + h/2*k1); k3 = f(Y + h/2*k2); k4 = f(Y + h*k3)
        Y = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
    return Y

print("=== A. exact mu2=-4 residual-gauge solution through the orchestrator's full_rhs")
phi_h = -0.5
b0 = bg.integrate(phi_h, 0.5, dv=1e-4)
rho, phi, s = b0
rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
e0, e1 = 0.7, -0.4
def gauge_rhs(Y):
    rho, phi, s, e, de = Y
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
    return np.array([rp, s, bg.U1(phi) - 4*Hh*s, de, -(4*rho*rp*de + 2*e)/rho**2])
psi_g = lambda rho, Hh, s, e, de: -Hh*rho*rho*de - e
dpsi_g = lambda rho, Hh, s, e, de: s*s*rho*rho*de/3 + 2*Hh*Hh*rho*rho*de + 2*Hh*e
Yg = rk4(gauge_rhs, np.array([rho, phi, s, e0, e1]), 0.5, 2.0, 6000)
Yo = rk4(lambda Y: orch.full_rhs(Y, -4.0), np.array([rho, phi, s, psi_g(rho, Hh, s, e0, e1), dpsi_g(rho, Hh, s, e0, e1)]), 0.5, 2.0, 6000)
r2, p2, s2, e2, de2 = Yg; rp2 = math.sqrt(1 + r2*r2*(s2*s2/12 - bg.U(p2)/6)); H2 = rp2/r2
print("  background at y=2: rho=%.5f phi=%.5f phi'=%.5f" % (r2, p2, s2))
print("  psi  at y=2: gauge construction %.12f   orchestrator master eq %.12f   diff %.2e" % (psi_g(r2, H2, s2, e2, de2), Yo[3], psi_g(r2, H2, s2, e2, de2) - Yo[3]))
print("  psi' at y=2: gauge construction %.12f   orchestrator master eq %.12f   diff %.2e" % (dpsi_g(r2, H2, s2, e2, de2), Yo[4], dpsi_g(r2, H2, s2, e2, de2) - Yo[4]))
chi_cod = -3*(Yo[4] + 2*H2*Yo[3])/s2
print("  chi  at y=2: gauge -phi' rho^2 e' = %.12f   Codazzi expression from orchestrator psi = %.12f" % (-s2*r2*r2*de2, chi_cod))
for mu2_wrong in (-3.9,):
    Yw = rk4(lambda Y: orch.full_rhs(Y, mu2_wrong), np.array([rho, phi, s, psi_g(rho, Hh, s, e0, e1), dpsi_g(rho, Hh, s, e0, e1)]), 0.5, 2.0, 6000)
    print("  sensitivity: same data with mu2 = %.1f in the master equation -> psi(2) differs by %.2e" % (mu2_wrong, Yw[3] - Yo[3]))

print("=== B. Rayleigh identity on the registered shell")
sol = json.load(open(os.path.join(ROOT, "background", "REGISTERED_SHELL_FLOAT_SOLUTION.json")))
t, c, ph, yb = sol["t"], sol["c"], sol["phi_h"], sol["v_b"]
mu0 = -7.717871623
n = 40000; y0 = 1e-3
alpha = 0.5 + math.sqrt(2.25 - mu0)
b = bg.series(ph, y0)
Y = np.array([b[0], b[1], b[2], 1.0, alpha/y0, 0.0, 0.0])      # + two quadrature accumulators
def rhsB(Y):
    d = orch.full_rhs(Y[:5], mu0)
    rho, phi, s, psi, dpsi = Y[:5]
    rp = d[0]; P = rho*rho*psi; dP = 2*rho*rp*psi + rho*rho*dpsi
    return np.concatenate([d, [dP*dP/(rho*rho*s*s) + 2*P*P/(3*rho*rho), P*P/(rho**4*s*s)]])
Y = rk4(rhsB, Y, y0, yb, n)
rho, phi, s, psi, dpsi, Q1, Q2 = Y
rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); spp = bg.U1(phi) - 4*rp/rho*s
beta = rho**4*s*(spp + 0.5*2*bg.W2(phi)*s)
P = rho*rho*psi
lam = mu0 + 4
print("  phi''+sigma''phi'/2 at shell = %.6e  (t = %.0e)   beta = %.6e  (negative => one mode with mu2 < -4 allowed)" % (spp + bg.W2(phi)*s, t, beta))
print("  LHS int[pP'^2+2P^2/(3rho^2)] = %.8e ;  lam*(int wP^2 + P_b^2/beta) = %.8e ; ratio %.8f" % (Q1, lam*(Q2 + P*P/beta), Q1/(lam*(Q2 + P*P/beta))))
print("  int w P^2 = %.6e,  P_b^2/beta = %.6e  => w-norm N = %.6e (negative), lam*N = %.6e (positive)" % (Q2, P*P/beta, Q2 + P*P/beta, lam*(Q2 + P*P/beta)))
# shell BC check in the pP' form
dP = 2*rho*rp*psi + rho*rho*dpsi
print("  shell condition p P' - lam P/beta = %.3e (relative %.1e)" % (dP/(rho*rho*s*s) - lam*P/beta, (dP/(rho*rho*s*s) - lam*P/beta)/abs(lam*P/beta)))

print("=== C. exact-solution toy (test scalar on dS-sliced AdS5 with Robin brane)")
def hyp2f1(a, b_, c_, z, tol=1e-17):
    term = 1.0; tot = 1.0; k = 0
    while abs(term) > tol*abs(tot) or k < 5:
        term *= (a + k)*(b_ + k)/((c_ + k)*(k + 1))*z; tot += term; k += 1
        if k > 200000: raise RuntimeError
    return tot
def chi_exact(y, mu2, M2):
    be = -1.5 + math.sqrt(2.25 - mu2); a1 = (be + 2 + math.sqrt(4 + M2))/2; a2 = (be + 2 - math.sqrt(4 + M2))/2; cc = be + 2.5
    # Pfaff: F(a1,a2;c;-sh^2) = ch^{-2 a1} F(a1, c-a2; c; th^2)
    return math.sinh(y)**be*math.cosh(y)**(-2*a1)*hyp2f1(a1, cc - a2, cc, math.tanh(y)**2)
def d_exact(y, mu2, M2, h=2e-3):
    f = lambda x: chi_exact(x, mu2, M2)
    D1 = (f(y + h) - f(y - h))/(2*h); D2 = (f(y + h/2) - f(y - h/2))/h
    return (4*D2 - D1)/3
M2 = 0.5; yb2 = 1.5
# verify ODE
mu_t = -1.3; yy = 0.9; h = 1e-3
f = lambda x: chi_exact(x, mu_t, M2)
d2 = (-f(yy + 2*h) + 16*f(yy + h) - 30*f(yy) + 16*f(yy - h) - f(yy - 2*h))/(12*h*h)
d1 = (-f(yy + 2*h) + 8*f(yy + h) - 8*f(yy - h) + f(yy - 2*h))/(12*h)
print("  ODE residual of the closed form at y=0.9, mu2=-1.3: %.2e (terms ~ %.2e)" % (d2 + 4/math.tanh(yy)*d1 + (mu_t/math.sinh(yy)**2 - M2)*f(yy), abs(d2)))
def shoot_toy(mu2, L, n=8000, y0=1e-3):
    be = -1.5 + math.sqrt(2.25 - mu2)
    def rhs(Y):
        return np.array([1.0, Y[2], -4/math.tanh(Y[0])*Y[2] - (mu2/math.sinh(Y[0])**2 - M2)*Y[1]])
    Y = rk4(rhs, np.array([y0, 1.0, be/y0]), y0, yb2, n)
    Bc = Y[2] + L*Y[1]
    return Bc/(abs(Y[2]) + abs(L*Y[1]) + 1e-300)
def exact_mis(mu2, L):
    return d_exact(yb2, mu2, M2) + L*chi_exact(yb2, mu2, M2)
def bisect(fn, a, b_, it=50):
    fa = fn(a)
    for _ in range(it):
        m = 0.5*(a + b_); fm = fn(m)
        if fa*fm <= 0: b_ = m
        else: a, fa = m, fm
    return 0.5*(a + b_)
for L in (0.0, -0.4, 0.8, -1.5):
    grid = np.linspace(-12, 2.2, 72)
    vs = [shoot_toy(m, L, n=2000) for m in grid]; ve = [exact_mis(m, L) for m in grid]
    rs = [bisect(lambda m: shoot_toy(m, L), grid[i], grid[i+1], 40) for i in range(71) if vs[i]*vs[i+1] < 0]
    re = [bisect(lambda m: exact_mis(m, L), grid[i], grid[i+1], 40) for i in range(71) if ve[i]*ve[i+1] < 0]
    print("  Robin L=%+.1f, M^2=%.1f, y_b=%.1f: shooting roots %s | exact (2F1) roots %s" % (L, M2, yb2, ["%.8f" % r for r in rs], ["%.8f" % r for r in re]))
# massless Neumann: exact zero mode mu2 = 0
M2 = 0.0
r0 = bisect(lambda m: shoot_toy(m, 0.0), -0.5, 0.5, 45)
print("  massless scalar, Neumann brane: shooting root %.3e (exact 0)" % r0)
print("=== pure AdS5 + dS brane, no bulk scalar: (C) => psi = C/rho^2, (H) => (mu2+4) psi = 0: no scalar mode except the mu2=-4 cone-singular gauge/radion mode (analytic).")
