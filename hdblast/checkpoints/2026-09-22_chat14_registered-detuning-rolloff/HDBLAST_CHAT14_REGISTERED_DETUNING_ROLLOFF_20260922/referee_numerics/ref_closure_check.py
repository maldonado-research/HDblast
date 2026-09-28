#!/usr/bin/env python3
"""Referee check of the v2 boundary closure, coordinate map and dissipation scaling (numpy only).
A. algebra of the 4th-order closure: exactness on polynomials, truncation order of D1/D2 at the last nodes for both closures
B. stability: eigenvalues of the semi-discrete operator for u_tt = u_zz (Neumann right / Dirichlet left) with each closure
C. chain-rule check on the actual stretched grid; grid-stretch ratio per cell; effective Kreiss-Oliger damping rate in v2 vs v1."""
import math, numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rolloff5d_v2 as V2

h = 1.0
def ghosts(F, gx, closure):
    if closure == 2: return F[-2] + 2*h*gx, F[-3] + 4*h*gx
    e1 = 4*h*gx - 2*F[-3] + 6*F[-2] - (10/3)*F[-1] + F[-4]/3
    e2 = 5*e1 - 10*F[-1] + 10*F[-2] - 5*F[-3] + F[-4]
    return e1, e2
def D1(E): return (E[:-4] - 8*E[1:-3] + 8*E[3:-1] - E[4:])/(12*h)
def D2(E): return (-E[:-4] + 16*E[1:-3] - 30*E[2:-2] + 16*E[3:-1] - E[4:])/(12*h*h)

print("=== A. closure algebra")
# A1: derive e1 symbolically-by-hand check: substitute quartic-extrapolated e2 into the centred D1 at node N and compare with the
#     standard 5-point one-sided formula f'(0) = (-F[-3] + 6F[-2] - 18F[-1] + 10F[0] + 3F[1])/(12h)
rng = np.random.default_rng(1)
for trial in range(3):
    F = rng.normal(size=8); gx = rng.normal()
    e1, e2 = ghosts(F, gx, 4)
    E = np.concatenate((F, [e1, e2]))
    d1N = D1(E)[-1]
    onesided = (-F[-4] + 6*F[-3] - 18*F[-2] + 10*F[-1] + 3*e1)/12
    print("trial %d: centred D1 at shell with ghosts = %.15f   target g = %.15f   one-sided 5-pt formula = %.15f" % (trial, d1N, gx, onesided))
# A2: exactness on polynomials of degree <= 4 (ghost values should equal the true polynomial values)
for deg in range(6):
    x = np.arange(-7, 1, 1.0)                       # nodes ..., -1, 0 (shell at 0)
    c = rng.normal(size=deg + 1); p = np.polynomial.Polynomial(c)
    F = p(x); gx = p.deriv()(0.0)
    for cl in (2, 4):
        e1, e2 = ghosts(F, gx, cl)
        print("deg %d closure %d: ghost errors e1-p(1)=%.2e e2-p(2)=%.2e" % (deg, cl, e1 - p(1.0), e2 - p(2.0)))
# A3: truncation order of D1 and D2 at nodes N and N-1 for a smooth function under grid refinement (in physical units)
print("--- truncation order (f = sin(2z+0.3) on z<=0, Neumann slope exact)")
for cl in (2, 4):
    prev = None
    for hh in (0.1, 0.05, 0.025, 0.0125):
        x = np.arange(-2.0, 1e-12, hh); f = np.sin(2*x + 0.3); g = 2*math.cos(0.3)
        e1, e2 = ghosts(f, g*hh, cl)          # ghosts in index space with slope g*hh
        E = np.concatenate((f, [e1, e2]))
        d1 = D1(E)/hh; d2 = D2(E)/hh**2
        errs = (abs(d1[-1] - g), abs(d1[-2] - 2*math.cos(2*x[-2] + 0.3)), abs(d2[-1] + 4*math.sin(0.3)), abs(d2[-2] + 4*math.sin(2*x[-2] + 0.3)))
        line = "closure %d h=%.4f  |D1 err| shell %.2e, N-1 %.2e | |D2 err| shell %.2e, N-1 %.2e" % (cl, hh, *errs)
        if prev is not None: line += "   orders: " + " ".join("%.2f" % (math.log(p_/e_)/math.log(2)) if e_ > 0 else "exact" for p_, e_ in zip(prev, errs))
        print(line); prev = errs

print("\n=== B. semi-discrete stability of u_tt = u_zz, Dirichlet left (two zero ghosts), Neumann right (u_z = 0) via each closure")
for cl in (2, 4):
    n = 60
    Mat = np.zeros((n, n))
    for j in range(n):
        F = np.zeros(n); F[j] = 1.0
        e1, e2 = ghosts(F, 0.0, cl)
        E = np.concatenate(([0.0, 0.0], F, [e1, e2]))
        Mat[:, j] = D2(E)
    ev = np.linalg.eigvals(Mat)
    # for u_tt = M u the modes are exp(±sqrt(lambda) t): need lambda real and <= 0
    print("closure %d: max Re(lambda) = %+.3e, max |Im(lambda)| = %.3e, min Re = %.3f (continuum: -pi^2/4/h^2 ... )" % (cl, ev.real.max(), np.abs(ev.imag).max(), ev.real.min()))
    grow = np.sqrt(ev.astype(complex)); print("   largest growth rate Re sqrt(lambda) (per unit h^-1 time) = %.3e" % grow.real.max())
# same for a first-order-in-time system with the actual RK4 + CFL 0.5 (amplification per step) including KO as coded in v2
print("--- fully discrete: RK4, dt = 0.5 h, system (u, v): u_t = v + KO(u), v_t = D2 u + KO(v); spectral radius of the step map")
def ko_v2(F, gx, cl, eps):
    e1, e2 = ghosts(F, gx, cl); E = np.concatenate(([0.0, 0.0], F, [e1, e2])); E = np.concatenate(([0.0], E, [2*E[-1] - E[-2]]))
    return eps*(E[:-6] - 6*E[1:-5] + 15*E[2:-4] - 20*E[3:-3] + 15*E[4:-2] - 6*E[5:-1] + E[6:])/(64*h)
for cl in (2, 4):
    for eps in (0.0, 0.05):
        n = 60; N2 = 2*n
        def rhs(y):
            u, v = y[:n], y[n:]
            e1, e2 = ghosts(u, 0.0, cl); E = np.concatenate(([0.0, 0.0], u, [e1, e2]))
            return np.concatenate((v + ko_v2(u, 0.0, cl, eps), D2(E) + ko_v2(v, 0.0, cl, eps)))
        dt = 0.5*h
        def step(y):
            k1 = rhs(y); k2 = rhs(y + dt/2*k1); k3 = rhs(y + dt/2*k2); k4 = rhs(y + dt*k3); return y + dt/6*(k1 + 2*k2 + 2*k3 + k4)
        S = np.zeros((N2, N2))
        for j in range(N2):
            y = np.zeros(N2); y[j] = 1.0; S[:, j] = step(y)
        ev = np.linalg.eigvals(S); rad = np.abs(ev).max()
        print("closure %d eps=%.2f: spectral radius of the RK4 step map = %.10f  (growth per unit time in h^-1 units: %+.3e)" % (cl, eps, rad, math.log(rad)/dt))

print("\n=== C. the stretched grid actually used (dzf = 1.5e-4, dzc = 8e-3, L = 10)")
z, zp, zpp = V2.build_grid(10.0, 1.5e-4, 8e-3, 0.06, 0.02)
print("nodes %d, dz at shell %.4e, dz far %.4e, z[0] = %.4f" % (len(z), zp[-1], zp[0], z[0]))
ratio = zp[1:]/zp[:-1]
i = np.argmax(np.abs(np.log(ratio)))
print("max per-cell stretch ratio dz_{i+1}/dz_i = %.4f at z = %.4f (i.e. %.1f%% per cell); ratio > 1.10 on %d cells, > 1.05 on %d cells" % (
    ratio[i], z[i], 100*(ratio[i] - 1), int((ratio > 1.10).sum()), int((ratio > 1.05).sum())))
# chain rule check: derivative of a smooth function through dz1/dz2 as coded (reimplemented) vs analytic
def dz1(F, zp): E = np.concatenate(([F[0], F[0]], F, [F[-1], F[-1]])); return D1(E)/zp        # interior only matters here
def dz2(F, zp, zpp): E = np.concatenate(([F[0], F[0]], F, [F[-1], F[-1]])); return D2(E)/zp**2 - D1(E)*zpp/zp**3
f = np.sin(3*z); f1 = 3*np.cos(3*z); f2 = -9*np.sin(3*z)
sl = slice(5, -5)
print("chain rule: max |f_z err| = %.2e, max |f_zz err| = %.2e (interior, f = sin 3z); at the transition region |z+0.06|<0.05: %.2e / %.2e" % (
    np.abs(dz1(f, zp) - f1)[sl].max(), np.abs(dz2(f, zp, zpp) - f2)[sl].max(),
    np.abs(dz1(f, zp) - f1)[np.abs(z + 0.06) < 0.05].max(), np.abs(dz2(f, zp, zpp) - f2)[np.abs(z + 0.06) < 0.05].max()))
# a wall-width feature (width 0.0127) crossing the coarse region: points per width
print("points across a 0.0127-wide feature: %.1f at the shell, %.1f at z = -0.1, %.1f at z = -1" % (0.0127/zp[-1], 0.0127/np.interp(-0.1, z, zp), 0.0127/np.interp(-1.0, z, zp)))
print("Kreiss-Oliger: v1 damping rate of the grid mode = eps/dz = %.1f (dz=4e-4) ; v2 as coded (h = 1 in xi) = eps/h = %.3f per unit t at every node;"
      " per time step dt = %.2e that is %.1e (v1 at CFL 0.5: 0.025 per step)" % (0.05/4e-4, 0.05, 0.5*zp.min(), 0.05*0.5*zp.min()))
z4, zp4, _ = V2.build_grid(10.0, 1.5e-4, 2e-3, 0.06, 0.02); r4 = zp4[1:]/zp4[:-1]
print("finecoarse grid (dzc = 2e-3): %d nodes, max stretch %.1f%% per cell" % (len(z4), 100*(r4.max() - 1)))
