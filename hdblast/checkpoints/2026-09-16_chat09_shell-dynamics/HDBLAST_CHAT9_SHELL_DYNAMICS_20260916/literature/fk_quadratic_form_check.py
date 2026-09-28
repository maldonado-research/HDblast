#!/usr/bin/env python3
"""Frolov-Kofman (hep-th/0309002) self-adjoint form (17a) with the SOFT boundary condition (14), for the registered shell.
FK (17a):  -(g Y_w)_w + f Y = lam g Y,  Y = a^2 Phi,  g = 3/(2 a phi_w^2),  f = 1/a,  lam = m^2 + 4H^2  (here H = 1, a = rho, dw = dy/rho).
FK (14) in the registered chart (bulk on y < y_b):  Y_y B = lam Y / rho^2  at y_b,  B = phi''/phi' + sigma_t''/2   (eigenvalue-dependent b.c.).
Multiplying (17a) by Y and integrating from the regular cone to the shell gives the identity
     Q[Y] := int 3 Y_y^2/(2 rho^2 phi'^2) dy + int Y^2/rho^2 dy  =  lam * N[Y],
     N[Y] := int 3 Y^2/(2 rho^4 phi'^2) dy  +  3 Y_b^2/(2 rho_b^4 phi'_b^2 B)        (boundary weight has the sign of B).
Q > 0 always, so lam < 0 (m^2 < -4H^2, impossible with FK's rigid b.c. (16)) requires N < 0, which requires B < 0.
This script evaluates Q, N at the numerically found eigenvalue and checks Q = lam N.  Floating point diagnostic only.
run: python3 fk_quadratic_form_check.py"""
import math, os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prior_art_crosscheck as pc
bg = pc.bg
I = bg.I_plus(); c = 2/I - 4/3; t = 1e-3; dv = 2e-4; v0 = 1e-2
sol = bg.solve_shell(t, c, guess=(math.log10(8.7855e-7), 8.6), dv=dv)
mu2 = pc.eigen(sol['phi_h'], sol['v_b'], t, c, -9.0, -6.0, dv)
lam = mu2 + 4
rho, phi, s = [float(x) for x in bg.series(sol['phi_h'], v0)]
p = (5 + math.sqrt(9 - 4*mu2))/2
Y = 1.0; Z = 1.5*p*Y/v0/(rho*rho*s*s)
n = max(1, int(math.ceil((sol['v_b'] - v0)/dv))); h = (sol['v_b'] - v0)/n
def f(rho, phi, s, Y, Z):
    rp = math.sqrt(1 + rho*rho*(s*s/12 - pc.U(phi)/6))
    return (rp, s, pc.U1(phi) - 4*(rp/rho)*s, (2/3)*rho*rho*s*s*Z, Y/(rho*rho) - 1.5*lam*Y/(rho**4*s*s))
def dens(rho, phi, s, Y, Z):
    Yy = (2/3)*rho*rho*s*s*Z
    return (1.5*Yy*Yy/(rho*rho*s*s), Y*Y/(rho*rho), 1.5*Y*Y/(rho**4*s*s))
acc = [0.0, 0.0, 0.0]; d0 = dens(rho, phi, s, Y, Z)
# analytic cone pieces (Y ~ y^p, phi' ~ k y, rho ~ y): tiny, added for completeness
k = s/v0
acc[0] += 1.5*p*p/(k*k)/(2*p-5)/v0**5
acc[1] += v0**(2*p-1)/(2*p-1)/v0**(2*p)
acc[2] += 1.5/(k*k)/(2*p-5)/v0**5
for _ in range(n):
    y0 = (rho, phi, s, Y, Z)
    k1 = f(*y0); k2 = f(*[a + h/2*b for a, b in zip(y0, k1)]); k3 = f(*[a + h/2*b for a, b in zip(y0, k2)]); k4 = f(*[a + h*b for a, b in zip(y0, k3)])
    rho, phi, s, Y, Z = [a + h/6*(b1 + 2*b2 + 2*b3 + b4) for a, b1, b2, b3, b4 in zip(y0, k1, k2, k3, k4)]
    d1 = dens(rho, phi, s, Y, Z)
    for i in range(3): acc[i] += 0.5*h*(d0[i] + d1[i])
    d0 = d1
rp = math.sqrt(1 + rho*rho*(s*s/12 - pc.U(phi)/6))
B = pc.U1(phi)/s - 4*rp/rho + 0.5*(2*pc.W2(phi))
Q = acc[0] + acc[1]; Nbulk = acc[2]; Nbdy = 1.5*Y*Y/(rho**4*s*s*B); N = Nbulk + Nbdy
nrm = Y*Y
out = dict(mu2=mu2, lam=lam, B=B, B_rho_b2=B*rho*rho, Q=Q/nrm, N_bulk=Nbulk/nrm, N_boundary=Nbdy/nrm, N=N/nrm, Q_over_lamN=Q/(lam*N))
for a, b in out.items(): print("%-12s = %.10g" % (a, b))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fk_quadratic_form_check_output.json"), "w"), indent=1)
