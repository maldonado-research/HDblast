#!/usr/bin/env python3
"""Independent re-derivation (audit) of the B1 shell-matter algebra, written from scratch. Output audit/AUDIT_ALGEBRA.json"""
import json, os, sympy as sp
t = sp.symbols('t'); G, ps, k, Gam, kap = sp.symbols('G phi_s k Gamma kappa', positive=True)
a = sp.Function('a')(t); phi = sp.Function('phi')(t); N = sp.Function('N')(t); R = sp.Function('R')(t)
H = a.diff(t)/a; v = phi.diff(t)
m = sp.sqrt(G**2*(phi-ps)**2); p = k/a; om = sp.sqrt(p**2 + G**2*(phi-ps)**2)
res = {}
# single node gas with decay
rho = N*om/a**3; P = N*p**2/(3*om)/a**3
j = N/a**3*sp.diff(om, phi)                      # n d omega/d phi
Ndot = -Gam*(m/om)*N                              # time-dilated decay
Q = Gam*m*N/a**3
lhs = sp.diff(rho, t).subs(sp.Derivative(N, t), Ndot) + 3*H*(rho + P)
res['ward_with_decay'] = sp.simplify(lhs - (j*v - Q))
lhs_nodil = sp.diff(rho, t).subs(sp.Derivative(N, t), -Gam*N) + 3*H*(rho + P)
res['control_no_time_dilation_nonzero'] = sp.simplify(lhs_nodil - (j*v - Q)) != 0
res['j_formula'] = sp.simplify(j - G**2*(phi-ps)*N/a**3/om)
# shell budget: lambda = sigma/kappa^2, n.phi = -(sigma' + kappa^2 J)/2, claim d(lam+rho_tot)/dt + 3H(rho+p) = -2 (n.phi) v/kappa^2
sig = sp.Function('sigma'); J = sp.symbols('J')
lam = sig(phi)/kap**2
rho_tot_dot_plus = J*v  # from Ward + radiation ledger with Q cancelling
nphi = -(sig(phi).diff(phi) + kap**2*J)/2
res['shell_budget'] = sp.simplify(sp.diff(lam, t) + rho_tot_dot_plus - (-2*nphi*v/kap**2))
nphi_wrong = +(sig(phi).diff(phi) + kap**2*J)/2
res['control_plus_sign_nonzero'] = sp.simplify(sp.diff(lam, t) + rho_tot_dot_plus - (-2*nphi_wrong*v/kap**2)) != 0
# production integrals (instant preheating, massless energy)
kk, q = sp.symbols('kk q', positive=True)
n_int = sp.integrate(4*sp.pi*kk**2*sp.exp(-sp.pi*kk**2/q), (kk, 0, sp.oo))/(2*sp.pi)**3
e_int = sp.integrate(4*sp.pi*kk**3*sp.exp(-sp.pi*kk**2/q), (kk, 0, sp.oo))/(2*sp.pi)**3
res['n_integral'] = sp.simplify(n_int - q**sp.Rational(3, 2)/(8*sp.pi**3))
res['e_integral'] = sp.simplify(e_int - q**2/(4*sp.pi**4))
# rad decomposition (s+R)^2/36 - s^2/36 = sR/18 + R^2/36
s_, R_ = sp.symbols('s R')
res['rad_decomposition'] = sp.expand((s_+R_)**2/36 - s_**2/36 - (s_*R_/18 + R_**2/36))
out = {k: str(v) for k, v in res.items()}
out['all_ok'] = all(out[k] in ('0', 'True') for k in out)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AUDIT_ALGEBRA.json'), 'w'), indent=1)
print(out)
