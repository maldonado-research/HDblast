#!/usr/bin/env python3
"""Audit script 3: exact (sympy) checks of the formulas the preheating workstream relies on.
Writes verify/SYMPY_CHECKS.json.  Each entry records the simplified residual (0 means exact agreement)."""
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
out = {}
s, k, m0, q, G, v, h0 = sp.symbols('s kappa mu0 q G v h0', positive=True)

# 1. X = a^{3/2} chi removes the friction term: chi'' + 3h chi' + (k^2/a^2 + m^2) chi = 0
a = sp.Function('a')(s); chi = sp.Function('chi')(s); msq = sp.Function('m2')(s)
h = sp.diff(a, s)/a
eq_chi = sp.diff(chi, s, 2) + 3*h*sp.diff(chi, s) + (k**2/a**2 + msq)*chi
X = a**sp.Rational(3, 2)*chi
Om2 = k**2/a**2 + msq - sp.Rational(9, 4)*h**2 - sp.Rational(3, 2)*sp.diff(h, s)
res = sp.simplify(sp.expand(sp.diff(X, s, 2) + Om2*X - a**sp.Rational(3, 2)*eq_chi))
out['mode_equation_Omega2'] = str(res)

# wrong-formula control: dropping the (3/2)h' term must leave a nonzero residual
Om2_wrong = k**2/a**2 + msq - sp.Rational(9, 4)*h**2
out['mode_equation_wrong_control_nonzero'] = str(sp.simplify(sp.expand(sp.diff(X, s, 2) + Om2_wrong*X - a**sp.Rational(3, 2)*eq_chi)) != 0)

# 2. instant-preheating number: (1/2pi^2) int k^2 exp(-pi(k^2+mu0^2)/q) dk = q^{3/2}/(8 pi^3) exp(-pi mu0^2/q)
Nint = sp.integrate(k**2*sp.exp(-sp.pi*(k**2 + m0**2)/q), (k, 0, sp.oo))/(2*sp.pi**2)
out['instant_number'] = str(sp.simplify(Nint - q**sp.Rational(3, 2)/(8*sp.pi**3)*sp.exp(-sp.pi*m0**2/q)))
# mass-shift part of D: mu0^2 -> -(9/4)h^2 - (3/2)h'  =>  N ~ N0 (1 + pi(9/4 h^2 + 3/2 h')/q)
h1 = sp.symbols('h1')
lead = sp.series(sp.exp(-sp.pi*(-sp.Rational(9, 4)*h0**2 - sp.Rational(3, 2)*h1)/q), q, sp.oo, 2).removeO()
out['mass_shift_D_coefficients'] = str(sp.expand(sp.simplify((lead - 1)*q)))

# 3. residual-vacuum threshold: H^2 = ((sigma + x)/6)^2 - mu^2, H_vac^2 = (sigma/6)^2 - mu^2, x = kappa5^2 rho
sig, x, Hv = sp.symbols('sigma x H_vac', positive=True)
H2 = Hv**2 + ((sig + x)**2 - sig**2)/36
out['H2_expansion'] = str(sp.expand(H2 - Hv**2 - sig*x/18 - x**2/36))
xc = sp.solve(sp.Eq(H2 - Hv**2, Hv**2), x)
out['rho_crit_roots'] = [str(r) for r in xc]
out['rho_crit_matches'] = any(sp.simplify(r - (sp.sqrt(sig**2 + 36*Hv**2) - sig)) == 0 for r in xc)
# deceleration threshold: radiation (w=1/3) + vacuum; qdec: d/dt(aH) = 0  <=>  H' + H^2 = 0 with rho' = -4 H rho
# H^2 = Hv^2 + (2 sig x + x^2)/36 ; d(H^2)/dN = (2 sig + 2x)/36 * dx/dN, dx/dN = -4x
Ndot = -4*x
dH2dN = sp.diff(H2, x)*Ndot
acc = dH2dN/2 + H2          # a''/a = H' + H^2 = (1/2) dH^2/dN + H^2
xd = sp.solve(sp.Eq(acc, 0), x)
out['decel_roots'] = [str(r) for r in xd]
out['decel_matches'] = any(sp.simplify(r - (sp.sqrt(sig**2 + 108*Hv**2) - sig)/3) == 0 for r in xd)

# 4. Ward identity for the late-time particle description (per mode): rho_k = n_k w/a^3, p_k = n_k k^2/(3 a^2 w a^3)
#    with m^2 = G^2 (phi - phi*)^2 and <chi^2>_k = n_k/(w a^3):  d rho/ds + 3 h (rho + p) = j v
t = sp.symbols('t'); A = sp.Function('A')(t); ph = sp.Function('phi')(t); nk, phs = sp.symbols('n_k phi_s')
w = sp.sqrt(k**2/A**2 + G**2*(ph - phs)**2)
rho = nk*w/A**3; p = nk*k**2/(3*A**2*w*A**3); chi2 = nk/(w*A**3)
hh = sp.diff(A, t)/A; j = G**2*(ph - phs)*chi2
out['ward_identity_particle'] = str(sp.simplify(sp.diff(rho, t) + 3*hh*(rho + p) - j*sp.diff(ph, t)))

# 5. large-G law: R = b rho/crit with b = (G Dmax)^-3, rho = N m_e/a_e^3, N = (G v)^{3/2}/(8 pi^3), m_e = G De
De, Dm, ae, cr = sp.symbols('De Dmax a_e crit', positive=True)
Rcut = (G*Dm)**-3*(G*v)**sp.Rational(3, 2)/(8*sp.pi**3)*G*De/ae**3/cr
K = v**sp.Rational(3, 2)*De/(8*sp.pi**3*ae**3*Dm**3*cr)
out['K_over_sqrtG_law'] = str(sp.simplify(Rcut - K/sp.sqrt(G)))

# 6. Coleman-Weinberg source without the log: V = m^4/(64 pi^2) [ln - c]; dV/dphi (no log) = m^2 dm^2/dphi/(32 pi^2)
phi = sp.symbols('phi')
m2 = G**2*(phi - phs)**2
out['CW_j_nolog'] = str(sp.simplify(m2*sp.diff(m2, phi)/(32*sp.pi**2) - G**2*(phi - phs)*m2/(16*sp.pi**2)))

# 7. registered tension derivative and c value
p_ = sp.symbols('p'); c_ = sp.Rational(2)/sp.Float('1.0357712571566784', 30) - sp.Rational(4, 3)
Wp = 1 - p_ + p_**3/3
d_ = sp.Rational(1, 1000)
sigma = 2*Wp + d_*(1 + c_*p_)
out['sigma_prime'] = str(sp.simplify(sp.diff(sigma, p_) - (2*(p_**2 - 1) + d_*c_)))
out['c_value'] = str(sp.N(c_, 17))
out['sigma_f_at_phi_b'] = str(sp.N(sigma.subs(p_, sp.Float('0.9999159473169134', 20)), 15))
(HERE/'SYMPY_CHECKS.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
