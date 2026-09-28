"""Auditor check of the analytic statements (README section 3-4), written independently.
  1. energy identity  [rho^2 Z X]' = 3 lam Z^2 - rho^2 X^2/3 - 2 rho^2 phi'^2 Z^2
  2. Wronskian        [rho^2 (Z2 X1 - X2 Z1)]' = 3 (lam1 - lam2) Z1 Z2
  3. shell boundary terms with X_b = -3 lam Z_b/(B rho_b^2)
  4. growth criterion Re p > 0  <=>  Re mu2 < (Im mu2)^2/9 (sampled test, 20000 random points)
  5. decoupled AdS limit: mu2=0 radial equation in AdS5 (L=9, rho = 9 sinh(y/9), m^2 = 28/9) is the Gegenbauer eq. C_14^(2)
  6. U''(1) = 28/9, U'''(1) = 100/9, Delta = 18 ; B -> 2 + (Delta-4)/L = 32/9
Output: identities.json
"""
import json, random, cmath
from pathlib import Path
import sympy as sp
y = sp.symbols('y'); l, l1, l2, Bs, rb = sp.symbols('lambda lambda1 lambda2 B rho_b')
rho = sp.Function('rho')(y); ph = sp.Function('phi')(y); U1 = sp.Function('U1')(y)
H = rho.diff(y)/rho; g = (U1 - 4*H*ph.diff(y))/ph.diff(y)
def rhs(X, Z, lam): return (g*X + (3*lam/rho**2 - 2*ph.diff(y)**2)*Z, -(2*H + g)*Z - X/3)
X1, Z1, X2, Z2 = [sp.Function(n)(y) for n in ('X1', 'Z1', 'X2', 'Z2')]
d1, d2 = rhs(X1, Z1, l1), rhs(X2, Z2, l2)
sub = {X1.diff(y): d1[0], Z1.diff(y): d1[1], X2.diff(y): d2[0], Z2.diff(y): d2[1]}
out = {}
out['energy_identity'] = sp.simplify((rho**2*Z1*X1).diff(y).subs(sub) - (3*l1*Z1**2 - rho**2*X1**2/3 - 2*rho**2*ph.diff(y)**2*Z1**2)) == 0
out['wronskian_identity'] = sp.simplify((rho**2*(Z2*X1 - X2*Z1)).diff(y).subs(sub) - 3*(l1 - l2)*Z1*Z2) == 0
Zb, Zb1, Zb2 = sp.symbols('Z_b Z_b1 Z_b2')
out['energy_boundary_term'] = sp.simplify(rb**2*Zb*(-3*l*Zb/(Bs*rb**2)) + 3*l*Zb**2/Bs) == 0
out['wronskian_boundary_term'] = sp.simplify(rb**2*(Zb2*(-3*l1*Zb1/(Bs*rb**2)) - (-3*l2*Zb2/(Bs*rb**2))*Zb1) + 3*(l1 - l2)*Zb1*Zb2/Bs) == 0
# wrong-sign control: flipping the sign of X/3 in Z' must break the energy identity
d1b = (d1[0], -(2*H + g)*Z1 + X1/3)
out['CONTROL_energy_identity_fails_with_wrong_Z_sign'] = sp.simplify((rho**2*Z1*X1).diff(y).subs({X1.diff(y): d1b[0], Z1.diff(y): d1b[1]}) - (3*l1*Z1**2 - rho**2*X1**2/3 - 2*rho**2*ph.diff(y)**2*Z1**2)) != 0
random.seed(1); bad = 0
for _ in range(20000):
    m = complex(random.uniform(-50, 50), random.uniform(-50, 50))
    pgrow = (-1.5 + cmath.sqrt(2.25 - m)).real > 0
    crit = m.real < m.imag**2/9
    bad += (pgrow != crit)
out['growth_criterion_mismatches_of_20000'] = bad
f = sp.symbols('f'); W = 1 - f + f**3/3; U = sp.diff(W, f)**2/2 - sp.Rational(2, 3)*W**2
out['U(1)'] = str(U.subs(f, 1)); out['U2(1)'] = str(sp.diff(U, f, 2).subs(f, 1)); out['U3(1)'] = str(sp.diff(U, f, 3).subs(f, 1))
L = 1/(W.subs(f, 1)/3); out['AdS_radius_L'] = str(L)
m2L2 = sp.diff(U, f, 2).subs(f, 1)*L**2; Delta = 2 + sp.sqrt(4 + m2L2)
out['m2L2'] = str(m2L2); out['Delta'] = str(Delta); out['B_limit_2+(Delta-4)/L'] = str(2 + (Delta - 4)/L)
chi = sp.Function('chi'); x = sp.symbols('x')
# AdS: rho = L sinh(y/L); mu2=0: chi'' + 4 H chi' - m^2 chi = 0 ; with x = cosh(y/L) compare with Gegenbauer (1-x^2)C'' - 5 x C' + n(n+4) C = 0
# with x = cosh(y/L): d/dy = sinh/L d/dx and 4H = 4 cosh/(L sinh)  =>  (x^2-1)C'' + 5x C' - m^2 L^2 C = 0  (derived by hand;
# checked here by exact substitution on the y-form using a rational y-sample would need exp(); use the x-form exactly)
def geg_res(n):
    Cn = sp.gegenbauer(n, 2, x)
    return sp.expand((1 - x**2)*Cn.diff(x, 2) - 5*x*Cn.diff(x) + m2L2*Cn)
out['gegenbauer_C14_2_solves_AdS_mu2_0_exact_xform'] = geg_res(14) == 0
out['CONTROL_gegenbauer_C13_2_fails'] = geg_res(13) != 0
# chain-rule step checked symbolically for a generic function
Fx = sp.Function('F')
yy = sp.cosh(y/L)
lhs = sp.diff(Fx(yy), y, 2) + 4*sp.cosh(y/L)/(L*sp.sinh(y/L))*sp.diff(Fx(yy), y)
rhs_x = ((yy**2 - 1)*sp.Subs(sp.Derivative(Fx(x), x, 2), x, yy) + 5*yy*sp.Subs(sp.Derivative(Fx(x), x), x, yy))/L**2
out['chain_rule_x_form'] = sp.simplify((lhs - rhs_x).doit().rewrite(sp.exp)) == 0
out = {k: (bool(v) if isinstance(v, bool) else v) for k, v in out.items()}
Path(__file__).with_name('identities.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=1))
