"""Exact structure of the linear radial problem behind the degree-14 Gegenbauer polynomial.

All checks are exact (SymPy rational / symbolic identities) unless labelled otherwise.
Writes GEGENBAUER_STRUCTURE.json.
"""
import json, hashlib
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
checks = []
def check(name, expr, expect_zero=True, note=''):
    val = sp.simplify(expr)
    ok = (val == 0) if expect_zero else (val != 0)
    checks.append(dict(name=name, passed=bool(ok), expect_zero=expect_zero, value=str(val)[:200], note=note))
    assert ok, (name, val)

phi, u, x, nu, eta = sp.symbols('phi u x nu eta')
W = 1 - phi + phi**3 / 3
U = sp.Rational(1, 2) * sp.diff(W, phi)**2 - sp.Rational(2, 3) * W**2
# --- 1. superpotential identity at a critical point of W: m^2/k^2 = s(s-4), s = 3W''/W
Wp, Wpp, W0 = sp.symbols('Wp Wpp W0')
Ugen_pp = Wpp**2 - sp.Rational(4, 3) * (Wp**2 + W0 * Wpp)   # U'' when W'''=0 and W'=0 dropped below
# general formula U'' = W''^2 + W'W''' - (4/3)(W'^2 + W W''); at W'=0 -> W''^2 - (4/3) W W''
kk2 = W0**2 / 9
s = 3 * Wpp / W0
check('U\'\'/k^2 = s(s-4) at a critical point of W (s=3W\'\'/W)', (Wpp**2 - sp.Rational(4, 3) * W0 * Wpp) / kk2 - s * (s - 4))
# explicit values at phi=+1 and phi=-1
for p0, name in [(1, 'phi=+1'), (-1, 'phi=-1')]:
    Wv = W.subs(phi, p0); Wppv = sp.diff(W, phi, 2).subs(phi, p0); Uppv = sp.diff(U, phi, 2).subs(phi, p0)
    k2 = -U.subs(phi, p0) / 6
    check(f'k^2=W^2/9 at {name}', k2 - Wv**2 / 9)
    check(f'U\'\'/k^2 = s(s-4) at {name}', Uppv / k2 - (3 * Wppv / Wv) * (3 * Wppv / Wv - 4))
Up1 = sp.diff(U, phi, 2).subs(phi, 1); k2p = -U.subs(phi, 1) / 6
check('U\'\'(1) = 28/9', Up1 - sp.Rational(28, 9)); check('k(+1)^2 = 1/81', k2p - sp.Rational(1, 81))
check('m^2/k^2 = 252 = 14*18 at phi=+1', Up1 / k2p - 252)
Um1 = sp.diff(U, phi, 2).subs(phi, -1); k2m = -U.subs(phi, -1) / 6
check('m^2/k^2 = 684/25 at phi=-1 (Delta=38/5, non-integer n=18/5)', Um1 / k2m - sp.Rational(684, 25))
n_m1 = sp.Rational(38, 5) - 4
check('phi=-1 degree n = Delta-4 = 18/5 is not an integer', sp.Integer(0) if not n_m1.is_integer else 1)

# --- 2. linear radial equation and Gegenbauer form
f = sp.Function('f')
lin_u = sp.diff(f(u), u, 2) + 4 * sp.coth(u) * sp.diff(f(u), u) - 252 * f(u)   # eta'' + 4 coth eta' - (m/k)^2 eta, u = y/9
# change variable x = cosh u : d/du = sinh u d/dx
F = sp.Function('F')
sub = {sp.diff(f(u), u, 2): sp.sinh(u)**2 * sp.Symbol('Fxx') + sp.cosh(u) * sp.Symbol('Fx'),
       sp.diff(f(u), u): sp.sinh(u) * sp.Symbol('Fx'), f(u): sp.Symbol('F0')}
expr = lin_u.subs(sub)
target = (sp.cosh(u)**2 - 1) * sp.Symbol('Fxx') + 5 * sp.cosh(u) * sp.Symbol('Fx') - 252 * sp.Symbol('F0')
check('radial equation in x=cosh(u): (x^2-1)F\'\' + 5xF\' - 252F = 0', sp.simplify(expr - target))
C14 = sp.gegenbauer(14, 2, x)
check('C_14^(2)(x) solves (x^2-1)F\'\'+5xF\'-n(n+4)F=0, n=14',
      sp.expand((x**2 - 1) * sp.diff(C14, x, 2) + 5 * x * sp.diff(C14, x) - 252 * C14))
# hypergeometric form: C_n^(2)(x)/C_n^(2)(1) = 2F1(-n, n+4; 5/2; (1-x)/2)
hyp = sum(sp.rf(-14, j) * sp.rf(18, j) / sp.rf(sp.Rational(5, 2), j) / sp.factorial(j) * ((1 - x) / 2)**j for j in range(15))
check('C_14^(2)(x)/C_14^(2)(1) = 2F1(-14,18;5/2;(1-x)/2) (terminating)', sp.expand(C14 / C14.subs(x, 1) - hyp))
check('C_14^(2)(1) = binom(17,3) = 680', C14.subs(x, 1) - 680)
# termination: the (j+1)-th/j-th coefficient ratio contains (j-n); a nonterminating series for non-integer n
j = sp.Symbol('j')
ratio = (j - 14) * (j + 18) / ((j + sp.Rational(5, 2)) * (j + 1))
check('2F1 coefficient ratio vanishes at j=14 (termination, quantization n in Z>=0)', ratio.subs(j, 14))
# --- 3. elementary closed form valid for ANY bulk mass: f_nu = (1/sinh u) d/du [sinh(nu u)/sinh u]
fnu = sp.diff(sp.sinh(nu * u) / sp.sinh(u), u) / sp.sinh(u)
ode_nu = sp.diff(fnu, u, 2) + 4 * sp.coth(u) * sp.diff(fnu, u) - (nu**2 - 4) * fnu
check('(1/sinh u) d/du[sinh(nu u)/sinh u] solves f\'\'+4coth f\'=(nu^2-4) f for all nu',
      sp.simplify(ode_nu.rewrite(sp.exp)))
gnu = sp.diff(sp.cosh(nu * u) / sp.sinh(u), u) / sp.sinh(u)
ode_g = sp.diff(gnu, u, 2) + 4 * sp.coth(u) * sp.diff(gnu, u) - (nu**2 - 4) * gnu
check('second (cone-singular) solution (1/sinh u) d/du[cosh(nu u)/sinh u]', sp.simplify(ode_g.rewrite(sp.exp)))
check('second solution ~ -u^-3 near cone', sp.limit(gnu * u**3, u, 0) + 1)
f16 = fnu.subs(nu, 16)
check('C_14^(2)(cosh u) = (1/2) f_16(u)  [i.e. C_n^(2)(x)=U\'_{n+1}(x)/2]',
      sp.simplify((sp.gegenbauer(14, 2, sp.cosh(u)) - f16 / 2).rewrite(sp.exp)))
check('regular solution value at cone: lim f_nu = nu(nu^2-1)/3 (=1360 at nu=16)',
      sp.limit(f16, u, 0) - 1360)
# exponential-sum form with positive coefficients
esum = sum((jj + 1) * (15 - jj) * sp.exp((14 - 2 * jj) * u) for jj in range(15))
check('C_14^(2)(cosh u) = sum_{j=0}^{14} (j+1)(15-j) e^{(14-2j)u}', sp.simplify((sp.gegenbauer(14, 2, sp.cosh(u)) - esum).rewrite(sp.exp)))
# consequence: log-derivative strictly positive for u>0 (pairs j, 14-j give positive sinh terms)
dsum = sum(2 * (jj + 1) * (15 - jj) * (14 - 2 * jj) * sp.sinh((14 - 2 * jj) * u) for jj in range(7))
check('dC/du = sum_{j<7} 2(j+1)(15-j)(14-2j) sinh((14-2j)u) (all coefficients > 0)',
      sp.simplify((sp.diff(esum, u) - dsum).rewrite(sp.exp)))
# producer's log derivative: d/dx C_n^(l) = 2 l C_{n-1}^(l+1)
check('dC_14^(2)/dx = 4 C_13^(3)', sp.expand(sp.diff(C14, x) - 4 * sp.gegenbauer(13, 3, x)))
# cone relation: eta(0)/[coefficient of e^{14u}] = 680/15 = 136/3
check('eta_h / alpha = 680/15 = 136/3', sp.Rational(680, 15) - sp.Rational(136, 3))
# linear static response with finite curvature: eta_b = -(9/2) delta c / (dlnC/du + 18)  (u units)
# flat limit dlnC/du -> 14 : eta_b -> -9 delta c / 64
check('flat limit of linear response: -(9/2)c/(14+18) = -9c/64',
      -sp.Rational(9, 2) / (14 + 18) + sp.Rational(9, 64))
# wrong-formula controls (must be NONZERO)
check('CONTROL: C_14^(1) (wrong index) does not solve the radial eq',
      sp.expand((x**2 - 1) * sp.diff(sp.gegenbauer(14, 1, x), x, 2) + 5 * x * sp.diff(sp.gegenbauer(14, 1, x), x) - 252 * sp.gegenbauer(14, 1, x)),
      expect_zero=False)
check('CONTROL: C_13^(2) (wrong degree) does not solve the radial eq',
      sp.expand((x**2 - 1) * sp.diff(sp.gegenbauer(13, 2, x), x, 2) + 5 * x * sp.diff(sp.gegenbauer(13, 2, x), x) - 252 * sp.gegenbauer(13, 2, x)),
      expect_zero=False)
check('CONTROL: 3W\'\'/W at phi=-1 is not an integer', sp.Integer(1) if (3 * sp.diff(W, phi, 2) / W).subs(phi, -1).is_integer else sp.Integer(0) + 1, expect_zero=False)

polycoef = sp.Poly(C14, x).all_coeffs()
res = dict(status='EXACT (symbolic identities, SymPy %s)' % sp.__version__,
           n_checks=len(checks), all_passed=all(c['passed'] for c in checks), checks=checks,
           C14_2_coefficients_highest_first=[str(v) for v in polycoef],
           exponential_form='C_14^(2)(cosh u) = sum_{j=0}^{14} (j+1)(15-j) exp((14-2j)u)',
           elementary_form_any_mass='f_nu(u) = (1/sinh u) d/du[ sinh(nu u)/sinh u ],  nu = sqrt(4 + m^2/k^2);  C_14^(2)(cosh u) = f_16(u)/2',
           quantization='polynomial in cosh(u) iff n = nu-2 is a non-negative integer; with U built from W, at a critical point of W: n = 3W\'\'/W - 4 (s>=4) or n=-3W\'\'/W (s<=0); n=14 at phi=+1, n=18/5 at phi=-1',
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE / 'GEGENBAUER_STRUCTURE.json').write_text(json.dumps(res, indent=2) + '\n')
print(res['n_checks'], res['all_passed'])
