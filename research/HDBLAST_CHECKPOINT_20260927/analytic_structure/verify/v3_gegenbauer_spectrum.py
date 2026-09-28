"""Auditor's independent checks of the Gegenbauer structure and the leading-order spectrum claims.

Exact parts use sympy; spectral parts use DIRECT ODE shooting (mpmath.odefun / scipy) and do not use the
Legendre-function formulas of the audited script, except where a formula itself is being tested.
Writes verify/V3_GEGENBAUER_SPECTRUM.json
"""
import json, math, time
from pathlib import Path
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
t0 = time.time()
res = {}
checks = []
def chk(name, ok, val=''):
    checks.append(dict(name=name, passed=bool(ok), value=str(val)[:160]))

# ---------------- exact: superpotential numbers
ph = sp.Symbol('phi')
W = 1 - ph + ph**3 / 3
U = sp.Rational(1, 2) * sp.diff(W, ph)**2 - sp.Rational(2, 3) * W**2
for p0 in (1, -1):
    k2 = -U.subs(ph, p0) / 6
    m2 = sp.diff(U, ph, 2).subs(ph, p0) / k2
    s = 3 * sp.diff(W, ph, 2).subs(ph, p0) / W.subs(ph, p0)
    Delta = 2 + sp.sqrt(4 + m2)
    chk(f'phi={p0}: k^2 = W^2/9', sp.simplify(k2 - W.subs(ph, p0)**2 / 9) == 0, k2)
    chk(f'phi={p0}: m^2/k^2 = s(s-4)', sp.simplify(m2 - s * (s - 4)) == 0, (m2, s))
    res[f'phi_{p0}'] = dict(m2_over_k2=str(m2), s=str(s), Delta=str(sp.nsimplify(Delta)), n=str(sp.nsimplify(Delta - 4)))
chk('Delta(+1)=18, n=14', res['phi_1']['Delta'] == '18' and res['phi_1']['n'] == '14', res['phi_1'])
chk('Delta(-1)=38/5, n=18/5 (s=-18/5 = 4-Delta there, i.e. s is the OTHER root)', res['phi_-1']['Delta'] == '38/5', res['phi_-1'])
# junction Robin coefficient b = sigma''/(2k) at leading order sigma = 2W: = W''/k = 3W''/W
b = sp.diff(2 * W, ph, 2).subs(ph, 1) / 2 * 9
chk('scalar Robin coefficient: f_u = -(9/2) sigma\'\' f = -18 f', b == 18, b)

# ---------------- exact: Gegenbauer closed forms
x, u, nu = sp.symbols('x u nu', positive=True)
C14 = sp.gegenbauer(14, 2, x)
expsum = sum((j + 1) * (15 - j) * sp.exp((14 - 2 * j) * u) for j in range(15))
chk('C14^(2)(cosh u) = sum (j+1)(15-j) e^{(14-2j)u}', sp.simplify(sp.expand((C14.subs(x, sp.cosh(u))).rewrite(sp.exp)) - expsum) == 0)
chk('C14^(2)(1) = 680', C14.subs(x, 1) == 680)
U15 = sp.chebyshevu(15, x)
chk('C14^(2) = (1/2) dU15/dx', sp.expand(C14 - sp.diff(U15, x) / 2) == 0)
chk('U15(cosh u) = sinh(16u)/sinh(u)', sp.simplify((U15.subs(x, sp.cosh(u)) * sp.sinh(u) - sp.sinh(16 * u)).rewrite(sp.exp)) == 0)
poly_claim = 245760*x**14 - 745472*x**12 + 878592*x**10 - 506880*x**8 + 147840*x**6 - 20160*x**4 + 1008*x**2 - 8
chk('explicit polynomial of C14^(2) as stated', sp.expand(C14 - poly_claim) == 0)
# generic nu: f_nu = (1/sinh u) d/du[sinh(nu u)/sinh u] solves f''+4coth f' = (nu^2-4) f
fnu = sp.diff(sp.sinh(nu * u) / sp.sinh(u), u) / sp.sinh(u)
r = sp.diff(fnu, u, 2) + 4 * sp.cosh(u) / sp.sinh(u) * sp.diff(fnu, u) - (nu**2 - 4) * fnu
vals = [abs(complex(sp.N(r.subs({nu: nv, u: uv}), 30))) for nv in (sp.Rational(7, 3), 16, sp.Rational(51, 10)) for uv in (sp.Rational(1, 3), 2)]
chk('f_nu solves the radial equation (numerical substitution at 6 points, 30 digits)', max(vals) < 1e-20, max(vals))
fsing = sp.diff(sp.cosh(nu * u) / sp.sinh(u), u) / sp.sinh(u)
r2 = sp.diff(fsing, u, 2) + 4 * sp.cosh(u) / sp.sinh(u) * sp.diff(fsing, u) - (nu**2 - 4) * fsing
vals2 = [abs(complex(sp.N(r2.subs({nu: nv, u: uv}), 30))) for nv in (sp.Rational(7, 3), 16) for uv in (sp.Rational(1, 3), 2)]
chk('singular partner also solves it', max(vals2) < 1e-20, max(vals2))
lead = sp.limit(fsing.subs(nu, 5) * u**3, u, 0)
chk('singular partner ~ -u^-3 at the cone', lead == -1, lead)
chk('f_16 = 2 C14^(2)(cosh u)', sp.simplify((fnu.subs(nu, 16) - 2 * C14.subs(x, sp.cosh(u))).rewrite(sp.exp)) == 0)
# wrong-formula control: nu=15 (n=13) does NOT give the m^2=252 solution
chk('control: f_15 does not solve the m^2/k^2=252 equation', abs(complex(sp.N((sp.diff(fnu, u, 2) + 4 * sp.cosh(u) / sp.sinh(u) * sp.diff(fnu, u) - 252 * fnu).subs({nu: 15, u: 1}), 20))) > 1)
# d ln C/du > 0 for u>0: C = sum over pairs 2 a_j cosh((14-2j)u) + a_7 with a_j>0
dC = sp.diff(expsum, u)
gvals = [float(sp.N((dC / expsum).subs(u, uv))) for uv in (0.01, 0.1, 0.5, 1, 3.36)]
chk('d ln C/du > 0 at sample points', min(gvals) > 0, gvals)
# linear static response eta_b = -(9/2) delta c/(g+18) -> -(9/64) c delta at large u (g -> 14)
chk('large-u limit of static response = -9c delta/64', sp.Rational(9, 2) / (14 + 18) == sp.Rational(9, 64))

# ---------------- tensor: direct ODE test of  f'(u_b) = K(u_b,mu) (mu - 3/2) P^{-mu}_{1/2}(cosh u_b)
mp.mp.dps = 25
def shoot_mp(M2, ub, m5sq, u0=mp.mpf('1e-3')):
    mu = mp.sqrt(mp.mpf(9) / 4 - M2)
    # Frobenius leading term u^{mu-3/2} (1 + a u^2): a from the indicial recursion of f''+4/u f'+(M2/u^2)f with coth, sinh corrections
    # use exact small-u data from the Legendre-free series: f = u^s (1 + a2 u^2), s = mu - 3/2
    s = mu - mp.mpf(3) / 2
    # coefficient a2: plug into f''+4coth(u)f'+(M2/sinh^2 u - m5sq)f = 0 with coth = 1/u + u/3, 1/sinh^2 = 1/u^2 - 1/3
    a2 = (m5sq + M2 / 3 - 4 * s / 3) / ((s + 2) * (s + 1) + 4 * (s + 2) + M2)
    f0 = u0**s * (1 + a2 * u0**2); f1 = s * u0**(s - 1) + a2 * (s + 2) * u0**(s + 1)
    sol = mp.odefun(lambda uu, y: [y[1], -4 * mp.coth(uu) * y[1] - (M2 / mp.sinh(uu)**2 - m5sq) * y[0]], u0, [f0, f1])
    y = sol(ub)
    return y[1] / y[0]
def Pm(nu_, mu_, xx): return mp.legenp(nu_, -mu_, xx, type=3)
trows = []
for ub in [mp.mpf('0.5'), mp.mpf('2.0'), mp.mpf('3.364')]:
    for M2 in [mp.mpf('0.3'), mp.mpf('1.2'), mp.mpf('2.0')]:
        mu = mp.sqrt(mp.mpf(9) / 4 - M2); xx = mp.cosh(ub)
        ld = shoot_mp(M2, ub, 0)
        # predicted: f'/f = (mu - 3/2) P_{1/2}/(sinh(u) P_{3/2})  [from (x^2-1)P'_nu = nu x P_nu - (nu-mu) P_{nu-1}]
        pred = (mu - mp.mpf(3) / 2) * Pm(mp.mpf(1) / 2, mu, xx) / (mp.sinh(ub) * Pm(mp.mpf(3) / 2, mu, xx))
        wrong = (mu - mp.mpf(1) / 2) * Pm(-mp.mpf(1) / 2, mu, xx) / (mp.sinh(ub) * Pm(mp.mpf(1) / 2, mu, xx))
        trows.append(dict(u_b=float(ub), M2=float(M2), shoot=float(ld), pred=float(pred), abs_err=float(abs(ld - pred)),
                          wrong_nu_err=float(abs(ld - wrong))))
res['tensor_condition_vs_direct_ODE'] = trows
chk('tensor: f\'/f from direct ODE equals (mu-3/2)P_{1/2}/(sinh P_{3/2})', max(r_['abs_err'] for r_ in trows) < 1e-6, max(r_['abs_err'] for r_ in trows))
chk('tensor control (wrong degree) fails', min(r_['wrong_nu_err'] for r_ in trows) > 1e-3, min(r_['wrong_nu_err'] for r_ in trows))
# positivity of P^{-mu}_{1/2}(x), x>1, independent wide sweep (mu up to 12, x up to cosh 12)
mn = None; bad = 0
for xx in [1.0001, 1.01, 1.2, 2, 5, 30, 1e3, 1e5, float(mp.cosh(12))]:
    for mu in np.linspace(0.005, 12, 120):
        v = Pm(mp.mpf(1) / 2, mp.mpf(mu), mp.mpf(xx))
        if v <= 0: bad += 1
        nv = v * mp.gamma(1 + mu) / ((mp.mpf(xx) - 1) / (mp.mpf(xx) + 1))**(mu / 2)
        mn = nv if mn is None else min(mn, nv)
res['P_minus_mu_half_positivity'] = dict(nonpositive_count=bad, min_normalised=float(mn))
chk('P^{-mu}_{1/2}(x) > 0 on sweep', bad == 0, (bad, float(mn)))
# direct-ODE sign scan of tensor f'(u_b) for M^2 in (0, 9/4) (scipy float, independent grid)
def shoot_f(M2, ub, m5sq, u0=1e-4):
    mu = math.sqrt(9 / 4 - M2); s = mu - 1.5
    a2 = (m5sq + M2 / 3 - 4 * s / 3) / ((s + 2) * (s + 1) + 4 * (s + 2) + M2)
    y0 = [1 + a2 * u0**2, s / u0 + a2 * (s + 2) * u0]   # f/u0^s normalisation (log-derivative unaffected)
    # integrate g = f / u0^s scaled; use log-derivative Riccati-free linear system
    sol = solve_ivp(lambda uu, y: [y[1], -4 / math.tanh(uu) * y[1] - (M2 / math.sinh(uu)**2 - m5sq) * y[0]],
                    (u0, ub), y0, method='DOP853', rtol=1e-11, atol=1e-14)
    return sol.y[1, -1] / sol.y[0, -1]
tscan = []
for ub in [0.02, 0.2, 1.0, 3.364, 6.0]:
    Fv = [shoot_f(m, ub, 0.0) for m in np.linspace(0.01, 2.24, 70)]
    tscan.append(dict(u_b=ub, sign_changes=int(sum(a * b_ < 0 for a, b_ in zip(Fv, Fv[1:]))), max=max(Fv), min=min(Fv)))
res['tensor_direct_scan'] = tscan
chk('tensor: no sign change of f\'(u_b) for 0<M^2<9/4 (direct ODE)', all(t['sign_changes'] == 0 for t in tscan), tscan)

# ---------------- decoupled scalar: threshold u_b where F + 18 = 0 has a root, by direct ODE (scipy)
def Fs(M2, ub): return shoot_f(M2, ub, 252.0) + 18
# at M^2 -> 9/4 (use 2.2499) find critical u_b
g = lambda ub: Fs(2.2499, ub)
grid = [0.03, 0.04, 0.05, 0.055, 0.06, 0.065, 0.07, 0.1]
gv = [g(v) for v in grid]
ucrit = None
for a, b_, ga, gb in zip(grid, grid[1:], gv, gv[1:]):
    if ga * gb < 0:
        ucrit = brentq(g, a, b_, xtol=1e-10); break
res['scalar_threshold'] = dict(grid=grid, F_plus_b=gv, u_b_crit=ucrit, H_over_k_crit=(1 / math.sinh(ucrit) if ucrit else None))
chk('scalar: critical H/k ~ 16.6 (audited: u_b=0.0600, H/k=16.65)', ucrit is not None and abs(1 / math.sinh(ucrit) - 16.65) < 0.2,
    res['scalar_threshold']['H_over_k_crit'])
# bound state at H/k = 50
ub50 = math.asinh(1 / 50)
m2grid = np.linspace(0.01, 2.24, 120); fv = [Fs(m, ub50) for m in m2grid]
roots = [brentq(lambda m: Fs(m, ub50), a, b_) for a, b_, fa, fb in zip(m2grid, m2grid[1:], fv, fv[1:]) if fa * fb < 0]
res['scalar_H_over_k_50'] = dict(roots_M2=roots)
chk('scalar: bound state M^2 ~ 1.004 at H/k=50', len(roots) == 1 and abs(roots[0] - 1.004) < 0.01, roots)
# registered shell u_b ~ 3.36 and u_b=1: F+b positive, no root
for ub in [1.0, 3.364]:
    fv = [Fs(m, ub) for m in np.linspace(0.01, 2.24, 60)]
    res[f'scalar_u_b_{ub}'] = dict(min=min(fv), max=max(fv))
    chk(f'scalar: no root and F+b>0 at u_b={ub}', min(fv) > 0, (min(fv), max(fv)))
# tachyonic side M^2 < 0: f'/f > 0
fv = [shoot_f(m, 1.0, 252.0) for m in [-0.5, -2, -10]]
chk('scalar: M^2<0 gives f\'/f>0 at u_b=1', min(fv) > 0, fv)

res['checks'] = checks
res['all_passed'] = all(c_['passed'] for c_ in checks)
res['runtime_s'] = time.time() - t0
(HERE / 'V3_GEGENBAUER_SPECTRUM.json').write_text(json.dumps(res, indent=1, default=str) + '\n')
for c_ in checks: print(c_['passed'], c_['name'], c_['value'][:100])
print('all passed', res['all_passed'], 'runtime', res['runtime_s'])
