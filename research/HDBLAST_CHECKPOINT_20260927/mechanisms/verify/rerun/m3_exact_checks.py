"""Mechanism screen, part 3: exact (symbolic) statements used by the screen.

 E1  The registered potential has the 'fake-supersymmetric' form: phi'=W', A'=-W/3 solves the flat
     static Einstein-scalar equations of the model (A'^2 = phi'^2/12 - U/6, phi''+4A'phi' = U').
 E2  The flat BPS wall between the bulk vacua phi=-1 (l=9/5) and phi=+1 (l=9) has effective tension
     Delta W = 4/3, exactly the critical value 3(1/l_- - 1/l_+) of a thin wall between them.
 E3  Thin-wall 'dark bubble' junction (non-Z2 wall, interior l_-, exterior l_+, normal pointing out):
       sqrt(1/l_-^2 + (H^2+k/a^2) - mu_-/a^4) - sqrt(1/l_+^2 + (H^2+k/a^2) - mu_+/a^4) = (lambda+rho)/3
     (kappa_5=1).  Expanded to first order in (H^2+k/a^2, mu, rho): coefficients of rho, mu_-, mu_+.
     Control: the Z2 limit reproduces the package's RS relation H^2 = (sigma+R)^2/36 - 1/l^2.
 E4  In the frozen-scalar phase (v=0, j=const) the package's Weyl identity gives Wy a^4 = const.
 E5  Sudden conversion at fixed sigma+R keeps nA continuous, hence H continuous; the Weyl term must
     then equal minus the radiation's own terms up to (H_before^2 - H_vac^2): exact algebra.
 E6  High-frequency (short-wave) partition for a dissipative scalar junction
       n.phi = -(sigma' + Y v)/2 ,  j = Y v / kappa_5^2 :
     an outgoing bulk wave carries 2/(2+Y) of the tension power, the shell matter Y/(2+Y).
"""
import sympy as sp
from common import dump

checks = {}; notes = {}
y = sp.symbols('y'); p = sp.Function('phi')(y)
Wf = lambda x: 1 - x + x**3/3
Uf = lambda x: sp.Rational(1, 2)*sp.diff(Wf(x), x)**2 - sp.Rational(2, 3)*Wf(x)**2
x = sp.symbols('x')
Wx = Wf(x); Ux = sp.expand(sp.Rational(1, 2)*sp.diff(Wx, x)**2 - sp.Rational(2, 3)*Wx**2)
# E1
phip = sp.diff(Wx, x); App = -Wx/3
constraint = sp.simplify(App**2 - (phip**2/12 - Ux/6))
# phi'' = d/dy W'(phi) = W''(phi) phi' = W'' W'
scalar = sp.simplify(sp.diff(phip, x)*phip + 4*App*phip - sp.diff(Ux, x))
checks['E1_BPS_constraint_identity'] = constraint == 0
checks['E1_BPS_scalar_identity'] = scalar == 0
# wrong-formula control: A' = -W/2 fails
checks['E1_control_wrong_A_fails'] = sp.simplify((-Wx/2)**2 - (phip**2/12 - Ux/6)) != 0
# E2
lm = 3/Wx.subs(x, -1); lp = 3/Wx.subs(x, 1)
checks['E2_l_minus_is_9_over_5'] = sp.nsimplify(lm) == sp.Rational(9, 5)
checks['E2_l_plus_is_9'] = sp.nsimplify(lp) == 9
DeltaW = Wx.subs(x, -1) - Wx.subs(x, 1)
crit = 3*(1/lm - 1/lp)
checks['E2_BPS_tension_equals_critical'] = sp.simplify(DeltaW - crit) == 0
notes['E2'] = dict(DeltaW=str(DeltaW), critical=str(sp.nsimplify(crit)),
                   U_minus=str(Ux.subs(x, -1)), U_plus=str(Ux.subs(x, 1)))
# E3
X, mum, mup, lam, rho, Lm, Lp, eps = sp.symbols('X mu_m mu_p lambda rho l_m l_p epsilon', real=True)
a = sp.symbols('a', positive=True)
lhs = sp.sqrt(1/Lm**2 + eps*X - eps*mum/a**4) - sp.sqrt(1/Lp**2 + eps*X - eps*mup/a**4)
ser = sp.series(lhs, eps, 0, 2).removeO()
eq = sp.Eq(ser.subs(eps, 1), (lam + rho)/3)
Xsol = sp.solve(eq, X)[0]
coef_rho = sp.simplify(sp.diff(Xsol, rho)); coef_mum = sp.simplify(sp.diff(Xsol, mum)); coef_mup = sp.simplify(sp.diff(Xsol, mup))
Lambda_term = sp.simplify(Xsol.subs({rho: 0, mum: 0, mup: 0}))
vals = {Lm: sp.Rational(9, 5), Lp: 9}
notes['E3'] = dict(X_first_order=str(sp.simplify(Xsol)), coef_rho=str(coef_rho), coef_mu_interior=str(coef_mum),
                   coef_mu_exterior=str(coef_mup), Lambda_term=str(Lambda_term),
                   model_values=dict(coef_rho=str(sp.nsimplify(coef_rho.subs(vals))), coef_mu_interior=str(coef_mum.subs(vals)),
                                     coef_mu_exterior=str(coef_mup.subs(vals)),
                                     Lambda_term=str(sp.simplify(Lambda_term.subs(vals))),
                                     lambda_critical=str(sp.solve(Lambda_term.subs(vals), lam)[0])))
checks['E3_brane_rho_coefficient_negative_for_l_plus_gt_l_minus'] = bool(coef_rho.subs(vals) < 0)
checks['E3_exterior_mass_positive_radiation'] = bool(coef_mup.subs({**vals, a: 1}) > 0)
checks['E3_interior_mass_negative_radiation'] = bool(coef_mum.subs({**vals, a: 1}) < 0)
checks['E3_critical_tension_equals_DeltaW'] = sp.solve(Lambda_term.subs(vals), lam)[0] == sp.Rational(4, 3)
# Z2 control: two identical interiors, normal toward the shell on each side -> 2 sqrt(1/l^2+X) = sigma_tot/3
L, sig, R = sp.symbols('l sigma R', positive=True)
Xz = sp.solve(sp.Eq(2*sp.sqrt(1/L**2 + X), (sig + R)/3), X)[0]
checks['E3_control_Z2_limit_matches_package'] = sp.simplify(Xz - ((sig + R)**2/36 - 1/L**2)) == 0
# E4 Weyl identity in frozen scalar phase
tau = sp.symbols('tau'); Hs = sp.Function('H')(tau); Wy = sp.Function('Wy')(tau)
ks, w = sp.symbols('k_s w'); v = 0
rhs = sp.Rational(1, 6)*(4*ks*w*v - 4*Hs*v**2 - 0 + w*0 - 0)
ode = sp.Eq(sp.diff(Wy, tau) + 4*Hs*Wy, rhs)
aa = sp.Function('a')(tau)
# check d(Wy a^4)/dtau = 0 using a' = H a
expr = sp.diff(Wy*aa**4, tau).subs(sp.diff(aa, tau), Hs*aa).subs(sp.diff(Wy, tau), -4*Hs*Wy)
checks['E4_Wy_a4_conserved_when_v_zero'] = sp.simplify(expr) == 0
# E5 sudden conversion
sf, Rr, Hb2, Hv2 = sp.symbols('sigma_f R H_b2 H_v2', positive=True)
Wy_req = Hb2 - ((sf + Rr)**2/36 - sf**2/36 + Hv2)     # uses H_vac^2 = sigma_f^2/36 - sigma_f'^2/48 + U_f/6
checks['E5_Wy_required_equals_minus_radiation_terms_plus_gap'] = sp.simplify(Wy_req - (-(sf*Rr/18 + Rr**2/36) + (Hb2 - Hv2))) == 0
# E6 impedance partition
Y, s1, vv = sp.symbols('Y sigma1 v', real=True)
# outgoing wave: n.phi = +v (phi = F(t+z), n = e^{-B} d_z, v = e^{-B} phi_t); junction n.phi = -(sigma1 + Y v)/2
vsol = sp.solve(sp.Eq(vv, -(s1 + Y*vv)/2), vv)[0]
tension_power = -s1*vsol           # -(d sigma/dtau) = -sigma' v  (kappa5=1)
matter_power = Y*vsol**2
bulk_power = 2*vsol*vsol           # both copies: -2 T_n tau = -2 (n.phi) v ... magnitude 2 v^2
checks['E6_power_balance'] = sp.simplify(tension_power - matter_power - bulk_power) == 0
frac = sp.simplify(matter_power/tension_power)
checks['E6_matter_fraction_is_Y_over_2_plus_Y'] = sp.simplify(frac - Y/(2 + Y)) == 0
notes['E6'] = dict(v=str(vsol), matter_fraction=str(frac),
                   validity='short-wavelength limit (frequency >> bulk curvature 1/l and >> H); the registered roll has '
                            'rate ~1.66 H0 << 1/l, so this is a known-limit calibration only')
dump('M3_EXACT_CHECKS.json', dict(status='exact-verified (SymPy %s)' % sp.__version__, checks=checks, notes=notes,
                                  all_passed=all(checks.values())))
for k, v in checks.items(): print(k, v)
print(notes)
assert all(checks.values())
