"""Audit V1: independent symbolic re-derivations of the formulas the mechanism screen relies on.

Method (independent of the workstream's m3 script): start from the Shiromizu-Maeda-Sasaki (SMS) projected
Einstein equation on a Z2 brane, kappa_5 = 1, signature (-,+,+,+,+), bulk action R/2 - (1/2)(d phi)^2 - U,
  G^mu_nu(4) = (2/3)[T^mu_nu + (T_nn - T/4) delta^mu_nu] + K K^mu_nu - K^mu_a K^a_nu - (1/2) delta (K^2 - K.K) - E^mu_nu,
with brane-frame components v = u.d phi, w = n.d phi, and mixed extrinsic curvature K^tau_tau = k_t, K^i_j = k_s delta.
Junction (Z2, S = -sigma q + radiation R, p = R/3):  K_mu nu = -(1/2)(S_mu nu - S q/3)
   -> k_s = -(sigma+R)/6, k_t = -(sigma-3R)/6  (package's nA = (sigma+R)/6, nB = (sigma-3R)/6 up to orientation sign).
Checks:
 D1  H^2 identity  H^2 = (sigma+R)^2/36 + v^2/12 - w^2/12 + U/6 + Wy,  Wy = -E^tau_tau... (sign fixed below)
 D2  effective pressure and the Weyl transport identity with shell radiation and friction source
       dWy/dtau + 4 H Wy = (1/6)[4 k w v - 4 H v^2 - v vdot + w wdot - U' v],  k = (sigma+R)/6, w = -(sigma'+Y v)/2,
     given dR/dtau + 4HR = Y v^2 and dsigma/dtau = sigma' v.  Control: with the wrong ledger sign it fails.
 D3  Codazzi energy balance: d(sigma+R)/dtau + 4HR = -2 w v  <=>  dR/dtau + 4HR = Y v^2  (for w = -(sigma'+Yv)/2).
 D4  BPS structure: for ANY superpotential with l = 3/W at the vacua, W(-1)-W(+1) = 3(1/l_- - 1/l_+)  (tautology check),
     and U''(-1) (bulk mass at phi=-1).
 D5  Dark-bubble junction: exact vs linearised Lambda_4 at lambda = 0 and near lambda_crit (validity of '10/81 - 5 lambda/54').
 D6  Leading H_vac^2 at the '+1' vacuum for sigma = 2W + delta*eps(phi): delta*eps(1)/27 + delta^2 eps(1)^2/36 at phi_b = 1
     (the -B^2/384 term needs the scalar profile and is not re-derived here).
 D7  At c = 0 the constant phi = -1 bulk solves both static junctions (U'(-1) = 0, sigma'(-1) = delta c = 0):
     the initial-shell family meets the trivial solution there, so a continuation parameterised by log10(phi_h + 1)
     cannot pass c = 0.
Output: V1_INDEPENDENT_DERIVATIONS.json
"""
import json, sympy as sp
from pathlib import Path

HERE = Path(__file__).resolve().parent
res = {}; notes = {}

v, w, U, sig, R, H, E = sp.symbols('v w U sigma R H E', real=True)
# 5D stress tensor components in the orthonormal (tau, n) frame on the brane; only tau and n gradients.
grad2 = -v**2 + w**2
T_tt_mixed = -(v**2) - (grad2/2 + U)        # T^tau_tau = g^tt T_tt, T_tt = v^2 + (grad2/2+U)
T_ii_mixed = -(grad2/2 + U)                  # T^i_i
T_nn = w**2 - (grad2/2 + U)
T = grad2 - 5*(grad2/2 + U)
kt = -(sig - 3*R)/6; ks = -(sig + R)/6
Ktr = kt + 3*ks; KK = kt**2 + 3*ks**2
G_tt_mixed = sp.Rational(2, 3)*(T_tt_mixed + T_nn - T/4) + Ktr*kt - kt**2 - (Ktr**2 - KK)/2 - E
G_ii_mixed = sp.Rational(2, 3)*(T_ii_mixed + T_nn - T/4) + Ktr*ks - ks**2 - (Ktr**2 - KK)/2 - E/3*(-1)  # E traceless: E^i_i = -E^t_t/3
# FRW: G^tau_tau = -3H^2, G^i_i = -(2 Hdot + 3 H^2)
H2 = sp.solve(sp.Eq(G_tt_mixed, -3*H**2), H**2)[0]
Wy = sp.symbols('Wy')
H2_pkg = (sig + R)**2/36 + v**2/12 - w**2/12 + U/6 + Wy
# identify Wy = E/3 (E = E^tau_tau)
res['D1_H2_identity_matches_package_with_Wy_eq_Etautau_over_3'] = sp.simplify(H2.subs(E, 3*Wy) - H2_pkg) == 0
# control: using the wrong junction normalisation K = -(S - S q/3) (factor 2) must fail
ks_w, kt_w = 2*ks, 2*kt
G_tt_wrong = sp.Rational(2, 3)*(T_tt_mixed + T_nn - T/4) + (kt_w+3*ks_w)*kt_w - kt_w**2 - ((kt_w+3*ks_w)**2 - (kt_w**2+3*ks_w**2))/2 - E
res['D1_control_wrong_junction_factor_fails'] = sp.simplify(sp.solve(sp.Eq(G_tt_wrong, -3*H**2), H**2)[0].subs(E, 3*Wy) - H2_pkg) != 0

# D2: effective 4D fluid (brane-local part) rho_eff = -G^t_t|_{E=0}, p_eff = G^i_i|_{E=0}; Weyl fluid rho_W = 3 Wy, p_W = Wy.
rho_eff = -(G_tt_mixed.subs(E, 0)); p_eff = G_ii_mixed.subs(E, 0)
notes['rho_eff'] = str(sp.factor(sp.expand(rho_eff))); notes['p_eff'] = str(sp.factor(sp.expand(p_eff)))
res['D2_rho_plus_p_eff'] = str(sp.factor(sp.expand(rho_eff + p_eff)))
tau = sp.symbols('tau'); Y, Hs = sp.symbols('Y H', real=True)
phi = sp.Function('phi')(tau); Rf = sp.Function('R')(tau); Wyf = sp.Function('Wy')(tau); wf = sp.Function('w')(tau)
Sg = sp.Function('Sigma'); Uf = sp.Function('Uf')
vf = sp.diff(phi, tau)
rho_eff_t = rho_eff.subs({sig: Sg(phi), R: Rf, v: vf, w: wf, U: Uf(phi)})
p_eff_t = p_eff.subs({sig: Sg(phi), R: Rf, v: vf, w: wf, U: Uf(phi)})
# total conservation (4D Bianchi, E traceless): d(rho_eff + 3Wy)/dtau + 3H(rho_eff + p_eff + 4Wy) = 0
cons = sp.diff(rho_eff_t, tau) + 3*Hs*(rho_eff_t + p_eff_t) + 3*(sp.diff(Wyf, tau) + 4*Hs*Wyf)
def check_transport(ledger_sign):
    Rdot = -4*Hs*Rf + ledger_sign*Y*vf**2
    dWy = sp.solve(cons.subs(sp.diff(Rf, tau), Rdot), sp.diff(Wyf, tau))[0]
    k = (Sg(phi) + Rf)/6
    wexpr = -(sp.diff(Sg(phi), phi) + Y*vf)/2
    rhs_pkg = (4*k*wf*vf - 4*Hs*vf**2 - vf*sp.diff(vf, tau) + wf*sp.diff(wf, tau) - sp.diff(Uf(phi), phi)*vf)/6 - 4*Hs*Wyf
    diff = sp.simplify((dWy - rhs_pkg).subs(wf, wexpr).doit())
    return diff
d_ok = check_transport(+1); d_bad = check_transport(-1)
res['D2_weyl_transport_with_matter_matches_package_form'] = d_ok == 0
res['D2_control_wrong_ledger_sign_fails'] = d_bad != 0
# D3 Codazzi: brane fluid rho_b = sigma + R, p_b = -sigma + R/3 ; energy balance d rho_b + 3H(rho_b+p_b) = -2 w v
sp1 = sp.symbols('sigma1', real=True)
Rdot_from_codazzi = sp.solve(sp.Eq(sp1*v + sp.Symbol('Rdot') + 4*H*R, -2*(-(sp1 + Y*v)/2)*v), sp.Symbol('Rdot'))[0]
res['D3_codazzi_gives_ledger_Rdot_plus_4HR_eq_Yv2'] = sp.simplify(Rdot_from_codazzi - (-4*H*R + Y*v**2)) == 0

# D4 BPS tautology and bulk mass at phi=-1
x = sp.symbols('x'); Wx = 1 - x + x**3/3; Ux = sp.Rational(1, 2)*sp.diff(Wx, x)**2 - sp.Rational(2, 3)*Wx**2
Wm, Wp = sp.symbols('W_m W_p', positive=True)
res['D4_DeltaW_eq_critical_is_identity_for_any_W'] = sp.simplify((Wm - Wp) - 3*(1/(3/Wm) - 1/(3/Wp))) == 0
res['D4_U_prime_at_minus1_is_zero'] = sp.simplify(sp.diff(Ux, x).subs(x, -1)) == 0
Upp = sp.nsimplify(sp.diff(Ux, x, 2).subs(x, -1)); notes['U_pp_at_minus1'] = str(Upp)
notes['m2_l2_at_minus1'] = str(sp.nsimplify(Upp*sp.Rational(9, 5)**2))

# D5 dark bubble exact vs linear
lam, X = sp.symbols('lambda X', real=True)
lm, lp = sp.Rational(9, 5), 9
def X_exact(lv):
    f = sp.sqrt(1/lm**2 + X) - sp.sqrt(sp.Rational(1, lp**2) + X) - sp.Rational(lv)/3
    try:
        return float(sp.nsolve(f, X, 0.001))
    except Exception as e:
        return None
lin = lambda lv: float(sp.Rational(10, 81) - sp.Rational(5, 54)*sp.Rational(lv))
d5 = {}
for lv in ['4/3', '13/10', '6/5', '1', '1/2']:
    d5[lv] = dict(exact=X_exact(lv), linear=lin(lv))
# large-X limit of LHS: sqrt(a+X)-sqrt(b+X) -> 0, so small lambda needs X -> infinity; exact solution exists for all 0<lambda<4/3
d5['note'] = 'X = H^2 (flat slicing). Linear formula valid only for |X| << 1/l_+^2 = 1/81, i.e. lambda close to 4/3.'
res['D5_linear_coefficient_matches_exact_near_critical'] = abs((d5['13/10']['exact'] - 0)/(4/3 - 1.3) - 5/54) / (5/54) < 0.2 if d5['13/10']['exact'] else False
notes['D5'] = d5

# D6 leading H_vac^2 at phi_b = 1 (no scalar gradient): sigma = 2/3 + delta*A
dl, A = sp.symbols('delta A', positive=True)
H2v = sp.expand(((sp.Rational(2, 3) + dl*A)**2/36 + Ux.subs(x, 1)/6))
res['D6_Hvac2_leading'] = sp.simplify(H2v - (dl*A/27 + dl**2*A**2/36)) == 0

# D7 trivial solution at c = 0
c = sp.symbols('c', real=True)
sig1_at_m1 = (2*sp.diff(Wx, x) + dl*c).subs(x, -1)
res['D7_sigma_prime_at_minus1_is_delta_c'] = sp.simplify(sig1_at_m1 - dl*c) == 0
notes['D7'] = ('phi == -1 is an exact bulk solution (U\'(-1)=0) and satisfies the scalar junction iff sigma\'(-1)=delta*c=0; '
               'so at c=0 the shell family meets the constant solution; the linear response phi_h+1 changes sign with c.')

out = dict(status='exact-verified (SymPy %s) unless marked' % sp.__version__, checks={k: (bool(v) if isinstance(v, (bool, sp.logic.boolalg.BooleanAtom)) else v) for k, v in res.items()}, notes=notes)
(HERE/'V1_INDEPENDENT_DERIVATIONS.json').write_text(json.dumps(out, indent=2, default=str) + '\n')
print(json.dumps(out, indent=1, default=str))
