"""Spherically symmetric (l = 0) sector of P_H and its explicit inversion.

Ansatz h0 = a(s) x x^T, psi0 = Psi(s).  The engine gives
   P_H(h0,psi0) = ca1 a' + ca0 a + cp1 Psi' + cp0 Psi     (' = d/ds).
Exact identities verified here (using the background equations):
 (I1) cp1 Psi' + cp0 Psi = -4 s^{-3/2} Omega^{-5} d/ds( s^{5/2} Omega^3 Phi' Psi )
 (I2) for every a:  s^{3/2} Omega^5 (ca1 a' + ca0 a)
          = d/ds( 3 s^{5/2}(Omega + 2 s Omega') a ) - (3/2) s^{3/2} Omega a
Consequently, for a compactly supported l = 0 source f0(s) on [s1,s2]:
   choose a = kappa*beta(s) with a bump beta >= 0, beta != 0, and
   kappa = -(2/3) int s^{3/2} Omega^5 f0 ds / int s^{3/2} Omega beta ds ;
   then Psi(s) = -(4 s^{5/2} Omega^3 Phi')^{-1} int_{s1}^{s} sig^{3/2} Omega^5
                 (f0 - kappa P_a beta) dsig
is smooth, vanishes outside [s1,s2], and P_H(h0,psi0) = f0 exactly.
"""
import json
import sympy as sp
from engine import X, O, F, register, reduce_bg, bg_rules, dS
from linear_ops import PH_PM

Aa = register('Aa'); Ps = register('Ps')
xv = sp.Matrix(X)
DH, _ = PH_PM(h=Aa[0] * (xv * xv.T), psi=Ps[0])
s = sp.Symbol('s_')
e = sp.together(reduce_bg(DH))
e = sp.simplify(e.subs({X[0]: 0, X[1]: 0, X[2]: 0, X[3]: sp.sqrt(s)}))
jv = (Aa[0], Aa[1], Aa[2], Ps[0], Ps[1], Ps[2])
coef = {str(v): sp.factor(sp.diff(e, v)) for v in jv}
lin_ok = sp.simplify(e - sum(sp.diff(e, v) * v for v in jv)) == 0
rules, _ = bg_rules()
O0, O1, F1 = O[0], O[1], F[1]
I1 = sp.simplify(coef['Ps1'] * Ps[1] + coef['Ps0'] * Ps[0]
                 + 4 * s**sp.Rational(-3, 2) / O0**5 * dS(s**sp.Rational(5, 2) * O0**3 * F1 * Ps[0], rules))
I2 = sp.simplify(s**sp.Rational(3, 2) * O0**5 * (coef['Aa1'] * Aa[1] + coef['Aa0'] * Aa[0])
                 - dS(3 * s**sp.Rational(5, 2) * (O0 + 2 * s * O1) * Aa[0], rules)
                 + sp.Rational(3, 2) * s**sp.Rational(3, 2) * O0 * Aa[0])
out = {k: str(v) for k, v in coef.items()}
out['linear_first_order_in_a_and_Psi'] = bool(lin_ok and coef['Aa2'] == 0 and coef['Ps2'] == 0)
out['I1_psi_total_derivative'] = bool(I1 == 0)
out['I2_a_bracket_minus_three_halves'] = bool(I2 == 0)
print(json.dumps(out, indent=2))
json.dump(out, open('../results/l0_sector.json', 'w'), indent=2)
