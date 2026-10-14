"""Full l = 0 sector of the linearized Hamiltonian constraint (original engine).
h0 = alpha(s) delta + a(s) x x^T,  psi0 = Psi(s).  Radial gauge xi = zeta(s) x gives
alpha_g = 2 Omega^2 zeta + 4 s Omega Omega' zeta, a_g = 4 Omega^2 zeta', Psi_g = 2 s Phi' zeta.
Checks: (1) gauge directions are annihilated; (2) the weighted functional
L = int s^{3/2} Omega^5 P_H ds reduces (mod total s-derivatives) to
int [c_alpha alpha + c_a a] ds with no Psi term; reports c_alpha, c_a."""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'basis', 'code'))
import sympy as sp
from engine import X, O, F, register, reduce_bg, bg_rules, dS
from linear_ops import PH_PM

Al = register('Al'); Aa = register('Aa'); Ps = register('Ps'); Z = register('Z')
s = sp.Symbol('s_')
xv = sp.Matrix(X)
rep = {X[0]: 0, X[1]: 0, X[2]: 0, X[3]: sp.sqrt(s)}
def radial(e):
    return sp.simplify(sp.together(reduce_bg(e)).subs(rep))
DH, _ = PH_PM(h=Al[0] * sp.eye(4) + Aa[0] * (xv * xv.T), psi=Ps[0])
e = radial(DH)
jv = [Al[0], Al[1], Al[2], Aa[0], Aa[1], Aa[2], Ps[0], Ps[1], Ps[2]]
coef = {str(v): sp.factor(sp.diff(e, v)) for v in jv}
lin = sp.simplify(e - sum(sp.diff(e, v) * v for v in jv)) == 0
rules, _ = bg_rules()
O0, O1, F1 = O[0], O[1], F[1]
# gauge check
zeta = Z[0]
gsub = {Al[0]: 2 * O0**2 * zeta + 4 * s * O0 * O1 * zeta,
        Aa[0]: 4 * O0**2 * Z[1], Ps[0]: 2 * s * F1 * zeta}
def jets_of(expr0):
    j0 = expr0; j1 = dS(j0, rules); j2 = dS(j1, rules)
    return j0, j1, j2
ga = jets_of(gsub[Al[0]]); gaa = jets_of(gsub[Aa[0]]); gp = jets_of(gsub[Ps[0]])
gsubs = {Al[0]: ga[0], Al[1]: ga[1], Al[2]: ga[2], Aa[0]: gaa[0], Aa[1]: gaa[1], Aa[2]: gaa[2],
         Ps[0]: gp[0], Ps[1]: gp[1], Ps[2]: gp[2]}
gauge_res = sp.simplify(sp.together(e.subs(gsubs, simultaneous=True).subs(rules)))
# weighted functional: W = s^{3/2} Omega^5 * e ; integrate by parts symbolically:
# write W = sum_k w_k(s) * f^{(k)} for f in {Al, Aa, Ps}; the reduced density is
# sum_k (-1)^k d^k w_k / ds^k
w = s**sp.Rational(3, 2) * O0**5
red = {}
for nm, J in (('alpha', Al), ('a', Aa), ('Psi', Ps)):
    acc = 0
    for k in range(3):
        wk = sp.together(w * coef[str(J[k])])
        dk = wk
        for _ in range(k):
            dk = dS(dk, rules)
        acc += (-1)**k * dk
    red[nm] = sp.factor(sp.simplify(sp.together(acc.subs(rules))))
out = {'coefficients': {k: str(v) for k, v in coef.items()},
       'linear': bool(lin), 'gauge_residual': str(gauge_res),
       'weighted_reduced_density': {k: str(v) for k, v in red.items()}}
print(json.dumps(out, indent=2))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'l0_full_sector.json'), 'w'), indent=2)
