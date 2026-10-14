"""Canonical form of the displacement (GG) part of the areal-radius-weighted
monopole density W = s^{3/2} Omega^5 f0 for A = diag(1,-1,0,0):
W_GG = sum_{j<=k} w_jk G_j G_k   ->  (mod total s-derivatives)
       c2 G''^2 + c1 G'^2 + c0 G^2.
Also confirms the velocity part W_KK = s^{7/2}[8 F1^2 s^2 K^2 + 3(2 s K' + 7 K)^2]/(8 F1^2 Omega)."""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'basis', 'code'))
import sympy as sp
from engine import O, F, register, bg_rules, dS
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(here, '..', 'basis', 'results', 'l0_weight.json')))
G = register('G'); K = register('K')
s = sp.Symbol('s_')
rules, _ = bg_rules()
WGG = sp.sympify(d['GG']); WKK = sp.sympify(d['KK'])
F1, O0 = F[1], O[0]
kk_hand = s**sp.Rational(7, 2) * (8 * F1**2 * s**2 * K[0]**2 + 3 * (2 * s * K[1] + 7 * K[0])**2) / (8 * F1**2 * O0)
kk_ok = sp.simplify(WKK - kk_hand) == 0
# quadratic-form coefficients
g = G[:3]
w = {}
for i in range(3):
    for j in range(i, 3):
        cij = sp.diff(WGG, g[i], g[j])
        w[(i, j)] = sp.together(cij if i != j else cij / 2)
# reduce: G0 G1 -> -(1/2) w01' G0^2 ; G1 G2 -> -(1/2) w12' G1^2 ;
#         G0 G2 -> -w02' G0 G1 - w02 G1^2 -> (1/2) w02'' G0^2 - w02 G1^2
D = lambda e, k=1: e if k == 0 else D(dS(e, rules), k - 1)
c2 = w[(2, 2)]
c1 = w[(1, 1)] - D(w[(1, 2)]) / 2 - w[(0, 2)]
c0 = w[(0, 0)] - D(w[(0, 1)]) / 2 + D(w[(0, 2)], 2) / 2
c2, c1, c0 = [sp.factor(sp.simplify(sp.together(c.subs(rules)))) for c in (c2, c1, c0)]
out = {'KK_matches_hand_derivation': bool(kk_ok), 'c2': str(c2), 'c1': str(c1), 'c0': str(c0)}
print(json.dumps(out, indent=2))
json.dump(out, open(os.path.join(here, 'wgg_canonical.json'), 'w'), indent=2)
