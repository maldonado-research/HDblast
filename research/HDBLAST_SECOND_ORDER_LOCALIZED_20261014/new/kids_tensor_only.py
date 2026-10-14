"""Kernel of the metric-only adjoint (psi frozen to zero) on an annulus.

P_g(h) = DH[h, psi=0] = (1/2) DR(h) + (1/2) h(dphi, dphi)  (non-densitized),
P_g^* N = (1/2)(Hess N - gamma Lap N - N Ric) + (1/2) N dphi (x) dphi.
This is a single-order (order 2) operator with injective symbol, so Delay's
Theorem 1.3 / Lemma 8.1 (arXiv:1003.0535) applies to P_g directly.
For N = n(s) Y_l(x) we compute the rank of the pointwise linear system in
(n, n', n'') with symbolic background jets.  Run from basis/code (uses the
original engine).
"""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'basis', 'code'))
import sympy as sp
from engine import X, O, F, U, J3, register, d, reduce_bg, bg_rules, christoffel, ricci
from fields import gamma0, phi0_grad

g0 = gamma0(); gi0 = g0.inv()
g = [[J3(g0[i, j]) for j in range(4)] for i in range(4)]
gi = [[J3(gi0[i, j]) for j in range(4)] for i in range(4)]
Gam = christoffel(g, gi)
Gam0 = [[[Gam[k][i][j].a for j in range(4)] for i in range(4)] for k in range(4)]
Ric0 = [[c.a for c in row] for row in ricci(Gam)]
gp = phi0_grad()

def Tstar(N):
    dN = [d(N, k) for k in range(4)]
    Hs = [[d(dN[j], i) - sum(Gam0[k][i][j] * dN[k] for k in range(4)) for j in range(4)] for i in range(4)]
    L = sum(gi0[i, j] * Hs[i][j] for i in range(4) for j in range(4))
    return [[sp.Rational(1, 2) * (Hs[i][j] - g0[i, j] * L - N * Ric0[i][j])
             + sp.Rational(1, 2) * N * gp[i] * gp[j] for j in range(4)] for i in range(4)]

Nj = register('N')
n0, n1, n2 = Nj[0], Nj[1], Nj[2]
cases = {'l0': sp.Integer(1), 'l1': X[0], 'l2': X[0]**2 - X[1]**2, 'l2b': X[0] * X[1],
         'l3': X[0]**3 - 3 * X[0] * X[1]**2, 'l4': X[0]**4 - 6 * X[0]**2 * X[1]**2 + X[1]**4}
pts = [{X[0]: sp.Rational(3, 7), X[1]: sp.Rational(-2, 5), X[2]: sp.Rational(1, 3), X[3]: sp.Rational(5, 11)},
       {X[0]: sp.Rational(-1, 3), X[1]: sp.Rational(2, 9), X[2]: sp.Rational(4, 7), X[3]: sp.Rational(1, 6)}]
rules, s = bg_rules()
res = {}
for name, Y in cases.items():
    T = Tstar(n0 * Y)
    eqs = [reduce_bg(T[i][j]) for i in range(4) for j in range(i, 4)]
    entry = {}
    for ip, pt in enumerate(pts):
        rows = []
        for e in eqs:
            e = sp.together(e.subs(pt))
            rows.append([sp.factor(sp.together(sp.diff(e, v))) for v in (n0, n1, n2)])
        Mx = sp.Matrix(rows)
        rk = Mx.rank(simplify=True)
        ns = Mx.nullspace(simplify=True)
        sval = sum(c**2 for c in pt.values())
        O2 = rules[O[2]].subs(s, sval)
        info = {'rank': int(rk), 'nullspace': [str([sp.factor(z) for z in v]) for v in ns]}
        if ns:
            v = ns[0]
            info['nullspace_is_Omega_jet'] = bool(sp.simplify(v[0] * O[1] - v[1] * O[0]) == 0 and
                                                  sp.simplify(v[0] * O2 - v[2] * O[0]) == 0)
        entry[f'point{ip}'] = info
    res[name] = entry
    print(name, entry, flush=True)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kids_tensor_only.json')
json.dump(res, open(out, 'w'), indent=2)
