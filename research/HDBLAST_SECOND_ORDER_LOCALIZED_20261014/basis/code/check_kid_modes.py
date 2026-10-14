"""Mode-by-mode classification of ker P_H^* on an annulus.

For N = n(s) Y_l(x) with Y_l a harmonic polynomial of degree l, all components
of P_H^* N are linear in the jets (n0, n1, n2).  We compute the rank of that
pointwise linear system (rational in background jets) and its null space.
  l = 0: rank 3  -> n = 0 (no spherically symmetric static potential)
  l = 1: rank 2, null space (Omega, Omega', Omega'') -> N = c Omega x_a
  l >= 2: rank 3 -> n = 0 (angular trace-free Hessian)
Since P_H^* is SO(4)-equivariant, a kernel element decomposes into such
components; vanishing of each component's jets on an open set and the
second-order ODE then gives vanishing on the whole connected annulus.
"""
import json
import sympy as sp
from engine import X, O, F, register, reduce_bg, S2
from check_kids import PHstar

Nj = register('N')
n0, n1, n2 = Nj[0], Nj[1], Nj[2]
cases = {'l0': sp.Integer(1), 'l1': X[0], 'l2': X[0]**2 - X[1]**2,
         'l2b': X[0] * X[1], 'l3': X[0]**3 - 3 * X[0] * X[1]**2,
         'l4': X[0]**4 - 6 * X[0]**2 * X[1]**2 + X[1]**4}
# generic rational evaluation point (s = |x|^2 is substituted by reduce_bg)
pt = {X[0]: sp.Rational(3, 7), X[1]: sp.Rational(-2, 5), X[2]: sp.Rational(1, 3),
      X[3]: sp.Rational(5, 11)}
res = {}
for name, Y in cases.items():
    T, sc = PHstar(n0 * Y)
    eqs = [reduce_bg(T[i][j]) for i in range(4) for j in range(i, 4)] + [reduce_bg(sc)]
    rows = []
    for e in eqs:
        e = sp.together(e.subs(pt))
        rows.append([sp.diff(e, v) for v in (n0, n1, n2)])
    Mx = sp.Matrix(rows)
    # symbolic rank: use fraction-free elimination on cleared denominators
    Mx = Mx.applyfunc(lambda z: sp.factor(sp.together(z)))
    rk = Mx.rank(simplify=True)
    ns = Mx.nullspace(simplify=True)
    entry = {'rank': int(rk), 'nullspace': [str([sp.factor(z) for z in v]) for v in ns]}
    # for l = 1 compare the null vector with (Omega, Omega', Omega'') at the point
    if ns:
        from engine import bg_rules
        rules, s = bg_rules()
        O2 = rules[O[2]].subs(s, sum(c**2 for c in pt.values()))
        v = ns[0]
        target = sp.Matrix([O[0], O[1], O2])
        entry['nullspace_is_Omega_jet'] = bool(sp.simplify(v[0] * target[1] - v[1] * target[0]) == 0
                                               and sp.simplify(v[0] * target[2] - v[2] * target[0]) == 0)
    res[name] = entry
    print(name, entry, flush=True)
json.dump(res, open('../results/kid_modes.json', 'w'), indent=2)
