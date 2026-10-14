"""Symbolic (all-s) solution of the metric-only adjoint equations P_g^* N = 0
for N = n(s) (l = 0) and N = n(s) x_1 (l = 1).  Reduces each system to ODEs and
shows: l = 1 kernel = span{Omega x_1}; l = 0 kernel is at most one-dimensional,
given by an explicit first-order ODE plus a compatibility condition."""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'basis', 'code'))
import sympy as sp
from engine import X, O, F, register, reduce_bg, bg_rules, dS
from kids_tensor_only import Tstar
rules, s = bg_rules()
Nj = register('M')   # fresh jets for n(s)
n0, n1, n2, n3 = Nj[0], Nj[1], Nj[2], Nj[3]
O0, O1, F1 = O[0], O[1], F[1]
def comp(e, rep):
    return sp.factor(sp.together(reduce_bg(e).subs(rep)))
out = {}
# ---- l = 0 : T_ij = cd delta_ij + cx x_i x_j
T = Tstar(n0)
cd = comp(T[0][0], {X[0]: 0, X[1]: 0, X[2]: 0, X[3]: sp.sqrt(s)})
cx = comp(T[0][1] / (X[0] * X[1]), {X[0]: 1, X[1]: 1, X[2]: 0, X[3]: sp.sqrt(s - 2)})
sol2 = sp.solve(cx, n2)[0]
first = sp.factor(sp.numer(sp.together(cd.subs(n2, sol2))))
lam = sp.solve(first, n1)[0]            # n' = lam n
lam_c = sp.factor(lam / n0)
# compatibility: n'' from differentiating n' = lam n must equal sol2
n2_from = sp.together((dS(lam_c, rules) + lam_c**2) * n0)
compat = sp.factor(sp.together((n2_from - sol2.subs(n1, lam)).subs(rules)))
out['l0'] = {'cd': str(cd), 'cx': str(cx), 'lambda': str(lam_c),
             'compatibility_residual': str(compat)}
print('l0 lambda =', lam_c); print('l0 compatibility residual =', compat, flush=True)
# ---- l = 1 : N = n(s) x1
T = Tstar(n0 * X[0])
x1, x2, x3, x4 = X
# T_ij = A x1 delta_ij + B x1 x_i x_j + C (delta_i1 x_j + delta_j1 x_i)
A_ = comp(T[1][1] / X[0], {X[1]: 0, X[2]: 0})
A_ = sp.factor(sp.together(reduce_bg(T[1][1] / X[0]).subs({X[1]: 0, X[2]: 0, X[0]: sp.Rational(1, 2), X[3]: sp.sqrt(s - sp.Rational(1, 4))})))
B_ = sp.factor(sp.together(reduce_bg(T[1][2] / (X[0] * X[1] * X[2])).subs({X[0]: 1, X[1]: 1, X[2]: 1, X[3]: sp.sqrt(s - 3)})))
C_ = sp.factor(sp.together(reduce_bg((T[0][1] - B_.subs(s, sum(v**2 for v in X)) * X[0]**2 * X[1]) / X[1]).subs({X[0]: 1, X[1]: 1, X[2]: 0, X[3]: sp.sqrt(s - 2)})))
eqs = [sp.numer(sp.together(e)) for e in (A_, B_, C_)]
Mx = sp.Matrix([[sp.diff(e, v) for v in (n0, n1, n2)] for e in eqs])
Mx = Mx.applyfunc(lambda z: sp.factor(sp.together(z.subs(rules))))
ns = Mx.nullspace(simplify=True)
target = sp.Matrix([O0, O1, rules[O[2]]])
isO = [bool(sp.simplify(v[0] * target[1] - v[1] * target[0]) == 0 and
            sp.simplify(v[0] * target[2] - v[2] * target[0]) == 0) for v in ns]
minors = [sp.factor(Mx.extract([i, j], [0, 1]).det()) for (i, j) in ((0, 1), (0, 2), (1, 2))]
out['l1'] = {'A': str(A_), 'B': str(B_), 'C': str(C_), 'rank': int(Mx.rank(simplify=True)),
             'nullspace_dim': len(ns), 'nullspace_is_Omega_jet': isO,
             'two_by_two_minors_in_(n,n1)': [str(m) for m in minors]}
print('l1', out['l1'], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kids_symbolic_l01.json'), 'w'), indent=2)
