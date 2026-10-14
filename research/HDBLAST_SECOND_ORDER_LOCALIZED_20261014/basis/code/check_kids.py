"""Local KIDs of the time-symmetric slice on an annulus.

P_H(h,psi) = DH[h,psi] (non-densitized, at pi = varpi = 0):
  DH = 1/2(-Lap tr h + div div h - <Ric,h>) + 1/2 h^{ij}phi_i phi_j
       - <D phi, D psi> - U' psi
Formal L^2(dmu_gamma) adjoint:
  P_H^* N = ( 1/2(Hess N - gamma Lap N - N Ric) + 1/2 N dphi dphi ,
              div(N D phi) - N U' ).
P_M(pi',varpi') = -2 div pi' + varpi' dphi;  P_M^* Y ~ (L_Y gamma, Y(phi)).

Checks: (i) the formal adjoint identity int N P_H(h,psi) = int <P_H^*N,(h,psi)>
up to a divergence, verified pointwise via an explicit divergence-free test
(skipped: standard) -- instead we verify directly that
(ii) boost lapses N = Omega x_i lie in ker P_H^*,
(iii) rotations Y = M x lie in ker P_M^*,
(iv) radial N(s): the three radial equations force N = 0 unless an explicit
background combination vanishes.
"""
import json, time
import sympy as sp
from engine import X, O, F, U, J3, register, d, reduce_bg, is_zero, S2, christoffel, ricci
from fields import gamma0, phi0_grad

out = {}
g0 = gamma0()
gi0 = g0.inv()
g = [[J3(g0[i, j]) for j in range(4)] for i in range(4)]
gi = [[J3(gi0[i, j]) for j in range(4)] for i in range(4)]
Gam = christoffel(g, gi)
Gam0 = [[[Gam[k][i][j].a for j in range(4)] for i in range(4)] for k in range(4)]
Ric0 = [[c.a for c in row] for row in ricci(Gam)]
gp = phi0_grad()
O0 = O[0]
sqrtg = O0**4

def hess(N):
    dN = [d(N, k) for k in range(4)]
    return [[d(dN[j], i) - sum(Gam0[k][i][j] * dN[k] for k in range(4))
             for j in range(4)] for i in range(4)]

def lap(N):
    Hs = hess(N)
    return sum(gi0[i, j] * Hs[i][j] for i in range(4) for j in range(4))

def PHstar(N):
    Hs = hess(N)
    L = lap(N)
    T = [[sp.Rational(1, 2) * (Hs[i][j] - g0[i, j] * L - N * Ric0[i][j])
          + sp.Rational(1, 2) * N * gp[i] * gp[j] for j in range(4)] for i in range(4)]
    # div(N D phi) with density: (1/sqrtg) d_i (sqrtg N gamma^{ij} phi_j)
    vec = [sqrtg * N * sum(gi0[i, j] * gp[j] for j in range(4)) for i in range(4)]
    sc = sum(d(vec[i], i) for i in range(4)) / sqrtg - N * U[1]
    return T, sc

# (ii) boost lapses
okb = True
for a in range(4):
    T, sc = PHstar(O0 * X[a])
    for i in range(4):
        for j in range(i, 4):
            okb &= is_zero(T[i][j])[0]
    okb &= is_zero(sc)[0]
out['boost_lapses_in_kernel'] = okb
print('boost lapse N = Omega x_a in ker P_H^*:', okb, flush=True)

# (iii) rotations: L_Y gamma0 and Y(phi0)
okr = True
for (a, b) in [(0, 1), (0, 2), (1, 3), (2, 3)]:
    Y = [0] * 4
    Y[a], Y[b] = -X[b], X[a]
    Yl = [g0[i, i] * Y[i] for i in range(4)]
    for i in range(4):
        for j in range(4):
            LY = d(Yl[j], i) + d(Yl[i], j) - 2 * sum(Gam0[k][i][j] * Yl[k] for k in range(4))
            okr &= is_zero(LY)[0]
    okr &= is_zero(sum(Y[i] * gp[i] for i in range(4)))[0]
out['rotations_in_kernel_PM_star'] = okr
print('rotations in ker P_M^*:', okr, flush=True)

# (iv) radial N(s)
Nj = register('N')
T, sc = PHstar(Nj[0])
# decompose T_ij = cd * delta_ij + cx * x_i x_j  (radial symmetry)
cd = sp.together(reduce_bg(T[0][0]).subs({X[0]: 0}))  # at x1 = 0: T_11 = cd
cx = sp.together(reduce_bg(T[0][1] / (X[0] * X[1])))
scr = sp.together(reduce_bg(sc))
s = sp.Symbol('s_')
# express in s: replace x-polynomial dependence (only through |x|^2) by s
def to_s(e):
    e = sp.together(e)
    num, den = sp.fraction(e)
    # substitute x4^2 = s - x1^2 - x2^2 - x3^2 by evaluating on x=(0,0,0,sqrt s)
    rep = {X[0]: 0, X[1]: 0, X[2]: 0, X[3]: sp.sqrt(s)}
    return sp.simplify(sp.together(num.subs(rep)) / sp.together(den.subs(rep)))
cd_s, cx_s, sc_s = to_s(cd), to_s(cx), to_s(scr)
# check consistency of the x-structure: T_ij - cd delta - cx x x == 0
okstruct = True
for i in range(4):
    for j in range(4):
        e = T[i][j] - (cd_s.subs(s, S2) if i == j else 0) - cx_s.subs(s, S2) * X[i] * X[j]
        okstruct &= is_zero(e)[0]
out['radial_tensor_structure'] = okstruct
lam = sp.Symbol('lam'); lamp = sp.Symbol('lamp')
# scalar equation is first order: solve for N1 = N'(s)
N0, N1, N2 = Nj[0], Nj[1], Nj[2]
solN1 = sp.solve(sp.Eq(sc_s, 0), N1)
print('scalar eq solved for N_s:', solN1, flush=True)
res = {}
if len(solN1) == 1:
    n1 = solN1[0]
    # second derivative from differentiating the first-order relation
    from engine import dS, bg_rules
    rules, _ = bg_rules()
    n2 = dS(n1, rules)
    n2 = n2.subs(N1, n1)
    for name, e in [('delta', cd_s), ('xx', cx_s)]:
        e2 = sp.factor(sp.together(e.subs(N2, n2).subs(N1, n1)))
        res[name] = str(e2)
        print(name, 'component after scalar eq:', e2, flush=True)
out['radial_reduced'] = res
out['radial_raw'] = {'cd': str(sp.factor(cd_s)), 'cx': str(sp.factor(cx_s)),
                     'scalar': str(sp.factor(sc_s))}
json.dump(out, open('../results/kids.json', 'w'), indent=2)
print(json.dumps(out, indent=2))
