"""Independent re-derivation of the small-delta series of the static +1 branch (auditor's code).

Written from scratch (does not import or read the audited series code).  Differences from the audited method:
  * bulk coefficients solved by direct residual division by the known linear operators
    (p-14)(p+18) for eta and 2p for P (no sympy.solve), in exact Fractions (the bulk problem has no c);
  * the ansatz is extended to total degree 9 INCLUDING the secular term at the (2,7) resonance:
        eta  ⊃  (nu + kappa*u) X^2 T^7 ,   nu arbitrary (free),  kappa fixed by the equations;
  * junctions solved with c, nu and U := u_b (log of T_b) kept as symbols, so that the dependence of the
    delta^9 shell coefficients on the free constant nu and on log(delta) can be read off exactly.
Model (package conventions): u = y/9, eta = phi-1, R = rho/9 = (e^u/2) P(X,T), X = alpha e^{14u}, T = e^{-2u}
  R_u^2 = 1 + R^2 (eta_u^2/12 + V),  V = 1 - 21 e^2 - 25 e^3 + 9/4 e^4 + 6 e^5 + e^6
  eta_uu + 4 (R_u/R) eta_u = F,       F = 252 e + 450 e^2 - 54 e^3 - 180 e^4 - 36 e^5
  shell: R_u/R = 3W + (3/2) delta (1 + c phi),  eta_u = -(9/2)(2 W_phi + c delta),  W(1+e) = 1/3 + e^2 + e^3/3
The potential expansions V, F are re-derived here from W with sympy (not copied).
Output: verify/V1_SERIES_INDEPENDENT.json
"""
import json, time
from fractions import Fraction as Fr
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
t0 = time.time()
out = {}

# ---------- 0. re-derive V and F from the superpotential
ph, e = sp.symbols('phi e')
W = 1 - ph + ph**3 / 3
U = sp.Rational(1, 2) * sp.diff(W, ph)**2 - sp.Rational(2, 3) * W**2
k2 = -U.subs(ph, 1) / 6
Vpoly = sp.Poly(sp.expand((-U / (6 * k2)).subs(ph, 1 + e)), e)
Fpoly = sp.Poly(sp.expand((sp.diff(U, ph) / k2).subs(ph, 1 + e)), e)
Wsh = sp.expand(W.subs(ph, 1 + e)); Wps = sp.expand(sp.diff(W, ph).subs(ph, 1 + e))
out['potentials'] = dict(k2=str(k2), V=str(Vpoly.as_expr()), F=str(Fpoly.as_expr()), W_at_1_plus_e=str(Wsh), Wphi=str(Wps))
Vc = {m[0]: Fr(int(sp.numer(v)), int(sp.denom(v))) for m, v in Vpoly.terms()}
Fc = {m[0]: Fr(int(sp.numer(v)), int(sp.denom(v))) for m, v in Fpoly.terms()}
assert k2 == sp.Rational(1, 81)

K = 9
# series elements: dict {(m, j, l): coeff}, l = power of u (0 or 1)
def tr(a):
    return {k: v for k, v in a.items() if k[0] + k[1] <= K and k[2] <= 1 and v != 0}
def add(*xs):
    o = {}
    for a in xs:
        for k, v in a.items():
            o[k] = o.get(k, 0) + v
    return tr(o)
def sc(s, a):
    return tr({k: s * v for k, v in a.items()})
def mul(a, b):
    o = {}
    for (m1, j1, l1), v in a.items():
        for (m2, j2, l2), w in b.items():
            if m1 + m2 + j1 + j2 <= K and l1 + l2 <= 1:
                kk = (m1 + m2, j1 + j2, l1 + l2)
                o[kk] = o.get(kk, 0) + v * w
    return tr(o)
def D(a):   # d/du of X^m T^j u^l = p X^m T^j u^l + l X^m T^j u^{l-1}
    o = {}
    for (m, j, l), v in a.items():
        p = 14 * m - 2 * j
        o[(m, j, l)] = o.get((m, j, l), 0) + p * v
        if l:
            o[(m, j, l - 1)] = o.get((m, j, l - 1), 0) + l * v
    return tr(o)
def poly_of(coeffs, x, maxdeg):
    pw = [{(0, 0, 0): Fr(1)}]
    for _ in range(maxdeg):
        pw.append(mul(pw[-1], x))
    return add(*[sc(cf, pw[d]) for d, cf in coeffs.items()])

def residuals(N, P):
    DN = D(N); DP = D(P)
    Fs = poly_of(Fc, N, 5)
    E1 = add(mul(P, D(DN)), sc(4, mul(add(P, DP), DN)), sc(-1, mul(P, Fs)))
    Vs = poly_of(Vc, N, 6)
    Q = add(sc(Fr(1, 12), mul(DN, DN)), Vs)
    PD = add(P, DP)
    E2 = add(mul(PD, PD), {(0, 1, 0): Fr(-4)}, sc(-1, mul(mul(P, P), Q)))
    return E1, E2

def solve_bulk(nu_free, with_secular=True):
    N = {(1, 0, 0): Fr(1)}
    P = {(0, 0, 0): Fr(1), (0, 1, 0): Fr(-1)}
    kappa = None; S27 = None
    for deg in range(1, K + 1):
        E1, E2 = residuals(N, P)
        for (m, j, l), v in list(E1.items()):
            if m + j != deg or l != 0 or m < 1:
                continue
            p = 14 * m - 2 * j
            op = (p - 14) * (p + 18)
            if op != 0:
                N[(m, j, 0)] = N.get((m, j, 0), 0) - v / op
            else:
                S27 = v
                if with_secular:
                    # L[kappa u X^m T^j] = kappa (2p+4) X^m T^j at leading order
                    kappa = -v / (2 * p + 4)
                    N[(m, j, 1)] = kappa
                    N[(m, j, 0)] = Fr(nu_free)
        for (m, j, l), v in list(E2.items()):
            if m + j != deg or l != 0 or m < 2:
                continue
            p = 14 * m - 2 * j
            P[(m, j, 0)] = P.get((m, j, 0), 0) - v / (2 * p)
        E1, E2 = residuals(N, P)
        bad = {k: v for k, v in list(E1.items()) + list(E2.items()) if k[0] + k[1] <= deg and v != 0}
        if bad and not (not with_secular and deg == 9):
            raise AssertionError(('residual', deg, bad))
    return N, P, kappa, S27

Nb, Pb, kappa, S27 = solve_bulk(0)
Nb5, Pb5, _, _ = solve_bulk(5)
out['bulk'] = dict(S_2_7=str(S27), kappa=str(kappa),
                   note='E1 residual at X^2T^7 with n_{2,7}=0 is S; secular coefficient kappa = -S/32; residuals vanish through degree 9 with the secular term')
# degree-9 check without the secular term: must fail (control)
try:
    solve_bulk(0, with_secular=False)
    E1, E2 = residuals(*solve_bulk(0, with_secular=False)[:2])
    out['bulk']['control_no_secular_residual_at_2_7'] = str(E1.get((2, 7, 0), 0))
except AssertionError as ex:
    out['bulk']['control_no_secular'] = 'fails: ' + str(ex)[:120]

# linear part vs Gegenbauer C_14^(2)(cosh u)/15 : sum_j (j+1)(15-j)/15 e^{(14-2j)u}
lin_ok = all(Nb.get((1, j, 0), 0) == Fr((j + 1) * (15 - j), 15) for j in range(0, K))
out['bulk']['linear_part_equals_C14_over_15'] = lin_ok
assert lin_ok

# ---------- junctions with symbols c, nu, U (U = u_b, enters only through the secular term)
c, d, nu, Ub = sp.symbols('c delta nu U')
L = K + 1
def zs(): return [sp.Integer(0)] * L
def smul(a, b):
    o = zs()
    for i in range(L):
        if a[i] == 0: continue
        for j in range(L - i):
            if b[j] != 0: o[i + j] += a[i] * b[j]
    return [sp.expand(x) for x in o]
def sadd(*xs):
    o = zs()
    for a in xs:
        for i in range(L): o[i] += a[i]
    return [sp.expand(x) for x in o]
def ssc(k, a): return [sp.expand(k * x) for x in a]
def cst(k): o = zs(); o[0] = sp.sympify(k); return o
def evalser(S, Xp, Tp, nuval):
    o = zs()
    for (m, j, l), v in S.items():
        coef = sp.Rational(v.numerator, v.denominator)
        if l == 1:
            coef = coef * Ub
        o = sadd(o, ssc(coef, smul(Xp[m], Tp[j])))
    return o
# series with symbolic nu: take Nb (nu=0) and replace (2,7,0) by nu; D acts on u-term: D(kappa u X^2T^7) = 14 kappa u X^2T^7 + kappa X^2T^7
NbS = {k: v for k, v in Nb.items() if k != (2, 7, 0)}   # nu enters only via the explicit symbolic term below
DNbS = D(NbS)   # (2,7,0) of DNbS contains kappa (from u-term) + 0; add 14*nu separately
Xs = zs(); Ts = zs(); dl = zs(); dl[1] = sp.Integer(1)
def junctions(Xs, Ts):
    Xp = [cst(1)]; Tp = [cst(1)]
    for _ in range(K):
        Xp.append(smul(Xp[-1], Xs)); Tp.append(smul(Tp[-1], Ts))
    eta = evalser(NbS, Xp, Tp, None)
    eta = sadd(eta, ssc(nu, smul(Xp[2], Tp[7])))
    etaD = evalser(DNbS, Xp, Tp, None)
    etaD = sadd(etaD, ssc(14 * nu, smul(Xp[2], Tp[7])))
    Pv = evalser(Pb, Xp, Tp, None); DPv = evalser(D(Pb), Xp, Tp, None)
    e2 = smul(eta, eta); e3 = smul(e2, eta)
    W = sadd(cst(sp.Rational(1, 3)), e2, ssc(sp.Rational(1, 3), e3))
    Wp = sadd(ssc(2, eta), e2)
    # R_u/R = 1 + DP/P ;  J1 := P*(R_u/R) - P*(3W + 3/2 delta (1 + c(1+eta)))
    J1 = sadd(Pv, DPv, ssc(-1, smul(Pv, sadd(ssc(3, W), ssc(sp.Rational(3, 2), smul(dl, sadd(cst(1 + c), ssc(c, eta))))))))
    J2 = sadd(etaD, ssc(sp.Rational(9, 2), sadd(ssc(2, Wp), ssc(c, dl))))
    return J1, J2, eta, Pv
for n in range(1, L):
    J1, J2, _, _ = junctions(Xs, Ts)
    # linear dependence on (t_n, x_n) at order n: J1 ~ 2 t_n,  J2 ~ 32 x_n  (checked below by vanishing)
    Ts[n] = sp.factor(-J1[n] / 2); Xs[n] = sp.factor(-J2[n] / 32)
J1, J2, eta, Pv = junctions(Xs, Ts)
assert all(sp.expand(x) == 0 for x in J1) and all(sp.expand(x) == 0 for x in J2)
# H^2 = 1/rho_b^2 = k^2/R_b^2 = (4/81) T_b / P_b^2
def sinv(a):
    o = zs(); o[0] = 1 / a[0]
    for n in range(1, L):
        o[n] = sp.expand(-sum(a[k] * o[n - k] for k in range(1, n + 1)) / a[0])
    return o
Pinv = sinv(Pv)
H2 = ssc(sp.Rational(4, 81), smul(Ts, smul(Pinv, Pinv)))
eta_c = [sp.factor(x) for x in eta]; H2_c = [sp.factor(x) for x in H2]
out['eta_b'] = [str(x) for x in eta_c]; out['H2'] = [str(x) for x in H2_c]
out['T_b'] = [str(sp.factor(x)) for x in Ts]; out['X_b'] = [str(sp.factor(x)) for x in Xs]
out['delta9_dependence'] = dict(
    eta_b_depends_on_nu=bool(sp.diff(eta_c[9], nu) != 0), eta_b_depends_on_U=bool(sp.diff(eta_c[9], Ub) != 0),
    H2_depends_on_nu=bool(sp.diff(H2_c[9], nu) != 0), H2_depends_on_U=bool(sp.diff(H2_c[9], Ub) != 0),
    X_b9_depends_on_nu=bool(sp.diff(Xs[9], nu) != 0), X_b9_depends_on_U=bool(sp.diff(Xs[9], Ub) != 0))

# compare with audited exact coefficients through delta^8
aud = json.loads((HERE.parent / 'SERIES_COEFFICIENTS.json').read_text())
cmp = {}
for name, mine in [('eta_b', eta_c), ('H2', H2_c), ('X_b', [sp.factor(x) for x in Xs]), ('T_b', [sp.factor(x) for x in Ts])]:
    theirs = [sp.sympify(s, locals={'c': c}) for s in aud[name]]
    cmp[name] = [bool(sp.simplify(mine[i] - theirs[i]) == 0) for i in range(len(theirs))]
out['match_audited_through_delta8'] = cmp
creg = sp.Rational(2) / sp.Rational('1.0357712571566784') - sp.Rational(4, 3)
vals = {}
for lab, cv in [('reg', creg), ('cm04', sp.Rational(-2, 5)), ('cp13', sp.Rational(13, 10))]:
    vals[lab] = dict(eta_b_delta9=str(sp.N(eta_c[9].subs(c, cv), 25)), H2_delta9=str(sp.N(H2_c[9].subs(c, cv), 25)),
                     eta_b=[str(sp.N(x.subs(c, cv), 30)) for x in eta_c], H2=[str(sp.N(x.subs(c, cv), 30)) for x in H2_c])
out['values'] = vals
out['eta_b_delta9_exact'] = str(eta_c[9]); out['H2_delta9_exact'] = str(H2_c[9])
out['runtime_s'] = time.time() - t0
(HERE / 'V1_SERIES_INDEPENDENT.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: out[k] for k in ['bulk', 'delta9_dependence', 'match_audited_through_delta8']}, indent=1))
for lab in vals:
    print(lab, vals[lab]['eta_b_delta9'], vals[lab]['H2_delta9'])
print('runtime', out['runtime_s'])
