"""Exact small-detuning expansion of the static "+1" branch (exact in c).

Method (see README.md, section 1):
  u = y/9 = k y  (k = 1/9 is the AdS curvature at phi = 1),  R = k*rho,  eta = phi - 1.
  Bulk equations (regular cone at u = 0):
     R_u^2 = 1 + R^2 (eta_u^2/12 + V(eta)),   V = -U/(6k^2) = 1 - 21 eta^2 - 25 eta^3 + 9/4 eta^4 + 6 eta^5 + eta^6
     eta_uu + 4 (R_u/R) eta_u = 252 eta + 450 eta^2 - 54 eta^3 - 180 eta^4 - 36 eta^5   (= U_eta/k^2)
  Junctions (single shell, doubled bulk):
     R_u/R = (3/2) sigma = 3 W + (3/2) delta (1 + c phi),      eta_u = -(9/2) sigma'
  Formal double series in X = alpha e^{14u} (regular growing mode) and T = e^{-2u}:
     eta = N(X,T) = X (1 + ...) + O(X^2),   R = (e^u/2) P(X,T),  P = 1 - T + O(X^2).
  D = d/du acts as D X^m T^j = (14 m - 2 j) X^m T^j.  Cone regularity only enters through
  homogeneous pieces suppressed by alpha ~ delta^8 (see README), and the first resonance
  (p = 14 m - 2 j = 14 with m = 2, j = 7) occurs at total degree 9.  The formal series is
  therefore trusted through total degree 8 in (X_b, T_b) ~ (delta, delta).
Outputs exact rational-in-c coefficients to JSON.
"""
import json, sys, time, hashlib
from pathlib import Path
import sympy as sp

K = int(sys.argv[1]) if len(sys.argv) > 1 else 8   # total degree kept (<= 8 before resonance)
assert K <= 8, "degree 9 hits the X^2 T^7 resonance (log terms); not implemented"
OUT = Path(__file__).resolve().parent

c, d = sp.symbols('c delta')

# ---- truncated bivariate polynomial arithmetic: dict {(m,j): coeff}
def trunc(p):
    return {k: v for k, v in p.items() if k[0] + k[1] <= K and v != 0}
def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return trunc(out)
def scal(s, p):
    return trunc({k: s * v for k, v in p.items()})
def mul(p, q):
    out = {}
    for (a, b), v in p.items():
        for (e, f), w in q.items():
            if a + b + e + f <= K:
                out[(a + e, b + f)] = out.get((a + e, b + f), 0) + v * w
    return trunc(out)
def D(p):
    return trunc({(m, j): (14 * m - 2 * j) * v for (m, j), v in p.items()})
def powers(p, n):
    res = [{(0, 0): sp.Integer(1)}]
    for _ in range(n):
        res.append(mul(res[-1], p))
    return res

# unknown coefficients
nsym = {}; psym = {}
for m in range(1, K + 1):
    for j in range(0, K + 1 - m):
        if (m, j) != (1, 0):
            nsym[(m, j)] = sp.Symbol(f'n_{m}_{j}')
for m in range(2, K + 1):
    for j in range(0, K + 1 - m):
        psym[(m, j)] = sp.Symbol(f'p_{m}_{j}')
N = {(1, 0): sp.Integer(1)}; N.update(nsym)
P = {(0, 0): sp.Integer(1), (0, 1): sp.Integer(-1)}; P.update(psym)

def residuals(N, P):
    Np = powers(N, 6); DN = D(N); DP = D(P)
    F = add(scal(252, Np[1]), scal(450, Np[2]), scal(-54, Np[3]), scal(-180, Np[4]), scal(-36, Np[5]))
    E1 = add(mul(P, D(DN)), scal(4, mul(add(P, DP), DN)), scal(-1, mul(P, F)))
    Q = add(scal(sp.Rational(1, 12), mul(DN, DN)), scal(-21, Np[2]), scal(-25, Np[3]),
            scal(sp.Rational(9, 4), Np[4]), scal(6, Np[5]), Np[6])
    PD = add(P, DP)
    E2 = add(mul(PD, PD), {(0, 1): sp.Integer(-4)}, scal(-1, mul(mul(P, P), add({(0, 0): sp.Integer(1)}, Q))))
    return E1, E2

t0 = time.time()
sol = {}
for deg in range(1, K + 1):
    Ns = {k: (v.subs(sol) if hasattr(v, 'subs') else v) for k, v in N.items()}
    Ps = {k: (v.subs(sol) if hasattr(v, 'subs') else v) for k, v in P.items()}
    E1, E2 = residuals(Ns, Ps)
    eqs = []; unk = []
    for (m, j), v in E1.items():
        if m + j == deg and m >= 1:
            eqs.append(sp.expand(v))
    for (m, j), v in E2.items():
        if m + j == deg:
            eqs.append(sp.expand(v))
    unk = [s for k, s in nsym.items() if k[0] + k[1] == deg] + [s for k, s in psym.items() if k[0] + k[1] == deg]
    eqs = [e for e in eqs if e != 0]
    if unk:
        s = sp.solve(eqs, unk, dict=True)
        assert len(s) == 1, (deg, eqs, unk)
        sol.update({k: sp.Rational(v) for k, v in s[0].items()})
    # all remaining equations of this degree must vanish identically
    E1, E2 = residuals({k: (v.subs(sol) if hasattr(v, 'subs') else v) for k, v in N.items()},
                       {k: (v.subs(sol) if hasattr(v, 'subs') else v) for k, v in P.items()})
    for (m, j), v in list(E1.items()) + list(E2.items()):
        if m + j <= deg:
            assert sp.simplify(v) == 0, ('residual', deg, (m, j), v)
Nf = {k: sp.Rational(v.subs(sol)) if hasattr(v, 'subs') else v for k, v in N.items()}
Pf = {k: sp.Rational(v.subs(sol)) if hasattr(v, 'subs') else v for k, v in P.items()}
print('bulk series solved', time.time() - t0, flush=True)

# check: linear part equals C_14^(2)(cosh u)/15 * alpha
u = sp.symbols('u', positive=True)
lin = sum(v * sp.exp((14 - 2 * j) * u) for (m, j), v in Nf.items() if m == 1)
w = sp.Symbol('w')
geg_poly = sp.Poly(sp.expand(sp.gegenbauer(14, 2, (w + 1 / w) / 2) / 15 * w ** 14), w)
gdict = {e: v for (e,), v in geg_poly.terms()}
lin_check = []
for (m, j), v in Nf.items():
    if m == 1:
        # term e^{(14-2j)u}; in geg*e^{14u} it is w^{28-2j}
        lin_check.append(sp.simplify(v - gdict.get(28 - 2 * j, 0)))
assert all(x == 0 for x in lin_check), lin_check

# ---- junctions: solve X_b, T_b as series in delta (truncated series arithmetic, coefficients polynomial in c)
L = K + 1
def sz(): return [sp.Integer(0)] * L
def smul(a, b):
    out = sz()
    for i in range(L):
        if a[i] == 0: continue
        for j in range(L - i):
            if b[j] != 0: out[i + j] += a[i] * b[j]
    return [sp.expand(x) for x in out]
def sadd(*xs):
    out = sz()
    for a in xs:
        for i in range(L): out[i] += a[i]
    return [sp.expand(x) for x in out]
def sscal(k, a): return [sp.expand(k * x) for x in a]
def const(k): o = sz(); o[0] = sp.sympify(k); return o
def sev(p, Xp, Tp):
    out = sz()
    for (m, j), v in p.items():
        term = smul(Xp[m], Tp[j])
        out = sadd(out, sscal(v, term))
    return out
def spows(a, n):
    r = [const(1)]
    for _ in range(n): r.append(smul(r[-1], a))
    return r
DNf = D(Nf); DPf = D(Pf)
Xs = sz(); Ts = sz(); deltaS = sz(); deltaS[1] = sp.Integer(1)
for order in range(1, K + 1):
    Xp = spows(Xs, K); Tp = spows(Ts, K)
    eta = sev(Nf, Xp, Tp); etaD = sev(DNf, Xp, Tp); Pb = sev(Pf, Xp, Tp); DPb = sev(DPf, Xp, Tp)
    e2 = smul(eta, eta); e3 = smul(e2, eta)
    br = sadd(sscal(3, e2), e3, sscal(sp.Rational(3, 2), smul(deltaS, sadd(const(1 + c), sscal(c, eta)))))
    J1 = sadd(DPb, sscal(-1, smul(Pb, br)))
    J2 = sadd(etaD, sscal(18, eta), sscal(9, e2), sscal(sp.Rational(9, 2) * c, deltaS))
    Ts[order] = sp.factor(-J1[order] / 2); Xs[order] = sp.factor(-J2[order] / 32)
    print('junction order', order, time.time() - t0, flush=True)
# verify junctions vanish through order K
Xp = spows(Xs, K); Tp = spows(Ts, K)
eta = sev(Nf, Xp, Tp); etaD = sev(DNf, Xp, Tp); Pb = sev(Pf, Xp, Tp); DPb = sev(DPf, Xp, Tp)
e2 = smul(eta, eta); e3 = smul(e2, eta)
br = sadd(sscal(3, e2), e3, sscal(sp.Rational(3, 2), smul(deltaS, sadd(const(1 + c), sscal(c, eta)))))
J1 = sadd(DPb, sscal(-1, smul(Pb, br))); J2 = sadd(etaD, sscal(18, eta), sscal(9, e2), sscal(sp.Rational(9, 2) * c, deltaS))
assert all(sp.expand(x) == 0 for x in J1) and all(sp.expand(x) == 0 for x in J2)
def sinv(a):
    # a[0] != 0
    out = sz(); out[0] = 1 / a[0]
    for n in range(1, L):
        out[n] = sp.expand(-sum(a[k] * out[n - k] for k in range(1, n + 1)) / a[0])
    return out
def spow_half_neg(a):
    # a[0] = 1 ; returns a^(-1/2) via recurrence for f = a^r: n a0 f_n = sum_{k=1}^n ((r+1)k - n) a_k f_{n-k}
    r = sp.Rational(-1, 2); f = sz(); f[0] = sp.Integer(1)
    for n in range(1, L):
        f[n] = sp.expand(sum(((r + 1) * k - n) * a[k] * f[n - k] for k in range(1, n + 1)) / n)
    return f
def slog1(a):
    # log(a), a[0]=1 : (log a)' = a'/a
    ia = sinv(a); da = [sp.expand((i + 1) * a[i + 1]) for i in range(L - 1)] + [0]
    q = smul(da, ia); out = sz()
    for n in range(1, L): out[n] = sp.expand(q[n - 1] / n)
    return out
eta_c = [sp.factor(x) for x in eta]
Pinv = sinv(Pb); H2s = sscal(sp.Rational(4, 81), smul(Ts, smul(Pinv, Pinv)))
H2_c = [sp.factor(x) for x in H2s]
# ratio = H2 * 27/(delta (1+c)): shift by one power of delta
ratio = [sp.expand(H2s[i + 1] * 27 / (1 + c)) for i in range(L - 1)] + [sp.Integer(0)]
Sser = spow_half_neg(ratio); S_c = [sp.factor(x) for x in Sser[:K]]
# eta_h/delta^8 = (136/3) (X_b/delta) (T_b/delta)^7
Xd = Xs[1:] + [sp.Integer(0)]; Td = Ts[1:] + [sp.Integer(0)]
tmp = Xd
for _ in range(7): tmp = smul(tmp, Td)
etah_c = [sp.factor(sp.Rational(136, 3) * x) for x in tmp[:K]]
t1 = Ts[1]
lnf = slog1([sp.expand(x / t1) for x in Td])
yb_c = [sp.factor(-sp.Rational(9, 2) * x) for x in lnf[:K]]
known_ok = dict(
    eta_first=sp.simplify(eta_c[1] + sp.Rational(9, 64) * c) == 0,
    H2_first=sp.simplify(H2_c[1] - (1 + c) / 27) == 0,
    H2_second=sp.simplify(H2_c[2] - ((1 + c) ** 2 / 36 - c ** 2 / 384)) == 0,
)
print(known_ok)
assert all(known_ok.values())

def s(x):
    return str(sp.factor(x))
res = dict(
    status='EXACT_FORMAL_SERIES (exact rational arithmetic in c); validity through total degree %d' % K,
    conventions='u=y/9, eta=phi-1, T_b=e^{-2u_b}, X_b=alpha e^{14 u_b}; eta_b = sum a_n delta^n; H^2 = sum h_n delta^n;'
                ' rho_b = sqrt(27/(delta(1+c))) * sum s_n delta^n; y_b = (9/2) ln(4/(3(1+c)delta)) + sum y_n delta^n;'
                ' eta_h = delta^8 * sum e_n delta^n (leading-order cone relation eta_h=(136/3) X_b T_b^7, valid up to relative O(delta^8)).',
    K=K,
    known_terms_verified=known_ok,
    eta_b=[s(x) for x in eta_c], H2=[s(x) for x in H2_c], rho_b_factor=[s(x) for x in S_c],
    y_b_correction=[s(x) for x in yb_c], eta_h_over_delta8=[s(x) for x in etah_c],
    X_b=[s(x) for x in Xs],
    T_b=[s(x) for x in Ts],
    bulk_N={f'{m},{j}': str(v) for (m, j), v in sorted(Nf.items())},
    bulk_P={f'{m},{j}': str(v) for (m, j), v in sorted(Pf.items())},
    linear_part_equals_gegenbauer_C14_2_over_15=True,
    runtime_seconds=time.time() - t0,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
)
(OUT / 'SERIES_COEFFICIENTS.json').write_text(json.dumps(res, indent=2) + '\n')
for name in ['eta_b', 'H2', 'rho_b_factor', 'eta_h_over_delta8']:
    print(name)
    for i, x in enumerate(res[name][:5]):
        print('  ', i, x)
