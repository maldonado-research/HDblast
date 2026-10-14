"""Independent numerical check of the second-order localized construction.

Shares no code with the original engine (engine.py / fields.py / linear_ops.py).
Method: explicit sympy expressions in x for every field on a concrete exact
test background, direct x-differentiation, then exact truncated power series
in eps evaluated in mpmath (40 digits) at sample points.  The constraint
maps are coded from the standard ADM formulas (4 spatial dimensions,
S = int sqrt(-g)[R/2 - (grad phi)^2/2 - U]):
  H   = R/2 + (2/det g)(tr(pi)^2/3 - |pi|^2) - varpi^2/(2 det g)
        - |D phi|^2/2 - U(phi)
  M_i = -2 D_j pi^j_i + varpi d_i phi        (pi, varpi densitized)
Test background: u = 1/Omega = 1 + S^2/5 + S^3/7 (any u with u'' > 0 works),
Phi'^2 = 3 u''/u, which solves the 5D equations (see bg5d.py).
Usage: python -I indep_check.py  -> writes indep_check.json
"""
import json, sys, time, pickle, itertools, hashlib, os
import sympy as sp
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bg5d import derive, Om as OmF, Ph as PhF, S

mp.mp.dps = 40
X = sp.symbols('x1:5', real=True)
s_x = sum(v**2 for v in X)
t0 = time.time()
out = {}

# ---------------------------------------------------------------- background
u = 1 + S**2 / 5 + S**3 / 7
Om_S = 1 / u
Php_S = sp.sqrt(3 * sp.diff(u, S, 2) / u)          # Phi'(S)
alpha, beta, box = derive()
def concrete(e):
    """Replace Om(S), Ph(S) and their derivatives by the test background."""
    reps = {}
    for k in range(4, 0, -1):
        reps[sp.Derivative(OmF(S), (S, k))] = sp.diff(Om_S, S, k)
        reps[sp.Derivative(PhF(S), (S, k))] = sp.diff(Php_S, S, k - 1)
    e = e.subs(reps)
    return e.subs(OmF(S), Om_S)
beta_c = concrete(beta)
out['background_beta_ODE_residual_samples'] = [
    float(abs(sp.N(beta_c.subs(S, sv), 30))) for sv in (sp.Rational(1, 2), 1, sp.Rational(3, 2))]
U_S = sp.Rational(3, 2) * concrete(alpha) / Om_S**2
Uphi_S = concrete(box)
Uphiphi_S = sp.diff(Uphi_S, S) / Php_S
consist = sp.diff(U_S, S) - Uphi_S * Php_S          # dU/dS = U_phi Phi'
out['background_dU_consistency_samples'] = [
    float(abs(sp.N(consist.subs(S, sv), 30))) for sv in (sp.Rational(1, 2), 1, sp.Rational(3, 2))]
print('background ok', out['background_beta_ODE_residual_samples'],
      out['background_dU_consistency_samples'], flush=True)

def onx(e):
    return e.subs(S, s_x)

Om_x = onx(Om_S); Php_x = onx(Php_S)

# ------------------------------------------------------------- profiles/data
Gf = (1 + S / 2) * sp.exp(-S)          # displacement profile (generic)
Kf = S**2 * sp.exp(-S / 3)             # velocity profile (generic)
H2 = X[0]**2 - X[1]**2
xv = sp.Matrix(X)
I4 = sp.eye(4)
P = s_x * I4 - xv * xv.T
def euler(e):
    return sum(X[i] * sp.diff(e, X[i]) for i in range(4))

def first_order(Gs, Ks):
    q = onx(Gs) * H2
    p = onx(Ks) * H2
    h = 2 * q * (s_x * I4 - 3 * xv * xv.T)
    dphi = -sp.Rational(3, 2) * (euler(q) + 4 * q) / (Om_x**2 * Php_x)
    dpi = -sp.Rational(3, 2) * p / Om_x * (xv * xv.T)
    dvarpi = -sp.Rational(3, 2) * Om_x * (euler(p) + 5 * p) / Php_x
    return h, dphi, dpi, dvarpi

# --------------------------------------- numeric radial machinery (1D only)
NJ = 9
def jets_fn(expr):
    """Return f(s) -> [expr, expr', ..., expr^(NJ-1)] (1D symbolic derivatives only)."""
    ds = [expr]
    for _ in range(NJ - 1):
        ds.append(sp.diff(ds[-1], S))
    fs = [sp.lambdify(S, e_, 'mpmath') for e_ in ds]
    return lambda sv: [fn(sv) for fn in fs]
OmJ = jets_fn(Om_S); PhpJ = jets_fn(Php_S)
U0f = sp.lambdify(S, U_S, 'mpmath'); U1f = sp.lambdify(S, Uphi_S, 'mpmath'); U2f = sp.lambdify(S, Uphiphi_S, 'mpmath')

here = os.path.dirname(os.path.abspath(__file__))
base = os.environ.get('SO_RESULTS', os.path.join(here, '..', 'basis', 'results'))
cm = json.load(open(os.path.join(base, 'construct_cc_M.json')))['nonzero_coefficients']
ch = json.load(open(os.path.join(base, 'construct_cc_H.json')))['nonzero_coefficients']
out['construction_inputs_sha256'] = {
    f_: hashlib.sha256(open(os.path.join(base, f_), 'rb').read()).hexdigest()
    for f_ in ('construct_cc_M.json', 'construct_cc_H.json')}
JET_SYMS = ([sp.Symbol(f'G{k}') for k in range(NJ)] + [sp.Symbol(f'K{k}') for k in range(NJ)] +
            [sp.Symbol(f'O{k}') for k in range(NJ)] + [sp.Symbol(f'F{k}') for k in range(1, NJ)] +
            [sp.Symbol('s_')])
def coef_fn(expr_str):
    return sp.lambdify(JET_SYMS, sp.sympify(expr_str), 'mpmath')
CM = {k: coef_fn(v) for k, v in cm.items()}
CH = {k: coef_fn(v) for k, v in ch.items()}

class Packet:
    """Numeric first- and second-order fields for profiles G, K (sympy in S)."""
    def __init__(self, Gs, Ks):
        self.GJ = jets_fn(sp.sympify(Gs)); self.KJ = jets_fn(sp.sympify(Ks))
    def jetvals(self, sv):
        o = OmJ(sv); f = PhpJ(sv)   # f[k] = Phi^{(k+1)}
        return self.GJ(sv) + self.KJ(sv) + o + f[:NJ - 1] + [sv]
    def radial(self, sv):
        j = self.jetvals(sv)
        return j
    # ---- first order (inherited reconstruction), with D = x.grad applied exactly:
    # q = G(s) H2(x): D q = (2 s G' + 2 G) H2 ;  similarly for p.
    def v_fields(self, x):
        sv = sum(c * c for c in x)
        G = self.GJ(sv); K = self.KJ(sv); Om = OmJ(sv)[0]; Pp = PhpJ(sv)[0]
        H2 = x[0]**2 - x[1]**2
        q = G[0] * H2; p = K[0] * H2
        Dq = (2 * sv * G[1] + 2 * G[0]) * H2; Dp = (2 * sv * K[1] + 2 * K[0]) * H2
        h = [[2 * q * (sv * (1 if i == j else 0) - 3 * x[i] * x[j]) for j in range(4)] for i in range(4)]
        dphi = -mp.mpf(3) / 2 * (Dq + 4 * q) / (Om**2 * Pp)
        dpi = [[-mp.mpf(3) / 2 * p / Om * x[i] * x[j] for j in range(4)] for i in range(4)]
        dvarpi = -mp.mpf(3) / 2 * Om * (Dp + 5 * p) / Pp
        return h, dphi, dpi, dvarpi
    def w_fields(self, x):
        sv = sum(c * c for c in x)
        jv = self.jetvals(sv)
        X1, X2 = x[0]**2, x[1]**2
        Q = (X1 - X2)**2
        Pm = [[sv * (1 if i == j else 0) - x[i] * x[j] for j in range(4)] for i in range(4)]
        beta = CM['b2'](*jv); nu = CM['n2'](*jv)
        pi2 = [[beta * Q * Pm[i][j] for j in range(4)] for i in range(4)]
        varpi2 = nu * Q
        C1 = CH['c1'](*jv); C2 = CH['c2'](*jv)
        # chi = C1 (x1^2 + x2^2) + C2 (x1^2 - x2^2)^2 with C's frozen (Hessian in x only)
        Hs = [[mp.mpf(0)] * 4 for _ in range(4)]
        # Hessian of (x1^2+x2^2): diag(2,2,0,0); of (x1^2-x2^2)^2:
        u_ = X1 - X2
        Hs[0][0] = 2 * C1 + C2 * (8 * X1 + 4 * u_)
        Hs[1][1] = 2 * C1 + C2 * (8 * X2 - 4 * u_)
        Hs[0][1] = Hs[1][0] = C2 * (-8 * x[0] * x[1])
        PHP = [[sum(Pm[i][k] * Hs[k][l] * Pm[l][j] for k in range(4) for l in range(4))
                for j in range(4)] for i in range(4)]
        tr = sum(PHP[i][i] for i in range(4))
        hT = [[sv * PHP[i][j] - tr * Pm[i][j] / 3 for j in range(4)] for i in range(4)]
        return hT, pi2, varpi2

# ------------------------------------------- truncated eps-series arithmetic
N = 3  # keep eps^0, eps^1, eps^2
def smul(a, b):
    return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(N)]
def sadd(*args):
    return [sum(a[k] for a in args) for k in range(N)]
def sscale(c, a):
    return [c * v for v in a]
def sinv(a):
    a0 = a[0]; r0 = 1 / a0
    return [r0, -a[1] * r0**2, (a[1]**2 / a0 - a[2]) * r0**2]
ZERO = [mp.mpf(0)] * N
PAIRS = [(i, j) for i in range(4) for j in range(i, 4)]

def get(D, key):
    v = D.get(key, 0)
    return mp.mpf(v) if not isinstance(v, mp.mpf) else v

def constraints_at(D, sval):
    n = 4
    def gser(i, j):
        i, j = min(i, j), max(i, j)
        return [get(D, ('g', i, j, o)) for o in range(N)]
    def dgser(i, j, k):
        i, j = min(i, j), max(i, j)
        return [get(D, ('dg', i, j, o, k)) for o in range(N)]
    def ddgser(i, j, k, l):
        i, j = min(i, j), max(i, j); k, l = min(k, l), max(k, l)
        return [get(D, ('ddg', i, j, o, k, l)) for o in range(N)]
    g = [[gser(i, j) for j in range(n)] for i in range(n)]
    dg = [[[dgser(i, j, k) for k in range(n)] for j in range(n)] for i in range(n)]
    ddg = [[[[ddgser(i, j, k, l) for l in range(n)] for k in range(n)] for j in range(n)] for i in range(n)]
    G0 = mp.matrix([[g[i][j][0] for j in range(n)] for i in range(n)])
    A = G0**-1
    G1 = mp.matrix([[g[i][j][1] for j in range(n)] for i in range(n)])
    G2 = mp.matrix([[g[i][j][2] for j in range(n)] for i in range(n)])
    B1 = -A * G1 * A
    B2 = A * G1 * A * G1 * A - A * G2 * A
    gi = [[[A[i, j], B1[i, j], B2[i, j]] for j in range(n)] for i in range(n)]
    def mm(Xs, Ys):
        return [[sadd(*[smul(Xs[i][k], Ys[k][j]) for k in range(n)]) for j in range(n)] for i in range(n)]
    dgi = []
    for k in range(n):
        dgk = [[dg[i][j][k] for j in range(n)] for i in range(n)]
        t1 = mm(mm(gi, dgk), gi)
        dgi.append([[sscale(-1, t1[i][j]) for j in range(n)] for i in range(n)])
    half = mp.mpf(1) / 2
    Lw = [[[sscale(half, sadd(dg[l][j][i], dg[l][i][j], sscale(-1, dg[i][j][l])))
            for j in range(n)] for i in range(n)] for l in range(n)]
    def dlow(m, l, i, j):
        return sscale(half, sadd(ddg[l][j][i][m], ddg[l][i][j][m], sscale(-1, ddg[i][j][l][m])))
    Gam = [[[sadd(*[smul(gi[k][l], Lw[l][i][j]) for l in range(n)]) for j in range(n)]
            for i in range(n)] for k in range(n)]
    def dGam(m, k, i, j):
        return sadd(*[sadd(smul(dgi[m][k][l], Lw[l][i][j]), smul(gi[k][l], dlow(m, l, i, j)))
                      for l in range(n)])
    Ric = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            acc = ZERO
            for k in range(n):
                acc = sadd(acc, dGam(k, k, i, j), sscale(-1, dGam(j, k, i, k)))
                for l in range(n):
                    acc = sadd(acc, smul(Gam[k][k][l], Gam[l][i][j]),
                               sscale(-1, smul(Gam[k][j][l], Gam[l][i][k])))
            Ric[i][j] = Ric[j][i] = acc
    R = sadd(*[smul(gi[i][j], Ric[i][j]) for i in range(n) for j in range(n)])
    d0 = mp.det(G0)
    X1m = A * G1; X2m = A * G2
    tr = lambda M_: sum(M_[i, i] for i in range(n))
    det = [d0, d0 * tr(X1m), d0 * (tr(X2m) + tr(X1m)**2 / 2 - tr(X1m * X1m) / 2)]
    idet = sinv(det)
    def piser(i, j):
        i, j = min(i, j), max(i, j)
        return [get(D, ('pi', i, j, o)) for o in range(N)]
    def dpiser(i, j, k):
        i, j = min(i, j), max(i, j)
        return [get(D, ('dpi', i, j, o, k)) for o in range(N)]
    pi = [[piser(i, j) for j in range(n)] for i in range(n)]
    trpi = sadd(*[smul(g[i][j], pi[i][j]) for i in range(n) for j in range(n)])
    sq = ZERO
    for i, j, k, l in itertools.product(range(n), repeat=4):
        sq = sadd(sq, smul(smul(g[i][k], g[j][l]), smul(pi[i][j], pi[k][l])))
    vp = [get(D, ('vp', o)) for o in range(N)]
    dphi = [[get(D, ('dphi', k, o)) for o in range(N)] for k in range(n)]
    grad2 = sadd(*[smul(gi[i][j], smul(dphi[i], dphi[j])) for i in range(n) for j in range(n)])
    U0 = UF[0](sval); U1 = UF[1](sval); U2 = UF[2](sval)
    dl = [mp.mpf(0), DPHI[1], DPHI[2]]
    Upot = sadd([U0, 0, 0], sscale(U1, dl), sscale(U2 / 2, smul(dl, dl)))
    H = sadd(sscale(half, R),
             sscale(2, smul(idet, sadd(sscale(mp.mpf(1) / 3, smul(trpi, trpi)), sscale(-1, sq)))),
             sscale(-half, smul(idet, smul(vp, vp))),
             sscale(-half, grad2), sscale(-1, Upot))
    M = []
    for i in range(n):
        acc = ZERO
        for j in range(n):
            for k in range(n):
                acc = sadd(acc, smul(dpiser(j, k, j), g[k][i]), smul(pi[j][k], dg[k][i][j]))
                for l in range(n):
                    acc = sadd(acc, sscale(-1, smul(Gam[k][j][i], smul(pi[j][l], g[l][k]))))
        M.append(sadd(sscale(-2, acc), smul(vp, dphi[i])))
    return H, M

DPHI = [0, 0, 0]
def UF0(sv): return U0f(sv)
UF = [U0f, U1f, U2f]

def numeric_D(field_fn, x0, need):
    """Build the derivative dictionary D by high-precision numerical
    differentiation of field_fn(x) -> (g_c, pi_c, phi_c, varpi_c), where
    each *_c entry is a list of eps-coefficients (numbers)."""
    D = {}
    def comp(xs, kind, i, j, o):
        g_c, pi_c, phi_c, vp_c = field_fn(list(xs))
        if kind == 'g': return g_c[i][j][o]
        if kind == 'pi': return pi_c[i][j][o]
        if kind == 'phi': return phi_c[o]
        if kind == 'vp': return vp_c[o]
    cache = {}
    _ff = field_fn
    def field_fn(xs):
        key = tuple(mp.mpf(c)._mpf_ for c in xs) + (mp.mp.prec,)
        if key not in cache:
            cache[key] = _ff(xs)
        return cache[key]
    base_ = field_fn(list(x0))
    for (i, j) in PAIRS:
        for o in range(N):
            D[('g', i, j, o)] = base_[0][i][j][o]
            D[('pi', i, j, o)] = base_[1][i][j][o]
    for o in range(N):
        D[('vp', o)] = base_[3][o]
    with mp.workdps(mp.mp.dps + 30):
        for (i, j) in PAIRS:
            for o in need['g']:
                f = lambda *xs, i=i, j=j, o=o: comp(xs, 'g', i, j, o)
                for k in range(4):
                    oc = [0] * 4; oc[k] = 1
                    D[('dg', i, j, o, k)] = mp.diff(f, tuple(x0), tuple(oc))
                    for l in range(k, 4):
                        oc2 = [0] * 4; oc2[k] += 1; oc2[l] += 1
                        D[('ddg', i, j, o, k, l)] = mp.diff(f, tuple(x0), tuple(oc2))
            for o in need['pi']:
                f = lambda *xs, i=i, j=j, o=o: comp(xs, 'pi', i, j, o)
                for k in range(4):
                    oc = [0] * 4; oc[k] = 1
                    D[('dpi', i, j, o, k)] = mp.diff(f, tuple(x0), tuple(oc))
        sv0 = sum(c * c for c in x0)
        for k in range(4):
            D[('dphi', k, 0)] = 2 * x0[k] * PhpJ(sv0)[0]   # exact gradient of phi0
            for o in (1, 2):
                f = lambda *xs, o=o: comp(xs, 'phi', 0, 0, o)
                oc = [0] * 4; oc[k] = 1
                D[('dphi', k, o)] = mp.diff(f, tuple(x0), tuple(oc))
    DPHI[1] = base_[2][1]; DPHI[2] = base_[2][2]
    return D

def phi0_val(x):
    # the background scalar enters only through its gradient; use the exact
    # radial antiderivative along s numerically: phi0 = int_0^s Phi'(S) dS
    sv = sum(c * c for c in x)
    return mp.quad(lambda t: PhpJ(t)[0], [0, sv])

def make_field(pk, mode):
    def field(x):
        x = [mp.mpf(c) for c in x]
        sv = sum(c * c for c in x)
        Om2 = OmJ(sv)[0]**2
        if mode == 'v':
            h, dphi, dpi, dvarpi = pk.v_fields(x)
            g_c = [[[Om2 * (1 if i == j else 0), h[i][j], 0] for j in range(4)] for i in range(4)]
            pi_c = [[[0, dpi[i][j], 0] for j in range(4)] for i in range(4)]
            phi_c = [0, dphi, 0]
            vp_c = [0, dvarpi, 0]
        else:
            hT, pi2, varpi2 = pk.w_fields(x)
            g_c = [[[Om2 * (1 if i == j else 0), hT[i][j], 0] for j in range(4)] for i in range(4)]
            pi_c = [[[0, pi2[i][j], 0] for j in range(4)] for i in range(4)]
            phi_c = [0, 0, 0]
            vp_c = [0, varpi2, 0]
        return g_c, pi_c, phi_c, vp_c
    return field

def evaluate(pk, mode, xs):
    need = {'g': [0, 1], 'pi': [1]}
    D = numeric_D(make_field(pk, mode), xs, need)
    return constraints_at(D, sum(mp.mpf(c)**2 for c in xs))

def rotation_matrix(seed):
    import random
    rnd = random.Random(seed)
    M_ = mp.matrix(4, 4)
    for i in range(4):
        for j in range(4):
            M_[i, j] = mp.mpf(rnd.uniform(-1, 1))
    Q, _ = mp.qr(M_)
    return Q

def run_case(label, Gs, Ks, points_per_sphere=3, radii=(mp.mpf('0.7'), mp.mpf('1.1'))):
    t1 = time.time()
    pk = Packet(Gs, Ks)
    res = {'spheres': []}
    for r in radii:
        base_pt = mp.matrix([mp.mpf('0.31'), mp.mpf('-0.52'), mp.mpf('0.44'), mp.mpf('0.17')])
        base_pt = base_pt / mp.norm(base_pt) * r
        sph = {'r': mp.nstr(r, 6), 'points': []}
        for kpt in range(points_per_sphere):
            Rm = rotation_matrix(1000 * kpt + 7) if kpt else mp.eye(4)
            v = Rm * base_pt
            xs = [v[i] for i in range(4)]
            Hv, Mv = evaluate(pk, 'v', xs)
            SH = -2 * Hv[2]; SM = [-2 * m[2] for m in Mv]
            fo = [Hv[1]] + [m[1] for m in Mv]
            Hw, Mw = evaluate(pk, 'w', xs)
            DH = Hw[1]; DM = [m[1] for m in Mw]
            scaleM = max(abs(c) for c in SM) + mp.mpf('1e-40')
            sph['points'].append({
                'x': [mp.nstr(c, 12) for c in xs],
                'background_H': mp.nstr(Hv[0], 5),
                'first_order_residual_max': mp.nstr(max(abs(c) for c in fo), 5),
                'S_H': mp.nstr(SH, 20),
                'S_M_scale': mp.nstr(scaleM, 8),
                'momentum_residual_rel': mp.nstr(max(abs(DM[i] - SM[i]) for i in range(4)) / scaleM, 5),
                'hamiltonian_l0_remainder': mp.nstr(SH - DH, 25),
            })
            print(label, 'r', mp.nstr(r, 3), 'pt', kpt, sph['points'][-1], '%.0fs' % (time.time() - t1), flush=True)
        rem = [mp.mpf(p_['hamiltonian_l0_remainder']) for p_ in sph['points']]
        mean = sum(rem) / len(rem)
        sph['remainder_spread_rel'] = mp.nstr((max(rem) - min(rem)) / (abs(mean) + mp.mpf('1e-40')), 5)
        sph['f0_from_remainder'] = mp.nstr(mean, 20)
        # independent hand formula for the velocity-only l=0 density
        sv = r**2
        o = OmJ(sv)[0]; f = PhpJ(sv)[0]; K = pk.KJ(sv)
        sph['f0_velocity_hand_formula'] = mp.nstr(
            sv**2 / 6 * (6 * K[0]**2 * sv**2 + mp.mpf(9) / 4 * (2 * sv * K[1] + 7 * K[0])**2 / f**2) / o**6, 20)
        res['spheres'].append(sph)
    return res

if __name__ == '__main__':
    out['cases'] = {}
    out['cases']['commuting_A_eq_B'] = run_case('cc', Gf, Kf)
    out['cases']['pure_velocity_G0'] = run_case('vel', sp.Integer(0), Kf, points_per_sphere=2)
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(here, 'indep_check.json'), 'w'), indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != 'cases'}, indent=2))
