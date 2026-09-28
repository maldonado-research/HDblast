#!/usr/bin/env python3
"""Independent symbolic audit of the continuum equations used by the workstream.

Written from scratch (does not import the workstream's controls.py):
 1. Einstein tensor of ds^2 = e^{2B}(-dt^2+dz^2) + e^{2A} dx_3^2 computed from the metric,
    with action normalisation R/2 - (1/2)(d phi)^2 - U(phi) (kappa_5^2 = 1).
    Checks that the workstream's H and M are fixed multiples of (G-T)_tt and (G-T)_tz,
    and that the evolution equations EA, EB, EP reproduce the remaining field equations
    modulo H (i.e. are a legitimate free-evolution form).
 2. Jet-space derivation (own total-derivative operators, recursive elimination of
    second time derivatives) of the damped constraint-transport law
    (d_t - d_z)[e^{3A}C+] = -kappa e^{3A}C+, (d_t + d_z)[e^{3A}C-] = -kappa e^{3A}C+(A_t-A_z)/(A_t+A_z),
    with kappa allowed to depend on z; and the undamped law as a known limit.
 3. Negative controls: wrong sign, wrong denominator, and the naive H-only damping.
 4. Constant-coefficient Fourier analysis of the constraint subsystem for the H-only and C+ damping.
 5. Perturbation-variable mapping: the lab's discrete H, M, damping denominators equal the
    continuum ones with A = t + ln rho + a, B = ln rho + b, background A_z = B_z = rho_z/rho = hc.
 6. Static junction data: the lab/frozen boundary slopes ga, gf equal the Israel-type conditions
    e^{-B} A_z = sigma(phi)/6, e^{-B} phi_z = -sigma'(phi)/2, expanded exactly in (b, f), and their
    time derivatives gpa, gpf.
 7. BPS consistency: U = W_phi^2/2 - 2W^2/3 with phi_y = W_phi, A_y = -W/3 solves the static
    equations in this normalisation (fine-tuned flat-brane limit), and sigma = 2W matches.
Writes sympy_independent.json.
"""
import json, time
from pathlib import Path
import sympy as sp

OUT = Path(__file__).resolve().parent / 'sympy_independent.json'
res = dict(status='PASS', checks=[])
def check(name, ok, **kw):
    res['checks'].append(dict(name=name, pass_=bool(ok), **{k: str(v) for k, v in kw.items()}))
    if not ok: res['status'] = 'FAIL'
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)

t0 = time.time()
# ------------------------------------------------------------------ 1. Einstein tensor
t, z, x1, x2, x3 = sp.symbols('t z x1 x2 x3', real=True)
A = sp.Function('A')(t, z); B = sp.Function('B')(t, z); p = sp.Function('phi')(t, z)
Uf = sp.Function('U')
X = [t, z, x1, x2, x3]
g = sp.diag(-sp.exp(2 * B), sp.exp(2 * B), sp.exp(2 * A), sp.exp(2 * A), sp.exp(2 * A))
gi = g.inv()
n = 5
Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n)) / 2
         for c in range(n)] for b in range(n)] for a in range(n)]
def ricci(b, c):
    return sp.simplify(sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
                           + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a] for d in range(n))
                           for a in range(n)))
Ric = sp.Matrix(n, n, lambda b, c: ricci(b, c) if (b == c or (b, c) in [(0, 1), (1, 0)]) else 0)
Rs = sp.simplify(sum(gi[a, a] * Ric[a, a] for a in range(n)))
G = (Ric - Rs * g / 2).applyfunc(sp.simplify)
dp = [sp.diff(p, xx) for xx in X]
kin = sum(gi[a, a] * dp[a] ** 2 for a in range(n))
T = sp.Matrix(n, n, lambda a, b: dp[a] * dp[b] - g[a, b] * (kin / 2 + Uf(p)))
E = (G - T).applyfunc(sp.expand)

At, Az, Bt, Bz, pt, pz = [sp.diff(f, v) for f in (A, B, p) for v in (t, z)]
Azz = sp.diff(A, z, 2)
H = -2 * sp.exp(2 * B) * Uf(p) + 6 * At ** 2 + 6 * At * Bt - 12 * Az ** 2 + 6 * Az * Bz - 6 * Azz - pt ** 2 - pz ** 2
M = -3 * sp.diff(A, t, z) - 3 * At * Az + 3 * At * Bz + 3 * Az * Bt - pt * pz
kH = sp.simplify(H / E[0, 0]); kM = sp.simplify(M / E[0, 1])
check('H is a constant multiple of (G-T)_tt', kH.is_number, factor=kH)
check('M is a constant multiple of (G-T)_tz', kM.is_number, factor=kM)
EA = Azz - 3 * At ** 2 + 3 * Az ** 2 + sp.Rational(2, 3) * sp.exp(2 * B) * Uf(p)
EB = sp.diff(B, z, 2) + 3 * At ** 2 - 3 * Az ** 2 - pt ** 2 / 2 + pz ** 2 / 2 - sp.exp(2 * B) * Uf(p) / 3
EP = sp.diff(p, z, 2) - 3 * At * pt + 3 * Az * pz - sp.exp(2 * B) * sp.diff(Uf(p), p)
sub2 = {sp.diff(A, t, 2): EA, sp.diff(B, t, 2): EB, sp.diff(p, t, 2): EP}
for (i, j), nm in [((1, 1), 'zz'), ((2, 2), 'xx')]:
    e = sp.expand(E[i, j].subs(sub2))
    # must be proportional to H with a coefficient free of second time derivatives
    c = sp.simplify(e / H)
    ok = sp.simplify(e - c * H) == 0 and not c.has(sp.Derivative)
    check(f'(G-T)_{nm} vanishes modulo H on the evolution equations', ok, coefficient=c)
box = sp.exp(-2 * B) * (-sp.diff(p, t, 2) - 3 * At * pt + sp.diff(p, z, 2) + 3 * Az * pz)
check('scalar equation box(phi)=U_phi reproduces EP', sp.simplify(sp.solve(sp.Eq(box, sp.diff(Uf(p), p)), sp.diff(p, t, 2))[0] - EP) == 0)
# wrong-normalisation control: a kinetic factor 2 in T breaks the proportionality of H
T2 = sp.Matrix(n, n, lambda a, b: 2 * dp[a] * dp[b] - g[a, b] * (kin + Uf(p)))
check('control: doubled scalar kinetic normalisation breaks H proportionality',
      not sp.simplify(H / sp.expand(G[0, 0] - T2[0, 0])).is_number)
res['einstein_factors'] = dict(H_over_Ett=str(kH), M_over_Etz=str(kM))

# ------------------------------------------------------------------ 2. jet-space transport law
N = 6
Aj = [[sp.Symbol(f'A_{i}_{j}') for j in range(N)] for i in range(N)]
Bj = [[sp.Symbol(f'B_{i}_{j}') for j in range(N)] for i in range(N)]
Pj = [[sp.Symbol(f'P_{i}_{j}') for j in range(N)] for i in range(N)]
Ufun = sp.Function('U'); P0 = Pj[0][0]
kz = [sp.Symbol(f'k_{j}') for j in range(N)]   # kappa(z) and its z-derivatives
fields = {}
for F, tag in ((Aj, 'A'), (Bj, 'B'), (Pj, 'P')):
    for i in range(N):
        for j in range(N): fields[F[i][j]] = (F, i, j)
def Dz(e):
    out = 0
    for s in e.free_symbols:
        if s in fields:
            F, i, j = fields[s]; out += sp.diff(e, s) * F[i][j + 1]
        elif s in kz:
            out += sp.diff(e, s) * kz[kz.index(s) + 1]
    for u in e.atoms(sp.Function):
        pass
    # U(P0) and its derivatives depend on z through P0 (chain rule handled by sympy via P0 symbol)
    return out
def Dt(e):
    out = 0
    for s in e.free_symbols:
        if s in fields:
            F, i, j = fields[s]; out += sp.diff(e, s) * F[i + 1][j]
    return out
def ev_rhs(damp):
    a = Aj; b = Bj; q = Pj
    U0 = Ufun(P0); U1 = sp.diff(U0, P0)
    eA = a[0][2] - 3 * a[1][0] ** 2 + 3 * a[0][1] ** 2 + sp.Rational(2, 3) * sp.exp(2 * b[0][0]) * U0
    eB = b[0][2] + 3 * a[1][0] ** 2 - 3 * a[0][1] ** 2 - q[1][0] ** 2 / 2 + q[0][1] ** 2 / 2 - sp.exp(2 * b[0][0]) * U0 / 3 + damp
    eP = q[0][2] - 3 * a[1][0] * q[1][0] + 3 * a[0][1] * q[0][1] - sp.exp(2 * b[0][0]) * U1
    return eA, eB, eP
def reduce_tt(e, rules):
    """Recursively eliminate every jet variable with >=2 time derivatives."""
    for _ in range(12):
        hit = False
        for s in sorted(e.free_symbols, key=str):
            if s in fields and fields[s][1] >= 2:
                F, i, j = fields[s]
                base = rules[id(F)]
                r = base
                for _k in range(i - 2): r = Dt(r)
                for _k in range(j): r = Dz(r)
                e = e.subs(s, r); hit = True
        if not hit: return sp.expand(e)
    raise RuntimeError('reduction did not terminate')
def reduce_full(e, rules):
    # Dt of a reduced expression can reintroduce second time derivatives; iterate
    return reduce_tt(sp.expand(e), rules)
a, b, q = Aj, Bj, Pj
U0 = Ufun(P0)
Hj = -2 * sp.exp(2 * b[0][0]) * U0 + 6 * a[1][0] ** 2 + 6 * a[1][0] * b[1][0] - 12 * a[0][1] ** 2 + 6 * a[0][1] * b[0][1] - 6 * a[0][2] - q[1][0] ** 2 - q[0][1] ** 2
Mj = -3 * a[1][1] - 3 * a[1][0] * a[0][1] + 3 * a[1][0] * b[0][1] + 3 * a[0][1] * b[1][0] - q[1][0] * q[0][1]
Cp, Cm = Hj + 2 * Mj, Hj - 2 * Mj
w = sp.exp(3 * a[0][0]); k = kz[0]
def laws(damp):
    eA, eB, eP = ev_rhs(damp)
    rules = {id(Aj): eA, id(Bj): eB, id(Pj): eP}
    # For the Dt-chain inside reduce_tt the rule for F[i][j] with i>2 needs Dt of the rule, which
    # itself may contain second time derivatives -> handled by the recursion loop.
    lp = reduce_full(Dt(w * Cp) - Dz(w * Cp) + k * w * Cp, rules)
    lm = reduce_full(Dt(w * Cm) + Dz(w * Cm) + k * w * Cp * (a[1][0] - a[0][1]) / (a[1][0] + a[0][1]), rules)
    return sp.simplify(lp), sp.simplify(lm)
lp0, lm0 = laws(0)
lp0k = sp.simplify(lp0.subs(k, 0)); lm0k = sp.simplify(lm0.subs(k, 0))
check('undamped limit: (d_t - d_z)[e^{3A}(H+2M)] = 0 on the evolution equations', lp0k == 0)
check('undamped limit: (d_t + d_z)[e^{3A}(H-2M)] = 0 on the evolution equations', lm0k == 0)
damp = -k * Cp / (6 * (a[1][0] + a[0][1]))
lp, lm = laws(damp)
check('damped law C+ with z-dependent kappa (jet-space derivation)', lp == 0)
check('damped law C- with z-dependent kappa (jet-space derivation)', lm == 0)
lpw, _ = laws(+k * Cp / (6 * (a[1][0] + a[0][1])))
check('negative control: wrong sign violates the C+ law', lpw != 0)
lpd, _ = laws(-k * Cp / (6 * (a[1][0] - a[0][1])))
check('negative control: denominator A_t-A_z violates the C+ law', lpd != 0)
# H-only damping: derive the source terms it generates
eA, eB, eP = ev_rhs(sp.Symbol('S'))
rules = {id(Aj): eA, id(Bj): eB, id(Pj): eP}
S = sp.Symbol('S')
dCp = reduce_full(Dt(w * Cp) - Dz(w * Cp), rules); dCm = reduce_full(Dt(w * Cm) + Dz(w * Cm), rules)
check('source S in B_tt enters (d_t-d_z)[e^3A C+] as 6 e^3A (A_t+A_z) S', sp.simplify(sp.diff(dCp, S) - 6 * w * (a[1][0] + a[0][1])) == 0)
check('source S in B_tt enters (d_t+d_z)[e^3A C-] as 6 e^3A (A_t-A_z) S', sp.simplify(sp.diff(dCm, S) - 6 * w * (a[1][0] - a[0][1])) == 0)

# ------------------------------------------------------------------ 4. Fourier analysis
kk, kap, r = sp.symbols('k kappa r', positive=True)
lam = sp.Symbol('lam')
# constraint subsystem for (C+, C-) with frozen coefficients, r = A_z/A_t, e^{3A} weight absorbed;
# (d_t - d_z) C+ = source_+, (d_t + d_z) C- = source_-; Fourier e^{ikz}: d_t C+ = ik C+ + ..., d_t C- = -ik C- + ...
def growth(src_p, src_m):
    Mx = sp.Matrix([[sp.I * kk + sp.diff(src_p, 'Cp'), sp.diff(src_p, 'Cm')],
                    [sp.diff(src_m, 'Cp'), -sp.I * kk + sp.diff(src_m, 'Cm')]])
    return Mx
Cps, Cms = sp.symbols('Cp Cm')
Hs = (Cps + Cms) / 2
# H-only: S = -kappa H/(6 A_t): source+ = 6(A_t+A_z)S = -kappa H (1+r), source- = -kappa H (1-r)
Mh = growth(-kap * Hs * (1 + r), -kap * Hs * (1 - r))
num = []
for rv in (0.5, 2.0, 17.0):
    Mn = sp.lambdify((kk,), Mh.subs({kap: 10, r: rv}), 'numpy')
    import numpy as _np
    eg = _np.linalg.eigvals(_np.array(Mn(1e5), dtype=complex))
    num.append(dict(r=rv, max_re=float(eg.real.max()), predicted=float(10 * (rv - 1) / 2) if rv > 1 else float(max(10 * (rv - 1) / 2, -10 * (1 + rv) / 2))))
res['honly_high_k_growth_rates'] = num
check('H-only damping: high-k max growth rate equals kappa(r-1)/2 for r = A_z/A_t > 1 (numeric eigenvalues, k=1e5, kappa=10)',
      all(abs(d['max_re'] - d['predicted']) < 1e-3 for d in num if d['r'] > 1) and num[0]['max_re'] < 0)
Mc = growth(-kap * Cps, -kap * Cps * (1 - r) / (1 + r))
evc = list(Mc.eigenvals())
check('C+ damping: frozen-coefficient eigenvalues ik-kappa and -ik (no growth)',
      set(sp.simplify(e) for e in evc) == {sp.I * kk - kap, -sp.I * kk})

# ------------------------------------------------------------------ 5. perturbation mapping
rho = sp.Function('rho')(z); a_, b_, f_ = [sp.Function(s_)(t, z) for s_ in ('a', 'b', 'f')]
ph0 = sp.Function('phi0')(z)
Ab = t + sp.log(rho) + a_; Bb = sp.log(rho) + b_; pb_ = ph0 + f_
hc = sp.diff(rho, z) / rho; phz = sp.diff(ph0, z)
rep = {A: Ab, B: Bb, p: pb_}
Hc = H.subs(rep).doit(); Mc_ = M.subs(rep).doit()
pa, pbv, pf = sp.diff(a_, t), sp.diff(b_, t), sp.diff(f_, t)
az, bz, fz = sp.diff(a_, z), sp.diff(b_, z), sp.diff(f_, z)
azz, paz = sp.diff(a_, z, 2), sp.diff(a_, t, z)
# lab formulae (su = rho^2 (e^{2b}U(phi0+f) - U(phi0)))
su = rho ** 2 * (sp.exp(2 * b_) * Uf(ph0 + f_) - Uf(ph0))
H_lab = (-2 * su + 12 * pa + 6 * pa * pa + 6 * (1 + pa) * pbv - 18 * hc * az + 6 * hc * bz
         - 12 * az * az + 6 * az * bz - 6 * azz - pf * pf - 2 * phz * fz - fz * fz)
M_lab = -3 * paz - 3 * (1 + pa) * (az - bz) + 3 * (hc + az) * pbv - pf * (phz + fz)
bgH = sp.simplify(Hc - H_lab)          # should be the background Hamiltonian constraint only
bgH0 = bgH.subs({a_: 0, b_: 0, f_: 0}).doit()
check('lab H = continuum H minus the static-background H (difference independent of a,b,f)',
      sp.simplify(bgH - bgH0) == 0, background_part=sp.simplify(bgH0))
check('lab M = continuum M exactly (background M vanishes)', sp.simplify(Mc_ - M_lab) == 0)
check('lab damping denominator 1+pa+hc+az equals A_t+A_z', sp.simplify((sp.diff(Ab, t) + sp.diff(Ab, z)) - (1 + pa + hc + az)) == 0)
check('lab weight (rho/rb)^3 e^{3(t+a)} equals e^{3A}/rb^3', sp.simplify(sp.exp(3 * Ab) - rho ** 3 * sp.exp(3 * (t + a_))) == 0)

# ------------------------------------------------------------------ 6. junction slopes
bb, ff, pbb, pff, rb, phb, dl, cc = sp.symbols('b f p_b p_f rho_b phi_b delta c', real=True)
Wf = lambda x: 1 - x + x ** 3 / 3
sig = lambda x: 2 * Wf(x) + dl * (1 + cc * x)
sigp = lambda x: sp.diff(sig(sp.Symbol('y')), sp.Symbol('y')).subs(sp.Symbol('y'), x)
hcb = rb * sig(phb) / 6          # static background: e^{-B}A_z = sigma/6 at B=ln rho_b
phzb = -rb * sigp(phb) / 2       # e^{-B}phi_z = -sigma'/2
ga_exact = rb * sp.exp(bb) * sig(phb + ff) / 6 - hcb
gf_exact = -rb * sp.exp(bb) * sigp(phb + ff) / 2 - phzb
s0 = 2 * Wf(phb) + dl * (1 + cc * phb); s10 = 2 * (phb ** 2 - 1) + dl * cc; s20 = 4 * phb
ds = s10 * ff + 2 * phb * ff ** 2 + sp.Rational(2, 3) * ff ** 3; ds1 = s20 * ff + 2 * ff ** 2
ga_lab = rb * (s0 * (sp.exp(bb) - 1) + sp.exp(bb) * ds) / 6
gf_lab = -rb * (s10 * (sp.exp(bb) - 1) + sp.exp(bb) * ds1) / 2
check('ga (lab/frozen) equals e^{-B}A_z = sigma(phi)/6 exactly', sp.simplify(sp.expand(ga_exact - ga_lab)) == 0)
check('gf (lab/frozen) equals e^{-B}phi_z = -sigma_phi/2 exactly', sp.simplify(sp.expand(gf_exact - gf_lab)) == 0)
tt = sp.Symbol('tt'); bt = sp.Function('bt')(tt); ft = sp.Function('ft')(tt)
gpa_exact = sp.diff(ga_exact.subs({bb: bt, ff: ft}), tt).subs({sp.diff(bt, tt): pbb, sp.diff(ft, tt): pff}).subs({bt: bb, ft: ff})
gpf_exact = sp.diff(gf_exact.subs({bb: bt, ff: ft}), tt).subs({sp.diff(bt, tt): pbb, sp.diff(ft, tt): pff}).subs({bt: bb, ft: ff})
gpa_lab = rb * sp.exp(bb) * ((s0 + ds) * pbb + (s10 + ds1) * pff) / 6
gpf_lab = -rb * sp.exp(bb) * ((s10 + ds1) * pbb + 4 * (phb + ff) * pff) / 2
check('gpa equals d/dt of the metric junction slope', sp.simplify(sp.expand(gpa_exact - gpa_lab)) == 0)
check('gpf equals d/dt of the scalar junction slope', sp.simplify(sp.expand(gpf_exact - gpf_lab)) == 0)
# control: using sigma/3 (two-sided jump instead of Z2 one-sided) fails
check('control: sigma/3 metric slope does not match ga', sp.simplify(sp.expand(rb * sp.exp(bb) * sig(phb + ff) / 3 - 2 * hcb - ga_lab)) != 0)

# ------------------------------------------------------------------ 7. BPS consistency in this normalisation
y = sp.Symbol('y'); ph = sp.Function('ph')(y); Ay = sp.Function('Ay')(y)
Wp = lambda x: 1 - x + x ** 3 / 3
UU = lambda x: sp.diff(Wp(sp.Symbol('u')), sp.Symbol('u')).subs(sp.Symbol('u'), x) ** 2 / 2 - sp.Rational(2, 3) * Wp(x) ** 2
# static flat slicing ds^2 = dy^2 + e^{2A(y)} eta_4: 6A'^2 = phi'^2/2 - U, phi'' + 4A'phi' = U_phi
Wph = sp.diff(Wp(ph), ph)
subsB = {sp.diff(ph, y): Wph, sp.diff(Ay, y): -Wp(ph) / 3}
c1 = sp.simplify(6 * (Wp(ph) / 3) ** 2 - (Wph ** 2 / 2 - UU(ph)))
phpp = sp.diff(Wph, ph) * Wph
c2 = sp.simplify(phpp + 4 * (-Wp(ph) / 3) * Wph - sp.diff(UU(ph), ph))
check('BPS: U = W_phi^2/2 - 2W^2/3 solves 6A_y^2 = phi_y^2/2 - U with phi_y=W_phi, A_y=-W/3', c1 == 0)
check('BPS: scalar equation phi_yy + 4A_y phi_y = U_phi holds on the flow', c2 == 0)
check('BPS: at delta=0 sigma=2W gives |A_y| = sigma/6 = W/3 and |phi_y| = sigma_phi/2 = |W_phi|',
      sp.simplify(sig(ph).subs(dl, 0) / 6 - Wp(ph) / 3) == 0 and sp.simplify(sigp(ph).subs(dl, 0) / 2 - sp.diff(Wp(ph), ph)) == 0)
res['seconds'] = time.time() - t0
OUT.write_text(json.dumps(res, indent=1) + '\n')
print(res['status'], res['seconds'])
