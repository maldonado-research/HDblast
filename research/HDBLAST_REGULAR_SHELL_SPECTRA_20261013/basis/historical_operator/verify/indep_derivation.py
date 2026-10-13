"""Independent auditor derivation (does NOT import or copy derive_linearized.py).

Method (different from the audited script):
  * builds the FULL nonlinear Einstein tensor G_AB = R_AB - R g_AB/2 of the perturbed metric (epsilon kept symbolic)
    and the full stress tensor T_AB = d_A f d_B f - g_AB((df)^2/2 + U), then takes d/d(epsilon) at 0;
  * works only in longitudinal gauge (g_yy = 1 + 2 eps N Y, g_mn = rho^2 (1 + 2 eps psi Y) gamma_mn), Y = e^{p t};
  * the displaced shell y = y_b + eps zeta Y is handled via the level-set normal and the explicit
    covariant derivative nabla_A n_B projected on the shell tangent vectors;
  * background relations are obtained by SOLVING the zeroth-order equations here (not typed in).
Checks: background eqs; trace-free part => N = -2 psi; yt constraint; Gauss constraint; that the claimed first-order
(X, Z) system implies every first-order field equation; junction trace-free part forces zeta = 0;
trace junction == momentum constraint; scalar junction == X' + sigma''/2 X + 2 phi' Psi == B X + 3 lam Z / rho^2;
tensor equation and Neumann condition.  Deliberate wrong-formula controls must fail.
Output: indep_derivation.json
"""
import json, time, hashlib
from pathlib import Path
import sympy as sp

T0 = time.time()
y, t, x1, x2, x3 = sp.symbols('y t x1 x2 x3', real=True)
eps, p, mu2, zeta = sp.symbols('epsilon p mu2 zeta')
C = [y, t, x1, x2, x3]
rho = sp.Function('rho')(y); ph = sp.Function('phi')(y)
N = sp.Function('N')(y); psi = sp.Function('psi')(y); chi = sp.Function('chi')(y); hT = sp.Function('hT')(y)
U0, U1, U2, U3 = [sp.Function(n)(y) for n in ('U0', 'U1', 'U2', 'U3')]
Y = sp.exp(p*t); a2 = sp.exp(2*t)
res = {}

def ok(name, expr):
    v = sp.simplify(expr) == 0
    res[name] = bool(v); print(('PASS ' if v else 'FAIL ') + name, flush=True)
    return v

def ricci(g):
    gi = g.inv()
    n = 5
    Gam = [[[sp.expand(sum(gi[a, d]*(sp.diff(g[d, b], C[c]) + sp.diff(g[d, c], C[b]) - sp.diff(g[b, c], C[d])) for d in range(n))/2)
             for c in range(n)] for b in range(n)] for a in range(n)]
    R = sp.zeros(n)
    for b in range(n):
        for c in range(b, n):
            r = sum(sp.diff(Gam[a][b][c], C[a]) - sp.diff(Gam[a][b][a], C[c]) for a in range(n))
            r += sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for a in range(n) for d in range(n))
            R[b, c] = R[c, b] = r
    return R, Gam, gi

def lin(e):
    """first-order coefficient in eps"""
    return sp.diff(e, eps).subs(eps, 0)

# U(phi + eps chi Y) to first order, with d/dy U_k = U_{k+1} phi'
Ufull = U0 + eps*U1*chi*Y
def dsubU(e):
    for _ in range(3):
        e = e.subs({sp.Derivative(U0, y): U1*ph.diff(y), sp.Derivative(U1, y): U2*ph.diff(y), sp.Derivative(U2, y): U3*ph.diff(y)})
    return e

# ------------------------------------------------------------------ metric and field equations
gam = sp.diag(-1, a2, a2, a2)
def metric(Nf, psif, hten=0):
    g = sp.zeros(5)
    g[0, 0] = 1 + 2*eps*Nf*Y
    for m in range(4):
        g[m+1, m+1] = rho**2*(1 + 2*eps*psif*Y)*gam[m, m]
    if hten != 0:
        g[2, 3] = g[3, 2] = rho**2*a2*eps*hten*Y
    return g

def field_eqs(g, fld, Uf):
    R, Gam, gi = ricci(g)
    Rs = sum(gi[a, b]*R[a, b] for a in range(5) for b in range(5))
    df = [sp.diff(fld, c) for c in C]
    kin = sum(gi[a, b]*df[a]*df[b] for a in range(5) for b in range(5))
    E = sp.zeros(5)
    for a in range(5):
        for b in range(5):
            E[a, b] = R[a, b] - Rs*g[a, b]/2 - (df[a]*df[b] - g[a, b]*(kin/2 + Uf))
    sqrtg = sp.sqrt(-g.det())
    box = sum(sp.diff(sqrtg*gi[a, b]*df[b], C[a]) for a in range(5) for b in range(5))/sqrtg
    return E, box, Gam, gi

# zeroth order
g0 = metric(0, 0)
E0, box0, _, _ = field_eqs(g0, ph, U0)
E0 = E0.applyfunc(lambda z: sp.simplify(dsubU(z)))
S0 = sp.simplify(box0 - U1)
# solve background: from yy constraint -> U0 ; from xx -> rho'' ; scalar -> phi''
Hs = rho.diff(y)/rho
U0sol = sp.solve(E0[0, 0], U0)[0]
rpp = sp.solve(sp.simplify((E0[2, 2]/a2).subs(U0, U0sol)), rho.diff(y, 2))[0]
fpp = sp.solve(S0, ph.diff(y, 2))[0]
print('U0 =', sp.simplify(U0sol)); print("rho'' =", sp.simplify(rpp)); print("phi'' =", sp.simplify(fpp))
ok('bg_first_integral_is_rho1^2=1+rho^2(phi1^2/12-U/6)', sp.expand((1 + rho**2*(ph.diff(y)**2/12 - U0sol/6)) - rho.diff(y)**2))
ok("bg_scalar_is_phi''=U1-4Hphi'", fpp - (U1 - 4*Hs*ph.diff(y)))
ok('bg_tt_consistent', sp.simplify((E0[1, 1]).subs(U0, U0sol).subs(rho.diff(y, 2), rpp)))

def bg(e):
    for _ in range(4):
        e = dsubU(e)
        e = e.subs({ph.diff(y, 3): sp.diff(fpp, y), rho.diff(y, 3): sp.diff(rpp, y)})
        e = dsubU(e)
        e = e.subs({ph.diff(y, 2): fpp, rho.diff(y, 2): rpp})
    e = dsubU(e).subs(U0, U0sol)
    return sp.simplify(e)
ok("bg_U1_consistency(d/dy of first integral)", bg(sp.diff(U0sol, y) - U1*ph.diff(y)))

# first order, longitudinal gauge
g = metric(N, psi)
E, box, Gam, gi = field_eqs(g, ph + eps*chi*Y, Ufull)
lam = mu2 + 4
def red(e):
    """remove harmonic factor, apply background, reduce p modulo p^2+3p+mu2"""
    e = bg(sp.expand(e))
    num, den = sp.fraction(sp.together(e))
    num = sp.rem(sp.Poly(sp.expand(num), p), sp.Poly(p**2 + 3*p + mu2, p)).as_expr()
    den = sp.rem(sp.Poly(sp.expand(den), p), sp.Poly(p**2 + 3*p + mu2, p)).as_expr()
    return sp.simplify(num/den)

Eyy = red(lin(E[0, 0])/Y); Eyt = red(lin(E[0, 1])/(p*Y)); Ett = red(lin(E[1, 1])/Y); Exx = red(lin(E[2, 2])/(a2*Y))
Ssc = red(lin(box)/Y - U2*chi)
ok('offdiag_yx_zero', lin(E[0, 2])); ok('offdiag_x1x2_zero', lin(E[2, 3]))
# E_mn = A gamma_mn Y + Bc nabla_m nabla_n Y  with nabla_t nabla_t Y = p^2 Y, nabla_i nabla_j Y = -p a2 Y delta_ij
A_, B_ = sp.symbols('A_ B_')
# before reducing mod p we need unreduced tt, xx:
Ett_u = bg(sp.expand(lin(E[1, 1])/Y)); Exx_u = bg(sp.expand(lin(E[2, 2])/(a2*Y)))
sol = sp.solve([sp.Eq(-A_ + B_*p**2, Ett_u), sp.Eq(A_ - B_*p, Exx_u)], [A_, B_], dict=True)[0]
tracefree = sp.factor(sp.simplify(sol[B_]))
print('trace-free coefficient:', tracefree)
ok('tracefree_proportional_to_N+2psi', sp.simplify(tracefree.subs(N, -2*psi)))
res['tracefree_nonzero_without_N=-2psi'] = bool(sp.simplify(tracefree.subs(N, 0)) != 0)

# momentum constraint
Eyt_l = sp.simplify(Eyt.subs(N, -2*psi).doit())
print('yt:', sp.factor(Eyt_l))
ok("yt_equals_const*(psi'+2Hpsi+phi'chi/3)", sp.simplify(sp.diff(sp.simplify(Eyt_l/(psi.diff(y) + 2*Hs*psi + ph.diff(y)*chi/3)), y)))
res['yt_ratio'] = str(sp.simplify(Eyt_l/(psi.diff(y) + 2*Hs*psi + ph.diff(y)*chi/3)))

# claimed first-order system in (chi, psi):  psi' = -2H psi - phi' chi/3 ; chi' = g chi + (3 lam/(rho^2 phi') - 2 phi') psi
gg = fpp/ph.diff(y)
dpsi = -2*Hs*psi - ph.diff(y)*chi/3
dchi = gg*chi + (3*lam/(rho**2*ph.diff(y)) - 2*ph.diff(y))*psi
first = {psi.diff(y): dpsi, chi.diff(y): dchi}
second = {psi.diff(y, 2): sp.diff(dpsi, y), chi.diff(y, 2): sp.diff(dchi, y)}
def on_system(e):
    e = e.subs(N, -2*psi).doit()
    e = e.subs(second).subs(first)
    e = e.subs(first)
    return bg(sp.expand(e))
for nm, e in [('yy', Eyy), ('yt', Eyt), ('tt', red(Ett_u)), ('xx', red(Exx_u)), ('scalar', Ssc)]:
    ok('claimed_system_implies_' + nm, on_system(e))
# control: wrong sign of the gravitational coupling term 3 lam/(rho^2 phi') must violate the equations
dchi_bad = gg*chi + (-3*lam/(rho**2*ph.diff(y)) - 2*ph.diff(y))*psi
e = Eyy.subs(N, -2*psi).doit().subs({chi.diff(y, 2): sp.diff(dchi_bad, y), psi.diff(y, 2): sp.diff(dpsi, y)}).subs({chi.diff(y): dchi_bad, psi.diff(y): dpsi})
res['CONTROL_wrong_sign_grav_term_detected'] = bool(bg(sp.expand(e)) != 0)
dpsi_bad = -2*Hs*psi + ph.diff(y)*chi/3
e = Eyt.subs(N, -2*psi).doit().subs({psi.diff(y): dpsi_bad})
res['CONTROL_wrong_sign_momentum_constraint_detected'] = bool(bg(sp.expand(e)) != 0)

# (X, Z) with psi = phi' Z, chi = X
Zf = sp.Function('Z')(y); Xf = sp.Function('X')(y)
Zp = sp.simplify(bg((dpsi.subs({psi: ph.diff(y)*Zf, chi: Xf}) - fpp*Zf)/ph.diff(y)))
Xp = sp.simplify(bg(dchi.subs({psi: ph.diff(y)*Zf, chi: Xf})))
ok("Z'=-(2H+g)Z-X/3", Zp - (-(2*Hs + gg)*Zf - Xf/3))
ok("X'=gX+(3lam/rho^2-2phi'^2)Z", bg(Xp - (gg*Xf + (3*lam/rho**2 - 2*ph.diff(y)**2)*Zf)))

# ------------------------------------------------------------------ junctions on the displaced shell
dl, cc, f = sp.symbols('delta c f')
sig = 2*(1 - f + f**3/3) + dl*(1 + cc*f)
s0 = lambda z: sig.subs(f, z); s1 = lambda z: sp.diff(sig, f).subs(f, z); s2 = lambda z: sp.diff(sig, f, 2).subs(f, z)
Fl = y - eps*zeta*Y
dF = [sp.diff(Fl, c) for c in C]
nrm = sp.sqrt(sum(gi[a, b]*dF[a]*dF[b] for a in range(5) for b in range(5)))
nl = [d/nrm for d in dF]
nu = [sum(gi[a, b]*nl[b] for b in range(5)) for a in range(5)]
tang = []
for m in range(4):
    v = [0]*5; v[m+1] = 1; v[0] = eps*zeta*sp.diff(Y, C[m+1]); tang.append(v)
def K(m, n_):
    return sum(tang[m][a]*tang[n_][b]*(sp.diff(nl[b], C[a]) - sum(Gam[c][a][b]*nl[c] for c in range(5))) for a in range(5) for b in range(5))
def h(m, n_):
    return sum(tang[m][a]*tang[n_][b]*g[a, b] for a in range(5) for b in range(5))
fld = ph + eps*chi*Y
Jtt = K(0, 0) - s0(fld)/6*h(0, 0); Jxx = K(1, 1) - s0(fld)/6*h(1, 1)
Jsc = sum(nu[a]*sp.diff(fld, C[a]) for a in range(5)) + s1(fld)/2
yb = sp.Symbol('y_b')
def on_shell(e):
    e0 = e.subs(eps, 0); e1 = lin(e)
    return e0, e1 + zeta*Y*sp.diff(e0, y)
bgj = {rho.diff(y): rho*s0(ph)/6, ph.diff(y): -s1(ph)/2}
def js(e):
    e = bg(sp.expand(e))
    e = e.subs(bgj)
    return sp.simplify(e)
jtt0, jtt1 = on_shell(Jtt); jxx0, jxx1 = on_shell(Jxx); jsc0, jsc1 = on_shell(Jsc)
ok('bg_metric_junction_is_rho1/rho=sigma/6', sp.simplify((jxx0/a2).subs(rho.diff(y), rho*s0(ph)/6)))
ok("bg_scalar_junction_is_phi1=-sigma1/2", sp.simplify(jsc0.subs(ph.diff(y), -s1(ph)/2)))
jtt1 = js(jtt1/Y); jxx1 = js(jxx1/(a2*Y)); jsc1 = js(jsc1/Y)
solJ = sp.solve([sp.Eq(-A_ + B_*p**2, jtt1), sp.Eq(A_ - B_*p, jxx1)], [A_, B_], dict=True)[0]
Jtf = sp.factor(sp.simplify(solJ[B_])); Jtr = sp.simplify(solJ[A_])
print('junction trace-free:', Jtf); print('junction trace:', sp.factor(Jtr)); print('junction scalar:', sp.factor(jsc1))
ok('junction_tracefree_depends_only_on_zeta', sp.simplify(Jtf.subs(zeta, 0)))
res['junction_tracefree_coefficient_of_zeta'] = str(sp.simplify(sp.diff(Jtf, zeta)))
res['junction_tracefree_forces_zeta_0'] = bool(sp.simplify(sp.diff(Jtf, zeta)) != 0)
# with zeta = 0, N = -2 psi
Jtr0 = js(Jtr.subs(zeta, 0).subs(N, -2*psi).doit())
ok('trace_junction_equals_momentum_constraint', js(Jtr0.subs(psi.diff(y), dpsi)))
res['trace_junction_zeta0'] = str(sp.factor(Jtr0))
jsc0z = js(jsc1.subs(zeta, 0).subs(N, -2*psi).doit())
Mform = chi.diff(y) + s2(ph)/2*chi + 2*ph.diff(y)*psi
r = sp.simplify(jsc0z/js(Mform))
print('scalar junction / (chi\' + s2/2 chi + 2 phi\' psi) =', r)
ok("scalar_junction_equals_chi'+sigma''/2 chi+2phi'psi", sp.simplify(r - 1))
# in (X,Z): B X + 3 lam Z/rho^2 with B = g + sigma''/2
MXZ = js(Mform.subs(chi.diff(y), dchi).subs({psi: ph.diff(y)*Zf, chi: Xf}))
Bexpr = (fpp/ph.diff(y)).subs(bgj) + s2(ph)/2
ok('M=BX+3lamZ/rho^2', js(sp.expand(MXZ - Bexpr*Xf - 3*lam*Zf/rho**2)))
# controls on the junction: wrong-sign sigma'' and missing 2phi'psi must NOT match
res['CONTROL_flip_sigma2_detected'] = bool(sp.simplify(jsc0z - js(chi.diff(y) - s2(ph)/2*chi + 2*ph.diff(y)*psi)) != 0)
res['CONTROL_drop_2phipsi_detected'] = bool(sp.simplify(jsc0z - js(chi.diff(y) + s2(ph)/2*chi)) != 0)

# ------------------------------------------------------------------ tensor sector
gt = metric(0, 0, hten=hT)
Et, boxt, Gt, git = field_eqs(gt, ph, U0)
e12 = bg(sp.expand(lin(Et[2, 3])/(a2*Y)))
num, den = sp.fraction(sp.together(e12))
e12 = sp.simplify(sp.rem(sp.Poly(sp.expand(num), p), sp.Poly(p**2 + 3*p + mu2, p)).as_expr()/den)
rt = sp.simplify(e12/(hT.diff(y, 2) + 4*Hs*hT.diff(y) + mu2*hT/rho**2))
print('tensor ratio', rt)
ok('tensor_eq_proportional(ratio free of hT)', 0 if not rt.has(hT) else 1)
res['tensor_ratio'] = str(rt)
others = [sp.simplify(bg(lin(Et[i, j]))) for i in range(5) for j in range(i, 5) if (i, j) != (2, 3)]
res['tensor_other_components_vanish'] = all(o == 0 for o in others)
nl_t = [1/sp.sqrt(git[0, 0]), 0, 0, 0, 0]
K12 = -sum(Gt[c][2][3]*nl_t[c] for c in range(5))
jt = js(lin(K12 - s0(ph)/6*gt[2, 3])/(a2*Y))
print('tensor junction', jt)
ok('tensor_junction_Neumann', sp.simplify(jt.subs(hT.diff(y), 0)))
res['tensor_junction_expr'] = str(jt)

res_all = {k: v for k, v in res.items() if isinstance(v, bool)}
out = dict(checks=res, n_bool=len(res_all), n_true=sum(res_all.values()),
           note='CONTROL_* entries must be True (= wrong formula detected); all other booleans must be True',
           runtime_s=time.time() - T0, sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('indep_derivation.json').write_text(json.dumps(out, indent=2) + '\n')
print(out['n_true'], '/', out['n_bool'], 'runtime', out['runtime_s'])
