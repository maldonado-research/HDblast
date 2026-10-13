"""Symbolic derivation of the linearized 5D Einstein-scalar equations about a
de Sitter-sliced background, with a Z2 shell and BOTH junction conditions.

Model (units kappa_5^2 = 1):
    S = int sqrt(-g) [ R/2 - (1/2)(d phi)^2 - U(phi) ]  (two Z2 copies)  - int_shell sqrt(-h) sigma(phi)
    field equations  E_AB := R_AB - d_A phi d_B phi - (2/3) U g_AB = 0 ,   Box phi - U_phi = 0
    background       ds^2 = dy^2 + rho(y)^2 gamma ,  gamma = unit dS4 in flat slicing
                     gamma = -dt^2 + e^{2t} dx^2 ,   phi = phi(y)
    shell            Z2, bulk on y < y_b, K_mn = (sigma/6) h_mn ,  n^A d_A phi = -sigma'(phi)/2

dS4 harmonics are realised concretely as Y = e^{p t} (homogeneous in x), for which
    Box_gamma Y = mu2 Y with mu2 = -(p^2 + 3 p),
    nabla_t nabla_t Y = p^2 Y , nabla_i nabla_j Y = -p gamma_ij Y .
The trace (gamma_mn Y) and trace-free (nabla_m nabla_n Y) structures are separable
whenever p not in {0, 1} (the Euclidean l=0, l=1 harmonics); those two cases are
treated separately below.  Every coefficient ODE is checked to depend on p only
through mu2 (reduction modulo p^2 + 3p + mu2).

Output: DERIVATION_RESULTS.json (all checks with pass/fail) and printed equations.
"""
import json, hashlib, time
from pathlib import Path
import sympy as sp

T0 = time.time()
y, t, x1, x2, x3 = sp.symbols('y t x1 x2 x3', real=True)
p, eps, mu2, zeta = sp.symbols('p epsilon mu2 zeta')
X = [y, t, x1, x2, x3]
rho = sp.Function('rho')(y); phi = sp.Function('phi')(y)
N, Bf, psi, E, chi = [sp.Function(n)(y) for n in ('N', 'B', 'psi', 'E', 'chi')]
T_g, L_g = sp.Function('T')(y), sp.Function('L')(y)
U0, U1, U2 = [sp.Function(n)(y) for n in ('U0', 'U1', 'U2')]   # U, U_phi, U_phiphi on background
a2 = sp.exp(2*t)
Y = sp.exp(p*t)
checks = {}

def chk(name, expr):
    ok = sp.simplify(expr) == 0
    checks[name] = bool(ok)
    print(('PASS ' if ok else 'FAIL ') + name, flush=True)
    return ok

# ---------------------------------------------------------------- background
g0 = sp.diag(1, -rho**2, rho**2*a2, rho**2*a2, rho**2*a2)
g0i = g0.inv()
rp = sp.Derivative(rho, y)
# background equations used for substitution
H = rho.diff(y)/rho
bgsub = {
    phi.diff(y, 2): U1 - 4*H*phi.diff(y),
    rho.diff(y, 2): -rho*(phi.diff(y)**2/4 + U0/6),
}
U0_constraint = 6/rho**2 + phi.diff(y)**2/2 - 6*rho.diff(y)**2/rho**2   # from rho'^2 = 1 + rho^2(phi'^2/12 - U/6)

def bgsimp(e):
    """Apply background equations (second derivatives, then the first integral)."""
    for _ in range(3):
        e = e.subs({phi.diff(y, 3): sp.diff(U1 - 4*H*phi.diff(y), y),
                    rho.diff(y, 3): sp.diff(-rho*(phi.diff(y)**2/4 + U0/6), y)})
        e = e.subs(sp.Derivative(U1, y), U2*phi.diff(y)).subs(sp.Derivative(U0, y), U1*phi.diff(y))
        e = e.subs(bgsub)
    e = e.subs(sp.Derivative(U1, y), U2*phi.diff(y)).subs(sp.Derivative(U0, y), U1*phi.diff(y))
    e = e.subs(bgsub)
    e = e.subs(U0, U0_constraint)
    return sp.simplify(e)

def christoffel(g, gi):
    n = 5
    return [[[sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                  for d in range(n))/2 for c in range(n)] for b in range(n)] for a in range(n)]

def ricci_from(G0, G1):
    """Zeroth- and first-order Ricci from Christoffels Gamma = G0 + eps G1."""
    n = 5
    R0 = sp.zeros(n); R1 = sp.zeros(n)
    for b in range(n):
        for c in range(n):
            r0 = sum(sp.diff(G0[a][b][c], X[a]) - sp.diff(G0[a][b][a], X[c]) for a in range(n))
            r0 += sum(G0[a][a][d]*G0[d][b][c] - G0[a][c][d]*G0[d][b][a] for a in range(n) for d in range(n))
            r1 = sum(sp.diff(G1[a][b][c], X[a]) - sp.diff(G1[a][b][a], X[c]) for a in range(n))
            r1 += sum(G1[a][a][d]*G0[d][b][c] + G0[a][a][d]*G1[d][b][c]
                      - G1[a][c][d]*G0[d][b][a] - G0[a][c][d]*G1[d][b][a] for a in range(n) for d in range(n))
            R0[b, c] = r0; R1[b, c] = r1
    return R0, R1

def first_order_gamma(h):
    """Gamma1 = d/deps Gamma(g0 + eps h) at eps=0."""
    n = 5; hi = -g0i*h*g0i
    G = [[[0]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                G[a][b][c] = sum(g0i[a, d]*(sp.diff(h[d, b], X[c]) + sp.diff(h[d, c], X[b]) - sp.diff(h[b, c], X[d]))
                                 + hi[a, d]*(sp.diff(g0[d, b], X[c]) + sp.diff(g0[d, c], X[b]) - sp.diff(g0[b, c], X[d]))
                                 for d in range(n))/2
    return G

G0 = christoffel(g0, g0i)

def linear_equations(h, dphi):
    """Return (E1 matrix, scalar-equation residual) at first order."""
    G1 = first_order_gamma(h)
    R0, R1 = ricci_from(G0, G1)
    dph0 = [sp.diff(phi, xx) for xx in X]
    dph1 = [sp.diff(dphi, xx) for xx in X]
    E1 = sp.zeros(5)
    for A in range(5):
        for Bq in range(5):
            E1[A, Bq] = R1[A, Bq] - (dph0[A]*dph1[Bq] + dph1[A]*dph0[Bq]) - sp.Rational(2, 3)*(U1*dphi*g0[A, Bq] + U0*h[A, Bq])
    hi = -g0i*h*g0i
    box1 = sum(hi[A, Bq]*(sp.diff(phi, X[A], X[Bq]) - sum(G0[C][A][Bq]*dph0[C] for C in range(5)))
               + g0i[A, Bq]*(sp.diff(dphi, X[A], X[Bq]) - sum(G0[C][A][Bq]*dph1[C] + G1[C][A][Bq]*dph0[C] for C in range(5)))
               for A in range(5) for Bq in range(5))
    S1 = box1 - U2*dphi
    return E1, S1, R0

# zeroth order check: background equations are consistent
_, _, R0 = linear_equations(sp.zeros(5), 0)
E0 = sp.zeros(5)
for A in range(5):
    for Bq in range(5):
        E0[A, Bq] = R0[A, Bq] - sp.diff(phi, X[A])*sp.diff(phi, X[Bq]) - sp.Rational(2, 3)*U0*g0[A, Bq]
box0 = sum(g0i[A, Bq]*(sp.diff(phi, X[A], X[Bq]) - sum(G0[C][A][Bq]*sp.diff(phi, X[C]) for C in range(5)))
           for A in range(5) for Bq in range(5))
chk('background_yy', bgsimp(E0[0, 0])); chk('background_tt', bgsimp(E0[1, 1])); chk('background_xx', bgsimp(E0[2, 2]/a2))
chk('background_scalar', bgsimp(box0 - U1))

# ---------------------------------------------------------------- scalar perturbations
def nnY_matrix(Yf, pp):
    M = sp.zeros(4)
    M[0, 0] = pp**2*Yf
    for i in range(1, 4): M[i, i] = -pp*a2*Yf
    return M
gam = sp.diag(-1, a2, a2, a2)
def scalar_h(Nf, Bff, psif, Ef):
    h = sp.zeros(5)
    h[0, 0] = 2*Nf*Y
    h[0, 1] = h[1, 0] = rho**2*Bff*p*Y          # rho^2 B nabla_m Y
    nn = nnY_matrix(Y, p)
    for m in range(4):
        for n_ in range(4):
            h[m+1, n_+1] = rho**2*(2*psif*gam[m, n_]*Y + 2*Ef*nn[m, n_])
    return h

h = scalar_h(N, Bf, psi, E)
E1, S1, _ = linear_equations(h, chi*Y)

def reduce_p(e):
    """Reduce polynomial dependence on p modulo p^2 + 3p + mu2 and report residual p-dependence."""
    e = sp.expand(sp.simplify(e))
    num, den = sp.fraction(sp.together(e))
    num = sp.rem(sp.Poly(sp.expand(num), p), sp.Poly(p**2 + 3*p + mu2, p)).as_expr()
    den = sp.rem(sp.Poly(sp.expand(den), p), sp.Poly(p**2 + 3*p + mu2, p)).as_expr()
    return sp.simplify(num/den)

# components (divide out harmonic dependence)
Eyy = bgsimp(sp.expand(E1[0, 0]/Y))
Eyt = bgsimp(sp.expand(E1[0, 1]/(p*Y)))
Ett = bgsimp(sp.expand(E1[1, 1]/Y))
Exx = bgsimp(sp.expand(E1[2, 2]/(a2*Y)))
Ssc = bgsimp(sp.expand(S1/Y))
chk('offdiag_yx_vanish', sp.simplify(E1[0, 2])); chk('offdiag_tx_vanish', sp.simplify(E1[1, 2])); chk('offdiag_xx_vanish', sp.simplify(E1[2, 3]))
# E_mn = a gamma_mn Y + b nabla_m nabla_n Y :  tt = -a + b p^2 ; xx = a - b p
asym, bsym = sp.symbols('a b')
sol = sp.solve([sp.Eq(-asym + bsym*p**2, Ett), sp.Eq(asym - bsym*p, Exx)], [asym, bsym], dict=True)[0]
Etr, Etf = sp.simplify(sol[asym]), sp.simplify(sol[bsym])
eqs = {}
for name, e in [('yy', Eyy), ('yt', Eyt), ('trace', Etr), ('tracefree', Etf), ('scalar', Ssc)]:
    r = reduce_p(e)
    checks['p_enters_only_via_mu2_' + name] = not r.has(p)
    eqs[name] = sp.simplify(r)
    print(name, ':', eqs[name], flush=True)

# ---------------------------------------------------------------- gauge transformations
xi = [T_g*Y] + [sum(gam.inv()[m, n_]*sp.diff(Y, X[n_+1]) for n_ in range(4))*L_g for m in range(4)]
def lie_g0(xiv):
    Lg = sp.zeros(5)
    for A in range(5):
        for Bq in range(5):
            Lg[A, Bq] = sum(xiv[C]*sp.diff(g0[A, Bq], X[C]) + g0[C, Bq]*sp.diff(xiv[C], X[A]) + g0[A, C]*sp.diff(xiv[C], X[Bq]) for C in range(5))
    return sp.simplify(Lg)
Lg = lie_g0(xi)
gauge_shift = {N: T_g.diff(y), Bf: L_g.diff(y) + T_g/rho**2, psi: H*T_g, E: L_g, chi: phi.diff(y)*T_g}
hg = scalar_h(gauge_shift[N], gauge_shift[Bf], gauge_shift[psi], gauge_shift[E])
checks['gauge_shifts_reproduce_Lie_derivative'] = all(sp.simplify(z) == 0 for z in (Lg - hg))
print('gauge shifts reproduce Lie derivative:', checks['gauge_shifts_reproduce_Lie_derivative'])
# invariance of all linear equations under the gauge shifts (on the background shell of solutions)
fsym = [N, Bf, psi, E, chi]
for name, e in eqs.items():
    eg = e.subs({fsym[i]: gauge_shift[fsym[i]] for i in range(5)}).doit()
    chk('gauge_invariance_' + name, bgsimp(sp.expand(eg)))

# ---------------------------------------------------------------- gauge-invariant variables & master system
# sigma_hat = B - E' shifts by T/rho^2; invariants:
#   Psi = psi - H rho^2 (B-E'),  Xgi = chi - phi' rho^2 (B-E'),  Ngi = N - (rho^2 (B-E'))'
# In longitudinal gauge (B=E=0) they coincide with (psi, chi, N).
long = {Bf: 0, E: 0}
Lyt = sp.simplify(eqs['yt'].subs(long).doit()); Ltf = sp.simplify(eqs['tracefree'].subs(long).doit())
print('longitudinal yt :', sp.factor(Lyt)); print('longitudinal tracefree :', sp.factor(Ltf))
Nsol = sp.solve(Ltf, N)[0]
chk('tracefree_gives_N_eq_minus2psi', Nsol + 2*psi)
# Hamiltonian (Gauss) constraint: C = E_yy - (1/2) g0^{CD} E_CD  -> (1/2)[E_yy - rho^-2 (4 a + mu2 b)]
Cgauss = sp.simplify((eqs['yy'] - (4*eqs['trace'] + mu2*eqs['tracefree'])/rho**2)/2)
Lc = sp.simplify(Cgauss.subs(long).doit().subs(N, -2*psi).doit())
checks['gauss_constraint_has_no_second_derivatives'] = not (Lc.has(psi.diff(y, 2)) or Lc.has(chi.diff(y, 2)))
Lyt2 = sp.simplify(Lyt.subs(N, -2*psi).doit())
dpsi = sp.solve(Lyt2, psi.diff(y))[0]
dchi = sp.solve(Lc.subs(psi.diff(y), dpsi), chi.diff(y))[0]
print('psi\' =', sp.simplify(dpsi)); print('chi\' =', sp.simplify(dchi))
# compare with the first-order system  Psi' = -2 H Psi - phi' X/3 ;  X' = (phi''/phi') X + (3 lam/(rho^2 phi') - 2 phi') Psi
lam = mu2 + 4
phipp = U1 - 4*H*phi.diff(y)
chk('momentum_constraint_form', dpsi - (-2*H*psi - phi.diff(y)*chi/3))
chk('gauss_constraint_form', bgsimp(dchi - (phipp/phi.diff(y)*chi + (3*lam/(rho**2*phi.diff(y)) - 2*phi.diff(y))*psi)))
# Bianchi-type consistency: remaining second-order equations hold identically on the first-order system
first = {psi.diff(y): dpsi, chi.diff(y): dchi}
second = {psi.diff(y, 2): sp.diff(dpsi, y), chi.diff(y, 2): sp.diff(dchi, y)}
for name in ('trace', 'scalar', 'yy'):
    e = eqs[name].subs(long).doit().subs(N, -2*psi).doit()
    e = e.subs(second).subs(first).doit()
    e = e.subs(first)
    chk('implied_by_first_order_system_' + name, bgsimp(sp.expand(e)))

# Z = Psi/phi' formulation (finite when phi' -> 0)
Zf = sp.Function('Z')(y); Xf = sp.Function('X')(y)
g_ = sp.Symbol('g')  # g = phi''/phi'
Zp = sp.simplify((dpsi.subs(psi, phi.diff(y)*Zf) - phi.diff(y, 1).diff(y)*Zf)/phi.diff(y))
Zp = sp.simplify(bgsimp(Zp.subs(chi, Xf)))
Xp = sp.simplify(bgsimp(dchi.subs({psi: phi.diff(y)*Zf, chi: Xf})))
print("Z' =", Zp); print("X' =", Xp)
chk('Z_equation', bgsimp(Zp - (-(2*H + phipp/phi.diff(y))*Zf - Xf/3)))
chk('X_equation', bgsimp(Xp - (phipp/phi.diff(y)*Xf + (3*lam/rho**2 - 2*phi.diff(y)**2)*Zf)))

# identities used for the analytic bounds (real mu2 solutions; two-solution form)
X1, Z1, X2, Z2 = [sp.Function(n)(y) for n in ('X1', 'Z1', 'X2', 'Z2')]
l1, l2 = sp.symbols('lambda1 lambda2')
def zx_rhs(Xv, Zv, l):
    return (phipp/phi.diff(y)*Xv + (3*l/rho**2 - 2*phi.diff(y)**2)*Zv, -(2*H + phipp/phi.diff(y))*Zv - Xv/3)
d1 = zx_rhs(X1, Z1, l1); d2 = zx_rhs(X2, Z2, l2)
sub12 = {X1.diff(y): d1[0], Z1.diff(y): d1[1], X2.diff(y): d2[0], Z2.diff(y): d2[1]}
# W12 = rho^2 (Z2 X1 - X2 Z1)   (i.e. (rho^2/phi')(Psi2 X1 - X2 Psi1))
W12 = rho**2*(Z2*X1 - X2*Z1)
chk('wronskian_identity_[rho^2(Z2X1-X2Z1)]\'=3(l1-l2)Z1Z2', bgsimp(sp.expand(W12.diff(y).subs(sub12) - 3*(l1 - l2)*Z1*Z2)))
Q11 = rho**2*Z1*X1
chk('energy_identity_[rho^2 Z X]\'=3 l Z^2 - rho^2 X^2/3 - 2 rho^2 phi^2 Z^2',
    bgsimp(sp.expand(Q11.diff(y).subs(sub12) - (3*l1*Z1**2 - rho**2*X1**2/3 - 2*rho**2*phi.diff(y)**2*Z1**2))))

# ---------------------------------------------------------------- junction conditions on the displaced shell
yb = sp.Symbol('y_b')
dl, cc = sp.symbols('delta c')
Wp = lambda f: 1 - f + f**3/3
sig = lambda f: 2*Wp(f) + dl*(1 + cc*f)
F = y - eps*zeta*Y                     # shell: F = y_b  (y_b constant)
gfull = g0 + eps*h
gfi = g0i - eps*g0i*h*g0i
dF = [sp.diff(F, xx) for xx in X]
norm2 = sum(gfi[A, Bq]*dF[A]*dF[Bq] for A in range(5) for Bq in range(5))
nlow = [dF[A]/sp.sqrt(norm2) for A in range(5)]
Gfull = [[[G0[a][b][c] + eps*g for c, g in enumerate(row)] for b, row in enumerate(plane)] for a, plane in enumerate(first_order_gamma(h))]
# tangent vectors e_m = d_m + eps zeta d_m Y d_y
emat = []
for m in range(4):
    v = [0]*5; v[m+1] = 1; v[0] = eps*zeta*sp.diff(Y, X[m+1]); emat.append(v)
def Kmn(m, n_):
    return sum(emat[m][A]*emat[n_][Bq]*(sp.diff(nlow[Bq], X[A]) - sum(Gfull[C][A][Bq]*nlow[C] for C in range(5)))
               for A in range(5) for Bq in range(5))
def hmn(m, n_):
    return sum(emat[m][A]*emat[n_][Bq]*gfull[A, Bq] for A in range(5) for Bq in range(5))
phifull = phi + eps*chi*Y
nup = [sum(gfi[A, Bq]*nlow[Bq] for Bq in range(5)) for A in range(5)]
Jphi = sum(nup[A]*sp.diff(phifull, X[A]) for A in range(5))
fsig = sp.Symbol('f')
sigp = sp.diff(sig(fsig), fsig)
Jtt = Kmn(0, 0) - sig(phifull)/6*hmn(0, 0)
Jxx = Kmn(1, 1) - sig(phifull)/6*hmn(1, 1)
Jsc = Jphi + sigp.subs(fsig, phifull)/2

def at_shell(expr):
    """First-order part of expr evaluated on the displaced shell y = y_b + eps zeta Y."""
    e0 = expr.subs(eps, 0)
    e1 = sp.diff(expr, eps).subs(eps, 0)
    tot = e1 + zeta*Y*sp.diff(e0, y)
    return tot, e0
bgjunc = {rho.diff(y): rho*sig(phi)/6, phi.diff(y): -sigp.subs(fsig, phi)/2}
def shell_simp(e):
    e = bgsimp(sp.expand(e))
    e = e.subs(bgjunc)
    return sp.simplify(e)
J = {}
for name, e in [('tt', Jtt), ('xx', Jxx), ('scalar', Jsc)]:
    e1, e0 = at_shell(e)
    J[name] = e1; J[name + '_bg'] = e0
# background junctions hold identically after substituting them
chk('background_metric_junction', shell_simp(J['xx_bg']/a2))
chk('background_scalar_junction', shell_simp(J['scalar_bg']))
jtt = shell_simp(J['tt']/Y); jxx = shell_simp(J['xx']/(a2*Y)); jsc = shell_simp(J['scalar']/Y)
solJ = sp.solve([sp.Eq(-asym + bsym*p**2, jtt), sp.Eq(asym - bsym*p, jxx)], [asym, bsym], dict=True)[0]
Jtr, Jtf, Jsc_r = reduce_p(solJ[asym]), reduce_p(solJ[bsym]), reduce_p(jsc)
checks['junction_p_only_via_mu2'] = not (Jtr.has(p) or Jtf.has(p) or Jsc_r.has(p))
print('junction trace     :', sp.factor(Jtr)); print('junction trace-free:', sp.factor(Jtf)); print('junction scalar    :', sp.factor(Jsc_r))
# gauge invariance of the junctions with zeta -> zeta - T(y_b)
gsub = {fsym[i]: gauge_shift[fsym[i]] for i in range(5)}
for name, e in [('trace', Jtr), ('tracefree', Jtf), ('scalar', Jsc_r)]:
    eg = e.subs(gsub).doit().subs(zeta, -T_g)   # pure gauge perturbation with zeta = -T : must vanish
    chk('junction_gauge_invariance_' + name, shell_simp(sp.expand(eg)))
# longitudinal gauge, N = -2 psi
Jl = {k: shell_simp(v.subs(long).doit().subs(N, -2*psi).doit()) for k, v in [('trace', Jtr), ('tracefree', Jtf), ('scalar', Jsc_r)]}
print('longitudinal junctions:', Jl)
chk('tracefree_junction_is_pure_zeta', sp.simplify(Jl['tracefree'].subs(zeta, 0)))
checks['tracefree_junction_forces_zeta_0'] = bool(sp.simplify(sp.diff(Jl['tracefree'], zeta)) != 0)
# trace junction with zeta=0 reduces to the bulk momentum constraint (at the shell)
jt0 = Jl['trace'].subs(zeta, 0)
chk('trace_junction_equals_momentum_constraint', shell_simp(jt0.subs(psi.diff(y), dpsi)))
# scalar junction with zeta=0 -> M = X' + (sigma''/2) X + 2 phi' Psi  (proportional)
js0 = sp.simplify(Jl['scalar'].subs(zeta, 0))
sig2 = sp.diff(sig(fsig), fsig, 2).subs(fsig, phi)
Mform = chi.diff(y) + sig2/2*chi + 2*phi.diff(y)*psi
ratio = sp.simplify(js0/shell_simp(Mform))
print('scalar junction / M =', ratio)
checks['scalar_junction_proportional_to_M'] = bool(sp.simplify(sp.diff(ratio, y)) == 0 and not ratio.has(chi) and not ratio.has(psi))
# M in (X, Z) variables: using X' from the system,  M = (phi''/phi' + sigma''/2) X + 3 lam Z / rho^2
Mz = shell_simp(Mform.subs(chi.diff(y), dchi)).subs({chi: Xf, psi: phi.diff(y)*Zf})
chk('M_in_XZ_form', shell_simp(sp.expand(Mz - ((phipp/phi.diff(y)).subs(bgjunc) + sig2/2)*Xf - 3*lam*Zf/rho**2)))

# ---------------------------------------------------------------- special harmonic p = 1 (mu2 = -4, l=1)
# nabla nabla Y = -gamma Y: E is not independent (only psi - E).  Test: pure gauge bulk (regular) plus a
# relative shell displacement Delta: junction residuals must be (0, B phi' Delta) with B = phi''/phi' + sigma''/2
Delta = sp.Symbol('Delta')
Jtr1 = sp.simplify(solJ[asym]).subs(p, 1)  # not separable at p=1; use components directly
e_tt = shell_simp(jtt.subs(p, 1).subs({N: 0, Bf: 0, psi: 0, E: 0, chi: 0}).subs(zeta, Delta))
e_xx = shell_simp(jxx.subs(p, 1).subs({N: 0, Bf: 0, psi: 0, E: 0, chi: 0}).subs(zeta, Delta))
e_sc = shell_simp(jsc.subs(p, 1).subs({N: 0, Bf: 0, psi: 0, E: 0, chi: 0}).subs(zeta, Delta))
chk('l1_displacement_metric_junction_tt_vanishes', e_tt); chk('l1_displacement_metric_junction_xx_vanishes', e_xx)
chk('l1_displacement_scalar_junction_equals_phi2+sig2phi1/2', shell_simp(e_sc - Delta*(phipp + sig2/2*phi.diff(y)).subs(bgjunc)))
# at p=1 the regular bulk solution space is pure gauge: verify the l=1 bulk equations in GN gauge (N=B=0)
# for fields (Pt = psi - E, chi) reduce to a first-order 2x2 system whose solutions are the 2-parameter
# gauge family  Pt = (H + F) e0 - L0 ,  chi = phi' e0  with F' = 1/rho^2.
Pt = sp.Function('Pt')(y); Ff = sp.Function('F')(y); e0s, L0s = sp.symbols('e0 L0')
eq_p1 = {}
for name, e in [('yy', Eyy), ('yt', Eyt), ('tt', Ett), ('xx', Exx), ('scalar', Ssc)]:
    ee = e.subs(p, 1).subs({N: 0, Bf: 0}).doit().subs(psi, Pt + E).doit()
    eq_p1[name] = bgsimp(sp.expand(ee))
checks['l1_E_drops_out'] = all(not v.has(E) for v in eq_p1.values())
gfam = {Pt: (H + Ff)*e0s - L0s, chi: phi.diff(y)*e0s}
for name, e in eq_p1.items():
    ee = e.subs(gfam).doit().subs({Ff.diff(y, 2): sp.diff(1/rho**2, y)}).subs({Ff.diff(y): 1/rho**2})
    chk('l1_gauge_family_solves_' + name, bgsimp(sp.expand(ee)))
Cg1 = bgsimp(sp.expand((eq_p1['yy'] - (-eq_p1['tt'] + 3*eq_p1['xx'])/rho**2)/2))
checks['l1_gauss_constraint_first_order'] = not (Cg1.has(Pt.diff(y, 2)) or Cg1.has(chi.diff(y, 2)))
sol_p1 = sp.solve([eq_p1['yt'], Cg1], [Pt.diff(y), chi.diff(y)], dict=True)
fam_first = [bgsimp(sp.expand((Pt.diff(y) - sol_p1[0][Pt.diff(y)]).subs(gfam).doit().subs({Ff.diff(y): 1/rho**2}))),
             bgsimp(sp.expand((chi.diff(y) - sol_p1[0][chi.diff(y)]).subs(gfam).doit().subs({Ff.diff(y): 1/rho**2})))]
chk('l1_gauge_family_spans_first_order_system', fam_first[0]**2 + fam_first[1]**2)
checks['l1_bulk_is_first_order_2x2'] = bool(len(sol_p1) == 1)
print('l=1 first-order system:', sol_p1)

# ---------------------------------------------------------------- tensor sector
hT = sp.Function('hT')(y)
ht = sp.zeros(5); ht[2, 3] = ht[3, 2] = rho**2*a2*hT*Y     # TT, spatially homogeneous: h_x1x2
E1t, S1t, _ = linear_equations(ht, 0)
comp = bgsimp(sp.expand(E1t[2, 3]/(a2*Y)))
others = [sp.simplify(E1t[i, j]) for i in range(5) for j in range(5) if (i, j) not in ((2, 3), (3, 2))]
checks['tensor_other_components_vanish'] = all(o == 0 for o in others) and sp.simplify(S1t) == 0
tensor_eq = reduce_p(comp)
print('tensor:', sp.factor(tensor_eq))
tensor_ratio = sp.simplify(tensor_eq/(hT.diff(y, 2) + 4*H*hT.diff(y) + mu2*hT/rho**2))
chk('tensor_equation_is_h2+4Hh1+mu2h/rho^2', sp.expand(tensor_eq*(-2/rho**2) - (hT.diff(y, 2) + 4*H*hT.diff(y) + mu2*hT/rho**2)))
print('tensor ratio:', tensor_ratio)
# tensor junction: K_12 - sigma/6 h_12  -> proportional to hT'
gfull_t = g0 + eps*ht; gfi_t = g0i - eps*g0i*ht*g0i
Gt = [[[G0[a][b][c] + eps*g for c, g in enumerate(row)] for b, row in enumerate(plane)] for a, plane in enumerate(first_order_gamma(ht))]
nl = [1/sp.sqrt(gfi_t[0, 0]), 0, 0, 0, 0]
K12 = -sum(Gt[C][2][3]*nl[C] for C in range(5))
jt = sp.diff(K12 - sig(phi)/6*gfull_t[2, 3], eps).subs(eps, 0)
jt = shell_simp(jt/(a2*Y))
print('tensor junction:', jt)
chk('tensor_junction_is_Neumann', sp.simplify(jt.subs(hT.diff(y), 0)))

# ---------------------------------------------------------------- deliberate wrong-formula controls (must FAIL)
controls = {}
wrong_shift = dict(gauge_shift); wrong_shift[psi] = -H*T_g
e_bad = eqs['yt'].subs({fsym[i]: wrong_shift[fsym[i]] for i in range(5)}).doit()
controls['wrong_gauge_shift_psi=-HT_breaks_yt_invariance'] = bool(bgsimp(sp.expand(e_bad)) != 0)
e_bad = Jsc_r.subs(gsub).doit().subs(zeta, +T_g)
controls['wrong_zeta_shift_sign_breaks_scalar_junction_invariance'] = bool(shell_simp(sp.expand(e_bad)) != 0)
M_bad = chi.diff(y) + sig2/2*chi          # omit the 2 phi' psi gravitational term
controls['scalar_junction_without_2phiPsi_term_is_not_the_junction'] = bool(sp.simplify(js0 - shell_simp(M_bad)) != 0)
controls['wrong_Z_equation_sign_detected'] = bool(bgsimp(Zp - (-(2*H + phipp/phi.diff(y))*Zf + Xf/3)) != 0)
for k, v in controls.items():
    print(('CONTROL-OK ' if v else 'CONTROL-FAILED ') + k)
out = dict(status='SYMBOLIC_DERIVATION', checks=checks, all_pass=all(checks.values()), wrong_formula_controls_detected=controls, all_controls_detected=all(controls.values()),
           n_checks=len(checks), n_pass=sum(checks.values()),
           equations=dict(
               linear_yy=str(eqs['yy']), linear_yt=str(eqs['yt']), linear_trace=str(eqs['trace']),
               linear_tracefree=str(eqs['tracefree']), linear_scalar=str(eqs['scalar']),
               gauge_shifts={str(k): str(v) for k, v in gauge_shift.items()},
               first_order_system="Z' = -(2H + phi''/phi') Z - X/3 ;  X' = (phi''/phi') X + (3(mu2+4)/rho^2 - 2 phi'^2) Z ;  Psi = phi' Z",
               shell_condition="M = (phi''/phi' + sigma''/2) X + 3(mu2+4) Z/rho^2 = 0 at y_b  (equivalently X' + sigma''/2 X + 2 phi' Psi = 0)",
               junction_trace_longitudinal=str(Jl['trace']), junction_tracefree_longitudinal=str(Jl['tracefree']),
               junction_scalar_longitudinal=str(Jl['scalar']),
               tensor="hT'' + 4 H hT' + mu2 hT/rho^2 = 0 , hT'(y_b) = 0",
               tensor_ratio=str(tensor_ratio)),
           runtime_s=time.time() - T0,
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('DERIVATION_RESULTS.json').write_text(json.dumps(out, indent=2) + '\n')
print('ALL PASS' if out['all_pass'] else 'SOME CHECKS FAILED', out['n_pass'], '/', out['n_checks'], 'time', out['runtime_s'])
