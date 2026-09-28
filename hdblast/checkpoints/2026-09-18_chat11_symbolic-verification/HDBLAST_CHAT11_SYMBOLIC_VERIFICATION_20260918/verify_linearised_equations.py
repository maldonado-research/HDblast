#!/usr/bin/env python3
"""HDBLAST Chat 11 - SYMBOLIC verification (SymPy, exact) of the linearised scalar-sector equations used by the
Chat 9 tachyon certificate and the Chat 10 stability certificate.

Set-up.  Coordinates (y, tau, x1, x2, x3).  Unit de Sitter_4 in flat slicing: gamma = -dtau^2 + e^{2 tau} dx^2  (R[gamma]_{mu nu} = 3 gamma).
  metric   ds^2 = (1 + 2 eps xi(y) Y) dy^2 + rho(y)^2 (1 + 2 eps psi(y) Y) gamma_{mu nu} dx^mu dx^nu      (longitudinal gauge)
  scalar   phi  = phi0(y) + eps chi(y) Y ,     Y = T(tau) E(x1),   E'' = -k^2 E,   T'' = -3 T' - (k^2 e^{-2 tau} + mu2) T   <=>  Box_gamma Y = mu2 Y
  field equations  EQ_AB := R_AB - d_A phi d_B phi - (2/3) U(phi) g_AB = 0,     SC := Box_g phi - U'(phi) = 0,   U arbitrary.
The script computes EQ_AB and SC exactly to first order in eps and proves, as polynomial identities:
  (0) the background equations;   (1) EQ_{tau x1}  is proportional to (xi + 2 psi);   (2) EQ_{y mu}  <=>  psi' - H xi = -phi0' chi/3  (Codazzi);
  (3) with xi = -2 psi, chi = -3(psi' + 2 H psi)/phi0' and the master equation
         psi'' + 2(H - phi0''/phi0') psi' + [ -(4/3) phi0'^2 - 4 H phi0''/phi0' + (2 + mu2)/rho^2 ] psi = 0 ,   H = rho'/rho,
      ALL fifteen components EQ_AB and SC vanish at first order;
  (4) the shell-condition identities  chi' + 2 phi0' psi + (s2/2) chi  ==  B chi + 3 lam psi/(rho^2 phi0')  ==  (3 psi/phi0') [ lam/rho^2 - B (psi'/psi + 2H) ]
      on solutions, with B = phi0''/phi0' + s2/2, lam = mu2 + 4  (s2 = sigma_t''(phi_b)).
Run with the isolated environment:   ../../.hdblast_venv/bin/python verify_linearised_equations.py
"""
import json, sys, time
import sympy as sp
t0 = time.time()
y, tau, x1, x2, x3, eps, k, mu2 = sp.symbols('y tau x1 x2 x3 epsilon k mu2')
X = [y, tau, x1, x2, x3]
rho, ph0, psi, xi, chi = [sp.Function(n)(y) for n in ('rho', 'phi0', 'psi', 'xi', 'chi')]
T = sp.Function('T')(tau); E = sp.Function('E')(x1); Y = T*E
U0, U1, U2 = sp.symbols('U0 U1 U2')            # U(phi0), U'(phi0), U''(phi0) as values; dU0/dy = U1 phi0', dU1/dy = U2 phi0'

def lin(e):                                     # truncate to first order in eps
    e = sp.expand(e); return sp.expand(e.coeff(eps, 0) + eps*e.coeff(eps, 1))
gam = [-1, sp.exp(2*tau), sp.exp(2*tau), sp.exp(2*tau)]
g = [1 + 2*eps*xi*Y] + [rho**2*(1 + 2*eps*psi*Y)*gm for gm in gam]                      # diagonal metric
ginv = [1 - 2*eps*xi*Y] + [(1 - 2*eps*psi*Y)/(rho**2*gm) for gm in gam]                 # inverse to O(eps)
n = 5
def dg(a, c): return sp.diff(g[a], X[c])
Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(b, n):
            s = 0
            if a == c: s += dg(a, b)
            if a == b: s += dg(a, c)
            if b == c: s -= dg(b, a)
            v = lin(sp.Rational(1, 2)*ginv[a]*s) if s != 0 else 0
            Gam[a][b][c] = Gam[a][c][b] = v
def ricci(b, c):
    r = 0
    for a in range(n):
        r += sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
        for d in range(n):
            r += Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a]
    return lin(r)
phi = ph0 + eps*chi*Y
dphi = [sp.diff(phi, v) for v in X]
Uphi = U0 + eps*U1*chi*Y; Uphi1 = U1 + eps*U2*chi*Y
EQ = {}
for b in range(n):
    for c in range(b, n):
        EQ[(b, c)] = lin(ricci(b, c) - dphi[b]*dphi[c] - sp.Rational(2, 3)*Uphi*(g[b] if b == c else 0))
sqrtg_lin = lin(rho**4*sp.exp(3*tau)*(1 + eps*xi*Y + 4*eps*psi*Y))                       # sqrt(-g) to O(eps)
box = sum(sp.diff(lin(sqrtg_lin*ginv[a]*dphi[a]), X[a]) for a in range(n))
SC = lin(lin(box*(1 - eps*xi*Y - 4*eps*psi*Y)/(rho**4*sp.exp(3*tau))) - Uphi1)
print("linearised Ricci tensor and scalar equation built in %.1f s" % (time.time() - t0), flush=True)

# ---- substitution rules
rp, pp = sp.symbols('rho1 phi1')                 # rho', phi0' as algebraic symbols after differentiation
H = rp/rho
bg_rho2 = rp**2/rho - 1/rho - rho*pp**2/3        # rho''  from (rho'/rho)' = -1/rho^2 - phi0'^2/3
bg_ph2 = U1 - 4*H*pp                             # phi0''
U0_con = 6*(1 - rp**2)/rho**2 + pp**2/2          # Hamiltonian constraint  rho'^2 = 1 + rho^2(phi0'^2/12 - U0/6)
def bg_reduce(e, extra=None):
    """replace derivatives by the background ODEs (highest order first), then symbols"""
    rules3 = {sp.diff(rho, y, 3): None, sp.diff(ph0, y, 3): None}
    # third derivatives: differentiate the second-derivative rules
    r2 = (sp.diff(rho, y)**2/rho - 1/rho - rho*sp.diff(ph0, y)**2/3)
    p2 = (U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))
    r3 = sp.diff(r2, y); p3 = sp.diff(p2, y) + U2*sp.diff(ph0, y)      # dU1/dy = U2 phi0'
    e = e.subs({sp.diff(rho, y, 3): r3, sp.diff(ph0, y, 3): p3})
    e = e.subs({sp.diff(rho, y, 2): r2, sp.diff(ph0, y, 2): p2})
    e = e.subs({sp.diff(rho, y, 2): r2, sp.diff(ph0, y, 2): p2})
    e = e.subs({sp.diff(rho, y): rp, sp.diff(ph0, y): pp})
    e = e.subs(U0, U0_con)
    return e
T1 = sp.symbols('T1'); E1 = sp.symbols('E1'); Ts, Es = sp.symbols('T0 E0')
def harm_reduce(e):
    e = e.subs(sp.diff(T, tau, 2), -3*sp.diff(T, tau) - (k**2*sp.exp(-2*tau) + mu2)*T)
    e = e.subs(sp.diff(E, x1, 2), -k**2*E)
    return e.subs({sp.diff(T, tau): T1, sp.diff(E, x1): E1}).subs({T: Ts, E: Es})
def zero(e): return sp.simplify(sp.together(sp.expand(e))) == 0

results = {}
def record(name, ok, note=""):
    results[name] = dict(passed=bool(ok), note=note); print("[%s] %s %s" % ("PASS" if ok else "FAIL", name, note), flush=True)

# (0) background
bgok = all(zero(bg_reduce(EQ[key].coeff(eps, 0))) for key in EQ) and zero(bg_reduce(SC.coeff(eps, 0)))
record("0_background_equations_all_components", bgok)

first = {key: harm_reduce(bg_reduce(EQ[key].coeff(eps, 1))) for key in EQ}
SC1 = harm_reduce(bg_reduce(SC.coeff(eps, 1)))

# (1) traceless constraint: EQ_{tau x1} proportional to (xi + 2 psi)
e = sp.factor(sp.simplify(first[(1, 2)]))
ratio = sp.simplify(first[(1, 2)]/(xi + 2*psi))
record("1_offdiagonal_tau_x1_proportional_to_xi_plus_2psi", (not ratio.has(xi)) and (not ratio.has(psi)) and ratio != 0, "EQ_{tau x1} = (%s) * (xi + 2 psi)" % sp.simplify(ratio))

# (2) Codazzi: EQ_{y tau} and EQ_{y x1}
cod = -3*(sp.diff(psi, y) - H*xi) - pp*chi
r_yt = sp.simplify(first[(0, 1)]/cod); r_yx = sp.simplify(first[(0, 2)]/cod)
okc = all((not r.has(psi)) and (not r.has(xi)) and (not r.has(chi)) and r != 0 for r in (r_yt, r_yx))
record("2_momentum_constraint_is_Codazzi_relation", okc, "EQ_{y tau} = (%s)*C, EQ_{y x1} = (%s)*C, C = -3(psi' - H xi) - phi0' chi" % (r_yt, r_yx))

# (3) everything vanishes on the master system
g_ = bg_ph2/pp
master_psi2 = -2*(H - g_)*sp.diff(psi, y) - (-sp.Rational(4, 3)*pp**2 - 4*H*g_ + (2 + mu2)/rho**2)*psi
chi_expr = -3*(sp.diff(psi, y) + 2*sp.diff(rho, y)/rho*psi)/sp.diff(ph0, y)
def on_shell(e):
    e = e.subs(xi, -2*psi)
    e = e.subs({sp.diff(chi, y, 2): sp.diff(chi_expr, y, 2), sp.diff(chi, y): sp.diff(chi_expr, y), chi: chi_expr}).doit()
    e = bg_reduce(e)
    # psi''' and psi'' from the master equation (master_psi2 is written with rp, pp symbols: rebuild its y-derivative via functions)
    m2f = (-2*(sp.diff(rho, y)/rho - (U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))/sp.diff(ph0, y))*sp.diff(psi, y)
           - (-sp.Rational(4, 3)*sp.diff(ph0, y)**2 - 4*sp.diff(rho, y)/rho*(U1 - 4*sp.diff(rho, y)/rho*sp.diff(ph0, y))/sp.diff(ph0, y) + (2 + mu2)/rho**2)*psi)
    m3f = sp.diff(m2f, y) + sp.diff(m2f, U1)*U2*sp.diff(ph0, y)
    e = e.subs(sp.diff(psi, y, 3), m3f).subs(sp.diff(psi, y, 2), m2f).subs(sp.diff(psi, y, 2), m2f)
    return bg_reduce(e)
names = {0: 'y', 1: 'tau', 2: 'x1', 3: 'x2', 4: 'x3'}
allok = True
for key in sorted(first):
    ok = zero(on_shell(first[key])); allok &= ok
    record("3_EQ_%s_%s_vanishes_on_master_system" % (names[key[0]], names[key[1]]), ok)
ok = zero(on_shell(SC1)); allok &= ok
record("3_scalar_field_equation_vanishes_on_master_system", ok)

# (4) shell-condition identities
s2, lam = sp.symbols('s2 lam')
Bq = bg_ph2/pp + s2/2
chi_s = -3*(sp.diff(psi, y) + 2*H*psi)/pp
lhs = on_shell(sp.diff(chi, y) + 2*sp.diff(ph0, y)*psi + (s2/2)*chi)
mid = on_shell(chi)*Bq + 3*(mu2 + 4)*psi/(rho**2*pp)
rhs = (3*psi/pp)*((mu2 + 4)/rho**2 - Bq*(-(sp.diff(psi, y)/psi + 2*H))*(-1))
record("4a_scalar_junction_equals_B_chi_plus_3lam_psi_form", zero(lhs - mid))
record("4b_mismatch_equals_(3psi/phi')[lam w + B R],_R=-(psi'/psi+2H)", zero(mid - (3*psi/pp)*((mu2 + 4)/rho**2 + Bq*(-(sp.diff(psi, y)/psi + 2*H)))))

# (5) NECESSITY: with the two constraints imposed (but NOT the master equation) every remaining equation is a multiple of the
#     master-equation residual  Mres := psi'' + 2(H - phi0''/phi0') psi' + [ -(4/3)phi0'^2 - 4H phi0''/phi0' + (2+mu2)/rho^2 ] psi
def on_constraints(e):
    e = e.subs(xi, -2*psi)
    e = e.subs({sp.diff(chi, y, 2): sp.diff(chi_expr, y, 2), sp.diff(chi, y): sp.diff(chi_expr, y), chi: chi_expr}).doit()
    return bg_reduce(e)
Mres = sp.diff(psi, y, 2) - master_psi2
ryy = on_constraints(first[(0, 0)])
ratio_yy = sp.simplify(ryy/Mres)
record("5_EQ_yy_on_constraints_is_a_nonzero_multiple_of_the_master_residual", (not ratio_yy.has(psi)) and ratio_yy != 0, "EQ_yy = (%s) * Mres" % ratio_yy)

ok_all = all(v["passed"] for v in results.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok_all, runtime_s=round(time.time() - t0, 1), checks=results), open("SYMBOLIC_VERIFICATION_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok_all else "SOME CHECKS FAILED", "(%.1f s)" % (time.time() - t0)); sys.exit(0 if ok_all else 1)
