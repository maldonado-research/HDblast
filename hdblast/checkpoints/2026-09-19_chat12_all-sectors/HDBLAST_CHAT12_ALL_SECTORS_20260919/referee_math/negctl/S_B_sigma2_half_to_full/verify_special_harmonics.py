#!/usr/bin/env python3
"""Chat 12 - the SPECIAL scalar harmonics (mu2 = -4:  nabla_mu nabla_nu Y = -gamma_{mu nu} Y, the five 'l = 1' functions), exact symbolic checks.
For these the trace-free (mu nu) equation is empty, so xi = -2 psi is not forced and the shell may bend; they were excluded in Chats 9-11.
 S1  Y = e^tau and Y = e^tau x1 are special harmonics (Hessian = -gamma Y, Box Y = -4 Y).
 S2  every bulk perturbation in this sector is a diffeomorphism.  With v^y = a(y)Y, v^mu = b(y) nabla^mu Y the background changes by
        xi_g = a',  psi_g = H a - b,  chi_g = phi0' a,  (dy dx^mu component) = a + rho^2 b'.
     Imposing (dy dx) = 0 and xi_g + 2 psi_g = 0 gives  a = -rho^2 b',   rho^2 b'' + 4 rho rho' b' + 2 b = 0   (a two-parameter family), and then
     psi_g solves the master equation at mu2 = -4 and chi_g obeys the Codazzi relation: the two-dimensional solution space of the
     master equation at mu2 = -4 consists of pure-gauge configurations (the map b -> psi_g is injective unless phi0' = 0).
     Near the cone b ~ 1/y, a -> const is the rigid translation of the apex (a regular diffeomorphism) and gives the regular branch psi ~ y^3.
 S3  what remains is the position of the shell relative to the bulk solution.  For the surface y = y_b + eps zeta Y in the UNPERTURBED bulk:
        K^mu_nu - sigma_t(phi)/6 delta^mu_nu = O(eps^2)            (trace Israel condition holds identically, using H' = -1/rho^2 - phi0'^2/3, sigma_t' = -2 phi0'),
        n.d phi + sigma_t'(phi)/2 = eps zeta Y phi0' B + O(eps^2),   B = phi0''/phi0' + sigma_t''/2.
     Hence B != 0 forces zeta = 0: no physical mode in the special sector.  (B = +5.3187e-4 is interval-certified for S_8/5; B = -2.68e-4 on the registered shell.)
Run:  <venv python> verify_special_harmonics.py"""
import json, sys
import sympy as sp
from lin_gr import *
res = {}
def record(n, ok, note=""): res[n] = dict(passed=bool(ok), note=note); print("[%s] %s %s" % ("PASS" if ok else "FAIL", n, note), flush=True)
xs = X[1:]; gi4 = gam.inv()
def hess(Yf): return sp.Matrix(4, 4, lambda m, n_: sp.diff(Yf, xs[m], xs[n_]) - sum(chr4(l, m, n_)*sp.diff(Yf, xs[l]) for l in range(4)))
specials = {"e^tau": sp.exp(tau), "e^tau*x1": sp.exp(tau)*x1}
ok = True
for nm, Yf in specials.items():
    Hs = hess(Yf); ok &= all(sp.simplify(Hs[m, n_] + gam[m, n_]*Yf) == 0 for m in range(4) for n_ in range(4))
    ok &= sp.simplify(sum(gi4[m, m]*Hs[m, m] for m in range(4)) + 4*Yf) == 0
record("S1_two_explicit_special_harmonics_Hess=-gamma_Y_and_Box_Y=-4Y", ok)

# S2: gauge family and the master equation at mu2 = -4
a, b = sp.Function('a')(y), sp.Function('b')(y)
g0 = background_metric(); ok = True
for nm, Yf in specials.items():
    v = [a*Yf] + [b*sum(gi4[m, n_]*sp.diff(Yf, xs[n_]) for n_ in range(4)) for m in range(4)]
    Lie = sp.Matrix(5, 5, lambda A, B_: sum(v[C]*sp.diff(g0[A, B_], X[C]) + g0[C, B_]*sp.diff(v[C], X[A]) + g0[A, C]*sp.diff(v[C], X[B_]) for C in range(5)))
    ok &= sp.simplify(Lie[0, 0] - 2*sp.diff(a, y)*Yf) == 0
    ok &= all(sp.simplify(Lie[0, m+1] - (a + rho**2*sp.diff(b, y))*sp.diff(Yf, xs[m])) == 0 for m in range(4))
    ok &= all(sp.simplify(Lie[m+1, n_+1] - 2*rho**2*(sp.diff(rho, y)/rho*a - b)*Yf*gam[m, n_]) == 0 for m in range(4) for n_ in range(4))
record("S2a_diffeomorphism_acts_as_xi=a'_psi=Ha-b_dydx=a+rho^2b'", ok)
Hf = sp.diff(rho, y)/rho
a_of_b = -rho**2*sp.diff(b, y)
gauge_cond = sp.diff(a_of_b, y) + 2*(Hf*a_of_b - b)                                  # xi_g + 2 psi_g = 0
b2 = sp.solve(sp.Eq(gauge_cond, 0), sp.diff(b, y, 2))[0]
record("S2b_residual_gauge_ODE_rho^2 b''+4 rho rho' b'+2b=0", sp.simplify(b2 - (-(4*rho*sp.diff(rho, y)*sp.diff(b, y) + 2*b)/rho**2)) == 0)
psi_g = Hf*a_of_b - b; chi_g = sp.diff(ph0, y)*a_of_b
r2 = (sp.diff(rho, y)**2/rho - 1/rho - rho*sp.diff(ph0, y)**2/3); p2 = (U1 - 4*Hf*sp.diff(ph0, y))
def red(e):
    b3 = sp.diff(b2, y)
    e = e.subs(sp.diff(b, y, 4), sp.diff(b3, y)).subs(sp.diff(b, y, 3), b3)
    for _ in range(3): e = e.subs(sp.diff(b, y, 3), b3).subs(sp.diff(b, y, 2), b2)
    for _ in range(3): e = e.subs({sp.diff(rho, y, 3): sp.diff(r2, y), sp.diff(ph0, y, 3): sp.diff(p2, y) + U2*sp.diff(ph0, y)}).subs({sp.diff(rho, y, 2): r2, sp.diff(ph0, y, 2): p2})
    for _ in range(2): e = e.subs(sp.diff(b, y, 2), b2)
    return sp.simplify(e)
mu2 = -4
master = sp.diff(psi_g, y, 2) + 2*(Hf - p2/sp.diff(ph0, y))*sp.diff(psi_g, y) + (-sp.Rational(4, 3)*sp.diff(ph0, y)**2 - 4*Hf*p2/sp.diff(ph0, y) + (2 + mu2)/rho**2)*psi_g
record("S2c_gauge_mode_psi_g_solves_the_master_equation_at_mu2=-4", red(master) == 0)
record("S2d_gauge_mode_obeys_Codazzi_chi_g=-3(psi_g'+2H psi_g)/phi0'", red(chi_g + 3*(sp.diff(psi_g, y) + 2*Hf*psi_g)/sp.diff(ph0, y)) == 0)
yy = sp.symbols('Y', positive=True)
ind = sp.simplify((yy**2*sp.diff(yy**sp.Symbol('q'), yy, 2) + 4*yy*sp.diff(yy**sp.Symbol('q'), yy) + 2*yy**sp.Symbol('q'))/yy**sp.Symbol('q'))
record("S2e_cone_exponents_of_b_are_-1_and_-2_(b~1/y,_a->const_=_rigid_translation_of_the_apex)", set(sp.solve(ind, sp.Symbol('q'))) == {-1, -2})

# S3: displaced shell in the unperturbed bulk
zeta, s0, s1, s2, yb = sp.symbols('zeta sigma0 sigma1 sigma2 y_b')
ok_tr = True; ok_sc = True
for nm, Yf in specials.items():
    n_low = [1] + [-eps*zeta*sp.diff(Yf, xs[m]) for m in range(4)]                     # unit normal to first order
    g0i = g0.inv()
    def chr5(A, B_, C): return sum(sp.Rational(1, 2)*g0i[A, D]*(sp.diff(g0[D, B_], X[C]) + sp.diff(g0[D, C], X[B_]) - sp.diff(g0[B_, C], X[D])) for D in range(5))
    Kmn = sp.Matrix(4, 4, lambda m, n_: sp.diff(n_low[n_+1], xs[m]) - sum(chr5(C, m+1, n_+1)*n_low[C] for C in range(5)))
    shift = lambda e: lin(e + eps*zeta*Yf*sp.diff(e, y))                               # evaluate at y_b + eps zeta Y
    hind_inv = sp.Matrix(4, 4, lambda m, n_: shift(gi4[m, n_]/rho**2))
    Kmix = (hind_inv*Kmn.applyfunc(shift)).applyfunc(lin)
    phi_at = ph0 + eps*zeta*Yf*sp.diff(ph0, y)
    sig = s0 + s1*(phi_at - ph0); dsig = s1 + s2*(phi_at - ph0)
    bgj = {sp.diff(rho, y, 2): r2, sp.diff(ph0, y, 2): p2}
    for m in range(4):
        for n_ in range(4):
            e = lin(Kmix[m, n_] - (sig/6 if m == n_ else 0)).subs(bgj)
            e = e.subs(s0, 6*sp.diff(rho, y)/rho).subs(s1, -2*sp.diff(ph0, y))          # background junctions
            ok_tr &= sp.simplify(e) == 0
    ndphi = shift(sp.diff(ph0, y))                                                      # n^A d_A phi = phi0'(y_b + eps zeta Y) + O(eps^2)
    mism = lin(ndphi + dsig/2).subs(s1, -2*sp.diff(ph0, y))
    Bq = sp.diff(ph0, y, 2)/sp.diff(ph0, y) + s2
    ok_sc &= sp.simplify(mism - eps*zeta*Yf*sp.diff(ph0, y)*Bq) == 0
record("S3a_trace_Israel_condition_holds_identically_for_the_displaced_shell", ok_tr)
record("S3b_scalar_junction_mismatch_=_eps_zeta_Y_phi0'_B", ok_sc)
ok = all(v_["passed"] for v_ in res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok, checks=res), open("SPECIAL_HARMONICS_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok else "SOME CHECKS FAILED"); sys.exit(0 if ok else 1)
