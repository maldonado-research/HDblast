#!/usr/bin/env python3
"""HDBLAST Chat 11 - SYMBOLIC verification of the shell boundary-condition derivation (SymPy, exact).
General scalar perturbation of  ds^2 = dy^2 + rho^2 gamma  (gamma = unit dS_4, flat slicing), harmonic Y = T(tau)E(x1), Box Y = mu2 Y:
   g_yy = 1 + 2 xi Y,  g_{y mu} = Bs d_mu Y,  g_{mu nu} = rho^2[(1 + 2 psi Y) gamma_{mu nu} + 2 Es Hess_{mu nu}Y],  phi = phi0 + chi Y.
 (J1) gauge rules: under  x^A -> x^A + v^A,  v^y = a(y) Y,  v^mu = b(y) gamma^{mu nu} d_nu Y, the background changes by the Lie derivative
        delta g_yy = 2 a' Y,  delta g_{y mu} = (a + rho^2 b') d_mu Y,  delta g_{mu nu} = 2 rho rho' a Y gamma_{mu nu} + 2 rho^2 b Hess_{mu nu} Y,  delta phi = phi0' a Y.
 (J2) Gaussian-normal gauge (xi = Bs = 0):  delta K^mu_nu = psi' Y delta^mu_nu + Es' Hess^mu_nu Y   (K_{mu nu} = (1/2) d_y g_{mu nu}).
 (J3) special harmonics: the gamma-traceless Hessian of a homogeneous harmonic vanishes iff T'' = T', which forces mu2 = -4 (or Y const):
        only for those can the shell bend in longitudinal gauge.
 (J4) algebra of the transformation longitudinal -> Gaussian-normal and the resulting shell conditions
        (a' = xi, b' = -a/rho^2, traceless Israel => a(y_b) = 0; trace Israel == Codazzi; scalar junction => chi' - phi0' xi = -(s2/2) chi).
Run:  ../../.hdblast_venv/bin/python verify_junction_logic.py"""
import json, sys, time
import sympy as sp
t0 = time.time()
y, tau, x1, x2, x3, k, mu2 = sp.symbols('y tau x1 x2 x3 k mu2'); X = [y, tau, x1, x2, x3]
rho, ph0, a, b, psi, Es = [sp.Function(n)(y) for n in ('rho', 'phi0', 'a', 'b', 'psi', 'Es')]
T = sp.Function('T')(tau); E = sp.Function('E')(x1); Y = T*E
gam = sp.diag(-1, sp.exp(2*tau), sp.exp(2*tau), sp.exp(2*tau)); gami = gam.inv(); xs = X[1:]
# Christoffels and Hessian on dS_4
def chr4(l, m, n_): return sum(sp.Rational(1, 2)*gami[l, s]*(sp.diff(gam[s, m], xs[n_]) + sp.diff(gam[s, n_], xs[m]) - sp.diff(gam[m, n_], xs[s])) for s in range(4))
Hess = sp.Matrix(4, 4, lambda m, n_: sp.diff(Y, xs[m], xs[n_]) - sum(chr4(l, m, n_)*sp.diff(Y, xs[l]) for l in range(4)))
boxY = sp.simplify(sum(gami[m, n_]*Hess[m, n_] for m in range(4) for n_ in range(4)))
res = {}
def record(name, ok, note=""):
    res[name] = dict(passed=bool(ok), note=note); print("[%s] %s %s" % ("PASS" if ok else "FAIL", name, note), flush=True)
def Z(e): return sp.simplify(e) == 0
# sanity: unit dS_4 and the harmonic rule
Ric4 = sp.Matrix(4, 4, lambda m, n_: sum(sp.diff(chr4(l, m, n_), xs[l]) - sp.diff(chr4(l, m, l), xs[n_]) + sum(chr4(l, l, s)*chr4(s, m, n_) - chr4(l, n_, s)*chr4(s, m, l) for s in range(4)) for l in range(4)))
record("S0_gamma_is_unit_de_Sitter_R=3gamma", all(Z(Ric4[m, n_] - 3*gam[m, n_]) for m in range(4) for n_ in range(4)))
record("S1_box_Y_formula", Z(boxY - (-sp.diff(T, tau, 2) - 3*sp.diff(T, tau))*E - sp.exp(-2*tau)*T*sp.diff(E, x1, 2)))

# (J1) Lie derivative of the background
g0 = sp.zeros(5, 5); g0[0, 0] = 1
for m in range(4):
    for n_ in range(4): g0[m+1, n_+1] = rho**2*gam[m, n_]
v = [a*Y] + [b*sum(gami[m, n_]*sp.diff(Y, xs[n_]) for n_ in range(4)) for m in range(4)]
Lie = sp.Matrix(5, 5, lambda A, B_: sum(v[C]*sp.diff(g0[A, B_], X[C]) + g0[C, B_]*sp.diff(v[C], X[A]) + g0[A, C]*sp.diff(v[C], X[B_]) for C in range(5)))
ok = Z(Lie[0, 0] - 2*sp.diff(a, y)*Y)
ok &= all(Z(Lie[0, m+1] - (a + rho**2*sp.diff(b, y))*sp.diff(Y, xs[m])) for m in range(4))
ok &= all(Z(Lie[m+1, n_+1] - (2*rho*sp.diff(rho, y)*a*Y*gam[m, n_] + 2*rho**2*b*Hess[m, n_])) for m in range(4) for n_ in range(4))
ok &= Z(v[0]*sp.diff(ph0, y) - sp.diff(ph0, y)*a*Y)
record("J1_gauge_transformation_rules_(Lie_derivative)", ok)

# (J2) first-order extrinsic curvature in Gaussian-normal gauge
eps = sp.symbols('epsilon')
h0 = rho**2*gam; dh = 2*rho**2*(psi*Y*gam + Es*Hess)
h0i = h0.inv(); dhi = -h0i*dh*h0i
dK = sp.Rational(1, 2)*(dhi*sp.diff(h0, y) + h0i*sp.diff(dh, y))                     # delta K^mu_nu
target = sp.diff(psi, y)*Y*sp.eye(4) + sp.diff(Es, y)*(gami*Hess)
record("J2_delta_K_mixed_=_psi'_Y_delta_+_Es'_Hess", all(Z(dK[m, n_] - target[m, n_]) for m in range(4) for n_ in range(4)))

# (J3) special harmonics: homogeneous Y = T(tau)
Yh = T; Hh = sp.Matrix(4, 4, lambda m, n_: sp.diff(Yh, xs[m], xs[n_]) - sum(chr4(l, m, n_)*sp.diff(Yh, xs[l]) for l in range(4)))
boxh = sum(gami[m, m]*Hh[m, m] for m in range(4)); TF = Hh - sp.Rational(1, 4)*gam*boxh
T1, T2 = sp.symbols('T1 T2'); TFs = TF.subs({sp.diff(T, tau, 2): T2, sp.diff(T, tau): T1})
prop = all(Z(TFs[m, m] - sp.Rational(3, 4)*(T2 - T1)*(1 if m == 0 else sp.exp(2*tau)/3)) for m in range(4))
mu_special = sp.simplify(boxh.subs(T, sp.exp(tau)).doit()/sp.exp(tau))
record("J3_tracefree_Hessian_prop_to_(T''-T')_and_T=e^tau_has_mu2=-4", prop and mu_special == -4, "Box e^tau = %s e^tau" % mu_special)

# (J4) longitudinal -> Gaussian-normal algebra (functions of y only)
xi, chi = sp.Function('xi')(y), sp.Function('chi')(y); s1, s2 = sp.symbols('sigma1 sigma2')     # sigma_t', sigma_t'' at the shell
H = sp.diff(rho, y)/rho
psib = psi - H*a; chib = chi - sp.diff(ph0, y)*a                                              # barred (GN) fields, with a' = xi
dpsib = sp.diff(psib, y).subs(sp.diff(a, y), xi); dchib = sp.diff(chib, y).subs(sp.diff(a, y), xi)
at0 = lambda e: e.subs(a, 0)                                                                  # no bending: a(y_b) = 0
trace_israel = at0(dpsib) - s1*at0(chib)/6
codazzi = (sp.diff(psi, y) - H*xi) + sp.diff(ph0, y)*chi/3
record("J4a_trace_Israel_condition_is_Codazzi_given_phi0'=-sigma'/2", Z((trace_israel - codazzi).subs(sp.diff(ph0, y), -s1/2)))
scalar_j = at0(dchib) + s2*at0(chib)/2
record("J4b_scalar_junction_is_chi'-phi0'xi+(s2/2)chi", Z(scalar_j - (sp.diff(chi, y) - sp.diff(ph0, y)*xi + s2*chi/2)))
record("J4c_with_xi=-2psi_gives_chi'+2phi0'psi+(s2/2)chi", Z(scalar_j.subs(xi, -2*psi) - (sp.diff(chi, y) + 2*sp.diff(ph0, y)*psi + s2*chi/2)))

ok_all = all(v_["passed"] for v_ in res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok_all, runtime_s=round(time.time() - t0, 1), checks=res), open("JUNCTION_LOGIC_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok_all else "SOME CHECKS FAILED", "(%.1f s)" % (time.time() - t0)); sys.exit(0 if ok_all else 1)
