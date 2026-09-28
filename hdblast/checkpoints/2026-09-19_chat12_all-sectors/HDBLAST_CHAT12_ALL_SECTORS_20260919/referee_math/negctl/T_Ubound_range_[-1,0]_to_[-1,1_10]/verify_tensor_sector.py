#!/usr/bin/env python3
"""Chat 12 - TENSOR sector (4D-covariant transverse-traceless modes), exact symbolic checks + exact polynomial root count.
 T1  a TT perturbation  delta g_{x2 x3} = rho^2 e^{2 tau} h(y) T(tau) E(x1)  (E''=-k^2E)  obeys, with the massive-graviton rule
     T'' + 3T' + (k^2 e^{-2tau} + m2) T = 0:   EQ_{x2x3} = -(rho^2 e^{2tau} T E/2) [ h'' + 4H h' + m2 h/rho^2 ],  all other components vanish.
 T2  in the conformal coordinate dz = dy/rho with u = rho^{3/2} h:  -u_zz + V1 u = m2 u,  V1 = (3/2) Hc_z + (9/4) Hc^2,  Hc = rho_z/rho = rho'
     and  -d_z^2 + V1 = Q^+ Q,  Q = d_z - (3/2)Hc  (so m2 >= 0: no tensor tachyon);  the Z2 shell condition h' = 0 is Q u = 0.
 T3  partner potential  V2 = -(3/2)Hc_z + (9/4)Hc^2  satisfies the exact identity   V2 - 9/4 = rho^2 ( 9 phi'^2/16 - U/8 ).
 T4  U(phi) <= -1/6 < 0 on [-1, 0] and U < 0 on [-1, 3/20] for the registered potential (exact Sturm root counts),
     hence V2 > 9/4 wherever phi in [-1, 3/20]: a massive tensor bound state below the continuum threshold 9/4 is impossible
     (a mode v = Q u of Q Q^+ with Dirichlet data at the shell would have m2 = <v,QQ^+v>/<v,v> >= min V2 > 9/4).
Run:  <venv python> verify_tensor_sector.py"""
import json, sys
import sympy as sp
from lin_gr import *
m2 = sp.symbols('m2'); h = sp.Function('h')(y); T = sp.Function('T')(tau); E = sp.Function('E')(x1)
res = {}
def record(n, ok, note=""): res[n] = dict(passed=bool(ok), note=note); print("[%s] %s %s" % ("PASS" if ok else "FAIL", n, note), flush=True)
dg = sp.zeros(5, 5); dg[3, 4] = dg[4, 3] = rho**2*sp.exp(2*tau)*h*T*E
eq = first_order_equations(dg)
def harm(e):
    e = e.subs(sp.diff(T, tau, 2), -3*sp.diff(T, tau) - (k**2*sp.exp(-2*tau) + m2)*T).subs(sp.diff(E, x1, 2), -k**2*E)
    return sp.simplify(bg_reduce(e))
H = rp/rho
target = -(rho**2*sp.exp(2*tau)*T*E/2)*(sp.diff(h, y, 2) + 4*H*sp.diff(h, y) + m2*h/rho**2)
ok23 = sp.simplify(harm(eq[(3, 4)]) - target) == 0
others = all(harm(v) == 0 for key, v in eq.items() if key != (3, 4))
record("T1_TT_mode_radial_equation_h''+4Hh'+m2 h/rho^2=0_and_all_other_components_vanish", ok23 and others)
# T2/T3 in conformal coordinate
u = sp.Function('u')(y); dz = lambda f: rho*sp.diff(f, y)
Hc = sp.diff(rho, y)
V1 = sp.Rational(3, 2)*dz(Hc) + sp.Rational(9, 4)*Hc**2; V2 = -sp.Rational(3, 2)*dz(Hc) + sp.Rational(9, 4)*Hc**2
hh = u/rho**sp.Rational(3, 2)
radial = sp.diff(hh, y, 2) + 4*sp.diff(rho, y)/rho*sp.diff(hh, y) + m2*hh/rho**2
record("T2a_Schroedinger_form_-u_zz+V1u=m2u", sp.simplify(radial*rho**sp.Rational(7, 2) - (dz(dz(u)) - V1*u + m2*u)) == 0)
Q = lambda f: dz(f) - sp.Rational(3, 2)*Hc*f; Qd = lambda f: -dz(f) - sp.Rational(3, 2)*Hc*f
record("T2b_factorisation_Q^+Q=-d_z^2+V1_and_QQ^+=-d_z^2+V2", sp.simplify(Qd(Q(u)) - (-dz(dz(u)) + V1*u)) == 0 and sp.simplify(Q(Qd(u)) - (-dz(dz(u)) + V2*u)) == 0)
record("T2c_shell_condition_h'=0_is_Qu=0", sp.simplify(Q(u) - rho**sp.Rational(5, 2)*sp.diff(hh, y)) == 0)
record("T3_identity_V2-9/4=rho^2(9phi'^2/16-U/8)", sp.simplify(bg_reduce(V2 - sp.Rational(9, 4)) - rho**2*(9*pp**2/16 - (6*(1 - rp**2)/rho**2 + pp**2/2)/8)) == 0)
# T4 exact root counts for the registered potential
p = sp.symbols('p'); W = 1 - p + p**3/3; U = sp.expand(sp.diff(W, p)**2/2 - sp.Rational(2, 3)*W**2)
PU = sp.Poly(U, p); PU6 = sp.Poly(U + sp.Rational(1, 6), p)
n_roots_U = PU.count_roots(-1, sp.Rational(3, 20)); val = U.subs(p, -sp.Rational(1, 2))
inner = sp.Poly(sp.cancel((U + sp.Rational(1, 6))/p), p)                    # U + 1/6 = p * inner(p)
record("T4a_U_has_no_root_on_[-1,3/20]_and_is_negative_there", n_roots_U == 0 and val < 0, "roots: %d, U(-1/2) = %s" % (n_roots_U, val))
record("T4b_U<=-1/6_on_[-1,0]  (U+1/6 = phi*P(phi), P>0 on [-1,0])", inner.count_roots(-1, sp.Rational(1,10)) == 0 and (U + sp.Rational(1, 6)).subs(p, sp.Rational(1, 10)) <= 0, "P roots on [-1,0]: %d" % inner.count_roots(-1, 0))
ok = all(v["passed"] for v in res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok, checks=res), open("TENSOR_SECTOR_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok else "SOME CHECKS FAILED"); sys.exit(0 if ok else 1)
