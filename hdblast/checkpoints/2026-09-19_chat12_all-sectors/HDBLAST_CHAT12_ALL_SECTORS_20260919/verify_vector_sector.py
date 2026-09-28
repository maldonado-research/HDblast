#!/usr/bin/env python3
"""Chat 12 - VECTOR sector (4D-covariant transverse vector harmonics V_mu on unit dS_4), exact symbolic checks.
General vector perturbation:  delta g_{y mu} = Bv(y) V_mu,   delta g_{mu nu} = 2 rho^2 Fv(y) nabla_(mu V_nu),   delta phi = 0.
 V1  gauge rule: under x^mu -> x^mu + zeta(y) V^mu:  Fv -> Fv + zeta,  Bv -> Bv + rho^2 zeta'   (so sigma_V := Bv - rho^2 Fv' is gauge invariant,
     and the gauge Fv = 0 is always reachable; the shell position is untouched because v^y = 0).
 V2  explicit transverse harmonic V = (0,0,T(tau)E(x1),0) (E'' = -k^2 E); in the gauge Fv = 0 the first-order field equations are
        EQ_{x1 x2} = +(E' T/2)      (Bv' + 2 H Bv),
        EQ_{tau x2} = +(E (T'-2T)/2)(Bv' + 2 H Bv),
        EQ_{y x2}  = -(Bv/(2 rho^2)) [ (Box + 3) V ]_{x2},         all other components vanish identically,
     where (Box + 3)V_nu = nabla^mu(nabla_mu V_nu + nabla_nu V_mu) for transverse V on unit dS_4.
 V3  consequences (elementary): if (Box+3)V != 0 then Bv = 0: no vector perturbation at all.  If (Box+3)V = 0 then Bv = C/rho^2, which
     is singular at the regular cone (rho -> 0) unless C = 0.  Either way the vector sector contains no regular mode, with or without the shell.
Run:  <venv python> verify_vector_sector.py"""
import json, sys
import sympy as sp
from lin_gr import *
Bv, Fv, zeta = [sp.Function(n)(y) for n in ('Bv', 'Fv', 'zeta')]
T = sp.Function('T')(tau); E = sp.Function('E')(x1)
res = {}
def record(n, ok, note=""): res[n] = dict(passed=bool(ok), note=note); print("[%s] %s %s" % ("PASS" if ok else "FAIL", n, note), flush=True)
xs = X[1:]; gi4 = gam.inv()
Vlow = [0, 0, T*E, 0]                                                 # V_mu (index down), only the x2 component
def cov_d(Vl):                                                         # nabla_mu V_nu
    return sp.Matrix(4, 4, lambda m, n_: sp.diff(Vl[n_], xs[m]) - sum(chr4(l, m, n_)*Vl[l] for l in range(4)))
DV = cov_d(Vlow); S = (DV + DV.T)                                      # 2 nabla_(mu V_nu)
div = sp.simplify(sum(gi4[m, m]*DV[m, m] for m in range(4)))
record("V0_harmonic_is_transverse", div == 0)
# V1 gauge rule via the Lie derivative of the background
g0 = background_metric(); Vup = [sum(gi4[m, n_]*Vlow[n_] for n_ in range(4)) for m in range(4)]
v = [0] + [zeta*c for c in Vup]
Lie = sp.Matrix(5, 5, lambda A, B_: sum(v[C]*sp.diff(g0[A, B_], X[C]) + g0[C, B_]*sp.diff(v[C], X[A]) + g0[A, C]*sp.diff(v[C], X[B_]) for C in range(5)))
ok = all(sp.simplify(Lie[0, m+1] - rho**2*sp.diff(zeta, y)*Vlow[m]) == 0 for m in range(4))
ok &= all(sp.simplify(Lie[m+1, n_+1] - rho**2*zeta*S[m, n_]) == 0 for m in range(4) for n_ in range(4)) and sp.simplify(Lie[0, 0]) == 0
record("V1_gauge_rule_Fv->Fv+zeta_Bv->Bv+rho^2 zeta'", ok)
# V2 field equations in the gauge Fv = 0
dg = sp.zeros(5, 5); dg[0, 3] = dg[3, 0] = Bv*T*E
eq = first_order_equations(dg)
red = lambda e: sp.simplify(bg_reduce(e))
H = rp/rho; C1 = sp.diff(Bv, y) + 2*H*Bv
ok12 = sp.simplify(red(eq[(2, 3)]) - sp.diff(E, x1)*T/2*C1) == 0
ok13 = sp.simplify(red(eq[(1, 3)]) - E*(sp.diff(T, tau) - 2*T)/2*C1) == 0
record("V2a_EQ_{x1x2}_and_EQ_{tau x2}_proportional_to_(Bv'+2H Bv)", ok12 and ok13)
# 4D operator (Box + 3) V on the x2 component:  nabla^mu (nabla_mu V_nu + nabla_nu V_mu)
def cov_d2(Tm):                                                        # nabla_l T_{m n}
    return [[[sp.diff(Tm[m, n_], xs[l]) - sum(chr4(s, l, m)*Tm[s, n_] + chr4(s, l, n_)*Tm[m, s] for s in range(4)) for n_ in range(4)] for m in range(4)] for l in range(4)]
DS = cov_d2(S); op = [sp.simplify(sum(gi4[l, l]*DS[l][l][n_] for l in range(4))) for n_ in range(4)]
record("V2b_EQ_{y x2}=-(Bv/(2rho^2))[(Box+3)V]_{x2}", sp.simplify(red(eq[(0, 3)]) + Bv/(2*rho**2)*op[2]) == 0, "[(Box+3)V]_x2 = %s" % sp.factor(op[2]))
record("V2c_other_components_of_(Box+3)V_vanish_for_this_harmonic", all(op[n_] == 0 for n_ in (0, 1, 3)))
record("V2d_all_other_field_equation_components_vanish", all(red(val) == 0 for key, val in eq.items() if key not in ((2, 3), (1, 3), (0, 3))))
# V3 the constraint integrates to Bv = C/rho^2
record("V3_(Bv'+2H Bv)=0_<=>_(rho^2 Bv)'=0", sp.simplify(sp.diff(rho**2*Bv, y) - rho**2*(sp.diff(Bv, y) + 2*sp.diff(rho, y)/rho*Bv)) == 0)
ok = all(v_["passed"] for v_ in res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok, checks=res), open("VECTOR_SECTOR_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok else "SOME CHECKS FAILED"); sys.exit(0 if ok else 1)
