#!/usr/bin/env python3
"""Referee (math) - closes the 'one explicit polarisation' gap of T1.
GENERAL symmetric gamma-traceless 4D tensor TT_{mu nu}(tau,x1,x2,x3) (nine arbitrary functions, NOT assumed transverse,
NOT assumed to be an eigenfunction), perturbation  delta g_{mu nu} = rho^2 h(y) TT_{mu nu},  delta g_{yA} = 0, delta phi = 0.
Claim checked (exact, all 15 components):
   EQ_{mu nu} = -(rho^2/2)(h'' + 4H h') TT_{mu nu} - (h/2) [ (Box - 2) TT ]_{mu nu} + h nabla_(mu D_nu),    D_nu := nabla^mu TT_{mu nu}
   EQ_{y nu}  =  c1 * h' D_nu   (c1 determined below),      EQ_{yy} = 0.
Hence for ANY transverse-traceless harmonic with (Box - 2) TT = m2 TT:  h'' + 4H h' + m2 h/rho^2 = 0   (all polarisations/helicities).
Also: the explicit harmonic of verify_tensor_sector.py satisfies (Box-2)TT = m2 TT under the script's T-rule (so the script's m2 is the
Fierz-Pauli/Lichnerowicz mass), and negative controls (4H->3H, Box-2 -> Box-3)."""
import json, sys, time
import sympy as sp
from lin_gr import *
t0 = time.time()
res = {}
def record(n, ok, note=""): res[n] = dict(passed=bool(ok), note=note); print("[%s] %s %s  (%.0fs)" % ("PASS" if ok else "FAIL", n, note, time.time()-t0), flush=True)
xs = X[1:]; gi4 = gam.inv()
h = sp.Function('h')(y)
names = {(0,0):'A00',(0,1):'A01',(0,2):'A02',(0,3):'A03',(1,1):'A11',(1,2):'A12',(1,3):'A13',(2,2):'A22',(2,3):'A23'}
TT = sp.zeros(4, 4)
for (m, n_), nm in names.items():
    f = sp.Function(nm)(*xs); TT[m, n_] = f; TT[n_, m] = f
# tracelessness: -T00 + e^{-2tau}(T11+T22+T33) = 0
TT[3, 3] = sp.exp(2*tau)*TT[0, 0] - TT[1, 1] - TT[2, 2]
G = [[[chr4(l, m, n_) for n_ in range(4)] for m in range(4)] for l in range(4)]
def cd2(Tm):   # nabla_l T_{mn}
    return [[[sp.diff(Tm[m, n_], xs[l]) - sum(G[s][l][m]*Tm[s, n_] + G[s][l][n_]*Tm[m, s] for s in range(4)) for n_ in range(4)] for m in range(4)] for l in range(4)]
def cd3(T3):   # nabla_k T_{l m n}
    return [[[[sp.diff(T3[l][m][n_], xs[kk]) - sum(G[s][kk][l]*T3[s][m][n_] + G[s][kk][m]*T3[l][s][n_] + G[s][kk][n_]*T3[l][m][s] for s in range(4)) for n_ in range(4)] for m in range(4)] for l in range(4)] for kk in range(4)]
D1 = cd2(TT)
Dv = [sp.expand(sum(gi4[l, l]*D1[l][l][n_] for l in range(4))) for n_ in range(4)]                  # D_nu
D2 = cd3(D1)
Box = sp.Matrix(4, 4, lambda m, n_: sp.expand(sum(gi4[l, l]*D2[l][l][m][n_] for l in range(4))))
cdD = sp.Matrix(4, 4, lambda m, n_: sp.diff(Dv[n_], xs[m]) - sum(G[s][m][n_]*Dv[s] for s in range(4)))
symD = (cdD + cdD.T)/2
dg = sp.zeros(5, 5)
for m in range(4):
    for n_ in range(4): dg[m+1, n_+1] = rho**2*h*TT[m, n_]
eq = first_order_equations(dg)
print("first-order equations done (%.0fs)" % (time.time()-t0), flush=True)
H = rp/rho
def Z(e): return sp.expand(bg_reduce(sp.expand(e))) == 0
def target(m, n_, cH=4, cL=2):
    return -(rho**2/2)*(sp.diff(h, y, 2) + cH*H*sp.diff(h, y))*TT[m, n_] - (h/2)*(Box[m, n_] - cL*TT[m, n_]) + h*symD[m, n_]
ok = all(Z(bg_reduce(eq[(m+1, n_+1)]) - target(m, n_)) for m in range(4) for n_ in range(m, 4))
record("G1_EQ_munu_=_-(rho^2/2)(h''+4Hh')TT-(h/2)(Box-2)TT+h*sym(nabla D)_for_GENERAL_traceless_TT", ok)
c1s = [c for c in (sp.Rational(1, 2), -sp.Rational(1, 2), 1, -1) if all(Z(bg_reduce(eq[(0, n_+1)]) - c*sp.diff(h, y)*Dv[n_]) for n_ in range(4))]
record("G2_EQ_ynu_proportional_to_h'_D_nu", len(c1s) == 1, "c1 = %s" % c1s)
record("G3_EQ_yy_vanishes_for_traceless_TT", Z(bg_reduce(eq[(0, 0)])))
record("NEG_4H_to_3H_rejected", not all(Z(bg_reduce(eq[(m+1, n_+1)]) - target(m, n_, cH=3)) for m in range(4) for n_ in range(m, 4)))
record("NEG_Box-2_to_Box-3_rejected", not all(Z(bg_reduce(eq[(m+1, n_+1)]) - target(m, n_, cL=3)) for m in range(4) for n_ in range(m, 4)))
# the script's explicit harmonic: (Box - 2) TT = m2 TT  under  T'' + 3T' + (k^2 e^{-2tau} + m2) T = 0,  and it is transverse
m2 = sp.symbols('m2'); T = sp.Function('T')(tau); E = sp.Function('E')(x1)
rep = {TT[m, n_]: 0 for (m, n_) in names}
def put(e):
    e = e.subs({sp.Function(nm)(*xs): (sp.exp(2*tau)*T*E if nm == 'A23' else 0) for nm in names.values()}).doit()
    e = e.subs(sp.diff(T, tau, 2), -3*sp.diff(T, tau) - (k**2*sp.exp(-2*tau) + m2)*T).subs(sp.diff(E, x1, 2), -k**2*E)
    return sp.simplify(e)
okL = all(put(Box[m, n_] - 2*TT[m, n_] - m2*TT[m, n_]) == 0 for m in range(4) for n_ in range(m, 4))
okD = all(put(Dv[n_]) == 0 for n_ in range(4))
record("G4_script_harmonic_is_transverse_and_(Box-2)TT=m2_TT_(m2_is_the_Fierz-Pauli_mass)", okL and okD)
ok = all(v["passed"] for v in res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok, checks=res, seconds=time.time()-t0), open("GENERAL_TENSOR_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok else "SOME CHECKS FAILED"); sys.exit(0 if ok else 1)
