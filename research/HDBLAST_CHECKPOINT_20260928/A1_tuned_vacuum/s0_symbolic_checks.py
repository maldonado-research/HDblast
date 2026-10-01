#!/usr/bin/env python3
"""A1 step 0: exact (symbolic) checks that the new solver relies on.  Output: S0_SYMBOLIC_CHECKS.json

 (1) Conformal-gauge Einstein-scalar system of the registered model (kappa_5^2 = 1):
        ds^2 = e^{2B}(-dT^2 + dZ^2) + e^{2A} dx_3^2 ,  R_ab = phi_a phi_b + (2/3) U g_ab ,  box phi = U'(phi).
     Checks: the three evolution equations used by the solver, and the two constraints H, M (Chat 13 definitions)
     are exact linear combinations of the field equations.
 (2) Wrong-formula controls: flipping one coefficient in each evolution equation must FAIL.
 (3) Reference solution in an arbitrary pair of null labels (the chart change of the proper-clock note,
     U -> U_K(U), V -> V_K(V) arbitrary):
        A = l(w)/2 + ln V_K(V),  B = l(w)/2 + ln(U_K'(U) V_K'(V))/2,  phi = Phi(w),  w = -U_K(U) V_K(V)
     solves the full PDE system iff (l, Phi) solve the static ODE in the invariant w
        2 w l'' = -8 l' - 3 w l'^2 - (2/3) e^l U(Phi),     4 w Phi'' = e^l U'(Phi) - 10 Phi' - 6 w l' Phi',
        2 l' + w l'^2 = w Phi'^2/3 - e^l U(Phi)/6         (first-order constraint).
     This is the O(4,1)-symmetric static shell family written in Kruskal-type null coordinates; it is regular across
     w = 0 (the light cone of the static solution's vertex), which is where the old conformal chart freezes.
 (4) Constraint transport: (d_T - d_Z)[e^{3A}(H+2M)] = 0 and (d_T + d_Z)[e^{3A}(H-2M)] = 0 on solutions (basis of Remedy B).
 (5) Delta N_eff <-> rho_DR/rho_SM mapping (exact rational arithmetic; see REGISTRATION.md section 2).
All checks are exact (sympy simplification to 0 or Fraction arithmetic); no floating point tolerance.
"""
import json, hashlib
from fractions import Fraction as Fr
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
out = {'checks': [], 'status': None}
def rec(name, ok, detail=''):
    out['checks'].append(dict(name=name, result='PASS' if ok else 'FAIL', detail=str(detail)))
    print(('PASS ' if ok else 'FAIL ') + name, detail)
    return ok

T, Z = sp.symbols('T Z', real=True)
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
A = sp.Function('A')(T, Z); B = sp.Function('B')(T, Z); ph = sp.Function('phi')(T, Z)
Uf = sp.Function('U'); U = Uf(ph); Up = sp.diff(Uf(sp.Symbol('p')), sp.Symbol('p')).subs(sp.Symbol('p'), ph)
X = [T, Z, x1, x2, x3]
g = sp.diag(-sp.exp(2*B), sp.exp(2*B), sp.exp(2*A), sp.exp(2*A), sp.exp(2*A))
gi = g.inv()
n = 5
Gam = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n))/2)
         for c in range(n)] for b in range(n)] for a in range(n)]
def Ric(b, c):
    r = 0
    for a in range(n):
        r += sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c])
        for d in range(n):
            r += Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a]
    return sp.simplify(r)
dph = [sp.diff(ph, x) for x in X]
E = {}
for (b, c) in [(0, 0), (0, 1), (1, 1), (2, 2)]:
    E[(b, c)] = sp.simplify(Ric(b, c) - dph[b]*dph[c] - sp.Rational(2, 3)*U*g[b, c])
sqrtg = sp.exp(B)*sp.exp(B)*sp.exp(3*A)
box = sum(sp.diff(sqrtg*gi[a, a]*sp.diff(ph, X[a]), X[a]) for a in range(2))/sqrtg
KG = sp.simplify(box - Up)

At, Az, Bt, Bz, ft, fz = [sp.diff(A, T), sp.diff(A, Z), sp.diff(B, T), sp.diff(B, Z), sp.diff(ph, T), sp.diff(ph, Z)]
Att, Azz, Btt, Bzz, ftt, fzz = [sp.diff(A, T, 2), sp.diff(A, Z, 2), sp.diff(B, T, 2), sp.diff(B, Z, 2), sp.diff(ph, T, 2), sp.diff(ph, Z, 2)]
e2B = sp.exp(2*B)
# solver's evolution equations
evoA = Att - (Azz - 3*At**2 + 3*Az**2 + sp.Rational(2, 3)*e2B*U)
evoB = Btt - (Bzz + 3*At**2 - 3*Az**2 - ft**2/2 + fz**2/2 - sp.Rational(1, 3)*e2B*U)
evoF = ftt - (fzz - 3*At*ft + 3*Az*fz - e2B*Up)
# constraints (Chat 13 definitions)
Mc = -3*sp.diff(A, T, Z) - 3*At*Az + 3*At*Bz + 3*Az*Bt - ft*fz
Hc = -2*U*e2B + 6*At**2 + 6*At*Bt - 12*Az**2 + 6*Az*Bz - 6*Azz - ft**2 - fz**2

# (1) express each as a combination of E_ab and KG: solve linear system for coefficients (rational functions)
a0, a1, a2, a3, a4 = sp.symbols('a0:5')
def combo_check(name, target):
    expr = a0*E[(0, 0)] + a1*E[(1, 1)] + a2*E[(2, 2)] + a3*E[(0, 1)] + a4*KG - target
    expr = sp.expand(sp.simplify(expr*sp.exp(-2*A)*sp.exp(2*A)))
    # collect in the independent derivative monomials
    derivs = [Att, Azz, Btt, Bzz, ftt, fzz, sp.diff(A, T, Z), sp.diff(B, T, Z), sp.diff(ph, T, Z), At, Az, Bt, Bz, ft, fz, U, Up]
    syms = sp.symbols('d0:%d' % len(derivs))
    e2 = expr.subs(dict(zip(derivs, syms)))
    poly = sp.Poly(sp.expand(e2), *syms)
    eqs = [sp.simplify(c) for c in poly.coeffs()]
    sol = sp.solve(eqs, [a0, a1, a2, a3, a4], dict=True)
    ok = len(sol) > 0
    return rec(name, ok, sol[0] if ok else 'no combination')

combo_check('evolution A_TT is a combination of field equations', evoA)
combo_check('evolution B_TT is a combination of field equations', evoB)
combo_check('evolution phi_TT is a combination of field equations', evoF)
combo_check('momentum constraint M is a combination of field equations', Mc)
combo_check('Hamiltonian constraint H is a combination of field equations', Hc)
# (2) wrong-formula controls
def control(name, target):
    expr = a0*E[(0, 0)] + a1*E[(1, 1)] + a2*E[(2, 2)] + a3*E[(0, 1)] + a4*KG - target
    derivs = [Att, Azz, Btt, Bzz, ftt, fzz, sp.diff(A, T, Z), sp.diff(B, T, Z), sp.diff(ph, T, Z), At, Az, Bt, Bz, ft, fz, U, Up]
    syms = sp.symbols('d0:%d' % len(derivs))
    poly = sp.Poly(sp.expand(sp.simplify(expr).subs(dict(zip(derivs, syms)))), *syms)
    sol = sp.solve([sp.simplify(c) for c in poly.coeffs()], [a0, a1, a2, a3, a4], dict=True)
    return rec('CONTROL (must have no combination): ' + name, len(sol) == 0, 'no combination found' if not sol else sol)
control('A_TT with 3A_Z^2 -> 2A_Z^2', Att - (Azz - 3*At**2 + 2*Az**2 + sp.Rational(2, 3)*e2B*U))
control('B_TT with -phi_T^2/2 -> -phi_T^2', Btt - (Bzz + 3*At**2 - 3*Az**2 - ft**2 + fz**2/2 - sp.Rational(1, 3)*e2B*U))
control('phi_TT with -3A_T phi_T -> -4A_T phi_T', ftt - (fzz - 4*At*ft + 3*Az*fz - e2B*Up))
control('H with -6A_ZZ -> -5A_ZZ', Hc + Azz)

# (1b) characteristic constraint transport (undamped): (d_T - d_Z)[e^{3A}(H+2M)] = 0 and (d_T + d_Z)[e^{3A}(H-2M)] = 0
#      on solutions of the evolution equations (Remedy B damps the outgoing combination C+ = H + 2M).
rules = {}
AttE = Azz - 3*At**2 + 3*Az**2 + sp.Rational(2, 3)*e2B*U
BttE = Bzz + 3*At**2 - 3*Az**2 - ft**2/2 + fz**2/2 - sp.Rational(1, 3)*e2B*U
fttE = fzz - 3*At*ft + 3*Az*fz - e2B*Up
def reduce_tt(expr):
    # replace every derivative with >= 2 time derivatives using the evolution equations (repeatedly)
    for _ in range(4):
        reps = {}
        for F, RHS in [(A, AttE), (B, BttE), (ph, fttE)]:
            for kz in range(0, 3):
                d = sp.diff(F, T, 2, Z, kz) if kz else sp.diff(F, T, 2)
                reps[d] = sp.diff(RHS, Z, kz) if kz else RHS
                for kt in range(3, 4):
                    d3 = sp.diff(F, T, kt, Z, kz) if kz else sp.diff(F, T, kt)
                    reps[d3] = sp.diff(RHS, T, kt - 2, Z, kz) if kz else sp.diff(RHS, T, kt - 2)
        expr = expr.subs(reps)
    return sp.simplify(sp.expand(expr))
for sgn, nm in [(+1, 'outgoing C+ = H+2M along d_T - d_Z'), (-1, 'ingoing C- = H-2M along d_T + d_Z')]:
    Cc = sp.exp(3*A)*(Hc + 2*sgn*Mc)
    tr = sp.diff(Cc, T) - sgn*sp.diff(Cc, Z)
    r = reduce_tt(tr)
    rec('constraint transport identity: ' + nm, r == 0, 'residual=' + str(r)[:80])

# (3) reference ansatz in arbitrary null labels
Uu, Vv = sp.symbols('Uu Vv', real=True)
UK = sp.Function('UK')(Uu); VK = sp.Function('VK')(Vv)
wsym = sp.Symbol('w')
l = sp.Function('l'); Phi = sp.Function('Phi')
w = -UK*VK
Aref = l(w)/2 + sp.log(VK)
Bref = l(w)/2 + sp.log(sp.diff(UK, Uu)*sp.diff(VK, Vv))/2
Fref = Phi(w)
# T = (U+V)/2, Z = (V-U)/2  =>  d_T = d_U + d_V, d_Z = d_V - d_U
def dT(e): return sp.diff(e, Uu) + sp.diff(e, Vv)
def dZ(e): return sp.diff(e, Vv) - sp.diff(e, Uu)
Us = Uf(Fref); Ups = sp.diff(Uf(sp.Symbol('p')), sp.Symbol('p')).subs(sp.Symbol('p'), Fref)
e2Bs = sp.exp(2*Bref)
resA = dT(dT(Aref)) - (dZ(dZ(Aref)) - 3*dT(Aref)**2 + 3*dZ(Aref)**2 + sp.Rational(2, 3)*e2Bs*Us)
resB = dT(dT(Bref)) - (dZ(dZ(Bref)) + 3*dT(Aref)**2 - 3*dZ(Aref)**2 - dT(Fref)**2/2 + dZ(Fref)**2/2 - sp.Rational(1, 3)*e2Bs*Us)
resF = dT(dT(Fref)) - (dZ(dZ(Fref)) - 3*dT(Aref)*dT(Fref) + 3*dZ(Aref)*dZ(Fref) - e2Bs*Ups)
resM = -3*dT(dZ(Aref)) - 3*dT(Aref)*dZ(Aref) + 3*dT(Aref)*dZ(Bref) + 3*dZ(Aref)*dT(Bref) - dT(Fref)*dZ(Fref)
resH = -2*Us*e2Bs + 6*dT(Aref)**2 + 6*dT(Aref)*dT(Bref) - 12*dZ(Aref)**2 + 6*dZ(Aref)*dZ(Bref) - 6*dZ(dZ(Aref)) - dT(Fref)**2 - dZ(Fref)**2
# substitute ODE: l'' and Phi'' from the second-order ODEs, and l'^... keep l' free; constraint used separately
lp, lpp, Pp, Ppp, lv, Pv = sp.symbols('lp lpp Pp Ppp lv Pv')
# generic substitution via replacing function by symbols
wz = sp.Symbol('wz')
def subs_funcs(e):
    e = e.doit()
    reps = {}
    for sb in e.atoms(sp.Subs):
        der = sb.expr
        if isinstance(der, sp.Derivative):
            fn = der.expr.func; order = sum(c for _, c in der.variable_count)
            key = {(l, 1): lp, (l, 2): lpp, (Phi, 1): Pp, (Phi, 2): Ppp}[(fn, order)]
            reps[sb] = key
    e = e.xreplace(reps)
    for d in list(e.atoms(sp.Derivative)):
        if d.expr in (l(w), Phi(w)):
            order = sum(c for _, c in d.variable_count)
            raise RuntimeError('unexpected derivative form')
    e = e.xreplace({l(w): lv, Phi(w): Pv})
    return e
lpp_ode = (-8*lp - 3*wsym*lp**2 - sp.Rational(2, 3)*sp.exp(lv)*Uf(Pv))/(2*wsym)
Upv = sp.diff(Uf(sp.Symbol('p')), sp.Symbol('p')).subs(sp.Symbol('p'), Pv)
Ppp_ode = (sp.exp(lv)*Upv - 10*Pp - 6*wsym*lp*Pp)/(4*wsym)
cons = 2*lp + wsym*lp**2 - (wsym*Pp**2/3 - sp.exp(lv)*Uf(Pv)/6)
uk, vk, u1, u2, u3, v1, v2, v3 = sp.symbols('uk vk u1 u2 u3 v1 v2 v3')
def jets(e):
    # replace the null-label maps and their derivatives by independent symbols (BEFORE eliminating vk)
    e = e.xreplace({sp.Derivative(UK, (Uu, 3)): u3, sp.Derivative(VK, (Vv, 3)): v3})
    e = e.xreplace({sp.Derivative(UK, (Uu, 2)): u2, sp.Derivative(VK, (Vv, 2)): v2})
    e = e.xreplace({sp.Derivative(UK, Uu): u1, sp.Derivative(VK, Vv): v1})
    e = e.xreplace({UK: uk, VK: vk})
    assert not e.has(Uu) and not e.has(Vv), 'unreplaced label dependence'
    return e.subs(vk, -wsym/uk)
def reduce_ref(rr, lpp_rule):
    e = jets(subs_funcs(sp.expand(rr)))
    return sp.simplify(e.subs({lpp: lpp_rule, Ppp: Ppp_ode}))
for nm, rr in [('A', resA), ('B', resB), ('phi', resF), ('M', resM), ('H', resH)]:
    e0 = jets(subs_funcs(sp.expand(rr)))
    e = reduce_ref(rr, lpp_ode)
    if e == 0:
        rec('reference ansatz solves %s-equation given the second-order w-ODE (depends on l\'\'/Phi\'\': %s)' % (nm, e0.has(lpp) or e0.has(Ppp)), True, '')
    else:
        q = sp.simplify(e/cons)
        ok2 = (not q.has(lp)) and (not q.has(Pp)) and sp.simplify(e - q*cons) == 0
        rec('reference ansatz solves %s-equation modulo the first-order constraint' % nm, bool(ok2), 'residual = (' + str(q)[:100] + ') x constraint')
# controls: wrong ODE / wrong constraint coefficient must fail
e = reduce_ref(resA, (-7*lp - 3*wsym*lp**2 - sp.Rational(2, 3)*sp.exp(lv)*Uf(Pv))/(2*wsym))
rec("CONTROL (must fail): w-ODE with -8 l' -> -7 l' does not solve the A-equation", e != 0, str(e)[:80])
e = reduce_ref(resH, lpp_ode)
cons_wrong = 2*lp + wsym*lp**2 - (wsym*Pp**2/4 - sp.exp(lv)*Uf(Pv)/6)
q = sp.simplify(e/cons_wrong)
rec("CONTROL (must fail): H-equation is not a multiple of a constraint with Phi'^2/3 -> Phi'^2/4", q.has(lp) or q.has(Pp) or q.has(lv), '')

# (5) Delta N_eff mapping, exact rationals
gBBN = Fr(43, 4)                   # 10.75: photons + e+- + 3 nu at T ~ 1 MeV
gprod = Fr(427, 4)                 # 106.75: full Standard Model
rho_nu1_over_rho_gamma_BBN = Fr(7, 8)       # one species (nu+nubar) at T_nu = T
rho_SM_over_rho_gamma_BBN = gBBN/2
ratio_per_dN_BBN = rho_nu1_over_rho_gamma_BBN/rho_SM_over_rho_gamma_BBN      # = 7/43
rec('rho_DR/rho_SM at BBN per unit Delta N_eff = 7/43', ratio_per_dN_BBN == Fr(7, 43), float(ratio_per_dN_BBN))
# production -> BBN: rho_DR a^4 const, rho_SM a^4 ∝ g* g_s^{-4/3} = g*^{-1/3} (g* = g_s); ratio(prod) = ratio(BBN) * (gprod/gBBN)^{1/3}
x = sp.Rational(427, 43)
fac = sp.nsimplify(x)**sp.Rational(1, 3)
out['neff_mapping'] = dict(
    rho_DR_over_rho_SM_BBN_per_dNeff='7/43', value_per_dNeff_BBN=float(ratio_per_dN_BBN),
    production_factor='(427/43)^(1/3) = (g*_prod/g*_BBN)^(1/3)', production_factor_value=float(fac),
    rho_DR_over_rho_SM_prod_per_dNeff=float(ratio_per_dN_BBN*float(fac)),
    rho_DR_over_rho_gamma_CMB_per_dNeff='(7/8)(4/11)^(4/3)', value_CMB=float(sp.Rational(7, 8)*sp.Rational(4, 11)**sp.Rational(4, 3)),
    thresholds={}
)
for r in [0.1, 0.03]:
    dN = r/(float(ratio_per_dN_BBN)*float(fac))
    out['neff_mapping']['thresholds'][str(r)] = dict(rho_DR_over_rho_SM_at_production=r, equivalent_Delta_N_eff=dN,
        rho_DR_over_rho_SM_at_BBN=dN*float(ratio_per_dN_BBN), rho_DR_over_rho_gamma_at_CMB=dN*float(sp.Rational(7, 8)*sp.Rational(4, 11)**sp.Rational(4, 3)))
rec('M5 cross-check: Delta N_eff=0.3 -> 0.10497 at production (g*=106.75)', abs(0.3*float(ratio_per_dN_BBN)*float(fac) - 0.1049713169074991) < 1e-12, 0.3*float(ratio_per_dN_BBN)*float(fac))

out['status'] = 'PASS' if all(c['result'] == 'PASS' for c in out['checks']) else 'FAIL'
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
out['label'] = 'exact-verified (sympy / exact rationals); the N_eff inputs themselves are snippet-level (not verified against papers)'
(HERE/'S0_SYMBOLIC_CHECKS.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
print('status', out['status'])
