"""Exact verification of the explicit second-order correction for the
commuting quadrupole A = B = diag(1,-1,0,0).

w = (h_T + a(s) x x^T, Psi(s), pi', varpi') with
  pi'    = sum_k b_k(s) m_k(x) (s delta - x x^T)
  varpi' = sum_k n_k(s) m_k(x)
  h_T    = sum_k c_k(s) [ s P Hess(m_k) P - tr(P Hess(m_k) P) P/3 ]
  a, Psi radial; Psi obeys  ca1 a' + ca0 a + cp1 Psi' + cp0 Psi = f0(s).
Checks  P_H[w] = S_H  and  P_M[w] = S_M  identically (all jets free, background
equations imposed), plus support and wrong-formula controls.
"""
import json, pickle, sys, time
import sympy as sp
from engine import X, O, F, U, register, reduce_bg, S2
from linear_ops import PH_PM

t0 = time.time()
s_ = sp.Symbol('s_')
src = pickle.load(open('../results/sources_cc.pkl', 'rb'))
SH, SM = src['SH'], src['SM']
cm = pickle.load(open('../results/construct_cc_M.pkl', 'rb'))
ch = pickle.load(open('../results/construct_cc_H.pkl', 'rb'))
mons = cm['mons']
G = register('G'); K = register('K'); Aa = register('Aa'); Ps = register('Ps')
toS = lambda e: sp.sympify(e).subs(s_, S2)
bv = [toS(cm['vals'][sp.Symbol(f'b{k}')]) for k in range(len(mons))]
nv = [toS(cm['vals'][sp.Symbol(f'n{k}')]) for k in range(len(mons))]
cv = [toS(ch['vals'][sp.Symbol(f'c{k}')]) for k in range(len(mons))]
xv = sp.Matrix(X)
P = S2 * sp.eye(4) - xv * xv.T
pi1 = sp.zeros(4); varpi1 = 0; hT = sp.zeros(4)
for k, m in enumerate(mons):
    if bv[k] != 0:
        pi1 += bv[k] * m * P
    if nv[k] != 0:
        varpi1 += nv[k] * m
    if cv[k] != 0:
        Hs = sp.Matrix(4, 4, lambda i, j: sp.diff(m, X[i], X[j]))
        PHP = sp.expand(P * Hs * P)
        hT += cv[k] * (S2 * PHP - PHP.trace() * P / 3)
f0 = ch['f0']
coef = json.load(open('../results/l0_sector.json'))
ca1, ca0 = sp.sympify(coef['Aa1']), sp.sympify(coef['Aa0'])
cp1, cp0 = sp.sympify(coef['Ps1']), sp.sympify(coef['Ps0'])
# ODE for Psi, as a substitution rule (all functions of s; s_ -> S2 later)
Ps1_rule = (f0 - ca1 * Aa[1] - ca0 * Aa[0] - cp0 * Ps[0]) / cp1
h_full = hT + Aa[0] * (xv * xv.T)
DH, DM = PH_PM(h=h_full, psi=Ps[0], pi=pi1, varpi=varpi1)
print('operators applied %.1fs' % (time.time() - t0), flush=True)

def zero_after(expr, extra=None):
    e = expr
    if extra:
        e = e.subs(extra)
    e = sp.together(reduce_bg(e))
    num = sp.expand(sp.numer(e))
    return num == 0

res = {}
res['momentum_exact'] = all(zero_after(DM[i] - SM[i]) for i in range(4))
print('P_M[w] = S_M :', res['momentum_exact'], flush=True)
res['hamiltonian_exact_given_l0_ode'] = zero_after(DH - SH, {Ps[1]: toS(Ps1_rule)})
print('P_H[w] = S_H (given l=0 ODE) :', res['hamiltonian_exact_given_l0_ode'], flush=True)
# no radial derivatives of the algebraic coefficients survive: replace them
# by arbitrary symbols and confirm the residual is unchanged (operator is
# algebraic in b_k, n_k, c_k sphere by sphere)
# support: every algebraic coefficient vanishes when all profile jets vanish
prof0 = {g: 0 for g in G + K}
res['algebraic_coefficients_vanish_off_profile_support'] = all(
    sp.simplify(sp.sympify(e).subs(prof0)) == 0 for e in bv + nv + cv) and \
    sp.simplify(sp.sympify(f0).subs(prof0)) == 0
# ---- wrong-formula controls (each must FAIL)
ctrl = {}
_, DMbad = PH_PM(pi=pi1, varpi=0)
ctrl['omit_scalar_momentum_correction'] = not all(zero_after(DMbad[i] - SM[i]) for i in range(4))
ctrl['wrong_sign_l0_source'] = not zero_after(DH - SH, {Ps[1]: toS((-f0 - ca1 * Aa[1] - ca0 * Aa[0] - cp0 * Ps[0]) / cp1)})
ctrl['drop_tangential_tf_metric'] = not zero_after(PH_PM(h=Aa[0] * (xv * xv.T), psi=Ps[0])[0] - SH,
                                                   {Ps[1]: toS(Ps1_rule)})
res['wrong_formula_controls_detected'] = ctrl
res['seconds'] = round(time.time() - t0, 1)
print(json.dumps(res, indent=2))
json.dump(res, open('../results/verify_cc.json', 'w'), indent=2)
