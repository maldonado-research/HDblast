"""Fast exact verification of the explicit second-order correction
(A = B = diag(1,-1,0,0)), alternative to verify_cc.py.

Step 1: apply the linearized operators to w with the radial coefficient
functions kept as FREE radial functions (fresh jets Bk, Nk, Ck); prove that no
s-derivative of any coefficient survives (so the sphere-by-sphere linear
algebra of construct_cc2 is an exact pointwise identity).
Step 2: substitute the solved coefficients and check P_M[w] = S_M as an exact
symbolic identity (after imposing the background equations).
The Hamiltonian identity P_H[h_T] = S_H - f0 is checked in
exact_random_point.py (exact rational arithmetic); a full symbolic expansion
of it was attempted here but did not finish in ~45 CPU minutes.
"""
import json, pickle, sys, os, time
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'basis', 'code'))
import sympy as sp
from engine import X, O, F, register, reduce_bg, S2
from linear_ops import PH_PM
from sphere import sphere_average
t0 = time.time()
res_dir = os.environ.get('SO_RESULTS', os.path.join(here, '..', 'basis', 'results'))
src = pickle.load(open(os.path.join(res_dir, 'sources_cc.pkl'), 'rb'))
SH, SM = src['SH'], src['SM']
cm = json.load(open(os.path.join(res_dir, 'construct_cc_M.json')))['nonzero_coefficients']
ch = json.load(open(os.path.join(res_dir, 'construct_cc_H.json')))['nonzero_coefficients']
s_ = sp.Symbol('s_')
X1, X2 = X[0]**2, X[1]**2
mons = [X1**a * X2**b for a in range(4) for b in range(4) if a + b <= 3]
G = register('G'); K = register('K')
xv = sp.Matrix(X)
P = S2 * sp.eye(4) - xv * xv.T
fresh = {}
pi2 = sp.zeros(4); varpi2 = 0; hT = sp.zeros(4)
for k, m in enumerate(mons):
    if f'b{k}' in cm:
        J = register(f'Bf{k}'); fresh[f'b{k}'] = J
        pi2 += J[0] * m * P
    if f'n{k}' in cm:
        J = register(f'Nf{k}'); fresh[f'n{k}'] = J
        varpi2 += J[0] * m
    if f'c{k}' in ch:
        J = register(f'Cf{k}'); fresh[f'c{k}'] = J
        Hs = sp.Matrix(4, 4, lambda i, j: sp.diff(m, X[i], X[j]))
        PHP = sp.expand(P * Hs * P)
        hT += J[0] * sp.expand(S2 * PHP - PHP.trace() * P / 3)
DH, _ = PH_PM(h=hT)
_, DM = PH_PM(pi=pi2, varpi=varpi2)
print('operators applied %.1fs' % (time.time() - t0), flush=True)
res = {}
# Step 1: derivative jets must drop out
deriv_ok = True
for key, J in fresh.items():
    for e in [DH] + DM:
        for jd in J[1:4]:
            c = sp.diff(e, jd)
            if c != 0:
                c = sp.expand(sp.numer(sp.together(reduce_bg(c))))
                if c != 0:
                    deriv_ok = False
                    print('derivative of', key, 'survives', flush=True)
res['no_radial_derivatives_of_coefficients'] = deriv_ok
print('step1 no radial derivatives:', deriv_ok, '%.1fs' % (time.time() - t0), flush=True)
# Step 2: substitute solved values
vals = {}
for key, J in fresh.items():
    expr = sp.sympify(cm.get(key, ch.get(key)))
    vals[J[0]] = expr.subs(s_, S2)
f0 = sphere_average(SH)
def zero(expr):
    e = sp.together(reduce_bg(expr))
    return sp.expand(sp.numer(e)) == 0
res['momentum_exact'] = all(zero(DM[i].subs(vals) - SM[i]) for i in range(4))
print('P_M[w] = S_M :', res['momentum_exact'], '%.1fs' % (time.time() - t0), flush=True)
# Step 3: wrong-formula controls
ctrl = {}
vals_bad = dict(vals); kk = list(fresh['n2'])[0]; vals_bad[kk] = 2 * vals[kk]
ctrl['doubled_scalar_momentum_coefficient'] = not all(zero(DM[i].subs(vals_bad) - SM[i]) for i in range(4))
res['wrong_formula_controls_detected'] = ctrl
res['seconds'] = round(time.time() - t0, 1)
print(json.dumps(res, indent=2))
json.dump(res, open(os.path.join(here, 'fast_verify.json'), 'w'), indent=2)
