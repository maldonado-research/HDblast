"""Finite exact-rational identity test of the Hamiltonian and
momentum identities for the explicit correction, A = B = diag(1,-1,0,0).

After the background rules are imposed, every quantity is a rational function
of independent symbols (x1..x4 and the free background/profile jets).  The
identities P_H[h_T] = S_H - f0 and P_M[w] = S_M are evaluated in EXACT
rational arithmetic at three seeded rational points. These finite tests are
not symbolic proofs, and no probability bound is claimed without a degree
bound and an analysis of the sampling distribution. Wrong-formula
controls must give nonzero values.  Step 1 (no radial derivatives of the
coefficients) is the separate symbolic check in fast_verify.py."""
import json, pickle, sys, os, time, random
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, '..', 'basis', 'code'))
import sympy as sp
from engine import X, O, F, register, reduce_bg, S2
from linear_ops import PH_PM
from sphere import sphere_average
t0 = time.time()
res_dir = os.path.join(here, '..', 'basis', 'results')
src = pickle.load(open(os.path.join(res_dir, 'sources_cc.pkl'), 'rb'))
SH, SM = src['SH'], src['SM']
cm = json.load(open(os.path.join(res_dir, 'construct_cc_M.json')))['nonzero_coefficients']
ch = json.load(open(os.path.join(res_dir, 'construct_cc_H.json')))['nonzero_coefficients']
s_ = sp.Symbol('s_')
X1, X2 = X[0]**2, X[1]**2
mons = [X1**a * X2**b for a in range(4) for b in range(4) if a + b <= 3]
G = register('G'); K = register('K')
xv = sp.Matrix(X); P = S2 * sp.eye(4) - xv * xv.T
fresh = {}; pi2 = sp.zeros(4); varpi2 = 0; hT = sp.zeros(4)
for k, m in enumerate(mons):
    if f'b{k}' in cm:
        J = register(f'Bf{k}'); fresh[f'b{k}'] = J; pi2 += J[0] * m * P
    if f'n{k}' in cm:
        J = register(f'Nf{k}'); fresh[f'n{k}'] = J; varpi2 += J[0] * m
    if f'c{k}' in ch:
        J = register(f'Cf{k}'); fresh[f'c{k}'] = J
        Hs = sp.Matrix(4, 4, lambda i, j: sp.diff(m, X[i], X[j]))
        PHP = sp.expand(P * Hs * P)
        hT += J[0] * sp.expand(S2 * PHP - PHP.trace() * P / 3)
DH, _ = PH_PM(h=hT)
_, DM = PH_PM(pi=pi2, varpi=varpi2)
# impose background rules on operators and sources (s_ -> |x|^2)
DHr = reduce_bg(DH); DMr = [reduce_bg(e) for e in DM]
SHr = reduce_bg(SH); SMr = [reduce_bg(e) for e in SM]
f0 = reduce_bg(sphere_average(SH).subs(s_, S2))
vals = {J[0]: reduce_bg(sp.sympify(cm.get(k, ch.get(k))).subs(s_, S2)) for k, J in fresh.items()}
print('prepared %.1fs' % (time.time() - t0), flush=True)
# derivative jets of the coefficients must be absent here (fast_verify step 1);
# set them to random values anyway: any dependence would show up as nonzero
syms = set()
for e in [DHr, SHr, f0] + DMr + SMr + list(vals.values()):
    syms |= e.free_symbols
coef_syms = set(J[0] for J in fresh.values())
deriv_syms = set(j for J in fresh.values() for j in J[1:])
free = sorted(syms - coef_syms, key=str)
rng = random.Random(20261014)
def rand_point():
    pt = {}
    for v in free:
        pt[v] = sp.Rational(rng.randint(-97, 97) or 1, rng.randint(1, 61))
    return pt
def ev(e, pt):
    return sp.nsimplify(e.subs(pt)) if False else e.xreplace(pt)
results = []
for trial in range(3):
    pt = rand_point()
    cv = {c: sp.simplify(vals[c].xreplace(pt)) for c in coef_syms}
    full = dict(pt); full.update(cv)
    rH = sp.simplify(DHr.xreplace(full) - (SHr.xreplace(pt) - f0.xreplace(pt)))
    rM = [sp.simplify(DMr[i].xreplace(full) - SMr[i].xreplace(pt)) for i in range(4)]
    scaleH = sp.Abs(SHr.xreplace(pt))
    # controls
    cbad = dict(full); kk = fresh['c1'][0]; cbad[kk] = -cv[kk]
    ctrlH = sp.simplify(DHr.xreplace(cbad) - (SHr.xreplace(pt) - f0.xreplace(pt))) != 0
    nbad = dict(full); kk = fresh['n2'][0]; nbad[kk] = 2 * cv[kk]
    ctrlM = any(sp.simplify(DMr[i].xreplace(nbad) - SMr[i].xreplace(pt)) != 0 for i in range(4))
    ctrl0 = sp.simplify(DHr.xreplace(full) - SHr.xreplace(pt)) != 0
    results.append({'H_residual_exact': str(rH), 'M_residual_exact': [str(r) for r in rM],
                    'S_H_value': str(sp.nsimplify(scaleH)) if scaleH.is_Rational else str(sp.N(scaleH, 15)),
                    'control_sign_flipped_tf_metric_detected': bool(ctrlH),
                    'control_doubled_scalar_momentum_detected': bool(ctrlM),
                    'control_omit_l0_subtraction_detected': bool(ctrl0)})
    print(trial, results[-1]['H_residual_exact'], results[-1]['M_residual_exact'],
          ctrlH, ctrlM, ctrl0, '%.1fs' % (time.time() - t0), flush=True)
out = {'method': 'finite exact rational evaluation at three seeded sample points',
       'free_symbols': [str(v) for v in free], 'trials': results,
       'all_identities_zero': all(r['H_residual_exact'] == '0' and all(m == '0' for m in r['M_residual_exact']) for r in results),
       'all_controls_detected': all(r['control_sign_flipped_tf_metric_detected'] and r['control_doubled_scalar_momentum_detected']
                                    and r['control_omit_l0_subtraction_detected'] for r in results),
       'seconds': round(time.time() - t0, 1)}
json.dump(out, open(os.path.join(here, 'exact_random_point.json'), 'w'), indent=2)
print(json.dumps({k: v for k, v in out.items() if k != 'trials'}, indent=2))
