"""Explicit second-order correction for the commuting state A = B = diag(1,-1,0,0).

Unknown polynomial coefficients are constants on each sphere |x|^2 = sigma and
are solved exactly by linear algebra.  Ansatz (all tangential, so radial
derivatives of the coefficients never appear):
  momentum   : pi'^{ij} = beta(x) (s delta_ij - x_i x_j),  varpi' = nu(x)
  hamiltonian: h = P Hess(chi) P - tr(P Hess(chi) P) P/(3s)  (tangential TF),
               psi = 0;  the l = 0 part of S_H is left to construct_l0.py.
beta, nu, chi are polynomials in x1^2, x2^2 (the symmetry of the source).
"""
import json, pickle, time
import sympy as sp
from engine import X, reduce_bg, S2
from sphere import sphere_average
from linear_ops import PH_PM

t0 = time.time()
d = pickle.load(open('../results/sources_cc.pkl', 'rb'))
SH, SM = d['SH'], d['SM']
sig = sp.Symbol('sigma')
X1, X2 = X[0]**2, X[1]**2
mons = [X1**a * X2**b for a in range(4) for b in range(4) if a + b <= 3]
bu = sp.symbols(f'b0:{len(mons)}'); nu_ = sp.symbols(f'n0:{len(mons)}')
cu = sp.symbols(f'c0:{len(mons)}')
beta = sum(c * m for c, m in zip(bu, mons))
nu = sum(c * m for c, m in zip(nu_, mons))
chi = sum(c * m for c, m in zip(cu, mons))
xv = sp.Matrix(X)
P = S2 * sp.eye(4) - xv * xv.T
pi1 = sp.expand(beta * P)
Hs = sp.Matrix(4, 4, lambda i, j: sp.diff(chi, X[i], X[j]))
PHP = sp.expand(P * Hs * P)
h1 = PHP - sp.expand(PHP.trace()) * P / (3 * S2)
DH, DM = PH_PM(h=h1, pi=pi1, varpi=nu)
print('operators applied %.1fs' % (time.time() - t0), flush=True)
f0 = sphere_average(SH)            # l = 0 part (function of s_ and jets)
s_ = sp.Symbol('s_')

def sphere_eqs(expr):
    e = sp.together(reduce_bg(expr))
    num, den = sp.fraction(e)
    num = sp.expand(num)
    # reduce modulo x4^2 = sigma - x1^2 - x2^2 - x3^2
    Pn = sp.Poly(num, X[3])
    red = 0
    for (k,), c in Pn.terms():
        red += c * (sig - X[0]**2 - X[1]**2 - X[2]**2)**(k // 2) * X[3]**(k % 2)
    red = sp.Poly(sp.expand(red), *X)
    return [sp.expand(c) for c in red.coeffs()]

eqs = []
for i in range(4):
    eqs += sphere_eqs(DM[i] - SM[i])
eqsM = len(eqs)
eqs += sphere_eqs(DH - (SH - f0.subs(s_, S2)))
print('equations: momentum %d, hamiltonian %d (%.1fs)' % (eqsM, len(eqs) - eqsM, time.time() - t0), flush=True)
unk = list(bu) + list(nu_) + list(cu)
sol = sp.linsolve(eqs, unk)
sol = list(sol)
print('solutions found:', len(sol), '%.1fs' % (time.time() - t0), flush=True)
if not sol:
    raise SystemExit('NO SOLUTION in ansatz')
sol = sol[0]
free = sorted(set().union(*[sp.sympify(v).free_symbols for v in sol]) & set(unk), key=str)
print('free parameters:', free, flush=True)
# fix free parameters to zero, convert sigma -> s
fix = {f: 0 for f in free}
vals = {u: sp.factor(sp.sympify(v).subs(fix).subs(sig, s_)) for u, v in zip(unk, sol)}
pickle.dump({'mons': mons, 'vals': vals, 'f0': f0}, open('../results/construct_cc.pkl', 'wb'))
nz = {str(k): str(v) for k, v in vals.items() if v != 0}
json.dump({'nonzero_coefficients': nz, 'free_parameters_set_to_zero': [str(f) for f in free],
           'f0': str(sp.factor(f0))}, open('../results/construct_cc.json', 'w'), indent=2)
print('nonzero coefficients:', len(nz), flush=True)
for k, v in nz.items():
    print(' ', k, '=', v[:200])
