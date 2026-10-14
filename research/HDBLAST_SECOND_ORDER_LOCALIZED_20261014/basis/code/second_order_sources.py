"""Second-order constraint sources S = -D^2 Phi[v,v] for quadrupole data.

Usage: py second_order_sources.py CASE
  CASE = cc   : A = B = diag(1,-1,0,0)           (commuting, the focused state)
         cd   : A = diag(1,-1,0,0), B = diag(0,0,1,-1)  (commuting, distinct)
         nc   : A = diag(1,-1,0,0), B = [A, M12] (non-commuting witness)
Data: q = G(s) x^T A x,  p = K(s) x^T B x  (G, K arbitrary radial profiles;
K = omega G recovers the frozen family).  Saves sympy srepr of S_H, S_M and
per-sphere angular moments.
"""
import json, sys, time, pickle
import sympy as sp
from engine import X, O, F, U, register, reduce_bg, S2
from fields import reconstruction, constraints, harmonic_quadratic

case = sys.argv[1]
A = sp.diag(1, -1, 0, 0)
if case == 'cc':
    B = A
elif case == 'cd':
    B = sp.diag(0, 0, 1, -1)
elif case == 'nc':
    M = sp.zeros(4); M[0, 1] = -1; M[1, 0] = 1
    B = A * M - M * A
else:
    raise SystemExit('unknown case')
G = register('G'); K = register('K')
q = G[0] * harmonic_quadratic(A)
p = K[0] * harmonic_quadratic(B)
v = reconstruction(q, p)
t0 = time.time()
Hm, Mm = constraints(v['h'], v['dphi'], v['dpi'], v['dvarpi'])
print('jets built %.1fs' % (time.time() - t0), flush=True)

def canon(e):
    e = sp.together(reduce_bg(e))
    n, dd = sp.fraction(e)
    return sp.expand(n) / sp.factor(dd)

SH = canon(-Hm.c)
print('S_H done %.1fs' % (time.time() - t0), flush=True)
SM = [canon(-Mm[i].c) for i in range(4)]
print('S_M done %.1fs' % (time.time() - t0), flush=True)
# first-order residuals as a sanity check
fo = [sp.expand(sp.numer(sp.together(reduce_bg(Hm.b)))) == 0] + \
     [sp.expand(sp.numer(sp.together(reduce_bg(Mm[i].b)))) == 0 for i in range(4)]
print('first-order constraints vanish:', fo, flush=True)
pickle.dump({'SH': SH, 'SM': SM, 'A': A, 'B': B, 'first_order_ok': all(fo)},
            open(f'../results/sources_{case}.pkl', 'wb'))
print('saved', flush=True)
