"""Exact S^3 moments and harmonic projections of polynomial fields.

Fields are rational in background/profile jets (constant on each sphere
|x|^2 = s) and polynomial in x.  For a monomial with even exponents 2a_i,
    <n^alpha> = Gamma(2) prod Gamma(a_i+1/2) / (pi^2 Gamma(2+sum a_i)),
the average over the unit S^3 (odd exponents give zero).
"""
import sympy as sp
from engine import X

s = sp.Symbol('s_')

def sphere_average(expr):
    """<expr> over |x|^2 = s, returned as expression in s_ and jets."""
    num, den = sp.fraction(sp.together(expr))
    if any(xx in den.free_symbols for xx in X):
        raise ValueError('denominator depends on x')
    P = sp.Poly(sp.expand(num), *X)
    acc = 0
    for mon, coeff in P.terms():
        if any(m % 2 for m in mon):
            continue
        a = [m // 2 for m in mon]
        avg = sp.gamma(2) * sp.prod([sp.gamma(ai + sp.Rational(1, 2)) for ai in a]) / (
            sp.pi**2 * sp.gamma(2 + sum(a)))
        acc += coeff * sp.simplify(avg) * s**sum(a)
    return sp.together(acc / den)

def parity(expr):
    """Return +1 if even, -1 if odd, 0 if mixed under x -> -x."""
    e = sp.expand(sp.numer(sp.together(expr)))
    m = sp.expand(e.subs({xx: -xx for xx in X}, simultaneous=True))
    if sp.expand(m - e) == 0:
        return 1
    if sp.expand(m + e) == 0:
        return -1
    return 0
