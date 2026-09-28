#!/usr/bin/env python3
"""Containment self-test of the rigorous arithmetic core against exact rational arithmetic (fractions).
Random Taylor models / intervals are combined with +, -, *, scale, addc, recip, Cauchy products; at random rational
parameter points u in [-1,1] and random selections inside the remainders the exact result must lie in the enclosure."""
import random
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import IV, TM, ONE
random.seed(20260916)
cc.set_degree(4)
D = 4
def rnd_tm(scale=1.0, spread=1e-3):
    c = [int(random.uniform(-1, 1) * scale * ONE)] + [int(random.uniform(-1, 1) * scale * spread ** k * ONE) for k in range(1, D + 1)]
    r = sorted([int(random.uniform(-1, 1) * 1e-20 * ONE), int(random.uniform(-1, 1) * 1e-20 * ONE)])
    return TM(c, r[0], r[1])
def sample(tm, u):
    """exact value of one member function of the TM at u (random point of the remainder)"""
    e = Fr(random.randint(tm.rlo, tm.rhi), ONE)
    return sum(Fr(c, ONE) * u ** k for k, c in enumerate(tm.c)) + e
def inside(tm, u, val):
    p = sum(Fr(c, ONE) * u ** k for k, c in enumerate(tm.c))
    return Fr(tm.rlo, ONE) <= val - p <= Fr(tm.rhi, ONE)
n = 0
for trial in range(300):
    a, b, c = rnd_tm(), rnd_tm(3.0), rnd_tm(0.5)
    if abs(b.c[0]) < ONE // 10:
        continue
    for _ in range(3):
        u = Fr(random.randint(-1000, 1000), 1000)
        va, vb, vc = sample(a, u), sample(b, u), sample(c, u)
        assert inside(a + b, u, va + vb) and inside(a - b, u, va - vb)
        assert inside(a * b, u, va * vb)
        assert inside(a.scale(-7, 3), u, va * Fr(-7, 3)) and inside(a.addc(Fr(5, 7)), u, va + Fr(5, 7))
        assert inside(b.recip(), u, 1 / vb)
        assert inside(TM.cauchy([a, b, c], [c, a, b], 2), u, va * vb + vb * va + vc * vc)
        bb = (a * b + c).bound()
        assert Fr(bb.lo, ONE) <= va * vb + vc <= Fr(bb.hi, ONE)
        n += 8
for trial in range(2000):
    x = sorted(random.uniform(-3, 3) for _ in range(2)); y = sorted(random.uniform(0.1, 4) for _ in range(2))
    X = IV.frac(Fr(x[0]), Fr(x[1])); Y = IV.frac(Fr(y[0]), Fr(y[1]))
    vx = Fr(x[0]) + (Fr(x[1]) - Fr(x[0])) * Fr(random.randint(0, 100), 100); vy = Fr(y[0]) + (Fr(y[1]) - Fr(y[0])) * Fr(random.randint(0, 100), 100)
    for iv, v in ((X * Y, vx * vy), (X - Y, vx - vy), (Y.recip(), 1 / vy), ((-Y).recip(), -1 / vy), (X.scale(-5, 9), vx * Fr(-5, 9)),
                  (IV.cauchy([X, Y], [Y, X], 1), vx * vx + vy * vy)):
        assert Fr(iv.lo, ONE) <= v <= Fr(iv.hi, ONE)
        n += 1
print("selftest_core: %d exact containment checks passed" % n)
