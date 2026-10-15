"""Outward fixed-point intervals for the explicit smooth comparison waveform.

All decisions use integers/rationals. No floating-point transcendental library
is used. This module does not prove the retained analytic PDE hypotheses.
"""
from fractions import Fraction as F
from math import factorial, isqrt

BITS = 96
SCALE = 1 << BITS


def ceil_div(n, d):
    return -((-n) // d)


def interval(lo, hi=None):
    lo, hi = F(lo), F(lo if hi is None else hi)
    if lo > hi:
        raise ValueError('Reversed interval')
    return (lo.numerator * SCALE // lo.denominator,
            ceil_div(hi.numerator * SCALE, hi.denominator))


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def neg(a):
    return -a[1], -a[0]


def mul(a, b):
    products = [x * y for x in a for y in b]
    return min(products) // SCALE, ceil_div(max(products), SCALE)


def div_int(a, n):
    if not isinstance(n, int) or n <= 0:
        raise ValueError('Positive integer divisor required')
    return a[0] // n, ceil_div(a[1], n)


def square(a):
    low = 0 if a[0] <= 0 <= a[1] else min(a[0] ** 2, a[1] ** 2)
    high = max(a[0] ** 2, a[1] ** 2)
    return low // SCALE, ceil_div(high, SCALE)


def exp_negative(a):
    """Enclose exp(-a), rational a>=0, including an explicit large-a tail."""
    a = F(a)
    if a < 0:
        raise ValueError('Nonnegative exponential argument required')
    if a == 0:
        return SCALE, SCALE
    if a >= 128:
        # e>2, a>=floor(a); truncating the exponent only widens this bound.
        exponent = min(a.numerator // a.denominator, BITS)
        return 0, 1 << (BITS - exponent)
    count = 0
    reduced = a
    while reduced > F(1, 8):
        reduced /= 2
        count += 1
    x = neg(interval(reduced))
    term = total = (SCALE, SCALE)
    for k in range(1, 25):
        term = div_int(mul(term, x), k)
        total = add(total, term)
    # |R_24| <= 2*(1/8)^25/25!, since exp(1/8)<2.
    radius = ceil_div(2 * SCALE, 8 ** 25 * factorial(25))
    total = total[0] - radius, total[1] + radius
    for _ in range(count):
        total = square(total)
    return max(0, total[0]), min(SCALE, total[1])


def q_of_v(v):
    if v is None:
        return (0, 0)
    if v < 4:
        raise ValueError('Bump reciprocal must be at least four')
    return exp_negative(v - 4)


def derivative_envelope(v):
    # v^2 exp(4-v) decreases for v>=4; its limit at infinity is zero.
    if v is None:
        return (0, 0)
    return mul(interval(v * v), q_of_v(v))


def profile_cell(left, right):
    """Continuous ranges of Q and Q' on a closed cell in [0,3]."""
    left, right = F(left), F(right)
    if not 0 <= left <= right <= 3:
        raise ValueError('Cell must lie in [0,3]')
    for boundary in [F(1), F(2)]:
        if left < boundary < right:
            a, b = profile_cell(left, boundary), profile_cell(boundary, right)
            return tuple((min(x[0], y[0]), max(x[1], y[1])) for x, y in zip(a, b))
    if right <= 1 or left >= 2:
        return (0, 0), (0, 0)
    wleft, wright = (left - 1) * (2 - left), (right - 1) * (2 - right)
    wmin = min(wleft, wright)
    wmax = F(1, 4) if left <= F(3, 2) <= right else max(wleft, wright)
    vmin = 1 / wmax
    vmax = None if wmin == 0 else 1 / wmin
    qlo, qhi = q_of_v(vmax)[0], q_of_v(vmin)[1]
    elo, ehi = derivative_envelope(vmax)[0], derivative_envelope(vmin)[1]
    derivative = mul(interval(3 - 2 * right, 3 - 2 * left), (elo, ehi))
    return (qlo, qhi), derivative


def comparison_cell(left, right, ell=F(1, 10**10)):
    q, derivative = profile_cell(left, right)
    exponential = (exp_negative(ell * F(right))[0], exp_negative(ell * F(left))[1])
    return add(derivative, mul(interval(ell), mul(exponential, q)))


def sharp_waveform_error():
    """Sol's positive boundary-kernel/characteristic majorant, all rational."""
    ell, M, m, b, gamma = F(1, 10**10), F(10**6), F(1, 1000), F(1000), F(250)
    T = 3 * ell
    U, r, s = 1 + T, M * T ** 2 / 2, b * T ** 2 / (2 * m)
    determinant = (1 - r) * (1 - s) - T * r / m
    if r >= 1 or s >= 1 or determinant <= 0:
        raise ValueError('Positive-majorant conditions failed')
    A = U * (1 - s + T / m) / determinant
    V = U / (m * determinant)
    eta = ell * (M * T * A + (1 + gamma * T) * V)
    V_tail = ell / m ** 2 + ell / m + r * A / m + s * V
    eta_tail = ell * (M * T * A + V_tail + gamma * T * V)
    return {'ell': ell, 'T': T, 'r': r, 's': s, 'determinant': determinant,
            'A_upper': A, 'V_upper': V, 'eta_upper': eta,
            'published_eta': F(10000004, 10**14), 'eta_tail_upper': eta_tail,
            'published_tail': F(4752, 10**17), 'norm_error_upper': F(10000005, 10**14)}


def sqrt_bounds(value):
    value = F(value)
    if value < 0:
        raise ValueError('Nonnegative square root required')
    k = isqrt(value.numerator * SCALE * SCALE // value.denominator)
    low = F(k, SCALE)
    return low, low if low * low == value else F(k + 1, SCALE)
