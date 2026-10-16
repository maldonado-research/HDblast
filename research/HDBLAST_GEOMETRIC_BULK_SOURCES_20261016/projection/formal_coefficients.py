"""Portable sparse Laurent ring for independent ordinary y/t jet symbols.

Monomials use sorted (name, integer exponent) pairs. No sampled substitution,
floating arithmetic, CAS, or assertion is used in equality checks.
"""
from fractions import Fraction as Q
from functools import lru_cache


@lru_cache(maxsize=131072)
def product_monomial(left, right):
    out = dict(left)
    for name, power in right:
        out[name] = out.get(name, 0) + power
        if not out[name]: del out[name]
    return tuple(sorted(out.items()))


class P:
    __slots__ = ('terms',)

    def __init__(self, value=0):
        if isinstance(value, P): self.terms = value.terms
        elif isinstance(value, dict):
            self.terms = {tuple(e): Q(v) for e, v in value.items() if v}
        else:
            value = Q(value)
            self.terms = {(): value} if value else {}

    def __add__(self, other):
        other = P(other)
        if not self: return other
        if not other: return self
        out = self.terms.copy()
        for e, value in other.terms.items():
            out[e] = out.get(e, Q(0)) + value
            if not out[e]: del out[e]
        return P(out)

    __radd__ = __add__
    def __neg__(self): return P({e: -v for e, v in self.terms.items()})
    def __sub__(self, other): return self + -P(other)
    def __rsub__(self, other): return P(other) + -self

    def __mul__(self, other):
        other = P(other)
        if not self or not other: return P()
        out = {}
        for e, value in self.terms.items():
            for f, weight in other.terms.items():
                g = product_monomial(e, f)
                out[g] = out.get(g, Q(0)) + value * weight
        return P(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = P(other)
        if len(other.terms) != 1:
            raise ValueError('FORMAL_DIVISOR_MUST_BE_NONZERO_MONOMIAL')
        monomial, weight = next(iter(other.terms.items()))
        inverse = tuple((name, -power) for name, power in monomial)
        return P({product_monomial(e, inverse): v / weight
                  for e, v in self.terms.items()})

    def __rtruediv__(self, other): return P(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError('FORMAL_EXPONENT_MUST_BE_INTEGER')
        if exponent < 0: return P(1) / (self ** -exponent)
        result, base = P(1), self
        while exponent:
            if exponent & 1: result *= base
            exponent //= 2
            if exponent: base *= base
        return result

    def __eq__(self, other):
        try: return self.terms == P(other).terms
        except (TypeError, ValueError): return False

    def __bool__(self): return bool(self.terms)
    def __repr__(self): return 'P(' + repr(self.terms) + ')'

    def diff(self, name):
        out = {}
        for monomial, value in self.terms.items():
            powers = dict(monomial)
            power = powers.get(name, 0)
            if power:
                powers[name] -= 1
                if not powers[name]: del powers[name]
                out[tuple(sorted(powers.items()))] = value * power
        return P(out)

    def evaluate(self, values):
        out = 0
        for monomial, value in self.terms.items():
            term = value
            for name, power in monomial:
                base = values[name]
                term *= (Q(base) if isinstance(base, int) else base) ** power
            out += term
        return out

    def encode(self):
        return [{'coefficient': str(value), 'powers': dict(monomial)}
                for monomial, value in sorted(self.terms.items())]


def symbol(name):
    if not isinstance(name, str) or not name: raise ValueError('FORMAL_SYMBOL_NAME')
    return P({((name, 1),): Q(1)})


def formal_inputs():
    """Independent off-shell ordinary derivatives through total order two.

    b has a=rho_y/rho, ay=a_y, theta=tanh(t), c=cosh(t), p=phi_y,
    py=phi_yy, Vp=V'(phi), Vpp=V''(phi), Vppp=V'''(phi).
    First harmonic constraints are deliberately NOT imposed here.
    """
    m = {name + suffix: symbol(name + suffix)
         for name in ('A', 'B', 'C', 'D', 'E', 'L', 'M')
         for suffix in ('', 'y', 't', 'yy', 'yt', 'tt')}
    s = {name: symbol(name) for name in ('s', 'sy', 'st', 'syy', 'syt', 'stt')}
    b = {name: symbol(name) for name in
         ('rho', 'c', 'theta', 'a', 'ay', 'p', 'py', 'Vp', 'Vpp', 'Vppp', 'ell')}
    return m, s, b
