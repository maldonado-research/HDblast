"""Exact sparse Laurent differential algebra; standard library only."""
from fractions import Fraction as F


class P:
    def __init__(self, terms=0):
        if isinstance(terms, P):
            self.terms = dict(terms.terms)
        elif isinstance(terms, dict):
            self.terms = {k: F(v) for k, v in terms.items() if v}
        else:
            self.terms = {(): F(terms)} if terms else {}

    @classmethod
    def var(cls, name):
        return cls({((name, 1),): F(1)})

    @classmethod
    def cv(cls, value):
        return value if isinstance(value, cls) else cls(value)

    def __add__(self, other):
        terms = dict(self.terms)
        for monomial, value in self.cv(other).terms.items():
            terms[monomial] = terms.get(monomial, F(0)) + value
        return P(terms)

    __radd__ = __add__

    def __neg__(self):
        return P({m: -v for m, v in self.terms.items()})

    def __sub__(self, other):
        return self + -self.cv(other)

    def __rsub__(self, other):
        return self.cv(other) + -self

    def __mul__(self, other):
        terms = {}
        for m, value in self.terms.items():
            for n, factor in self.cv(other).terms.items():
                powers = dict(m)
                for name, exponent in n:
                    powers[name] = powers.get(name, 0) + exponent
                key = tuple(sorted((name, e) for name, e in powers.items() if e))
                terms[key] = terms.get(key, F(0)) + value * factor
        return P(terms)

    __rmul__ = __mul__

    def __pow__(self, power):
        if not isinstance(power, int):
            raise ValueError('integer exponent required')
        if power < 0:
            if len(self.terms) != 1:
                raise ValueError('nonmonomial inverse')
            (m, value), = self.terms.items()
            return P({tuple((name, e * power) for name, e in m): value ** power})
        result = P(1)
        factor = self
        while power:
            if power & 1:
                result = result * factor
            factor = factor * factor
            power //= 2
        return result

    def __truediv__(self, other):
        return self * self.cv(other) ** -1

    def __rtruediv__(self, other):
        return self.cv(other) * self ** -1

    def __bool__(self):
        return bool(self.terms)

    def __eq__(self, other):
        return self.terms == self.cv(other).terms

    def diff(self, name):
        terms = {}
        for m, value in self.terms.items():
            powers = dict(m)
            exponent = powers.get(name, 0)
            if exponent:
                powers[name] = exponent - 1
                key = tuple(sorted((s, e) for s, e in powers.items() if e))
                terms[key] = value * exponent
        return P(terms)

    def names(self):
        return sorted({name for m in self.terms for name, _ in m})

    def encode(self):
        return [[[[name, e] for name, e in m], str(value)]
                for m, value in sorted(self.terms.items())]

    def substitute(self, values):
        result = P()
        for m, value in self.terms.items():
            term = P(value)
            for name, exponent in m:
                term = term * values.get(name, P.var(name)) ** exponent
            result = result + term
        return result

    def dr(self):
        return master_derivative(self, 'r')

    def dx(self):
        return master_derivative(self, 'x')


def derivation(poly, rules):
    poly = P.cv(poly)
    return sum((poly.diff(name) * rules.get(name, P()) for name in poly.names()), P())


def master_derivative(poly, axis):
    rules = {}
    for name in P.cv(poly).names():
        if name.startswith('W'):
            suffix = name[1:] + axis
            rules[name] = P.var('W' + ''.join(sorted(suffix)))
    return derivation(poly, rules)
