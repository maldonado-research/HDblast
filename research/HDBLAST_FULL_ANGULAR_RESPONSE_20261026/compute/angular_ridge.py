#!/usr/bin/env python3
"""Exact first reflected feedback for the complete 23-amplitude angular space.

All returned numbers are fractions. This evaluates a specified approximation,
not the actual HDBLAST field or its error without a separate native ledger.
"""
from fractions import Fraction as Q
from math import factorial

LABELS = tuple(f'{field}:{degree}' for field in ('K','L','N','J','w') for degree in (0,2,4)) + tuple(f'{field}:{degree}' for field in ('v','u') for degree in (2,4)) + ('H:2','H:4','T:+','T:-')
DIM = len(LABELS)
PARITY = tuple(-1 if name.startswith(('L:','v:')) else 1 for name in LABELS)
MEANS = LABELS

def angular_gram():
    scalar_norms = {0: Q(1), 2: Q(1,12), 4: Q(19,720)}
    weights = []
    for name in LABELS:
        if name.startswith('T:'):
            value = Q(1,6)
        else:
            field, degree = name.split(':')
            degree = int(degree)
            value = scalar_norms[degree]
            if field == 'L': value *= 2
            if field == 'J': value *= 3
            if field in ('v','u'): value *= 2*degree*(degree+2)
            if field == 'H': value *= Q(2,3)*degree*(degree+2)*(degree*(degree+2)-3)
        weights.append(value)
    return tuple(weights)

GRAM = angular_gram()


def need(ok, message):
    if not ok:
        raise ValueError(message)


def vector(values):
    out = tuple(map(Q, values))
    need(len(out) == DIM, "full23 amplitudes required")
    return out


ZERO = vector([0] * DIM)


def va(a, b):
    return tuple(x + y for x, y in zip(a, b))


def vs(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def matrix(rows):
    out = tuple(vector(row) for row in rows)
    need(len(out) == DIM, "full23 by23 matrix required")
    return out


def mv(a, v):
    return tuple(sum((aij * vj for aij, vj in zip(row, v)), Q(0)) for row in a)


def ma(a, b):
    return tuple(va(x, y) for x, y in zip(a, b))


def ms(a, b):
    return tuple(vs(x, y) for x, y in zip(a, b))


class PiecewiseLinear:
    """One continuous rational profile, with a specified reflection.

    primitive(t,k) is the k-fold primitive normalized at 0. Negative arguments
    use (-1)^k times the component reflection. All joins use shared node values.
    """

    def __init__(self, nodes, values, parity=PARITY):
        self.nodes = tuple(map(Q, nodes))
        self.values = tuple(vector(v) for v in values)
        self.parity = tuple(parity)
        need(len(self.nodes) == len(self.values) >= 2, "profile shape")
        need(self.nodes[0] == 0, "profile starts at zero")
        need(all(a < b for a, b in zip(self.nodes, self.nodes[1:])), "ordered nodes")
        need(len(self.parity) == DIM and set(self.parity) <= {-1, 1}, "reflection signs")
        need(all(s == 1 or v == 0 for s, v in zip(self.parity, self.values[0])),
             "odd profile must be zero at the origin")
        slopes = [scale(1 / (b - a), vs(vb, va0))
                  for a, b, va0, vb in
                  zip(self.nodes, self.nodes[1:], self.values, self.values[1:])]
        self.initial_slope = slopes[0]
        self.hinges = tuple((self.nodes[i], vs(slopes[i], slopes[i-1]))
                            for i in range(1, len(slopes)))

    def primitive(self, t, order=0):
        t = Q(t)
        need(isinstance(order, int) and 0 <= order <= 5, "primitive order 0..5")
        need(abs(t) <= self.nodes[-1], "outside registered profile footprint")
        negative = t < 0
        t = abs(t)
        out = va(scale(t**order / factorial(order), self.values[0]),
                 scale(t**(order+1) / factorial(order+1), self.initial_slope))
        for knot, change in self.hinges:
            if t > knot:
                out = va(out, scale((t-knot)**(order+1) / factorial(order+1), change))
        if negative:
            out = tuple(((-1)**order) * sign * v for sign, v in zip(self.parity, out))
        return out


def reflected_kernel(F, G, x, r, parity=PARITY):
    """Return Kf, (Kf)_x and (Kf)_r for physical f(s,v)=F(v+s)+G(v-s).

    F(t,k), G(t,k) return a vector profile or its normalized primitive k=0,1,2.
    The complete physical source is reflected with parity, including couplings.
    """
    original_F, original_G = F, G
    F = lambda y, k=0: vector(original_F(y, k))
    G = lambda y, k=0: vector(original_G(y, k))
    need(all(F(0, k) == G(0, k) == ZERO for k in (1, 2)),
         "primitives must be normalized at zero")
    x, r = Q(x), Q(r)
    need(x >= 0 and r >= 0, "physical half plane")
    need(len(parity) == DIM and set(parity) <= {-1, 1}, "reflection signs")
    a, q = r+x, r-x
    pf_a, qf_a, pg_a, qg_a = F(a, 1), F(a, 2), G(a, 1), G(a, 2)
    p = va(scale(x, F(a, 0)), scale(Q(1, 2), vs(pg_a, G(q, 1))))
    if q >= 0:
        value = scale(Q(1, 2), vs(va(vs(scale(x, pf_a),
                          scale(Q(1, 2), vs(qf_a, F(q, 2)))),
                          scale(Q(1, 2), vs(qg_a, G(q, 2)))), scale(x, G(q, 1))))
        m = va(scale(Q(1, 2), vs(pf_a, F(q, 1))), scale(x, G(q, 0)))
    else:
        t = -q
        physical = vs(va(vs(scale(x, pf_a),
                            scale(Q(1, 2), va(qf_a, F(t, 2)))),
                        scale(Q(1, 2), va(qg_a, G(-t, 2)))), scale(r, G(-t, 1)))
        reflected = va(vs(scale(t, F(t, 1)), F(t, 2)),
                       scale(Q(1, 2), va(G(t, 2), G(-t, 2))))
        value = scale(Q(1, 2), va(physical,
                         tuple(sign*v for sign, v in zip(parity, reflected))))
        reflected_m = va(scale(t, F(t, 0)),
                         scale(Q(1, 2), vs(G(t, 1), G(-t, 1))))
        m = va(va(scale(Q(1, 2), vs(pf_a, F(t, 1))), scale(r, G(-t, 0))),
               tuple(sign*v for sign, v in zip(parity, reflected_m)))
    return {"value": value, "dx": scale(Q(1, 2), va(p, m)),
            "dr": scale(Q(1, 2), vs(p, m))}


class CoupledRidges:
    """All 23 amplitudes, with independent forcing approximants u and du.

    du is not claimed to be the derivative of the piecewise-linear u. These
    two approximants appear only inside the forcing. The native certificate
    compares each with its own registered rational initial function.
    """

    def __init__(self, A, B, D, u, du, v):
        self.A, self.B, self.D = matrix(A), matrix(B), matrix(D)
        self.u, self.du, self.v = u, du, v
        need(u.parity == PARITY and v.parity == PARITY, "initial value/velocity parity")
        need(du.parity == tuple(-s for s in PARITY), "radial derivative parity")
        need(u.nodes == du.nodes == v.nodes, "one declared profile grid")
        need(u.values[0] == v.values[0] == ZERO, "zero initial collar at origin")
        self.plus, self.minus = ma(self.A, self.B), ms(self.B, self.A)

    def F(self, y, order=0):
        return scale(Q(1, 2),
             va(mv(self.plus, va(self.du.primitive(y, order),
                                self.v.primitive(y, order))),
                mv(self.D, va(self.u.primitive(y, order),
                              self.v.primitive(y, order+1)))))

    def G(self, y, order=0):
        return scale(Q(1, 2),
             va(va(mv(self.minus, self.du.primitive(y, order)),
                   scale(-1, mv(self.minus, self.v.primitive(y, order)))),
                mv(self.D, vs(self.u.primitive(y, order),
                              self.v.primitive(y, order+1)))))

    def physical_source(self, x, r):
        x, r = Q(x), Q(r)
        need(x >= 0 and r >= 0, "physical source domain")
        return va(self.F(r+x), self.G(r-x))

    def feedback(self, x, r, M):
        x, r, M = Q(x), Q(r), matrix(M)
        need(0 <= x <= Q(3, 2) and r >= max(Q(0), x-Q(5, 4))
             and x+r <= Q(7, 4), "outside retained causal cone")
        need(all(M[i][j] == 0 for i, sign in enumerate(PARITY) if sign == -1 for j in range(DIM)),
             "shell lift has no Dirichlet forcing in any angular row")
        out = reflected_kernel(self.F, self.G, x, r)
        t = x-r
        if t > 0:
            project = lambda row: tuple(v if s == 1 else Q(0)
                                        for s, v in zip(PARITY, row))
            g = mv(M, project(va(self.u.primitive(t), self.v.primitive(t, 1))))
            pg = mv(M, project(va(self.u.primitive(t, 1), self.v.primitive(t, 2))))
            out = {"value": vs(out["value"], pg), "dx": vs(out["dx"], g),
                   "dr": va(out["dr"], g)}
        return out


def from_data(data):
    p = data["profiles"]
    u = PiecewiseLinear(data["nodes"], p["u"])
    du = PiecewiseLinear(data["nodes"], p["du"], tuple(-s for s in PARITY))
    v = PiecewiseLinear(data["nodes"], p["v"])
    m = data["matrices"]
    return CoupledRidges(m["A"], m["B"], m["D"], u, du, v)


def encoded(out):
    return {key: [str(v) for v in row] for key, row in out.items()}
