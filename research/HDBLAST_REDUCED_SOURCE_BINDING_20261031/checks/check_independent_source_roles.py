#!/usr/bin/env python3
"""Independent exact map/frame joins. Standard library only; no source imports.

The physical first field is the registered seven-profile l=2 ansatz, with free
ordinary y,t derivatives. Full quotient-polynomial components are compared,
not angular point samples. This checker does not derive the PDE existence
theorem or claim generic spacetime metric completeness beyond that ansatz.
"""
import argparse
import ast
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
from itertools import product, permutations
import json
from pathlib import Path
from math import comb
import sys


def need(ok, why):
    if not ok:
        raise ValueError(why)


@lru_cache(maxsize=262144)
def monomial_product(a, b):
    powers = dict(a)
    for name, power in b:
        powers[name] = powers.get(name, 0) + power
    return tuple(sorted((n, e) for n, e in powers.items() if e))


class P:
    def __init__(self, value=0):
        if isinstance(value, P):
            self.d = value.d
        elif isinstance(value, dict):
            self.d = {tuple(k): Q(v) for k, v in value.items() if v}
        else:
            self.d = {(): Q(value)} if value else {}

    @staticmethod
    def var(name):
        return P({((name, 1),): 1})

    def __bool__(self):
        return bool(self.d)

    def __add__(self, other):
        if not isinstance(other, (P, int, Q)):
            return NotImplemented
        other = P(other)
        if not self:
            return other
        if not other:
            return self
        d = self.d.copy()
        for k, v in other.d.items():
            d[k] = d.get(k, 0) + v
            if not d[k]:
                del d[k]
        return P(d)

    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -P(other)

    def __rsub__(self, other):
        return P(other) + -self

    def __mul__(self, other):
        if not isinstance(other, (P, int, Q)):
            return NotImplemented
        other = P(other)
        d = {}
        for k, v in self.d.items():
            for j, w in other.d.items():
                h = monomial_product(k, j)
                d[h] = d.get(h, 0) + v * w
        return P(d)

    __rmul__ = __mul__

    def __pow__(self, n):
        need(isinstance(n, int), 'POWER_TYPE')
        if n < 0:
            need(len(self.d) == 1, 'NONMONOMIAL_INVERSE')
            (k, v), = self.d.items()
            return P({tuple((s, e * n) for s, e in k): v ** n})
        out, base = P(1), self
        while n:
            if n & 1:
                out = out * base
            n //= 2
            if n:
                base = base * base
        return out

    def __truediv__(self, other):
        return self * P(other) ** -1

    def __rtruediv__(self, other):
        return P(other) * self ** -1

    def diff(self, name):
        d = {}
        for k, v in self.d.items():
            powers = dict(k)
            e = powers.get(name, 0)
            if e:
                powers[name] -= 1
                d[tuple(sorted((s, t) for s, t in powers.items() if t))] = v * e
        return P(d)

    def names(self):
        return {n for k in self.d for n, _ in k}

    def encode(self):
        return [{'coefficient': str(v), 'powers': dict(k)} for k, v in sorted(self.d.items())]

    def substitute(self, values):
        out = P()
        for k, v in self.d.items():
            term = P(v)
            for name, power in k:
                term = term * P(values.get(name, P.var(name))) ** power
            out += term
        return out


RHO, A, AY, C, TH = [P.var(n) for n in ('rho', 'a', 'ay', 'c', 'theta')]
I, J = range(6), range(4)
SIG = (1, -1, 1, 1, 1, 1)
EZ = (0, 0, 0, 0)


def ordinary(p, axis):
    out = P()
    for name in p.names():
        if name.startswith('z:'):
            sector, degree, suffix = name[2:].split(':')
            suffix += 'y' if axis == 0 else 't'
            derivative = P.var('z:' + sector + ':' + degree + ':' + 'y' * suffix.count('y') + 't' * suffix.count('t'))
        elif name[0] in 'ABCDELM' and set(name[1:]) <= {'y', 't'}:
            suffix = name[1:] + ('y' if axis == 0 else 't')
            derivative = P.var(name[0] + 'y' * suffix.count('y') + 't' * suffix.count('t'))
        elif name.startswith('s') and set(name[1:]) <= {'y', 't'}:
            suffix = name[1:] + ('y' if axis == 0 else 't')
            derivative = P.var('s' + 'y' * suffix.count('y') + 't' * suffix.count('t'))
        elif axis == 0:
            derivative = {'rho': A * RHO, 'a': AY, 'ay': P.var('ayy')}.get(name, P())
        else:
            derivative = {'c': C * TH, 'theta': 1 - TH * TH}.get(name, P())
        out += p.diff(name) * derivative
    return out


@lru_cache(None)
def reduce_monomial(e):
    if e[3] < 2:
        return ((e, Q(1)),)
    base = list(e)
    base[3] -= 2
    d = {}
    for slot in (-1, 0, 1, 2):
        f = base.copy()
        if slot >= 0:
            f[slot] += 2
        for k, v in reduce_monomial(tuple(f)):
            d[k] = d.get(k, 0) + (v if slot == -1 else -v)
    return tuple((k, v) for k, v in d.items() if v)


class S:
    """S3 quotient polynomial with independent Laurent coefficients."""
    def __init__(self, value=0):
        if isinstance(value, S):
            self.d = value.d
        elif isinstance(value, dict):
            self.d = {k: P(v) for k, v in value.items() if P(v)}
        else:
            self.d = {EZ: P(value)} if P(value) else {}

    def __bool__(self):
        return bool(self.d)

    def __add__(self, other):
        other = S(other)
        d = self.d.copy()
        for k, v in other.d.items():
            d[k] = d.get(k, P()) + v
            if not d[k]:
                del d[k]
        return S(d)

    __radd__ = __add__

    def __neg__(self):
        return S({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -S(other)

    def __rsub__(self, other):
        return S(other) + -self

    def __mul__(self, other):
        other = S(other)
        d = {}
        for e, v in self.d.items():
            for f, w in other.d.items():
                for k, q in reduce_monomial(tuple(x + y for x, y in zip(e, f))):
                    d[k] = d.get(k, P()) + q * v * w
        return S(d)

    __rmul__ = __mul__

    def derivative(self, axis):
        return S({e: ordinary(v, axis) for e, v in self.d.items()})

    def partial(self, axis):
        d = {}
        for e, v in self.d.items():
            if e[axis]:
                f = list(e)
                f[axis] -= 1
                d[tuple(f)] = v * e[axis]
        return S(d)

    def encode(self):
        return [{'angular': list(e), 'coefficient': v.encode()} for e, v in sorted(self.d.items())]


X = [S({tuple(int(i == j) for j in J): 1}) for i in J]
K = {(i, j): S(int(i == j)) - X[i] * X[j] for i, j in product(J, repeat=2)}


def intrinsic(t):
    out = {}
    for idx, value in t.items():
        for d in J:
            result = sum((K[d, k] * value.partial(k) for k in J), S())
            for slot, v in enumerate(idx):
                replaced = list(idx)
                replaced[slot] = d
                result += X[v] * t[tuple(replaced)]
            out[(d,) + idx] = result
    return out


def make_bases():
    y = X[0] * X[0] - X[1] * X[1]
    radius = X[0] * X[0] + X[1] * X[1]
    z = radius - Q(1, 2)
    f = y * y - radius * Q(1, 2) + Q(1, 12)
    scalars = [S(1), z, f]
    gradients = [intrinsic({(): v}) for v in (z, f)]
    hessians = []
    for grad in gradients:
        h = intrinsic(grad)
        trace = sum((h[i, i] for i in J), S())
        hessians.append({(i, j): h[i, j] - trace * K[i, j] * Q(1, 3)
                         for i, j in product(J, repeat=2)})
    # Two parity partners from the fixed normalized algebraic Weyl tensor.
    weyl = {k: Q(0) for k in product(J, repeat=4)}
    for (i, j), value in { (0, 1): Q(-2, 3), (2, 3): Q(-2, 3),
                          (0, 2): Q(1, 3), (0, 3): Q(1, 3),
                          (1, 2): Q(1, 3), (1, 3): Q(1, 3)}.items():
        weyl[i, j, i, j] = weyl[j, i, j, i] = value
        weyl[i, j, j, i] = weyl[j, i, i, j] = -value
    eps = {p: (-1) ** sum(p[i] > p[j] for i in J for j in range(i + 1, 4))
           for p in permutations(J)}
    star = {k: sum((Q(eps.get((k[0], k[1], a, b), 0), 2) * weyl[a, b, k[2], k[3]]
                    for a, b in product(J, repeat=2)), Q(0)) for k in weyl}
    tt = [{(i, j): sum(((weyl[i, a, j, b] + sign * star[i, a, j, b]) * Q(1, 2) * X[a] * X[b]
                        for a, b in product(J, repeat=2)), S())
           for i, j in product(J, repeat=2)} for sign in (1, -1)]
    tensors = [{(i, j): value * K[i, j] for i, j in product(J, repeat=2)}
               for value in scalars] + hessians + tt
    return y, scalars, gradients, tensors


Y, SCALARS, GRADS, TENSORS = make_bases()


def connection():
    gamma = {}
    gamma[1, 1, 0] = S(A)
    gamma[0, 1, 1] = S(A)
    for i, j in product(J, repeat=2):
        gamma[j + 2, i + 2, 0] = A * K[i, j]
        gamma[0, i + 2, j + 2] = -A * K[i, j]
        gamma[j + 2, i + 2, 1] = (TH / RHO) * K[i, j]
        gamma[1, i + 2, j + 2] = (TH / RHO) * K[i, j]
    return gamma


GAMMA = connection()


def covariant(t, mutation=None):
    out = {}
    for idx, value in t.items():
        for d in I:
            if d == 0:
                result = value.derivative(0)
            elif d == 1:
                result = (1 if mutation == 'omit_physical_rho' else 1 / RHO) * value.derivative(1)
            else:
                i = d - 2
                result = sum((K[i, k] * value.partial(k) for k in J), S())
                for slot, v in enumerate(idx):
                    if v >= 2:
                        repl = list(idx)
                        repl[slot] = d
                        result += X[v - 2] * t[tuple(repl)]
                result = (1 / C if mutation == 'omit_physical_rho' else 1 / (RHO * C)) * result
            for slot, v in enumerate(idx):
                if mutation == 'omit_derivative_slot' and len(idx) == 2 and slot == 0:
                    continue
                for u in I:
                    gamma = GAMMA.get((u, d, v))
                    if mutation == 'omit_radial_connection' and (u == 0 or v == 0):
                        continue
                    if mutation == 'omit_FRW_connection' and d >= 2 and (u == 1 or v == 1):
                        continue
                    if gamma:
                        repl = list(idx)
                        repl[slot] = u
                        result -= gamma * t[tuple(repl)]
            out[(d,) + idx] = result
    return out


def first_field():
    q = {idx: S() for idx in product(I, repeat=2)}
    q[0, 0] = P.var('A') * Y
    q[0, 1] = q[1, 0] = P.var('B') * Y
    q[1, 1] = P.var('D') * Y
    grad = intrinsic({(): Y})
    hess = intrinsic(grad)
    for i in J:
        q[0, i + 2] = q[i + 2, 0] = P.var('C') * grad[i,]
        q[1, i + 2] = q[i + 2, 1] = P.var('E') * grad[i,]
        for j in J:
            q[i + 2, j + 2] = P.var('L') * Y * K[i, j] + P.var('M') * hess[i, j]
    return q


def decode(row):
    d = {}
    for term in row:
        key = tuple(sorted(term['powers'].items()))
        need(key not in d, 'MAP_DUPLICATE_MONOMIAL')
        d[key] = Q(term['coefficient'])
    return P(d)


def decode_maps(raw):
    return {name: ({field: [decode(row) for row in rows] for field, rows in value.items()}
                   if isinstance(value, dict) else [decode(row) for row in value])
            for name, value in raw['maps'].items()}


def reconstruct_covector(rows):
    out = {(i,): S() for i in I}
    for field, slot in (('y', 0), ('0', 1)):
        out[slot,] = sum((v * b for v, b in zip(rows[field], SCALARS)), S())
    for i in J:
        out[i + 2,] = sum((v * b[i,] for v, b in zip(rows['a'], GRADS)), S())
    return out


def reconstruct_tensor(rows):
    out = {idx: S() for idx in product(I, repeat=2)}
    for field, idx in (('yy', (0, 0)), ('y0', (0, 1)), ('00', (1, 1))):
        out[idx] = sum((v * b for v, b in zip(rows[field], SCALARS)), S())
    out[1, 0] = out[0, 1]
    for i in J:
        out[0, i + 2] = out[i + 2, 0] = sum((v * b[i,] for v, b in zip(rows['ya'], GRADS)), S())
        out[1, i + 2] = out[i + 2, 1] = sum((v * b[i,] for v, b in zip(rows['0a'], GRADS)), S())
        for j in J:
            out[i + 2, j + 2] = sum((v * b[i, j] for v, b in zip(rows['ab'], TENSORS)), S())
    return out


def equal_tensor(name, left, right):
    need(set(left) == set(right), name + ':SHAPE')
    for idx in sorted(left):
        remainder = left[idx] - right[idx]
        if remainder:
            witness = remainder.encode() if isinstance(remainder, S) else remainder.encode()
            raise ValueError(name + ':NONZERO:' + json.dumps({'component': idx, 'witness': witness[:1]}, sort_keys=True, separators=(',', ':')))


def equal_rows(name, left, right):
    need(set(left) == set(right), name + ':SHAPE')
    for field in left:
        need(len(left[field]) == len(right[field]), name + ':LENGTH')
        for i, (v, w) in enumerate(zip(left[field], right[field])):
            if v - w:
                raise ValueError(name + ':NONZERO:' + json.dumps({'field': field, 'mode': i, 'witness': (v - w).encode()[:1]}, sort_keys=True, separators=(',', ':')))


def sphere_integral(poly):
    out = P()
    for e, value in poly.d.items():
        if any(n % 2 for n in e):
            continue
        moment = Q(1)
        for n in e:
            for k in range(n // 2):
                moment *= 2 * k + 1
        for k in range(sum(e) // 2):
            moment /= 4 + 2 * k
        out += value * moment
    return out


def angular_gram_norms():
    norms = {}
    for label, basis in (('scalar', [{(): s} for s in SCALARS]),
                         ('gradient', GRADS), ('tensor', TENSORS)):
        diagonal = []
        for i, j in product(range(len(basis)), repeat=2):
            coefficient = sum((sphere_integral(basis[i][k] * basis[j][k]) for k in basis[i]), P())
            need(set(coefficient.d) <= {()}, 'GRAM_NONCONSTANT_COEFFICIENT')
            value = coefficient.d.get((), Q(0))
            need(bool(value) if i == j else not value, 'GRAM_RANK:' + label + ':' + str((i, j)))
            if i == j:
                diagonal.append(str(value))
        norms[label] = diagonal
    return norms


def coherent_lowering_witness(lowering):
    # A(y)Y dy^2 with all other metric profiles identically zero. Retain
    # free ordinary A jets through differentiation, then specialize at a
    # point. Thus no covariant prolongation is independently invented.
    other = {name: 0 for value in lowering.values() for coefficient in value.d.values()
             for name in coefficient.names() if name[0] in 'BCDELM' and set(name[1:]) <= {'y', 't'}}
    z = {k: S({e: v.substitute(other) for e, v in value.d.items()}) for k, value in lowering.items()}
    dz = covariant(z)
    ddz = covariant(dz)
    # Ricci multiplies the vanishing point value, but is included explicitly.
    ricci = Q(2, 3) * P.var('V')
    wave = {(b,): sum((SIG[d] * ddz[d, d, b] for d in I), S())
            + (ricci + (P.var('p') ** 2 if b == 0 else P())) * z[b,] for b in I}
    values = {'A': 0, 'Ay': 1, 'Ayy': 0, 'Ayyy': 0,
              'At': 0, 'Ayt': 0, 'Att': 0, 'Ayyt': 0, 'Aytt': 0, 'Attt': 0,
              'rho': 1, 'a': Q(1, 2), 'ay': Q(-4, 3), 'ayy': Q(7, 3),
              'theta': 0, 'c': 1, 'p': 1, 'V': 5}
    at_point = lambda tensor: {k: S({e: v.substitute(values) for e, v in value.d.items()}) for k, value in tensor.items()}
    equal_tensor('COHERENT_QC_VALUES_AGREE', at_point(z), {(b,): S() for b in I})
    result = at_point(wave)
    expected = {(b,): 10 * Y * Y if b == 0 else S() for b in I}
    equal_tensor('COHERENT_VALUE_EQUALITY_WAVE_DIFFERENCE', result, expected)
    return {'q': 'A(y) Y dy^2; other six metric profiles identically zero',
            'point_ordinary_jets': {k: str(v) for k, v in values.items()},
            'background_consistency': 'rho=p=1,a=1/2,ay=-4/3,V=5,Vp=0 gives ayy=7/3; c=1,theta=0',
            'QC_actual_minus_QC_background_value': 0,
            'L0_difference_radial_component': result[0,].encode(),
            'L0_difference_formula_at_point': '10 Y^2 = 20 a Y^2',
            'conclusion': 'Equality of QC values at one point does not permit exchange inside L0.'}


def background_reduce(poly):
    # Einstein constraint is algebraic in addition to the derivative rules.
    pp = P.var('p')
    poly = poly.substitute({'ay': -RHO ** -2 - pp * pp / 3,
                            'V': 6 * RHO ** -2 + pp * pp / 2 - 6 * A * A})
    out = P()
    for key, value in poly.d.items():
        powers = dict(key)
        theta_power = powers.pop('theta', 0)
        need(theta_power >= 0, 'NEGATIVE_THETA_POWER')
        for k in range(theta_power // 2 + 1):
            mon = powers.copy()
            mon['theta'] = theta_power % 2
            mon['c'] = mon.get('c', 0) - 2 * k
            out += P({tuple(sorted((n, e) for n, e in mon.items() if e)):
                      value * comb(theta_power // 2, k) * (-1) ** k})
    return out


def arbitrary_covector_operator(mutation=None):
    """Direct full-frame rough wave, independent of every producer formula."""
    z = {(b,): S() for b in I}
    for sector, b in (('radial', 0), ('temporal', 1)):
        z[b,] = sum((P.var('z:' + sector + ':' + str(l) + ':') * base
                     for l, base in zip((0, 2, 4), SCALARS)), S())
    for i in J:
        z[i + 2,] = sum((P.var('z:angular:' + str(l) + ':') * base[i,]
                        for l, base in zip((2, 4), GRADS)), S())
    dz = covariant(z)
    symgrad = {(a, b): (dz[a, b] + dz[b, a]) * Q(1, 2) for a, b in product(I, repeat=2)}
    dsym = covariant(symgrad)
    b0sym = {(b,): sum((SIG[d] * (dsym[d, d, b] - dsym[b, d, d] * Q(1, 2)) for d in I), S()) for b in I}
    altered = mutation if mutation in ('omit_derivative_slot', 'omit_radial_connection', 'omit_FRW_connection', 'omit_physical_rho') else None
    wave_dz = covariant(z, altered)
    ddz = covariant(wave_dz, altered)
    wave_sign = (1,) * 6 if mutation == 'positive_time_wave' else SIG
    ricci = Q(2, 3) * P.var('V')
    wave = {}
    for b in I:
        endomorphism = -ricci if b == 1 and mutation == 'ricci_time_minus' else ricci
        if b == 0 and mutation != 'omit_radial_p2':
            endomorphism += P.var('p') ** 2
        wave[b,] = sum((wave_sign[d] * ddz[d, d, b] for d in I), S()) + endomorphism * z[b,]
    factor = 1 if mutation == 'missing_half_operator' else Q(1, 2)
    for b in I:
        remainder = b0sym[b,] - factor * wave[b,]
        reduced = S({e: background_reduce(v) for e, v in remainder.d.items()})
        if reduced:
            raise ValueError('ARBITRARY_COVECTOR_OPERATOR:NONZERO:' + json.dumps({'component': b, 'witness': reduced.encode()[:1]}, sort_keys=True, separators=(',', ':')))
    names = {name for value in ddz.values() for coefficient in value.d.values()
             for name in coefficient.names() if name.startswith('z:')}
    need(len(names) == 48, 'COVECTOR_JET_COUNT:' + str(len(names)))
    rows = {}
    fields = {'y': {(): wave[0,]}, '0': {(): wave[1,]},
              'a': {(i,): wave[i + 2,] for i in J}}
    for label, basis in (('y', [{(): s} for s in SCALARS]),
                         ('0', [{(): s} for s in SCALARS]), ('a', GRADS)):
        coefficients = []
        for base in basis:
            numerator = sum((sphere_integral(fields[label][idx] * v) for idx, v in base.items()), P())
            denominator = sum((sphere_integral(v * v) for v in base.values()), P())
            coefficients.append(background_reduce(numerator / denominator))
        reconstructed = {idx: sum((co * base[idx] for co, base in zip(coefficients, basis)), S()) for idx in fields[label]}
        for idx in fields[label]:
            residual = fields[label][idx] - reconstructed[idx]
            need(not S({e: background_reduce(v) for e, v in residual.d.items()}), 'WAVE_MODAL_RECONSTRUCTION:' + label + ':' + str(idx))
        rows[label] = [value.encode() for value in coefficients]
    return {'identity': 'B0(sym nabla Z)=(Box0+Ric0)Z/2',
            'domain': 'Arbitrary eight retained covector modes, each with six free ordinary jets through order two',
            'free_ordinary_jets': len(names), 'full_ambient_covector_components': 6,
            'angular_coefficient_remainders': 0,
            'Ricci_endomorphism': '(2/3)V Z_b+p^2 delta_b^y Z_y',
            'independent_eight_modal_wave_rows': rows,
            'method': 'Own explicit connection, full first and second covariant derivatives including derivative slot, direct Lorentz contraction; no modal wave expressions imported.'}


MAP_REL = 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
MAP_SHA = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
NAT_REL = 'research/HDBLAST_NATURAL_FIELD_BINDING_20261030'
NAT_MANIFEST_SHA = '5121e4304fad7ab5856977b38c6ebd9fdc205900b543bc1c3d97384d246f378f'
ALLOWED_IMPORTS = {'argparse', 'ast', 'fractions', 'functools', 'hashlib', 'itertools', 'json', 'pathlib', 'sys', 'math'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(root):
    map_path = root / MAP_REL
    need(map_path.is_file() and not map_path.is_symlink() and sha(map_path) == MAP_SHA, 'MAP_PIN')
    manifest_path = root / NAT_REL / 'MANIFEST.json'
    need(sha(manifest_path) == NAT_MANIFEST_SHA, 'NATURAL_MANIFEST_PIN')
    closure = {MAP_REL: MAP_SHA, NAT_REL + '/MANIFEST.json': NAT_MANIFEST_SHA}
    manifest = json.loads(manifest_path.read_text())
    need(len(manifest['files']) == manifest['file_count_excluding_this_manifest'], 'NATURAL_MANIFEST_COUNT')
    for row in manifest['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts, 'NATURAL_MANIFEST_PATH')
        path = root / NAT_REL / row['path']
        need(path.is_file() and not path.is_symlink() and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'NATURAL_CLOSURE:' + row['path'])
        closure[NAT_REL + '/' + row['path']] = row['sha256']
    raw = json.loads(map_path.read_text())
    geo = root / 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    for name, pin in raw['engine_pins'].items():
        path = geo / 'geometry' / name
        need(path.is_file() and not path.is_symlink() and sha(path) == pin, 'ENGINE_CLOSURE:' + name)
        closure[str(path.relative_to(root))] = pin
    for name, pin in raw['parents'].items():
        path = geo / 'projection/inputs' / name
        need(path.is_file() and not path.is_symlink() and sha(path) == pin, 'MAP_PARENT:' + name)
        closure[str(path.relative_to(root))] = pin
    # Exactly one non-stdlib Python execution unit. Authenticate imports and
    # absence of optimization-sensitive checks before and after computation.
    this = Path(__file__).resolve()
    tree = ast.parse(this.read_text())
    imports = set()
    for node in ast.walk(tree):
        need(not isinstance(node, ast.Assert), 'ASSERT_IN_EXECUTED_CLOSURE')
        if isinstance(node, ast.Import):
            imports.update(a.name.split('.')[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.add(node.module.split('.')[0])
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            need(node.func.id not in {'eval', 'exec', '__import__'}, 'DYNAMIC_CODE_IN_EXECUTED_CLOSURE')
    need(imports == ALLOWED_IMPORTS, 'IMPORT_ALLOWLIST:' + str(sorted(imports)))
    closure['checker_sha256'] = sha(this)
    return {'files': closure, 'executed_nonstdlib_closure': [this.name], 'imports': sorted(imports), 'assertions': 0}


CONTROLS = ('lowering_factor_one', 'lowering_euclidean', 'pc1_factor_one',
            'pc1_euclidean', 'premature_background', 'unnormalized_sym',
            'reduced_gradient_sign', 'drop_PC1', 'delete_TT', 'natural_for_general',
            'principal_for_general', 'harmonic_scalar', 'omit_radial_p2',
            'ricci_time_minus', 'missing_half_operator', 'omit_derivative_slot',
            'omit_radial_connection', 'omit_FRW_connection', 'positive_time_wave', 'omit_physical_rho')


def run(root, mutation=None):
    before = authenticate(root)
    raw = json.loads((root / MAP_REL).read_text())
    mapped = decode_maps(raw)
    q = first_field()
    dq = covariant(q)
    lowerP = {(d, a, b): (dq[a, b, d] + dq[b, a, d] - dq[d, a, b]) * Q(1, 2)
              for d, a, b in product(I, repeat=3)}
    c1 = {(d,): sum((SIG[a] * lowerP[d, a, a] for a in I), S()) for d in I}
    # C1 is independently calculated by the trace-gradient formula too.
    c1_second = {(d,): sum((SIG[a] * (dq[a, a, d] - dq[d, a, a] * Q(1, 2)) for a in I), S()) for d in I}
    equal_tensor('C1_LINEAR_CONNECTION_CONTRACTION', c1, c1_second)
    actual = reconstruct_covector(mapped['QC_actual'])
    background = reconstruct_covector(mapped['QC_background'])
    lowering_sign = (1,) * 6 if mutation == 'lowering_euclidean' else SIG
    lowering_factor = 1 if mutation == 'lowering_factor_one' else 2
    lowering = {(b,): sum((lowering_factor * lowering_sign[d] * q[b, d] * c1[d,] for d in I), S()) for b in I}
    equal_tensor('LOWERING_FULL_ANGULAR', {k: actual[k] - background[k] for k in actual}, lowering)
    # Independently form the full-lowered Hessian rather than only its difference.
    direct_actual = {(d,): -2 * sum((SIG[a] * SIG[b] * q[a, b] * lowerP[d, a, b] for a, b in product(I, repeat=2)), S()) for d in I}
    equal_tensor('QC_ACTUAL_DIRECT_FULL_ANGULAR', actual, direct_actual)
    equal_rows('QC_ALIAS', mapped['QC'], mapped['QC_actual'])
    pc_sign = (1,) * 6 if mutation == 'pc1_euclidean' else SIG
    pc_factor = 1 if mutation == 'pc1_factor_one' else 2
    pc1 = {(a, b): pc_factor * sum((pc_sign[d] * lowerP[d, a, b] * c1[d,] for d in I), S()) for a, b in product(I, repeat=2)}
    equal_tensor('PC1_FULL_ANGULAR', reconstruct_tensor(mapped['PC1']), pc1)
    selected_qc = background if mutation == 'premature_background' else actual
    dqc = covariant(selected_qc)
    sym_factor = 1 if mutation == 'unnormalized_sym' else Q(1, 2)
    gradient = {(a, b): (dqc[a, b] + dqc[b, a]) * sym_factor for a, b in product(I, repeat=2)}
    equal_tensor('GAUGE_GRADIENT_FULL_ANGULAR', reconstruct_tensor(mapped['gauge_gradient']), gradient)
    if mutation == 'delete_TT':
        mapped['QEred_general']['ab'][5:] = [P(), P()]
    if mutation == 'natural_for_general':
        mapped['QEred_general'] = mapped['QE']
    if mutation == 'principal_for_general':
        mapped['QEred_general'] = mapped['metric_forcing']
    sign = 1 if mutation == 'reduced_gradient_sign' else -1
    expected = {f: [v + sign * w for v, w in zip(mapped['QE'][f], mapped['gauge_gradient'][f])] for f in mapped['QE']}
    equal_rows('REDUCED_HARMONIC_ALL20', mapped['QEred'], expected)
    expected_general = {f: [v + (P() if mutation == 'drop_PC1' else w) for v, w in zip(expected[f], mapped['PC1'][f])] for f in expected}
    equal_rows('REDUCED_GENERAL_ALL20', mapped['QEred_general'], expected_general)
    # Scalar general correction directly from the independently built C1.
    eta = {(): P.var('s') * Y}
    deta = covariant(eta)
    scalar_correction = 2 * sum(((P.var('p') * SIG[d] * q[0, d] - SIG[d] * deta[d,]) * c1[d,] for d in I), S())
    qphi_general = mapped['Qphi'] if mutation == 'harmonic_scalar' else mapped['Qphi_general']
    equal_tensor('SCALAR_GENERAL_CORRECTION_FULL_ANGULAR', {(): sum(((a - b) * base for a, b, base in zip(qphi_general, mapped['Qphi'], SCALARS)), S())}, {(): scalar_correction})
    norms = angular_gram_norms()
    witness = coherent_lowering_witness(lowering)
    operator = arbitrary_covector_operator(mutation)
    need(authenticate(root) == before, 'SOURCE_CLOSURE_CHANGED')
    return {'status': 'PASS_INDEPENDENT_EXACT_REGISTERED_SOURCE_ROLES',
            'authentication': before,
            'map_sha256': MAP_SHA,
            'independence': 'Own Laurent ring, quotient sphere calculus, basis reconstruction and explicit upper connection table; no scientific source is imported.',
            'first_field_domain': 'Registered seven-profile l=2 metric ansatz; all 42 ordinary profile jets through order two are free. No harmonic or first-natural equation imposed.',
            'angular_comparisons': {'lowering_covector_modes': 8, 'QC_actual_covector_modes': 8, 'PC1_tensor_modes': 20, 'gauge_gradient_tensor_modes': 20, 'scalar_general_modes': 3},
            'modal_joins': {'QEred': 20, 'QEred_general': 20, 'QC_alias': 8},
            'TT_source_terms_retained': {name: [len(mapped[name]['ab'][i].d) for i in (5, 6)] for name in ('QE', 'QEred', 'QEred_general', 'PC1', 'gauge_gradient')},
            'background_field_equations_needed_for_these_roles': False,
            'kinematic_derivative_relations_for_roles': 'rho_y=a rho; a_y=ay; c_t=c theta; theta_t=1-theta^2',
            'full_angular_remainders': 0,
            'independent_angular_Gram_diagonal_nonzero': norms,
            'coherent_pointwise_lowering_counterexample': witness,
            'independent_arbitrary_covector_operator': operator,
            'scope': 'Exact source-role/operator joins only. Generic reduced Bianchi and exported modal wave formulas are separately produced. No new PDE existence, external peer review, accuracy improvement or cosmological origin claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mutation', choices=CONTROLS)
    args = parser.parse_args()
    need(not args.output.exists(), 'OUTPUT_MUST_BE_ABSENT')
    report = run(args.root.resolve(), args.mutation)
    with args.output.open('x') as output:
        json.dump(report, output, sort_keys=True, indent=2)
        output.write('\n')
    print(report['status'])


if __name__ == '__main__':
    main()
