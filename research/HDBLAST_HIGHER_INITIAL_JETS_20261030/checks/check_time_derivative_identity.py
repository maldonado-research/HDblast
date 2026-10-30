#!/usr/bin/env python3
"""Exact independent AD of the authenticated original full23 assembly.

The manufactured rational backgrounds validate an algebraic derivative
identity. They are not enclosures of the actual background or source F_x.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
OP_PIN = 'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


class Dual:
    """A first-order rational jet, independent of the interval producer."""
    def __init__(self, value=0, derivative=0):
        self.v, self.d = Q(value), Q(derivative)

    @staticmethod
    def cv(value):
        return value if isinstance(value, Dual) else Dual(value)

    def __add__(self, other):
        other = self.cv(other)
        return Dual(self.v+other.v, self.d+other.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, other):
        return self + (-self.cv(other))

    def __rsub__(self, other):
        return self.cv(other) + (-self)

    def __mul__(self, other):
        other = self.cv(other)
        return Dual(self.v*other.v, self.d*other.v+self.v*other.d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.cv(other)
        need(other.v != 0, 'nonzero dual denominator')
        return Dual(self.v/other.v, (self.d*other.v-self.v*other.d)/other.v**2)

    def __rtruediv__(self, other):
        return self.cv(other)/self

    def __pow__(self, n):
        need(isinstance(n, int), 'integer dual power')
        if n == 0:
            return Dual(1)
        return Dual(self.v**n, n*self.v**(n-1)*self.d)


def derivative(matrix):
    return [[Dual.cv(z).d for z in row] for row in matrix]


def difference(a, b):
    return [[x-y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def nonzero(matrix):
    return any(z for row in matrix for z in row)


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
args = ap.parse_args()
path = args.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py'
before = path.read_bytes()
need(hashlib.sha256(before).hexdigest() == OP_PIN, 'INPUT_HASH:original linear operator')
spec = importlib.util.spec_from_file_location('time_AD_original_operator', path)
op = importlib.util.module_from_spec(spec)
spec.loader.exec_module(op)
m = op.matrices()
# The original assemble only reads the primitive matrices. Reusing the exact
# authenticated primitive object removes repeated construction from fixtures.
op.matrices = lambda: m
groups, controls = [], []
matrix_entries = cases = 0
last = None


def expected_terms(ell, rho, a, p, tau, s):
    return [(-2*ell**3*s*s*tau, m['Delta']),
            (-2*ell**3*rho*a*s*tau, m['Ca']),
            (2*ell**3*s*(s*s-tau*tau), m['Ch']),
            (ell**3*rho*a*s*s, op.add(m['Sah'], op.scale(-3, m['O']))),
            (2*ell**3*tau*s*s, m['Sh']),
            (3*ell**3*rho*p*s*s, m['W_L']),
            (-ell**3*rho*p*s*tau, m['Ephi'])]


def manufactured(seed, time_scale=True, omit_h_derivative=False, reverse_sech=False):
    u = Q(seed, seed+7)
    tau, s = 2*u/(1+u*u), (1-u*u)/(1+u*u)
    need(tau*tau+s*s == 1 and s > 0, 'exact rational hyperbolic value relation')
    ell, rho, a, p = Q(seed+1, 1000), Q(seed+5, 3), Q(seed-2, 7), -Q(seed+4, 11)
    factor = ell if time_scale else Q(1)
    tjet = Dual(tau, factor*s*s)
    sjet = Dual(s, (1 if reverse_sech else -1)*factor*s*tau)
    hjet = tjet/rho
    if omit_h_derivative:
        hjet = Dual(hjet.v)
    # Independent original assembly: curvature/potential/normal derivatives
    # are radial-background data and therefore have zero x jet.
    result = op.assemble(ell, rho, a, hjet, 1/sjet, tjet, p,
        Q(seed+2, 13), Q(seed-3, 17), Q(seed+3, 19),
        Q(seed+7, 23), Q(seed-5, 29), Q(seed+11, 31),
        Q(seed-7, 37), Q(seed+13, 41))
    return result, (ell, rho, a, p, tau, s)


for seed in (1, 3, 5, 11):
    result, scalars = manufactured(seed)
    ell, rho, a, p, tau, s = scalars
    D_expected = op.add(*(op.scale(c, matrix) for c, matrix in expected_terms(*scalars)))
    A_expected = op.scale(-3*ell*ell*s*s, m['I'])
    targets = {'A': A_expected, 'B': op.zeros(), 'D': D_expected, 'M': op.zeros()}
    for key, expected in targets.items():
        got = derivative(result[key])
        need(got == expected, 'original rational AD identity:'+key+':'+str(seed))
        groups.append('original_AD_'+key+'_seed_'+str(seed))
        matrix_entries += 23*23
    cases += 1
    last = scalars

# Genuine wrong-formula controls change the original AD result or drop a
# nonzero analytic product term; equality is checked on complete matrices.
correct, scalars = manufactured(3)
correct_D = derivative(correct['D'])
for name, kwargs in [('physical_t_instead_of_scaled_x', {'time_scale': False}),
                     ('omitted_h_x_product_rule', {'omit_h_derivative': True}),
                     ('wrong_sech_derivative_sign', {'reverse_sech': True})]:
    wrong, _ = manufactured(3, **kwargs)
    need(derivative(wrong['D']) != correct_D, 'reject wrong derivative convention:'+name)
    controls.append(name)
terms = expected_terms(*scalars)
for j, name in [(0, 'dropped_s_squared_Delta_derivative'),
                (2, 'dropped_cosmological_connection_product'),
                (3, 'dropped_Sah_minus_3O_term'),
                (5, 'dropped_scalar_metric_time_coupling'),
                (6, 'dropped_scalar_angular_sech_derivative')]:
    wrong = op.add(*(op.scale(c, matrix) for i, (c, matrix) in enumerate(terms) if i != j))
    need(wrong != correct_D, 'reject omitted actual coefficient term:'+name)
    controls.append(name)
need(nonzero(derivative(correct['A'])) and not nonzero(derivative(correct['B'])) and
     not nonzero(derivative(correct['M'])), 'exact time-independent B and M')
groups.append('B_and_M_x_exactly_zero')
need(path.read_bytes() == before, 'original operator readback unchanged')
print(json.dumps({'status': 'PASS_EXACT_ORIGINAL_COEFFICIENT_TIME_AD',
    'input_sha256': {'original_linear_operator': OP_PIN},
    'named_groups': groups, 'group_count': len(groups),
    'manufactured_background_cases': cases, 'matrix_entry_equalities': matrix_entries,
    'wrong_formula_controls': controls, 'wrong_formula_control_count': len(controls),
    'derivative_convention': 'partial_x at fixed r, t=ell*x; radial background data have zero x jet',
    'scope': 'exact independent automatic differentiation of the authenticated original full23 coefficient assembly on rational manufactured backgrounds; no actual bulk-source derivative, C2 error or Einstein residual bound'}, sort_keys=True))
