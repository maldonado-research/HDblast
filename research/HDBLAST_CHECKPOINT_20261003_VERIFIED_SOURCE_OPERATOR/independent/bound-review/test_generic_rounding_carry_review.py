#!/usr/bin/env python3
"""Independent entire-kernel reference for fabricated rounded/carry fixtures."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path


BASE = Path(__file__).resolve().parent.parent
SPEC = importlib.util.spec_from_file_location('reviewed_generic_operator_definitions', BASE / 'generic_operator.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def need(test, name):
    if not test:
        raise RuntimeError(name)


def mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def scale(a, q):
    return (a[0] * q, a[1] * q)


def abs_upper(a):
    return abs(a[0]) + abs(a[1])


def kernels(omega, h, count=200):
    power = (Q(1), Q(0))
    out = [(Q(0), Q(0)) for _ in range(3)]
    for n in range(count + 1):
        for offset in range(3):
            out[offset] = add(out[offset], scale(power, h ** (n + offset) / factorial(n + offset)))
        power = mul(power, (Q(0), omega))
    x = abs(omega * h)
    radii = tuple(h ** j * x ** (count + 1) / factorial(count + 1 + j)
                  for j in range(3))
    return out, radii


def contained(truth, truth_radius, candidate, radius):
    delta = add(truth, scale(candidate, Q(-1)))
    margin = radius - truth_radius
    return margin >= 0 and delta[0] ** 2 + delta[1] ** 2 <= margin ** 2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = 0
    h = Q(1, 64)
    gg = Q(3, 7)
    w0 = (Q(2, 11), Q(-1, 13))
    u0 = (Q(-1, 17), Q(3, 19))
    source_error = Q(1, 1000000)
    dw = Q(1, 700000)
    du = Q(1, 900000)
    for omega in (Q(0), Q(1, 2 ** 39), Q(32), Q(512)):
        certificate = MODULE.defect_certificate((gg,), source_error, h, omega, 96,
                     coefficient_bits=16, incoming_w=w0, incoming_u=u0,
                     incoming_w_error=dw, incoming_u_error=du)
        need(any(a != (Q(0), Q(0)) for a in certificate['residual_coefficients']),
             'deliberate coarse coefficient rounding is present in actual ODE defect')
        checks += 1
        (e0, e1, e2), (r0, r1, r2) = kernels(omega, h)
        for sign in (-1, 1):
            actual_g = gg + sign * source_error
            actual_w0 = add(w0, (sign * dw, Q(0)))
            actual_u0 = add(u0, (Q(0), -sign * du))
            me = add(scale(e1, actual_g), scale(mul(e0, actual_w0), Q(-1)))
            mer = abs(actual_g) * r1 + abs_upper(actual_w0) * r0
            dd = add(add(scale(e2, actual_g), scale(actual_u0, Q(-1))),
                     scale(mul(e1, actual_w0), Q(-1)))
            ddr = abs(actual_g) * r2 + abs_upper(actual_w0) * r1
            for field, truth, tr in (('Mexp', me, mer), ('Duhamel_u', dd, ddr)):
                center, radius = certificate[field]
                need(contained(truth, tr, center, radius),
                     'rounded recurrence plus inherited uncertainty encloses entire-kernel reference')
                checks += 1
        need(certificate['Mexp'][1] >= dw + h * source_error,
             'incoming-w and forcing radius retained')
        need(certificate['Duhamel_u'][1] >= du + h * dw + h * h * source_error / 2,
             'incoming-u/w and nested forcing radius retained')
        checks += 2
    for radius in (Q(0), Q(1, 7), Q(17, 13 ** 600), Q(1, 15 * 2 ** 90)):
        rounded = MODULE.outward_radius(radius, 512)
        need(radius <= rounded < radius + Q(1, 2 ** 512),
             'actual outward helper is exact dyadic ceiling')
        need(rounded.denominator & (rounded.denominator - 1) == 0 and
             rounded.denominator <= 2 ** 512, 'radius denominator growth is bounded')
        checks += 2
    for bad in (-Q(1, 7), 0.0):
        try:
            MODULE.outward_radius(bad)
        except ValueError:
            checks += 1
        else:
            raise RuntimeError('negative or nonexact radius accepted')
    receipt = {'status': 'PASS_FABRICATED_ROUNDING_AND_CARRY_REVIEW',
               'checks': checks, 'coarse_chosen_coefficient_bits': 16, 'kernel_degree': 96,
               'registered_source_model_calls': 0, 'physical_source_evaluations': 0,
               'scientific_input_arrays': 0,
               'generic_operator_sha256': hashlib.sha256((BASE / 'generic_operator.py').read_bytes()).hexdigest()}
    encoded = json.dumps(receipt, sort_keys=True, indent=2) + '\n'
    if args.output:
        args.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
