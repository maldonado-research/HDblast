#!/usr/bin/env python3
"""Exact algebra/continuous-bound constants only; no background or bump query.

The analytic proof supplies the meaning of the formal differential identities,
flatness, norm changes of variables and Sobolev estimates. This small verifier
checks the nontrivial coefficient algebra and rational constants they use.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def need(condition, message):
    if not condition:
        raise ValueError(message)


# Sparse rational polynomials with named formal variables.
VARS = ('A', 'Ar', 'b', 'br', 'brr', 'g', 'gr', 'u', 'w')
ZERO = (0,) * len(VARS)


def const(value):
    return {ZERO: F(value)} if value else {}


def var(name):
    power = list(ZERO)
    power[VARS.index(name)] = 1
    return {tuple(power): F(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for power, coefficient in poly.items():
            out[power] = out.get(power, F()) + coefficient
    return {power: coefficient for power, coefficient in out.items() if coefficient}


def mul(left, right):
    out = {}
    for p, a in left.items():
        for q, b in right.items():
            power = tuple(x + y for x, y in zip(p, q))
            out[power] = out.get(power, F()) + a * b
    return {power: coefficient for power, coefficient in out.items() if coefficient}


def scale(poly, scalar):
    return mul(poly, const(scalar))


def partial(poly, name):
    index = VARS.index(name)
    out = {}
    for power, coefficient in poly.items():
        if power[index]:
            fresh = list(power)
            fresh[index] -= 1
            out[tuple(fresh)] = coefficient * power[index]
    return out


def radial(poly):
    return add(*(mul(partial(poly, old), var(new))
                 for old, new in [('A', 'Ar'), ('b', 'br'), ('br', 'brr'), ('g', 'gr')]))


def encode(poly):
    return [{'coefficient': str(coefficient),
             'powers': {name: exponent for name, exponent in zip(VARS, power) if exponent}}
            for power, coefficient in sorted(poly.items())]


class Rejected(ValueError):
    pass


class QuietParser(argparse.ArgumentParser):
    def error(self, message):
        raise Rejected('invalid arguments')


def validate_output(path):
    if not path.is_absolute():
        raise Rejected('absolute output path')
    if any(part.is_symlink() for part in [path, *path.parents]):
        raise Rejected('output path contains symlink')
    if path != path.resolve():
        raise Rejected('normalized output path')
    if not path.parent.is_dir():
        raise Rejected('output parent exists')
    if path.exists():
        raise Rejected('output absent')
    if path.is_relative_to(Path(__file__).resolve().parent):
        raise Rejected('output outside component')


def main():
    parser = QuietParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    validate_output(args.output)
    ell = F(1, 10**6)
    A, Ar, b, br, brr, g, gr, u, w = map(var, VARS)
    one = const(1)
    # Direct raw substitution, in the z,z_r,z_rr coefficient basis.
    h = [scale(b, 3), const(-1), {}]
    hr = [scale(br, 3), scale(b, 3), const(-1)]
    raw = [mul(A, add(scale(one, 3 * ell**2) if j == 0 else {},
                       mul(add(scale(b, -1), g), h[j]), scale(hr[j], -1)))
           for j in range(3)]
    c = add(scale(b, -2), scale(g, -1))
    d = scale(add(const(ell**2), scale(mul(b, b), -1), mul(b, g), scale(br, -1)), 3)
    expected = [mul(A, d), mul(A, c), A]
    need(raw == expected, 'raw-to-profile differential operator identity')
    direct_derivative = [add(radial(raw[j]) if j < 3 else {}, raw[j - 1] if j else {})
                         for j in range(4)]
    expected_derivative = [add(mul(Ar, d), mul(A, radial(d))),
                           add(mul(Ar, c), mul(A, radial(c)), mul(A, d)),
                           add(Ar, mul(A, c)), A]
    need(direct_derivative == expected_derivative, 'profile derivative coefficient identity')

    # For E=exp(1-w), w=(1-u^2)^-1, d(E P)/du=E*(P_u+2u w^2 P_w-2u w^2 P).
    uw2 = mul(u, mul(w, w))
    bump = [mul(u, mul(u, u))]
    for unused in range(3):
        P = bump[-1]
        bump.append(add(partial(P, 'u'), scale(mul(uw2, partial(P, 'w')), 2),
                        scale(mul(uw2, P), -2)))
    jets = []
    caps = []
    for P in bump:
        jets.append(sum((coefficient for powers, coefficient in P.items()
                         if powers[VARS.index('u')] == 0), F()))
        cap = F()
        for powers, coefficient in P.items():
            need(all(exponent == 0 for name, exponent in zip(VARS, powers) if name not in ('u', 'w')),
                 'bump derivative polynomial contains only u,w')
            exponent = powers[VARS.index('w')]
            cap += abs(coefficient) * (exponent**exponent if exponent else 1)
        caps.append(cap)
    need(jets == [0, 0, 0, 6], 'exact centered bump jets')
    need(caps == [1, 11, 1302, 550638], 'whole-support rational bump derivative caps')
    # Integral scaling: n^(j-5/2) pointwise, n^(-1/2) for the L2 support width.
    l2_n = [F(j) - F(5, 2) - F(1, 2) for j in range(4)]
    l2_m = [2 * exponent for exponent in l2_n]
    need(l2_m == [-6, -4, -2, 0], 'localized L2 scaling exponents')
    r0, first_m, T = F(9, 8), 4, F(5, 4)
    need(1 < r0 - F(1, first_m**2) < r0 + F(1, first_m**2) < T,
         'every m>=4 support is strictly interior')
    need(F(6) * F(3, 2) == 9 and -F(2, 3) * F(1, 2) * 9 == -3,
         'exact D_r and conditioned scalar derivative center factors')
    length = T - 1
    need(length == F(1, 4) and length**3 / 3 < length and
         1 / length + length == F(17, 4), 'H2 bridge rational constants')
    need(F(1, 2) * F(4 * 10**33) == F(2 * 10**33), 'anchored H2 sufficient-class budget')
    result = {
        'status': 'PASS_EXACT_FORMAL_ALGEBRA_AND_CONTINUOUS_BOUND_CONSTANTS',
        'scope': 'Formal identities and rational constants only; no actual-background evaluation, oracle, PDE solve or numerical bump evaluation.',
        'profile_operator_coefficients_z_zr_zrr': list(map(encode, raw)),
        'profile_derivative_coefficients_z_zr_zrr_zrrr': list(map(encode, direct_derivative)),
        'bump': 'u^3 exp(1-1/(1-u^2)) for |u|<1, zero otherwise',
        'bump_derivative_polynomials_after_factoring_E': list(map(encode, bump)),
        'exact_bump_jets_at_zero': list(map(str, jets)),
        'continuous_bump_sup_caps_j0_to_3': list(map(str, caps)),
        'rational_safe_bump_L2_caps': [str(2 * cap) for cap in caps],
        'm_subsequence_L2_exponents_j0_to_3': list(map(str, l2_m)),
        'support_center': str(r0), 'minimum_integer_m': first_m,
        'exact_center_formulas': {'D_r': '9*m/g(r0)', 's_r': '-3*p(r0)*m/(a(r0)*g(r0))'},
        'H2_bridge': {'interval_length': str(length), 'generic_constant_squared': '17/4',
                      'anchored_constant': '1/2', 'anchored_sufficient_budget': '4000000000000000000000000000000000'},
        'analytic_arguments_in_proof': ['flatness of E at support endpoints',
            'max exp(1-w) w^k <= k^k for w>=1', 'L2 change of variables',
            'linear profile-image admissibility and fixed canonical sectors',
            'H1 budget scaling and unbounded derivative supremum', 'one-dimensional anchored H2 estimate'],
        'new_activity': {'actual_background_queries': 0, 'oracle_calls': 0, 'PDE_solves': 0,
                         'numerical_bump_queries': 0},
    }
    payload = (json.dumps(result, indent=2, sort_keys=True) + '\n').encode()
    with args.output.open('xb') as stream:
        stream.write(payload)
    print(json.dumps({'status': result['status'], 'sha256': hashlib.sha256(payload).hexdigest()}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (Rejected, ValueError) as exc:
        print(json.dumps({'status': 'REJECTED', 'witness': str(exc)}, sort_keys=True))
        raise SystemExit(2)
    except OSError:
        print(json.dumps({'status': 'REJECTED', 'witness': 'filesystem error'}, sort_keys=True))
        raise SystemExit(2)
