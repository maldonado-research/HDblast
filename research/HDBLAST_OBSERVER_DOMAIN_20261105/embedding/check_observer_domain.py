#!/usr/bin/env python3
"""Exact finite companions to OBSERVER_DOMAIN.md; no analytic oracle."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
EXPECTED = {
    'README.md', 'OBSERVER_DOMAIN.md', 'SOURCE_PINS.json',
    'check_observer_domain.py', 'inputs/README.source.md',
    'inputs/OBSERVER_FLOW_BUDGET.source.md',
    'inputs/AREA_DENSITY_BOUND.source.md', 'inputs/SOURCE_DEPENDENCIES.source.md',
}


def canonical(obj):
    return (json.dumps(obj, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True) + '\n').encode('ascii')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode('ascii') + b'\0' + data).hexdigest()


def authenticate():
    mp = ROOT / 'MANIFEST.json'
    require(not mp.is_symlink(), 'manifest symlink')
    raw = mp.read_bytes()
    manifest = json.loads(raw)
    require(raw == canonical(manifest), 'noncanonical manifest')
    entries = manifest['files']
    require(set(entries) == EXPECTED, 'manifest inventory mismatch')
    observed = set()
    for p in ROOT.rglob('*'):
        require(not p.is_symlink(), 'payload symlink')
        if p.is_file():
            observed.add(p.relative_to(ROOT).as_posix())
    require(observed == EXPECTED | {'MANIFEST.json'}, 'closed inventory mismatch')
    for name in sorted(entries):
        data = (ROOT / name).read_bytes()
        require(entries[name] == {'bytes': len(data), 'sha256': sha(data)},
                'manifest bytes mismatch: ' + name)
    pins_raw = (ROOT / 'SOURCE_PINS.json').read_bytes()
    pins = json.loads(pins_raw)
    require(pins_raw == canonical(pins), 'noncanonical source pins')
    require(pins['verified_local_head'] == '21b652a69f5b2c49b24cec17971a6d6b44d4afdd',
            'unexpected source head')
    require(pins['verified_local_tree'] == '1f710b94410a19151df6e02a462ebc341ffeb635',
            'unexpected source tree')
    require(len(pins['sources']) == 4, 'source count')
    require({s['copy'] for s in pins['sources']} ==
            {p for p in EXPECTED if p.startswith('inputs/')}, 'source inventory')
    for source in pins['sources']:
        data = (ROOT / source['copy']).read_bytes()
        require(len(data) == source['bytes'] and sha(data) == source['sha256']
                and blob(data) == source['git_blob_sha1'],
                'source bytes mismatch: ' + source['copy'])
    return {'manifest_sha256': sha(raw), 'source_pins_sha256': sha(pins_raw),
            'authenticated_payload_files': len(entries), 'authenticated_sources': 4}


# Sparse exact polynomials in b, eta, r, p, alpha, v, H, j.
N = 8


def const(x):
    x = F(x)
    return {(0,) * N: x} if x else {}


def variable(i):
    e = [0] * N
    e[i] = 1
    return {tuple(e): F(1)}


def add(*polys):
    out = {}
    for poly in polys:
        for e, x in poly.items():
            out[e] = out.get(e, F(0)) + x
    return {e: x for e, x in out.items() if x}


def scale(poly, c):
    return {e: x * c for e, x in poly.items() if x * c}


def mul(a, b):
    out = {}
    for ea, xa in a.items():
        for eb, xb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, F(0)) + xa * xb
    return {e: x for e, x in out.items() if x}


def sq(a):
    return mul(a, a)


def diff(poly, i):
    out = {}
    for e, x in poly.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = x * e[i]
    return out


def evaluate(poly, values):
    total = F(0)
    for e, x in poly.items():
        term = x
        for value, power in zip(values, e):
            term *= value ** power
        total += term
    return total


def frac(q):
    return str(q.numerator) + '/' + str(q.denominator)


def algebra():
    checks = []

    def same(label, a, b):
        require(a == b, label)
        checks.append(label)

    def positive_coefficients(label, poly):
        require(all(x >= 0 for x in poly.values()), label)
        checks.append(label)

    b, eta, r, p, alpha, v, H, j = [variable(i) for i in range(N)]
    one = const(1)
    same('block_cauchy_schwarz_identity',
         add(mul(add(sq(r), sq(p)), add(sq(alpha), sq(v))),
             scale(sq(add(mul(r, alpha), mul(p, v))), -1)),
         sq(add(mul(r, v), scale(mul(p, alpha), -1))))
    same('l1_dominates_l2_squared',
         add(sq(add(r, p)), scale(add(sq(r), sq(p)), -1)),
         scale(mul(r, p), 2))
    same('time_reflection_isometry', add(sq(scale(alpha, -1)), sq(v)),
         add(sq(alpha), sq(v)))
    same('cut_lorentzian_cancellation',
         add(sq(add(one, scale(b, -1))), scale(sq(b), -1)),
         add(one, scale(b, -2)))
    error = add(sq(add(one, b)), sq(b))
    same('cut_auxiliary_norm_expansion', error,
         add(one, scale(b, 2), scale(sq(b), 2)))
    lower = add(one, scale(b, -2), scale(mul(eta, error), -1))
    upper = add(mul(H, sq(add(one, b))), mul(eta, error))
    positive_coefficients('lower_decreases_with_b', scale(diff(lower, 0), -1))
    positive_coefficients('lower_decreases_with_eta', scale(diff(lower, 1), -1))
    positive_coefficients('upper_increases_with_b', diff(upper, 0))
    positive_coefficients('upper_increases_with_eta', diff(upper, 1))
    positive_coefficients('upper_increases_with_cosh_squared', diff(upper, 6))
    factorial_ratio = add(mul(add(scale(j, 2), const(2)),
                             add(scale(j, 2), one)), const(-2))
    same('factorial_majorant_induction_identity', factorial_ratio,
         add(scale(sq(j), 4), scale(j, 6)))
    positive_coefficients('factorial_majorant_induction_nonnegative', factorial_ratio)
    geometric_denominator = F(1) - F(1, 4) ** 2 / 2
    same('cosh_quarter_majorant', 1 / geometric_denominator, F(32, 31))
    cosh2 = F(32, 31) ** 2
    same('cosh_squared_majorant', cosh2, F(1024, 961))
    values = [F(0)] * N
    values[0], values[1], values[6] = F(1, 8), F(1, 8), cosh2
    lcut = evaluate(lower, values)
    ucut = evaluate(upper, values)
    same('lower_cut_endpoint', lcut, F(151, 256))
    lower_g = lcut / cosh2
    same('lower_normalized_endpoint', lower_g, F(145111, 262144))
    same('upper_normalized_endpoint', ucut, F(371177, 246016))
    same('upper_split_expression', ucut, F(1296, 961) + F(41, 256))
    require(lower_g > F(1, 2), 'strict lower spectrum')
    checks.append('strict_lower_spectrum')
    require(ucut < 2, 'strict upper spectrum')
    checks.append('strict_upper_spectrum')
    same('lower_spectrum_margin', lower_g - F(1, 2), F(14039, 262144))
    same('upper_spectrum_margin', 2 - ucut, F(120855, 246016))
    # The following are exact algebra cases, not physical input values.
    guard_cases = [
        ((F(0), F(1), F(0), F(0)), True),
        ((F(1), F(1), F(1, 4), F(1, 8)), True),
        ((F(1), F(1), F(1, 4), F(1, 4)), False),
        ((F(1), F(1), F(1, 2), F(1, 8)), False),
        ((F(-1), F(1), F(0), F(0)), False),
        ((F(1), F(0), F(0), F(0)), False),
        ((F(1), F(1), F(-1), F(0)), False),
    ]
    for index, (data, expected) in enumerate(guard_cases):
        a, rho, u, pp = data
        actual = a >= 0 and rho > 0 and u >= 0 and pp >= 0 and a * (rho*u+pp) < F(1, 2)
        same('rational_strict_guard_case_' + str(index), actual, expected)
    for index in range(21):
        same('factorial_coefficient_case_' + str(index),
             math.factorial(2*index) >= 2**index, True)
    return {
        'check_count': len(checks), 'checks': checks,
        'exact_constants': {
            'lower_cut': frac(lcut), 'cosh_squared_upper': frac(cosh2),
            'lower_G': frac(lower_g), 'upper_G': frac(ucut),
            'lower_G_minus_half': frac(lower_g-F(1, 2)),
            'two_minus_upper_G': frac(2-ucut),
        },
        'continuous_proof_required': [
            'full_cylinder_flow_existence_and_retention',
            'physical_fixed_gauge_and_complete_shell_fields',
            'uniform_velocity_and_Jacobi_parameter_bounds',
            'sine_concavity_and_product_minimizing_geodesic',
            'integration_on_cylinder_and_compact_embedding',
            'metric_gap_and_tensor_norm_hypotheses',
            'infinite_cosh_series_majorization',
        ],
        'scope': 'finite_exact_algebra_not_PDE_or_effective_amplitude',
    }


def output_path(raw):
    p = Path(raw)
    require(p.is_absolute(), 'output must be absolute')
    require('..' not in p.parts, 'output must not contain dot-dot')
    for ancestor in [p] + list(p.parents):
        require(not ancestor.is_symlink(), 'output path contains symlink')
    require(p.parent.is_dir(), 'output parent must exist')
    require(not p.exists(), 'output must be fresh')
    require(p != ROOT and ROOT not in p.parents, 'output must be outside closed capsule')
    require(p.name not in ('', '.', '..'), 'invalid output filename')
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        target = output_path(args.output)
        auth = authenticate()
        report = {'schema': 'hdblast-observer-domain-exact-v1', 'status': 'pass',
                  'authentication': auth, 'algebra': algebra()}
        data = canonical(report)
        # Exclusive, no-follow final creation; validated ancestors are not created.
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, 'O_NOFOLLOW', 0)
        fd = os.open(str(target), flags, 0o600)
        with os.fdopen(fd, 'wb') as out:
            out.write(data)
        sys.stdout.buffer.write(data)
        return 0
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
        print('observer-domain check rejected: ' + str(exc), file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
