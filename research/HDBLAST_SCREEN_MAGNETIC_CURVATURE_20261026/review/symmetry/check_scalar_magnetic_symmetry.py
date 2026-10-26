#!/usr/bin/env python3
"""Independent exact S3 scalar-square magnetic-TT witness; no frozen imports.

This checks angular admissibility, not the coefficient of physical Q_H.
All sphere polynomials are reduced in Q[n0,n1,n2,n3]/(sum n_i^2-1).
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations, product
from pathlib import Path

I = range(4)
O = (0, 0, 0, 0)
HERE = Path(__file__).resolve().parent
CAPSULE = 'HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026'

def guard_output(output):
    """Validate a fresh external destination before the proof calculation."""
    output = Path(output)
    if output.is_symlink():
        raise ValueError('OUTPUT_SYMLINK')
    if output.exists():
        raise ValueError('OUTPUT_EXISTS')
    absolute = output.absolute()
    if any(parent.is_symlink() for parent in absolute.parents):
        raise ValueError('OUTPUT_ANCESTOR_SYMLINK')
    target = output.resolve()
    if target == HERE or HERE in target.parents:
        raise ValueError('OUTPUT_INSIDE_SOURCE_TREE')
    roots = [p for p in (target, *target.parents)
             if p.name == CAPSULE and p.parent.name == 'research']
    roots += [p for p in (HERE, *HERE.parents)
              if p.name == CAPSULE and p.parent.name == 'research'
              and (target == p or p in target.parents)]
    if roots:
        raise ValueError('OUTPUT_INSIDE_PACKAGED_CAPSULE')
    if not target.parent.is_dir():
        raise ValueError('OUTPUT_PARENT_MISSING')
    return output

def add(*terms):
    result = {}
    for term in terms:
        for key, value in term.items():
            result[key] = result.get(key, F(0)) + value
    return {key: value for key, value in result.items() if value}

def scale(value, poly):
    return {key: value * c for key, c in poly.items() if value * c}

@lru_cache(None)
def reduce_monomial(key):
    if key[3] < 2:
        return {key: F(1)}
    base = list(key)
    base[3] -= 2
    result = reduce_monomial(tuple(base)).copy()
    for i in range(3):
        other = base.copy()
        other[i] += 2
        result = add(result, scale(-1, reduce_monomial(tuple(other))))
    return result

def reduce(poly):
    return add(*(scale(value, reduce_monomial(key)) for key, value in poly.items()))

def mul(left, right):
    result = {}
    for a, v in left.items():
        for b, w in right.items():
            key = tuple(x + y for x, y in zip(a, b))
            result[key] = result.get(key, F(0)) + v * w
    return reduce(result)

def partial(poly, index):
    result = {}
    for key, value in poly.items():
        if key[index]:
            other = list(key)
            other[index] -= 1
            result[tuple(other)] = value * key[index]
    return result

ONE = {O: F(1)}
N = []
for i in I:
    key = [0] * 4
    key[i] = 1
    N.append({tuple(key): F(1)})
P = {(i, j): add(ONE if i == j else {}, scale(-1, mul(N[i], N[j])))
     for i, j in product(I, repeat=2)}

def covariant_derivative(tensor):
    # For a tangent covariant tensor this is the fully projected ambient
    # derivative. The n_i terms follow by differentiating tangency.
    result = {}
    for key, poly in tensor.items():
        for k in I:
            terms = [mul(P[k, a], partial(poly, a)) for a in I]
            for slot, i in enumerate(key):
                other = list(key)
                other[slot] = k
                terms.append(mul(N[i], tensor[tuple(other)]))
            result[(k,) + key] = add(*terms)
    return result

def tadd(*terms):
    return {key: add(*(term[key] for term in terms)) for key in terms[0]}

def tscale(value, tensor):
    return {key: scale(value, poly) for key, poly in tensor.items()}

def equal(left, right):
    return all(not reduce(add(left[key], scale(-1, right[key]))) for key in left)

def trace(tensor):
    return add(*(tensor[i, i] for i in I))

def divergence(tensor):
    derivative = covariant_derivative(tensor)
    return {(j,): add(*(derivative[i, i, j] for i in I)) for j in I}

def sign(sequence):
    return (-1) ** sum(sequence[i] > sequence[j] for i in I for j in range(i + 1, 4))

EPS = {key: sign(key) for key in permutations(I)}
SCREEN_EPS = {(i, k, l): add(*(scale(EPS.get((m, i, k, l), 0), N[m]) for m in I))
              for i, k, l in product(I, repeat=3)}

def curl(tensor):
    derivative = covariant_derivative(tensor)
    return {(i, j): scale(F(1, 2), add(*(
        add(mul(SCREEN_EPS[i, k, l], derivative[k, j, l]),
            mul(SCREEN_EPS[j, k, l], derivative[k, i, l]))
        for k, l in product(I, repeat=2))))
        for i, j in product(I, repeat=2)}

def scalar_gradient_square(diagonal):
    y = add(*(scale(diagonal[i], mul(N[i], N[i])) for i in I))
    gradient = {(i,): scale(2, add(scale(diagonal[i], N[i]), scale(-1, mul(y, N[i]))))
                for i in I}
    return y, {(i, j): mul(gradient[i,], gradient[j,]) for i, j in product(I, repeat=2)}

def make_frozen_tt_definitions():
    # This C is the displayed angular_exact.py definition, transcribed as
    # rational Weyl data. Self-duality labels are not a curl-eigenvalue premise.
    diagonal = {(0, 1): F(-2, 3), (2, 3): F(-2, 3), (0, 2): F(1, 3),
                (0, 3): F(1, 3), (1, 2): F(1, 3), (1, 3): F(1, 3)}
    C = {key: F(0) for key in product(I, repeat=4)}
    for (i, j), value in diagonal.items():
        C[i, j, i, j] = C[j, i, j, i] = value
        C[i, j, j, i] = C[j, i, i, j] = -value
    star = {key: sum((F(EPS.get((key[0], key[1], a, b), 0), 2) * C[a, b, key[2], key[3]]
                     for a, b in product(I, repeat=2)), F(0)) for key in C}
    result = []
    for orientation in (1, -1):
        W = {key: (C[key] + orientation * star[key]) / 2 for key in C}
        result.append({(i, j): add(*(scale(W[i, a, j, b], mul(N[a], N[b]))
                                     for a, b in product(I, repeat=2)))
                       for i, j in product(I, repeat=2)})
    return result

def integrate(poly):
    result = F(0)
    for key, value in poly.items():
        if any(power % 2 for power in key):
            continue
        moment = F(1)
        for power in key:
            for j in range(power // 2):
                moment *= 2 * j + 1
        for j in range(sum(key) // 2):
            moment /= 4 + 2 * j
        result += value * moment
    return result

def inner(left, right):
    return sum((integrate(mul(left[key], right[key])) for key in left), F(0))

def reflection(tensor):
    return {key: {powers: value * (-1) ** (powers[3] + key.count(3))
                  for powers, value in poly.items()} for key, poly in tensor.items()}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mutation', choices=['magnetic_sign', 'electric_half', 'reflection_slots'])
    args = parser.parse_args()
    destination = guard_output(args.output)
    plus, minus = make_frozen_tt_definitions()
    y, source = scalar_gradient_square((1, -1, 0, 0))
    magnetic = curl(source)
    expected = tadd(tscale(-2, plus), tscale(2, minus))
    if args.mutation == 'magnetic_sign':
        expected = tscale(-1, expected)
    if not equal(magnetic, expected):
        raise ValueError('COMPLETE_MAGNETIC_TT_IDENTITY')
    norms = [inner(tensor, tensor) for tensor in (plus, minus)]
    electric = [inner(source, tensor) / norm for tensor, norm in zip((plus, minus), norms)]
    if args.mutation == 'electric_half':
        electric[0] /= 2
    if electric != [F(-2, 3), F(-2, 3)]:
        raise ValueError('ELECTRIC_TT_PROJECTION')
    if trace(magnetic) or any(divergence(magnetic).values()):
        raise ValueError('MAGNETIC_TRACE_DIVERGENCE')
    reflected = reflection(plus)
    if args.mutation == 'reflection_slots':
        reflected = {key: {powers: value * (-1) ** powers[3] for powers, value in poly.items()}
                     for key, poly in plus.items()}
    if not equal(reflected, minus):
        raise ValueError('REFLECTION_EXCHANGE_INCLUDING_SLOTS')
    if not equal(reflection(magnetic), tscale(-1, magnetic)):
        raise ValueError('MAGNETIC_PSEUDOTENSOR_REFLECTION')
    _, axisymmetric_source = scalar_gradient_square((3, -1, -1, -1))
    if any(curl(axisymmetric_source).values()):
        raise ValueError('AXISYMMETRIC_COMPARISON')
    report = {
        'schema': 'HDBLAST_INDEPENDENT_SCALAR_MAGNETIC_SYMMETRY_V1',
        'status': 'PASS_EXACT_ALL_POINT_ANGULAR_IDENTITIES',
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'arithmetic': 'independent Python Fraction quotient-polynomial arithmetic; no frozen imports; no sampling',
        'actual_seed': 'Y=n0^2-n1^2; A=diag(1,-1,0,0)',
        'SO3_axisymmetric_seed_comparison': 'A=diag(3,-1,-1,-1)',
        'spin4_symmetric_square': '(0,0)+(1,1)+(2,2)+(2,0)+(0,2); dimensions 1+9+25+5+5=45',
        'TT_norms_squared': list(map(str, norms)),
        'electric_TT_projection_DY_outer_DY': list(map(str, electric)),
        'magnetic_full_tensor_curl_DY_outer_DY': ['-2', '2'],
        'magnetic_trace_and_divergence': 'exact zero',
        'orientation_convention': 'epsilon_ikl=sum_m epsilon_mikl n_m; epsilon_0123=+1',
        'reflection': 'diag(1,1,1,-1) fixes Y, exchanges Tplus/Tminus, and sends magnetic witness to its negative',
        'axisymmetric_witness': 'exact zero',
        'limits': ['This proves a quadratic scalar angular jet may generate magnetic TT.',
                   'This does not compute the geometric Q_H coefficient for the actual first metric.',
                   'A nonzero allowed witness does not exclude cancellations among the full physical terms.']
    }
    with destination.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'output_name': destination.name}))

if __name__ == '__main__':
    main()
