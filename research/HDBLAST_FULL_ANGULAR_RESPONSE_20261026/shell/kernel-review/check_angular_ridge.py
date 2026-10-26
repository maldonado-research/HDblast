#!/usr/bin/env python3
"""Independent exact polygon/endpoint review of the new full23 ridge kernel.

The oracle integrates monomials on clipped physical retarded diamonds. It
never invokes a producer kernel to obtain expected integral values.
"""
import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import sys
from types import ModuleType

sys.dont_write_bytecode = True
DIM = 23
ZERO = (Q(0),)*DIM
KERNEL_PIN = 'd6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f'


def need(ok, label):
    if not ok:
        raise ValueError(label)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add(*rows):
    return tuple(sum(items, Q(0)) for items in zip(*rows))


def scale(c, row):
    return tuple(c*v for v in row)


def mv(matrix, vector):
    return tuple(sum((a*b for a, b in zip(row, vector)), Q(0)) for row in matrix)


def dot(a, b):
    return sum((u*v for u, v in zip(a, b)), Q(0))


def polyval(coeffs, t):
    out = ZERO
    for row in reversed(coeffs):
        out = add(scale(t, out), row)
    return out


def integrate_poly(coeffs, boundary=Q(0), boundary_value=ZERO):
    out = [ZERO]+[scale(Q(1, j+1), row) for j, row in enumerate(coeffs)]
    out[0] = add(boundary_value, scale(-1, polyval(out, boundary)))
    return tuple(out)


class Profiles:
    """Independent piecewise polynomials, with normalized iterated primitives."""
    def __init__(self, pieces):
        self.pieces = tuple((Q(a), Q(b), tuple(tuple(map(Q, row)) for row in poly))
                            for a, b, poly in pieces)
        self.cache = {0: self.pieces}

    def primitive_pieces(self, order):
        if order not in self.cache:
            previous = self.primitive_pieces(order-1)
            positive = [piece for piece in previous if piece[0] >= 0]
            negative = list(reversed([piece for piece in previous if piece[1] <= 0]))
            need(len(positive)+len(negative) == len(previous), 'ORACLE_SPLIT_AT_ORIGIN')
            pieces = []
            for sign, rows in ((1, positive), (-1, negative)):
                value = ZERO
                for a, b, poly in rows:
                    anchor, end = (a, b) if sign == 1 else (b, a)
                    primitive = integrate_poly(poly, anchor, value)
                    value = polyval(primitive, end)
                    pieces.append((a, b, primitive))
            self.cache[order] = tuple(sorted(pieces))
        return self.cache[order]

    def __call__(self, t, order=0):
        t = Q(t)
        for a, b, poly in self.primitive_pieces(order):
            if a <= t <= b:
                return polyval(poly, t)
        raise ValueError('ORACLE_PROFILE_FOOTPRINT')


def p1(nodes, values, parity):
    positive = []
    for a, b, va, vb in zip(nodes, nodes[1:], values, values[1:]):
        slope = scale(1/(b-a), add(vb, scale(-1, va)))
        positive.append((a, b, (add(va, scale(-a, slope)), slope)))
    negative = [(-b, -a, tuple(tuple(sign*((-1)**j)*v for sign, v in zip(parity, row))
                               for j, row in enumerate(poly)))
                for a, b, poly in positive]
    return Profiles(sorted(negative+positive))


def conv(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            out[i+j] += u*v
    return out


def power(poly, n):
    out = [Q(1)]
    for _ in range(n):
        out = conv(out, poly)
    return out


@lru_cache(maxsize=None)
def moment(polygon, i, j):
    """Green's theorem: integral s^i*v^j = contour s^(i+1)*v^j dv/(i+1)."""
    if len(polygon) < 3:
        return Q(0)
    result = Q(0)
    for a, b in zip(polygon, polygon[1:]+polygon[:1]):
        ds, dv = b[0]-a[0], b[1]-a[1]
        poly = conv(power([a[0], ds], i+1), power([a[1], dv], j))
        result += dv*sum((value/Q(k+1) for k, value in enumerate(poly)), Q(0))/(i+1)
    return result


def clip(polygon, cs, cv, constant):
    out = []
    if not polygon:
        return ()
    for a, b in zip(polygon, polygon[1:]+polygon[:1]):
        u, v = cs*a[0]+cv*a[1]+constant, cs*b[0]+cv*b[1]+constant
        inside_a, inside_b = u >= 0, v >= 0
        if inside_a:
            out.append(a)
        if inside_a != inside_b:
            t = u/(u-v)
            out.append((a[0]+t*(b[0]-a[0]), a[1]+t*(b[1]-a[1])))
    return tuple(out)


def polygon_profile(profile, cs, cv, polygon):
    out = ZERO
    for a, b, coefficients in profile.pieces:
        cell = clip(clip(polygon, cs, cv, -a), -cs, -cv, b)
        for n, row in enumerate(coefficients):
            coefficient = sum((Q(comb(n, i))*cs**i*cv**(n-i)*moment(cell, i, n-i)
                               for i in range(n+1)), Q(0))
            out = add(out, scale(coefficient, row))
    return out


def line_profile(profile, alpha, beta, left, right):
    """Integrate profile(alpha*s+beta) ds directly on each polynomial segment."""
    if right <= left:
        return ZERO
    if alpha == 0:
        return scale(right-left, profile(beta))
    out = ZERO
    for a, b, coefficients in profile.pieces:
        sa, sb = (a-beta)/alpha, (b-beta)/alpha
        l, u = max(left, min(sa, sb)), min(right, max(sa, sb))
        if u <= l:
            continue
        for n, row in enumerate(coefficients):
            coefficient = sum((Q(comb(n, i))*alpha**i*beta**(n-i)
                               *(u**(i+1)-l**(i+1))/(i+1) for i in range(n+1)), Q(0))
            out = add(out, scale(coefficient, row))
    return out


def direct_kernel(F, G, x, r, parity, blind=False):
    """Physical source first, then output reflection; no producer primitives."""
    x, r = Q(x), Q(r)
    polygon = ((Q(0), r-x), (x, r), (Q(0), r+x))
    positive = clip(polygon, 0, 1, 0)
    negative = clip(polygon, 0, -1, 0)
    physical = add(polygon_profile(F, 1, 1, positive), polygon_profile(G, -1, 1, positive))
    if blind:
        reflected = add(polygon_profile(F, 1, 1, negative), polygon_profile(G, -1, 1, negative))
    else:
        z = add(polygon_profile(F, 1, -1, negative), polygon_profile(G, -1, -1, negative))
        reflected = tuple(sign*v for sign, v in zip(parity, z))
    value = scale(Q(1, 2), add(physical, reflected))
    a, q = r+x, r-x
    upper = add(line_profile(F, 0, a, 0, x), line_profile(G, -2, a, 0, x))
    if q >= 0 or blind:
        lower = add(line_profile(F, 2, q, 0, x), line_profile(G, 0, q, 0, x))
    else:
        t = -q
        raw = add(line_profile(F, 0, t, 0, t), line_profile(G, -2, t, 0, t))
        lower = add(tuple(sign*v for sign, v in zip(parity, raw)),
                    line_profile(F, 2, q, t, x), line_profile(G, 0, q, t, x))
    return {'value': value, 'dx': scale(Q(1, 2), add(upper, lower)),
            'dr': scale(Q(1, 2), add(upper, scale(-1, lower)))}


def combine_profiles(A, B, D, u, du, v):
    F, G = [], []
    vp = v.primitive_pieces(1)
    for up, dp, vv, pv in zip(u.pieces, du.pieces, v.pieces, vp):
        a, b = up[:2]
        need(all(piece[:2] == (a, b) for piece in (dp, vv, pv)), 'ORACLE_ALIGNED_CELLS')
        fcoeffs, gcoeffs = [], []
        for j in range(max(len(piece[2]) for piece in (up, dp, vv, pv))):
            U, DU, V, PV = [piece[2][j] if j < len(piece[2]) else ZERO for piece in (up, dp, vv, pv)]
            # Expand from A*free_x+B*free_r+D*free, independently of CoupledRidges.F/G.
            fcoeffs.append(scale(Q(1, 2), add(mv(A, add(DU, V)), mv(B, add(DU, V)),
                                                  mv(D, add(U, PV)))))
            gcoeffs.append(scale(Q(1, 2), add(mv(A, add(scale(-1, DU), V)),
                                                  mv(B, add(DU, scale(-1, V))),
                                                  mv(D, add(U, scale(-1, PV))))))
        F.append((a, b, tuple(fcoeffs)))
        G.append((a, b, tuple(gcoeffs)))
    return Profiles(F), Profiles(G)


def free(u, du, v, x, r):
    a, b = r+x, r-x
    return {'value': scale(Q(1, 2), add(u(a), u(b), v(a, 1), scale(-1, v(b, 1)))),
            'dx': scale(Q(1, 2), add(du(a), scale(-1, du(b)), v(a), v(b))),
            'dr': scale(Q(1, 2), add(du(a), du(b), v(a), scale(-1, v(b))))}


def direct_feedback(F, G, u, v, x, r, M, parity, blind=False):
    out = direct_kernel(F, G, x, r, parity, blind)
    t = x-r
    if t > 0:
        projection = lambda row: tuple(value if sign == 1 else Q(0) for sign, value in zip(parity, row))
        g = mv(M, projection(add(u(t), v(t, 1))))
        pg = mv(M, projection(add(u(t, 1), v(t, 2))))
        out = {'value': add(out['value'], scale(-1, pg)), 'dx': add(out['dx'], scale(-1, g)),
               'dr': add(out['dr'], g)}
    return out


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--work-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-error', choices=('wrong-reflection', 'missing-one-D-row',
                                                  'commute-matrices', 'omit-PN-projection'))
    args = parser.parse_args()
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    auth = load('new_shell_authentication', args.work_root/'shell/check_full23_shell.py')
    roots = {key: args.repo_root.resolve()/'research'/label for key, label in auth.LABELS.items()}
    need(all(root not in args.output.resolve().parents for root in roots.values()), 'EXTERNAL_OUTPUT')
    before = {key: auth.authenticate(root, auth.PINS[key]) for key, root in roots.items()}
    need(sha(roots['latest']/'balanced/exact_boundary.py') == auth.HELPER_PIN, 'HELPER_PIN')
    old_text = (roots['current']/'reflected_ridge.py').read_text()
    kernel_path = args.work_root/'compute/angular_ridge.py'
    need(sha(kernel_path) == KERNEL_PIN, 'NEW_KERNEL_PIN')
    transformation = json.loads((args.work_root/'compute/KERNEL_TRANSFORMATION.json').read_text())
    need(sha(roots['current']/'reflected_ridge.py') == transformation['source_sha256'], 'OLD_SOURCE_PIN')
    reconstructed = old_text
    for row in transformation['changes']:
        need(reconstructed.count(row['old']) == row['occurrences'], 'DIMENSION_TRANSFORM_COUNT')
        reconstructed = reconstructed.replace(row['old'], row['new'])
    need(reconstructed == kernel_path.read_text() and len(transformation['changes']) == 10,
         'DIMENSION_ONLY_TRANSFORMATION')
    text = kernel_path.read_text()
    if args.inject_error == 'wrong-reflection':
        text = text.replace('    original_F, original_G = F, G',
                            '    parity = tuple(1 for _ in PARITY)\n    original_F, original_G = F, G')
    elif args.inject_error == 'missing-one-D-row':
        text = text.replace('if sign == -1 for j in range(DIM)',
                            'if sign == -1 and i != 16 for j in range(DIM)')
    elif args.inject_error == 'commute-matrices':
        text = text.replace('out = reflected_kernel(self.F, self.G, x, r)',
                            'out = _blind_reflected_kernel(self.F, self.G, x, r)')
    elif args.inject_error == 'omit-PN-projection':
        text = text.replace('project = lambda row: tuple(v if s == 1 else Q(0)\n                                        for s, v in zip(PARITY, row))',
                            'project = lambda row: tuple(row)')
    need(not args.inject_error or text != kernel_path.read_text(), 'REVIEW_MUTATION_NOT_APPLIED')
    k = ModuleType('reviewed_new_full23_kernel')
    exec(compile(text, str(kernel_path)+':review-mutant' if args.inject_error else str(kernel_path), 'exec'), k.__dict__)
    op = load('authenticated_frozen_full23_operator', roots['G']/'projection/inputs/linear_operator.py')
    expected_parity = tuple(-1 if name.startswith(('L:', 'v:')) else 1 for name in op.LABELS)
    need(k.LABELS == tuple(op.LABELS) and k.GRAM == tuple(op.gram()) and k.PARITY == expected_parity,
         'FULL23_LABEL_GRAM_PARITY')
    parity = expected_parity
    counts = {}
    samples = []

    def equal(actual, expected, label):
        need(actual == expected, label)
        counts[label.split(':')[0]] = counts.get(label.split(':')[0], 0)+1

    # Full23 polynomial forcing, all slots active and both reflection classes.
    Fcoeff = [tuple(Q((-1)**(i+j)*(i+2)*(j+1), 7+j) for i in range(DIM)) for j in range(5)]
    Gcoeff = [tuple(Q((i-9)*(j+2), 11+j) for i in range(DIM)) for j in range(5)]
    F = Profiles([(-4, 0, Fcoeff), (0, 4, Fcoeff)])
    G = Profiles([(-4, 0, Gcoeff), (0, 4, Gcoeff)])
    points = [(Q(0), Q(0)), (Q(0), Q(7, 4)), (Q(1, 2), Q(0)),
              (Q(5, 4), Q(0)), (Q(3, 2), Q(1, 4)), (Q(1, 4), Q(9, 8)),
              (Q(1, 2), Q(1, 2)), (Q(5, 7), Q(2, 9)), (Q(17, 24), Q(17, 24)),
              (Q(1, 8), Q(3, 2)), (Q(11, 8), Q(3, 8))]
    for x, r in points:
        actual, expected = k.reflected_kernel(F, G, x, r), direct_kernel(F, G, x, r, parity)
        for key in ('value', 'dx', 'dr'):
            equal(actual[key], expected[key], 'POLYNOMIAL_DIAMOND:'+key)
        if x == 0:
            equal(actual, dict.fromkeys(('value', 'dx', 'dr'), ZERO), 'POLYNOMIAL_INITIAL_ZERO')
        if r == 0:
            equal(tuple(actual[key][j] for key in ('value', 'dx') for j in range(DIM) if parity[j] == -1),
                  (Q(0),)*10, 'ALL_FIVE_POLYNOMIAL_DIRICHLET')
            equal(tuple(actual['dr'][j] for j in range(DIM) if parity[j] == 1),
                  (Q(0),)*18, 'ALL_EIGHTEEN_POLYNOMIAL_NEUMANN')
    # Exact one-sided polynomial continuations to r=x certify all C1 joins.
    join = Q(17, 24)
    step = Q(1, 4096)
    for sign in (-1, 1):
        values = [k.reflected_kernel(F, G, join, join+sign*j*step) for j in range(1, 8)]
        weights = [Q((-1)**(j-1)*comb(7, j)) for j in range(1, 8)]
        target = k.reflected_kernel(F, G, join, join)
        for key in ('value', 'dx', 'dr'):
            limit = add(*(scale(w, row[key]) for w, row in zip(weights, values)))
            equal(limit, target[key], 'EXACT_C1_RIDGE_JOIN:'+str(sign)+':'+key)

    nodes = tuple(map(Q, (0, 1, Q(5, 4), Q(3, 2), Q(7, 4))))
    values = {}
    for name, phase in (('u', 1), ('du', 2), ('v', 3)):
        values[name] = [ZERO, ZERO]+[tuple(Q((-1)**(i+j+phase)*(i+phase)*(j+2), 13+phase)
                                         for i in range(DIM)) for j in range(3)]
    A = tuple(tuple(Q((-1)**(i+j)*(i+1), 19+j) for j in range(DIM)) for i in range(DIM))
    B = tuple(tuple(Q((i+2)*(j+1), 53) for j in range(DIM)) for i in range(DIM))
    D = tuple(tuple(Q(i-j+1, 71) for j in range(DIM)) for i in range(DIM))
    M = tuple(tuple(Q((i+1)*(j+1), 101) if parity[i] == 1 else Q(0)
                    for j in range(DIM)) for i in range(DIM))
    u, du, v = [p1(nodes, values[name], tuple(-s for s in parity) if name == 'du' else parity)
                for name in ('u', 'du', 'v')]
    ku, kd, kv = [k.PiecewiseLinear(nodes, values[name], tuple(-s for s in parity) if name == 'du' else parity)
                  for name in ('u', 'du', 'v')]
    for own, producer in ((u, ku), (du, kd), (v, kv)):
        for t in tuple(nodes)+tuple((a+b)/2 for a, b in zip(nodes, nodes[1:])):
            for sign in (-1, 1):
                for order in range(6):
                    equal(producer.primitive(sign*t, order), own(sign*t, order), 'P1_NORMALIZED_PRIMITIVES')
    source = k.CoupledRidges(A, B, D, ku, kd, kv)
    cF, cG = combine_profiles(A, B, D, u, du, v)
    k._blind_reflected_kernel = lambda f, g, x, r: direct_kernel(cF, cG, x, r, parity, blind=True)
    for t in tuple(nodes)+tuple(-x for x in nodes)+tuple((a+b)/2 for a, b in zip(nodes, nodes[1:])):
        for order in range(5):
            equal(source.F(t, order), cF(t, order), 'COUPLED_F_PRIMITIVES')
            equal(source.G(t, order), cG(t, order), 'COUPLED_G_PRIMITIVES')
    witnesses = []
    physical_points = [(x, r) for x, r in points
                       if x <= Q(3, 2) and r >= max(Q(0), x-Q(5, 4)) and x+r <= Q(7, 4)]
    physical_points += [(Q(9, 8), Q(0)), (Q(3, 4), Q(3, 4)), (Q(1, 2), Q(5, 4)),
                        (Q(5, 4), Q(1, 4)), (Q(1, 8), Q(1, 2)),
                        (Q(9, 8), Q(1, 32)), (Q(5, 4), Q(1, 16)),
                        (Q(11, 8), Q(1, 8))]
    for x, r in physical_points:
        field = free(u, du, v, x, r)
        forcing = add(mv(A, field['dx']), mv(B, field['dr']), mv(D, field['value']))
        equal(source.physical_source(x, r), forcing, 'COMPLETE_COUPLED_PHYSICAL_SOURCE')
        expected = direct_feedback(cF, cG, u, v, x, r, M, parity)
        actual = source.feedback(x, r, M)
        for key in ('value', 'dx', 'dr'):
            equal(actual[key], expected[key], 'COUPLED_PHYSICAL_POLYGON:'+key)
        if r == 0:
            h = tuple(z if sign == 1 else Q(0) for sign, z in zip(parity, add(u(x), v(x, 1))))
            equal(tuple(actual[key][i] for key in ('value', 'dx') for i in range(DIM) if parity[i] == -1),
                  (Q(0),)*10, 'ALL_FIVE_COUPLED_DIRICHLET')
            equal(tuple(actual['dr'][i] for i in range(DIM) if parity[i] == 1),
                  tuple(mv(M, h)[i] for i in range(DIM) if parity[i] == 1), 'ALL_EIGHTEEN_ROBIN_LIFT')
            total_value = add(field['value'], actual['value'])
            total_dr = add(field['dr'], actual['dr'])
            defect = add(total_dr, scale(-1, mv(M, total_value)))
            # PN only: Dirichlet normal derivatives are intentionally unconstrained.
            equal(tuple(defect[i] for i in range(DIM) if parity[i] == 1),
                  tuple(-z for i, z in enumerate(mv(M, actual['value'])) if parity[i] == 1),
                  'FINITE_CANDIDATE_ROBIN_DEFECT')
            # Change the exact displacement while retaining the P1 feedback.
            delta = tuple(Q(i+1, 37)*x*x if sign == 1 else Q(0) for i, sign in enumerate(parity))
            actual_exact_displacement = add(total_value, delta)
            mismatch_defect = add(total_dr, scale(-1, mv(M, actual_exact_displacement)))
            expected_defect = scale(-1, mv(M, add(actual['value'], delta)))
            equal(tuple(mismatch_defect[i] for i in range(DIM) if parity[i] == 1),
                  tuple(expected_defect[i] for i in range(DIM) if parity[i] == 1),
                  'RATIONAL_VERSUS_P1_DISPLACEMENT_ROBIN_DEFECT')
        if x == 0 or x+r <= 1:
            equal(actual, dict.fromkeys(('value', 'dx', 'dr'), ZERO), 'COUPLED_INITIAL_ZERO_WEDGE')
        if r > 0 and x > r:
            wrong = direct_feedback(cF, cG, u, v, x, r, M, parity, blind=True)
            if wrong != expected:
                witnesses.append({'x': x, 'r': r,
                                  'correct_value_K0': expected['value'][0],
                                  'blind_commuted_value_K0': wrong['value'][0]})
        samples.append({'x': x, 'r': r, 'oracle_agrees_all23': True})
    need(witnesses, 'NONCOMMUTING_REFLECTION_MUST_HAVE_NONZERO_WITNESS')
    # Polynomial continuation across physical ridge/node/front interfaces.
    # Each one-sided value is a polynomial of degree <=4; seven samples are
    # deliberately more than the required interpolation degree.
    joins = [('r=x', Q(17, 24), Q(17, 24)),
             ('r+x=1', Q(1, 4), Q(3, 4)),
             ('r+x=5/4', Q(1, 4), Q(1)),
             ('r-x=1', Q(1, 6), Q(7, 6)),
             ('x-r=1', Q(7, 6), Q(1, 6))]
    for name, x, r in joins:
        target = source.feedback(x, r, M)
        for sign in (-1, 1):
            values_at_side = [source.feedback(x, r+sign*j*step, M) for j in range(1, 8)]
            weights = [Q((-1)**(j-1)*comb(7, j)) for j in range(1, 8)]
            for key in ('value', 'dx', 'dr'):
                continuation = add(*(scale(w, row[key]) for w, row in zip(weights, values_at_side)))
                equal(continuation, target[key], 'COUPLED_EXACT_C1_JOIN:'+name+':'+str(sign)+':'+key)
    # Every individual D row is protected, including D-column entries that the
    # old mean-only guard did not inspect. Five independent hard cases.
    rejected = []
    bad_matrix_cases = []
    for i, sign in enumerate(parity):
        if sign != -1:
            continue
        for column_class, j in (('PN-input', 0), ('PD-input', i)):
            bad = [list(row) for row in M]
            bad[i][j] = Q(1)
            try:
                source.feedback(Q(9, 8), Q(0), bad)
            except ValueError as error:
                need('no Dirichlet forcing' in str(error), 'DIRICHLET_ROW_REJECTION_REASON')
                bad_matrix_cases.append({'row': op.LABELS[i], 'column_class': column_class,
                                         'rejected': True})
            else:
                raise ValueError('MISSING_ONE_DIRICHLET_ROW_GUARD:'+op.LABELS[i])
        rejected.append(op.LABELS[i])
    equal(len(rejected), 5, 'FIVE_INDIVIDUAL_DIRICHLET_GUARDS')
    equal(len(bad_matrix_cases), 10, 'TEN_DIRICHLET_ROW_INPUT_CLASS_CONTROLS')
    # Embed the authenticated old mean matrices/profiles in the complete space.
    old = json.loads((roots['current']/'results/PROFILES.json').read_text())
    golden = json.loads((roots['current']/'results/EVALUATED_CANDIDATE.json').read_text())
    slots = [op.IDX[f'{field}:0'] for field in ('K', 'L', 'N', 'J', 'w')]
    def embed_vector(row):
        out = [Q(0)]*DIM
        for i, value in zip(slots, row):
            out[i] = Q(value)
        return out
    def embed_matrix(rows):
        out = [[Q(0)]*DIM for _ in range(DIM)]
        for i, row in zip(slots, rows):
            for j, value in zip(slots, row):
                out[i][j] = Q(value)
        return out
    embedded = {'nodes': old['nodes'],
                'profiles': {name: [embed_vector(row) for row in rows] for name, rows in old['profiles'].items()},
                'matrices': {name: embed_matrix(rows) for name, rows in old['matrices'].items()}}
    retained = k.from_data(embedded)
    for row in golden['samples']:
        actual = retained.feedback(Q(row['x']), Q(row['r']), embedded['matrices']['M'])
        for key in ('value', 'dx', 'dr'):
            equal(actual[key], tuple(embed_vector(row['feedback_exact'][key])), 'OLD_SIX_SAMPLE_MEAN_EMBEDDING')
    equal(len(golden['samples']), 6, 'OLD_MEAN_SAMPLE_COUNT')
    # Exact vector identity behind the two-term complete error ledger. This
    # uses an arbitrary dense T and makes no commutation or diagonal claim.
    T = tuple(tuple(Q((-1)**(i+j)*(i+j+1), 10**6) for j in range(DIM)) for i in range(DIM))
    U = tuple(Q((-1)**i*(i+1), 7) for i in range(DIM))
    error0 = tuple(Q(i-3, 101) for i in range(DIM))
    error1 = tuple(Q(i+2, 307) for i in range(DIM))
    b = add(U, scale(-1, mv(T, U)))
    B0 = add(b, scale(-1, error0))
    B1 = add(mv(T, B0), scale(-1, error1))
    equal(add(U, scale(-1, B0), scale(-1, B1)),
          add(mv(T, mv(T, U)), error0, mv(T, error0), error1), 'COMPLETE_TWO_TERM_ERROR_IDENTITY')
    need({key: auth.authenticate(root, auth.PINS[key]) for key, root in roots.items()} == before,
         'FROZEN_CLOSURES_CHANGED')
    need(sha(kernel_path) == KERNEL_PIN, 'NEW_KERNEL_CHANGED')
    report = {'schema': 'HDBLAST_FULL23_INDEPENDENT_RIDGE_REVIEW_V1', 'status': 'PASS',
              'kernel_sha256': KERNEL_PIN, 'authentication_before_and_after': before,
              'dimension_only_changes_verified': 10, 'oracle': 'Exact independent clipped-polygon Green moments plus direct endpoint integrals; normalized piecewise primitives integrated segment by segment.',
              'counts_by_exact_comparison_family': counts,
              'exact_comparison_groups': sum(counts.values()),
              'polynomial_degree': 4, 'coupled_matrix_dimension': DIM,
              'dirichlet_individual_bad_rows_rejected': rejected,
              'bad_matrix_rows_both_input_classes': bad_matrix_cases,
              'noncommuting_wrong_reflection_witnesses': witnesses,
              'coupled_sample_points': samples,
              'old_mean_embedding_samples': len(golden['samples']),
              'exact_C1_join_interfaces': [row[0] for row in joins],
              'finite_candidate_Robin_identity': 'PN V_r-M PN V=M PN(uP1-u0-B1) at the shell; for coherent uP1=u0 this is -M PN B1.',
              'complete_error_identity': 'U-B0-B1=T^2 U+(I+T)(b-B0)+(T B0-B1).',
              'complete_error_bound': 'If ||T||<=q, ||U||<=Ucap, ||b-B0||<=e0 and ||T B0-B1||<=e1 in the same full23 X norm, then ||U-B0-B1||<=q^2*Ucap+(1+q)*e0+e1.',
              'limits': ['Finite exact kernel/representation controls, not a new actual HDBLAST source or continuous response-error certificate.',
                         'The oracle does not call producer kernels to obtain expected integrals.',
                         'All 23 amplitudes, five D rows, eighteen N rows, and noncommuting coupling/reflection cases are retained.',
                         'No old campaign replay or frozen file edit.']}
    with args.output.open('x') as handle:
        json.dump(encode(report), handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': 'PASS', 'checks': sum(counts.values()), 'old_samples': 6,
                      'bad_D_rows_rejected': len(rejected), 'reflection_witnesses': len(witnesses),
                      'report_sha256': sha(args.output)}, sort_keys=True))


if __name__ == '__main__':
    main()
