#!/usr/bin/env python3
"""Validated causal panel quadrature of the first actual-bulk correction.

Every panel uses a whole-panel source interval and an exact positive weight
mass.  The affine source is enclosed by its genuine component caps.  The
actual Robin history is folded into moments, with a rigorous factorial tail.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
EVALUATOR_PIN = '65879ec37d509f750d70c82a77c40d73ac9d34e0477caef399b71d99fa466e02'
LOCAL_PIN = '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899'
BOUND_PIN = '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f'


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def iv(value):
    need(len(value) == 2, 'interval shape')
    lo, hi = map(Q, value)
    need(lo <= hi, 'interval order')
    return lo, hi


def add(left, right):
    return left[0] + right[0], left[1] + right[1]


def neg(value):
    return -value[1], -value[0]


def sub(left, right):
    return add(left, neg(right))


def mul(left, right):
    values = [x * y for x in left for y in right]
    return min(values), max(values)


def scale(value, scalar):
    return mul(value, (scalar, scalar))


def total(values):
    out = (Q(0), Q(0))
    for value in values:
        out = add(out, value)
    return out


def mv(matrix, vector):
    need(len(matrix) == len(vector) and all(len(row) == len(vector) for row in matrix),
         'matrix/vector dimensions')
    return [total(mul(entry, value) for entry, value in zip(row, vector)) for row in matrix]


def vector_add(left, right):
    return [add(x, y) for x, y in zip(left, right)]


def wide(value, radius):
    need(radius >= 0, 'nonnegative outward radius')
    return value[0] - radius, value[1] + radius


def intersect(left, right):
    out = max(left[0], right[0]), min(left[1], right[1])
    need(out[0] <= out[1], 'independent correction enclosures must overlap')
    return out


def encode(value):
    return list(map(str, value))


def sqrt_upper(value):
    need(value >= 0, 'nonnegative squared norm')
    grid = 1 << 112
    integer = math.isqrt(value.numerator * grid * grid // value.denominator)
    if integer * integer * value.denominator < value.numerator * grid * grid:
        integer += 1
    result = Q(integer, grid)
    need(result * result >= value, 'outward rational square root')
    return result


def load_pinned(path, digest):
    data = path.read_bytes()
    need(hashlib.sha256(data).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(data)


def area_moment(a, j):
    """Integral over s,y>=0, 1<s+y<a of (a-s-y)^j/j!."""
    d = a - 1
    return d ** (j + 1) / math.factorial(j) * (a / (j + 1) - d / (j + 2))


def panel_mass(a, j, alpha, theta):
    lo, hi = alpha
    need(1 <= lo <= hi <= a and -1 <= theta[0] <= theta[1] <= 1,
         'physical active alpha/theta panel')
    # Positive primitive in u=a-alpha; the alpha/2 Jacobian is included.
    def primitive(u):
        return a * u ** (j + 1) / (j + 1) - u ** (j + 2) / (j + 2)
    mass = (theta[1] - theta[0]) / (2 * math.factorial(j)) * (
        primitive(a - lo) - primitive(a - hi))
    need(mass >= 0, 'positive exact panel mass')
    return mass


def partition(lo, hi, panels, extras=()):
    need(panels > 0 and lo < hi, 'positive quadrature partition')
    points = {lo + (hi - lo) * Q(j, panels) for j in range(panels + 1)}
    points.update(x for x in extras if lo < x < hi)
    points = sorted(points)
    need(points[0] == lo and points[-1] == hi and
         all(x < y for x, y in zip(points, points[1:])), 'exact partition cover')
    return list(zip(points, points[1:]))


def geometry_checks():
    a = Q(9, 8)
    aa = partition(Q(1), a, 7, [Q(65, 64), Q(66, 64)])
    tt = partition(Q(-1), Q(1), 11, [Q(-8, 9), Q(8, 9)])
    groups = 0
    for j in range(7):
        mass = sum(panel_mass(a, j, x, t) for x in aa for t in tt)
        need(mass == area_moment(a, j), 'exact weighted causal triangle cover')
        groups += 1
    need(area_moment(a, 0) == Q(17, 128), 'causal area is not a thin d*d triangle')
    need(sum(a * (t[1] - t[0]) / 2 for t in tt) == a,
         'complete physical null-edge measure')
    # Physical-coordinate polynomial fixtures integrate alpha^p and theta^q.
    # These are independent of interval source evaluation.
    fixtures = 0
    for p in range(5):
        for q in range(5):
            alpha_integral = (a ** (p + 2) - 1) / (p + 2)
            theta_integral = Q(0) if q % 2 else Q(2, q + 1)
            transformed = alpha_integral * theta_integral / 2
            physical = Q(0)
            # (s+y)^p ((y-s)/(s+y))^q dsdy =
            # alpha^(p+1)*theta^q/2; rectangular exact polynomial integral.
            for ar in aa:
                for tr in tt:
                    physical += ((ar[1] ** (p + 2) - ar[0] ** (p + 2)) / (p + 2)
                                 * (tr[1] ** (q + 1) - tr[0] ** (q + 1)) / (q + 1) / 2)
            need(physical == transformed, 'exact transformed polynomial fixture')
            fixtures += 1
    d, m = a - 1, Q(1, 1000)
    for n in range(6):
        tail = m ** (n + 1) * area_moment(a, n + 1) / (1 - m * d / (n + 2))
        partial = sum(m ** j * area_moment(a, j) for j in range(n + 1, n + 12))
        need(partial < tail, 'factorial tail exceeds positive finite tail fixtures')
    # Test matrix order: an output reflection cannot be moved through coupling.
    matrix = [[(Q(0), Q(0)), (Q(2), Q(2))],
              [(Q(3), Q(3)), (Q(0), Q(0))]]
    source = [(Q(5), Q(5)), (Q(7), Q(7))]
    parity = [1, -1]
    correct = [scale(z, p) for z, p in zip(mv(matrix, source), parity)]
    wrong = mv(matrix, [scale(z, p) for z, p in zip(source, parity)])
    need(correct != wrong, 'noncommuting source/reflection corruption rejected')
    return {'named_groups': groups + 5, 'polynomial_fixtures': fixtures,
            'tail_fixtures': 6, 'noncommuting_reflection_controls': 1}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path)
    parser.add_argument('--contact-capsule', type=Path)
    parser.add_argument('--bulk-capsule', type=Path)
    parser.add_argument('--evaluator', type=Path, default=Path(__file__).with_name('interior_evaluator.py'))
    parser.add_argument('--alpha-panels', type=int, default=16)
    parser.add_argument('--theta-panels', type=int, default=32)
    parser.add_argument('--primitive-panels', type=int, default=128)
    parser.add_argument('--history-degree', type=int, default=3)
    parser.add_argument('--self-test', action='store_true')
    args = parser.parse_args()
    controls = geometry_checks()
    if args.self_test:
        print(json.dumps({'status': 'PASS', 'geometry_controls': controls}, sort_keys=True))
        return
    need(args.repo_root is not None, 'repo root is required')
    need(args.history_degree >= 0, 'nonnegative history degree')
    need(args.alpha_panels > 0 and args.theta_panels > 0 and args.primitive_panels > 0,
         'positive quadrature panel counts')
    root = args.repo_root.resolve()
    bulk = args.bulk_capsule or root / 'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
    local = load_pinned(bulk / 'reference/LOCAL.json', LOCAL_PIN)
    bound = load_pinned(bulk / 'reference/BULK_BOUND.json', BOUND_PIN)
    need(local['status'] == bound['status'] == 'PASS' and local['point'] == ['9/8', '0'],
         'authenticated prior local correction certificate')
    data = args.evaluator.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    need(EVALUATOR_PIN is not None and digest == EVALUATOR_PIN, 'INPUT_HASH:interior evaluator')
    spec = importlib.util.spec_from_file_location('verified_interior_bulk_evaluator', args.evaluator)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    evaluator = module.InteriorEvaluator(root=root, primitive_panels=args.primitive_panels,
                                         contact_capsule=args.contact_capsule, bulk_capsule=bulk)
    evaluator.authenticate()
    labels, weights, parity = evaluator.labels, list(map(Q, evaluator.weights)), evaluator.parity
    matrix = [[iv(z) for z in row] for row in evaluator.M]
    need(len(labels) == len(weights) == len(parity) == 23 and min(weights) > 0,
         'complete positive full23 modal space')
    need([z['label'] for z in local['full23_point_table']] == labels,
         'same prior modal row order')
    m, a = Q(evaluator.m), Q(9, 8)
    d, degree = a - 1, args.history_degree
    need(0 < m * d < 1 and m == Q(bound['actual_coefficient_H_norm_caps']['M']),
         'same actual constant Robin norm cap')
    for i in range(23):
        for j in range(23):
            if parity[i] < 0 or parity[j] < 0:
                need(matrix[i][j] == (0, 0), 'Robin/Dirichlet closure')
    Rlocal = Q(local['actual_reduced_bulk_residual_local_H_cap'])
    need(Q(evaluator.FH) == Q(local['genuine_pointwise_affine_F_H_cap']),
         'same genuine affine-source units')
    fcap = list(map(Q, evaluator.Fcap))
    need(len(fcap) == 23 and min(fcap) >= 0, 'actual affine source component caps')
    # Splitting alpha at profile joins improves ranges without omitting any
    # transition. The evaluator also handles intervals crossing these joins.
    aa = partition(Q(1), a, args.alpha_panels, [Q(65, 64), Q(66, 64)])
    tt = partition(Q(-1), Q(1), args.theta_panels, [Q(-8, 9), Q(8, 9)])
    moments = [[(Q(0), Q(0)) for _ in range(23)] for _ in range(degree + 1)]
    masses = [Q(0)] * (degree + 1)
    lower_panels = 0
    for alpha in aa:
        for theta in tt:
            values = evaluator.evaluate(alpha, theta)
            residual = list(map(iv, values['r0_lower']))
            need(len(residual) == 23, 'complete lower-order source panel')
            for j in range(degree + 1):
                mass = panel_mass(a, j, alpha, theta)
                moments[j] = vector_add(moments[j], [scale(z, mass) for z in residual])
                masses[j] += mass
            lower_panels += 1
    lower_moments = [[encode(z) for z in row] for row in moments]
    for j in range(degree + 1):
        need(masses[j] == area_moment(a, j), 'whole weighted source triangle covered')
        moments[j] = [wide(z, c * masses[j]) for z, c in zip(moments[j], fcap)]
    # Applying successive powers of the same actual M is enclosed by interval
    # matrix products; independence broadens the range, never freezes M.
    trace = [(Q(0), Q(0))] * 23
    history_terms = []
    for j, moment in enumerate(moments):
        term = moment
        for _ in range(j):
            term = [neg(z) for z in mv(matrix, term)]
        trace = vector_add(trace, term)
        history_terms.append([encode(z) for z in term])
    tail_H = (Rlocal * m ** (degree + 1) * area_moment(a, degree + 1)
              / (1 - m * d / (degree + 2)))
    tail_components = [Q(0) if parity[i] < 0 or all(z == (0, 0) for z in matrix[i])
                       else tail_H * sqrt_upper(1 / weights[i]) for i in range(23)]
    trace = [wide(z, radius) if p > 0 else (Q(0), Q(0))
             for z, radius, p in zip(trace, tail_components, parity)]
    p_edge = [(Q(0), Q(0))] * 23
    edge_mass = Q(0)
    for theta in tt:
        values = evaluator.evaluate((a, a), theta)
        residual = list(map(iv, values['r0_lower']))
        mass = a * (theta[1] - theta[0]) / 2
        p_edge = vector_add(p_edge, [scale(z, mass) for z in residual])
        edge_mass += mass
    need(edge_mass == a, 'complete physical edge covered')
    lower_edge = [encode(z) for z in p_edge]
    p_edge = [wide(z, a * c) for z, c in zip(p_edge, fcap)]
    mtrace = mv(matrix, trace)
    Ctime = [sub(z, bz) if p > 0 else (Q(0), Q(0))
             for z, bz, p in zip(p_edge, mtrace, parity)]
    Cradial = [bz if p > 0 else z for z, bz, p in zip(p_edge, mtrace, parity)]
    analytic_trace = Q(local['bulk_update_correction_local_value_and_trace_H_cap'])
    analytic_time = Q(local['bulk_update_correction_shell_time_derivative_H_cap'])
    table, nonzero = [], {'value': 0, 'dx': 0, 'dr': 0}
    strict = 0
    for i, label in enumerate(labels):
        rowdual = sqrt_upper(sum(max(abs(lo), abs(hi)) ** 2 / weights[j]
                                 for j, (lo, hi) in enumerate(matrix[i])))
        factor = sqrt_upper(1 / weights[i])
        if parity[i] > 0:
            radii = {'value': factor * analytic_trace, 'dx': factor * analytic_time,
                     'dr': rowdual * analytic_trace}
        else:
            radii = {'value': Q(0), 'dx': Q(0), 'dr': factor * a * Rlocal}
        raw = {'value': trace[i], 'dx': Ctime[i], 'dr': Cradial[i]}
        correction = {key: intersect(value, (-radii[key], radii[key]))
                      for key, value in raw.items()}
        for key, value in correction.items():
            nonzero[key] += int(value[0] > 0 or value[1] < 0)
            strict += int(value[1] - value[0] < 2 * radii[key])
        table.append({'label': label, 'parity': parity[i],
                      'raw_quadrature_correction': {key: encode(z) for key, z in raw.items()},
                      'correction': {key: encode(z) for key, z in correction.items()},
                      'prior_analytic_symmetric_correction_radii': {key: str(z) for key, z in radii.items()}})
    evaluator.authenticate()
    need(hashlib.sha256(args.evaluator.read_bytes()).hexdigest() == digest,
         'interior evaluator unchanged during execution')
    result = {'status': 'PASS', 'method': 'exact positive-weight whole-panel interval causal quadrature',
              'point': ['9/8', '0'], 'labels': labels,
              'input_sha256': {'interior_evaluator': digest, 'prior_local': LOCAL_PIN,
                               'prior_global': BOUND_PIN},
              'geometry_controls': controls, 'alpha_panels': len(aa), 'theta_panels': len(tt),
              'source_panels': lower_panels, 'edge_panels': len(tt),
              'primitive_panels_requested': args.primitive_panels, 'history_degree': degree,
              'moment_positive_weight_masses': list(map(str, masses)),
              'lower_order_moments': lower_moments,
              'complete_affine_source_moments': [[encode(z) for z in row] for row in moments],
              'history_terms': history_terms, 'history_tail_H_cap': str(tail_H),
              'history_tail_component_caps': list(map(str, tail_components)),
              'lower_order_physical_edge': lower_edge, 'complete_physical_edge': [encode(z) for z in p_edge],
              'actual_affine_source_component_caps': list(map(str, fcap)),
              'strictly_narrower_correction_boxes': strict, 'nonzero_correction_intervals': nonzero,
              'full23_point_table': table, 'interior_evaluator_metadata': evaluator.metadata(),
              'source_scope': 'same fixed actual affine source enclosed by authenticated uniform component caps; no midpoint or zero source substitution',
              'correction_definition': '(I-L M PN)^-1 K[(La-L0)Bs+La(W0-Bs)+F]',
              'scope': 'conditional reduced-system first actual-bulk correction at one shell point; no self-consistent fixed point, quantitative C2 or Einstein residual, nonlinear spacetime or physical observable'}
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
