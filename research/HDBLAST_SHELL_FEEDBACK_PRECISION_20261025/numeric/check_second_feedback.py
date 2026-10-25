#!/usr/bin/env python3
"""Independent rational controls and interval join for the B2 component ledger."""
import argparse
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('reviewed_second_feedback', HERE/'second_feedback_ledger.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
core = module.core


def product(a, b):
    # Independent rectangular row/column multiplication.
    columns = list(zip(*b))
    return [[sum((x*y for x, y in zip(row, column)), F(0)) for column in columns] for row in a]


def plus(*matrices):
    return [[sum(values, F(0)) for values in zip(*rows)] for rows in zip(*matrices)]


def scale(matrix, value):
    return [[value*x for x in row] for row in matrix]


def zero(rows, columns):
    return [[F(0)]*columns for _ in range(rows)]


def h_at(m, a):
    c, k = a-1, (a*a-1)/2
    C = plus(m['A'], m['B'])
    upper = [plus(scale(m['D'], k), scale(m['M'], c))[i]+scale(C, k)[i] for i in range(5)]
    lower = [plus(scale(m['D'], a), m['M'])[i]+scale(C, a)[i] for i in range(5)]
    return upper+lower


def lambda_at(m, a):
    C = plus(m['A'], m['B'])
    return [plus(scale(m['D'], a), m['M'])[i]+scale(C, a)[i] for i in range(5)]


def second_at(outer, inner, a, time, wrong_shell_weight=False):
    C = plus(outer['A'], outer['B'])
    bulk = product([outer['D'][i]+C[i] for i in range(5)], h_at(inner, a))
    weight = a-1 if wrong_shell_weight else time-a
    shell = product(outer['M'], lambda_at(inner, a))
    return plus(scale(bulk, a), scale(shell, weight))


def boole(function, left, right):
    # Closed five-node Boole rule is exact for every polynomial of degree<=5.
    values = [function(left+(right-left)*F(i, 4)) for i in range(5)]
    weighted = plus(*(scale(value, F(weight)) for value, weight in zip(values, (7, 32, 12, 32, 7))))
    return scale(weighted, (right-left)/90)


def fixture(seed):
    result = {}
    for kind, name in enumerate(('A', 'B', 'D', 'M')):
        matrix = [[F(1+((seed+13*i+7*j+5*kind+3*i*j) % 19), 10000*(kind+1))
                   for j in range(5)] for i in range(5)]
        if name == 'M':
            for i in range(5):
                matrix[1][i] = matrix[i][1] = F(0)
        result[name] = matrix
    return result


def outward_endpoint(value, upper):
    if value == 0:
        return '0'
    magnitude = abs(value)
    exponent = len(str(magnitude.numerator))-len(str(magnitude.denominator))
    while magnitude < F(10)**exponent:
        exponent -= 1
    while magnitude >= F(10)**(exponent+1):
        exponent += 1
    unit = F(10)**(exponent-8)
    count = (value/unit).__ceil__() if upper else (value/unit).__floor__()
    return str(count)+'e'+str(exponent-8)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    roots = {key: args.repo_root.resolve()/'research'/name for key, name in core.LABELS.items()}
    before = core.authenticate(roots)
    checks = []
    T, E, S = fixture(3), fixture(11), fixture(17)
    sets = {'T': T, 'Delta': E, 'S': S}
    for name, matrices in sets.items():
        core.need(any(matrices['A'][1][j] or matrices['A'][j][1] or matrices['B'][1][j] or matrices['B'][j][1]
                      for j in range(5) if j != 1), 'INTERIOR_L0_COUPLING_FIXTURE')
        checks.append({'control': 'retained_L0_interior_couplings_'+name, 'passed': True})
    for time, left in ((F(9, 8), F(1)), (F(9, 8), F(17, 16)), (F(5, 4), F(1))):
        for outer_name, inner_name in (('T', 'T'), ('Delta', 'T'), ('S', 'Delta'), ('T', 'Delta'), ('Delta', 'S')):
            outer, inner = sets[outer_name], sets[inner_name]
            polynomial = module.second_shell_kernel(module.blocks(outer), module.blocks(inner), time)
            expected = boole(lambda a: second_at(outer, inner, a, time), left, time)
            measured = module.integral(polynomial, left, time)
            core.need(measured == expected, 'INDEPENDENT_BOOLE_ORDERED_SECOND_KERNEL')
            checks.append({'control': 'ordered_second_kernel_Boole', 'outer': outer_name,
                           'inner': inner_name, 'time': str(time), 'left': str(left), 'passed': True})
        third = module.polynomial_multiply(module.second_shell_kernel(module.blocks(T), module.blocks(T), time),
                                           module.blocks(T)['H'])
        expected_third = boole(lambda a: product(second_at(T, T, a, time), h_at(T, a)), left, time)
        measured_third = module.integral(third, left, time)
        core.need(measured_third == expected_third, 'INDEPENDENT_BOOLE_FULL_THIRD_KERNEL')
        checks.append({'control': 'full_third_kernel_Boole', 'time': str(time), 'left': str(left), 'passed': True})
        wrong = boole(lambda a: product(second_at(T, T, a, time, True), h_at(T, a)), left, time)
        core.need(wrong != expected_third, 'WRONG_SHELL_WEIGHT_NOT_DETECTED')
        checks.append({'control': 'wrong_a_minus_one_shell_weight_rejected', 'time': str(time), 'left': str(left), 'passed': True})
        quadratic = boole(lambda a: second_at(T, T, a, time), left, time)
        core.need(quadratic != expected_third, 'MISSING_THIRD_OPERATOR_NOT_DETECTED')
        checks.append({'control': 'missing_third_operator_rejected', 'time': str(time), 'left': str(left), 'passed': True})
    ordered = boole(lambda a: second_at(E, T, a, F(9, 8)), F(1), F(9, 8))
    reversed_order = boole(lambda a: second_at(T, E, a, F(9, 8)), F(1), F(9, 8))
    core.need(ordered != reversed_order, 'NONCOMMUTING_ORDER_FIXTURE_DEGENERATE')
    checks.append({'control': 'noncommuting_outer_inner_swap_rejected', 'passed': True})
    # Direct finite-matrix identity, independent of comparison block construction.
    t = [[F(2, 100), F(3, 100)], [F(5, 100), F(7, 100)]]
    s = [[F(1, 100), F(4, 100)], [F(2, 100), F(3, 100)]]
    d = plus(t, scale(s, -1))
    target = plus(product(t, t), scale(product(s, s), -1))
    core.need(target == plus(product(d, t), product(s, d)), 'ORDERED_DELTA_IDENTITY_A')
    core.need(target == plus(product(t, d), product(d, s)), 'ORDERED_DELTA_IDENTITY_B')
    core.need(target != plus(product(t, d), product(s, d)), 'WRONG_NONCOMMUTING_IDENTITY_NOT_DETECTED')
    checks.append({'control': 'both_exact_noncommuting_identities_and_wrong_order', 'passed': True})
    H = h_at(T, F(9, 8))
    R, norm = core.positive_resolvent(H)
    beta = list(map(F, range(1, 11)))
    rb = [sum((a*b for a, b in zip(row, beta)), F(0)) for row in R]
    hb = core.matvec(H, beta)
    h2b = core.matvec(product(H, H), beta)
    exact_tail = [a-b-c-d for a, b, c, d in zip(rb, beta, hb, h2b)]
    third_R = core.matvec(product(product(H, H), H), rb)
    third_only = core.matvec(product(product(H, H), H), beta)
    core.need(exact_tail == third_R, 'FULL_THIRD_RESOLVENT_REMAINDER')
    core.need(all(a > b for a, b in zip(third_R, third_only)), 'MISSING_RESOLVENT_NOT_DETECTED')
    checks.append({'control': 'full_third_resolvent_remainder_and_missing_R_rejected',
                   'comparison_max_row_sum': str(norm), 'passed': True})

    numeric_normal = HERE/'SECOND_FEEDBACK_NORMAL.json'
    numeric_optimized = HERE/'SECOND_FEEDBACK_OPTIMIZED.json'
    balanced = HERE.parent/'balanced/b2'
    balanced_normal = balanced/'EXACT_B2_NORMAL.json'
    balanced_optimized = balanced/'EXACT_B2_OPTIMIZED.json'
    core.need(core.sha(numeric_normal) == core.sha(numeric_optimized)
              == '23883937b7e5a55f1155e8b5c795644253ed5cf0abd51e0c40a09ece50220974', 'NUMERIC_REPORT_PIN_AND_MODES')
    core.need(core.sha(balanced_normal) == core.sha(balanced_optimized)
              == '3db587c30dca20f01e4964ac87f8f5ba457d55d15f7720d09eb840ed3e0e4935', 'CANDIDATE_REPORT_PIN_AND_MODES')
    numeric = json.loads(numeric_normal.read_text())
    candidate = json.loads(balanced_normal.read_text())
    core.need(numeric['complete_closures_authenticated_before_and_after']
              and candidate['complete_closures_authenticated_before_and_after'], 'AUTHENTICATION_FLAGS')
    core.need(numeric['no_old_u_du_p1_defect'], 'CANDIDATE_DEFINITION')
    joined = []
    for row in numeric['results']:
        time = F(row['time'])
        sample = next(sample for sample in candidate['samples'] if F(sample['x']) == time and F(sample['r']) == 0)
        bounds = list(map(F, row['complete_shell_value_error_upper']))
        generic_bounds = list(map(F, row['generic_value_derivative_error'][:5]))
        intervals = []
        generic_intervals = []
        for i, pair in enumerate(sample['evaluation']['three_term']):
            candidate_lo, candidate_hi = map(F, pair)
            lo, hi = candidate_lo-bounds[i], candidate_hi+bounds[i]
            intervals.append({'mean': core.MEANS[i], 'conditional_exact_response_interval': [lo, hi],
                              'outward_interval': [outward_endpoint(lo, False), outward_endpoint(hi, True)],
                              'strictly_negative': hi < 0, 'strictly_positive': lo > 0})
            lo, hi = candidate_lo-generic_bounds[i], candidate_hi+generic_bounds[i]
            generic_intervals.append({'mean': core.MEANS[i], 'conditional_exact_response_interval': [lo, hi],
                                      'outward_interval': [outward_endpoint(lo, False), outward_endpoint(hi, True)],
                                      'strictly_negative': hi < 0, 'strictly_positive': lo > 0})
        core.need(intervals[3]['strictly_negative'], 'JOINED_J_INTERVAL_INCLUDES_ZERO')
        core.need(generic_intervals[3]['strictly_negative'], 'GENERIC_JOINED_J_INTERVAL_INCLUDES_ZERO')
        checks.append({'control': 'joined_conditional_J_interval_strictly_negative', 'time': str(time), 'passed': True})
        checks.append({'control': 'simpler_generic_comparison_J_interval_strictly_negative', 'time': str(time), 'passed': True})
        joined.append({'x': time, 'r': F(0), 'components': intervals,
                       'simpler_generic_comparison_components': generic_intervals,
                       'candidate_interval_already_contains_quadrature_and_arithmetic': True,
                       'actual_response_error_added_once': True})
    core.need(core.authenticate(roots) == before, 'CLOSURE_CHANGED_AFTER_CONTROLS')
    pins = {path.name: core.sha(path) for path in
            (HERE/'component_ledger.py', HERE/'second_feedback_ledger.py', Path(__file__),
             numeric_normal, numeric_optimized, balanced_normal, balanced_optimized,
             balanced/'B2_INDEPENDENT_REVIEW.json')}
    result = {'schema': 'HDBLAST_SECOND_FEEDBACK_INDEPENDENT_KERNEL_AND_JOIN_REVIEW_V1',
              'status': 'PASS_EXACT_NONCOMMUTING_CONTROLS_AND_CONDITIONAL_J_SAMPLE_JOIN',
              'source_and_report_sha256': pins, 'authentication': before,
              'complete_closures_authenticated_before_and_after': True,
              'controls': checks, 'controls_count': len(checks), 'joined_samples': joined,
              'independent_analytic_review': {
                  'reviewer': '/root/coupled_response_astra', 'status': 'APPROVED',
                  'scope': 'Complete noncommutative B2 identity; Lambda/W third tail; both coefficient and input routes.',
                  'premises_checked': ['Retained shell M and nominal Mstar are time-independent.',
                                      'Full interior L0 couplings remain in all matrix products.',
                                      'All local envelopes are frozen at the cell right endpoint.',
                                      'Integration begins at actual zero-wedge boundary a=1.'],
                  'software_controls_are_separate_from_analytic_review': True},
              'limits': ['Two retained-shell samples only; no whole-history sign.',
                         'No intrinsic-observable conclusion without separate geometric tube proof.',
                         'Internal analytic/software review is not external peer review.']}
    with args.output.open('x') as handle:
        json.dump(core.encode(result), handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': core.sha(args.output),
                      'controls_count': len(checks), 'J_intervals_outward': {str(s['x']): s['components'][3]['outward_interval']
                                                                          for s in joined}}, sort_keys=True))


if __name__ == '__main__':
    main()
