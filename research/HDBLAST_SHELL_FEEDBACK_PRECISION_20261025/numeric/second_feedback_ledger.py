#!/usr/bin/env python3
"""Complete conditional error of B0+S B0+S^2 B0, with exact nominal u0.

S is the full frozen constant-matrix wave/lift operator. Frozen P1 errors
bound norms of nominal u0 and its derivative only; they are not defects in
this candidate. The frozen P1 velocity still defines B0 exactly.
"""
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('component_ledger_core', HERE/'component_ledger.py')
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)
add, matvec, matmul = core.add, core.matvec, core.matmul


def zeros(rows, columns):
    return [[F(0)]*columns for _ in range(rows)]


def matrix_add(*matrices):
    return [[sum((matrix[i][j] for matrix in matrices), F(0))
             for j in range(len(matrices[0][0]))] for i in range(len(matrices[0]))]


def matrix_scale(matrix, coefficient):
    return [[coefficient*entry for entry in row] for row in matrix]


def polynomial_add(*polynomials):
    result = []
    rows, columns = len(polynomials[0][0]), len(polynomials[0][0][0])
    for k in range(max(map(len, polynomials))):
        result.append(matrix_add(*(p[k] if k < len(p) else zeros(rows, columns) for p in polynomials)))
    return result


def polynomial_multiply(left, right):
    rows, columns = len(left[0]), len(right[0][0])
    result = [zeros(rows, columns) for _ in range(len(left)+len(right)-1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i+j] = matrix_add(result[i+j], matmul(a, b))
    return result


def integral(polynomial, left, right):
    result = zeros(len(polynomial[0]), len(polynomial[0][0]))
    for k, matrix in enumerate(polynomial):
        result = matrix_add(result, matrix_scale(matrix, (right**(k+1)-left**(k+1))/(k+1)))
    core.need(all(value >= 0 for row in result for value in row), 'NONNEGATIVE_INTEGRATED_KERNEL')
    return result


def evaluate(polynomial, argument):
    result = zeros(len(polynomial[0]), len(polynomial[0][0]))
    for k, matrix in enumerate(polynomial):
        result = matrix_add(result, matrix_scale(matrix, argument**k))
    core.need(all(value >= 0 for row in result for value in row), 'NONNEGATIVE_COMPARISON_AT_ARGUMENT')
    return result


def blocks(matrices):
    C = matrix_add(matrices['A'], matrices['B'])
    D, M = matrices['D'], matrices['M']
    z = zeros(5, 5)
    H0 = [matrix_add(matrix_scale(D, -F(1, 2)), matrix_scale(M, -1))[i]
          +matrix_scale(C, -F(1, 2))[i] for i in range(5)]
    H0 += [M[i][:]+z[i][:] for i in range(5)]
    H1 = [M[i][:]+z[i][:] for i in range(5)]+[D[i][:]+C[i][:] for i in range(5)]
    H2 = [matrix_scale(D, F(1, 2))[i]+matrix_scale(C, F(1, 2))[i] for i in range(5)]
    H2 += [z[i][:]+z[i][:] for i in range(5)]
    Lambda = [[M[i][:]+z[i][:] for i in range(5)],
              [D[i][:]+C[i][:] for i in range(5)]]
    outer = [D[i][:]+C[i][:] for i in range(5)]
    return {'H': [H0, H1, H2], 'Lambda': Lambda, 'outer': outer, 'M': M}


def second_shell_kernel(outer, inner, time):
    # K_outer(inner z) uses the full interior ten-state comparison.
    # L_outer(inner z) uses the nested PN shell integral; no L0 interior
    # row/column is removed from either matrix product.
    bulk = polynomial_multiply([zeros(5, 10), outer['outer']], inner['H'])
    shell = polynomial_multiply([matrix_scale(outer['M'], time), matrix_scale(outer['M'], -1)],
                                inner['Lambda'])
    return polynomial_add(bulk, shell)


def upward(values):
    grid = 1 << 112
    return [F((value*grid).__ceil__(), grid) for value in values]


def run(context):
    roots, frozen, profiles = context['roots'], context['frozen'], context['profiles']
    selected = context['selected']
    scope, ten_from = context['scope'], context['ten_from']
    before, dependency_before = context['before'], context['dependency_before']
    matrices = context['matrices']
    nominal_abs = {name: [[abs(value) for value in row] for row in matrix]
                   for name, matrix in matrices.items()}
    star = blocks(nominal_abs)
    cache = {row['time_limit']: row for row in context['ledgers'][1:]}
    ten_cache = {}

    def get(time):
        if time not in cache:
            cache[time] = scope('second_feedback_interval_upper', time, time, time-1)
        ledger = cache[time]
        if time not in ten_cache:
            ten_cache[time] = ten_from(ledger)
        ten = ten_cache[time]
        entry = ledger['matrix_entries']
        actual = blocks({name: entry[name+'_abs'] for name in ('A', 'B', 'D', 'M')})
        delta = blocks(entry['coefficient_entry_errors'])
        H = evaluate(actual['H'], time)
        D = evaluate(delta['H'], time)
        Hstar = evaluate(star['H'], time)
        R, contraction = core.positive_resolvent(H)
        beta, e0 = ten['beta_value_derivative'], ten['e0_value_derivative']
        Rbeta, Re0 = upward(matvec(R, beta)), upward(matvec(R, e0))
        return ledger, ten, actual, delta, H, D, Hstar, Rbeta, Re0, contraction

    def velocity_integral(time):
        result = [F(0)]*5
        for rows, initial in ((frozen['rows'], True), (selected['cell_certificates'], False)):
            for row in rows:
                left, right = map(F, row['r'])
                width = max(F(0), min(time, right)-left)
                for i in range(5):
                    cap = (sum((F(row['channel_component_error_caps'][channel][2][i])
                                for channel in range(3)), F(0)) if initial
                           else F(row['component_error_caps'][2][i]))
                    result[i] += width*cap
        return result

    results = []
    for time in (F(9, 8), F(5, 4)):
        ledger, ten, actual, delta, H, D, Hstar, Rbeta, Re0, contraction = get(time)
        H2, H3 = matmul(H, H), matmul(matmul(H, H), H)
        coefficient_left = matrix_add(D, matmul(D, H), matmul(Hstar, D))
        coefficient_right = matrix_add(D, matmul(H, D), matmul(D, Hstar))
        generic_tail = upward(matvec(H3, Rbeta))
        generic_coeff_left = upward(matvec(coefficient_left, ten['beta_value_derivative']))
        generic_coeff_right = upward(matvec(coefficient_right, ten['beta_value_derivative']))
        generic_error = upward(add(Re0, generic_tail,
                                  [min(a, b) for a, b in zip(generic_coeff_left, generic_coeff_right)]))
        boundaries = [F(1)]+[F(r) for r in profiles['nodes'] if F(1) < F(r) <= time]
        tail, input_route_a, input_first, input_higher = ([F(0)]*5 for _ in range(4))
        coefficient_a, coefficient_b = ([F(0)]*5 for _ in range(2))
        coefficient_first, coefficient_delta_T, coefficient_star_delta = ([F(0)]*5 for _ in range(3))
        intervals = []
        for left, right in zip(boundaries, boundaries[1:]):
            local, local_ten, T, E, h, d, hs, rb, re, norm = get(right)
            beta, e0 = local_ten['beta_value_derivative'], local_ten['e0_value_derivative']
            Wtt = second_shell_kernel(T, T, time)
            third_kernel = polynomial_multiply(Wtt, T['H'])
            tail_part = matvec(integral(third_kernel, left, right), rb)
            tail = upward(add(tail, tail_part))
            input_route_a = upward(add(input_route_a,
                                       matvec(integral(T['Lambda'], left, right), re)))
            initial_cell_u = [max(sum((F(row['channel_component_error_caps'][channel][0][i])
                                      for channel in range(3)), F(0))
                                  for row in frozen['rows'] if F(row['r'][0]) < right and F(row['r'][1]) > left)
                              for i in range(5)]
            shell_input = add([(right-left)*(a+b) for a, b in zip(initial_cell_u, velocity_integral(right))],
                              [(right-left)*(right+left-2)/2*core.sqrt_upper(core.SOURCE**2/F(weight))
                               for weight in core.WEIGHTS])
            input_first = upward(add(input_first,
                [(right*right-left*left)/2*x for x in matvec(T['outer'], e0)],
                matvec(T['M'], shell_input)))
            input_higher = upward(add(input_higher,
                                     matvec(integral(Wtt, left, right), re)))
            first_part = matvec(integral(E['Lambda'], left, right), beta)
            delta_T_part = matvec(integral(second_shell_kernel(E, T, time), left, right), beta)
            star_delta_part = matvec(integral(second_shell_kernel(star, E, time), left, right), beta)
            coefficient_first = upward(add(coefficient_first, first_part))
            coefficient_delta_T = upward(add(coefficient_delta_T, delta_T_part))
            coefficient_star_delta = upward(add(coefficient_star_delta, star_delta_part))
            coefficient_a = upward(add(coefficient_a, first_part, delta_T_part, star_delta_part))
            coefficient_b = upward(add(coefficient_b, first_part,
                matvec(integral(second_shell_kernel(T, E, time), left, right), beta),
                matvec(integral(second_shell_kernel(E, star, time), left, right), beta)))
            intervals.append({'a': [left, right], 'comparison_max_row_sum': norm,
                              'third_tail': upward(tail_part)})
        initial_at_time = [min(sum((F(row['channel_component_error_caps'][channel][0][i])
                                   for channel in range(3)), F(0))
                               for row in frozen['rows'] if F(row['r'][0]) <= time <= F(row['r'][1]))
                           for i in range(5)]
        direct_input = add(initial_at_time, velocity_integral(time),
                           [(time-1)*core.sqrt_upper(core.SOURCE**2/F(weight)) for weight in core.WEIGHTS])
        input_a = upward(add(direct_input, input_route_a))
        input_b = upward(add(direct_input, input_first, input_higher))
        chosen_input = [min(a, b) for a, b in zip(input_a, input_b)]
        chosen_coefficient = [min(a, b) for a, b in zip(coefficient_a, coefficient_b)]
        error = upward(add(chosen_input, tail, chosen_coefficient))
        results.append({'time': time, 'integration_cells': len(boundaries)-1,
                        'state_order': list(core.MEANS)+list(core.MEANS),
                        'generic_comparison_max_row_sum': contraction,
                        'generic_input_R_e0': Re0, 'generic_third_tail_H3_R_beta': generic_tail,
                        'generic_coefficient_D_DH_HstarD': generic_coeff_left,
                        'generic_coefficient_D_HD_DHstar': generic_coeff_right,
                        'generic_value_derivative_error': generic_error,
                        'generic_error_outward': list(map(core.outward, generic_error)),
                        'direct_input_shell': direct_input,
                        'input_route_direct_Lambda_R_e0': input_a,
                        'input_route_first_plus_nested_resolvent': input_b,
                        'chosen_input_shell': chosen_input,
                        'integrated_full_third_tail': tail,
                        'coefficient_first_Delta_B0': coefficient_first,
                        'coefficient_Delta_T_B0': coefficient_delta_T,
                        'coefficient_S_Delta_B0': coefficient_star_delta,
                        'coefficient_route_Delta_DeltaT_SDelta': coefficient_a,
                        'coefficient_route_Delta_TDelta_DeltaS': coefficient_b,
                        'chosen_coefficient_shell': chosen_coefficient,
                        'complete_shell_value_error_upper': error,
                        'complete_shell_value_error_outward': list(map(core.outward, error)),
                        'interval_tail_certificates': intervals})
    core.need(core.authenticate(roots) == before, 'FINAL_CLOSURE_CHANGED')
    core.need(context['direct'].authenticate(context['depmap']) == dependency_before,
              'FINAL_UPSTREAM_AUTHENTICATION_CHANGED')
    result = {'schema': 'HDBLAST_EXACT_U_SECOND_FEEDBACK_COMPONENT_LEDGER_V1',
              'candidate': 'B0+S B0+S^2 B0; exact registered nominal u0/du0, frozen vP1; full constant matrices S=Tstar',
              'full_error_identity': 'R_T(B-B0)+T^3 R_T B0+[Delta+Delta T+S Delta]B0',
              'generic_positive_bound': 'R e0+H^3 R beta+[D+D H+Hstar D]beta',
              'coefficient_alternative': 'D+H D+D Hstar',
              'no_old_u_du_p1_defect': True,
              'p1_envelope_scope': 'Frozen u/du P1 envelopes bound exact nominal profile norms beta only.',
              'input_scope': 'Actual-minus-nominal displacement/velocity, frozen vP1 reconstruction and omitted forcing charged once in resolvent input.',
              'authentication': before,
              'complete_closures_authenticated_before_and_after': True,
              'arithmetic': 'Exact fractions; positive inverse identities exact; all accumulation roundings upward to 2^-112.',
              'results': results,
              'limits': ['Conditional retained five-mean response only.',
                         'Candidate evaluation quadrature/arithmetic errors must be added separately.',
                         'No sign claim without joining candidate intervals.',
                         'No frozen source edits, profile campaign rerun, new mesh or branch gate.']}
    output = context['args'].output
    with output.open('x') as handle:
        json.dump(core.encode(result), handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'output': str(output), 'sha256': core.sha(output),
                      'complete_shell_value_error_outward': {str(row['time']): row['complete_shell_value_error_outward']
                                                             for row in results}}, sort_keys=True))


if __name__ == '__main__':
    core.main(context_handler=run)
