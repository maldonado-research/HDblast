#!/usr/bin/env python3
"""Join certified actual initial profiles to the retained second-response inverse.

This evaluates upper bounds, not the second-response waveform. All decisions
use exact rational arithmetic and input hashes supplied separately in PINS.json.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


def authenticate(root):
    pins = json.loads((root / 'PINS.json').read_text())
    data = {}
    for name, pin in pins.items():
        raw = (root / pin['path']).read_bytes()
        need(hashlib.sha256(raw).hexdigest() == pin['sha256'], 'input_hash_' + name)
        need(len(raw) == pin['bytes'], 'input_length_' + name)
        data[name] = json.loads(raw)
    return data, pins


def least_terms(b, q, factor, tolerance):
    n = 0
    tail = factor * b / (1-q)
    while tail > tolerance:
        n += 1
        tail *= q
    need(tail <= tolerance, 'tail_upper')
    need(n == 0 or tail / q > tolerance, 'tail_minimality')
    return {'terms': n, 'tail_upper_exact': str(tail),
            'previous_tail_exceeds_target': n == 0 or tail / q > tolerance}


def run(root):
    data, pins = authenticate(root)
    old, initial, evaluation = (data[k] for k in ('parent', 'initial', 'evaluation'))
    need(old['status'] == 'PASS_EXACT_SECOND_RESPONSE_MAJORANTS', 'parent_status')
    need(initial['status'] == 'PASS_CONTINUOUS_ACTUAL_SECOND_INITIAL_PROFILE_ENCLOSURES', 'initial_status')
    need(evaluation['status'] == 'PASS_UNPUBLISHED_SELECTED23_EVALUATION_BOUND', 'evaluation_status')
    need(initial['domain']['time'] == '0' and initial['domain']['whole_sphere'], 'initial_domain')
    need(initial['domain']['smooth_zero_extension_outside_packet'], 'initial_support')
    need(initial['domain']['ell'] == '1/1000000', 'initial_scale')
    need(initial['scope']['mixed_shift_velocity_retained'] and
         initial['scope']['scalar_initial_and_velocity_retained'] and
         initial['scope']['initial_lapse_square_retained'], 'initial_conventions')
    cap = lambda k: Q(initial['component_caps'][k]['exact_upper'])
    # The positive Frobenius norm of a symmetric 5x5 metric tensor is at
    # most 5 times a bound on each frame component, including mixed entries.
    # Add the scalar norm. These bounds hold over the whole initial annulus.
    u0 = 5*cap('metric_value') + cap('scalar_value')
    ur = 5*cap('metric_scaled_radial_partial') + cap('scalar_scaled_radial_partial')
    ux = 5*cap('metric_scaled_time_partial') + cap('scalar_scaled_time_partial')
    free = max(u0/3 + ux, ur + ux)
    q = Q(old['contraction']['registered_upper'])
    need(q == Q(2109, 1000000), 'retained_contraction')
    force = Q(old['source_join']['F_H_upper'])
    boundary = Q(old['source_join']['G_H_upper'])
    b = free + 3*force + boundary
    response = b/(1-q)
    simple = Q('1.03e33')
    need(response < simple, 'refined_response_simple_bound')
    sharp_squared = Q(evaluation['squared_constant'])
    need(sharp_squared == Q(321, 19), 'selected_evaluation_squared')
    point_factor = Q(4111, 1000)
    need(point_factor**2 > sharp_squared, 'rounded_point_factor')
    need(point_factor < 17, 'improved_selected_point_bound')
    targets = [('1e30', Q(10**30)), ('1e12', Q(10**12)),
               ('1', Q(1)), ('1e-6', Q(1,10**6)), ('1e-12', Q(1,10**12))]
    table = {label: {
        'X_norm': least_terms(b, q, Q(1), tol),
        'selected_pointwise_value': least_terms(b, q, 3*point_factor, tol),
        'selected_scaled_first_derivative': least_terms(b, q, point_factor, tol)
        } for label,tol in targets}
    return {
        'status': 'PASS_REFINED_FORMAL_SECOND_RESPONSE_MAJORANTS',
        'input_sha256': {k:v['sha256'] for k,v in pins.items()},
        'conditional_premises': [
            'Retained smooth background, actual first master/gauge solutions, actual second initial prescription and accepted full causal inverse.',
            'Sources, selected initial data, full linear operator and boundary rows preserve the selected quadratic-feature space E_v; see the separate angular closure proof and review.',
            'Inherited bulk/shell source and operator bounds remain unchanged and are authenticated premises.'
        ],
        'norm': 'max(sup|U|_H/3, sup|U_x|_H, sup|U_r|_H), probability-normalized S3 and positive Frobenius metric norm',
        'initial_H_bounds': dict(zip(['U0','U0_r','U0_x','free_X'], map(str,[u0,ur,ux,free]))),
        'retained_join': {'q':str(q),'F_H':str(force),'G_H':str(boundary)},
        'b_X_upper': str(b),
        'response': {'X_exact_upper':str(response),'X_simple_upper':str(simple),
                     'pointwise_factor_upper':str(point_factor),
                     'pointwise_value_upper':str(3*point_factor*simple),
                     'pointwise_scaled_first_derivative_upper':str(point_factor*simple),
                     'old_registered_X_upper':old['response']['X_simple_upper'],
                     'old_to_new_registered_X_bound_ratio':str(Q(old['response']['X_simple_upper'])/simple)},
        'least_exact_Neumann_tail_terms':table,
        'scope':{
            'improved_upper_bound_not_computed_amplitude':True,
            'selected_solution_sharp_point_factor_squared':str(sharp_squared),
            'full_287_space_point_factor_unchanged':17,
            'actual_second_response_numerically_evaluated':False,
            'background_approximation_error_budget_complete':False,
            'quadratic_source_approximation_errors_complete':False,
            'continuous_second_response_residual_certified':False,
            'tail_counts_exclude_input_quadrature_interpolation_and_rounding_errors':True,
            'exact_nonlinear_family_or_invariant_force':False
        },
        'arithmetic':'Exact fractions.Fraction acceptance; no floating-point comparisons.'
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse_existing_output')
    report = run(Path(__file__).resolve().parent)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
