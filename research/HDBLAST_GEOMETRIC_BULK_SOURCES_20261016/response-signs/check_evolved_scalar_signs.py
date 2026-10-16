#!/usr/bin/env python3
"""Exact quantitative persistence of two accepted initial scalar-mean signs.

This joins authenticated premises. It does not recompute the parent initial
profiles, solve the evolution, or claim a finite-amplitude nonlinear result.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


MUTATIONS = {
    'wrong_mean_norm': 'mean_projection_Cauchy_Schwarz',
    'wrong_derivative_weight': 'unit_derivative_weight_in_X',
    'omit_evolution_drift': 'continuous_time_integral_enclosure',
    'wrong_time_scale': 'physical_time_equals_ell_times_x',
    'invent_initial_margin': 'initial_mean_interval_authenticated',
    'expand_radius_cell': 'entire_radius_cell_in_initial_witness',
    'lose_sign_room': 'strict_evolved_sign_room',
}


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def authenticate(here):
    pins = json.loads((here / 'INPUT_PINS.json').read_text())
    data = {}
    for name, pin in pins.items():
        path = Path(pin['path'])
        need(not path.is_absolute() and '..' not in path.parts, 'portable_input_' + name)
        raw = (here / path).read_bytes()
        need(len(raw) == pin['bytes'], 'input_length_' + name)
        need(hashlib.sha256(raw).hexdigest() == pin['sha256'], 'input_hash_' + name)
        if path.suffix == '.json':
            data[name] = json.loads(raw)
    return data, pins


def run(here, mutation=None):
    data, pins = authenticate(here)
    initial, witness, refined, parent = (data[k] for k in
                                         ('initial', 'witnesses', 'refined', 'parent'))
    need(initial['status'] == 'PASS_CONTINUOUS_ACTUAL_SECOND_INITIAL_PROFILE_ENCLOSURES',
         'accepted_initial_status')
    need(witness['status'] == 'PASS_ACTUAL_SECOND_INITIAL_SCALAR_MEAN_SIGNED_WITNESSES',
         'accepted_witness_status')
    need(refined['status'] == 'PASS_REFINED_FORMAL_SECOND_RESPONSE_MAJORANTS',
         'accepted_refined_status')
    need(parent['status'] == 'PASS_EXACT_SECOND_RESPONSE_MAJORANTS', 'accepted_parent_status')
    for report in (witness, refined):
        digest = report.get('initial_sha256', report.get('input_sha256', {}).get('initial'))
        need(digest == pins['initial']['sha256'], 'same_actual_initial_data')
    need(refined['input_sha256']['parent'] == pins['parent']['sha256'], 'same_causal_inverse')
    need(initial['conventions']['scalar_mean'] ==
         'Probability-normalized sphere mean, W0=1; zeta0', 'initial_mean_normalization')
    need(initial['domain']['time'] == '0' and initial['domain']['whole_sphere'],
         'initial_whole_sphere_at_time_zero')
    need(initial['domain']['smooth_zero_extension_outside_packet'], 'initial_zero_extension')
    need(initial['scope']['initial_lapse_square_retained'] and
         initial['scope']['mixed_shift_velocity_retained'] and
         initial['scope']['scalar_initial_and_velocity_retained'],
         'same_complete_initial_prescription')

    # Reconstruct only the norm join, not the expensive parent profile proof.
    cap = lambda name: Q(initial['component_caps'][name]['exact_upper'])
    u0 = 5 * cap('metric_value') + cap('scalar_value')
    ur = 5 * cap('metric_scaled_radial_partial') + cap('scalar_scaled_radial_partial')
    ux = 5 * cap('metric_scaled_time_partial') + cap('scalar_scaled_time_partial')
    free = max(u0 / 3 + ux, ur + ux)
    for name, value in (('U0', u0), ('U0_r', ur), ('U0_x', ux), ('free_X', free)):
        need(value == Q(refined['initial_H_bounds'][name]), 'refined_norm_join_' + name)
    q = Q(parent['contraction']['registered_upper'])
    need(q == Q(2109, 10**6) == Q(refined['retained_join']['q']), 'retained_full_contraction')
    bulk = Q(parent['source_join']['F_H_upper'])
    shell = Q(parent['source_join']['G_H_upper'])
    need(bulk == Q(refined['retained_join']['F_H']) and
         shell == Q(refined['retained_join']['G_H']), 'retained_bulk_and_shell_join')
    b = free + 3 * bulk + shell
    exact_response = b / (1 - q)
    simple_response = Q('1.03e33')
    need(b == Q(refined['b_X_upper']) and
         exact_response == Q(refined['response']['X_exact_upper']) and
         simple_response == Q(refined['response']['X_simple_upper']) and
         exact_response < simple_response, 'refined_actual_response_upper')
    need(refined['norm'] ==
         'max(sup|U|_H/3, sup|U_x|_H, sup|U_r|_H), probability-normalized S3 and positive Frobenius metric norm',
         'accepted_X_norm_definition')

    # The scalar mean is pairing with (0,W0); ||(0,W0)||_H^2=integral 1 dmu=1.
    # Its norm is exactly one, witnessed by the constant scalar itself.
    mean_norm_squared = Q(1)
    mean_factor = Q(1, 2) if mutation == 'wrong_mean_norm' else Q(1)
    need(mean_factor >= 0 and mean_factor**2 >= mean_norm_squared,
         'mean_projection_Cauchy_Schwarz')
    # The derivative entries of X have unit weight; the value entry has weight 1/3.
    derivative_factor = Q(1, 3) if mutation == 'wrong_derivative_weight' else Q(1)
    need(derivative_factor >= 1, 'unit_derivative_weight_in_X')
    mean_dx_upper = mean_factor * derivative_factor * simple_response

    ell = Q(initial['domain']['ell'])
    need(ell == Q(parent['domain']['ell']) == Q(1, 10**6), 'same_scaled_coordinates')
    x_interval = [Q(1, 1000), Q(1, 100) if mutation == 'lose_sign_room' else Q(1, 200)]
    need(0 < x_interval[0] < x_interval[1] <= 3, 'positive_closed_time_interval')
    t_interval = [ell * x for x in x_interval]
    if mutation == 'wrong_time_scale':
        t_interval = x_interval[:]
    need(all(t == ell * x for t, x in zip(t_interval, x_interval)),
         'physical_time_equals_ell_times_x')
    # Fundamental theorem of calculus in the accepted C1(D;H) solution class.
    drift_required = x_interval[1] * mean_dx_upper
    drift = Q(0) if mutation == 'omit_evolution_drift' else drift_required
    need(drift >= drift_required, 'continuous_time_integral_enclosure')

    rows = initial['profile_cells']
    cells = []
    expected = [('positive', ['73/64', '585/512'], 1),
                ('negative', ['951/512', '119/64'], -1)]
    for name, ends, sign in expected:
        a, z = map(Q, ends)
        selected = [row for row in rows if list(map(Q, row['r'])) == [a, z]]
        need(len(selected) == 1, 'unique_initial_witness_' + name)
        row = selected[0]
        published = [v for v in witness['witnesses'] if v['sign'] == name]
        need(len(published) == 1 and published[0]['r'] == ends and
             published[0]['certified_scalar_mean_interval'] == row['scalar_mean'],
             'published_witness_matches_actual_initial_cell_' + name)
        lo, hi = map(Q, row['scalar_mean'])
        need(lo <= hi and (lo > 7 * 10**30 if sign > 0 else hi < -7 * 10**30),
             'strict_initial_witness_sign_' + name)
        candidate_initial = [lo, hi]
        if mutation == 'invent_initial_margin' and sign > 0:
            candidate_initial[0] = Q('7.1e30')
        need(candidate_initial[0] <= lo and candidate_initial[1] >= hi,
             'initial_mean_interval_authenticated')
        r_interval = [a, z]
        if mutation == 'expand_radius_cell':
            r_interval[0] -= Q(1, 512)
        need(a <= r_interval[0] < r_interval[1] <= z,
             'entire_radius_cell_in_initial_witness')
        need(r_interval[0] >= 0 and x_interval[1] + r_interval[1] <= 8,
             'whole_time_radius_rectangle_in_causal_domain')
        evolved = [candidate_initial[0] - drift, candidate_initial[1] + drift]
        signed_room = evolved[0] if sign > 0 else -evolved[1]
        need(signed_room > 0, 'strict_evolved_sign_room')
        common_room = 7 * 10**30 - drift
        need(common_room == Q('1.85e30') and signed_room > common_room,
             'rounded_common_strict_sign_margin')
        cells.append({
            'sign': name,
            'scaled_r_closed_interval': list(map(str, r_interval)),
            'initial_scalar_mean_interval': list(map(str, candidate_initial)),
            'evolved_scalar_mean_enclosure': list(map(str, evolved)),
            'signed_magnitude_strict_lower_exact': str(signed_room),
            'signed_magnitude_exceeds_common': str(common_room),
        })

    # A zero exists at each time by continuity in r. No zero branch is selected.
    zero_bracket = [cells[0]['scaled_r_closed_interval'][1],
                    cells[1]['scaled_r_closed_interval'][0]]
    need(Q(zero_bracket[0]) < Q(zero_bracket[1]), 'nonempty_open_zero_bracket')
    return {
        'status': 'PASS_QUANTITATIVE_EVOLVED_FORMAL_SCALAR_MEAN_SIGNS',
        'arithmetic': 'Exact fractions.Fraction decisions; no floating-point acceptance.',
        'input_sha256': {name: pin['sha256'] for name, pin in pins.items()},
        'premise_authentication': 'Frozen accepted reports authenticated; norm join reconstructed. Parent all-point proofs remain premises.',
        'mean_projection': {
            'measure_mass': '1', 'W0': '1', 'functional_norm_exact': '1',
            'rule': '|integral w dmu| <= ||w||_L2(dmu) <= |(k,w)|_H',
            'selected_point_evaluation_factor_needed': False,
        },
        'response_join': {'q': str(q), 'X_exact_upper': str(exact_response),
                          'X_simple_strict_upper': str(simple_response),
                          'scalar_mean_x_derivative_strict_upper': str(mean_dx_upper)},
        'domain': {'ell': str(ell), 'scaled_x_closed_interval': list(map(str, x_interval)),
                   'physical_t_closed_interval': list(map(str, t_interval)),
                   'time_conversion': 't=ell*x', 'all_rectangles_in_D': True},
        'persistence': {'rule': '|mean_w U(x,r)-mean_w U(0,r)| <= x sup|mean_w U_x|',
                        'uniform_drift_strict_upper': str(drift), 'cells': cells},
        'zero_consequence': {'at_every_time_in_closed_interval': True,
                             'open_scaled_radius_bracket': zero_bracket,
                             'distinct_zeros_at_least_at_each_time': 1,
                             'continuous_zero_branch_claimed': False},
        'scope': {
            'actual_integral_defined_selected_formal_second_coefficient': True,
            'registered_gauge_and_initial_prescription_required': True,
            'quantitative_persistence_corollary_of_accepted_premises': True,
            'actual_U2_waveform_numerically_evaluated': False,
            'complete_numerical_residual_certified': False,
            'shell_boundary_sign_or_zero': False,
            'gauge_invariant_observable_or_force': False,
            'finite_amplitude_nonlinear_family_or_remainder': False,
            'external_novelty_claim': False,
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject', choices=sorted(MUTATIONS))
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse_existing_output')
    report = run(Path(__file__).resolve().parent, args.inject)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(report['status'])


if __name__ == '__main__':
    main()
