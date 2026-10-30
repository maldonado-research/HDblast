#!/usr/bin/env python3
"""Conditional global coordinate second-derivative seminorm normalization.

No scientific imports. The positive denominator is an authenticated actual
initial endpoint derivative; it is never the old affine candidate's error.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

PINS = {
    'actual_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'independent_endpoint': 'e21d2cb172893495b557bbcdcee198d33d21a7ba287c5215d7319eff3d9dc9cc',
    'prior_initial': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
    'prior_c2': 'fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0',
    'retained_join_errors': 'e1d0d6e3d8650cbf402ccb4f5a87fd20c20db651e3bd323c68f0bdea0219dd35',
}
checks = []
controls = []
hashes = {}


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)
    checks.append(reason)


def load(path, role, digest=None):
    raw = Path(path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    need(actual == (digest or PINS[role]), 'authenticate '+role)
    hashes[role] = actual
    return json.loads(raw)


def interval(row):
    a, b = map(Q, row)
    need(a <= b, 'ordered signed derivative enclosure')
    return a, b


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def upper4(value):
    need(value > 0, 'positive coordinate quotient upper')
    unit = Q(1)
    while value >= 10000*unit:
        unit *= 10
    while value < 1000*unit:
        unit /= 10
    return (value//unit+1)*unit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--combined-report', type=Path, required=True)
    ap.add_argument('--combined-sha256', required=True)
    ap.add_argument('--c2-report', type=Path, required=True)
    ap.add_argument('--c2-sha256', required=True)
    ap.add_argument('--inject-error', choices=('affine-error-denominator', 'wrong-Gram', 'drop-trace',
                    'pointwise-relative', 'separate-xx-denominator', 'full-C2-norm', 'physical-accuracy'))
    a = ap.parse_args()
    root = a.repo_root
    initial = load(root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/INITIAL.json', 'actual_initial')
    independent = load(root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/ENDPOINT.json', 'independent_endpoint')
    prior_i = load(root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json', 'prior_initial')
    prior_c2 = load(root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/C2.json', 'prior_c2')
    joins = load(root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/JOINS.json', 'retained_join_errors')
    combined = load(a.combined_report, 'combined', a.combined_sha256)
    c2 = load(a.c2_report, 'c2', a.c2_sha256)
    need(combined['status'] == 'PASS_COMBINED_UNCHANGED_INITIAL_CURVATURE_REFINEMENT' and
         c2['status'] == 'PASS_COMBINED_INITIAL_ORIGINAL_TARGET_C2_ERROR', 'successful uniform error inputs')
    need(combined['same_candidate_W1_quadratic'] is True and combined['candidate_changed'] is False and
         c2['candidate_changed'] is False and c2['actual_sources_changed'] is False, 'same original target and W1_quadratic')
    need(combined['refined_u_rr_H_cap_is_upper_bound'] is True and
         combined['bound_quantity'] == 'sup_H_norm_actual_minus_smoothed_initial_u_rr', 'uniform initial error numerator role')
    labels, gram = combined['labels'], combined['Gram']
    need(labels == c2['labels'] == initial['labels'] == prior_i['labels'] == prior_c2['labels'] and
         gram == c2['Gram'] == initial['gram'] == prior_i['Gram'] == prior_c2['Gram'], 'same all23 positive norm across initial and evolved estimates')
    weights = list(map(Q, gram))
    if a.inject_error == 'wrong-Gram':
        weights[12] *= 2
    need(len(labels) == len(weights) == 23 and min(weights) > 0, 'complete positive angular H')
    need(labels[12] == 'w:0' and weights[12] == 1, 'actual scalar mean has unit Gram weight')
    need(Q(initial['ell']) == Q(combined['ell']) == Q(1, 10**6), 'same scaled radial coordinate convention')
    actual_row = initial['endpoint']['rows'][12]
    direct_row = independent['rows'][12]
    need(actual_row['label'] == direct_row['label'] == 'w:0' and Q(direct_row['input_Gram']) == 1,
         'paired and independent actual endpoint component')
    quantity = 'actual_u_rr'
    if a.inject_error == 'affine-error-denominator':
        quantity = 'actual_minus_smoothed_u_rr_at_fixed_width'
    need(quantity == 'actual_u_rr', 'denominator uses actual target derivative rather than affine approximation error')
    paired_iv = interval(actual_row[quantity])
    independent_iv = interval(direct_row['actual_u_rr'])
    need(max(paired_iv[0], independent_iv[0]) <= min(paired_iv[1], independent_iv[1]),
         'independent and paired actual endpoint intervals overlap')
    scale = Q('3.249e33')
    need(paired_iv[1] < -scale and independent_iv[1] < -scale,
         'strict positive target endpoint H-scale in both actual-jet routes')
    R, x = Q(7, 4), Q(0)
    need(Q(initial['endpoint']['r']) == R and len(independent['endpoint']) == 2 and
         list(map(Q, independent['endpoint'])) == [0, R], 'same actual outer initial endpoint')
    need(0 <= x <= Q(3, 2) and max(Q(0), x-Q(5, 4)) <= R <= Q(7, 4)-x,
         'actual initial corner belongs to declared closed Dplus')
    need(c2['domain'] == prior_c2['domain'] == 'D+:0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x',
         'same evolved physical-side domain')
    trace = a.inject_error != 'drop-trace'
    need(trace, 'retained classical target initial rr trace is an explicit required premise')
    allowed_claim = 'global_rr_sup_and_max_second_derivative_seminorm'
    requested = {
        'pointwise-relative': 'pointwise_relative',
        'separate-xx-denominator': 'relative_xx_with_rr_denominator',
        'full-C2-norm': 'full_C2_norm',
        'physical-accuracy': 'physical_accuracy',
    }.get(a.inject_error, allowed_claim)
    need(requested == allowed_claim, 'only declared global coordinate rr and max-second seminorm ratios justified')
    initial_cap = Q(combined['refined_uniform_initial_u_rr_H_upper'])
    need(initial_cap == Q(c2['initial_error_H_caps']['e0_rr']) > 0, 'same uniform initial error numerator')
    caps = {key: Q(c2['original_target_second_coordinate_error_H_caps'][key]) for key in ('xx', 'xr', 'rr')}
    old_caps = {key: Q(prior_c2['original_target_second_coordinate_error_H_caps'][key]) for key in ('xx', 'xr', 'rr')}
    need(min(caps.values()) > 0 and all(caps[k] < old_caps[k] for k in caps), 'all three complete coordinate error numerators improve')
    ratios = {'initial_rr_H_sup': initial_cap/scale,
              'evolved_rr_H_sup': caps['rr']/scale,
              'evolved_max_second_derivative_H_sup_seminorm': max(caps.values())/scale}
    baseline = {'initial_rr_H_sup': Q(prior_i['refined_uniform_initial_u_rr_H_upper'])/scale,
                'evolved_rr_H_sup': old_caps['rr']/scale,
                'evolved_max_second_derivative_H_sup_seminorm': max(old_caps.values())/scale}
    need(baseline['initial_rr_H_sup'] < Q('3.648e-6') and
         baseline['evolved_max_second_derivative_H_sup_seminorm'] < Q('3.654e-6'),
         'reproduce safe PR103 coordinate normalization baseline')
    need(all(0 < ratios[k] < baseline[k] for k in ratios), 'strictly tighter conditional global coordinate quotient bounds')
    need(joins['status'] == 'PASS_INNER_JOIN_ACTUAL_INITIAL_CURVATURE_ERROR' and
         joins['initial_candidate_unchanged'] is True and joins['actual_target_unchanged'] is True and
         joins['labels'] == labels and joins['positive_Gram'] == gram, 'same actual join target candidate and full norm')
    first_join = joins['joins'][0]
    local = first_join['rows'][12]
    need(Q(first_join['r']) == Q(65, 64) and local['label'] == 'w:0', 'registered local scalar-mean calibration point')
    local_target = interval(local['actual_u_rr'])
    local_error = interval(local['actual_minus_smoothed_u_rr'])
    need(local_target[0] > 0 and local_error[0] > 0, 'nonzero actual initial component and positive approximation error')
    local_ratio_lower = local_error[0]/local_target[1]
    need(local_ratio_lower > Q(999, 1000), 'local initial component relative error exceeds99.9percent despite small global norm ratio')
    actual_H_square_upper = Q(0)
    error_H_square_lower = Q(0)
    for i, row in enumerate(first_join['rows']):
        need(row['label'] == labels[i] and Q(row['positive_Gram_weight']) == weights[i],
             'same full23 H calibration row '+labels[i])
        au, er = interval(row['actual_u_rr']), interval(row['actual_minus_smoothed_u_rr'])
        actual_H_square_upper += weights[i]*max(abs(v) for v in au)**2
        minimum = Q(0) if er[0] <= 0 <= er[1] else min(abs(v) for v in er)
        error_H_square_lower += weights[i]*minimum**2
    H_error_lower = Q(first_join['fixed_width_error_H_norm_strict_lower'])
    need(H_error_lower > 0 and H_error_lower**2 < error_H_square_lower and actual_H_square_upper > 0,
         'independent full23 weighted actual/error point norm bounds')
    grid = 1 << 112
    H_actual_upper = Q(math.isqrt((actual_H_square_upper*grid*grid).__ceil__())+1, grid)
    need(H_actual_upper**2 >= actual_H_square_upper, 'outward pointwise actual H norm upper')
    local_H_ratio_lower = H_error_lower/H_actual_upper
    need(local_H_ratio_lower > Q(999, 1000),
         'same-H pointwise initial relative error exceeds99.9percent')
    calibration = {'x': Q(0), 'r': Q(65, 64), 'component': 'w:0',
        'quantity': 'initial radial second derivative', 'actual_target_interval': local_target,
        'actual_minus_candidate_interval': local_error,
        'component_relative_absolute_error_lower_bound': local_ratio_lower,
        'component_relative_absolute_error_readable_strict_lower': Q(999, 1000),
        'scope': 'One initial coordinate component ratio; no physical model defect. It disproves interpreting a small ratio of global suprema as a uniform pointwise component-relative error bound.'}
    full_H_calibration = {'x': Q(0), 'r': Q(65, 64),
        'quantity': 'initial radial second derivative in the same complete positive23 H norm',
        'actual_target_H_norm_upper': H_actual_upper,
        'actual_error_H_norm_strict_lower': H_error_lower,
        'pointwise_H_relative_error_strict_lower': local_H_ratio_lower,
        'pointwise_H_relative_error_readable_strict_lower': Q(999, 1000),
        'actual_H_squared_upper_from_all23_intervals': actual_H_square_upper,
        'error_H_squared_lower_from_all23_intervals': error_H_square_lower,
        'nonzero_target_denominator': 'actual scalar mean radial second derivative is positive',
        'scope': 'Concrete lower bound on the pointwise initial coordinate H-norm relative error, using the same H norm as the global quotient. It excludes a uniform pointwise99.9percent accuracy claim despite a small global quotient; no physical or nonlinear obstruction is implied.'}
    # Exact polynomial counterexamples distinguish these mathematical claims.
    # For p(x,r)=1+2x+3r all second derivatives vanish, while p(0,0)=1.
    affine = {(0, 0): Q(1), (1, 0): Q(2), (0, 1): Q(3)}
    def derivative(p, axis):
        result = {}
        for powers, coefficient in p.items():
            n = powers[axis]
            if n:
                powers = list(powers)
                powers[axis] -= 1
                result[tuple(powers)] = n*coefficient
        return result
    need(all(not derivative(derivative(affine, i), j) for i, j in ((0, 0), (0, 1), (1, 1))) and
         affine[(0, 0)] == 1, 'affine polynomial has zero second seminorm and nonzero value')
    controls.append('full C2 relabeling excluded by seminorm definition')
    # U=r^3/6 and e=epsilon*r^2/2 have U_rr=r and e_rr=epsilon:
    # the global quotient is epsilon on [0,1], but the target vanishes at0.
    eps = Q(1, 10**6)
    target_poly, error_poly = {(0, 3): Q(1, 6)}, {(0, 2): eps/2}
    target_rr = derivative(derivative(target_poly, 1), 1)
    error_rr = derivative(derivative(error_poly, 1), 1)
    need(target_rr == {(0, 1): Q(1)} and error_rr == {(0, 0): eps},
         'global derivative quotient does not ensure local nonzero target')
    controls.append('pointwise relabeling excluded by exact zero-target fixture')
    for name, wrong, correct in (
        ('inverted target scale', initial_cap*scale, ratios['initial_rr_H_sup']),
        ('replace complete second error by minimum', min(caps.values())/scale, max(caps.values())/scale),
        ('use old affine error lower as target scale', caps['rr']/Q('1.6245e33'), ratios['evolved_rr_H_sup']),
    ):
        need(wrong != correct, 'reject '+name)
        controls.append(name)
    report = {'status': 'PASS_GLOBAL_COORDINATE_SECOND_SEMINORM_NORMALIZATION',
        'candidate': c2['candidate'], 'candidate_changed': False, 'actual_sources_changed': False,
        'target': c2['target'], 'domain': c2['domain'], 'labels': labels, 'Gram': gram,
        'input_sha256': hashes,
        'normalization_quantity': allowed_claim,
        'denominator_role': 'actual target coordinate derivative supremum; lower bounded by actual initial scalar-mean rr endpoint',
        'target_scale_strict_lower': scale,
        'actual_initial_endpoint_r': R,
        'paired_actual_endpoint_u_rr_interval': paired_iv,
        'independent_actual_endpoint_u_rr_interval': independent_iv,
        'global_coordinate_ratios': ratios,
        'global_coordinate_ratio_readable_strict_uppers': {k: upper4(v) for k, v in ratios.items()},
        'baseline_PR103_global_coordinate_ratios': baseline,
        'ratio_bound_improvement_factors': {k: baseline[k]/ratios[k] for k in ratios},
        'local_component_calibration': calibration,
        'local_same_H_calibration': full_H_calibration,
        'initial_denominator_definition': 'sup_[0,7/4] ||actual U0_rr||_H',
        'evolved_rr_denominator_definition': 'sup_Dplus ||actual U_rr||_H',
        'maximum_second_seminorm_definition': 'max_ab_in_xx_xr_rr sup_Dplus ||V_ab||_H',
        'required_premises': {'same_actual_target_initial_trace_Urr_equals_U0rr': True,
                              'classical_physical_side_C2_with_one_sided_corner_trace': True,
                              'registered_coordinate_and_positive_full23_H_norm': True},
        'claimed_conclusions': {'global_coordinate_rr_sup_ratio': True,
                                'global_coordinate_maximum_second_derivative_seminorm_ratio': True,
                                'uniform_pointwise_relative_error_upper_bound': False,
                                'one_initial_point_same_H_relative_error_lower_bound': True,
                                'one_initial_component_relative_error_lower_bound': True,
                                'separate_relative_xx_or_xr_sup_ratios': False,
                                'full_C2_norm_relative_error': False,
                                'physical_covariant_or_observational_relative_accuracy': False,
                                'finite_amplitude_nonlinear_or_cosmological_realization': False},
        'denominator_limitation': 'The global scale can be dominated by the large initial outer-endpoint derivative and can mask large relative errors in smaller components or local regions.',
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Conditional ratio of global coordinate H-sup seminorms for the unchanged formal second coefficient, with a concrete local initial error lower bound. No uniform pointwise, invariant, covariant, physical, finite-amplitude nonlinear or cosmological accuracy claim.'}
    print(json.dumps(encode(report), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
