#!/usr/bin/env python3
"""Reassemble the unchanged-candidate C2 ledger with a refined initial cap.

All numerical decisions use exact rational arithmetic. No scientific source
code is imported. The complete retained ledger is reconstructed and checked
against the authenticated PR102 report before substituting the one improved
initial curvature error bound.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

PINS = {
    'initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'repair': '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
    'transfer': '4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3',
    'time': '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
    'response': 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
    'candidate': '44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f',
    'original_c2': '64967ab9afda32150153c934522bfb3d9784a6a96515fec6ef34343057d71a58',
    'original_fx': '63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60',
    'prior_c2': 'fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0',
    'prior_initial': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
}
checks = []
controls = []
hashes = {}


def need(ok, label):
    if not ok:
        raise RuntimeError(label)
    checks.append(label)


def load(path, role, digest=None):
    raw = Path(path).read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    need(sha == (digest or PINS[role]), 'authenticate ' + role)
    hashes[role] = sha
    return json.loads(raw)


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def reject(name, wrong, correct):
    need(wrong != correct, 'reject ' + name)
    controls.append(name)


def upper4(value):
    need(value > 0, 'positive error upper bound')
    unit = Q(1)
    while value >= 10000*unit:
        unit *= 10
    while value < 1000*unit:
        unit /= 10
    return (value // unit + 1)*unit


def validate_refinement(refined, labels, gram):
    need(refined['status'].startswith('PASS'), 'successful refined initial report')
    need(refined['same_candidate_W1_quadratic'] is True and refined['candidate_changed'] is False,
         'refinement changes estimates only, not candidate')
    need(refined['refined_u_rr_H_cap_is_upper_bound'] is True and
         refined['bound_quantity'] == 'sup_H_norm_actual_minus_smoothed_initial_u_rr',
         'refined input is an upper bound on the specified initial error')
    need(refined['labels'] == labels and refined['Gram'] == gram and len(labels) == 23 and
         min(map(Q, gram)) > 0, 'refined input has the same complete positive H norm')
    need(Q(refined['refined_uniform_initial_u_rr_H_upper']) >= 0,
         'refined initial upper bound is nonnegative')


def report_scope_controls(refined, labels, gram):
    variants = [
        ('NOT_PASS status', lambda d: d.update(status='NOT_PASS'), 'successful refined initial report'),
        ('changed candidate', lambda d: d.update(candidate_changed=True), 'refinement changes estimates only, not candidate'),
        ('different candidate identity', lambda d: d.update(same_candidate_W1_quadratic=False), 'refinement changes estimates only, not candidate'),
        ('lower bound input role', lambda d: d.update(refined_u_rr_H_cap_is_upper_bound=False), 'refined input is an upper bound'),
        ('point curvature quantity', lambda d: d.update(bound_quantity='point_signed_actual_minus_smoothed_initial_u_rr'), 'refined input is an upper bound'),
        ('incomplete all23', lambda d: d.update(labels=d['labels'][:-1]), 'refined input has the same complete positive H norm'),
        ('wrong positive Gram', lambda d: d['Gram'].__setitem__(0, '2'), 'refined input has the same complete positive H norm'),
        ('negative upper bound', lambda d: d.update(refined_uniform_initial_u_rr_H_upper='-1'), 'refined initial upper bound is nonnegative'),
    ]
    for name, mutate, expected in variants:
        changed = copy.deepcopy(refined)
        mutate(changed)
        start = len(checks)
        caught = False
        try:
            validate_refinement(changed, labels, gram)
        except RuntimeError as exc:
            caught = expected in str(exc)
        del checks[start:]
        need(caught, 'reject mutated report ' + name)
        controls.append('mutated report ' + name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--initial-refinement', type=Path, required=True)
    ap.add_argument('--initial-refinement-sha256', required=True)
    args = ap.parse_args()
    root = args.repo_root
    hi = root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference'
    co = root/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference'
    initial = load(hi/'INITIAL.json', 'initial')
    repair = load(hi/'REPAIR.json', 'repair')
    transfer = load(hi/'TRANSFER.json', 'transfer')
    time = load(hi/'TIME.json', 'time')
    response = load(root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json', 'response')
    candidate = load(root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json', 'candidate')
    original = load(co/'C2.json', 'original_c2')
    fx = load(co/'FX.json', 'original_fx')
    prior_c2 = load(root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/C2.json', 'prior_c2')
    prior_initial = load(root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json', 'prior_initial')
    refined = load(args.initial_refinement, 'initial_refinement', args.initial_refinement_sha256)
    need(original['status'] == 'PASS_CONDITIONAL_ORIGINAL_TARGET_C2_ERROR', 'closed retained original C2 gate')
    need(fx['status'].startswith('PASS') and fx['scope']['actual_original_full23_F_x_bound'] is True and
         fx['source_coordinate'] == 'x=t/ell', 'original unfiltered scaled-time source derivative')
    need(original['candidate'].startswith('W1_quadratic'), 'retained W1_quadratic candidate')
    need(original['boundary_error'] == {'g': '0', 'g_x': '0', 'M_x': '0'}, 'unchanged homogeneous error boundary')
    labels = transfer['labels']
    gram = transfer['Gram']
    validate_refinement(refined, labels, gram)
    report_scope_controls(refined, labels, gram)
    need(labels == initial['labels'] == repair['labels'] == time['labels'] == refined['labels'] and
         gram == initial['gram'] == repair['gram'] == time['Gram'] == refined['Gram'],
         'same all23 labels and positive Gram')
    need(len(labels) == 23 and min(map(Q, gram)) > 0, 'complete positive H norm')
    act = {k: Q(v) for k, v in response['actual_coefficient_H_norm_caps'].items()}
    frozen = {k: Q(v) for k, v in transfer['frozen_H_operator_caps'].items()}
    difference = {k: Q(v) for k, v in response['actual_coefficient_error_H_norm_caps'].items()}
    der = {k: Q(v) for k, v in time['global']['H_operator_caps'].items()}
    need(all(act[k] == frozen[k]+difference[k] for k in act), 'retained actual coefficient split')
    need(der['B_x'] == der['M_x'] == 0, 'actual B and boundary M time independent')
    span = Q(3, 2)
    c = act['A']+act['B']+Q(3, 4)*act['D']
    comm = der['A_x']+der['B_x']+Q(3, 4)*der['D_x']
    b = act['M']/4
    inverse = 1/(1-b)
    q = Q(response['retained_coupled_X_operator_cap'])
    need(span*c == Q(transfer['k_actual_bulk_X_operator_cap']) and
         b == Q(transfer['b_actual_boundary_X_operator_cap']) and
         inverse == Q(transfer['P_boundary_inverse_X_cap']), 'same Robin inverse and bulk operator')
    need(0 < b < 1 and span*c+b <= q == Q(3409, 8000000) < 1, 'same full contraction')
    need(comm == Q(time['coefficient_commutator_X_to_H_upper']), 'same complete time commutator')
    absolute = {k: Q(v) for k, v in candidate['initial_absolute_H_caps'].items()}
    direct = transfer['direct_reproof']
    Bs = Q(direct['new_smooth_Bs_X_cap'])
    Dfield = Q(direct['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
    C = Q(direct['new_actual_bulk_correction_X_cap'])
    E = Q(transfer['best_reproved_global_C1_error_X_cap'])
    F0 = Q(direct['new_reduced_forcing_defect_H_cap'])
    need(Bs == Q(candidate['free_wave_X_cap']) and Bs+Dfield == Q(candidate['W0_X_cap']),
         'unchanged absolute free-wave and W0 caps')
    need(F0 == c*C, 'unchanged original reduced forcing defect')
    free_second = absolute['u_rr']+absolute['v_r']
    f00 = sum(frozen[k]*absolute[field] for k, field in [('A', 'v'), ('B', 'u_r'), ('D', 'u')])
    f0x = (frozen['A']+frozen['B'])*free_second+frozen['D']*Bs
    need(f00 == Q(candidate['W0_initial_forcing_H_cap']) and
         f0x == Q(candidate['frozen_forcing_x_H_cap']), 'retained absolute second-jet forcing ledger')
    D_x = inverse*(f00+span*f0x+act['M']*Bs+Q(candidate['source_caps']['G_x_H']))
    r0_x_constant = (act['A']+act['B'])*D_x+(act['D']+comm)*Dfield+(
        difference['A']+difference['B'])*free_second+(difference['D']+comm)*Bs
    r0_initial = sum(difference[k]*absolute[field] for k, field in [('A', 'v'), ('B', 'u_r'), ('D', 'u')])+Q(candidate['source_caps']['F_H'])
    C_x_constant = inverse*(r0_initial+span*r0_x_constant)
    C_x_slope = inverse*span
    F1_constant = (act['A']+act['B'])*C_x_constant+(act['D']+comm)*C
    F1_slope = (act['A']+act['B'])*C_x_slope
    Fx = Q(fx['genuine_scaled_F_x_H_cap'])
    need(Fx == Q(original['actual_original_bulk_Fx_H_cap']) >= 0, 'same actual original Fx bound')
    F1 = F1_constant+F1_slope*Fx
    eu, eur, ev = map(Q, response['initial_smoothed_data_H_error_caps'])
    evr = Q(initial['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
    old_eurr = Q(repair['new_uniform_initial_u_rr_H_upper'])
    eurr = Q(refined['refined_uniform_initial_u_rr_H_upper'])
    need(0 <= eurr < old_eurr, 'strictly improved nonnegative initial curvature cap')
    need(max((eu+span*ev)/3, eur+ev) == Q(direct['new_initial_free_wave_X_error_cap']),
         'unchanged initial C1 comparison')
    common_initial = act['A']*ev+act['B']*eur+act['D']*eu
    slope = span*F1_slope/(1-q)
    def evaluate(curvature):
        acceleration = curvature+common_initial
        J = max((ev+span*acceleration)/3, evr+acceleration)
        z0 = (J+span*(F1_constant+comm*E))/(1-q)
        xx = z0+slope*Fx
        return {'acceleration': acceleration, 'J': J, 'z0': z0,
                'xx': xx, 'xr': xx, 'rr': xx+c*E+F0}
    old = evaluate(old_eurr)
    prior = evaluate(Q(prior_initial['refined_uniform_initial_u_rr_H_upper']))
    new = evaluate(eurr)
    need(old['acceleration'] == Q(original['initial_differentiated_acceleration_H_cap']) and
         old['J'] == Q(original['differentiated_initial_free_wave_X_cap']), 'reproduce retained initial C2 ledger')
    need(old['z0'] == Q(original['error_xx_xr_H_affine_in_Fx']['constant']) and
         slope == Q(original['error_xx_xr_H_affine_in_Fx']['Fx_coefficient']), 'reproduce retained affine C2 ledger')
    need(F1 == Q(original['original_reduced_defect_time_derivative_H_cap']) and
         D_x == Q(original['W0_minus_Bs_time_derivative_X_cap']), 'reproduce retained differentiated correction ledger')
    need(all(old[k] == Q(original['original_target_second_coordinate_error_H_caps'][k]) for k in ('xx', 'xr', 'rr')),
         'reproduce each retained C2 numeric cap exactly')
    need(all(prior[k] == Q(prior_c2['original_target_second_coordinate_error_H_caps'][k]) for k in ('xx', 'xr', 'rr')) and
         prior['J'] == Q(prior_c2['differentiated_initial_free_wave_X_cap']),
         'reproduce PR103 exact majorant and all three second-error caps')
    need(all(0 < new[k] < old[k] for k in ('xx', 'xr', 'rr')), 'all three C2 caps strictly improve')
    need(all(new[k] < prior[k] for k in ('xx', 'xr', 'rr')), 'all three C2 caps strictly improve from PR103')
    need(new['J'] == evr+new['acceleration'], 'differentiated radial initial branch still dominates')
    improvement = {k: old[k]/new[k] for k in ('xx', 'xr', 'rr')}
    inherited_floor = evaluate(Q(0))
    reject('replace absolute candidate u_rr by error u_rr', free_second-absolute['u_rr']+eurr, free_second)
    reject('omit unchanged initial v_r error', new['J']-evr, new['J'])
    reject('omit complete commutator', new['z0']-span*comm*E/(1-q), new['z0'])
    reject('omit complete feedback inverse', new['z0']*(1-q), new['z0'])
    reject('omit original Fx', new['xx']-slope*Fx, new['xx'])
    reject('reuse old initial curvature cap', old['xx'], new['xx'])
    reject('interpret error refinement as zero actual error', Q(0), new['xx'])
    report = {
        'status': 'PASS_COMBINED_INITIAL_ORIGINAL_TARGET_C2_ERROR',
        'candidate': original['candidate'], 'target': original['target'],
        'domain': original['domain'], 'norm': original['norm'],
        'input_sha256': hashes, 'labels': labels, 'Gram': gram,
        'candidate_changed': False, 'actual_sources_changed': False, 'actual_target_unchanged': True,
        'only_replaced_bound': 'initial error e0_rr; all absolute candidate, C1, v0_r, source and operator caps retained',
        'initial_error_H_caps': {'e0': eu, 'e0_r': eur, 'v0': ev, 'e0_rr': eurr, 'v0_r': evr},
        'retained_initial_curvature_error_H_cap': old_eurr,
        'initial_differentiated_acceleration_H_cap': new['acceleration'],
        'differentiated_initial_free_wave_X_cap': new['J'],
        'complete_original_C1_error_X_cap': E,
        'complete_original_reduced_defect_H_cap': F0,
        'original_reduced_defect_time_derivative_H_cap': F1,
        'actual_original_bulk_Fx_H_cap': Fx,
        'actual_Fx_contribution_to_second_error_upper': slope*Fx,
        'boundary_error': original['boundary_error'],
        'error_xx_xr_H_affine_in_Fx': {'constant': new['z0'], 'Fx_coefficient': slope},
        'error_rr_H_affine_in_Fx': {'constant': new['z0']+c*E+F0, 'Fx_coefficient': slope},
        'original_target_second_coordinate_error_H_caps': {k: new[k] for k in ('xx', 'xr', 'rr')},
        'original_target_second_coordinate_error_H_readable_strict_uppers': {k: upper4(new[k]) for k in ('xx', 'xr', 'rr')},
        'retained_to_refined_upper_bound_ratios': improvement,
        'PR103_to_refined_upper_bound_ratios': {k: prior[k]/new[k] for k in ('xx', 'xr', 'rr')},
        'PR103_second_coordinate_error_H_caps': {k: prior[k] for k in ('xx', 'xr', 'rr')},
        'initial_free_wave_fraction_of_differentiated_majorant_numerator': new['J']/(new['J']+span*(F1+comm*E)),
        'counterfactual_zero_e0_rr_ledger': {k: inherited_floor[k] for k in ('xx', 'xr', 'rr')},
        'counterfactual_scope': 'Sensitivity of this upper-bound formula with all other caps frozen; not an attainable error or a true-error lower bound',
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Conditional continuous physical-side coordinate C2 error bounds for the unchanged original-source W1_quadratic. Refining a certified bound does not change the candidate or true error. No relative physical accuracy, complete constraint propagation, nonlinear realization or cosmological evidence.'
    }
    print(json.dumps(encode(report), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
