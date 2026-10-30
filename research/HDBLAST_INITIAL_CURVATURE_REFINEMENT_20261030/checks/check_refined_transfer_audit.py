#!/usr/bin/env python3
"""Independent delta identity audit of the unchanged-candidate C2 transfer.

Uses the authenticated old theorem and its affine dependence on the initial
free-wave majorant; imports no producer or transfer checker code.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def load(path, digest):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'authenticated report')
    return json.loads(raw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--refinement', type=Path, required=True)
    ap.add_argument('--refinement-sha256', required=True)
    ap.add_argument('--transfer', type=Path, required=True)
    ap.add_argument('--transfer-sha256', required=True)
    a = ap.parse_args()
    old = load(a.repo_root/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/C2.json',
               '64967ab9afda32150153c934522bfb3d9784a6a96515fec6ef34343057d71a58')
    response = load(a.repo_root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json',
                    'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84')
    r = load(a.refinement, a.refinement_sha256)
    t = load(a.transfer, a.transfer_sha256)
    checks = []

    def ck(ok, reason):
        need(ok, reason)
        checks.append(reason)

    ck(old['status'] == 'PASS_CONDITIONAL_ORIGINAL_TARGET_C2_ERROR', 'retained closed C2 theorem')
    ck(r['status'].startswith('PASS') and t['status'] == 'PASS_REFINED_INITIAL_ORIGINAL_TARGET_C2_ERROR', 'successful new inputs')
    ck(r['same_candidate_W1_quadratic'] and not r['candidate_changed'] and
       not t['candidate_changed'] and not t['actual_sources_changed'], 'same candidate and original sources')
    ck(t['candidate'] == old['candidate'] and t['target'] == old['target'] and
       t['domain'] == old['domain'] and t['norm'] == old['norm'], 'same target norm and physical domain')
    ck(t['labels'] == r['labels'] and t['Gram'] == r['Gram'] and
       len(t['labels']) == 23 and min(map(Q, t['Gram'])) > 0, 'complete same positive angular H norm')
    olddata = {k: Q(v) for k, v in old['initial_error_H_caps'].items()}
    data = {k: Q(v) for k, v in t['initial_error_H_caps'].items()}
    new_curvature = Q(r['refined_uniform_initial_u_rr_H_upper'])
    ck(r['refined_u_rr_H_cap_is_upper_bound'] is True and
       r['bound_quantity'] == 'sup_H_norm_actual_minus_smoothed_initial_u_rr', 'upper bound has the correct error role')
    ck(data['e0_rr'] == new_curvature and 0 <= new_curvature < olddata['e0_rr'], 'one improved initial error bound')
    ck(new_curvature < Q('1.185e28'), 'public readable initial bound')
    ck(all(data[k] == v for k, v in olddata.items() if k != 'e0_rr'), 'other initial error bounds retained')
    q = Q(response['retained_coupled_X_operator_cap'])
    ck(q == Q(3409, 8000000), 'same feedback contraction')
    old_accel = Q(old['initial_differentiated_acceleration_H_cap'])
    new_accel = old_accel - olddata['e0_rr'] + new_curvature
    J = max((data['v0'] + Q(3, 2)*new_accel)/3, data['v0_r'] + new_accel)
    old_J = Q(old['differentiated_initial_free_wave_X_cap'])
    ck(new_accel == Q(t['initial_differentiated_acceleration_H_cap']), 'independent acceleration delta identity')
    ck(J == Q(t['differentiated_initial_free_wave_X_cap']) and J < old_J, 'independent free-wave max branch')
    shift = (J-old_J)/(1-q)
    results = {}
    for key in ('xx', 'xr', 'rr'):
        previous = Q(old['original_target_second_coordinate_error_H_caps'][key])
        actual = Q(t['original_target_second_coordinate_error_H_caps'][key])
        expected = previous + shift
        ck(0 < actual == expected < previous, 'exact C2 delta identity '+key)
        ck(Q(t['retained_to_refined_upper_bound_ratios'][key]) == previous/actual,
           'exact upper-bound ratio '+key)
        ck(actual < Q('1.187e28') and previous > 524000*actual,
           'public readable C2 bound and upper-certificate reduction '+key)
        ck(actual < Q(t['original_target_second_coordinate_error_H_readable_strict_uppers'][key]),
           'outward readable upper '+key)
        results[key] = str(actual)
    retained = ('actual_original_bulk_Fx_H_cap', 'actual_Fx_contribution_to_second_error_upper',
                'complete_original_C1_error_X_cap', 'complete_original_reduced_defect_H_cap',
                'original_reduced_defect_time_derivative_H_cap', 'boundary_error')
    ck(all(t[k] == old[k] for k in retained), 'source and correction derivative ledgers retained exactly')
    ck(t['boundary_error'] == {'g': '0', 'g_x': '0', 'M_x': '0'}, 'homogeneous unchanged error boundary')
    controls = []
    for name, wrong in (
        ('omit_feedback_inverse_in_delta', J-old_J),
        ('subtract_initial_bound_without_max_branch', new_curvature-olddata['e0_rr']),
        ('reuse_old_bound', Q(0)),
    ):
        ck(wrong != shift, 'reject '+name)
        controls.append(name)
    result = {'status': 'PASS_INDEPENDENT_REFINED_C2_DELTA_AUDIT',
        'method': 'old authenticated exact theorem plus independently evaluated initial acceleration/max and affine delta identity',
        'initial_acceleration': str(new_accel), 'initial_free_wave_majorant': str(J),
        'each_second_error_upper_change': str(shift), 'second_error_H_caps': results,
        'refinement_sha256': a.refinement_sha256, 'transfer_sha256': a.transfer_sha256,
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Internal audit of conditional absolute coordinate error bounds; an upper-bound reduction is not a reduction or a lower bound of actual error.'}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
