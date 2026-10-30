#!/usr/bin/env python3
"""Combine independently refined paired and nominal-convolution upper caps."""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

PINS = {
    'old_higher': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'old_refinement': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
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
    sha = hashlib.sha256(raw).hexdigest()
    need(sha == (digest or PINS[role]), 'authenticate '+role)
    hashes[role] = sha
    return json.loads(raw)


def norm(caps, weights):
    square = sum((w*z*z for w, z in zip(weights, caps)), Q(0))
    need(square >= 0, 'nonnegative positive-Gram square')
    if not square:
        return Q(0)
    grid = 1 << 112
    k = math.isqrt((square*grid*grid).__ceil__())+1
    result = Q(k, grid)
    need(result*result >= square, 'outward complete H square root')
    return result


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def upper4(value):
    need(value > 0, 'positive continuous error bound for outward display')
    unit = Q(1)
    while value >= 10000*unit:
        unit *= 10
    while value < 1000*unit:
        unit /= 10
    return (value//unit+1)*unit


def common(d, labels, gram, role):
    need(d['status'].startswith('PASS'), 'successful '+role+' report')
    need(d['same_candidate_W1_quadratic'] is True and d['candidate_changed'] is False,
         role+' changes estimates only, not the candidate')
    need(d['labels'] == labels and d['Gram'] == gram, role+' same complete23 Gram')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--paired-report', type=Path, required=True)
    ap.add_argument('--paired-sha256', required=True)
    ap.add_argument('--nominal-report', type=Path, required=True)
    ap.add_argument('--nominal-sha256', required=True)
    a = ap.parse_args()
    hi = load(a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/INITIAL.json', 'old_higher')
    old = load(a.repo_root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json', 'old_refinement')
    paired = load(a.paired_report, 'paired', a.paired_sha256)
    nominal = load(a.nominal_report, 'nominal', a.nominal_sha256)
    labels, gram = old['labels'], old['Gram']
    weights = list(map(Q, gram))
    need(len(labels) == len(weights) == 23 and min(weights) > 0, 'complete positive full23 H')
    need(labels == hi['labels'] and gram == hi['gram'], 'same original initial norm')
    common(paired, labels, gram, 'paired')
    common(nominal, labels, gram, 'nominal')
    need(paired['paired_u_rr_H_cap_is_upper_bound'] is True and
         paired['bound_quantity'] == 'sup_H_norm_paired_actual_minus_nominal_initial_u_rr',
         'paired complete initial-error upper role')
    need(nominal['nominal_M3_upper'] is True and nominal['actual_order_six_bound_claimed'] is False and
         nominal['bound_quantity'] == 'piecewise_nominal_u_rrr_component_upper',
         'explicit nominal third derivative upper role, no actual sixth-order assumption')
    m3 = list(map(Q, nominal['refined_nominal_M3u_component_caps']))
    old_m3 = list(map(Q, hi['piecewise_nominal_M3u']))
    need(len(m3) == 23 and all(0 <= n <= o for n, o in zip(m3, old_m3)),
         'componentwise nonnegative no-larger nominal M3 bounds')
    du, dv = Q(1, 2**30), Q(1, 2048)
    need(Q(paired['fixed_displacement_kernel_width']) == Q(nominal['fixed_displacement_kernel_width']) == du and
         Q(paired['fixed_velocity_kernel_width']) == Q(nominal['fixed_velocity_kernel_width']) == dv,
         'same fixed convolution widths')
    jumps = hi['curvature_jumps']
    need([Q(row['r']) for row in jumps] == [Q(65, 64), Q(66, 64), Q(7, 4)] and
         all(Q(v) == 0 for v in jumps[0]['jump_abs_component_caps']), 'same isolated joins and zero first jump')
    inner = [max(Q(row['jump_abs_component_caps'][i]) for row in jumps[:2]) for i in range(23)]
    need(inner == list(map(Q, old['unchanged_inner_jump_component_caps'])), 'same quadratic exterior removes only outer jump')
    smooth = [Q(35, 128)*du*m+j/2 for m, j in zip(m3, inner)]
    smoothing_h = norm(smooth, weights)
    need(smooth == list(map(Q, nominal['refined_curvature_smoothing_component_caps'])) and
         smoothing_h == Q(nominal['refined_curvature_smoothing_H_cap']), 'independently rebuild entire smoothing charge')
    channels = [[Q(v) for v in row] for row in paired['paired_channel_H_caps']]
    need(len(channels) == 3 and all(len(row) == 2 and min(row) >= 0 for row in channels), 'three complete paired channels')
    paired_h = sum((row[0] for row in channels), Q(0))
    need(paired_h == Q(paired['paired_u_rr_H_upper']), 'sum all paired channel suprema')
    need(paired_h <= Q(old['paired_u_rr_H_upper']) and
         smoothing_h <= Q(old['unchanged_curvature_smoothing_H_cap']), 'both majorant parts no larger than PR103')
    new = paired_h+smoothing_h
    previous = Q(old['refined_uniform_initial_u_rr_H_upper'])
    need(0 < new < previous, 'strictly tighter combined uniform initial-curvature error cap')
    velocity = Q(paired['unchanged_full_initial_v_r_error_H_upper'])
    need(velocity == Q(old['unchanged_full_initial_v_r_error_H_upper']), 'complete velocity error retained')
    # Each mutation goes through the actual acceptance predicate, rather than
    # relying on a status string or a coincidental numerical inequality.
    for name, key, value in (('NOT_PASS status', 'status', 'NOT_PASS'),
                             ('different candidate', 'candidate_changed', True),
                             ('wrong full23 norm', 'Gram', ['2']+gram[1:])):
        bad = copy.deepcopy(paired)
        bad[key] = value
        begin = len(checks)
        caught = False
        try:
            common(bad, labels, gram, 'paired')
        except RuntimeError:
            caught = True
        del checks[begin:]
        need(caught, 'reject '+name)
        controls.append(name)
    for name, wrong in (('omit nominal smoothing charge', paired_h),
                        ('omit paired initial comparison', smoothing_h),
                        ('maximum replaces sum of channels', max(row[0] for row in channels)+smoothing_h)):
        need(wrong != new, 'reject '+name)
        controls.append(name)
    report = {'status': 'PASS_COMBINED_UNCHANGED_INITIAL_CURVATURE_REFINEMENT',
        'candidate': old['candidate'], 'same_candidate_W1_quadratic': True, 'candidate_changed': False,
        'target_unchanged': True, 'initial_data_unchanged': True, 'F_G_and_coefficients_unchanged': True,
        'refined_u_rr_H_cap_is_upper_bound': True,
        'bound_quantity': 'sup_H_norm_actual_minus_smoothed_initial_u_rr',
        'labels': labels, 'Gram': gram, 'gram': gram, 'ell': old['ell'],
        'fixed_displacement_kernel_width': du, 'fixed_velocity_kernel_width': dv,
        'paired_channel_H_caps': channels, 'paired_u_rr_H_upper': paired_h,
        'refined_nominal_M3u_component_caps': m3, 'retained_inner_jump_component_caps': inner,
        'refined_curvature_smoothing_component_caps': smooth, 'refined_curvature_smoothing_H_cap': smoothing_h,
        'refined_uniform_initial_u_rr_H_upper': new,
        'refined_uniform_initial_u_rr_H_readable_strict_upper': upper4(new),
        'previous_PR103_initial_u_rr_H_upper': previous,
        'upper_certificate_reduction_factor_from_PR103': previous/new,
        'unchanged_full_initial_v_r_error_H_upper': velocity,
        'input_sha256': hashes, 'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Continuous initial error upper bound for unchanged W1_quadratic. Only interval majorants change; no smaller actual error, new candidate/source, physical accuracy or nonlinear/cosmological realization.'}
    print(json.dumps(encode(report), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
