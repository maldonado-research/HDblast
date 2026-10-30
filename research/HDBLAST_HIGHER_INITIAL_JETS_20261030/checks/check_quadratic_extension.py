#!/usr/bin/env python3
"""A distinct curvature-matched initial extension, with exact error charges.

The old fixed candidate is preserved. The quadratic extension coefficient is
the same exact analytic nominal left curvature, enclosed by pinned reports.
No midpoint is chosen as a replacement for that function value.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path

INITIAL_PIN = 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207'
ENDPOINT_PIN = 'e21d2cb172893495b557bbcdcee198d33d21a7ba287c5215d7319eff3d9dc9cc'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def load(path, digest):
    raw = path.read_bytes()
    need(digest is not None and hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return json.loads(raw)


def iv(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(Q, z)
    need(a <= b, 'ordered interval')
    return a, b


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def wide(a, e):
    need(e >= 0, 'nonnegative remainder')
    return a[0]-e, a[1]+e


def scale(a, c):
    return (a[0]*c, a[1]*c) if c >= 0 else (a[1]*c, a[0]*c)


def norm(caps, weights):
    q = sum(w*z*z for w, z in zip(weights, caps))
    grid = 1 << 112
    k = math.isqrt(q.numerator*grid*grid//q.denominator)
    if Q(k, grid)**2 < q:
        k += 1
    result = Q(k, grid)
    need(result*result >= q, 'outward rational norm')
    return result


def enc(z):
    return list(map(str, z))


ap = argparse.ArgumentParser()
ap.add_argument('--initial-report', type=Path, required=True)
ap.add_argument('--endpoint-report', type=Path, required=True)
args = ap.parse_args()
d = load(args.initial_report, INITIAL_PIN)
a = load(args.endpoint_report, ENDPOINT_PIN)
need(d['status'] == a['status'] == 'PASS' and d['initial_data_unchanged'] is True, 'successful old-candidate inputs')
labels, weights = d['labels'], list(map(Q, d['gram']))
delta, R = Q(d['fixed_displacement_kernel_width']), Q(d['endpoint']['r'])
need(len(labels) == len(weights) == 23 and min(weights) > 0 and delta == Q(1, 2**30) and R == Q(7, 4), 'same full23 state and scales')
rows = d['endpoint']['rows']
independent = {r['label']: r for r in a['rows']}
need(len(rows) == len(independent) == 23 and set(independent) == set(labels), 'complete independent endpoint data')

# Exact moments of the retained positive even compact kernel.
polynomial = [(0, Q(35, 32)), (2, -Q(105, 32)), (4, Q(105, 32)), (6, -Q(35, 32))]
half_moment = lambda p: sum(c/Q(n+p+1) for n, c in polynomial)
need(half_moment(0) == Q(1, 2) and half_moment(1) == Q(35, 256) and
     half_moment(2)/2 == Q(1, 36), 'exact half-kernel moments')
# A polynomial with constant interior curvature supplies an exact independent
# fixture: affine continuation loses half its curvature; matching quadratic
# continuation restores it. The linear and constant parts contribute zero.
for curvature in (Q(-7, 3), Q(0), Q(11, 5)):
    affine = curvature*half_moment(0)
    repaired = affine+curvature*half_moment(0)
    need(repaired == curvature and affine == curvature/2, 'exact quadratic reproduction at the extension join')
need(half_moment(1) != Q(35, 128), 'doubled first half-moment corruption rejected')

kappa_caps, out = [], []
for i, row in enumerate(rows):
    need(row['label'] == labels[i], 'same modal order')
    kappa = iv(row['nominal_left_u_rr'])
    paired_error = iv(row['paired_actual_minus_nominal_u_rr'])
    remainder = Q(row['fixed_width_half_kernel_remainder_abs_upper'])
    need(remainder == Q(d['piecewise_nominal_M3u'][i])*delta*Q(35, 256), 'half-kernel Lipschitz remainder')
    old_error = iv(row['actual_minus_smoothed_u_rr_at_fixed_width'])
    audited_old = iv(independent[labels[i]]['actual_minus_affine_smoothed_u_rr'])
    need(max(old_error[0], audited_old[0]) <= min(old_error[1], audited_old[1]), 'independent old endpoint errors overlap')
    # Adding (kappa/2) to the smoothed curvature cancels the old half-curvature
    # term algebraically, avoiding subtraction of two enormous intervals.
    new_error = wide(paired_error, remainder)
    kappa_caps.append(max(abs(z) for z in kappa))
    out.append({'label': labels[i], 'exact_nominal_curvature_interval': enc(kappa),
        'old_actual_minus_smoothed_u_rr': enc(old_error),
        'new_actual_minus_quadratic_smoothed_u_rr': enc(new_error),
        'change_initial_value_at_R': enc(scale(kappa, delta**2/36)),
        'change_initial_radial_derivative_at_R': enc(scale(kappa, Q(35, 256)*delta)),
        'change_initial_curvature_at_R': enc(scale(kappa, Q(1, 2)))})
need(labels[12] == 'w:0', 'registered scalar mean')
old_mean = iv(out[12]['old_actual_minus_smoothed_u_rr'])
new_mean = iv(out[12]['new_actual_minus_quadratic_smoothed_u_rr'])
old_lower = Q('1.6245e33')
new_upper = Q('1.820e27')
need(old_mean[1] < -old_lower and max(abs(z) for z in new_mean) < new_upper and
     old_lower > 800000*new_upper, 'genuine endpoint initial-error improvement, not comparison of two upper bounds')
kappa_H = norm(kappa_caps, weights)
value_change = kappa_H*delta**2/36
derivative_change = kappa_H*Q(35, 256)*delta
free_X_change = max(value_change/3, derivative_change)
need(free_X_change == derivative_change and Q(9, 8) < R-delta and Q(5, 4) < R-delta, 'free norm charge and protected causal regions')

# The nominal second derivative now matches across R. Only the two inner
# cutoff jumps remain. Its third derivative is zero on the quadratic exterior.
jumps = d['curvature_jumps']
need([Q(row['r']) for row in jumps] == [Q(65, 64), Q(66, 64), R], 'all old curvature joins identified')
inner_caps = [max(Q(row['jump_abs_component_caps'][i]) for row in jumps[:2]) for i in range(23)]
new_smoothing_caps = [Q(m)*Q(35, 128)*delta+j/2 for m, j in zip(d['piecewise_nominal_M3u'], inner_caps)]
channels = d['continuous']['paired_channel_H_caps']
need(len(channels) == 3 and d['continuous']['continuous_cells'] == 48, 'same continuous paired comparison')
new_uniform_initial_u2 = sum(Q(row[0]) for row in channels)+norm(new_smoothing_caps, weights)
old_uniform_initial_u2 = Q(d['continuous']['actual_minus_smoothed_u_rr_H_upper'])
need(new_uniform_initial_u2 <= old_uniform_initial_u2, 'removing the exterior jump cannot enlarge this initial-data upper certificate')
need(new_uniform_initial_u2 < Q('6.219e33'), 'strict repaired uniform initial-curvature error upper bound')
print(json.dumps({'status': 'PASS', 'new_candidate': 'W1_quadratic, a distinct auxiliary finite update with curvature-matched exterior displacement extension before the same convolution',
    'definition': 'u_ext_q(R+h)=u_nom(R)+u_nom_r(R) h+(1/2)u_nom_rr(R-) h^2 for h>=0; same nominal function on [0,R], same parity, velocity, delta_u, actual coefficients and fixed F/G',
    'old_candidate_preserved': True, 'underlying_actual_target_unchanged': True,
    'R': str(R), 'delta_u': str(delta), 'labels': labels, 'gram': list(map(str, weights)),
    'curvature_H_upper': str(kappa_H), 'initial_value_change_H_upper': str(value_change),
    'initial_radial_derivative_change_H_upper': str(derivative_change),
    'free_wave_X_change_upper': str(free_X_change),
    'difference_support': 'initial difference zero for r<=R-delta_u; reflected full source/Volterra differences vanish when x+r<R-delta_u',
    'old_selected_point_and_retained_shell_values_unchanged': True,
    'endpoint_initial_scalar_mean_old_abs_strict_lower': str(old_lower),
    'endpoint_initial_scalar_mean_new_abs_strict_upper': str(new_upper),
    'endpoint_scalar_mean_actual_error_reduction_factor_strict_lower': '800000',
    'new_uniform_initial_u_rr_H_upper': str(new_uniform_initial_u2),
    'new_uniform_initial_u_rr_H_readable_strict_upper': '6.219e33',
    'old_uniform_initial_u_rr_H_upper': str(old_uniform_initial_u2),
    'velocity_and_its_initial_error_bound_unchanged': True,
    'remaining_cutoff_curvature_jump_locations': ['65/64', '66/64'],
    'rows': out, 'input_sha256': {'initial_report': INITIAL_PIN, 'independent_endpoint': ENDPOINT_PIN},
    'scope': 'explicit change of numerical initial extension and candidate with rigorous C1 perturbation input; genuine improvement of one endpoint initial curvature error, not a global evolved C2 improvement, full differentiated-source bound, Einstein residual, physical-model obstruction or nonlinear spacetime'}, sort_keys=True))
