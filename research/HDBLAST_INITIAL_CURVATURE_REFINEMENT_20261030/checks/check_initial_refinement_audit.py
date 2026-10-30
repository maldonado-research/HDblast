#!/usr/bin/env python3
"""Independent exact-cover/leaf/aggregation audit for unchanged initial data.

No refinement-producer imports. Original authenticated order-five algebra is
reused to check new leaf evaluations; aggregate and convolution arithmetic
are independently reconstructed here.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

CORE_PIN = 'f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782'
INITIAL_PIN = 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207'
REPAIR_PIN = '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002'
checks = []
controls = []


def need(ok, message):
    if not ok:
        raise RuntimeError(message)
    checks.append(message)


def checked_bytes(path, digest, role):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'authenticate ' + role)
    return raw


def read_report(path, digest, role):
    return json.loads(checked_bytes(path, digest, role))


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def norm(caps, weights):
    square = sum((w*z*z for w, z in zip(weights, caps)), Q(0))
    grid = 1 << 112
    k = math.isqrt((square*grid*grid).__ceil__()) + 1
    value = Q(k, grid) if square else Q(0)
    need(value*value >= square, 'positive weighted rational norm')
    return value


def reject(name, fn):
    try:
        fn()
    except (RuntimeError, ValueError, KeyError, IndexError):
        controls.append(name)
    else:
        raise RuntimeError('negative control was accepted: ' + name)


def matrix_q(m):
    return [[Q(v) for v in row] for row in m]


def structural(d, original, labels, weights):
    need(d['status'] == 'PASS_UNCHANGED_QUADRATIC_INITIAL_CURVATURE_REFINEMENT' and
         d['F_G_and_coefficients_unchanged'] is True and
         d['refined_u_rr_H_cap_is_upper_bound'] is True and
         d['actual_error_reduction_claimed'] is False,
         'successful unchanged-source upper-certificate scope')
    need(d['labels'] == labels and list(map(Q, d['gram'])) == weights,
         'same full23 labels and positive H weights')
    need(d['same_candidate_W1_quadratic'] is True and d['candidate_changed'] is False and
         d['target_unchanged'] is True and d['initial_data_unchanged'] is True,
         'unchanged actual target and W1_quadratic candidate')
    need(Q(d['fixed_displacement_kernel_width']) == Q(1, 2**30) and
         Q(d['fixed_velocity_kernel_width']) == Q(1, 2048),
         'unchanged displacement and velocity convolution widths')
    leaves = d['refinement']['leaves']
    need(len(leaves) >= 2, 'nontrivial cutoff subdivision')
    ranges = [list(map(Q, leaf['r'])) for leaf in leaves]
    need(ranges[0][0] == Q(65, 64) and ranges[-1][1] == Q(66, 64),
         'refined cutoff endpoints exact')
    need(all(a < b for a, b in ranges) and
         all(a[1] == b[0] for a, b in zip(ranges, ranges[1:])),
         'exact ordered cutoff leaf cover without gaps or overlaps')
    retained = [cell for cell in d['complete_cover'] if cell['role'] == 'retained']
    expected = [dict(cell, role='retained') for cell in original['continuous']['partition']
                if list(map(Q, cell['r'])) != [Q(65, 64), Q(66, 64)]]
    need(retained == expected and len(retained) == 47,
         'all47 unchanged old cells and authenticated per-cell H caps')
    expected_refined = [{'r': leaf['r'], 'channel_H_caps': leaf['channel_H_caps'], 'role': 'refined'} for leaf in leaves]
    need([cell for cell in d['complete_cover'] if cell['role'] == 'refined'] == expected_refined,
         'complete cover uses exactly audited refined leaves')
    whole = d['complete_cover']
    wr = [list(map(Q, row['r'])) for row in whole]
    need(wr[0][0] == 1 and wr[-1][1] == Q(7, 4) and
         all(a[1] == b[0] for a, b in zip(wr, wr[1:])),
         'complete active-domain cover from1 to7/4')
    return leaves, whole


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--higher-capsule', type=Path)
    ap.add_argument('--refined-report', type=Path, required=True)
    ap.add_argument('--refined-sha256', required=True)
    args = ap.parse_args()
    hi = args.higher_capsule or args.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
    original = read_report(hi/'reference/INITIAL.json', INITIAL_PIN, 'old higher initial report')
    repair = read_report(hi/'reference/REPAIR.json', REPAIR_PIN, 'same quadratic extension report')
    d = read_report(args.refined_report, args.refined_sha256, 'refinement report')
    labels = original['labels']
    weights = list(map(Q, original['gram']))
    need(len(labels) == len(weights) == 23 and min(weights) > 0,
         'complete positive full23 norm')
    leaves, whole = structural(d, original, labels, weights)
    core_path = hi/'checks/initial_second_core.py'
    checked_bytes(core_path, CORE_PIN, 'original order-five core before import')
    spec = importlib.util.spec_from_file_location('independent_initial_refinement_core', core_path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    h = mod.HigherInitial(args.repo_root)
    ell = h.ell
    need(ell == Q(1, 10**6) and list(h.ctx.op.LABELS) == labels and
         h.ctx.weights == weights, 'same registered algebra and derivative scale')

    for k, leaf in enumerate(leaves):
        a, b = map(Q, leaf['r'])
        displacement, velocity = h.paired(h.I(a, b))
        actual = []
        for ch in range(3):
            actual.append([
                [(2*ell**2*row.d[ch].c[2]).mag() for row in displacement],
                [(ell*row.d[ch].c[1]).mag() for row in velocity],
            ])
        claimed = [[[Q(z) for z in jet] for jet in channel]
                   for channel in leaf['paired_channel_component_caps']]
        need(actual == claimed, 'independent original-core component evaluation leaf' + str(k))
        caps_h = matrix_q(leaf['channel_H_caps'])
        for ch in range(3):
            for jet in range(2):
                need(norm(actual[ch][jet], weights) == caps_h[ch][jet],
                     'independent full23 norm leaf' + str(k) + ':' + str(ch) + ':' + str(jet))
    maxima = [[max(Q(cell['channel_H_caps'][ch][jet]) for cell in whole)
               for jet in range(2)] for ch in range(3)]
    need(maxima == matrix_q(d['paired_channel_H_caps']),
         'global channel H maxima over complete cover')
    need(all(max(Q(ch[0]) for ch in leaf['channel_H_caps']) <= Q('1e28') for leaf in leaves),
         'each cutoff leaf meets displacement-channel threshold')

    # Independent integration of k(s)=35/32*(1-s^2)^3 on the half interval.
    polynomial = [(0, Q(35, 32)), (2, -Q(105, 32)),
                  (4, Q(105, 32)), (6, -Q(35, 32))]
    m0 = sum(c/Q(n+1) for n, c in polynomial)
    m1 = sum(c/Q(n+2) for n, c in polynomial)
    need(m0 == Q(1, 2) and m1 == Q(35, 256), 'exact half-kernel mass and first moment')
    du, dv, R = Q(1, 2**30), Q(1, 2048), Q(7, 4)
    need(2*du < Q(1, 64), 'kernel footprint meets at most one isolated curvature join')
    jumps = original['curvature_jumps']
    need([Q(row['r']) for row in jumps] == [Q(65, 64), Q(66, 64), R],
         'retained curvature joins complete')
    need(all(Q(v) == 0 for v in jumps[0]['jump_abs_component_caps']),
         'first inner nominal curvature jump exactly zero')
    inner = [max(Q(row['jump_abs_component_caps'][i]) for row in jumps[:2])
             for i in range(23)]
    smooth_u = [Q(m)*2*m1*du+m0*j
                for m, j in zip(original['piecewise_nominal_M3u'], inner)]
    need(smooth_u == list(map(Q, d['unchanged_curvature_smoothing_component_caps'])),
         'same inner-jump and global third-derivative smoothing charges')
    need(repair['velocity_and_its_initial_error_bound_unchanged'] is True,
         'same rounded-node velocity and endpoint extension')
    smooth_v = [Q(m)*2*m1*dv for m in original['piecewise_nominal_M2v']]
    need(smooth_v == list(map(Q, original['velocity_gradient_convolution_component_caps'])),
         'retained velocity convolution charge')
    # Recheck the P1 secant/nodal charges from authenticated original sources.
    base = h.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results'
    ip = read_report(base/'INITIAL_REPORT.json', mod.INITIAL_PIN, 'original nominal interpolation report')
    profiles = read_report(base/'INITIAL_PROFILES.json', mod.PROFILES_PIN, 'original P1 profile nodes')
    nodes = list(map(Q, profiles['nodes']))
    need(len(nodes)-1 == len(ip['interpolation_cells']) and nodes[-1] == R,
         'same complete P1 node partition')
    p1 = [Q(0)]*23
    for i, cell in enumerate(ip['interpolation_cells']):
        left, right = nodes[i:i+2]
        step = right-left
        need(step > 0 and list(map(Q, cell['r'])) == [left, right], 'exact P1 interpolation cell')
        for j in range(23):
            value = step*Q(cell['second_derivative_component_caps'][2][j])/2
            value += 2*Q(cell['component_error_caps'][2][j])/step
            p1[j] = max(p1[j], value)
    need(p1 == list(map(Q, original['P1_velocity_derivative_component_error_caps'])),
         'independent original P1 secant and nodal-radius charge')
    smooth_u_h, smooth_v_h, p1_h = norm(smooth_u, weights), norm(smooth_v, weights), norm(p1, weights)
    u2 = sum((row[0] for row in maxima), Q(0))+smooth_u_h
    v1 = sum((Q(row[1]) for row in original['continuous']['paired_channel_H_caps']), Q(0))+smooth_v_h+p1_h
    tighter_paired_v_diagnostic = sum((row[1] for row in maxima), Q(0))+smooth_v_h+p1_h
    need(u2 == Q(d['new_uniform_initial_u_rr_H_upper']) and
         v1 == Q(d['unchanged_full_initial_v_r_error_H_upper']),
         'independent complete continuous higher initial-error ledger')
    old_u, old_v = Q(repair['new_uniform_initial_u_rr_H_upper']), Q(original['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
    need(u2 < old_u and v1 == old_v and tighter_paired_v_diagnostic <= v1,
         'uniform displacement upper certificate improves with retained full velocity cap')
    need(u2 < Q('1.185e28') and old_u > 500000*u2,
         'strict readable initial upper-cap improvement without changing actual data')

    # Structural controls are run on copies; each must fail explicitly under -O too.
    bad = copy.deepcopy(d)
    del bad['refinement']['leaves'][len(bad['refinement']['leaves'])//2]
    reject('omitted cutoff leaf', lambda: structural(bad, original, labels, weights))
    bad = copy.deepcopy(d)
    next(cell for cell in bad['complete_cover'] if cell['role'] == 'retained')['channel_H_caps'][0][0] = '1'
    reject('changed retained cell', lambda: structural(bad, original, labels, weights))
    bad = copy.deepcopy(d)
    bad['fixed_displacement_kernel_width'] = str(2*du)
    reject('changed displacement width', lambda: structural(bad, original, labels, weights))
    bad = copy.deepcopy(d)
    bad['gram'][0] = str(2*weights[0])
    reject('wrong Gram weight', lambda: structural(bad, original, labels, weights))
    reject('omitted displacement smoothing charge', lambda: need(u2 == sum(row[0] for row in maxima), 'require positive smoothing charge'))
    reject('omitted P1 velocity-gradient charge', lambda: need(v1 == sum(row[1] for row in maxima)+smooth_v_h, 'require positive P1 charge'))
    nonzero = next((Q(z) for row in leaves for ch in row['paired_channel_component_caps']
                    for jet in ch for z in jet if Q(z)), None)
    need(nonzero is not None, 'nonzero derivative-scale fixture')
    reject('missing u_rr factorial', lambda: need(nonzero == nonzero/2, 'ordinary radial derivative needs2!'))
    h.authenticate()
    checked_bytes(core_path, CORE_PIN, 'original core readback')
    checked_bytes(args.refined_report, args.refined_sha256, 'refinement report readback')
    print(json.dumps(encode({
        'status': 'PASS_INDEPENDENT_UNCHANGED_INITIAL_REFINEMENT_AUDIT',
        'input_sha256': {'refined_report': args.refined_sha256, 'original_initial_report': INITIAL_PIN,
                         'same_quadratic_extension_report': REPAIR_PIN, 'original_order_five_core': CORE_PIN},
        'candidate': 'unchanged W1_quadratic / W1_dagger',
        'labels': labels, 'gram': weights,
        'new_cutoff_leaves_recomputed': len(leaves), 'old_cells_preserved': 47,
        'complete_active_cover_cells': len(whole),
        'paired_channel_H_caps': maxima,
        'displacement_curvature_smoothing_H_upper': smooth_u_h,
        'velocity_gradient_convolution_H_upper': smooth_v_h,
        'P1_velocity_derivative_H_upper': p1_h,
        'actual_minus_smoothed_u_rr_H_upper': u2,
        'actual_minus_smoothed_P1_v_r_H_upper': v1,
        'old_u_rr_H_upper': old_u, 'old_v_r_H_upper': old_v,
        'normalization': 'u_rr=2ell^2*c2; v_r=ell*c1; positive full23 H',
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'producer_code_imported': False,
        'scope': 'Independent source-core leaf recomputation and exact global initial-error ledger audit for same candidate; absolute coordinate data bounds only, no physical relative accuracy or nonlinear/cosmological conclusion.',
    }), sort_keys=True))


if __name__ == '__main__':
    main()
