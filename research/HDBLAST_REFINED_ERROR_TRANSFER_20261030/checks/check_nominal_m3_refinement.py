#!/usr/bin/env python3
"""Refine explicit nominal u_rrr bounds, preserving the original candidate.

Only intervals used to prove a piecewise derivative upper bound change.
No actual sixth derivative or derivative of an unknown remainder is used.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

PINS = {
    'core': 'f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782',
    'nominal_initial': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'higher_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'repair': '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
    'previous_refinement': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
}


def enc(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [enc(v) for v in x]
    return x


def cover(rows, left, right):
    return bool(rows) and rows[0][0] == left and rows[-1][1] == right and all(
        a < b for a, b in rows) and all(a[1] == b[0] for a, b in zip(rows, rows[1:]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--higher-capsule', type=Path)
    ap.add_argument('--core-file', type=Path)
    ap.add_argument('--inject-error', choices=('missing-cover-cell', 'wrong-branch',
        'missing-factorial', 'missing-ell-cube', 'omit-inner-jump', 'parent-maximum'))
    a = ap.parse_args()
    groups, controls = [], []

    def need(ok, label):
        if not ok:
            raise RuntimeError(label)
        groups.append(label)

    def load(path, role):
        raw = path.read_bytes()
        need(hashlib.sha256(raw).hexdigest() == PINS[role], 'authenticate '+role)
        return raw

    root = a.repo_root.resolve()
    hi = a.higher_capsule or root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
    core = a.core_file or hi/'checks/initial_second_core.py'
    core_raw = load(core, 'core')
    nominal_path = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_REPORT.json'
    nominal = json.loads(load(nominal_path, 'nominal_initial'))
    higher = json.loads(load(hi/'reference/INITIAL.json', 'higher_initial'))
    repair = json.loads(load(hi/'reference/REPAIR.json', 'repair'))
    previous_path = root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json'
    previous = json.loads(load(previous_path, 'previous_refinement'))
    spec = importlib.util.spec_from_file_location('authenticated_nominal_m3_core', core)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    h = mod.HigherInitial(root)
    c, I, ell = h.ctx, h.I, h.ell
    labels, weights = list(c.op.LABELS), c.weights
    need(labels == nominal['labels'] == higher['labels'] == repair['labels'] == previous['labels'] and
         weights == list(map(Q, nominal['gram'])) == list(map(Q, previous['Gram'])), 'same full23 positive Gram')
    need(len(labels) == len(weights) == 23 and min(weights) > 0, 'complete positive H norm')
    need(labels[12:15] == ['w:0', 'w:2', 'w:4'] and ell == Q(1, 10**6), 'registered scalar order and coordinate scale')
    need(previous['same_candidate_W1_quadratic'] is True and previous['candidate_changed'] is False,
         'same unchanged W1_quadratic candidate')
    rows = nominal['interpolation_cells']
    old_cover = [tuple(map(Q, row['r'])) for row in rows]
    if a.inject_error == 'missing-cover-cell':
        old_cover = old_cover[1:]
    need(len(rows) == 189 and cover(old_cover, Q(0), Q(7, 4)), 'complete nominal parent cover')
    old = list(map(Q, nominal['nominal_second_derivative_component_caps'][1]))
    old_cells = [list(map(Q, row['second_derivative_component_caps'][1])) for row in rows]
    need([max(row[i] for row in old_cells) for i in range(23)] == old ==
         list(map(Q, higher['piecewise_nominal_M3u'])) == list(map(Q, previous['unchanged_nominal_M3u_component_caps'])),
         'old nominal M3 maxima reconstructed exactly')
    denominator = 10
    selected = [i for i, row in enumerate(old_cells) if any(row[j] > old[j]/denominator for j in (12, 13, 14))]
    need(len(selected) == 35 and old_cover[selected[0]][0] == Q(273, 256) and
         old_cover[selected[-1]][1] == Q(77, 64), 'declared tenth-threshold selected work list')
    need(all(old_cover[i][0] >= Q(66, 64) for i in selected), 'all selected parents are plain branch')
    branch = 'cutoff' if a.inject_error == 'wrong-branch' else 'plain'
    need(branch == 'plain', 'forced plain nominal branch')
    factor = 1 if a.inject_error == 'missing-factorial' else math.factorial(3)
    scale = Q(1) if a.inject_error == 'missing-ell-cube' else ell**3
    need(factor == 6, 'ordinary third derivative factorial')
    need(scale == Q(1, 10**18), 'scaled radial third derivative ell cube')
    choose = max if a.inject_error == 'parent-maximum' else min
    need(choose(Q(7), Q(11)) == Q(7), 'parent cap intersection uses minimum')
    delta = Q(higher['fixed_displacement_kernel_width'])
    jumps = higher['curvature_jumps']
    need(delta == Q(1, 2**30) and [Q(row['r']) for row in jumps] == [Q(65, 64), Q(66, 64), Q(7, 4)],
         'same smoothing width and complete old joins')
    need(all(Q(v) == 0 for v in jumps[0]['jump_abs_component_caps']), 'zero first inner nominal curvature jump')
    inner = [max(Q(row['jump_abs_component_caps'][i]) for row in jumps[:2]) for i in range(23)]
    if a.inject_error == 'omit-inner-jump':
        inner = [Q(0)]*23
    need(inner == list(map(Q, previous['unchanged_inner_jump_component_caps'])), 'retained inner jump allowance')
    need(repair['velocity_and_its_initial_error_bound_unchanged'] is True, 'same velocity and quadratic outer continuation')

    leaves, splits = [], []
    maximum_depth = 4
    evaluations = 0
    for index in selected:
        left, right = old_cover[index]
        stack = [(left, right, 0)]
        parent = old_cells[index]
        while stack:
            left, right, depth = stack.pop()
            u, v = h.nominal(I(left, right), branch)
            need(len(u) == 23 and all(z.n >= 3 for z in u), 'nominal leaf has all third output jets')
            raw = [(factor*scale*z.c[3]).mag() for z in u]
            cap = [min(x, y) for x, y in zip(raw, parent)]
            passed = all(cap[j] <= old[j]/denominator for j in (12, 13, 14))
            evaluations += 1
            if not passed and depth < maximum_depth:
                mid = (left+right)/2
                splits.append({'parent_index': index, 'r': [left, right], 'midpoint': mid, 'depth': depth})
                stack.append((mid, right, depth+1))
                stack.append((left, mid, depth+1))
            else:
                leaves.append({'parent_index': index, 'r': [left, right], 'depth': depth,
                    'raw_component_caps': raw, 'effective_component_caps': cap,
                    'scalar_tenth_target_met': passed})
    for index in selected:
        ranges = [tuple(row['r']) for row in leaves if row['parent_index'] == index]
        need(cover(ranges, *old_cover[index]), 'exact selected parent leaf cover '+str(index))
    retained = [{'parent_index': i, 'r': old_cover[i], 'component_caps': old_cells[i], 'role': 'retained'}
                for i in range(len(rows)) if i not in selected]
    complete = retained + [{'parent_index': row['parent_index'], 'r': row['r'],
                'component_caps': row['effective_component_caps'], 'role': 'refined'} for row in leaves]
    complete.sort(key=lambda row: row['r'][0])
    need(cover([tuple(row['r']) for row in complete], Q(0), Q(7, 4)), 'complete refined nominal domain cover')
    need(len(retained) == 154 and evaluations == len(leaves)+len(splits) and len(leaves) == len(splits)+len(selected),
         'finite binary trees and unchanged old cells')
    new = [max(row['component_caps'][i] for row in complete) for i in range(23)]
    need(all(Q(0) <= x <= y for x, y in zip(new, old)), 'each refined nominal derivative bound cannot increase')
    need(any(new[i] < old[i] for i in range(23)), 'strict nominal M3 improvement')
    mu = Q(35, 128)*delta
    smooth = [mu*z+j/2 for z, j in zip(new, inner)]
    old_smooth = [mu*z+j/2 for z, j in zip(old, inner)]
    old_H, new_H = c.norm(old_smooth), c.norm(smooth)
    need(old_H == Q(previous['unchanged_curvature_smoothing_H_cap']), 'old complete smoothing ledger reproduced exactly')
    need(new_H < old_H, 'strict complete smoothing upper bound improvement')
    misses = [i for i, row in enumerate(leaves) if not row['scalar_tenth_target_met']]
    need(all(leaves[i]['depth'] == maximum_depth for i in misses), 'threshold misses retained honestly at finite depth limit')
    # Explicit invalid cover and scaling fixtures: optimization cannot remove decisions.
    fixtures = [('omitted initial parent', not cover(old_cover[1:], Q(0), Q(7, 4))),
                ('duplicated first parent', not cover([old_cover[0]]+old_cover, Q(0), Q(7, 4))),
                ('ordinary third derivative is not normalized coefficient', math.factorial(3) != 1),
                ('radial scale cannot be omitted', ell**3 != 1),
                ('sum retains inner jump', smooth != [mu*z for z in new]),
                ('parent maximum is not a valid tightening', max(Q(7), Q(11)) != min(Q(7), Q(11)))]
    for name, okay in fixtures:
        need(okay, 'reject '+name)
        controls.append(name)
    h.authenticate()
    load(core, 'core'); load(nominal_path, 'nominal_initial')
    load(hi/'reference/INITIAL.json', 'higher_initial'); load(hi/'reference/REPAIR.json', 'repair')
    load(previous_path, 'previous_refinement')
    need(hashlib.sha256(core.read_bytes()).hexdigest() == hashlib.sha256(core_raw).hexdigest(), 'unchanged local core readback')
    report = {'status': 'PASS_NOMINAL_THIRD_DERIVATIVE_REFINEMENT', 'candidate': previous['candidate'],
        'same_candidate_W1_quadratic': True, 'candidate_changed': False,
        'underlying_actual_target_unchanged': True, 'F_G_and_coefficients_unchanged': True,
        'labels': labels, 'Gram': weights, 'ell': ell,
        'bound_quantity': 'piecewise_nominal_u_rrr_component_upper', 'nominal_M3_upper': True,
        'actual_order_six_bound_claimed': False, 'branch': branch,
        'ordinary_scaled_third_derivative': '6 ell^3 c[3]',
        'original_nominal_M3u_component_caps': old, 'refined_nominal_M3u_component_caps': new,
        'original_nominal_M3u_H_cap': c.norm(old), 'refined_nominal_M3u_H_cap': c.norm(new),
        'fixed_displacement_kernel_width': delta, 'fixed_velocity_kernel_width': Q(1, 2048),
        'unchanged_inner_jump_component_caps': inner,
        'original_curvature_smoothing_component_caps': old_smooth,
        'refined_curvature_smoothing_component_caps': smooth,
        'original_curvature_smoothing_H_cap': old_H, 'refined_curvature_smoothing_H_cap': new_H,
        'smoothing_upper_certificate_reduction_factor': old_H/new_H,
        'selected_parent_indices': selected, 'retained_parent_count': len(retained),
        'refinement': {'scalar_target_fraction_of_old_global': Q(1, denominator),
            'maximum_depth': maximum_depth, 'number_of_leaves': len(leaves), 'number_of_splits': len(splits),
            'number_of_evaluations': evaluations, 'threshold_miss_leaf_indices': misses,
            'leaves': leaves, 'splits': splits}, 'complete_cover': complete,
        'input_sha256': PINS, 'named_checks': groups, 'wrong_formula_controls_rejected': controls,
        'scope': 'Piecewise third derivatives of explicit nominal initial displacement and its unchanged convolution allowance only. No actual sixth derivative, new candidate, true error reduction, physical accuracy, nonlinear family or cosmological result.'}
    print(json.dumps(enc(report), sort_keys=True))


if __name__ == '__main__':
    main()
