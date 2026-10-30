#!/usr/bin/env python3
"""Further exact-cover refinement of the unchanged W1_quadratic paired error.

Only interval enclosures change. All nominal/actual functions, smoothing
widths, outer continuation, original source and finite candidate are fixed.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True

PINS = {
    'core': 'f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782',
    'higher_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'initial_repair': '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
    'previous_refinement': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
}
checks = []
controls = []


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)
    checks.append(reason)


def load(path, role):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == PINS[role], 'authenticate '+role)
    return raw


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def tensor_q(x):
    return [[[Q(v) for v in row] for row in channel] for channel in x]


def cover(rows, lo, hi):
    return (bool(rows) and rows[0][0] == lo and rows[-1][1] == hi and
            all(a < b for a, b in rows) and all(a[1] == b[0] for a, b in zip(rows, rows[1:])))


def reject(name, predicate):
    need(not predicate, 'reject '+name)
    controls.append(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--higher-capsule', type=Path)
    ap.add_argument('--previous-refinement', type=Path)
    ap.add_argument('--core-file', type=Path)
    ap.add_argument('--inject-error', choices=('omit-curvature-factor', 'omit-radial-scale', 'changed-width', 'wrong-parent-cover'))
    a = ap.parse_args()
    hi = a.higher_capsule or a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
    core = a.core_file or hi/'checks/initial_second_core.py'
    prev_path = a.previous_refinement or a.repo_root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json'
    raw = load(core, 'core')
    old = json.loads(load(hi/'reference/INITIAL.json', 'higher_initial'))
    repair = json.loads(load(hi/'reference/REPAIR.json', 'initial_repair'))
    previous = json.loads(load(prev_path, 'previous_refinement'))
    spec = importlib.util.spec_from_file_location('paired_refinement_authenticated_original_core', core)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    h = module.HigherInitial(a.repo_root)
    c, I, ell = h.ctx, h.I, h.ell
    labels, weights = list(c.op.LABELS), list(c.weights)
    need(labels == old['labels'] == previous['labels'] and weights == list(map(Q, previous['Gram'])), 'same full23 labels and Gram')
    need(len(labels) == 23 and min(weights) > 0 and weights == list(map(Q, old['gram'])), 'complete positive H')
    need(previous['status'].startswith('PASS') and previous['same_candidate_W1_quadratic'] is True and
         previous['candidate_changed'] is False and previous['target_unchanged'] is True and
         previous['initial_data_unchanged'] is True and previous['F_G_and_coefficients_unchanged'] is True,
         'unchanged successful original-source candidate and target')
    du, dv = Q(previous['fixed_displacement_kernel_width']), Q(previous['fixed_velocity_kernel_width'])
    if a.inject_error == 'changed-width':
        du *= 2
    need(du == Q(1, 2**30) and dv == Q(1, 2048), 'same fixed convolution widths')
    scale = 2*ell**2
    if a.inject_error == 'omit-curvature-factor':
        scale /= 2
    if a.inject_error == 'omit-radial-scale':
        scale = Q(2)
    need(scale == 2*ell**2, 'ordinary scaled radial curvature factor')
    seeds = previous['refinement']['leaves']
    left, right, R = Q(65, 64), Q(66, 64), Q(7, 4)
    seed_cover = [tuple(map(Q, z['r'])) for z in seeds]
    if a.inject_error == 'wrong-parent-cover':
        seed_cover = seed_cover[1:]
    need(len(seeds) == 10 and cover(seed_cover, left, right), 'authenticated exact ten-leaf parent cover')
    old_cells = old['continuous']['partition']
    retained = [z for z in previous['complete_cover'] if z['role'] == 'retained']
    expected = [dict(z, role='retained') for z in old_cells if list(map(Q, z['r'])) != [left, right]]
    need(retained == expected and len(retained) == 47, 'all47 untouched old rows retained exactly')
    threshold, depth_budget, leaf_budget = Q('1e26'), 16, 128
    leaves, splits, evaluations = [], [], []
    # Every child inherits only the authenticated PR103 seed component caps.
    # No unpublished internal-parent enclosure is used in a leaf certificate.
    stack = [(Q(seed['r'][0]), Q(seed['r'][1]), int(seed['depth']), i,
              tensor_q(seed['paired_channel_component_caps'])) for i, seed in reversed(list(enumerate(seeds)))]
    while stack:
        lo, up, depth, seed_index, inherited = stack.pop()
        u, v = h.paired(I(lo, up))
        raw_caps = [[[ (scale*z.d[ch].c[2]).mag() for z in u],
                     [ (ell*z.d[ch].c[1]).mag() for z in v]] for ch in range(3)]
        effective = [[[min(x, y) for x, y in zip(raw_row, old_row)]
                      for raw_row, old_row in zip(raw_channel, old_channel)]
                     for raw_channel, old_channel in zip(raw_caps, inherited)]
        norms = [[c.norm(row) for row in channel] for channel in effective]
        entry = {'r': [lo, up], 'depth': depth, 'seed_index': seed_index,
                 'raw_paired_channel_component_caps': raw_caps,
                 'inherited_paired_channel_component_caps': inherited,
                 'paired_channel_component_caps': effective, 'channel_H_caps': norms}
        evaluations.append({'r': [lo, up], 'depth': depth, 'seed_index': seed_index,
                            'channel_H_caps': norms})
        missed = max(z[0] for z in norms) > threshold
        can_split = depth < depth_budget and len(seeds)+len(splits) < leaf_budget
        if missed and can_split:
            mid = (lo+up)/2
            splits.append({'r': [lo, up], 'depth': depth, 'seed_index': seed_index, 'midpoint': mid,
                           'paired_channel_component_caps': effective, 'channel_H_caps': norms})
            stack.append((mid, up, depth+1, seed_index, inherited))
            stack.append((lo, mid, depth+1, seed_index, inherited))
        else:
            entry['threshold_satisfied'] = not missed
            entry['stop_reason'] = 'threshold' if not missed else ('depth_budget' if depth == depth_budget else 'leaf_budget')
            leaves.append(entry)
    need(cover([tuple(z['r']) for z in leaves], left, right), 'gap-free exact cutoff leaf cover')
    need(len(leaves) == len(seeds)+len(splits) and len(evaluations) == len(leaves)+len(splits), 'complete retained binary refinement forest')
    need(len(leaves) <= leaf_budget and max(z['depth'] for z in leaves) <= depth_budget, 'declared finite evaluation budget')
    need(all(seed_cover[z['seed_index']][0] <= z['r'][0] < z['r'][1] <= seed_cover[z['seed_index']][1] for z in leaves), 'every child lies in its authenticated seed')
    need(all(z['inherited_paired_channel_component_caps'] == tensor_q(seeds[z['seed_index']]['paired_channel_component_caps']) for z in leaves), 'all leaf inherited caps are authenticated PR103 seed caps')
    whole = retained+[{'r': z['r'], 'channel_H_caps': z['channel_H_caps'], 'role': 'refined'} for z in leaves]
    whole.sort(key=lambda z: Q(z['r'][0]))
    need(cover([tuple(map(Q, z['r'])) for z in whole], Q(1), R), 'complete active initial-domain cover')
    channel_caps = [[max(Q(z['channel_H_caps'][ch][jet]) for z in whole) for jet in range(2)] for ch in range(3)]
    paired = sum(z[0] for z in channel_caps)
    previous_paired = Q(previous['paired_u_rr_H_upper'])
    need(0 < paired < previous_paired, 'strictly improved paired comparison upper certificate')
    smoothing = Q(previous['unchanged_curvature_smoothing_H_cap'])
    fallback = paired+smoothing
    need(fallback < Q(previous['refined_uniform_initial_u_rr_H_upper']), 'improves even with old smoothing charge')
    need(Q(old['continuous']['actual_minus_smoothed_P1_v_r_H_upper']) == Q(previous['unchanged_full_initial_v_r_error_H_upper']), 'full velocity-gradient error cap retained')
    ranges = [tuple(z['r']) for z in leaves]
    reject('omitted cutoff leaf', cover(ranges[1:], left, right))
    reject('overlapping cutoff leaf', cover([ranges[0]]+ranges, left, right))
    reject('reversed cutoff leaf', cover([(ranges[0][1], ranges[0][0])]+ranges[1:], left, right))
    reject('maximum replaces telescoping sum', paired == max(z[0] for z in channel_caps))
    reject('paired upper replaces complete smoothed error upper', paired == fallback)
    reject('unresolved old cutoff passes new threshold', max(Q(z[0]) for z in seeds[-1]['channel_H_caps']) <= threshold)
    h.authenticate()
    load(core, 'core')
    load(hi/'reference/INITIAL.json', 'higher_initial')
    load(hi/'reference/REPAIR.json', 'initial_repair')
    load(prev_path, 'previous_refinement')
    report = {
        'status': 'PASS_UNCHANGED_QUADRATIC_PAIRED_CUTOFF_REFINEMENT',
        'candidate': previous['candidate'], 'same_candidate_W1_quadratic': True, 'candidate_changed': False,
        'target_unchanged': True, 'initial_data_unchanged': True, 'F_G_and_coefficients_unchanged': True,
        'labels': labels, 'Gram': weights, 'gram': weights, 'ell': ell,
        'fixed_displacement_kernel_width': du, 'fixed_velocity_kernel_width': dv,
        'bound_quantity': 'sup_H_norm_paired_actual_minus_nominal_initial_u_rr',
        'paired_u_rr_H_cap_is_upper_bound': True,
        'input_sha256': PINS,
        'refinement': {'replaced_old_cell': [left, right], 'retained_old_cells': 47,
                       'previous_cutoff_leaves': 10, 'leaf_threshold_H': threshold,
                       'maximum_depth': depth_budget, 'leaf_budget': leaf_budget,
                       'number_of_leaves': len(leaves), 'number_of_splits': len(splits),
                       'number_of_evaluations': len(evaluations), 'threshold_missed_leaves': sum(not z['threshold_satisfied'] for z in leaves),
                       'leaves': leaves, 'splits': splits, 'evaluations': evaluations},
        'complete_cover': whole, 'paired_channel_H_caps': channel_caps,
        'paired_u_rr_H_upper': paired, 'previous_paired_u_rr_H_upper': previous_paired,
        'paired_upper_certificate_reduction_factor': previous_paired/paired,
        'unchanged_full_initial_v_r_error_H_upper': Q(previous['unchanged_full_initial_v_r_error_H_upper']),
        'fallback_old_smoothing_H_cap': smoothing,
        'fallback_with_old_smoothing_uniform_initial_u_rr_H_upper': fallback,
        'new_nominal_M3_combination_performed_by_this_report': False,
        'actual_error_reduction_claimed': False,
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Tighter paired interval certificate for the unchanged auxiliary candidate. A complete smoothed initial error requires an independently valid convolution charge. No new candidate, lower true error, physical/covariant accuracy, nonlinear or cosmological result.'}
    print(json.dumps(encode(report), sort_keys=True))


if __name__ == '__main__':
    main()
