#!/usr/bin/env python3
"""Independent core replay of nominal-M3 and paired-cutoff refinements.

No new refinement producer is imported. All nominal and paired leaf values
are recomputed from the authenticated retained source core. The two complete
covers and the unchanged convolution law give one combined initial error cap.
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

PINS = {
    'core': 'f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782',
    'nominal_initial': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'initial_profiles': '555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1',
    'higher_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'repair': '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
    'previous_refinement': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
}
checks = []
controls = []


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)
    checks.append(reason)


def bytes_checked(path, digest, role):
    data = Path(path).read_bytes()
    need(hashlib.sha256(data).hexdigest() == digest, 'authenticate '+role)
    return data


def read(path, digest, role):
    return json.loads(bytes_checked(path, digest, role))


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def qvector(row, n=23):
    need(len(row) == n, 'complete component vector')
    out = list(map(Q, row))
    need(min(out) >= 0, 'nonnegative upper component vector')
    return out


def qtensor(tensor):
    need(len(tensor) == 3 and all(len(channel) == 2 for channel in tensor), 'three paired channels and two jet slots')
    return [[qvector(row) for row in channel] for channel in tensor]


def qmatrix(row):
    need(len(row) == 3 and all(len(z) == 2 for z in row), 'three-channel H matrix')
    return [list(map(Q, channel)) for channel in row]


def norm(caps, weights):
    square = sum((w*v*v for w, v in zip(weights, caps)), Q(0))
    need(len(caps) == len(weights) == 23 and min(weights) > 0 and min(caps) >= 0,
         'complete positive weighted H input')
    if not square:
        return Q(0)
    grid = 1 << 112
    q = square*grid*grid
    ceiling = (q.numerator+q.denominator-1)//q.denominator
    answer = Q(math.isqrt(ceiling)+1, grid)
    need(answer*answer >= square, 'outward H square-root inequality')
    # Independently check closeness to the least dyadic outward square root.
    k = math.isqrt(q.numerator//q.denominator)
    if Q(k, grid)**2 < square:
        k += 1
    need(Q(k, grid) <= answer <= Q(k+1, grid), 'retained H rounding protocol within one outward grid unit')
    return answer


def exact_cover(ranges, lower, upper):
    need(bool(ranges) and ranges[0][0] == lower and ranges[-1][1] == upper,
         'cover exact outer endpoints')
    need(all(a < b for a, b in ranges) and
         all(a[1] == b[0] for a, b in zip(ranges, ranges[1:])),
         'ordered gap-free nonoverlapping cover')


def validate_tree(leaves, splits, roots, index_key):
    def key(row):
        a, b = map(Q, row['r'])
        return int(row[index_key]), a, b, int(row['depth'])
    leaf_map = {key(row): row for row in leaves}
    split_map = {key(row): row for row in splits}
    need(len(leaf_map) == len(leaves) and len(split_map) == len(splits) and
         not set(leaf_map).intersection(split_map), 'unique disjoint tree nodes')
    seen = set()
    stack = list(roots)
    while stack:
        item = stack.pop()
        need(item not in seen, 'tree node reached once')
        seen.add(item)
        index, lower, upper, depth = item
        need(lower < upper and depth >= 0, 'positive tree interval and depth')
        if item in split_map:
            midpoint = Q(split_map[item]['midpoint'])
            need(midpoint == (lower+upper)/2, 'exact binary midpoint')
            stack.extend(((index, lower, midpoint, depth+1), (index, midpoint, upper, depth+1)))
        else:
            need(item in leaf_map, 'every terminal tree node is a recorded leaf')
    need(seen == set(leaf_map).union(split_map), 'all declared tree nodes reached from authenticated roots')
    need(len(leaves) == len(splits)+len(roots), 'complete finite binary forest count')


def common(d, labels, weights):
    need(d['status'].startswith('PASS') and d['same_candidate_W1_quadratic'] is True and
         d['candidate_changed'] is False and d['F_G_and_coefficients_unchanged'] is True,
         'successful unchanged candidate and original sources')
    need(d['labels'] == labels and list(map(Q, d['Gram'])) == weights,
         'same complete labels and positive Gram')
    need(Q(d['ell']) == Q(1, 10**6) and
         Q(d['fixed_displacement_kernel_width']) == Q(1, 2**30) and
         Q(d['fixed_velocity_kernel_width']) == Q(1, 2048), 'same derivative scale and both smoothing widths')


def structural_nominal(d, nominal, labels, weights):
    common(d, labels, weights)
    need(d['status'] == 'PASS_NOMINAL_THIRD_DERIVATIVE_REFINEMENT' and
         d['underlying_actual_target_unchanged'] is True and
         d['nominal_M3_upper'] is True and d['actual_order_six_bound_claimed'] is False and
         d['bound_quantity'] == 'piecewise_nominal_u_rrr_component_upper' and
         d['branch'] == 'plain' and d['ordinary_scaled_third_derivative'] == '6 ell^3 c[3]',
         'explicit nominal third derivative only with forced plain branch')
    old_rows = nominal['interpolation_cells']
    ranges = [tuple(map(Q, row['r'])) for row in old_rows]
    need(len(ranges) == 189, 'all189 authenticated nominal parents')
    exact_cover(ranges, Q(0), Q(7, 4))
    parent_caps = [qvector(row['second_derivative_component_caps'][1]) for row in old_rows]
    old = [max(row[j] for row in parent_caps) for j in range(23)]
    need(old == list(map(Q, nominal['nominal_second_derivative_component_caps'][1])) ==
         list(map(Q, d['original_nominal_M3u_component_caps'])), 'old nominal M3 vector reproduced from all189 parents')
    selected = [i for i, row in enumerate(parent_caps) if any(row[j] > old[j]/10 for j in (12, 13, 14))]
    need(d['selected_parent_indices'] == selected and len(selected) == 35,
         'exact tenth-threshold selected35 work list')
    need(all(ranges[i][0] >= Q(66, 64) for i in selected), 'selected parents remain wholly on plain branch')
    leaves = d['refinement']['leaves']
    maximum_depth = int(d['refinement']['maximum_depth'])
    need(maximum_depth == 4 and Q(d['refinement']['scalar_target_fraction_of_old_global']) == Q(1, 10),
         'declared nominal finite-depth and threshold budget')
    for leaf in leaves:
        i = int(leaf['parent_index'])
        need(i in selected, 'nominal leaf has selected authenticated parent')
        a, b = map(Q, leaf['r'])
        lower, upper = ranges[i]
        depth = int(leaf['depth'])
        need(0 <= depth <= maximum_depth and lower <= a < b <= upper,
             'nominal child within authenticated parent and finite depth')
        width = (upper-lower)/2**depth
        need(b-a == width and ((a-lower)/width).denominator == 1, 'nominal leaf is exact dyadic parent subdivision')
        raw = qvector(leaf['raw_component_caps'])
        effective = qvector(leaf['effective_component_caps'])
        need(effective == [min(x, y) for x, y in zip(raw, parent_caps[i])], 'nominal child cap uses authenticated parent minimum')
        target = all(effective[j] <= old[j]/10 for j in (12, 13, 14))
        need(leaf['scalar_tenth_target_met'] == target and (target or depth == maximum_depth),
             'nominal threshold misses retained honestly at finite budget')
    for i in selected:
        exact_cover([tuple(map(Q, row['r'])) for row in leaves if int(row['parent_index']) == i], *ranges[i])
    validate_tree(leaves, d['refinement']['splits'], [(i, *ranges[i], 0) for i in selected], 'parent_index')
    statistics = d['refinement']
    need(statistics['number_of_leaves'] == len(leaves) and
         statistics['number_of_splits'] == len(statistics['splits']) and
         statistics['number_of_evaluations'] == len(leaves)+len(statistics['splits']) and
         statistics['threshold_miss_leaf_indices'] == [i for i, row in enumerate(leaves) if not row['scalar_tenth_target_met']],
         'honest nominal evaluation and threshold-miss counts')
    expected = [{'parent_index': i, 'r': list(old_rows[i]['r']), 'component_caps': list(old_rows[i]['second_derivative_component_caps'][1]),
                 'role': 'retained'} for i in range(189) if i not in selected]
    expected += [{'parent_index': row['parent_index'], 'r': row['r'], 'component_caps': row['effective_component_caps'],
                  'role': 'refined'} for row in leaves]
    expected.sort(key=lambda row: Q(row['r'][0]))
    need(d['complete_cover'] == expected and d['retained_parent_count'] == 154,
         'nominal complete cover preserves exactly154 unselected old cells')
    exact_cover([tuple(map(Q, row['r'])) for row in expected], Q(0), Q(7, 4))
    return old, parent_caps, leaves, expected


def structural_paired(d, previous, higher, labels, weights):
    common(d, labels, weights)
    need(d['status'] == 'PASS_UNCHANGED_QUADRATIC_PAIRED_CUTOFF_REFINEMENT' and
         d['target_unchanged'] is True and d['initial_data_unchanged'] is True and
         d['paired_u_rr_H_cap_is_upper_bound'] is True and d['actual_error_reduction_claimed'] is False and
         d['bound_quantity'] == 'sup_H_norm_paired_actual_minus_nominal_initial_u_rr',
         'paired result is unchanged-data upper certificate')
    seeds = previous['refinement']['leaves']
    need(len(seeds) == 10, 'ten authenticated PR103 paired seeds')
    exact_cover([tuple(map(Q, row['r'])) for row in seeds], Q(65, 64), Q(66, 64))
    leaves = d['refinement']['leaves']
    maximum_depth, budget = int(d['refinement']['maximum_depth']), int(d['refinement']['leaf_budget'])
    need(maximum_depth == 16 and budget == 128 and len(leaves) <= budget,
         'declared finite paired depth and leaf budget')
    need(Q(d['refinement']['leaf_threshold_H']) == Q('1e26'), 'paired refinement declared threshold')
    for row in leaves:
        i = int(row['seed_index'])
        need(0 <= i < len(seeds), 'paired leaf has authenticated seed')
        lower, upper = map(Q, seeds[i]['r'])
        a, b = map(Q, row['r'])
        depth = int(row['depth'])
        base_depth = int(seeds[i]['depth'])
        need(base_depth <= depth <= maximum_depth and lower <= a < b <= upper,
             'paired leaf lies in authenticated seed at finite depth')
        width = (upper-lower)/2**(depth-base_depth)
        need(b-a == width and ((a-lower)/width).denominator == 1, 'paired leaf exact dyadic seed subdivision')
        inherited = qtensor(row['inherited_paired_channel_component_caps'])
        seed = qtensor(seeds[i]['paired_channel_component_caps'])
        need(inherited == seed, 'paired inheritance exactly pinned PR103 seed cap')
        raw = qtensor(row['raw_paired_channel_component_caps'])
        effective = qtensor(row['paired_channel_component_caps'])
        correct = [[[min(x, y) for x, y in zip(fresh, old)] for fresh, old in zip(fc, oc)]
                   for fc, oc in zip(raw, seed)]
        need(effective == correct, 'paired effective cap is fresh-child authenticated-seed minimum')
    for i, seed in enumerate(seeds):
        exact_cover([tuple(map(Q, row['r'])) for row in leaves if int(row['seed_index']) == i], *map(Q, seed['r']))
    roots = [(i, *map(Q, row['r']), int(row['depth'])) for i, row in enumerate(seeds)]
    validate_tree(leaves, d['refinement']['splits'], roots, 'seed_index')
    statistics = d['refinement']
    need(statistics['number_of_leaves'] == len(leaves) and
         statistics['number_of_splits'] == len(statistics['splits']) and
         statistics['number_of_evaluations'] == len(leaves)+len(statistics['splits']) and
         statistics['threshold_missed_leaves'] == sum(not row['threshold_satisfied'] for row in leaves) and
         statistics['retained_old_cells'] == 47 and statistics['previous_cutoff_leaves'] == 10,
         'honest paired forest and threshold-miss counts')
    retained = [dict(row, role='retained') for row in higher['continuous']['partition']
                if list(map(Q, row['r'])) != [Q(65, 64), Q(66, 64)]]
    expected = retained+[{'r': row['r'], 'channel_H_caps': row['channel_H_caps'], 'role': 'refined'} for row in leaves]
    expected.sort(key=lambda row: Q(row['r'][0]))
    need(len(retained) == 47 and d['complete_cover'] == expected, 'all47 old paired cells retained in complete new cover')
    exact_cover([tuple(map(Q, row['r'])) for row in expected], Q(1), Q(7, 4))
    return seeds, leaves, expected


def reject(name, callback):
    begin = len(checks)
    try:
        callback()
    except (RuntimeError, ValueError, KeyError, IndexError, TypeError):
        del checks[begin:]
        controls.append(name)
        checks.append('reject '+name)
    else:
        raise RuntimeError('negative control accepted: '+name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--higher-capsule', type=Path)
    ap.add_argument('--nominal-report', type=Path, required=True)
    ap.add_argument('--nominal-sha256', required=True)
    ap.add_argument('--paired-report', type=Path, required=True)
    ap.add_argument('--paired-sha256', required=True)
    ap.add_argument('--combined-report', type=Path)
    ap.add_argument('--combined-sha256')
    a = ap.parse_args()
    root = a.repo_root.resolve()
    hi = a.higher_capsule or root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
    core = hi/'checks/initial_second_core.py'
    nominal_path = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_REPORT.json'
    profiles_path = nominal_path.with_name('INITIAL_PROFILES.json')
    previous_path = root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/INITIAL.json'
    raw_core = bytes_checked(core, PINS['core'], 'retained core before import')
    nominal = read(nominal_path, PINS['nominal_initial'], 'old nominal initial report')
    profiles = read(profiles_path, PINS['initial_profiles'], 'unchanged nominal profiles and parity')
    higher = read(hi/'reference/INITIAL.json', PINS['higher_initial'], 'old higher initial report')
    repair = read(hi/'reference/REPAIR.json', PINS['repair'], 'same quadratic extension')
    previous = read(previous_path, PINS['previous_refinement'], 'PR103 initial refinement')
    new_nominal = read(a.nominal_report, a.nominal_sha256, 'new nominal M3 report')
    new_paired = read(a.paired_report, a.paired_sha256, 'new paired refinement report')
    labels, weights = nominal['labels'], list(map(Q, nominal['gram']))
    need(len(labels) == len(weights) == 23 and min(weights) > 0 and labels[12:15] == ['w:0', 'w:2', 'w:4'],
         'complete registered full23 angular H norm')
    need(labels == higher['labels'] == repair['labels'] == previous['labels'] and
         weights == list(map(Q, higher['gram'])) == list(map(Q, repair['gram'])) == list(map(Q, previous['Gram'])),
         'same all23 norm in all retained sources')
    parity = [-1 if label.startswith(('L:', 'v:')) else 1 for label in labels]
    need(profiles['labels'] == labels and list(map(Q, profiles['gram'])) == weights and
         profiles['parity'] == {'u': parity, 'du': [-p for p in parity], 'v': parity},
         'same original reflected profile parity without changing initial functions')
    old_m3, parent_caps, nominal_leaves, nominal_cover = structural_nominal(new_nominal, nominal, labels, weights)
    seeds, paired_leaves, paired_cover = structural_paired(new_paired, previous, higher, labels, weights)
    spec = importlib.util.spec_from_file_location('independent_combined_initial_retained_core', core)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    h = module.HigherInitial(root)
    ell = h.ell
    need(ell == Q(1, 10**6) and list(h.ctx.op.LABELS) == labels and h.ctx.weights == weights,
         'same authenticated core labels Gram and coordinate scale')

    # Explicit-polynomial ordinary derivative fixture, independent of interval
    # producer code: D_r^3(r^3)=6 under the normalized physical-z jet convention.
    cubic = h.J([h.I(0), 1/ell, 0, 0, 0, 0, 0])**3
    ordinary_cubic = 6*ell**3*cubic.c[3]
    need(ordinary_cubic.lo <= 6 <= ordinary_cubic.hi and
         ordinary_cubic.hi-ordinary_cubic.lo < Q('1e-12') and
         not (ordinary_cubic.lo <= 1 <= ordinary_cubic.hi),
         'ordinary third derivative enclosure requires3! and ell^3')
    nominal_raw = []
    for j, leaf in enumerate(nominal_leaves):
        lower, upper = map(Q, leaf['r'])
        displacement, _ = h.nominal(h.I(lower, upper), 'plain')
        need(len(displacement) == 23 and all(row.n >= 3 for row in displacement), 'all explicit nominal output third jets '+str(j))
        actual = [(6*ell**3*row.c[3]).mag() for row in displacement]
        need(actual == qvector(leaf['raw_component_caps']), 'independent nominal leaf component evaluation '+str(j))
        nominal_raw.append(actual)
        expected = [min(x, y) for x, y in zip(actual, parent_caps[int(leaf['parent_index'])])]
        need(expected == qvector(leaf['effective_component_caps']), 'independent nominal parent cap intersection '+str(j))
    new_m3 = [max(Q(row['component_caps'][j]) for row in nominal_cover) for j in range(23)]
    need(new_m3 == list(map(Q, new_nominal['refined_nominal_M3u_component_caps'])) and
         all(0 <= x <= y for x, y in zip(new_m3, old_m3)), 'independent complete nominal M3 maxima')
    need(norm(new_m3, weights) == Q(new_nominal['refined_nominal_M3u_H_cap']) and
         norm(old_m3, weights) == Q(new_nominal['original_nominal_M3u_H_cap']), 'nominal positive H packing')

    paired_raw = []
    for j, leaf in enumerate(paired_leaves):
        lower, upper = map(Q, leaf['r'])
        displacement, velocity = h.paired(h.I(lower, upper))
        actual = [[[(2*ell**2*row.d[channel].c[2]).mag() for row in displacement],
                   [(ell*row.d[channel].c[1]).mag() for row in velocity]] for channel in range(3)]
        need(actual == qtensor(leaf['raw_paired_channel_component_caps']), 'independent paired leaf component evaluation '+str(j))
        paired_raw.append(actual)
        seed = qtensor(seeds[int(leaf['seed_index'])]['paired_channel_component_caps'])
        effective = [[[min(x, y) for x, y in zip(fr, old)] for fr, old in zip(fc, oc)]
                     for fc, oc in zip(actual, seed)]
        need(effective == qtensor(leaf['paired_channel_component_caps']), 'independent paired leaf inherited-cap intersection '+str(j))
        hc = [[norm(row, weights) for row in channel] for channel in effective]
        need(hc == qmatrix(leaf['channel_H_caps']), 'independent paired leaf H caps '+str(j))
        met = max(row[0] for row in hc) <= Q('1e26')
        need(leaf['threshold_satisfied'] == met, 'honest paired leaf threshold metadata '+str(j))
    channel_maxima = [[max(Q(row['channel_H_caps'][channel][jet]) for row in paired_cover)
                       for jet in range(2)] for channel in range(3)]
    paired_cap = sum((row[0] for row in channel_maxima), Q(0))
    need(channel_maxima == qmatrix(new_paired['paired_channel_H_caps']) and
         paired_cap == Q(new_paired['paired_u_rr_H_upper']), 'independent full paired comparison over complete active cover')

    # Integrate the unchanged normalized positive kernel directly on [0,1].
    coefficients = [Q(35, 32), 0, -Q(105, 32), 0, Q(105, 32), 0, -Q(35, 32)]
    half_mass = sum((value/Q(j+1) for j, value in enumerate(coefficients)), Q(0))
    half_first = sum((value/Q(j+2) for j, value in enumerate(coefficients)), Q(0))
    need(half_mass == Q(1, 2) and half_first == Q(35, 256), 'independent exact half kernel moments')
    delta = Q(1, 2**30)
    need(2*delta < Q(1, 64), 'each kernel footprint meets at most one isolated inner join')
    jumps = higher['curvature_jumps']
    need([Q(row['r']) for row in jumps] == [Q(65, 64), Q(66, 64), Q(7, 4)], 'retained complete curvature joins')
    need(all(Q(z) == 0 for z in jumps[0]['jump_abs_component_caps']), 'exact zero first inner nominal curvature jump')
    inner_jump = [max(Q(row['jump_abs_component_caps'][j]) for row in jumps[:2]) for j in range(23)]
    need(inner_jump == list(map(Q, previous['unchanged_inner_jump_component_caps'])) ==
         list(map(Q, new_nominal['unchanged_inner_jump_component_caps'])), 'all second inner jump allowances retained')
    need(repair['velocity_and_its_initial_error_bound_unchanged'] is True and
         repair['definition'].startswith('u_ext_q(R+h)=u_nom(R)+u_nom_r(R) h+(1/2)u_nom_rr(R-) h^2'),
         'same curvature-matched quadratic exterior and P1 velocity')
    old_smooth = [2*half_first*delta*z+half_mass*j for z, j in zip(old_m3, inner_jump)]
    new_smooth = [2*half_first*delta*z+half_mass*j for z, j in zip(new_m3, inner_jump)]
    old_smooth_h, new_smooth_h = norm(old_smooth, weights), norm(new_smooth, weights)
    need(old_smooth_h == Q(previous['unchanged_curvature_smoothing_H_cap']) ==
         Q(new_nominal['original_curvature_smoothing_H_cap']), 'old complete convolution charge reproduced exactly')
    need(new_smooth == list(map(Q, new_nominal['refined_curvature_smoothing_component_caps'])) and
         new_smooth_h == Q(new_nominal['refined_curvature_smoothing_H_cap']), 'independent new complete convolution charge')
    combined = paired_cap+new_smooth_h
    old_combined = Q(previous['refined_uniform_initial_u_rr_H_upper'])
    need(0 < combined < old_combined and new_smooth_h < old_smooth_h and
         paired_cap < Q(previous['paired_u_rr_H_upper']), 'both refined majorants and combined initial cap improve')
    vr = Q(higher['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
    need(vr == Q(previous['unchanged_full_initial_v_r_error_H_upper']) ==
         Q(new_paired['unchanged_full_initial_v_r_error_H_upper']), 'full P1 velocity gradient certificate retained')
    if a.combined_report is not None:
        need(a.combined_sha256 is not None, 'explicit combined report digest')
        combined_claim = read(a.combined_report, a.combined_sha256, 'combined initial producer report')
        common(combined_claim, labels, weights)
        need(combined_claim['status'] == 'PASS_COMBINED_UNCHANGED_INITIAL_CURVATURE_REFINEMENT' and
             combined_claim['refined_u_rr_H_cap_is_upper_bound'] is True and
             combined_claim['bound_quantity'] == 'sup_H_norm_actual_minus_smoothed_initial_u_rr',
             'combined report has complete uniform error upper role')
        need(combined_claim['input_sha256']['nominal'] == a.nominal_sha256 and
             combined_claim['input_sha256']['paired'] == a.paired_sha256,
             'combined report uses these same two authenticated refinements')
        need(Q(combined_claim['refined_uniform_initial_u_rr_H_upper']) == combined and
             Q(combined_claim['refined_curvature_smoothing_H_cap']) == new_smooth_h and
             Q(combined_claim['paired_u_rr_H_upper']) == paired_cap,
             'independent complete initial cap matches separate combined implementation')
        need(list(map(Q, combined_claim['refined_nominal_M3u_component_caps'])) == new_m3 and
             list(map(Q, combined_claim['refined_curvature_smoothing_component_caps'])) == new_smooth and
             qmatrix(combined_claim['paired_channel_H_caps']) == channel_maxima and
             Q(combined_claim['unchanged_full_initial_v_r_error_H_upper']) == vr,
             'all combined majorant components and retained velocity agree')

    # Corrupted structural/numerical copies must fail under normal and -O.
    bad = copy.deepcopy(new_nominal)
    del bad['refinement']['leaves'][len(bad['refinement']['leaves'])//2]
    reject('omitted nominal leaf', lambda: structural_nominal(bad, nominal, labels, weights))
    bad = copy.deepcopy(new_nominal)
    bad['branch'] = 'cutoff'
    reject('wrong nominal branch', lambda: structural_nominal(bad, nominal, labels, weights))
    bad = copy.deepcopy(new_nominal)
    bad['actual_order_six_bound_claimed'] = True
    reject('substitute actual sixth derivative claim', lambda: structural_nominal(bad, nominal, labels, weights))
    bad = copy.deepcopy(new_nominal)
    cell = next(row for row in bad['complete_cover'] if row['role'] == 'retained' and Q(row['component_caps'][12]) > 0)
    cell['component_caps'][12] = '0'
    reject('changed retained nominal cell', lambda: structural_nominal(bad, nominal, labels, weights))
    bad = copy.deepcopy(new_nominal)
    bad['nominal_M3_upper'] = False
    reject('nominal bound has wrong role', lambda: structural_nominal(bad, nominal, labels, weights))
    bad = copy.deepcopy(new_paired)
    del bad['refinement']['leaves'][len(bad['refinement']['leaves'])//2]
    reject('omitted paired leaf', lambda: structural_paired(bad, previous, higher, labels, weights))
    bad = copy.deepcopy(new_paired)
    ch, jet, component = next((ch, jet, j) for ch in range(3) for jet in range(2) for j in range(23)
                              if Q(bad['refinement']['leaves'][0]['inherited_paired_channel_component_caps'][ch][jet][j]) > 0)
    bad['refinement']['leaves'][0]['inherited_paired_channel_component_caps'][ch][jet][component] = '0'
    reject('arbitrary inherited paired cap', lambda: structural_paired(bad, previous, higher, labels, weights))
    bad = copy.deepcopy(new_paired)
    bad['Gram'][0] = str(2*weights[0])
    reject('wrong Gram weight', lambda: structural_paired(bad, previous, higher, labels, weights))
    bad = copy.deepcopy(new_paired)
    bad['candidate_changed'] = True
    reject('changed candidate identity', lambda: structural_paired(bad, previous, higher, labels, weights))
    bad = copy.deepcopy(new_paired)
    bad['status'] = 'NOT_PASS'
    reject('NOT_PASS status', lambda: structural_paired(bad, previous, higher, labels, weights))
    first_nonzero = next(z for row in nominal_raw for z in row if z)
    reject('omitted ordinary third factorial', lambda: need(first_nonzero == first_nonzero/6, 'nominal raw jet retains3!'))
    reject('omitted radial ell cube', lambda: need(first_nonzero == first_nonzero/ell**3, 'nominal raw jet retainsell^3'))
    reject('omitted inner curvature jump allowance', lambda: need(new_smooth == [2*half_first*delta*z for z in new_m3], 'convolution retainsinnerjump'))
    reject('maximum replaces telescoping sum', lambda: need(paired_cap == max(row[0] for row in channel_maxima), 'paired triangle sum required'))
    reject('paired cap replaces complete smoothed error', lambda: need(combined == paired_cap, 'smoothing allowance required'))
    h.authenticate()
    need(core.read_bytes() == raw_core, 'unchanged retained core readback')
    bytes_checked(nominal_path, PINS['nominal_initial'], 'nominal source readback')
    bytes_checked(profiles_path, PINS['initial_profiles'], 'profile parity source readback')
    bytes_checked(hi/'reference/INITIAL.json', PINS['higher_initial'], 'higher initial source readback')
    bytes_checked(hi/'reference/REPAIR.json', PINS['repair'], 'quadratic extension source readback')
    bytes_checked(previous_path, PINS['previous_refinement'], 'PR103 source readback')
    bytes_checked(a.nominal_report, a.nominal_sha256, 'new nominal report readback')
    bytes_checked(a.paired_report, a.paired_sha256, 'new paired report readback')
    if a.combined_report is not None:
        bytes_checked(a.combined_report, a.combined_sha256, 'combined report readback')
    result = {
        'status': 'PASS_INDEPENDENT_COMBINED_INITIAL_CURVATURE_AUDIT',
        'candidate': previous['candidate'], 'same_candidate_W1_quadratic': True, 'candidate_changed': False,
        'target_unchanged': True, 'initial_data_unchanged': True, 'F_G_and_coefficients_unchanged': True,
        'labels': labels, 'Gram': weights, 'gram': weights, 'ell': ell,
        'fixed_displacement_kernel_width': delta, 'fixed_velocity_kernel_width': Q(1, 2048),
        'bound_quantity': 'sup_H_norm_actual_minus_smoothed_initial_u_rr',
        'refined_u_rr_H_cap_is_upper_bound': True, 'actual_error_reduction_claimed': False,
        'input_sha256': dict(PINS, nominal_report=a.nominal_sha256, paired_report=a.paired_sha256,
                             combined_report=a.combined_sha256),
        'nominal_leaves_independently_recomputed': len(nominal_leaves),
        'nominal_old_parents': 189, 'nominal_untouched_parents': 154,
        'nominal_complete_cover_cells': len(nominal_cover),
        'paired_leaves_independently_recomputed': len(paired_leaves),
        'paired_untouched_old_cells': 47, 'paired_complete_cover_cells': len(paired_cover),
        'refined_nominal_M3u_component_caps': new_m3,
        'refined_curvature_smoothing_component_caps': new_smooth,
        'refined_curvature_smoothing_H_cap': new_smooth_h,
        'paired_channel_H_caps': channel_maxima, 'paired_u_rr_H_upper': paired_cap,
        'refined_uniform_initial_u_rr_H_upper': combined,
        'new_uniform_initial_u_rr_H_upper': combined,
        'previous_uniform_initial_u_rr_H_upper': old_combined,
        'initial_upper_certificate_reduction_factor': old_combined/combined,
        'unchanged_full_initial_v_r_error_H_upper': vr,
        'producer_code_imported': False, 'actual_order_six_bound_used': False,
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Independent replay using the retained authenticated analytic core plus full-cover/convolution aggregation. Conditional absolute scaled-coordinate initial-error bound for unchanged W1_quadratic; no true error reduction, pointwise relative or physical accuracy, nonlinear/cosmological realization.'}
    print(json.dumps(encode(result), sort_keys=True))


if __name__ == '__main__':
    main()
