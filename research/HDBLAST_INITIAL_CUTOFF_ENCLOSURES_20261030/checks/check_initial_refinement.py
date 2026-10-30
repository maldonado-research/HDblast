#!/usr/bin/env python3
"""Refine one dominant initial-error interval without changing any candidate.

All other interval caps are inherited from the pinned full23 initial report.
The same authenticated paired evaluator is used on exact dyadic subintervals.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path


PRODUCER_PIN = '2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1'
REPORT_PIN = 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'
TARGET = F(2*10**24)
ROW = 7
MAX_LEAVES = 128
MAX_DEPTH = 8


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def caps(value):
    need(len(value) == 3 and all(len(channel) == 3 for channel in value)
         and all(len(row) == 23 for channel in value for row in channel),
         'three channels by three kinds by23 components')
    out = [[[F(v) for v in row] for row in channel] for channel in value]
    need(all(v >= 0 for channel in out for row in channel for v in row),
         'nonnegative component caps')
    return out


def sqrtup(value):
    need(value >= 0, 'nonnegative norm square')
    if value == 0:
        return F(0)
    grid = 1 << 112
    ceiling = (value.numerator*grid*grid + value.denominator-1)//value.denominator
    out = F(math.isqrt(ceiling)+1, grid)
    need(out*out >= value, 'outward exact square root')
    return out


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', required=True, type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
producer = base/'full23_initial.py'
pin(producer, PRODUCER_PIN)
old = json.loads(pin(base/'results/INITIAL_REPORT.json', REPORT_PIN))
need(old['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER',
     'inherited initial scope')
need(old['fixed_mesh_denominator'] == 256 and old['continuous_error_cells'] == 192,
     'inherited cell roster')
need(old['footprint']['r'] == ['0', '7/4'] and
     old['footprint']['actual_zero_on'] == ['0', '1'] and
     old['footprint']['nominal_zero_on'] == ['0', '65/64'], 'initial footprint')
need(old['derivative_convention']['du'] == 'partial_r U' and
     old['derivative_convention']['v'] == 'partial_x U; x=t/ell',
     'registered derivative conventions')
labels = old['labels']
weights = list(map(F, old['gram']))
need(len(labels) == len(weights) == 23 and min(weights) > 0, 'positive modal norm')
oldrows = old['rows']
need(len(oldrows) == 192, 'complete inherited error rows')
rowcaps = [caps(row['channel_component_error_caps']) for row in oldrows]
old_channel = caps(old['channel_component_error_caps'])
previous = F(1)
for i, row in enumerate(oldrows):
    left, right = map(F, row['r'])
    need(left == previous and right-left == F(1, 256), 'old continuous cover')
    need(not any(left < join < right for join in (F(65, 64), F(66, 64))),
         'old analytic branch split')
    previous = right
need(previous == F(7, 4), 'old cover endpoint')
left0, right0 = map(F, oldrows[ROW]['r'])
need((left0, right0) == (F(263, 256), F(33, 32)), 'fixed dominant transition cell')


def aggregate(leaves):
    result = [[[F(0)]*23 for kind in range(3)] for channel in range(3)]
    allcaps = [row for i, row in enumerate(rowcaps) if i != ROW]
    allcaps += [leaf['channel_component_error_caps'] for leaf in leaves]
    for table in allcaps:
        for c in range(3):
            for k in range(3):
                for i in range(23):
                    result[c][k][i] = max(result[c][k][i], table[c][k][i])
    norms = [[sqrtup(sum(weights[i]*z*z for i, z in enumerate(row)))
              for row in channel] for channel in result]
    total = [sum(norms[c][kind] for c in range(3)) for kind in range(3)]
    return result, norms, total


oldleaf = {'path': '', 'r': [left0, right0],
           'channel_component_error_caps': rowcaps[ROW]}
reconstructed_old, oldnorms, oldtotal = aggregate([oldleaf])
need(reconstructed_old == old_channel and
     oldnorms == [[F(z) for z in row] for row in old['channel_H_caps']] and
     oldtotal == list(map(F, old['total_H_caps'])), 'reconstruct inherited aggregation')
need(oldtotal[1] > TARGET, 'material refinement target')
spec = importlib.util.spec_from_file_location('initial_refinement_pinned', producer)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
ctx = module.Context(root)  # authenticates all seven roots before science imports
need(ctx.op.LABELS == labels and ctx.weights == weights and
     list(module.CHANNELS) == old['channels'], 'same mathematical decomposition')


def evaluate(path, left, right):
    paired = ctx.paired_profiles(ctx.I(left, right))
    raw = [[[item.d[c].val().mag() for item in row] for row in paired]
           for c in range(3)]
    raw = caps(raw)
    # Both the new subinterval cap and inherited parent cap enclose the same
    # fixed difference. Their intersection is another valid magnitude cap.
    effective = [[[min(raw[c][k][i], rowcaps[ROW][c][k][i])
                   for i in range(23)] for k in range(3)] for c in range(3)]
    return {'path': path, 'r': [left, right],
            'raw_channel_component_error_caps': raw,
            'channel_component_error_caps': effective}


leaves = [oldleaf]
history = []
while True:
    component, norms, total = aggregate(leaves)
    if total[1] < TARGET:
        break
    need(len(leaves) < MAX_LEAVES, 'adaptive leaf budget exhausted')
    # Highest channel1 radial H-square first; lexicographic path breaks ties.
    score = lambda leaf: sum(weights[i]*z*z for i, z in enumerate(
        leaf['channel_component_error_caps'][1][1]))
    selected = sorted(leaves, key=lambda leaf: (-score(leaf), leaf['path']))[0]
    need(len(selected['path']) < MAX_DEPTH, 'adaptive depth budget exhausted')
    left, right = selected['r']
    mid = (left+right)/2
    children = [evaluate(selected['path']+'0', left, mid),
                evaluate(selected['path']+'1', mid, right)]
    leaves.remove(selected)
    leaves.extend(children)
    updated = aggregate(leaves)[2]
    need(all(a <= b for a, b in zip(updated, total)), 'refinement monotonicity')
    history.append({'split_path': selected['path'], 'r': [left, right],
                    'selected_channel1_du_H_square': score(selected),
                    'total_H_caps_after_split': updated})

leaves.sort(key=lambda leaf: leaf['r'][0])
endpoint = left0
for leaf in leaves:
    path = leaf['path']
    need(path and set(path) <= {'0', '1'}, 'nonempty binary leaf path')
    index, depth = int(path, 2), len(path)
    left = left0+(right0-left0)*F(index, 2**depth)
    right = left0+(right0-left0)*F(index+1, 2**depth)
    need(leaf['r'] == [left, right] and left == endpoint and left < right,
         'dyadic leaf cover')
    need('raw_channel_component_error_caps' in leaf, 'evaluated final leaf')
    endpoint = right
need(endpoint == right0, 'refined parent fully covered')
component, norms, total = aggregate(leaves)
need(total[1] < TARGET and total[1] < oldtotal[1]/100,
     'published radial improvement gate')
need(all(component[c][k][i] <= old_channel[c][k][i]
         for c in range(3) for k in range(3) for i in range(23)),
     'all componentwise caps nonincreasing')
coverage = [{'old_row_index': i, 'r': list(map(F, row['r'])),
             'refined': i == ROW, 'leaves': leaves if i == ROW else []}
            for i, row in enumerate(oldrows)]
ctx.reauthenticate()
pin(producer, PRODUCER_PIN)
pin(base/'results/INITIAL_REPORT.json', REPORT_PIN)
print(json.dumps(encode({
    'status': 'PASS',
    'method': 'exact dyadic adaptive subdivision using unchanged pinned paired initial evaluator',
    'input_sha256': {'initial_producer': PRODUCER_PIN, 'initial_report': REPORT_PIN},
    'labels': labels, 'gram': weights, 'channels': old['channels'],
    'derivative_convention': old['derivative_convention'],
    'same_candidate': True,
    'candidate_background_coefficients_sha256': old['candidate_background_coefficients_sha256'],
    'fixed_original_mesh_denominator': 256,
    'footprint': old['footprint'],
    'refined_old_row_index': ROW,
    'inherited_untouched_cell_count': 191,
    'refined_parent_leaf_count': len(leaves),
    'refined_parent_depth_max': max(len(leaf['path']) for leaf in leaves),
    'total_continuous_leaf_count': 191+len(leaves),
    'radial_total_H_target_strict_upper': TARGET,
    'adaptive_budget': {'max_leaves': MAX_LEAVES, 'max_depth': MAX_DEPTH},
    'adaptive_selection': 'maximum channel1 du Gram-square; binary path lexical tie-break',
    'adaptive_split_history': history,
    'coverage_cells': coverage,
    'channel_component_error_caps': component,
    'channel_H_caps': norms,
    'total_H_caps': total,
    'old_total_H_caps': oldtotal,
    'old_to_new_radial_upper_certificate_ratio': oldtotal[1]/total[1],
    'complete_closures_authenticated_before_and_after': True,
    'scope': 'conditional continuous initial-difference bounds for exactly the same frozen candidate on the retained footprint; tighter interval upper certificates, not measured changes in actual data or error; no nonlinear family, cosmological claim or numerical C2 error'
}), sort_keys=True))
