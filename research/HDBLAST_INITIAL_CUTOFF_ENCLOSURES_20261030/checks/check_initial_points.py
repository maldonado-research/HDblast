#!/usr/bin/env python3
"""Exact interval diagnostic of actual-minus-nominal initial differences.

Point intervals diagnose looseness of continuous cell upper bounds. They do
not replace continuous certificates or prove where the global maximum lies.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib.util, json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root / 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
producer = base / 'full23_initial.py'
pin(producer, '2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1')
old = json.loads(pin(base / 'results/INITIAL_REPORT.json',
    'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
spec = importlib.util.spec_from_file_location('initial_diagnostic_producer', producer)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
ctx = mod.Context(root)
need(old['labels'] == list(ctx.op.LABELS), 'modal labels')
labels = old['labels']
weights = list(map(F, old['gram']))
need(labels[12] == 'w:0' and weights == list(ctx.weights), 'registered scalar and Gram')
left, right = F(263, 256), F(33, 32)
oldrow = old['rows'][7]
need(list(map(F, oldrow['r'])) == [left, right], 'old last cutoff cell')
points = [left + (right-left)*F(k, 8) for k in range(9)]
points += [F(573, 512)]  # midpoint of second-ranked old continuous cell
rows = []
for r in points:
    paired = ctx.paired_profiles(r)
    channel_boxes = [[[(item.d[c].val().lo, item.d[c].val().hi)
                       for item in row] for row in paired] for c in range(3)]
    for c in range(3):
        for kind in range(3):
            need(all(lo <= hi for lo, hi in channel_boxes[c][kind]), 'ordered channel box')
    radial = channel_boxes[1][1]
    radial_caps = [max(abs(lo), abs(hi)) for lo, hi in radial]
    channel_H = [[ctx.norm([max(abs(lo), abs(hi)) for lo, hi in row])
                   for row in channel] for channel in channel_boxes]
    lower_norm2 = sum(w*(min(abs(lo), abs(hi)) if lo*hi > 0 else 0)**2
                      for w, (lo, hi) in zip(weights, radial))
    # The three paired channels telescope to the actual-minus-nominal data.
    # Only their interval sum, rather than one channel by itself, can support
    # a lower bound for the complete initial-data difference at a point.
    total_boxes = [[(sum(channel_boxes[c][kind][i][0] for c in range(3)),
                     sum(channel_boxes[c][kind][i][1] for c in range(3)))
                    for i in range(23)] for kind in range(3)]
    total_upper = [ctx.norm([max(abs(lo), abs(hi)) for lo, hi in row])
                   for row in total_boxes]
    total_lower2 = [sum(w*(min(abs(lo), abs(hi)) if lo*hi > 0 else 0)**2
                        for w, (lo, hi) in zip(weights, row)) for row in total_boxes]
    if r <= right:
        need(all(max(abs(lo), abs(hi)) <= F(oldrow['channel_component_error_caps'][1][1][i])
                 for i, (lo, hi) in enumerate(radial)), 'point inside old cell cap')
    rows.append({'r': str(r),
        'rational_packet_channel_radial_component_intervals':
            [[str(lo), str(hi)] for lo, hi in radial],
        'rational_packet_channel_radial_H_lower_square': str(lower_norm2),
        'channel_H_upper': [[str(z) for z in row] for row in channel_H],
        'total_actual_minus_nominal_initial_component_intervals':
            [[[str(lo), str(hi)] for lo, hi in row] for row in total_boxes],
        'total_actual_minus_nominal_initial_H_upper': list(map(str, total_upper)),
        'total_actual_minus_nominal_initial_H_lower_square': list(map(str, total_lower2))})
peak = next(row for row in rows if row['r'] == '1055/1024')
peak_box = list(map(F, peak['total_actual_minus_nominal_initial_component_intervals'][1][12]))
need(F(-35495*10**12) < peak_box[0] <= peak_box[1] < F(-35494*10**12),
     'readable nonzero total radial scalar-mean difference sample')
need(F(403*10**14)**2 < F(peak['total_actual_minus_nominal_initial_H_lower_square'][1])
     and F(peak['total_actual_minus_nominal_initial_H_upper'][1]) < F(4031*10**13),
     'readable nonzero total radial norm sample')
ctx.reauthenticate()
print(json.dumps({'status': 'PASS', 'labels': labels, 'points': rows,
    'dominant_old_cell': ['263/256', '33/32'],
    'dominant_old_radial_channel_H_cap': old['channel_H_caps'][1][1],
    'nonzero_total_radial_difference_sample': {
        'r': '1055/1024', 'label': 'w:0',
        'strict_component_bounds': ['-35495*10^12', '-35494*10^12'],
        'strict_H_norm_bounds': ['403*10^14', '4031*10^13']},
    'point_diagnostic_scope': 'authenticated pointwise difference intervals only; not a replacement continuous initial certificate or proof of a global maximum',
    'candidate_changed': False}, sort_keys=True))
