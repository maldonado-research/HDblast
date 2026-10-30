#!/usr/bin/env python3
"""Independently verify the released numerical statements from pinned reports."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pinned(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def iv(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(Q, z)
    need(a <= b, 'interval order')
    return a, b


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--capsule', type=Path)
ap.add_argument('--bulk-capsule', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
cap = args.capsule or root/'research/HDBLAST_CAUSAL_BULK_QUADRATURE_20261030'
bulk = args.bulk_capsule or root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
q = pinned(cap/'reference/QUADRATURE.json', '5bd5ce2d40fdf0805d06351160fe7d2e323324983b69aeee5732c1e28013c6f8')
t = pinned(cap/'reference/TRANSFER.json', 'ab80240cedec3807d40ac960457b384e68857acafa78ee9792730eb75a3dda99')
old = pinned(bulk/'reference/LOCAL.json', '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899')
r0 = pinned(bulk/'reference/POINT_RESIDUAL.json', '30dc43f1ee1f1830976fefffc589beed42431d8013a1218094f9ecde0327ec69')
p = pinned(root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/compute/response-normal/RESPONSE_PROFILES.json', '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
need(all(z['status'] == 'PASS' and z['point'] == ['9/8', '0'] for z in (q, t, old, r0)), 'same successful point')
labels = p['labels']
weights = list(map(Q, p['gram']))
need(len(labels) == len(weights) == 23 and min(weights) > 0, 'registered positive Gram')
signs = {k: [] for k in ('value', 'dx', 'dr')}
norm_squares = {k: Q(0) for k in signs}
containments = strict = 0
for i, (qr, tr, pr) in enumerate(zip(q['full23_point_table'], t['full23_point_table'], old['full23_point_table'])):
    need(qr['label'] == tr['label'] == pr['label'] == labels[i], 'registered modal order')
    for k in signs:
        c = iv(qr['correction'][k])
        need(c == iv(tr['correction'][k]), 'same correction transferred')
        if c[0] > 0 or c[1] < 0:
            signs[k].append(labels[i])
        norm_squares[k] += weights[i]*max(abs(z) for z in c)**2
        new = iv(tr['conditional_actual_target_U_best_point_boxes'][k])
        previous = iv(pr['conditional_actual_target_U_best_point_boxes'][k])
        need(previous[0] <= new[0] <= new[1] <= previous[1], 'same-target containment')
        containments += 1
        strict += new[1]-new[0] < previous[1]-previous[0]
need(containments == 69 and strict == 53 == t['strictly_narrower_best_target_boxes'], 'target box claims')
need({k: len(v) for k, v in signs.items()} == {'value': 16, 'dx': 16, 'dr': 17} == q['nonzero_correction_intervals'], 'correction sign claims')
need(signs == t['sign_determined_correction_rows'], 'same sign classifications')
need(labels[12] == 'w:0', 'registered scalar mean')
readable = {'value': ('1.596e21', '6.733e21'),
            'dx': ('8.355e22', '1.058e23'),
            'dr': ('4.147e17', '1.750e18')}
for k, (a, b) in readable.items():
    c = iv(q['full23_point_table'][12]['correction'][k])
    need(0 < Q(a) < c[0] <= c[1] < Q(b), 'strict scalar-mean correction inequality')
old_residual = iv(r0['full23_point_table'][12]['complete_reduced_bulk_residual'])
need(old_residual[1] < 0, 'opposite local residual sign')
Hcap = {'value': Q('7.646e21'), 'dx': Q('1.202e23'), 'dr': Q('1.061e21')}
for k, c in Hcap.items():
    need(norm_squares[k] < c*c, 'readable whole23 point correction norm cap')
need(Q(old['bulk_update_correction_local_value_and_trace_H_cap']) > 330*Hcap['value'] and
     Q(old['bulk_update_correction_shell_time_derivative_H_cap']) > 178*Hcap['dx'], 'comparison of upper certificates')
ratio = {k: Q(v) for k, v in t['scalar_mean_previous_to_best_width_ratios'].items()}
need(ratio['value'] > Q('1.125') and ratio['dr'] > Q('1.125') and ratio['dx'] > Q('1.330'), 'readable scalar-mean target width ratios')
need(t['unchanged_global_error_readable_strict_upper'] == old['global_one_update_error_readable_strict_upper'] == '35410000000000000000000000', 'global error theorem unchanged')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'correction_scalar_mean_strict_bounds': readable,
    'whole23_correction_point_H_strict_upper': {k: str(v) for k, v in Hcap.items()},
    'correction_sign_determined_component_counts': {k: len(v) for k, v in signs.items()},
    'target_point_box_containments': containments, 'strictly_narrower_target_boxes': strict,
    'scalar_mean_target_width_ratio_strict_lower': {'value': '1.125', 'dx': '1.330', 'dr': '1.125'},
    'interpretation': 'the negative local r0 scalar mean and positive causal C scalar mean show that the endpoint forcing-defect sign does not determine this accumulated correction sign',
    'scope': 'conditional finite-update correction and same-target point boxes; upper-certificate comparisons, not true-error reduction measurements or physical forces; global C1 theorem unchanged; no C2/full Einstein constraint certificate, nonlinear family or Big Bang evidence'}, sort_keys=True))
