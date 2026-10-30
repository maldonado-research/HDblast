#!/usr/bin/env python3
"""Cross-report arithmetic and a uniform forcing-defect bound for W1."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def interval(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(F, row)
    need(lo <= hi, 'interval order')
    return lo, hi


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', required=True, type=Path)
for name in ('bound', 'point', 'local'):
    ap.add_argument('--' + name + '-report', type=Path)
args = ap.parse_args()
ref = Path(__file__).resolve().parents[1] / 'reference'
b = read(args.bound_report or ref/'BULK_BOUND.json',
         '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f')
p = read(args.point_report or ref/'POINT_RESIDUAL.json',
         '30dc43f1ee1f1830976fefffc589beed42431d8013a1218094f9ecde0327ec69')
l = read(args.local_report or ref/'LOCAL.json',
         '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899')
profile = read(args.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/compute/response-normal/RESPONSE_PROFILES.json',
               '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
need(b['status'] == p['status'] == l['status'] == 'PASS', 'successful core reports')
need(p['point'] == l['point'] == ['9/8', '0'], 'same shell point')
labels, weights = profile['labels'], list(map(F, profile['gram']))
need(labels == b['labels'] and len(labels) == len(weights) == 23 and min(weights) > 0,
     'registered positive comparison norm')
need(len(p['full23_point_table']) == len(l['full23_point_table']) == 23,
     'complete point tables')
lower_square = F(0)
nonzero = []
for i, row in enumerate(p['full23_point_table']):
    need(row['label'] == labels[i], 'old defect row identity')
    lo, hi = interval(row['complete_reduced_bulk_residual'])
    minimum = lo if lo > 0 else -hi if hi < 0 else F(0)
    lower_square += weights[i]*minimum**2
    if minimum > 0:
        nonzero.append(labels[i])
need(lower_square == F(p['complete_residual_H_norm_squared_lower']) and
     lower_square > F(9*10**22)**2 and nonzero == p['nonzero_residual_rows'] and
     len(nonzero) == 17, 'old complete forcing defect lower norm and signs')

# With the actual-RHS-minus-principal convention, r1=La(W1-W0).
# The zero-wedge lemma gives |W1-W0|<=3 X(W1-W0)/4.
ops = {key: F(value) for key, value in b['actual_coefficient_H_norm_caps'].items()}
coefficient = ops['A']+ops['B']+F(3, 4)*ops['D']
need(F(3, 2)*coefficient == F(b['bulk_K_La_X_operator_cap']),
     'raw forcing and propagated response units agree')
new_uniform_defect_H = coefficient*F(b['first_actual_bulk_correction_X_cap'])
need(new_uniform_defect_H < F(3574*10**19), 'readable next forcing-defect cap')
need(lower_square > (F(5, 2)*new_uniform_defect_H)**2,
     'next uniform defect cap is smaller than old point lower norm by over 2.5')
need(F(l['global_one_update_error_exact_retained_cap']) ==
     F(b['one_update_retained_target_error_X_cap']) and
     F(l['global_one_update_error_readable_strict_upper']) ==
     F(b['one_update_readable_error_strict_upper']), 'same global W1 error certificate')
need(l['updated_best_point_box_containment_checks'] == 69 and
     l['strictly_narrower_best_point_boxes'] == 53,
     'local target point comparison results')
for i, row in enumerate(l['full23_point_table']):
    need(row['label'] == labels[i], 'local target row identity')
    for key in ('value', 'dx', 'dr'):
        interval(row['conditional_actual_target_U_best_point_boxes'][key])
mean = l['full23_point_table'][12]
need(mean['label'] == 'w:0', 'registered scalar mean')
for key, low, high in [('value', F(829036*10**25), F(829042*10**25)),
                      ('dx', F(1974621*10**26), F(1974623*10**26)),
                      ('dr', F(2154067*10**21), F(2154080*10**21))]:
    lo, hi = interval(mean['conditional_actual_target_U_best_point_boxes'][key])
    need(low < lo <= hi < high, 'strict readable target ' + key + ' bounds')
print(json.dumps({'status': 'PASS',
    'old_W0_forcing_defect_point_H_norm_strict_lower': str(9*10**22),
    'next_W1_forcing_defect_uniform_H_cap': str(new_uniform_defect_H),
    'next_W1_forcing_defect_uniform_H_readable_strict_upper': '3.574e22',
    'uniform_new_cap_vs_old_point_lower_factor_strict_lower': '2.5',
    'forcing_defect_convention': 'r1=La W1+F-(La W0+F)=La(W1-W0)',
    'global_W1_error_readable_strict_upper': b['one_update_readable_error_strict_upper'],
    'local_target_point_boxes_compared': 69, 'strictly_narrower_local_target_boxes': 53,
    'scope': 'conditional reduced forcing-defect majorant and registered target enclosures; no evaluation of the bulk convolution, quantitative C2 error or pointwise Einstein/constraint residual'}, sort_keys=True))
