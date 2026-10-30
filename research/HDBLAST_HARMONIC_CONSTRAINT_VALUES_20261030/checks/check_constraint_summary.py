#!/usr/bin/env python3
"""Exact report audit and all-angular-point conversion for harmonic values."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


PINS = {
    'FORMULAS': 'f3f037f890655ab05f936491d378a494b56e10f0c8a3d32dda0f1b8b2ea7556f',
    'TRANSFER': 'f809e0b5016f610bc717648a5feb58019d8b7a335d094ec5647740e15213dd73',
    'DIRECT': '08c2849ea32071821a1b4cd5550c20ca3031e2ed538148e22d9a59f17379c872',
}


def iv(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(Q, z)
    need(a <= b, 'ordered interval')
    return a, b


ap = argparse.ArgumentParser()
ap.add_argument('--capsule', type=Path, default=Path(__file__).resolve().parent.parent)
args = ap.parse_args()
reports = {}
for name, digest in PINS.items():
    raw = (args.capsule/'reference'/(name+'.json')).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+name)
    reports[name] = json.loads(raw)
f, t, d = [reports[k] for k in ('FORMULAS', 'TRANSFER', 'DIRECT')]
need(f['status'] == 'PASS_EXACT_HARMONIC_CONSTRAINT_FORMULAS' and
     t['status'] == d['status'] == 'PASS', 'successful independent reports')
need(d['point'] == t['shell_point']['point'] == ['9/8', '0'], 'same selected shell point')
need(f['input_labels'] == t['input_labels'] and f['input_Gram'] == t['input_Gram'] and
     f['output_labels'] == t['output_labels'] and f['output_Gram'] == t['output_Gram'], 'independent operator/Gram agreement')
need(f['group_count'] == 27 and f['control_count'] == 7, 'independent algebra verification count')
weights = dict(zip(t['output_labels'], map(Q, t['output_Gram'])))
need(len(weights) == 8 and min(weights.values()) > 0, 'complete positive covector Gram')
rows = d['rows']
need(len(rows) == 8 and {r['label'] for r in rows} == set(weights), 'all eight modal constraint rows')
raw_square = refined_square = Q(0)
for row in rows:
    raw, refined = iv(row['raw_interval']), iv(row['refined_interval'])
    need(raw[0] <= 0 <= raw[1] and refined[0] <= 0 <= refined[1], 'zero remains possible in every component')
    need(raw[0] <= refined[0] <= refined[1] <= raw[1], 'valid normal boundary refinement')
    weight = weights[row['label']]
    need(weight == Q(row['output_probability_sphere_Gram_weight']), 'independent covector output weight')
    raw_square += weight*max(abs(z) for z in raw)**2
    refined_square += weight*max(abs(z) for z in refined)**2
    if row['label'].startswith('R:'):
        need(refined == (0, 0), 'exact normal component at the shell')
    else:
        need(refined == raw, 'tangential components retain direct enclosures')
need(raw_square == Q(d['raw_direct_covector_norm_squared_upper']) and
     refined_square == Q(d['boundary_refined_covector_norm_squared_upper']), 'independent full covector interval norm sum')
need(d['nonzero_constraint_intervals'] == [] and d['Dirichlet_normal_derivatives_retained'] is True,
     'no nonzero sign inferred; normal derivatives retained')
need(t['shell_point']['exact_radial_component_error'] == ['0', '0'], 'shared-boundary normal error identity')

# Scalar hull dimensions 1+9+25 and gradient trace 9+25. Homogeneity and
# isotropy prove the kernel form in the companion mathematical derivation.
scalar_kernel = sum((l+1)**2 for l in (0, 2, 4))
gradient_kernel = Q(sum((l+1)**2 for l in (2, 4)), 3)
need(scalar_kernel == 35 and gradient_kernel == Q(34, 3) < scalar_kernel and
     Q(f['pointwise_from_sphere_L2_cap_squared']) == scalar_kernel, 'pointwise five-covector addition-kernel maximum')
direct_H = Q('1.319e24')
direct_point = Q('7.801e24')
need(refined_square < direct_H**2 and 35*refined_square < direct_point**2,
     'strict direct probability-sphere and every-angular-point bounds')
global_point, shell_transfer_point = Q('3.017e30'), Q('2.943e30')
need(35*Q(t['global']['frame_constraint_H_cap'])**2 < global_point**2 and
     35*Q(t['shell_point']['frame_constraint_H_cap'])**2 < shell_transfer_point**2,
     'strict conditional C1-transfer every-angular-point bounds')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'normal_harmonic_component_at_shell': 'exactly zero under same actual boundary/source identities',
    'raw_modal_intervals_containing_zero': 8,
    'nonzero_tangential_constraint_signs_proved': 0,
    'direct_positive_sphere_L2_strict_upper': str(direct_H),
    'direct_every_angular_point_positive_frame_strict_upper': str(direct_point),
    'conditional_global_every_angular_point_positive_frame_strict_upper': str(global_point),
    'conditional_shell_every_angular_point_positive_frame_strict_upper': str(shell_transfer_point),
    'angular_evaluation_kernel': {'normal': '35', 'time': '35', 'sphere_vector': '34/3'},
    'direct_requires_formal_target_zero_constraint': False,
    'C1_error_transfer_requires_formal_target_zero_constraint': True,
    'input_sha256': PINS,
    'scope': 'harmonic covector VALUES in a declared positive background frame norm; no exact tangential zero, relative smallness, constraint derivative, full Hamiltonian/momentum or Einstein residual, quantitative C2 error, nonlinear spacetime or Big Bang evidence'}, sort_keys=True))
