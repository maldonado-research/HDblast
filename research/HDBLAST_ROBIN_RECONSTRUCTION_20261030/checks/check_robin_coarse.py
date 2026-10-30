#!/usr/bin/env python3
"""Independent local boundary reconstruction bound from the old candidate residual.

No new trace-moment or boundary-resummation implementation is imported.
The authenticated previous candidate point and local cone certificates are inputs.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def sqrt_upper(q):
    need(q >= 0, 'nonnegative squared norm')
    n = math.isqrt(q.numerator // q.denominator)
    if n*n*q.denominator < q.numerator:
        n += 1
    need(F(n*n) >= q, 'outward square root')
    return F(n)


def interval(row):
    lo, hi = map(F, row)
    need(lo <= hi, 'ordered interval')
    return lo, hi


def widen(row, error):
    need(error >= 0, 'positive interval widening')
    return row[0] - error, row[1] + error


def encode(row):
    return list(map(str, row))


parser = argparse.ArgumentParser()
parser.add_argument('--repo-root', type=Path, required=True)
parser.add_argument('--contact-capsule', type=Path)
args = parser.parse_args()
root = args.repo_root.resolve()
contact = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
old = read(contact/'reference/SHELL.json',
           '9de23a17dfb96361e4d5ad46f31b60cb0da458492d6e21d26c9dd6890e0b1d64')
coarse = read(contact/'reference/COARSE.json',
              '2b8daad3410478d180ba46cffd07609ff57fce8f8e3fd2f77256ea4a0fe2214c')
source = read(contact/'reference/AFFINE_SOURCE.json',
              'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e')
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
profiles = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
                '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = read(base/'compute/response-normal/RESPONSE_REPORT.json',
                'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
initial = read(base/'initial/results/INITIAL_REPORT.json',
               'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
need(all(d['status'] == 'PASS' and d['point'] == ['9/8', '0']
         for d in (old, coarse, source)), 'previous successful sample')
need(response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'formal response status')
need('zero-wedge factors3/4 bulk and1/4 Robin' in response['full_H_source_and_contraction_scope']
     and old['zero_wedge_value_error_factor'] == coarse['zero_wedge_value_error_factor'] == '1/8',
     'inherited actual zero-wedge scope')
need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER',
     'initial certificate status')
labels = profiles['labels']
weights = list(map(F, profiles['gram']))
need(len(labels) == len(weights) == 23 and labels == old['labels'] == response['labels']
     == initial['labels'] == [z['label'] for z in source['rows']], 'registered space')
need(weights == list(map(F, response['gram'])) == list(map(F, initial['gram']))
     and min(weights) > 0, 'positive common Gram weights')
parity = [-1 if name.startswith(('L:', 'v:')) else 1 for name in labels]
need(parity == old['parity'] and parity.count(-1) == 5, 'mixed parity')
blocks = [[i for i, name in enumerate(labels) if name.endswith(':' + degree)]
          for degree in ('0', '2', '4', '+', '-')]
matrices = {name: [list(map(F, row)) for row in profiles['matrices'][name]]
            for name in ('A', 'B', 'D', 'M')}
norms = {name: F(response['constant_matrix_H_operator_upper'][name])
         for name in matrices}
for name, matrix in matrices.items():
    for i in range(23):
        for j in range(23):
            need(any(i in block and j in block for block in blocks) or matrix[i][j] == 0,
                 'off-block frozen coefficient')
    for block in blocks:
        need(sum(weights[i]*matrix[i][j]**2/weights[j]
                 for i in block for j in block) <= norms[name]**2,
             'frozen block operator cap')
M0 = matrices['M']
deltaM = F(response['coefficient_error_H_operator_upper']['M'])
entries = [list(map(F, row)) for row in response['coefficient_entry_error_caps']['M']]
actualM = [[interval(z) for z in row] for row in response['actual_coefficient_boxes']['M']]
for i in range(23):
    for j in range(23):
        lo, hi = actualM[i][j]
        need(max(abs(lo-M0[i][j]), abs(hi-M0[i][j])) <= entries[i][j],
             'actual coefficient error entry')
        need(any(i in block and j in block for block in blocks) or entries[i][j] == 0,
             'actual off-block coefficient error')
        if parity[i] < 0 or parity[j] < 0:
            need(M0[i][j] == 0 and (lo, hi) == (0, 0), 'actual/frozen Dirichlet rows/columns')
for block in blocks:
    need(sum(weights[i]*entries[i][j]**2/weights[j]
             for i in block for j in block) <= deltaM**2, 'coefficient error H cap')
nM = norms['M'] + deltaM
need(nM == F(coarse['actual_shell_M_H_cap']) and deltaM > 0,
     'actual versus frozen Robin uncertainty')
x, length = F(9, 8), F(1, 8)
du, dv = F(1, 2**30), F(1, 2048)
nodes = list(map(F, profiles['nodes']))
tabs = {name: [list(map(F, row)) for row in profiles['profiles'][name]]
        for name in ('u', 'du', 'v')}
need(nodes[1] - max(du, dv) > 1 and
     all(z == 0 for table in tabs.values() for row in table[:2] for z in row),
     'inherited zero collar')
# Recompute local history caps directly from nodal convexity and old continuous
# interpolation certificates. These cover every shell time 0<=t<=x, not merely x.
caps = {name: [F(0)]*23 for name in tabs}
for k, cell in enumerate(initial['interpolation_cells']):
    need(list(map(F, cell['r'])) == nodes[k:k+2], 'cell geometry')
    if nodes[k] > x + max(du, dv):
        break
    for j in range(23):
        for name, error_index in (('u', 0), ('du', 1)):
            bound = max(abs(tabs[name][k][j]), abs(tabs[name][k+1][j]))
            caps[name][j] = max(caps[name][j], bound + F(cell['component_error_caps'][error_index][j]))
        caps['v'][j] = max(caps['v'][j], abs(tabs['v'][k][j]), abs(tabs['v'][k+1][j]))
U0, U1, Vcap = [sqrt_upper(sum(weights[i]*caps[name][i]**2 for i in range(23)))
               for name in ('u', 'du', 'v')]
need([U0, U1, Vcap] == [F(coarse[key]) for key in ('local_U0_H_cap', 'local_U1_H_cap', 'local_V_H_cap')],
     'independently recomputed local profile caps')
trace0 = U0 + length*Vcap
Fcap = (norms['A'] + norms['B'])*(U1 + Vcap) + norms['D']*(U0 + x*Vcap)
Kcap = (x*x - 1)*Fcap/2
B1cap = Kcap + length*norms['M']*trace0
need(B1cap == F(coarse['B1_value_H_cap']), 'uniform history first-feedback cap')
need(Kcap > length*length*Fcap/2, 'active cone area is not a short triangle')
Gcap = F(source['source_H_norm_strict_integer_upper'])
need(F(source['source_H_norm_square_upper']) < Gcap*Gcap and Gcap == 487122,
     'inherited uniform actual source cap')
# R_V=-M0 B1-(Mactual-M0)V-G. All caps below are uniform over 0<=t<=x.
# The actual source G=0 on t<=1 is inherited from the actual C1 formal
# response's zero wedge and its exact Robin law, not from the magnitude cap.
Rcap = nM*B1cap + deltaM*trace0 + Gcap
contraction = nM*length
need(0 < contraction < F(1, 10000), 'small but nonzero Robin-history contraction')
qcap = Rcap/(1-contraction)
need(qcap >= Rcap and qcap == Rcap + contraction*qcap, 'resolvent bound arithmetic')
# q+M Jq=-R_V, C=Lq. The corresponding exact reconstruction has the same
# initial data and principal bulk forcing, and satisfies the actual Robin law.
# q has no Dirichlet output, hence all three PD corrections are identically zero.
old_error = F(old['C1_error_readable_upper'])
need(old_error == F(2184*10**25), 'old smooth candidate error bound')
new_error_local = old_error + qcap
rows = []
for i, name in enumerate(labels):
    previous = old['full23_point_table'][i]
    need(previous['label'] == name and previous['parity'] == parity[i], 'old point row identity')
    factor = sqrt_upper(1/weights[i])
    zero_Robin_row = parity[i] > 0 and all(z == (0, 0) for z in actualM[i])
    if zero_Robin_row:
        # Here R_V=-G and q=G exactly. This is a componentwise source bound,
        # rather than a conversion of the full-H correction norm.
        derivative_radius = F(source['rows'][i]['abs_upper'])
        value_radius = length*derivative_radius
    else:
        value_radius = length*qcap*factor if parity[i] > 0 else F(0)
        derivative_radius = qcap*factor if parity[i] > 0 else F(0)
    reconstructed = {key: widen(interval(previous['candidate'][key]),
                                value_radius if key == 'value' else derivative_radius)
                     for key in ('value', 'dx', 'dr')}
    actual = {key: widen(reconstructed[key], factor*new_error_local*(length if key == 'value' else 1))
              for key in ('value', 'dx', 'dr')}
    # The generic triangle transfer is conservative. On Robin rows the shared
    # exact law for U and W cancels G, so the radial error is M(U-W).
    radial_shared = (reconstructed['dr'] if zero_Robin_row else
                    widen(reconstructed['dr'], factor*nM*length*new_error_local)) if parity[i] > 0 else None
    if parity[i] < 0:
        reconstructed['value'] = reconstructed['dx'] = (F(0), F(0))
        actual['value'] = actual['dx'] = (F(0), F(0))
    rows.append({'label': name, 'parity': parity[i],
                 'zero_Robin_row_source_only_correction': zero_Robin_row,
                 'old_candidate': previous['candidate'],
                 'correction_H_component_value_radius': str(value_radius),
                 'correction_H_component_derivative_radius': str(derivative_radius),
                 'reconstructed_boundary_candidate': {key: encode(value) for key, value in reconstructed.items()},
                 'conditional_actual_formal_response': {key: encode(value) for key, value in actual.items()},
                 'conditional_actual_formal_shared_Robin_radial': encode(radial_shared) if radial_shared else None})
mean = rows[12]
need(mean['label'] == 'w:0', 'primary scalar mean')
need(all(F(mean['reconstructed_boundary_candidate'][key][0]) > 0 for key in ('value', 'dx', 'dr')),
     'reconstructed scalar mean jet positive')
need(all(F(mean['conditional_actual_formal_response'][key][0]) > 0 for key in ('value', 'dx'))
     and F(mean['conditional_actual_formal_shared_Robin_radial'][0]) > 0,
     'conditional actual scalar mean jet positive')
for i in (21, 22):
    need(rows[i]['zero_Robin_row_source_only_correction'] and
         rows[i]['old_candidate'] == {key: ['0', '0'] for key in ('value', 'dx', 'dr')},
         'TT zero-old-candidate source-only correction')
    cap = F(source['rows'][i]['abs_upper'])
    need(interval(rows[i]['reconstructed_boundary_candidate']['dr']) == (-cap, cap)
         and interval(rows[i]['conditional_actual_formal_shared_Robin_radial']) == (-cap, cap),
         'TT actual inward derivative is the actual source')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'method': 'independent residual-resolvent norm enclosure around authenticated previous point; no new trace-moment evaluator import',
    'history_domain': '0<=t<=9/8; zero for t<=1',
    'history_length': str(length), 'local_U0_H_cap': str(U0), 'local_U1_H_cap': str(U1),
    'local_V_H_cap': str(Vcap), 'B0_history_trace_H_cap': str(trace0),
    'B1_history_value_H_cap': str(B1cap), 'actual_M_H_cap': str(nM),
    'actual_minus_frozen_M_H_cap': str(deltaM), 'G_H_cap': str(Gcap),
    'old_candidate_Robin_residual_H_cap': str(Rcap),
    'Robin_history_contraction': str(contraction), 'correction_q_H_cap': str(qcap),
    'point_correction_value_H_cap': str(length*qcap),
    'point_correction_derivative_H_cap': str(qcap),
    'local_C1_actual_minus_reconstructed_bound': str(new_error_local),
    'full23_point_table': rows,
    'boundary_identity': 'PN W_r=Mactual PN W+G exactly for the integral-defined reconstruction; interval arithmetic alone is not an independently evaluated zero residual',
    'principal_bulk_identity': 'W_xx-W_rr=f0, where f0 uses B0 and frozen bulk matrices',
    'scope': 'conditional boundary reconstruction and one-point enclosure sharing prior certified candidate/source inputs; not full bulk self-consistency, numerical C2 error, nonlinear spacetime or physical observable'}, sort_keys=True))
