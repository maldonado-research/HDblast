#!/usr/bin/env python3
"""Exact rational majorant for a boundary-preserving actual-bulk update.

This verifies a conditional norm certificate for an integral-defined candidate.
It does not integrate the actual bulk correction or evaluate an affine source.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def read(path, digest):
    return json.loads(pin(path, digest))


def sqrtup(value):
    need(value >= 0, 'nonnegative norm square')
    n = math.isqrt(value.numerator // value.denominator)
    if n*n < value:
        n += 1
    need(Q(n*n) >= value, 'outward square root')
    return Q(n)


def interval(value):
    need(len(value) == 2, 'interval shape')
    lo, hi = map(Q, value)
    need(lo <= hi, 'interval order')
    return lo, hi


def caps(value):
    need(len(value) == 3 and all(len(c) == 3 for c in value)
         and all(len(k) == 23 for c in value for k in c), 'channel/kind/component shape')
    result = [[[Q(z) for z in k] for k in c] for c in value]
    need(all(z >= 0 for c in result for k in c for z in k), 'positive component caps')
    return result


def verify_initial_cover(initial, refined, labels, weights):
    need(refined['same_candidate'] is True and refined['status'] == 'PASS'
         and refined['labels'] == labels and list(map(Q, refined['gram'])) == weights,
         'same refined initial candidate and modal space')
    need(refined['channels'] == initial['channels']
         and refined['derivative_convention'] == initial['derivative_convention']
         and refined['footprint'] == initial['footprint']
         and refined['candidate_background_coefficients_sha256'] ==
             initial['candidate_background_coefficients_sha256'], 'same initial functions')
    oldrows, cover = initial['rows'], refined['coverage_cells']
    need(len(oldrows) == len(cover) == 192, 'complete inherited cover')
    maximum = [[[Q(0)]*23 for _ in range(3)] for _ in range(3)]
    untouched = leaves_count = parents_count = 0
    for k, (old, parent) in enumerate(zip(oldrows, cover)):
        lo, hi = interval(old['r'])
        need((lo, hi) == (Q(256+k, 256), Q(257+k, 256))
             and parent['old_row_index'] == k and interval(parent['r']) == (lo, hi),
             'inherited parent identity and footprint')
        inherited = caps(old['channel_component_error_caps'])
        tables = []
        if not parent['refined']:
            need(parent['leaves'] == [], 'untouched leaf list')
            tables.append(inherited)
            untouched += 1
        else:
            parents_count += 1
            cursor, paths = lo, set()
            for leaf in sorted(parent['leaves'], key=lambda row: Q(row['r'][0])):
                path = leaf['path']
                need(isinstance(path, str) and path and set(path) <= {'0', '1'}
                     and path not in paths, 'unique nonempty binary path')
                paths.add(path)
                step = (hi-lo)/2**len(path)
                a, b = interval(leaf['r'])
                need(a == lo+int(path, 2)*step and b == a+step and a == cursor
                     and lo <= a < b <= hi, 'exact gap-free dyadic cover')
                need(not any(a < j < b for j in (Q(65,64), Q(66,64))),
                     'single analytic branch')
                raw, used = (caps(leaf[key]) for key in
                             ('raw_channel_component_error_caps', 'channel_component_error_caps'))
                need(all(used[c][j][i] == min(raw[c][j][i], inherited[c][j][i])
                         for c in range(3) for j in range(3) for i in range(23)),
                     'intersection with valid inherited cap')
                tables.append(used)
                cursor = b
                leaves_count += 1
            need(cursor == hi, 'refined parent endpoint')
        for table in tables:
            for c in range(3):
                for j in range(3):
                    for i in range(23):
                        maximum[c][j][i] = max(maximum[c][j][i], table[c][j][i])
    need((untouched, parents_count, leaves_count) == (191, 1, 10), 'refined cover counts')
    need(maximum == caps(refined['channel_component_error_caps']), 'component maxima')
    H = [[Q(z) for z in row] for row in refined['channel_H_caps']]
    need(len(H) == 3 and all(len(row) == 3 for row in H), 'channel norm shape')
    for c in range(3):
        for j in range(3):
            need(H[c][j] >= 0 and sum(weights[i]*maximum[c][j][i]**2 for i in range(23))
                 <= H[c][j]**2, 'positive Gram norm cap')
    total = [sum(H[c][j] for c in range(3)) for j in range(3)]
    need(total == list(map(Q, refined['total_H_caps'])), 'channel triangle totals')
    return total


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--initial-capsule', type=Path)
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--source-report', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
initial_capsule = args.initial_capsule or root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030'
robin_capsule = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
pins = {
    'prior_response': 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
    'refined_initial': 'd97f7393cb9467e78387368754981e81cf763da043cab0a4df65978f52b72122',
    'initial_report': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'profiles': '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc',
    'response': 'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd',
    'Robin_transfer': '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52',
    'affine_shell': 'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e',
    'source_response_review': '69a7a4741d599aad07d952a9b48423bf9b06f18c2b2431f3de8e821134db610c',
    'common_shell_report': '559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a',
    'full23_shell_packing': '0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe'}
prior = read(initial_capsule/'reference/RESPONSE.json', pins['prior_response'])
refined = read(initial_capsule/'reference/REFINED.json', pins['refined_initial'])
initial = read(base/'initial/results/INITIAL_REPORT.json', pins['initial_report'])
profiles = read(base/'compute/response-normal/RESPONSE_PROFILES.json', pins['profiles'])
response = read(base/'compute/response-normal/RESPONSE_REPORT.json', pins['response'])
transfer = read(robin_capsule/'reference/TRANSFER.json', pins['Robin_transfer'])
source_path = args.source_report or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030/reference/AFFINE_SOURCE.json'
source = read(source_path, pins['affine_shell'])
pin(root/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md',
    pins['source_response_review'])
shell = read(root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json',
             pins['common_shell_report'])
packing = read(base/'shell/verified/FULL23_SHELL_NORMAL.json', pins['full23_shell_packing'])
need(all(d['status'] == 'PASS' for d in (prior, refined, transfer, source))
     and initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER'
     and response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'inherited scope')
need(shell['status'] == 'PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING'
     and packing['status'] == 'PASS' and source['input_sha256'] == {
        'COMMON_HARMONIC_SHELL_NORMAL.json': pins['common_shell_report'],
        'FULL23_SHELL_NORMAL.json': pins['full23_shell_packing']}, 'affine source provenance')
labels, weights = profiles['labels'], list(map(Q, profiles['gram']))
need(len(labels) == len(weights) == 23 and min(weights) > 0 and
     labels == initial['labels'] == response['labels'] == transfer['labels'] and
     weights == list(map(Q, initial['gram'])) == list(map(Q, response['gram'])),
     'same positive registered modal space')
need(response['domain'] == {'r_min': 'max(0,x-5/4)', 'x': ['0','3/2'],
                            'x_plus_r_max': '7/4'}, 'retained cone')
need(initial['footprint']['actual_zero_on'] == ['0','1'] and
     initial['footprint']['nominal_zero_on'] == ['0','65/64'], 'inherited zero collars')
need(initial['derivative_convention']['du'] == 'partial_r U' and
     initial['derivative_convention']['v'] == 'partial_x U; x=t/ell' and
     initial['derivative_convention']['second_derivatives'][0] == 'partial_r^2 U'
     and 'U and V are C1' in initial['regularity'], 'initial smoothing and derivative premises')
parity = [-1 if label.startswith(('L:', 'v:')) else 1 for label in labels]
need(profiles['parity']['u'] == profiles['parity']['v'] == parity and
     profiles['parity']['du'] == [-s for s in parity], 'compatible derivative reflection')
blocks = [[i for i, label in enumerate(labels) if label.endswith(':'+degree)]
          for degree in ('0','2','4','+','-')]
need(sum(map(len, blocks)) == 23, 'complete orthogonal blocks')
sameblock = lambda i, j: any(i in b and j in b for b in blocks)
op = {k: Q(response['constant_matrix_H_operator_upper'][k]) for k in ('A','B','D','M')}
error = {k: Q(response['coefficient_error_H_operator_upper'][k]) for k in op}
for key in op:
    frozen = [list(map(Q, row)) for row in profiles['matrices'][key]]
    entry = [list(map(Q, row)) for row in response['coefficient_entry_error_caps'][key]]
    actual_boxes = [[interval(z) for z in row] for row in response['actual_coefficient_boxes'][key]]
    need(profiles['matrices'][key] == response['constant_matrices'][key], 'frozen matrix agreement')
    need(len(frozen) == len(entry) == len(actual_boxes) == 23 and
         all(len(row) == 23 for matrix in (frozen, entry, actual_boxes) for row in matrix),
         'operator matrix shapes')
    for i in range(23):
        for j in range(23):
            lo, hi = actual_boxes[i][j]
            need(entry[i][j] >= 0 and max(abs(lo-frozen[i][j]), abs(hi-frozen[i][j]))
                 <= entry[i][j], 'actual-minus-frozen entry cap')
            need(sameblock(i, j) or frozen[i][j] == entry[i][j] == 0, 'off-block closure')
            if key == 'M' and (parity[i] < 0 or parity[j] < 0):
                need(frozen[i][j] == 0 and (lo, hi) == (0,0), 'Robin projector closure')
    for block in blocks:
        need(sum(weights[i]*frozen[i][j]**2/weights[j] for i in block for j in block)
             <= op[key]**2, 'frozen operator Gram cap')
        need(sum(weights[i]*entry[i][j]**2/weights[j] for i in block for j in block)
             <= error[key]**2, 'coefficient-error Gram cap')
actual = {k: op[k]+error[k] for k in op}
need({k: str(v) for k,v in actual.items()} == prior['actual_coefficient_H_norm_caps'],
     'recomputed actual coefficient norms')

# Reconstruct the same smoothed initial free-wave mismatch and size.
delta_u, delta_v = Q(1,2**30), Q(1,2048)
nodes = list(map(Q, profiles['nodes']))
vtab = [list(map(Q, row)) for row in profiles['profiles']['v']]
need(nodes[0] == 0 and nodes[-1] == Q(7,4) and
     min(b-a for a,b in zip(nodes,nodes[1:])) >= 2*delta_v and
     nodes[1]-delta_v > 1 and all(v == 0 for row in vtab[:2] for v in row),
     'velocity mesh and compatible smoothing bands')
slopes = [[(vtab[k+1][i]-vtab[k][i])/(nodes[k+1]-nodes[k]) for i in range(23)]
          for k in range(len(nodes)-1)]
jumps = [[b-a for a,b in zip(left,right)] for left,right in zip(slopes,slopes[1:])]
av = Q(35,256)*delta_v*sqrtup(max(sum(weights[i]*row[i]**2 for i in range(23))
                               for row in jumps))
M2 = Q(initial['nominal_second_derivative_H_caps'][0])
au, bu = M2*delta_u**2/18, Q(35,128)*M2*delta_u
nu, nur, nv = map(Q, initial['nominal_H_caps'])
eu, eur, ev = map(Q, initial['P1_H_error_caps'])
du, dur, dv = verify_initial_cover(initial, refined, labels, weights)
initial_errors = [du+au, dur+bu, dv+ev+av]
einitial = max((initial_errors[0]+Q(3,2)*initial_errors[2])/3,
               initial_errors[1]+initial_errors[2])
B0cap = max((nu+au+Q(3,2)*(nv+ev+av))/3, nur+bu+nv+ev+av)
need(einitial == Q(prior['initial_free_wave_X_error_cap']) and
     B0cap == Q(prior['smooth_B0_X_cap']), 'same smoothed initial free wave')

# Reconstruct the inherited global boundary correction and full bulk residual.
history = transfer['global_retained_history']
need(Q(history['shell_history_end']) == Q(5,4) and
     Q(history['active_history_length']) == Q(1,4), 'whole-cone boundary history')
Gcap = Q(source['source_H_norm_strict_integer_upper'])
need(source['uniform_domain'] == '0<=x<=3,r=0' and Gcap <= 487122 and
     Q(source['source_H_norm_square_upper']) < Gcap**2, 'actual affine shell cap')
Rboundary = (op['M']*Q(history['B1_shell_trace_H_cap']) +
             error['M']*Q(history['V_shell_trace_H_cap']) + Gcap)
Ccap = Rboundary/(1-actual['M']/4)
need(Ccap == Q(history['correction_C1_X_cap']), 'global boundary-correction arithmetic')
q0 = Q(3,2)*(op['A']+op['B']+Q(3,4)*op['D'])+op['M']/4
b = actual['M']/4
k = Q(3,2)*(actual['A']+actual['B']+Q(3,4)*actual['D'])
q = k+b
qretained = Q(response['ledger']['retained_full_H_contraction'])
need(0 <= b < q < qretained == Q(3409,8000000) < 1, 'verified full contraction')
B1cap = q0*B0cap
feedbackcap = B1cap+Ccap
Rcoefficient = (error['A']+error['B']+Q(3,4)*error['D'])*B0cap
Rfeedback = (actual['A']+actual['B']+Q(3,4)*actual['D'])*feedbackcap
SOURCE = Q(response['ledger']['omitted_full_H_forcing_X_upper'])
need(SOURCE == Q(52162547172519,2), 'inherited affine response units')
KFcap = SOURCE+Gcap
R = Q(3,2)*(Rcoefficient+Rfeedback)+KFcap
same = {'boundary_correction_whole_cone_X_cap': Ccap,
        'frozen_coupled_X_operator_cap': q0, 'actual_coupled_X_operator_cap': q,
        'retained_coupled_X_operator_cap': qretained,
        'frozen_first_feedback_whole_cone_X_cap': B1cap,
        'W_minus_B0_whole_cone_X_cap': feedbackcap,
        'lower_order_coefficient_residual_H_cap': Rcoefficient,
        'lower_order_feedback_residual_H_cap': Rfeedback,
        'affine_bulk_response_KF_X_cap': KFcap,
        'bulk_residual_response_X_cap': R}
need(all(value == Q(prior[name]) for name,value in same.items()),
     'independent reconstruction of previous residual ledger')

# B=L(M PN .) is inverted before each actual-bulk update. No matrices commute.
# R already bounds X(K r0); it must not be treated as a raw forcing H cap.
D0 = R/(1-b)
beta = k/(1-b)
betaretained = (qretained-b)/(1-b)
initial_majorant = einitial/(1-q)
initial_majorant_retained = einitial/(1-qretained)
need(0 < beta <= betaretained < 1 and beta < q and betaretained < qretained,
     'boundary-inverted actual-bulk contraction')
E0 = (einitial+R)/(1-q)
E0retained = (einitial+R)/(1-qretained)
need(E0 == Q(prior['direct_C1_X_error_cap']) and
     E0retained == Q(prior['direct_C1_X_error_using_retained_contraction_cap']),
     'n0 reproduces previous direct certificate')
table = []
for n in range(5):
    tail = beta**n*D0/(1-beta)
    tailretained = betaretained**n*D0/(1-betaretained)
    En = initial_majorant+tail
    Enretained = initial_majorant_retained+tailretained
    need(En <= Enretained and En > initial_majorant, 'positive update-tail majorants')
    if n:
        previous = table[-1]
        need(En == einitial/(1-b)+beta*Q(previous['actual_target_error_X_cap']) and
             Enretained == einitial/(1-b)+betaretained*Q(previous['retained_target_error_X_cap']),
             'independent initial-mismatch recurrence')
    table.append({'bulk_update_count': n,
        'smoothed_actual_bulk_fixed_point_tail_X_cap': str(tail),
        'retained_smoothed_actual_bulk_fixed_point_tail_X_cap': str(tailretained),
        'actual_target_error_X_cap': str(En),
        'retained_target_error_X_cap': str(Enretained),
        'next_update_correction_X_cap': str(beta**n*D0),
        'retained_next_update_correction_X_cap': str(betaretained**n*D0)})
E1, E1retained = (Q(table[1][name]) for name in
                  ('actual_target_error_X_cap','retained_target_error_X_cap'))
need(E1 == (einitial+beta*R)/(1-q) and
     E1retained == (einitial+betaretained*R)/(1-qretained), 'one-update closed-form bound')
readable = Q(3541*10**22)
need(E1 <= E1retained < readable < E0retained and
     E0retained/readable > Q(57,10), 'readable one-update improvement gate')
need(D0+E1 == E0 and D0+E1retained == E0retained,
     'norm-only target enclosure conservation')
need(all(Q(table[n+1]['retained_target_error_X_cap']) <
         Q(table[n]['retained_target_error_X_cap']) for n in range(4)), 'decreasing upper majorants')

print(json.dumps({
    'status': 'PASS', 'method': 'exact rational boundary-inverted actual-bulk update majorant',
    'input_sha256': pins, 'labels': labels,
    'norm': 'X=max(sup|e|H/3,sup|e_x|H,sup|e_r|H) on retained Dplus',
    'same_smoothed_nominal_initial_data': True,
    'actual_coefficient_H_norm_caps': {name:str(value) for name,value in actual.items()},
    'bulk_K_La_X_operator_cap': str(k),
    'Robin_L_M_trace_X_operator_cap': str(b),
    'actual_coupled_X_operator_cap': str(q),
    'retained_coupled_X_operator_cap': str(qretained),
    'boundary_inverse_X_operator_cap': str(1/(1-b)),
    'actual_bulk_update_contraction_cap': str(beta),
    'retained_actual_bulk_update_contraction_cap': str(betaretained),
    'initial_free_wave_X_error_cap': str(einitial),
    'bulk_residual_response_X_cap': str(R),
    'first_actual_bulk_correction_X_cap': str(D0),
    'initial_data_limiting_upper_majorant': str(initial_majorant),
    'retained_initial_data_limiting_upper_majorant': str(initial_majorant_retained),
    'iteration_table': table,
    'one_update_actual_target_error_X_cap': str(E1),
    'one_update_retained_target_error_X_cap': str(E1retained),
    'one_update_readable_error_strict_upper': str(readable),
    'previous_to_one_update_readable_certificate_ratio': str(E0retained/readable),
    'norm_only_target_enclosure_gain': False,
    'norm_only_target_enclosure_identity': 'D0+E1=E0: unresolved correction widening plus the new error reproduces the old target norm radius',
    'candidate_scope': 'W1=P[Bs+K(La W0+F)+LG] is integral-defined with fixed actual coefficients and actual sources; the actual bulk correction has not been quadrature-evaluated here',
    'boundary_identity': 'each integral-defined update retains exact actual affine Robin and homogeneous Dirichlet laws',
    'source_units': 'R bounds X(K r0); SOURCE bounds X(KF+LG); no pointwise F or reduced residual evaluation is asserted',
    'limitations': ['The initial-data limiting majorant is not a lower bound on true error.',
        'A finite actual-bulk update is not a full self-consistent bulk solution.',
        'Inherited initial/operator/source/zero-wedge and formal completion premises remain conditional.',
        'No numerical C2 error, pointwise Einstein/constraint residual, nonlinear family or physical observable.']
}, sort_keys=True))
