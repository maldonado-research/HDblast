#!/usr/bin/env python3
"""Independent retained-majorant audit for one actual-bulk Robin update.

No new bulk-update or point-source implementation is imported. All numbers
are exact rational outer bounds from authenticated reports and recomputed
Gram operator inequalities. Candidate point boxes use an unresolved norm
bound for the correction; they are not evaluations of the bulk convolution.
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


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def interval(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(Q, row)
    need(lo <= hi, 'interval order')
    return lo, hi


def sqrt_integer(q):
    need(q >= 0, 'nonnegative squared norm')
    k = math.isqrt(q.numerator // q.denominator)
    if k*k < q:
        k += 1
    need(k*k >= q, 'outward integer root')
    return Q(k)


def sqrt_grid(q):
    need(q >= 0, 'nonnegative squared row norm')
    grid = 1 << 112
    k = math.isqrt(q.numerator*grid*grid // q.denominator)
    if k*k*q.denominator < q.numerator*grid*grid:
        k += 1
    value = Q(k, grid)
    need(value*value >= q, 'outward rational root')
    return value


def widen(row, radius):
    need(radius >= 0, 'nonnegative radius')
    return row[0]-radius, row[1]+radius


def enc(row):
    return list(map(str, row))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--initial-capsule', type=Path)
ap.add_argument('--robin-capsule', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
initial_cap = args.initial_capsule or root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030'
robin_cap = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
prior = read(initial_cap/'reference/RESPONSE.json',
             'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84')
refined = read(initial_cap/'reference/REFINED.json',
               'd97f7393cb9467e78387368754981e81cf763da043cab0a4df65978f52b72122')
point = read(robin_cap/'reference/RESUMMED.json',
             'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69')
transfer = read(robin_cap/'reference/TRANSFER.json',
                '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52')
profile = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
               '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = read(base/'compute/response-normal/RESPONSE_REPORT.json',
                'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
initial = read(base/'initial/results/INITIAL_REPORT.json',
               'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
need(all(d['status'] == 'PASS' for d in (prior, refined, point, transfer)),
     'successful prior certificates')
need(response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR' and
     initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER',
     'inherited formal-response scope')
need(prior['point'] == point['point'] == ['9/8', '0'], 'same shell sample')
labels = profile['labels']
w = list(map(Q, profile['gram']))
need(len(labels) == len(w) == 23 and min(w) > 0 and
     labels == [row['label'] for row in prior['full23_point_table']],
     'positive registered dimension and prior row identities')
need(labels == point['labels'] == refined['labels'] == response['labels'] == initial['labels']
     and w == list(map(Q, response['gram'])) == list(map(Q, initial['gram'])),
     'common full23 space')
need(response['domain'] == {'r_min': 'max(0,x-5/4)', 'x': ['0','3/2'],
                           'x_plus_r_max': '7/4'}, 'retained cone')
need(initial['footprint']['actual_zero_on'] == ['0','1'] and
     refined['footprint'] == initial['footprint'], 'same initial zero collar')
parity = [-1 if z.startswith(('L:', 'v:')) else 1 for z in labels]
need(profile['parity']['u'] == profile['parity']['v'] == parity and
     profile['parity']['du'] == [-p for p in parity], 'reflection and derivatives')
blocks = [[i for i, name in enumerate(labels) if name.endswith(':'+suffix)]
          for suffix in ('0','2','4','+','-')]
need(sum(map(len, blocks)) == 23, 'complete orthogonal blocks')
sameblock = lambda i, j: any(i in b and j in b for b in blocks)
op = {key: Q(response['constant_matrix_H_operator_upper'][key])
      for key in ('A','B','D','M')}
err = {key: Q(response['coefficient_error_H_operator_upper'][key]) for key in op}
actual = {key: op[key]+err[key] for key in op}
mat = {key: [list(map(Q, row)) for row in profile['matrices'][key]] for key in op}
actual_M = None
for key in op:
    em = [list(map(Q, row)) for row in response['coefficient_entry_error_caps'][key]]
    boxes = [[interval(z) for z in row] for row in response['actual_coefficient_boxes'][key]]
    need(len(mat[key]) == len(em) == len(boxes) == 23 and
         all(len(row) == 23 for matrix in (mat[key], em, boxes) for row in matrix),
         'matrix shapes')
    need(profile['matrices'][key] == response['constant_matrices'][key], 'frozen matrix agreement')
    for i in range(23):
        for j in range(23):
            lo, hi = boxes[i][j]
            need(em[i][j] >= max(abs(lo-mat[key][i][j]), abs(hi-mat[key][i][j])),
                 'actual coefficient error entry')
            need(sameblock(i,j) or mat[key][i][j] == em[i][j] == 0,
                 'orthogonal block closure')
            if key == 'M' and (parity[i] < 0 or parity[j] < 0):
                need((lo,hi) == (0,0) and mat[key][i][j] == 0,
                     'actual Robin rows and columns')
    for block in blocks:
        need(sum(w[i]*mat[key][i][j]**2/w[j] for i in block for j in block) <= op[key]**2,
             'frozen Gram operator cap')
        need(sum(w[i]*em[i][j]**2/w[j] for i in block for j in block) <= err[key]**2,
             'actual error Gram operator cap')
    need(actual[key] == Q(prior['actual_coefficient_H_norm_caps'][key]),
         'same actual operator norm')
    if key == 'M':
        actual_M = boxes

# Independently reconstruct the initial free-wave error.
du, dur, dv = map(Q, refined['total_H_caps'])
need([du,dur,dv] == list(map(Q, prior['refined_actual_minus_nominal_initial_H_caps'])),
     'same continuous initial certificate')
nodes = list(map(Q, profile['nodes']))
vtab = [list(map(Q, row)) for row in profile['profiles']['v']]
slopes = [[(vtab[k+1][i]-vtab[k][i])/(nodes[k+1]-nodes[k]) for i in range(23)]
          for k in range(len(nodes)-1)]
jumps = [[b-a for a,b in zip(left,right)] for left,right in zip(slopes,slopes[1:])]
av = Q(35,256*2048)*sqrt_integer(max(sum(w[i]*row[i]**2 for i in range(23)) for row in jumps))
M2 = Q(initial['nominal_second_derivative_H_caps'][0])
au = M2*Q(1,2**30)**2/18
bu = Q(35,128)*M2*Q(1,2**30)
ev = Q(initial['P1_H_error_caps'][2])
einitial = max((du+au+Q(3,2)*(dv+ev+av))/3, dur+bu+dv+ev+av)
need(einitial == Q(prior['initial_free_wave_X_error_cap']), 'independent initial free-wave error')

# Rebuild the previous bulk residual from norms independently. The affine
# response bound is in X units and must not be used as a raw pointwise F cap.
nu, nur, nv = map(Q, initial['nominal_H_caps'])
B0cap = max((nu+au+Q(3,2)*(nv+ev+av))/3, nur+bu+nv+ev+av)
history = transfer['global_retained_history']
need(Q(history['shell_history_end']) == Q(5,4) and
     Q(history['active_history_length']) == Q(1,4), 'global boundary history')
Gcap = Q(prior['actual_affine_boundary_G_H_cap'])
need(Gcap == Q(transfer['affine_shell_G_H_cap']) == 487122,
     'same inherited uniform affine shell cap')
boundary_residual = (op['M']*Q(history['B1_shell_trace_H_cap'])+
                     err['M']*Q(history['V_shell_trace_H_cap'])+Gcap)
Cglobal = boundary_residual/(1-actual['M']/4)
q0 = Q(3,2)*(op['A']+op['B']+Q(3,4)*op['D'])+op['M']/4
feedback = q0*B0cap+Cglobal
coefficient_residual = (err['A']+err['B']+Q(3,4)*err['D'])*B0cap
feedback_residual = (actual['A']+actual['B']+Q(3,4)*actual['D'])*feedback
SOURCE = Q(response['ledger']['omitted_full_H_forcing_X_upper'])
need(SOURCE == Q(52162547172519,2), 'same inherited affine response bound')
R = Q(3,2)*(coefficient_residual+feedback_residual)+SOURCE+Gcap
need(B0cap == Q(prior['smooth_B0_X_cap']) and
     Cglobal == Q(prior['boundary_correction_whole_cone_X_cap']) and
     R == Q(prior['bulk_residual_response_X_cap']), 'independent bulk residual arithmetic')
need(R > 0 and prior['source_norm_convention'] ==
     'SOURCE bounds K F+L G in X; X(K F)<=SOURCE+Gcap; SOURCE is not a pointwise F cap',
     'residual response source units')
need(Q(prior['affine_bulk_response_KF_X_cap']) ==
     Q(prior['inherited_affine_response_S_X_cap']) + Q(prior['actual_affine_boundary_G_H_cap']),
     'affine source response triangle')
k = Q(3,2)*(actual['A']+actual['B']+Q(3,4)*actual['D'])
b = actual['M']/4
qretained = Q(response['ledger']['retained_full_H_contraction'])
need(0 < b and 0 < k and k+b <= qretained == Q(3409,8000000) < 1,
     'raw coefficient caps fit retained contraction')
beta_actual = k/(1-b)
# This intentionally uses the older, larger retained majorant independently
# of the sharper bulk-update checker: k <= qretained-b is verified above.
beta = (qretained-b)/(1-b)
need(0 < beta_actual <= beta < 1, 'boundary-preconditioned contraction')
Eold = (einitial+R)/(1-qretained)
need(Eold == Q(prior['direct_C1_X_error_using_retained_contraction_cap']),
     'independent old residual certificate')
D = R/(1-b)
floor = einitial/(1-qretained)
Enew = einitial/(1-b)+beta*Eold
need(Enew == floor+beta*D/(1-beta) == (einitial+beta*R)/(1-qretained),
     'independent recurrence and floor-tail forms agree')
need(Enew < Q(3541*10**22) and Eold/Enew > Q(57,10),
     'readable one-update conditional error gate')
need(D+Enew == Eold, 'unresolved norm correction exhausts the old error radius')

oldrows = point['full23_point_table']
need(len(oldrows) == 23, 'full prior candidate sample')
rows = []
for i, row in enumerate(oldrows):
    need(row['label'] == labels[i] and row['parity'] == parity[i], 'sample row identity')
    factor = sqrt_grid(1/w[i])
    W0 = {key: interval(row['reconstructed_boundary_candidate'][key])
          for key in ('value','dx','dr')}
    rownorm = sqrt_grid(sum(max(abs(lo),abs(hi))**2/w[j]
                           for j,(lo,hi) in enumerate(actual_M[i])))
    if parity[i] > 0:
        radcorr, raderror = rownorm*D/8, rownorm*Enew/8
        W1 = {'value': widen(W0['value'], factor*D/8),
              'dx': widen(W0['dx'], factor*D),
              'dr': widen(W0['dr'], radcorr)}
        shared_actual = widen(W1['dr'], raderror)
        same_old_radial = widen(W0['dr'], rownorm*Eold/8)
        need(shared_actual == same_old_radial, 'shared Robin target-width identity')
    else:
        radcorr, raderror = factor*D, factor*Enew
        W1 = {'value': (Q(0),Q(0)), 'dx': (Q(0),Q(0)),
              'dr': widen(W0['dr'], radcorr)}
        shared_actual = widen(W1['dr'], raderror)
    Ubox = {key: widen(z, factor*Enew*(Q(1,8) if key == 'value' else 1))
            for key,z in W1.items()}
    if parity[i] < 0:
        Ubox['value'] = Ubox['dx'] = (Q(0),Q(0))
    for key in ('value','dx'):
        if parity[i] > 0:
            need(Ubox[key] == widen(W0[key], factor*Eold*(Q(1,8) if key=='value' else 1)),
                 'norm-only target point width unchanged')
    if i in (21,22):
        need(rownorm == 0 and W1['dr'] == W0['dr'], 'zero-M tensor radial correction')
    rows.append({'label':labels[i], 'parity':parity[i],
        'W0':{key:enc(z) for key,z in W0.items()},
        'norm_enclosed_actual_bulk_update_W1':{key:enc(z) for key,z in W1.items()},
        'norm_only_target_U_enclosure':{key:enc(z) for key,z in Ubox.items()},
        'shared_Robin_target_radial':enc(shared_actual),
        'W1_minus_W0_radial_correction_abs_upper':str(radcorr),
        'U_minus_W1_shared_Robin_radial_error_abs_upper':str(raderror)})
need(all(interval(rows[12]['norm_enclosed_actual_bulk_update_W1'][key])[0] > 0
         for key in ('value','dx','dr')), 'positive reconstructed mean first-jet boxes')
print(json.dumps({'status':'PASS',
    'method':'independent retained-contraction majorant; no new bulk evaluator imported',
    'point':['9/8','0'], 'actual_bulk_KLa_X_operator_cap':str(k),
    'actual_Robin_boundary_X_operator_cap':str(b),
    'retained_total_X_operator_cap':str(qretained),
    'actual_boundary_preconditioned_bulk_contraction_cap':str(beta_actual),
    'independent_retained_boundary_preconditioned_contraction_cap':str(beta),
    'initial_free_wave_error_X_cap':str(einitial),
    'prior_bulk_residual_response_X_cap':str(R),
    'prior_actual_target_error_X_cap':str(Eold),
    'first_actual_bulk_correction_X_cap':str(D),
    'one_update_actual_target_error_X_cap':str(Enew),
    'readable_one_update_error_strict_upper':str(Q(3541*10**22)),
    'asymptotic_initial_data_upper_certificate_X_cap':str(floor),
    'remaining_smooth_data_fixed_point_tail_after_one_update_X_cap':str(beta*D/(1-beta)),
    'prior_to_one_update_upper_certificate_ratio':str(Eold/Enew),
    'norm_width_identity_D_plus_Enew_equals_Eold':True,
    'full23_point_table':rows,
    'numerical_resolution_limit':'bulk correction is enclosed by a global norm majorant, not evaluated as a pointwise causal convolution; target value/time and best shared-Robin radial boxes retain the preceding norm resolution; a generic radial box may shrink but remains weaker than the shared-Robin box',
    'scope':'conditional actual-coefficient/affine-source integral bulk update with exact inherited Robin law; first update not a converged coupled solution, numerical C2 error, Einstein/constraint residual, nonlinear spacetime or physical observable'}, sort_keys=True))
