#!/usr/bin/env python3
"""Direct C1 residual bound for the integral-defined Robin reconstruction.

This is a JSON-only arithmetic check of a conditional analytic residual bound.
It does not import the new point evaluator or use the old P1 bulk forcing error.
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


def sqrtup(q):
    need(q >= 0, 'nonnegative norm square')
    k = math.isqrt(q.numerator // q.denominator)
    if k*k < q:
        k += 1
    need(k*k >= q, 'outward square root')
    return Q(k)


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--transfer-report', type=Path)
ap.add_argument('--source-report', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
initial = read(base/'initial/results/INITIAL_REPORT.json',
               'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
profiles = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
                '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = read(base/'compute/response-normal/RESPONSE_REPORT.json',
                'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
transfer_path = args.transfer_report or Path(__file__).resolve().parents[1]/'reference/TRANSFER.json'
transfer = read(transfer_path,
                '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52')
source_path = args.source_report or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030/reference/AFFINE_SOURCE.json'
source = read(source_path,
              'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e')
review_path = root/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md'
pin(review_path, '69a7a4741d599aad07d952a9b48423bf9b06f18c2b2431f3de8e821134db610c')
need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER'
     and response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR'
     and transfer['status'] == source['status'] == 'PASS', 'inherited scope/status')
labels = profiles['labels']
w = list(map(Q, profiles['gram']))
need(len(labels) == len(w) == 23 and min(w) > 0 and
     labels == response['labels'] == initial['labels'] == transfer['labels'] and
     w == list(map(Q, response['gram'])) == list(map(Q, initial['gram'])),
     'common full23 Hilbert space')
need(response['domain'] == {'r_min': 'max(0,x-5/4)', 'x': ['0','3/2'],
                            'x_plus_r_max': '7/4'}, 'retained cone')
need(initial['footprint']['actual_zero_on'] == ['0', '1'] and
     initial['footprint']['nominal_zero_on'] == ['0', '65/64'] and
     initial['footprint']['r'] == ['0', '7/4'], 'inherited initial zero collars')
need('U and V are C1' in initial['regularity'] and
     initial['derivative_convention']['second_derivatives'][0] == 'partial_r^2 U',
     'nominal C1,1 convolution premise')
parity = [-1 if name.startswith(('L:', 'v:')) else 1 for name in labels]
need(profiles['parity']['u'] == profiles['parity']['v'] == parity and
     profiles['parity']['du'] == [-s for s in parity], 'derivative-consistent reflection')
blocks = [[i for i, label in enumerate(labels) if label.endswith(':' + degree)]
          for degree in ('0', '2', '4', '+', '-')]
need(sum(map(len, blocks)) == 23, 'complete block decomposition')
sameblock = lambda i, j: any(i in b and j in b for b in blocks)
op = {key: Q(response['constant_matrix_H_operator_upper'][key])
      for key in ('A', 'B', 'D', 'M')}
error = {key: Q(response['coefficient_error_H_operator_upper'][key]) for key in op}
for key in op:
    frozen = [list(map(Q, row)) for row in profiles['matrices'][key]]
    entry = [list(map(Q, row)) for row in response['coefficient_entry_error_caps'][key]]
    boxes = [[tuple(map(Q, z)) for z in row]
             for row in response['actual_coefficient_boxes'][key]]
    need(profiles['matrices'][key] == response['constant_matrices'][key],
         'frozen matrix agreement')
    need(len(frozen) == len(entry) == len(boxes) == 23 and
         all(len(row) == 23 for matrix in (frozen, entry, boxes) for row in matrix),
         'matrix shapes')
    for i in range(23):
        for j in range(23):
            lo, hi = boxes[i][j]
            need(lo <= hi and entry[i][j] >= 0 and
                 max(abs(lo-frozen[i][j]), abs(hi-frozen[i][j])) <= entry[i][j],
                 'actual coefficient entry error cap')
            need(sameblock(i, j) or frozen[i][j] == entry[i][j] == 0,
                 'off-block closure')
            if key == 'M' and (parity[i] < 0 or parity[j] < 0):
                need(frozen[i][j] == 0 and (lo, hi) == (0,0),
                     'actual Robin split')
    for block in blocks:
        need(sum(w[i]*frozen[i][j]**2/w[j] for i in block for j in block)
             <= op[key]**2, 'frozen H operator cap')
        need(sum(w[i]*entry[i][j]**2/w[j] for i in block for j in block)
             <= error[key]**2, 'actual coefficient error H cap')
actual = {key: op[key] + error[key] for key in op}

# Reconstruct all smoothing and initial/free-wave error contributions directly.
nodes = list(map(Q, profiles['nodes']))
vtab = [list(map(Q, row)) for row in profiles['profiles']['v']]
delta_u, delta_v = Q(1,2**30), Q(1,2048)
need(nodes[0] == 0 and nodes[-1] == Q(7,4) and
     min(b-a for a,b in zip(nodes,nodes[1:])) >= 2*delta_v,
     'registered mesh and disjoint smoothing bands')
need(Q(65,64)-delta_u > 1 and nodes[1]-delta_v > 1 and
     all(z == 0 for row in vtab[:2] for z in row), 'smoothed zero collar')
slopes = [[(vtab[k+1][i]-vtab[k][i])/(nodes[k+1]-nodes[k]) for i in range(23)]
          for k in range(len(nodes)-1)]
jumps = [[b-a for a,b in zip(left,right)] for left,right in zip(slopes,slopes[1:])]
jumpnorm2 = max(sum(w[i]*row[i]**2 for i in range(23)) for row in jumps)
av = Q(35,256)*delta_v*sqrtup(jumpnorm2)
M2H = Q(initial['nominal_second_derivative_H_caps'][0])
au = M2H*delta_u**2/18
bu = Q(35,128)*M2H*delta_u
nu, nur, nv = map(Q, initial['nominal_H_caps'])
eu, eur, ev = map(Q, initial['P1_H_error_caps'])
du, dur, dv = map(Q, initial['total_H_caps'])
need([du,dur,dv] == [Q(response['ledger'][name]) for name in
                     ('delta_u_H','delta_du_H','delta_v_H')], 'actual initial difference caps')
need([eu,eur,ev] == [Q(response['ledger'][name]) for name in
                     ('epsilon_u_H','epsilon_du_H','epsilon_v_H')], 'P1 difference caps')
initial_errors = [du+au, dur+bu, dv+ev+av]
einitial = max((initial_errors[0]+Q(3,2)*initial_errors[2])/3,
               initial_errors[1]+initial_errors[2])
B0cap = max((nu+au+Q(3,2)*(nv+ev+av))/3, nur+bu+nv+ev+av)

# Independently recheck the new authenticated boundary-correction majorant.
history = transfer['global_retained_history']
need(Q(history['shell_history_end']) == Q(5,4) and
     Q(history['active_history_length']) == Q(1,4), 'whole-domain boundary history')
Gcap = Q(source['source_H_norm_strict_integer_upper'])
need(source['uniform_domain'] == '0<=x<=3,r=0' and Gcap <= 487122 and
     Q(source['source_H_norm_square_upper']) < Gcap**2, 'actual source norm bound')
B1trace, Vtrace = (Q(history[key]) for key in ('B1_shell_trace_H_cap','V_shell_trace_H_cap'))
Rboundary = op['M']*B1trace + error['M']*Vtrace + Gcap
Ccap = Rboundary/(1-actual['M']/4)
need(actual['M']/4 < 1 and Ccap == Q(history['correction_C1_X_cap']),
     'independent exact boundary majorant arithmetic')

# Compatible zero-wedge bounds are |B0|<=3 X(B0)/4 and
# |B0|shell<=X(B0)/4. Thus this is the WHOLE-CONE B1 cap, not a shell cap.
q0 = Q(3,2)*(op['A']+op['B']+Q(3,4)*op['D'])+op['M']/4
qactual = Q(3,2)*(actual['A']+actual['B']+Q(3,4)*actual['D'])+actual['M']/4
qretained = Q(response['ledger']['retained_full_H_contraction'])
need(0 <= q0 < qactual <= qretained == Q(3409,8000000) < 1,
     'verified actual coupled contraction')
B1cap = q0*B0cap
WminusB0cap = B1cap+Ccap
coefficient_residual = (error['A']+error['B']+Q(3,4)*error['D'])*B0cap
feedback_residual = (actual['A']+actual['B']+Q(3,4)*actual['D'])*WminusB0cap
residualH = coefficient_residual+feedback_residual

# SOURCE is an X bound for S=K F+L G, NOT a pointwise F cap.
# Triangle inequality gives X(K F)<=X(S)+X(L G)<=SOURCE+Gcap.
SOURCE = Q(response['ledger']['omitted_full_H_forcing_X_upper'])
need(SOURCE == Q(52162547172519,2), 'inherited full-H affine response source bound')
KFcap = SOURCE+Gcap
forcing_residual_X = Q(3,2)*residualH+KFcap
direct_residual_X = einitial+forcing_residual_X
E = direct_residual_X/(1-qactual)
Eretained = direct_residual_X/(1-qretained)
Eprior = Q(transfer['inherited_prior_exact_C1_X_error'])
need(E < Eretained < Q(1735*10**24) < 10**29, 'readable direct C1 target')
need(Eretained < Eprior/12, 'strict twelvefold improvement in the upper certificate')
need(E >= einitial and residualH > 0, 'initial error and nonzero bulk uncertainty retained')
print(json.dumps({
    'status':'PASS',
    'method':'direct initial and coupled-bulk residual estimate for exact actual-Robin auxiliary reconstruction',
    'initial_smoothed_data_H_error_caps':list(map(str,initial_errors)),
    'initial_free_wave_X_error_cap':str(einitial),
    'displacement_smoothing_value_H_error':str(au),
    'displacement_smoothing_derivative_H_error':str(bu),
    'velocity_smoothing_H_error':str(av),
    'smooth_B0_X_cap':str(B0cap),
    'frozen_first_feedback_whole_cone_X_cap':str(B1cap),
    'boundary_correction_whole_cone_X_cap':str(Ccap),
    'W_minus_B0_whole_cone_X_cap':str(WminusB0cap),
    'actual_coefficient_H_norm_caps':{k:str(v) for k,v in actual.items()},
    'actual_coefficient_error_H_norm_caps':{k:str(v) for k,v in error.items()},
    'frozen_coupled_X_operator_cap':str(q0),
    'actual_coupled_X_operator_cap':str(qactual),
    'retained_coupled_X_operator_cap':str(qretained),
    'lower_order_coefficient_residual_H_cap':str(coefficient_residual),
    'lower_order_feedback_residual_H_cap':str(feedback_residual),
    'lower_order_bulk_residual_H_cap_excluding_affine_F':str(residualH),
    'inherited_affine_response_S_X_cap':str(SOURCE),
    'actual_affine_boundary_G_H_cap':str(Gcap),
    'affine_bulk_response_KF_X_cap':str(KFcap),
    'bulk_residual_response_X_cap':str(forcing_residual_X),
    'total_direct_residual_X_cap':str(direct_residual_X),
    'direct_C1_X_error_cap':str(E),
    'direct_C1_X_error_using_retained_contraction_cap':str(Eretained),
    'readable_direct_C1_X_error_strict_upper':'1735*10^24',
    'inherited_prior_C1_X_error':str(Eprior),
    'prior_to_new_certificate_ratio_lower':str(Eprior/Eretained),
    'strict_twelvefold_certificate_improvement':True,
    'full_boundary_defect':'exactly zero for the integral-defined reconstruction, by analytic Robin trace equation; not an interval sum identity',
    'source_norm_convention':'SOURCE bounds K F+L G in X; X(K F)<=SOURCE+Gcap; SOURCE is not a pointwise F cap',
    'dependency_transfer_report_sha256':'7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52',
    'scope':'conditional C1 error upper bound against the retained formal full23 response, using inherited initial/operator/source/zero-wedge premises and exact Robin reconstruction; not measured error, numerical C2 error, pointwise Einstein residual, nonlinear family or physical observable'
},sort_keys=True))
