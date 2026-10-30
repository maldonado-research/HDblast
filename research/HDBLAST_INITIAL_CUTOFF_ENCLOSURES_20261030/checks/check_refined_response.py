#!/usr/bin/env python3
"""Transfer a same-candidate initial-cell refinement into the Robin C1 bound.

This JSON-only audit verifies exact dyadic coverage and max/Gram aggregation,
then recomputes the direct response ledger and all23 shell target intervals.
It changes no candidate, nominal profile, smoothing, operator or bulk source.
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


def interval(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(Q, row)
    need(lo <= hi, 'interval order')
    return lo, hi


def sqrtgrid(q):
    need(q >= 0, 'nonnegative squared dual norm')
    grid = 1 << 112
    n = math.isqrt(q.numerator*grid*grid//q.denominator)
    if n*n*q.denominator < q.numerator*grid*grid:
        n += 1
    result = Q(n, grid)
    need(result*result >= q, 'outward fractional square root')
    return result


def caps_array(value):
    need(len(value) == 3 and all(len(channel) == 3 for channel in value)
         and all(len(row) == 23 for channel in value for row in channel),
         'three-channel three-kind component-cap shape')
    out = [[[Q(z) for z in row] for row in channel] for channel in value]
    need(all(z >= 0 for channel in out for row in channel for z in row),
         'nonnegative component magnitude caps')
    return out


def aggregate(target, caps):
    for channel in range(3):
        for kind in range(3):
            for i in range(23):
                target[channel][kind][i] = max(target[channel][kind][i],
                                               caps[channel][kind][i])


def verify_refinement(refined, initial, labels, weights):
    """Independent coverage/max/norm audit; scientific leaf evaluation is replayed
    by the separately authenticated producer, not imported into this checker.
    """
    need(refined['status'] == 'PASS' and refined['input_sha256'] == {
        'initial_producer':'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1',
        'initial_report':'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'},
        'same authenticated initial producer and inherited report')
    need(refined['labels'] == labels and list(map(Q, refined['gram'])) == weights
         and refined['channels'] == initial['channels'], 'same registered initial space')
    need(refined['same_candidate'] is True and
         refined['candidate_background_coefficients_sha256'] ==
             initial['candidate_background_coefficients_sha256'] and
         refined['derivative_convention'] == initial['derivative_convention'] and
         refined['footprint'] == initial['footprint'],
         'candidate coefficients, derivative convention and footprint unchanged')
    oldrows, cover = initial['rows'], refined['coverage_cells']
    need(len(oldrows) == len(cover) == initial['continuous_error_cells'] == 192,
         'complete inherited-cell coverage roster')
    maxima = [[[Q(0)]*23 for _ in range(3)] for _ in range(3)]
    oldmaxima = [[[Q(0)]*23 for _ in range(3)] for _ in range(3)]
    leafcount = refinedcount = 0
    for k, (old, entry) in enumerate(zip(oldrows, cover)):
        left, right = interval(old['r'])
        need((left, right) == (Q(256+k, 256), Q(257+k, 256)),
             'inherited continuous footprint partition')
        need(entry['old_row_index'] == k and interval(entry['r']) == (left, right)
             and isinstance(entry['refined'], bool), 'coverage parent identity')
        parent = caps_array(old['channel_component_error_caps'])
        aggregate(oldmaxima, parent)
        leaves = entry['leaves']
        if not entry['refined']:
            need(leaves == [], 'untouched cell must retain inherited cap')
            aggregate(maxima, parent)
            continue
        refinedcount += 1
        need(len(leaves) > 0, 'refined parent has leaves')
        cursor, seen = left, set()
        ordered = sorted(leaves, key=lambda row: Q(row['r'][0]))
        for leaf in ordered:
            path = leaf['path']
            need(isinstance(path, str) and all(z in '01' for z in path)
                 and path not in seen, 'unique binary refinement path')
            seen.add(path)
            unit = (right-left) / (1 << len(path))
            position = int(path, 2) if path else 0
            a, b = interval(leaf['r'])
            need(a == left+position*unit and b == a+unit and a == cursor
                 and b > a and b <= right, 'exact dyadic leaf coverage without gaps/overlap')
            need(not (a < Q(65,64) < b) and not (a < Q(66,64) < b),
                 'leaf stays in one analytic cutoff branch')
            raw = caps_array(leaf['raw_channel_component_error_caps'])
            used = caps_array(leaf['channel_component_error_caps'])
            need(all(used[c][kind][i] == min(parent[c][kind][i], raw[c][kind][i])
                     for c in range(3) for kind in range(3) for i in range(23)),
                 'leaf cap is intersection of two valid bounds on the same quantity')
            aggregate(maxima, used)
            cursor = b
            leafcount += 1
        need(cursor == right, 'refined parent fully covered')
    need(oldmaxima == caps_array(initial['channel_component_error_caps']),
         'independent inherited component maxima')
    need(maxima == caps_array(refined['channel_component_error_caps']),
         'independent refined component maxima')
    Hcaps = [[Q(z) for z in row] for row in refined['channel_H_caps']]
    need(len(Hcaps) == 3 and all(len(row) == 3 for row in Hcaps), 'channel H shape')
    for c in range(3):
        for kind in range(3):
            square = sum(weights[i]*maxima[c][kind][i]**2 for i in range(23))
            need(Hcaps[c][kind] >= 0 and square <= Hcaps[c][kind]**2,
                 'positive Gram norm encloses continuous refined channel')
            need(Hcaps[c][kind] <= Q(initial['channel_H_caps'][c][kind]),
                 'refined channel does not enlarge inherited certificate')
    totals = [sum(Hcaps[c][kind] for c in range(3)) for kind in range(3)]
    need(totals == list(map(Q, refined['total_H_caps']))
         and all(x <= Q(y) for x,y in zip(totals, initial['total_H_caps'])),
         'independent three-channel triangle aggregation')
    return totals, {'inherited_cells': len(oldrows), 'refined_parent_cells': refinedcount,
                    'verified_refined_leaves': leafcount,
                    'component_caps_aggregated_per_leaf': 207}


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--transfer-report', type=Path)
ap.add_argument('--source-report', type=Path)
ap.add_argument('--refined-report', type=Path)
ap.add_argument('--robin-capsule', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
initial = read(base/'initial/results/INITIAL_REPORT.json',
               'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
profiles = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
                '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = read(base/'compute/response-normal/RESPONSE_REPORT.json',
                'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
transfer_path = args.transfer_report or robin/'reference/TRANSFER.json'
transfer = read(transfer_path,
                '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52')
source_path = args.source_report or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030/reference/AFFINE_SOURCE.json'
source = read(source_path,
              'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e')
refined_path = args.refined_report or Path(__file__).resolve().parents[1]/'reference/REFINED.json'
REFINED_PIN = 'd97f7393cb9467e78387368754981e81cf763da043cab0a4df65978f52b72122'
refined = read(refined_path, REFINED_PIN)
old_direct = read(robin/'reference/DIRECT.json',
                  '4506b4104f6d80e7b46b527b968ab7046765f6235cf605a1f4643d799f884b54')
old_sample = read(robin/'reference/RESUMMED.json',
                  'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69')
old_posterior = read(robin/'reference/POSTERIOR.json',
                     '41857abcc6f77d4fc09c986e0b39a0b5836ed4d16591fc75306348cb40621012')
review_path = root/'research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md'
pin(review_path, '69a7a4741d599aad07d952a9b48423bf9b06f18c2b2431f3de8e821134db610c')
need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER'
     and response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR'
     and transfer['status'] == source['status'] == 'PASS', 'inherited scope/status')
need(all(z['status'] == 'PASS' for z in (old_direct, old_sample, old_posterior)),
     'same-candidate direct and sample reference reports')
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
new_totals, coverage_audit = verify_refinement(refined, initial, labels, w)
new_du, new_dur, new_dv = new_totals
need(new_dur < dur, 'strictly improved same-candidate initial radial certificate')
initial_errors = [new_du+au, new_dur+bu, new_dv+ev+av]
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
unchanged = {
    'smooth_B0_X_cap': B0cap,
    'boundary_correction_whole_cone_X_cap': Ccap,
    'frozen_coupled_X_operator_cap': q0,
    'actual_coupled_X_operator_cap': qactual,
    'retained_coupled_X_operator_cap': qretained,
    'frozen_first_feedback_whole_cone_X_cap': B1cap,
    'W_minus_B0_whole_cone_X_cap': WminusB0cap,
    'lower_order_coefficient_residual_H_cap': coefficient_residual,
    'lower_order_feedback_residual_H_cap': feedback_residual,
    'lower_order_bulk_residual_H_cap_excluding_affine_F': residualH,
    'affine_bulk_response_KF_X_cap': KFcap,
    'bulk_residual_response_X_cap': forcing_residual_X}
need(all(value == Q(old_direct[key]) for key,value in unchanged.items()),
     'candidate, source and every bulk/norm residual contribution remain unchanged')
direct_residual_X = einitial+forcing_residual_X
E = direct_residual_X/(1-qactual)
Eretained = direct_residual_X/(1-qretained)
Eprior = Q(transfer['inherited_prior_exact_C1_X_error'])
Eprior_direct = Q(old_direct['direct_C1_X_error_using_retained_contraction_cap'])
readable_step = Q(10**22)
Ereadable = (Eretained//readable_step+1)*readable_step
need(E < Eretained < Ereadable < Eprior_direct < Q(1735*10**24),
     'strictly improved conservative response certificate')
need(E >= einitial and residualH > 0, 'initial error and nonzero bulk uncertainty retained')
need(old_sample['point'] == old_posterior['point'] == ['9/8','0']
     and old_sample['labels'] == labels and
     len(old_sample['full23_point_table']) == len(old_posterior['full23_point_table']) == 23,
     'same retained reconstructed sample')
point_rows=[]
point_checks=0
for i,(sample, prior) in enumerate(zip(old_sample['full23_point_table'],
                                      old_posterior['full23_point_table'])):
    need(sample['label'] == prior['label'] == labels[i] and
         sample['parity'] == prior['parity'], 'same registered sample row')
    factor = sqrtgrid(1/w[i])
    candidate = {key:interval(sample['reconstructed_boundary_candidate'][key])
                 for key in ('value','dx','dr')}
    need(all(candidate[key] == interval(prior['reconstructed_candidate'][key])
             for key in candidate), 'unchanged reconstructed candidate point')
    actual_point={key:(value[0]-factor*Ereadable*(Q(1,8) if key=='value' else 1),
                 value[1]+factor*Ereadable*(Q(1,8) if key=='value' else 1))
            for key,value in candidate.items()}
    if sample['parity'] < 0:
        need(candidate['value'] == candidate['dx'] == (0,0), 'exact Dirichlet trace')
        actual_point['value'] = actual_point['dx'] = (Q(0),Q(0))
        shared_radial = actual_point['dr']
        radial_error = None
    else:
        rowboxes=[interval(z) for z in response['actual_coefficient_boxes']['M'][i]]
        rownorm2=sum(max(abs(lo),abs(hi))**2/w[j]
                     for j,(lo,hi) in enumerate(rowboxes))
        radial_error=sqrtgrid(rownorm2)*Ereadable/8
        boundary=(candidate['dr'][0]-radial_error,candidate['dr'][1]+radial_error)
        shared_radial=(max(actual_point['dr'][0],boundary[0]),min(actual_point['dr'][1],boundary[1]))
        need(shared_radial[0] <= shared_radial[1], 'compatible shared Robin radial interval')
    for key,interval_new in actual_point.items():
        interval_old=interval(prior['conditional_actual_formal_response'][key])
        need(interval_old[0] <= interval_new[0] <= interval_new[1] <= interval_old[1],
             'refined response interval contained in old conditional interval')
        point_checks+=1
    oldshared=interval(prior['shared_Robin_actual_radial'])
    need(oldshared[0] <= shared_radial[0] <= shared_radial[1] <= oldshared[1],
         'shared Robin refined interval contained in old conditional interval')
    point_checks+=1
    if i in (21,22):need(radial_error==0, 'zero-M tensor radial error remains zero')
    point_rows.append({'label':labels[i],'parity':sample['parity'],
        'reconstructed_candidate':{key:list(map(str,value)) for key,value in candidate.items()},
        'conditional_actual_formal_response':{key:list(map(str,value)) for key,value in actual_point.items()},
        'shared_Robin_actual_radial':list(map(str,shared_radial)),
        'shared_Robin_radial_error_abs_upper':str(radial_error) if radial_error is not None else None})
mean=point_rows[12]
need(mean['label']=='w:0' and all(Q(mean['conditional_actual_formal_response'][key][0])>0
                                 for key in ('value','dx','dr'))
     and Q(mean['shared_Robin_actual_radial'][0])>0, 'positive refined formal scalar-mean first jet')
print(json.dumps({
    'status':'PASS',
    'method':'independent exact coverage/Gram audit and same-candidate direct C1 transfer',
    'refined_initial_report_sha256':REFINED_PIN,
    'coverage_audit':coverage_audit,
    'refined_actual_minus_nominal_initial_H_caps':list(map(str,new_totals)),
    'inherited_actual_minus_nominal_initial_H_caps':initial['total_H_caps'],
    'initial_error_certificate_reduction_ratios':[str(Q(a)/b) for a,b in zip(initial['total_H_caps'],new_totals)],
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
    'readable_direct_C1_X_error_strict_upper':str(Ereadable),
    'previous_PR95_direct_C1_X_error_using_retained_contraction_cap':str(Eprior_direct),
    'PR95_to_refined_error_certificate_ratio':str(Eprior_direct/Eretained),
    'PR95_to_refined_readable_error_certificate_ratio':str(Eprior_direct/Ereadable),
    'initial_term_fraction_of_direct_numerator':str(einitial/direct_residual_X),
    'inherited_prior_C1_X_error':str(Eprior),
    'prior_to_new_certificate_ratio_lower':str(Eprior/Eretained),
    'unchanged_bulk_and_operator_quantities_verified':list(unchanged),
    'point':['9/8','0'],
    'full23_point_table':point_rows,
    'refined_point_interval_containment_checks':point_checks,
    'full_boundary_defect':'exactly zero for the integral-defined reconstruction, by analytic Robin trace equation; not an interval sum identity',
    'source_norm_convention':'SOURCE bounds K F+L G in X; X(K F)<=SOURCE+Gcap; SOURCE is not a pointwise F cap',
    'dependency_transfer_report_sha256':'7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52',
    'scope':'same candidate W and target U; only inherited initial-cell enclosures refined; conditional C1 upper certificate and modal point intervals, not measured error, numerical C2 error, pointwise Einstein residual, nonlinear family or physical observable'
},sort_keys=True))
