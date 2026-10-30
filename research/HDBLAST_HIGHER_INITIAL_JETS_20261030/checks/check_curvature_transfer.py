#!/usr/bin/env python3
"""Charge a curvature-matched outer extension to a new finite candidate.

The old W1 is preserved. This propagates only the controlled change of the
smoothed nominal displacement; it does not certify second-derivative error.
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


INITIAL_PIN = 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207'
BULK_PIN = '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f'
PROFILE_PIN = '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'
RESPONSE_PIN = 'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd'
OLD_INITIAL_PIN = 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'
REFINED_RESPONSE_PIN = 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84'
ROBIN_TRANSFER_PIN = '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52'
SHELL_SOURCE_PIN = 'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e'


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return json.loads(raw)


def iv(z):
    need(len(z) == 2, 'interval shape')
    a, b = map(Q, z)
    need(a <= b, 'ordered interval')
    return a, b


def sqrt_up(q):
    need(q >= 0, 'nonnegative curvature norm square')
    n = math.isqrt(q.numerator//q.denominator)
    if n*n < q:
        n += 1
    need(n*n >= q, 'outward integer square root')
    return Q(n)


def encode(z):
    if isinstance(z, Q):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (tuple, list)):
        return [encode(v) for v in z]
    return z


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--initial-second-report', type=Path, required=True)
ap.add_argument('--initial-report-sha256')
ap.add_argument('--bulk-capsule', type=Path)
ap.add_argument('--initial-cutoff-capsule', type=Path)
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--contact-capsule', type=Path)
args = ap.parse_args()
digest = INITIAL_PIN or args.initial_report_sha256
need(digest is not None and (INITIAL_PIN is None or args.initial_report_sha256 in (None, INITIAL_PIN)), 'registered initial-second input digest')
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
bulk = args.bulk_capsule or root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
cutoff = args.initial_cutoff_capsule or root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030'
robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
contact = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
initial = pin(args.initial_second_report, digest)
old = pin(bulk/'reference/BULK_BOUND.json', BULK_PIN)
profile = pin(base/'compute/response-normal/RESPONSE_PROFILES.json', PROFILE_PIN)
response = pin(base/'compute/response-normal/RESPONSE_REPORT.json', RESPONSE_PIN)
oldinitial = pin(base/'initial/results/INITIAL_REPORT.json', OLD_INITIAL_PIN)
prior = pin(cutoff/'reference/RESPONSE.json', REFINED_RESPONSE_PIN)
robin_transfer = pin(robin/'reference/TRANSFER.json', ROBIN_TRANSFER_PIN)
shellsource = pin(contact/'reference/AFFINE_SOURCE.json', SHELL_SOURCE_PIN)
labels = profile['labels']
gram = list(map(Q, profile['gram']))
need(labels == initial['labels'] == old['labels'] == response['labels'] and
     gram == list(map(Q, initial['gram'])) == list(map(Q, response['gram'])) and
     len(labels) == len(gram) == 23 and min(gram) > 0, 'same complete registered full23 space')
need(initial['status'] == old['status'] == prior['status'] == robin_transfer['status'] == shellsource['status'] == 'PASS' and
     response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR' and
     oldinitial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER', 'successful authenticated inputs')
R, delta = Q(7, 4), Q(1, 2**30)
need(list(map(Q, initial['footprint'])) == [0, R] and Q(initial['fixed_displacement_kernel_width']) == delta,
     'same initial footprint and fixed displacement kernel width')
rows = initial['endpoint']['rows']
need(len(rows) == 23 and [r['label'] for r in rows] == labels, 'complete outer curvature rows')
kappa = [iv(row['nominal_left_u_rr']) for row in rows]
square = sum(w*max(abs(z) for z in c)**2 for w, c in zip(gram, kappa))
kappaH = sqrt_up(square)

# Exact half-kernel moments of 35/32*(1-s^2)^3 on s in [0,1].
def half_moment(j):
    return Q(35, 32)*sum(Q(c, j+2*k+1) for k, c in enumerate((1, -3, 3, -1)))
need(half_moment(0) == Q(1, 2) and half_moment(1) == Q(35, 256) and
     half_moment(2) == Q(1, 18), 'exact normalized half-kernel moments')
displacement_factor = delta**2/36
radial_derivative_factor = 35*delta/256
free_factor = max(displacement_factor/3, radial_derivative_factor)
need(free_factor == radial_derivative_factor, 'derivative controls fixed-width free-wave X change')

# Verify the recorded frozen homogeneous operator caps in the same H norm.
blocks = [[i for i, z in enumerate(labels) if z.endswith(':'+l)] for l in ('0', '2', '4', '+', '-')]
need(sorted(i for block in blocks for i in block) == list(range(23)), 'orthogonal frozen block partition')
which = {i: b for b, block in enumerate(blocks) for i in block}
frozen_caps = {}
block_checks = 0
for key in ('A', 'B', 'D', 'M'):
    need(profile['matrices'][key] == response['constant_matrices'][key], 'same frozen source matrix')
    matrix = [list(map(Q, row)) for row in profile['matrices'][key]]
    need(all(matrix[i][j] == 0 for i in range(23) for j in range(23) if which[i] != which[j]), 'frozen degree-block closure')
    cap = Q(response['constant_matrix_H_operator_upper'][key])
    for block in blocks:
        norm_square = sum(gram[i]/gram[j]*matrix[i][j]**2 for i in block for j in block)
        need(norm_square <= cap**2, 'frozen H operator cap')
        block_checks += 1
    frozen_caps[key] = cap
k0 = Q(3, 2)*(frozen_caps['A']+frozen_caps['B']+Q(3, 4)*frozen_caps['D'])
k = Q(old['bulk_K_La_X_operator_cap'])
b = Q(old['Robin_L_M_trace_X_operator_cap'])
actual = {key: Q(value) for key, value in old['actual_coefficient_H_norm_caps'].items()}
need(k == Q(3, 2)*(actual['A']+actual['B']+Q(3, 4)*actual['D']) and b == actual['M']/4,
     'same inherited actual operator and boundary inverse norm constants')
need(0 < b < 1 and k > 0 and k0 > 0 and k+b <= Q(old['retained_coupled_X_operator_cap']) < 1,
     'retained causal operator bounds')
Pcap = 1/(1-b)
W0factor = Pcap*(1+k0)*free_factor
W1factor = Pcap*(free_factor+k*W0factor)
need(W1factor == (1+k*(1+k0)/(1-b))/(1-b)*free_factor, 'ordered two-stage finite-update comparison')
oldE = Q(old['one_update_retained_target_error_X_cap'])
old_readable = Q(old['one_update_readable_error_strict_upper'])
need(0 < oldE < old_readable == Q('3.541e25'), 'retained old finite-candidate error certificate')

# Directly reprove the old numerical upper certificate for the new candidate.
# Unlike a triangle charge, this uses the unchanged global C1,1 majorants.
M2 = list(map(Q, oldinitial['nominal_second_derivative_component_caps'][0]))
M2H = Q(oldinitial['nominal_second_derivative_H_caps'][0])
need(oldinitial['labels'] == labels and list(map(Q, oldinitial['gram'])) == gram and
     len(M2) == 23 and all(max(map(abs, c)) <= m for c, m in zip(kappa, M2)),
     'quadratic extension curvature obeys every old component M2 cap')
need(sum(w*m*m for w, m in zip(gram, M2)) <= M2H*M2H,
     'same global H-Lipschitz derivative cap for both extensions')
au, bu = M2H*delta**2/18, Q(35,128)*M2H*delta
nodes = list(map(Q, profile['nodes']))
vtab = [list(map(Q, row)) for row in profile['profiles']['v']]
dv = Q(1,2048)
need(nodes[0] == 0 and nodes[-1] == R and nodes[1]-dv > 1 and
     min(b-a for a,b in zip(nodes,nodes[1:])) >= 2*dv and
     all(z == 0 for row in vtab[:2] for z in row), 'unchanged velocity and flat initial collar')
slopes = [[(right-left)/(b-a) for left,right in zip(vl,vr)]
          for a,b,vl,vr in zip(nodes,nodes[1:],vtab,vtab[1:])]
jumps = [[b-a for a,b in zip(left,right)] for left,right in zip(slopes,slopes[1:])]
av = Q(35,256)*dv*sqrt_up(max(sum(w*z*z for w,z in zip(gram,row)) for row in jumps))
nu,nur,nv = map(Q,oldinitial['nominal_H_caps'])
eu,eur,ev = map(Q,oldinitial['P1_H_error_caps'])
du,dur,dvv = map(Q,prior['refined_actual_minus_nominal_initial_H_caps'])
initial_errors = [du+au,dur+bu,dvv+ev+av]
einitial = max((initial_errors[0]+Q(3,2)*initial_errors[2])/3,initial_errors[1]+initial_errors[2])
Bs_cap = max((nu+au+Q(3,2)*(nv+ev+av))/3,nur+bu+nv+ev+av)
need(initial_errors == list(map(Q,prior['initial_smoothed_data_H_error_caps'])) and
     einitial == Q(old['initial_free_wave_X_error_cap']) == Q(prior['initial_free_wave_X_error_cap']) and
     Bs_cap == Q(prior['smooth_B0_X_cap']), 'same smoothing, actual-minus-interior-nominal, initial mismatch and free-wave size majorants')
need(au == Q(prior['displacement_smoothing_value_H_error']) and
     bu == Q(prior['displacement_smoothing_derivative_H_error']) and
     av == Q(prior['velocity_smoothing_H_error']), 'every smoothing contribution reconstructed')
history = robin_transfer['global_retained_history']
need(Q(history['shell_history_end']) == Q(5,4) and Q(history['active_history_length']) == Q(1,4) and
     Q(5,4)+max(delta,dv) < R-delta, 'entire shell-history initial footprint precedes the changed collar')
Gcap = Q(shellsource['source_H_norm_strict_integer_upper'])
need(shellsource['uniform_domain'] == '0<=x<=3,r=0' and Gcap == 487122, 'same fixed actual G cap')
errors = {key:Q(response['coefficient_error_H_operator_upper'][key]) for key in ('A','B','D','M')}
need(all(actual[key] == frozen_caps[key]+errors[key] for key in actual), 'same inherited actual coefficient caps')
C_boundary = (frozen_caps['M']*Q(history['B1_shell_trace_H_cap'])+
              errors['M']*Q(history['V_shell_trace_H_cap'])+Gcap)/(1-b)
need(C_boundary == Q(history['correction_C1_X_cap']) == Q(prior['boundary_correction_whole_cone_X_cap']), 'same unchanged boundary-history correction majorant')
q0 = k0+frozen_caps['M']/4
B1_cap = q0*Bs_cap
feedback_cap = B1_cap+C_boundary
Rcoef = (errors['A']+errors['B']+Q(3,4)*errors['D'])*Bs_cap
Rfeedback = (actual['A']+actual['B']+Q(3,4)*actual['D'])*feedback_cap
SOURCE = Q(response['ledger']['omitted_full_H_forcing_X_upper'])
need(SOURCE == Q(52162547172519,2), 'fixed affine response units, not a raw source derivative cap')
KF = SOURCE+Gcap
Rmajorant = Q(3,2)*(Rcoef+Rfeedback)+KF
same = {'frozen_coupled_X_operator_cap':q0,'frozen_first_feedback_whole_cone_X_cap':B1_cap,
        'W_minus_B0_whole_cone_X_cap':feedback_cap,'lower_order_coefficient_residual_H_cap':Rcoef,
        'lower_order_feedback_residual_H_cap':Rfeedback,'affine_bulk_response_KF_X_cap':KF,
        'bulk_residual_response_X_cap':Rmajorant}
need(all(v == Q(prior[name]) for name,v in same.items()) and
     Rmajorant == Q(old['bulk_residual_response_X_cap']), 'every bulk-residual majorant input preserved under direct reproof')
qret = Q(old['retained_coupled_X_operator_cap'])
beta_ret = (qret-b)/(1-b)
direct_E = (einitial+beta_ret*Rmajorant)/(1-qret)
need(direct_E == oldE and direct_E < old_readable and
     beta_ret == Q(old['retained_actual_bulk_update_contraction_cap']), 'same numerical one-update C1 upper certificate rederived for the new candidate')
direct_correction_cap = Rmajorant/(1-b)
need(direct_correction_cap == Q(old['first_actual_bulk_correction_X_cap']), 'same correction majorant re-proved from new-candidate residual bound')
delta_Bs = free_factor*kappaH
delta_W0, delta_W1 = W0factor*kappaH, W1factor*kappaH
newE = oldE+delta_W1
readable_newE = Q('3.588e25')
need(newE < readable_newE, 'strict readable new-candidate C1 target error')

# Charge the changed reduced forcing defect instead of silently reusing its
# former exact cap. Fixed F cancels, and r1=La(W1-W0) for either construction.
delta_correction = Pcap*(k*delta_W0+k0*delta_Bs)
La_first_order = k/Q(3, 2)
old_r1 = La_first_order*Q(old['first_actual_bulk_correction_X_cap'])
delta_r1 = La_first_order*delta_correction
new_r1 = old_r1+delta_r1
need(new_r1 < Q('3.574e22'), 'same readable reduced forcing-defect cap after explicit charge')

# The perturbation cannot reach any admitted boundary event. Keeping Pcap
# above is conservative; this support-refined comparison is recorded too.
threshold = R-delta
need(Q(9, 8) < Q(5, 4) < threshold < R and threshold > 1, 'selected event and complete admitted shell precede changed-data arrival')
W0_no_boundary = (1+k0)*free_factor*kappaH
W1_no_boundary = (1+k*(1+k0))*free_factor*kappaH
need(0 < W1_no_boundary < delta_W1, 'causal support sharpening of the general inverse bound')
print(json.dumps(encode({'status': 'PASS',
    'new_candidate': 'W1_dagger, a new finite update; old W1 and its archived results are preserved',
    'change': 'replace only u_nom outer affine extension by its curvature-matched quadratic extension before the same compact even mollifier; interior nominal data and P1 velocity unchanged',
    'fixed_functions': 'same frozen L0, actual La and constant M, same actual first field and affine sources F,G; no source midpoint or zero substitution',
    'input_sha256': {'initial_second_report': digest, 'old_bulk_bound': BULK_PIN, 'profiles': PROFILE_PIN, 'response_report': RESPONSE_PIN,
        'old_initial_report':OLD_INITIAL_PIN,'refined_response':REFINED_RESPONSE_PIN,
        'Robin_transfer':ROBIN_TRANSFER_PIN,'affine_shell_source':SHELL_SOURCE_PIN},
    'labels': labels, 'Gram': gram, 'R': R, 'delta_u': delta,
    'nominal_left_curvature_intervals': kappa, 'kappa_H_squared_upper': square, 'kappa_H_upper': kappaH,
    'half_kernel_moments': {str(j): half_moment(j) for j in range(3)},
    'delta_u_s_H_cap': displacement_factor*kappaH,
    'delta_u_s_radial_derivative_H_cap': radial_derivative_factor*kappaH,
    'delta_Bs_X_cap': delta_Bs,
    'frozen_H_operator_caps': frozen_caps, 'frozen_block_cap_checks': block_checks,
    'k0_frozen_bulk_X_operator_cap': k0, 'k_actual_bulk_X_operator_cap': k,
    'b_actual_boundary_X_operator_cap': b, 'P_boundary_inverse_X_cap': Pcap,
    'delta_W0_X_cap': delta_W0, 'delta_W1_X_cap': delta_W1,
    'old_target_W1_error_exact_upper': oldE,
    'new_target_W1_dagger_error_X_upper': newE,
    'new_target_W1_dagger_error_readable_strict_upper': readable_newE,
    'best_reproved_global_C1_error_X_cap':direct_E,
    'best_reproved_global_C1_readable_strict_upper':old_readable,
    'direct_reproof': {'extension_component_curvature_caps_unchanged':True,
        'M2_H':M2H,'displacement_smoothing_value_H_cap':au,
        'displacement_smoothing_radial_derivative_H_cap':bu,'velocity_smoothing_H_cap':av,
        'new_initial_free_wave_X_error_cap':einitial,'new_smooth_Bs_X_cap':Bs_cap,
        'boundary_correction_X_cap':C_boundary,'residual_majorants':same,
        'new_actual_bulk_correction_X_cap':direct_correction_cap,
        'new_reduced_forcing_defect_H_cap':(k/Q(3,2))*direct_correction_cap,
        'new_target_W1_dagger_error_X_upper':direct_E,
        'new_target_W1_dagger_error_readable_strict_upper':old_readable,
        'interpretation':'same numerical upper certificate re-proved for a different auxiliary candidate; not equality of actual errors'},
    'forcing_defect_change': {'definition': 'r1_dagger-r1=La[(W1_dagger-W0_dagger)-(W1-W0)]',
        'delta_W1_minus_W0_X_cap': delta_correction,
        'delta_reduced_forcing_defect_H_cap': delta_r1,
        'old_reduced_forcing_defect_H_cap': old_r1,
        'new_reduced_forcing_defect_H_cap': new_r1,
        'new_reduced_forcing_defect_readable_strict_upper': Q('3.574e22')},
    'per_kappa_H_factors': {'free_wave': free_factor, 'W0': W0factor, 'W1': W1factor},
    'causal_support': {'within_retained_footprint_initial_change': [threshold, R],
        'all_finite_update_changes_vanish_when': 'x+r<R-delta_u',
        'selected_point': ['9/8', '0'], 'selected_point_and_neighborhood_unchanged': True,
        'all_admitted_shell_events_unchanged': True,
        'PR98_selected_candidate_first_jets_unchanged': True,
        'PR99_selected_direct_constraint_values_unchanged': True},
    'support_refined_bounds': {'boundary_inverse_is_identity_on_differences': True,
        'delta_W0_X_cap': W0_no_boundary, 'delta_W1_X_cap': W1_no_boundary,
        'new_target_W1_dagger_error_X_upper': oldE+W1_no_boundary},
    'scope': 'conditional C1 comparison for a newly defined auxiliary finite candidate and causal inheritance of earlier selected-point values; does not repair every higher-initial-jet error, supply complete F_x, certify C2 accuracy or full Einstein constraints, solve the self-consistent second-order problem, construct a nonlinear family, or establish a Big Bang origin'},), sort_keys=True))
