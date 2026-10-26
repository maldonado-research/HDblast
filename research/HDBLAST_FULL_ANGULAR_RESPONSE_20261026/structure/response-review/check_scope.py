#!/usr/bin/env python3
"""Independent bounded arithmetic/support checks; no project imports."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
def need(ok, message):
    if not ok:
        raise ValueError(message)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--capsule-root', type=Path, default=BASE)
parser.add_argument('--normal-results', type=Path)
parser.add_argument('--optimized-results', type=Path)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
ROOT = args.capsule_root.resolve()
normal_results = (args.normal_results or ROOT/'compute/response-normal').resolve()
optimized_results = args.optimized_results.resolve() if args.optimized_results else None
target = args.output.resolve()
need(not target.exists(), 'FRESH_OUTPUT_REQUIRED')
need(target != ROOT and ROOT not in target.parents, 'OUTPUT_INSIDE_CAPSULE')
report_path = normal_results/'RESPONSE_REPORT.json'
profile_path = normal_results/'RESPONSE_PROFILES.json'
report = json.loads(report_path.read_text())
profile = json.loads(profile_path.read_text())
ledger = {k:F(v) for k,v in report['ledger'].items()}
q = F(3,2)*(2*F(1,10000)+F(3,4)*F(1,1000000))+F(1,4)*F(1,2000)
need(q == F(3409,8000000) == ledger['retained_full_H_contraction'], 'CONTRACTION_ARITHMETIC')
source = F(3,2)*17387513156573+3851400
need(source == F(52162547172519,2) == ledger['omitted_full_H_forcing_X_upper'], 'SOURCE_ARITHMETIC')
need(ledger['e0'] == ledger['free_input_error_X_upper']+ledger['epsilon_v_H']+source, 'E0')
need(ledger['e1'] == F(3,2)*(ledger['delta_f_coefficient']+ledger['delta_f_profile'])+ledger['delta_g_coefficient']+ledger['delta_g_profile'], 'E1')
need(ledger['actual_response_X_upper'] == (ledger['B0_X_upper']+ledger['e0'])/(1-q), 'UCAP')
error = q*q*ledger['actual_response_X_upper']+(1+q)*ledger['e0']+ledger['e1']
need(error == ledger['total_X_error_upper'] and error < F(10**29), 'TOTAL_ERROR')
residual = ledger['e0']+(1+q)*ledger['e1']+q*q*ledger['B0_X_upper']
resolvent_error = residual/(1-q)
residual_minus_two_term = 2*q*ledger['e1']/(1-q)
need(resolvent_error-error == residual_minus_two_term, 'RESOLVENT_TWO_TERM_IDENTITY')
need(0 <= residual_minus_two_term and resolvent_error < F(10**29), 'RESOLVENT_TARGET')
angular_factor = F(255,4)
need(5+8*9+8*25+5+5 == 287 and angular_factor**2 > F(9,16)*287*24,
     'INTRINSIC_ANGULAR_FACTOR')
intrinsic_angular_error = angular_factor*error
cstar = ledger['selected_evaluation_factor_upper']
need(cstar > 0 and cstar*cstar >= F(321,19), 'CSTAR_UPPER')
need(ledger['zero_wedge_pointwise_value_error_upper'] == F(3,4)*cstar*error, 'POINTWISE_VALUE')
need(ledger['pointwise_scaled_derivative_error_upper'] == cstar*error, 'POINTWISE_DERIVATIVE')

# Cutoff polynomial is C4 when joined to zero/one; sufficient for the
# registered displacement, which uses at most three packet derivatives.
cutoff = [F(0)]*5+[F(126),F(-420),F(540),F(-315),F(70)]
for order in range(5):
    at0 = cutoff[0]
    at1 = sum(cutoff)
    need(at0 == 0 and at1 == (1 if order == 0 else 0), 'CUTOFF_ENDPOINT_JET')
    cutoff = [F(i)*cutoff[i] for i in range(1,len(cutoff))]
nodes = list(map(F,profile['nodes']))
need(nodes[0] == 0 and nodes[1] == F(65,64) and nodes[-1] == F(7,4), 'INITIAL_GRID_DOMAIN')
need(all(a < b for a,b in zip(nodes,nodes[1:])), 'GRID_MONOTONE')
for key in ('u','du','v'):
    rows = profile['profiles'][key]
    need(len(rows) == len(nodes) and all(len(row)==23 for row in rows), 'PROFILE_SHAPES')
    need(all(F(v)==0 for row in rows[:2] for v in row), 'ZERO_INITIAL_COLLAR')
    need(all(F(row[i])==0 for row in rows for i in (21,22)), 'ZERO_TT_INITIAL_TABLES')
if optimized_results is not None:
    need(report_path.read_bytes() == (optimized_results/'RESPONSE_REPORT.json').read_bytes(), 'NORMAL_OPTIMIZED_REPORT')

out = {
    'status':'PASS_RESPONSE_SCOPE_ARITHMETIC_AND_SUPPORT',
    'producer_sha256':sha(ROOT/'compute/full23_response.py'),
    'kernel_sha256':sha(ROOT/'compute/angular_ridge.py'),
    'report_sha256':sha(report_path),
    'profiles_sha256':sha(profile_path),
    'full_H_contraction':str(q),
    'source_omission_X_upper':str(source),
    'candidate_zero_on':'x+r<65/64',
    'actual_and_error_zero_on':'x+r<1, using retained exact-source premise',
    'cutoff_endpoint_derivatives_through_order4':True,
    'total_X_error_upper':str(error),
    'complete_integral_residual_X_upper':str(residual),
    'resolvent_error_X_upper':str(resolvent_error),
    'resolvent_minus_two_term_error':str(residual_minus_two_term),
    'identity_2q_e1_over_1_minus_q':True,
    'passes_two_term_and_resolvent_1e29':True,
    'angular_equivariant_section_dimension':287,
    'derivative_image_dimension_upper':287,
    'angular_rough_spectral_cap':24,
    'intrinsic_angular_derivative_error_factor':str(angular_factor),
    'intrinsic_angular_derivative_error_upper':str(intrinsic_angular_error),
    'intrinsic_angular_scope':'Intrinsic D_kappa on the unit S3, with the positive Frobenius norm summed over angular derivative directions and tensor slots; not the physical spacetime covariant derivative.',
    'pointwise_value_error_upper':str(F(3,4)*cstar*error),
    'pointwise_scaled_derivative_error_upper':str(cstar*error),
    'Robin_residual_H_upper_from_exact_solution_error':str((1+F(1,8000))*error),
    'interpretation':'The Robin residual bound is a conditional consequence of the exact response boundary law, full X error, zero-wedge trace bound, and retained ||M||<=1/2000. It is not a new direct boundary residual evaluation.',
    'limits':['No frozen producer/import executed','No new coefficient or initial interval verification','No full response re-evaluation','No C2 Einstein residual']}
with target.open('x') as handle:
    handle.write(json.dumps(out,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':out['status'],'error_display':float(error),'Robin_residual_H_upper_display':float((1+F(1,8000))*error)}))
