#!/usr/bin/env python3
"""Independent C1 transfer to an exact actual-Robin reconstruction.

Uses authenticated JSON only. It imports no scientific project code and no new
point evaluator. The inherited analytic PDE and zero-collar premises remain
conditional; hashes authenticate the recorded inputs, not those proofs.
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


def rootup(q):
    need(q >= 0, 'nonnegative squared norm')
    n = math.isqrt(q.numerator // q.denominator)
    if n*n < q:
        n += 1
    need(n*n >= q, 'outward square root')
    return Q(n)


def norm2(row, gram):
    return sum((a*a*w for a, w in zip(row, gram)), Q(0))


parser = argparse.ArgumentParser()
parser.add_argument('--repo-root', required=True, type=Path)
parser.add_argument('--source-report', type=Path)
args = parser.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
profile = read(base/'compute/response-normal/RESPONSE_PROFILES.json',
               '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
response = read(base/'compute/response-normal/RESPONSE_REPORT.json',
                'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
initial = read(base/'initial/results/INITIAL_REPORT.json',
               'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
source_path = args.source_report or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030/reference/AFFINE_SOURCE.json'
source = read(source_path,
              'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e')
# Authenticate the two reports underlying the inherited uniform source bound.
source_input_pins = {
    'COMMON_HARMONIC_SHELL_NORMAL.json': '559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a',
    'FULL23_SHELL_NORMAL.json': '0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe'}
shell = read(root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json',
             source_input_pins['COMMON_HARMONIC_SHELL_NORMAL.json'])
packing = read(base/'shell/verified/FULL23_SHELL_NORMAL.json',
               source_input_pins['FULL23_SHELL_NORMAL.json'])
need(source['input_sha256'] == source_input_pins and source['status'] == 'PASS',
     'source provenance and status')
need(shell['status'] == 'PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING'
     and packing['status'] == 'PASS', 'inherited source scope')
need(source['uniform_domain'] == '0<=x<=3,r=0', 'uniform source coverage')
need(response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR' and
     initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER',
     'inherited response/initial scope')
need(response['domain'] == {'r_min': 'max(0,x-5/4)', 'x': ['0', '3/2'],
                           'x_plus_r_max': '7/4'}, 'registered retained cone')
need('zero-wedge factors3/4 bulk and1/4 Robin' in
     response['full_H_source_and_contraction_scope'], 'inherited zero-wedge shell history')
labels = profile['labels']
gram = list(map(Q, profile['gram']))
need(len(labels) == len(gram) == 23 and min(gram) > 0 and
     labels == initial['labels'] == response['labels'] == packing['labels'] and
     gram == list(map(Q, initial['gram'])) == list(map(Q, response['gram'])),
     'common positive modal space')
parity = [-1 if label.startswith(('L:', 'v:')) else 1 for label in labels]
need(parity.count(-1) == 5 and profile['parity']['u'] == profile['parity']['v'] == parity
     and profile['parity']['du'] == [-z for z in parity], 'mixed boundary split')
blocks = [[i for i, label in enumerate(labels) if label.endswith(':' + suffix)]
          for suffix in ('0', '2', '4', '+', '-')]
need(sum(map(len, blocks)) == 23, 'complete orthogonal block partition')
sameblock = lambda i, j: any(i in b and j in b for b in blocks)
ops = {key: Q(response['constant_matrix_H_operator_upper'][key])
       for key in ('A', 'B', 'D', 'M')}
matrices = {key: [list(map(Q, row)) for row in profile['matrices'][key]]
            for key in ops}
for key, matrix in matrices.items():
    need(len(matrix) == 23 and all(len(row) == 23 for row in matrix), 'matrix shape')
    for i in range(23):
        for j in range(23):
            need(sameblock(i, j) or matrix[i][j] == 0, 'off-block frozen entry')
    for block in blocks:
        need(sum(gram[i]*matrix[i][j]**2/gram[j] for i in block for j in block)
             <= ops[key]**2, 'frozen block operator cap')
M0 = matrices['M']
entry_error = [list(map(Q, row)) for row in response['coefficient_entry_error_caps']['M']]
actualM = [[tuple(map(Q, box)) for box in row]
           for row in response['actual_coefficient_boxes']['M']]
deltaM = Q(response['coefficient_error_H_operator_upper']['M'])
for i in range(23):
    for j in range(23):
        lo, hi = actualM[i][j]
        need(lo <= hi and max(abs(lo-M0[i][j]), abs(hi-M0[i][j])) <= entry_error[i][j],
             'actual-vs-frozen Robin entry enclosure')
        need(sameblock(i, j) or entry_error[i][j] == 0, 'off-block actual error')
        if parity[i] < 0 or parity[j] < 0:
            need(M0[i][j] == 0 and (lo, hi) == (0, 0), 'Robin only rows/columns')
for block in blocks:
    need(sum(gram[i]*entry_error[i][j]**2/gram[j] for i in block for j in block)
         <= deltaM**2, 'actual-vs-frozen Robin H norm cap')
need(deltaM > 0, 'nonzero actual coefficient uncertainty retained')
m = ops['M'] + deltaM
actual_bulk_ops = {}
for key in ('A', 'B', 'D'):
    errs = [list(map(Q, row)) for row in response['coefficient_entry_error_caps'][key]]
    delta = Q(response['coefficient_error_H_operator_upper'][key])
    for i in range(23):
        for j in range(23):
            lo, hi = map(Q, response['actual_coefficient_boxes'][key][i][j])
            frozen = matrices[key][i][j]
            need(lo <= hi and max(abs(lo-frozen), abs(hi-frozen)) <= errs[i][j],
                 'actual bulk coefficient entry error')
            need(sameblock(i, j) or errs[i][j] == 0, 'actual bulk off-block error')
    for block in blocks:
        need(sum(gram[i]*errs[i][j]**2/gram[j] for i in block for j in block)
             <= delta**2, 'actual bulk coefficient error H cap')
    actual_bulk_ops[key] = ops[key] + delta
Gcap = Q(source['source_H_norm_strict_integer_upper'])
need(Q(source['source_H_norm_square_upper']) < Gcap**2 and Gcap <= 487122,
     'strict inherited full source norm cap')
need([row['label'] for row in source['rows']] == labels and
     all(Q(source['rows'][i]['abs_upper']) == 0 for i in range(23) if parity[i] < 0),
     'source modal space and exact Dirichlet rows')

# Independently recompute the tighter nominal-convolution transfer ledger.
nodes = list(map(Q, profile['nodes']))
tab = {key: [list(map(Q, row)) for row in profile['profiles'][key]]
       for key in ('u', 'du', 'v')}
delta_u, delta_v = Q(1, 2**30), Q(1, 2048)
need(nodes[0] == 0 and nodes[-1] == Q(7, 4) and
     all(b > a for a, b in zip(nodes, nodes[1:])), 'initial mesh')
need(nodes[1] - max(delta_u, delta_v) > 1 and
     all(z == 0 for table in tab.values() for row in table[:2] for z in row),
     'smoothed nominal flat shell collar')
need(min(b-a for a, b in zip(nodes, nodes[1:])) >= 2*delta_v,
     'nonoverlapping velocity smoothing bands')
slopes = [[(tab['v'][k+1][i]-tab['v'][k][i])/(nodes[k+1]-nodes[k])
           for i in range(23)] for k in range(len(nodes)-1)]
jumps = [[b-a for a, b in zip(left, right)] for left, right in zip(slopes, slopes[1:])]
av = Q(35, 256)*delta_v*rootup(max(norm2(row, gram) for row in jumps))
M2H = Q(initial['nominal_second_derivative_H_caps'][0])
need(initial['derivative_convention']['second_derivatives'][0] == 'partial_r^2 U'
     and 'U and V are C1' in initial['regularity'], 'nominal C1,1 premise')
au, bu = M2H*delta_u**2/18, 35*M2H*delta_u/128
eU, eUr, eV = map(Q, initial['P1_H_error_caps'])
need([eU, eUr, eV] == [Q(response['ledger'][key]) for key in
                     ('epsilon_u_H', 'epsilon_du_H', 'epsilon_v_H')], 'old P1 ledger')
Eold = Q(response['ledger']['total_X_error_upper'])
free_transfer = max((au+Q(3, 2)*av)/3, bu+av)
feedback_transfer = Q(3, 2)*((ops['A']+ops['B'])*(eUr+bu+av) +
                              ops['D']*(eU+au+Q(3, 2)*av)) + ops['M']*(eU+au+av/4)
Eprior = Eold + free_transfer + feedback_transfer
need(Eprior < Q(2184*10**25), 'inherited tighter candidate transfer')


def bound_history(end):
    """Bound V, its actual-law residual and the exact boundary correction.

    Only initial data through end+max(delta) enter the shell history up to end.
    Positive convolution transfers cellwise continuous caps to smooth profiles.
    """
    end = Q(end)
    need(1 < end <= Q(5, 4), 'retained active shell history')
    uc, pc, vc = ([Q(0)]*23 for _ in range(3))
    for k, cell in enumerate(initial['interpolation_cells']):
        need(list(map(Q, cell['r'])) == nodes[k:k+2], 'cell geometry')
        if nodes[k] > end + max(delta_u, delta_v):
            break
        for i in range(23):
            uc[i] = max(uc[i], max(abs(tab['u'][k][i]), abs(tab['u'][k+1][i])) +
                        Q(cell['component_error_caps'][0][i]))
            pc[i] = max(pc[i], max(abs(tab['du'][k][i]), abs(tab['du'][k+1][i])) +
                        Q(cell['component_error_caps'][1][i]))
            vc[i] = max(vc[i], abs(tab['v'][k][i]), abs(tab['v'][k+1][i]))
    U0, U1, V = (rootup(norm2(row, gram)) for row in (uc, pc, vc))
    length = end - 1
    fcap = (ops['A']+ops['B'])*(U1+V) + ops['D']*(U0+end*V)
    free_trace = U0 + length*V
    # At the shell, the active PHYSICAL triangle has area (end^2-1)/2.
    # Its reflected complete output is an H-isometric parity extension.
    bulk_trace = (end**2-1)*fcap/2
    lift_trace = length*ops['M']*free_trace
    B1_trace = bulk_trace + lift_trace
    V_trace = free_trace + B1_trace
    Rcap = ops['M']*B1_trace + deltaM*V_trace + Gcap
    contraction = m*length
    need(0 <= contraction < 1, 'boundary correction contraction')
    qcap = Rcap/(1-contraction)
    # q+M Jq=-R; C=Lq, C_x=-q, C_r=q and |C|<=length*qcap.
    C_X = qcap
    C_value = length*qcap
    total = Eprior+C_X
    # The wave-principal part of C=Lq vanishes. For the reduced coupled
    # residual principal-minus-lower-minus-F, only its CHANGE is enclosed:
    # R[W]-R[V]=(A_actual-B_actual)q+D_actual Jq.
    residual_change = (actual_bulk_ops['A']+actual_bulk_ops['B'] +
                       length*actual_bulk_ops['D'])*qcap
    need(qcap >= Rcap and qcap == Rcap + contraction*qcap,
         'resolvent bound solves majorant equation')
    return {key: str(value) for key, value in {
        'shell_history_end': end, 'active_history_length': length,
        'nominal_u_H_cap': U0, 'nominal_du_H_cap': U1, 'nominal_v_H_cap': V,
        'frozen_physical_bulk_source_H_cap': fcap,
        'free_shell_trace_H_cap': free_trace,
        'reflected_bulk_shell_trace_H_cap': bulk_trace,
        'frozen_Robin_lift_shell_trace_H_cap': lift_trace,
        'B1_shell_trace_H_cap': B1_trace,
        'V_shell_trace_H_cap': V_trace,
        'actual_law_candidate_boundary_residual_H_cap': Rcap,
        'boundary_correction_contraction': contraction,
        'correction_q_H_cap': qcap,
        'correction_C1_X_cap': C_X,
        'correction_value_H_cap': C_value,
        'reduced_coupled_bulk_residual_change_H_cap': residual_change,
        'new_conditional_C1_X_error_cap': total}.items()}


global_history = bound_history(Q(5, 4))
sample_history = bound_history(Q(9, 8))
need(Q(global_history['correction_C1_X_cap']) < Q(6344*10**21),
     'readable global correction cap')
need(Q(global_history['new_conditional_C1_X_error_cap']) < Q(2185*10**25) < 10**29,
     'new global target retained')
need(Q(sample_history['correction_C1_X_cap']) <= Q(global_history['correction_C1_X_cap']),
     'shorter shell history does not enlarge correction cap')
print(json.dumps({
    'status': 'PASS',
    'method': 'independent JSON-only shell residual and Volterra contraction majorant',
    'labels': labels,
    'inherited_prior_exact_C1_X_error': str(Eprior),
    'frozen_Robin_M0_H_cap': str(ops['M']),
    'actual_Robin_M_minus_M0_H_cap': str(deltaM),
    'actual_Robin_M_H_cap': str(m),
    'actual_bulk_operator_H_caps': {key: str(value) for key, value in actual_bulk_ops.items()},
    'affine_shell_G_H_cap': str(Gcap),
    'global_retained_history': global_history,
    'sample_history': sample_history,
    'global_readable_correction_X_cap': '6344*10^21',
    'global_readable_new_C1_X_error_cap': '2185*10^25',
    'new_C1_target_1e29_met': True,
    'correction_identity': 'C=Lq; q+M Jq=-R_V; R_V=-M0 PN B1-(M-M0)PN V-G',
    'scope': 'conditional C1 transfer to integral-defined exact actual-Robin subproblem reconstruction; frozen bulk forcing is not the self-consistent coupled solution, no numerical C2 error or physical observable'
}, sort_keys=True))
