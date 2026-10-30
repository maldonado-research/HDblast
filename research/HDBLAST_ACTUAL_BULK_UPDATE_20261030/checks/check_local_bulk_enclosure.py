#!/usr/bin/env python3
"""Local causal majorant for the first actual-bulk correction at the shell.

This uses genuine uniform pointwise affine-source caps and local initial
profile caps. It does not quadrature-evaluate the actual bulk correction.
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


def iv(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(Q, row)
    need(lo <= hi, 'interval order')
    return lo, hi


def sqrt_integer(value):
    need(value >= 0, 'nonnegative squared norm')
    n = math.isqrt(value.numerator // value.denominator)
    if n*n < value:
        n += 1
    need(n*n >= value, 'outward integer square root')
    return Q(n)


def sqrt_grid(value):
    need(value >= 0, 'nonnegative squared row norm')
    grid = 1 << 112
    n = math.isqrt(value.numerator*grid*grid // value.denominator)
    if n*n*value.denominator < value.numerator*grid*grid:
        n += 1
    q = Q(n, grid)
    need(q*q >= value, 'outward rational square root')
    return q


def wide(row, radius):
    need(radius >= 0, 'nonnegative radius')
    return row[0]-radius, row[1]+radius


def enc(row):
    return list(map(str, row))


def inside(inner, outer):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--initial-capsule', type=Path)
ap.add_argument('--robin-capsule', type=Path)
ap.add_argument('--contact-capsule', type=Path)
ap.add_argument('--bound-report', type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
initial_capsule = args.initial_capsule or root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030'
robin_capsule = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
contact_capsule = args.contact_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
pins = {
    'prior_response': 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
    'prior_point': 'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69',
    'initial_report': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'profiles': '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc',
    'response': 'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd',
    'affine_shell': 'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e',
    'actual_affine_bulk': '45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f',
    'complete_geometric_source_map': '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'}
pins['bulk_update_majorant'] = '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f'
prior = read(initial_capsule/'reference/RESPONSE.json', pins['prior_response'])
point = read(robin_capsule/'reference/RESUMMED.json', pins['prior_point'])
initial = read(base/'initial/results/INITIAL_REPORT.json', pins['initial_report'])
profile = read(base/'compute/response-normal/RESPONSE_PROFILES.json', pins['profiles'])
response = read(base/'compute/response-normal/RESPONSE_REPORT.json', pins['response'])
bound_path = args.bound_report or Path(__file__).resolve().parents[1]/'reference/BULK_BOUND.json'
bound = read(bound_path, pins['bulk_update_majorant'])
shell = read(contact_capsule/'reference/AFFINE_SOURCE.json', pins['affine_shell'])
source = read(root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json',
              pins['actual_affine_bulk'])
map_path = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
map_raw = map_path.read_bytes()
need(hashlib.sha256(map_raw).hexdigest() == pins['complete_geometric_source_map'],
     'complete source map hash')
need(all(d['status'] == 'PASS' for d in (prior, point, shell, bound)), 'successful inherited certificates')
need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER'
     and response['status'] == 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'formal response scope')
need(source['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES'
     and source['scope']['actual_evolved_common_gauge_bulk_sources_evaluated'] is True,
     'actual common-gauge affine source scope')
need(source['source_map'] == {'bytes': len(map_raw),
     'sha256': pins['complete_geometric_source_map'],
     'status': 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'},
     'complete geometric source premise')
need(source['domain'] == {'ell': '1/1000000',
     'causal_triangle': '0<=x<=3,r>=0,x+r<=8',
     'angular': 'entire unit S3; exact unnormalized registered scalar/vector/tensor coefficient basis'},
     'uniform actual bulk source domain')
need(source['normalization'] ==
     'Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.',
     'actual affine bulk source normalization')
labels = profile['labels']
weights = list(map(Q, profile['gram']))
need(len(labels) == len(weights) == 23 and min(weights) > 0 and
     labels == point['labels'] == initial['labels'] == response['labels'] and
     weights == list(map(Q, initial['gram'])) == list(map(Q, response['gram'])),
     'registered positive modal norm')
need([r['label'] for r in prior['full23_point_table']] == labels
     and [r['label'] for r in shell['rows']] == labels, 'prior row order')
need(response['domain'] == {'r_min': 'max(0,x-5/4)', 'x':['0','3/2'],
                           'x_plus_r_max':'7/4'}, 'retained cone')
need(initial['derivative_convention']['du'] == 'partial_r U' and
     initial['derivative_convention']['v'] == 'partial_x U; x=t/ell',
     'initial derivative conventions')
parity = [-1 if z.startswith(('L:', 'v:')) else 1 for z in labels]
need(profile['parity']['u'] == profile['parity']['v'] == parity and
     profile['parity']['du'] == [-s for s in parity], 'compatible initial reflection')
blocks = [[i for i,z in enumerate(labels) if z.endswith(':'+suffix)]
          for suffix in ('0','2','4','+','-')]
need(sum(map(len, blocks)) == 23, 'complete modal blocks')
sameblock = lambda i,j: any(i in block and j in block for block in blocks)
op = {k: Q(response['constant_matrix_H_operator_upper'][k]) for k in ('A','B','D','M')}
error = {k: Q(response['coefficient_error_H_operator_upper'][k]) for k in op}
actual = {k: op[k]+error[k] for k in op}
M = None
for key in op:
    frozen = [list(map(Q,row)) for row in profile['matrices'][key]]
    errs = [list(map(Q,row)) for row in response['coefficient_entry_error_caps'][key]]
    boxes = [[iv(z) for z in row] for row in response['actual_coefficient_boxes'][key]]
    need(profile['matrices'][key] == response['constant_matrices'][key], 'frozen matrix identity')
    need(len(frozen) == len(errs) == len(boxes) == 23 and
         all(len(row) == 23 for matrix in (frozen,errs,boxes) for row in matrix), 'matrix shapes')
    for i in range(23):
        for j in range(23):
            lo,hi = boxes[i][j]
            need(max(abs(lo-frozen[i][j]),abs(hi-frozen[i][j])) <= errs[i][j],
                 'actual coefficient entry uncertainty')
            need(sameblock(i,j) or frozen[i][j] == errs[i][j] == 0, 'off-block closure')
            if key == 'M' and (parity[i] < 0 or parity[j] < 0):
                need((lo,hi) == (0,0) and frozen[i][j] == 0, 'mixed Robin closure')
    for block in blocks:
        need(sum(weights[i]*frozen[i][j]**2/weights[j] for i in block for j in block)
             <= op[key]**2, 'frozen H operator norm')
        need(sum(weights[i]*errs[i][j]**2/weights[j] for i in block for j in block)
             <= error[key]**2, 'coefficient-error H operator norm')
    need(actual[key] == Q(prior['actual_coefficient_H_norm_caps'][key]), 'same actual coefficient cap')
    if key == 'M':
        M = boxes

# Recover the genuine pointwise source norm in the current H normalization.
summary = source['global_cover']['summary']
component = summary['absolute_coefficient_caps']
bundle = summary['positive_bundle_bulk_forcing_H_bound']
mapping = [
    ('yy','metric_forcing/yy',('K:0','K:2','K:4'),1),
    ('y0','metric_forcing/y0',('L:0','L:2','L:4'),2),
    ('00','metric_forcing/00',('N:0','N:2','N:4'),1),
    ('ya','metric_forcing/ya',('v:2','v:4'),2),
    ('0a','metric_forcing/0a',('u:2','u:4'),2),
    ('ab','metric_forcing/ab',('J:0','J:2','J:4','H:2','H:4','T:+','T:-'),1),
    ('scalar','scalar_forcing',('w:0','w:2','w:4'),1)]
Fcomponents = [Q(0)]*23
source_square = Q(0)
gram_entries = 0
for field,prefix,fieldlabels,multiplicity in mapping:
    entry = bundle['block_Gram_certificates'][field]
    gram = [list(map(Q,row)) for row in entry['gram']]
    n = len(fieldlabels)
    need(entry['multiplicity'] == multiplicity and len(gram) == n and
         all(len(row) == n for row in gram), 'source Gram shape')
    cc = [Q(component[prefix+'/'+str(i)]) for i in range(n)]
    need(min(cc) >= 0, 'positive affine source caps')
    for i,name in enumerate(fieldlabels):
        Fcomponents[labels.index(name)] = cc[i]
        for j in range(n):
            need(multiplicity*gram[i][j] == (weights[labels.index(name)] if i == j else 0),
                 'source/current exact Gram normalization')
            gram_entries += 1
    square = multiplicity*sum(abs(gram[i][j])*cc[i]*cc[j] for i in range(n) for j in range(n))
    need(square == Q(entry['squared_norm_upper']), 'source block squared norm')
    source_square += square
need(source_square == Q(bundle['squared_upper']) ==
     sum(w*c*c for w,c in zip(weights,Fcomponents)), 'whole source H norm')
FH = Q(bundle['F_H_upper'])
need(FH == 17387513156573 and (FH-1)**2 <= source_square < FH**2,
     'genuine uniform affine F cap')
Gcap = Q(shell['source_H_norm_strict_integer_upper'])
need(shell['uniform_domain'] == '0<=x<=3,r=0' and Gcap == 487122 and
     Q(shell['source_H_norm_square_upper']) < Gcap**2, 'uniform affine shell cap')

# Positive convolution transfers local nominal/P1 caps to the smoothed data.
a,d = Q(9,8),Q(1,8)
udelta,vdelta = Q(1,2**30),Q(1,2048)
nodes = list(map(Q,profile['nodes']))
tab = {name: [list(map(Q,row)) for row in profile['profiles'][name]]
       for name in ('u','du','v')}
need(nodes[0] == 0 and nodes[-1] == Q(7,4) and
     nodes[1]-max(udelta,vdelta) > 1 and
     all(z == 0 for table in tab.values() for row in table[:2] for z in row),
     'smoothed zero collar and local footprint')
need(len(initial['interpolation_cells'])+1 == len(nodes), 'initial interpolation partition')
uc,pc,vc = ([Q(0)]*23 for _ in range(3))
used_cells = 0
for n,cell in enumerate(initial['interpolation_cells']):
    need(list(map(Q,cell['r'])) == nodes[n:n+2], 'initial cell endpoints')
    if nodes[n] > a+max(udelta,vdelta):
        break
    used_cells += 1
    for i in range(23):
        uc[i] = max(uc[i],max(abs(tab['u'][n][i]),abs(tab['u'][n+1][i]))+
                    Q(cell['component_error_caps'][0][i]))
        pc[i] = max(pc[i],max(abs(tab['du'][n][i]),abs(tab['du'][n+1][i]))+
                    Q(cell['component_error_caps'][1][i]))
        vc[i] = max(vc[i],abs(tab['v'][n][i]),abs(tab['v'][n+1][i]))
need(used_cells and nodes[used_cells] >= a+max(udelta,vdelta), 'local convolution coverage')
U0,U1,V = (sqrt_integer(sum(w*z*z for w,z in zip(weights,c))) for c in (uc,pc,vc))
Bsvalue = U0+a*V
Bsderivative = U1+V
Bstrace = U0+d*V
f0cap = (op['A']+op['B'])*Bsderivative+op['D']*Bsvalue
area = (a*a-1)/2
need(area == Q(17,128) and area > d*d/2, 'complete active physical triangle geometry')
m = actual['M']
need(0 < m*d < 1, 'local boundary trace contraction')
hcap = Bstrace+area*f0cap
Zcap = (hcap+d*Gcap)/(1-m*d)
gcap = m*Zcap+Gcap
need(gcap == (m*hcap+Gcap)/(1-m*d), 'equivalent total boundary-load majorants')
WminusBs_value = area*f0cap+d*gcap
WminusBs_derivative = a*f0cap+gcap
Rlocal = ((error['A']+error['B'])*Bsderivative+error['D']*Bsvalue+
          (actual['A']+actual['B'])*WminusBs_derivative+
          actual['D']*WminusBs_value+FH)
Dtrace = area*Rlocal/(1-m*d)
Ctime = a*Rlocal+m*Dtrace
need(Dtrace == area*Rlocal+d*m*Dtrace, 'local homogeneous Robin trace majorant')

# Independent conservative one-update global error; no new bulk code import.
k = Q(3,2)*(actual['A']+actual['B']+Q(3,4)*actual['D'])
b = m/4
q = Q(response['ledger']['retained_full_H_contraction'])
need(k+b <= q == Q(3409,8000000) < 1, 'retained global operator cap')
beta = (q-b)/(1-b)
ei = Q(prior['initial_free_wave_X_error_cap'])
Rglobal = Q(prior['bulk_residual_response_X_cap'])
Eold = (ei+Rglobal)/(1-q)
need(Eold == Q(prior['direct_C1_X_error_using_retained_contraction_cap']), 'same prior error norm')
E1 = (ei+beta*Rglobal)/(1-q)
E = Q(3541*10**22)
need(0 < E1 < E and E1 == Q(bound['one_update_retained_target_error_X_cap'])
     and E == Q(bound['one_update_readable_error_strict_upper']),
     'independently recomputed readable global one-update error')
need(Dtrace+d*E < d*Eold and Ctime+E < Eold,
     'local value/time target radii strictly improve global norm-only radii')

rows = []
checks = strict = 0
widths = {}
zero_Robin_rows = []
for i,(oldpoint,oldtarget) in enumerate(zip(point['full23_point_table'],prior['full23_point_table'])):
    need(oldpoint['label'] == oldtarget['label'] == labels[i] and
         oldpoint['parity'] == oldtarget['parity'] == parity[i], 'point row identity')
    W0 = {name: iv(oldpoint['reconstructed_boundary_candidate'][name]) for name in ('value','dx','dr')}
    factor = sqrt_grid(1/weights[i])
    rowdual = sqrt_grid(sum(max(abs(lo),abs(hi))**2/weights[j] for j,(lo,hi) in enumerate(M[i])))
    if parity[i] > 0:
        W1 = {'value':wide(W0['value'],factor*Dtrace),
              'dx':wide(W0['dx'],factor*Ctime),
              'dr':wide(W0['dr'],rowdual*Dtrace)}
        U = {'value':wide(W0['value'],factor*(Dtrace+d*E)),
             'dx':wide(W0['dx'],factor*(Ctime+E)),
             'dr':wide(W0['dr'],rowdual*(Dtrace+d*E))}
    else:
        need(W0['value'] == W0['dx'] == (0,0), 'Dirichlet old trace')
        W1 = {'value':(Q(0),Q(0)), 'dx':(Q(0),Q(0)),
              'dr':wide(W0['dr'],factor*a*Rlocal)}
        U = {'value':(Q(0),Q(0)), 'dx':(Q(0),Q(0)),
             'dr':wide(W0['dr'],factor*(a*Rlocal+E))}
    oldbest = {'value':iv(oldtarget['conditional_actual_formal_response']['value']),
               'dx':iv(oldtarget['conditional_actual_formal_response']['dx']),
               'dr':iv(oldtarget['shared_Robin_actual_radial'])}
    ratios = {}
    for name in U:
        need(inside(U[name],oldbest[name]), 'new local target box outside old best box: '+labels[i]+' '+name)
        checks += 1
        newwidth,oldwidth = U[name][1]-U[name][0],oldbest[name][1]-oldbest[name][0]
        if newwidth < oldwidth:
            strict += 1
        ratios[name] = str(oldwidth/newwidth) if newwidth else None
    if parity[i] > 0 and rowdual == 0:
        need(rowdual == 0 and W1['dr'] == U['dr'] == W0['dr'], 'zero-M tensor radial identity')
        zero_Robin_rows.append(labels[i])
    rows.append({'label':labels[i], 'parity':parity[i],
        'locally_enclosed_actual_bulk_candidate_W1':{name:enc(z) for name,z in W1.items()},
        'conditional_actual_target_U_best_point_boxes':{name:enc(z) for name,z in U.items()},
        'old_best_target_to_new_width_ratios':ratios})
    if i == 12:
        widths = ratios
need(labels[12] == 'w:0' and all(iv(rows[12]['conditional_actual_target_U_best_point_boxes'][name])[0] > 0
                              for name in ('value','dx','dr')), 'positive formal scalar mean first jet')
need(checks == 69 and strict > 0, 'complete local point containment audit')
mean = rows[12]['conditional_actual_target_U_best_point_boxes']
for key, lower, upper in [('value',Q(829036*10**25),Q(829042*10**25)),
                           ('dx',Q(1974621*10**26),Q(1974623*10**26)),
                           ('dr',Q(2154067*10**21),Q(2154080*10**21))]:
    lo,hi = iv(mean[key])
    need(lower < lo <= hi < upper, 'readable strict local scalar-mean '+key+' bounds')
print(json.dumps({'status':'PASS', 'method':'local causal-cone norm enclosure, no bulk quadrature',
    'input_sha256':pins, 'point':['9/8','0'],
    'local_domain':'s>=0,y>=0,s+y<=9/8; inherited zero wedge s+y<=1',
    'active_physical_source_area':str(area), 'history_length':str(d),
    'local_convolution_cells_used':used_cells,
    'nominal_u_local_H_cap':str(U0), 'nominal_ur_local_H_cap':str(U1),
    'P1_velocity_local_H_cap':str(V),
    'Bs_value_local_H_cap':str(Bsvalue), 'Bs_first_derivative_local_H_cap':str(Bsderivative),
    'Bs_shell_trace_local_H_cap':str(Bstrace),
    'frozen_f0_local_H_cap':str(f0cap), 'auxiliary_trace_h_local_H_cap':str(hcap),
    'W0_trace_Z_local_H_cap':str(Zcap), 'W0_boundary_load_MZ_plus_G_H_cap':str(gcap),
    'W0_minus_Bs_value_local_H_cap':str(WminusBs_value),
    'W0_minus_Bs_first_derivative_local_H_cap':str(WminusBs_derivative),
    'genuine_pointwise_affine_F_H_cap':str(FH), 'source_Gram_entries_compared':gram_entries,
    'source_H_norm_square':str(source_square),
    'actual_reduced_bulk_residual_local_H_cap':str(Rlocal),
    'bulk_update_correction_local_value_and_trace_H_cap':str(Dtrace),
    'bulk_update_correction_shell_time_derivative_H_cap':str(Ctime),
    'global_one_update_error_exact_retained_cap':str(E1),
    'global_one_update_error_readable_strict_upper':str(E),
    'prior_to_new_scalar_mean_point_width_ratios':widths,
    'updated_best_point_box_containment_checks':checks,
    'strictly_narrower_best_point_boxes':strict,
    'zero_Robin_rows_preserving_source_only_radial_box':zero_Robin_rows,
    'full23_point_table':rows,
    'support_premise':'zero wedge and source support are inherited separately from magnitude caps',
    'source_premise':'uniform actual common-gauge source projection and exact H normalization; not SOURCE=X(KF+LG) reinterpreted as a pointwise bound',
    'resolution':'new point enclosures use local analytic majorants, not evaluation of the actual bulk convolution or its sign',
    'scope':'conditional first actual-bulk integral update and registered modal target sample; no full self-consistency, quantitative C2 error, Einstein/constraint residual, nonlinear spacetime or physical observable'}, sort_keys=True))
