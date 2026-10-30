#!/usr/bin/env python3
"""Exact audit of a past-time filtered auxiliary forcing, not actual F_x.

All arithmetic is rational. No project source is imported or modified.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def encode(z):
    if isinstance(z, Q):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (list, tuple)):
        return [encode(v) for v in z]
    return z


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    pinned.append((path, digest))
    return json.loads(raw)


def value(p, x):
    return sum(a*x**i for i, a in enumerate(p))


def derivative(p):
    return [i*a for i, a in enumerate(p) if i]


def integral(p, a=Q(0), b=Q(1)):
    return sum(c*(b**(i+1)-a**(i+1))/(i+1) for i,c in enumerate(p))


def negative_control(name, predicate):
    try:
        need(predicate, name)
    except RuntimeError:
        controls.append(name)
    else:
        raise RuntimeError('wrong formula accepted:'+name)


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--higher-capsule', type=Path)
ap.add_argument('--bulk-capsule', type=Path)
ap.add_argument('--epsilon', default='1/1024')
args = ap.parse_args()
root = args.repo_root.resolve()
higher = args.higher_capsule or root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
bulk = args.bulk_capsule or root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
pinned = []
controls = []
groups = []
transfer = pin(higher/'reference/TRANSFER.json', '4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3')
local = pin(bulk/'reference/LOCAL.json', '1832fba309ceaf0a138bbddc0aaf14602f81796bc1a7aa8cd938a74470e1b899')
source = pin(root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json',
             '45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f')
map_path = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
map_digest = '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'
map_raw = map_path.read_bytes()
need(hashlib.sha256(map_raw).hexdigest() == map_digest, 'complete map hash')
pinned.append((map_path,map_digest))
need(transfer['status'] == local['status'] == 'PASS' and
     source['status'] == 'PASS_CONTINUOUS_ACTUAL_EVOLVED_GEOMETRIC_SOURCE_CELL_ENCLOSURES',
     'successful inherited actual-source and candidate reports')
need(source['source_map'] == {'bytes':len(map_raw), 'sha256':map_digest,
     'status':'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION'}, 'complete source-map premise')
need(source['scope']['actual_evolved_common_gauge_bulk_sources_evaluated'] is True and
     source['domain']['causal_triangle'] == '0<=x<=3,r>=0,x+r<=8', 'actual source domain and scope')
groups.append('authenticated actual full-source premise')

# Independently recover the complete source norm, in the same positive H.
labels = transfer['labels']
gram = list(map(Q, transfer['Gram']))
need(len(labels) == len(gram) == 23 and min(gram) > 0, 'complete positive full23 Gram')
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
caps = [Q(0)]*23
seen = []
norm_square = Q(0)
gram_entries = 0
for field,prefix,names,multiplicity in mapping:
    block = bundle['block_Gram_certificates'][field]
    bg = [list(map(Q,row)) for row in block['gram']]
    cc = [Q(component[prefix+'/'+str(i)]) for i in range(len(names))]
    need(block['multiplicity'] == multiplicity and len(bg) == len(names) and
         all(len(row) == len(names) for row in bg) and min(cc) >= 0, 'source block shape')
    for i,name in enumerate(names):
        caps[labels.index(name)] = cc[i]
        seen.append(name)
        for j in range(len(names)):
            need(multiplicity*bg[i][j] == (gram[labels.index(name)] if i == j else 0),
                 'same full23 source Gram normalization')
            gram_entries += 1
    square = multiplicity*sum(abs(bg[i][j])*cc[i]*cc[j]
                             for i in range(len(names)) for j in range(len(names)))
    need(square == Q(block['squared_norm_upper']), 'source block norm square')
    norm_square += square
FH = Q(bundle['F_H_upper'])
need(sorted(seen) == sorted(labels) and len(seen) == 23, 'all23 physical source components exactly once')
need(norm_square == Q(bundle['squared_upper']) == Q(local['source_H_norm_square']) ==
     sum(w*c*c for w,c in zip(gram,caps)), 'complete actual raw source norm')
need(FH == Q(local['genuine_pointwise_affine_F_H_cap']) == 17387513156573 and
     (FH-1)**2 <= norm_square < FH**2, 'genuine FH, not response SOURCE')
groups.extend(['seven physical source blocks', '93 source Gram normalization entries',
               'complete positive source norm', 'genuine FH integer enclosure'])

# Unit-width kernel and its exact moments.
k = list(map(Q,(0,0,30,-60,30)))
dk = derivative(k)
ddk = derivative(dk)
need(integral(k) == 1 and integral([Q(0)]+k) == Q(1,2) and
     integral([Q(0),Q(0)]+k) == Q(2,7), 'kernel exact mass, first and second moments')
need(k == [Q(0),Q(0),Q(30),Q(-60),Q(30)] and
     dk == [Q(0),Q(60),Q(-180),Q(120)], 'positive kernel and derivative factorizations')
need(all(value(p,x) == 0 for p in (k,dk) for x in (Q(0),Q(1))),
     'compact C1 kernel matches zero value and slope')
peak = value(k,Q(1,2))
variation = integral(dk,Q(0),Q(1,2))-integral(dk,Q(1,2),Q(1))
need(peak == Q(15,8) and variation == 2*peak == Q(15,4), 'sharp first derivative L1 norm')
# k''=60(1-6u+6u^2), roots (3+-sqrt(3))/6. The extrema
# of k' have modulus 10 sqrt(3)/3, so its total variation is
# 40 sqrt(3)/3 < 24, checked without floating-point square roots.
need(ddk == [Q(60),Q(-360),Q(360)] and Q(3) < Q(9,5)**2 and
     Q(40,3)*Q(9,5) == 24, 'outward second derivative L1 cap')
groups.extend(['exact kernel moments', 'nonnegative kernel factorizations',
               'kernel endpoint C1 matching', 'sharp first derivative variation',
               'outward second derivative variation'])

eps = Q(args.epsilon)
need(eps > 0, 'positive filter width')
epsilon_rows = []
for e in sorted(set((eps,Q(1,64),Q(1,1024),Q(1,65536))),reverse=True):
    epsilon_rows.append({'epsilon':e,'F_epsilon_H_cap':FH,
        'source_discrepancy_H_cap':2*FH,
        'F_epsilon_x_H_cap':variation*FH/e,
        'F_epsilon_xx_H_strict_upper':24*FH/e**2,
        'source_discrepancy_initial_H_cap':FH})
groups.append('positive-width scaled filter derivative caps')

# Symbolic slacks certify all past vertical points, not only mesh points.
# For x'=x-s >=0: R-x'-r=(R-x-r)+s,
# r-x'+5/4=(r-x+5/4)+s. All summands are nonnegative on D+.
T,R,H = Q(3,2),Q(7,4),Q(5,4)
domain_slack_coefficients = {'outer_past_slack':[Q(1),Q(1)],
                           'inner_past_slack':[Q(1),Q(1)],
                           'support_past_sum_slack':[Q(1),Q(1)]}
need(all(min(row) >= 0 for row in domain_slack_coefficients.values()) and
     0 < H < T < R and T <= 3 and R <= 8,
     'past closure, support closure, inherited source domain')
groups.extend(['entire D+ past vertical closure', 'inherited exact zero wedge preservation'])

b = Q(transfer['b_actual_boundary_X_operator_cap'])
Pcap = Q(transfer['P_boundary_inverse_X_cap'])
oldE = Q(transfer['best_reproved_global_C1_error_X_cap'])
readableE = Q(transfer['best_reproved_global_C1_readable_strict_upper'])
need(0 < b < 1 and Pcap == 1/(1-b), 'same actual Robin inverse')
need(oldE == Q(local['global_one_update_error_exact_retained_cap']) and
     readableE == Q(local['global_one_update_error_readable_strict_upper']) == Q('3.541e25'),
     'same exact re-established W1_quadratic majorant')
delta_X = T*2*FH*Pcap
newE = oldE+delta_X
need(delta_X == 3*FH/(1-b) and delta_X < Q('5.217e13') and newE < readableE,
     'filtered new candidate retains same strict readable C1 upper certificate')
groups.extend(['actual Robin inverse unchanged', 'exact new-candidate C1 perturbation charge',
               'same readable C1 certificate with explicit charge'])

# At x=0 the causal integral and its first time derivative vanish.
# k(s)=O(s^2), k'(s)=O(s), k'' locally bounded imply F_eps,F_eps,x,
# F_eps,xx ->0 as x->0+. The zero-data principal response has
# Delta_xx(0)=F_eps(0)-F(0)=-F(0), while its other second initial
# jets Delta_xr and Delta_rr vanish. This is distinct from a target
# C2 error estimate.
initial_delta_coefficients = {'value':Q(0),'time':Q(0),'radial':Q(0),
    'time_radial':Q(0),'radial_radial':Q(0),'time_time_per_F0':Q(-1)}
need(initial_delta_coefficients['time_time_per_F0'] == -1 and
     all(initial_delta_coefficients[name] == 0 for name in
         ('value','time','radial','time_radial','radial_radial')),
     'new initial acceleration discrepancy retained')
groups.extend(['zero initial filtered source and time jets',
               'new nonzero-permitted initial acceleration discrepancy',
               'flat physical shell corner compatibility'])

negative_control('kernel factor15 misnormalization', integral([z/2 for z in k]) == 1)
negative_control('peak substituted for derivative variation', variation <= peak)
negative_control('wrong first moment1/3', integral([Q(0)]+k) == Q(1,3))
negative_control('response SOURCE substituted for raw FH', FH == Q(52162547172519,2))
negative_control('omitted actual Robin inverse', delta_X <= 3*FH)
negative_control('old zero acceleration discrepancy reused',
                 initial_delta_coefficients['time_time_per_F0'] == 0)
# A future sample fails at the outer face and can activate F inside the
# original support wedge. These are exact geometry counterexamples.
x,r,s = Q(1,2),R-Q(1,2),Q(1,1024)
negative_control('future fixed-radius sample declared in D+', r <= R-(x+s))
support_sum, future_shift = 1-eps/4,eps/2
negative_control('future filter declared zero-wedge preserving', support_sum+future_shift <= 1)

for path,digest in pinned:
    need(hashlib.sha256(path.read_bytes()).hexdigest() == digest, 'input unchanged:'+str(path))
groups.append('all authenticated inputs unchanged on readback')
report = {'status':'PASS','new_candidate':'Z_epsilon, distinct from W1_quadratic and archived W1',
    'definition':'Z_epsilon=P[Bs_quadratic+K(La W0_quadratic+F_epsilon)+LG]',
    'fixed_target':'same actual U, F, G, coefficients, actual M, and first field',
    'epsilon':eps,'kernel_unit_coefficients':k,'kernel_mass':integral(k),
    'kernel_first_moment_per_epsilon':Q(1,2),'kernel_second_moment_per_epsilon_squared':Q(2,7),
    'kernel_derivative_L1_per_inverse_epsilon':variation,
    'kernel_second_derivative_L1_strict_upper_per_inverse_epsilon_squared':Q(24),
    'zero_extension':'F(x,r)=0 for x<0 before past-time convolution',
    'coordinate_scope':'component smoothing in registered scaled time and frame; not covariant averaging',
    'actual_full23_F_H_cap':FH,'complete_source_norm_square':norm_square,
    'source_Gram_entries_checked':gram_entries,'epsilon_comparison':epsilon_rows,
    'domain':{'T':T,'R':R,'inner_offset':H,
        'description':'0<=x<=T, max(0,x-H)<=r<=R-x',
        'past_nonnegative_time_available':True,'source_zero_wedge_preserved':True,
        'slack_sum_coefficients':domain_slack_coefficients},
    'actual_Robin_b':b,'actual_boundary_inverse_X_cap':Pcap,
    'new_candidate_minus_W1_quadratic_X_cap':delta_X,
    'new_candidate_minus_W1_quadratic_X_readable_strict_upper':Q('5.217e13'),
    'old_exact_target_error_X_cap':oldE,'new_exact_target_error_X_cap':newE,
    'new_target_error_X_readable_strict_upper':readableE,
    'initial_delta_jets_per_F0':initial_delta_coefficients,
    'target_error_forcing_initial':'r_epsilon(0,r)=F(0,r), not assumed zero',
    'complete_coupling_before_reflection':True,
    'old_point_values_exactly_preserved':False,
    'scope':'conditional auxiliary candidate regularization and C1 perturbation bound; no quantitative source approximation rate, actual F_x bound, evolved C2 target-error certificate, candidate second-jet theorem by this checker alone, natural Einstein residual theorem, exact nonlinear family or Big Bang-origin evidence',
    'input_sha256':{'higher_transfer':pinned[0][1],'local_bulk':pinned[1][1],
        'actual_full_source':pinned[2][1],'complete_geometric_map':pinned[3][1]},
    'named_check_groups':groups,'named_check_group_count':len(groups),
    'wrong_formula_controls_rejected':controls,'wrong_formula_control_count':len(controls)}
print(json.dumps(encode(report),sort_keys=True))
