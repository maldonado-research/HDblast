#!/usr/bin/env python3
"""Inherited actual affine full23 shell G enclosure at x=9/8, r=0.
A uniform shell amplitude cap is used; the source is not evaluated pointwise.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest,
         'INPUT_HASH:' + str(path))
    return json.loads(raw)


def weight(label):
    if label.startswith('T:'):
        return F(1, 6)
    field, degree = label.split(':')
    degree = int(degree)
    value = {0: F(1), 2: F(1, 12), 4: F(19, 720)}[degree]
    if field == 'L': value *= 2
    if field == 'J': value *= 3
    if field in ('v', 'u'): value *= 2*degree*(degree+2)
    if field == 'H': value *= F(2, 3)*degree*(degree+2)*(degree*(degree+2)-3)
    return value


parser = argparse.ArgumentParser()
parser.add_argument('--repo-root', type=Path, required=True)
args = parser.parse_args()
root = args.repo_root.resolve()/'research'
shell_path = root/'HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json'
full_path = root/'HDBLAST_FULL_ANGULAR_RESPONSE_20261026/shell/verified/FULL23_SHELL_NORMAL.json'
shell_hash = '559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a'
full_hash = '0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe'
shell = read(shell_path, shell_hash)
full = read(full_path, full_hash)
need(shell['status'] == 'PASS_CONTINUOUS_FULL_COMMON_FIRST_HARMONIC_PHYSICAL_SHELL_FORCING', 'source scope')
need(shell['domain']['scaled_shell'] == '0<=x<=3,r=0; entire unit S3', 'uniform source domain')
need(full['status'] == 'PASS' and full['robin_dimension'] == 18
     and full['universal_boundary_source_rank'] == 8, 'full23 packing scope')
labels = [f'{field}:{degree}' for field in ('K','L','N','J','w')
          for degree in (0,2,4)] + [f'{field}:{degree}' for field in ('v','u','H')
                                  for degree in (2,4)] + ['T:+','T:-']
need(full['labels'] == labels, 'registered labels')
B = [[F(value) for value in row] for row in full['source_basis_B_rows_in_LABELS_order']]
S = [[F(value) for value in row] for row in full['source_gram_Btranspose_H_B']]
need(len(B) == 23 and all(len(row) == 8 for row in B), 'source matrix shape')
need(len(S) == 8 and all(len(row) == 8 for row in S), 'source Gram shape')
weights = list(map(weight, labels))
for i in range(8):
    for j in range(8):
        need(S[i][j] == sum(weights[k]*B[k][i]*B[k][j] for k in range(23)),
             'source Gram entry:' + str((i,j)))
need(S[7][7] != sum((3*weights[k] if labels[k]=='T:+' else weights[k])
                    *B[k][7]*B[k][7] for k in range(23)), 'wrong TT weight control')
expected = {
    'alpha': 'a*h*h/4-sigma1*h*v/6-sigma2*v*v/6',
    'beta': 'a*h+sigma1*v/3',
    'gn_Y2': '4*sigma2*v*v/3+4*sigma1*h*v-10*a*h*h+8*a*h*tau-4*a*(4*X*X+12*L*L-64*L*d)',
    'gn_grad2': '8*a*V*V',
    'gn_hess2': '-16*a*d*d',
    'pi': 'p*h*h/4+sigma2*h*v/2+sigma3*v*v/2',
    'tau': '2*X+6*L-16*d'
}
need(full['source_profile_definitions'] == expected, 'registered source profiles')
need(full['source_profile_names'] == [
    '-ell*rho*gn_Y2','-ell*rho*gn_grad2','-ell*rho*gn_hess2',
    '2*ell*rho*(-alpha+2*beta*X)','ell*rho*pi','-2*ell*rho*beta*V',
    '2*ell*rho*(alpha-2*beta*L)','-4*ell*rho*beta*d'], 'registered source adapters')
P = shell['bounds']['parameters']
A = shell['bounds']['actual_field_amplitude_caps']
parameter = lambda key: F(P[key])
amplitude = lambda key: F(A[key]['exact_upper'])
ell = parameter('ell')
rho = max(map(F, P['rho_b_interval']))
a, sigma1, sigma2, sigma3, p = map(parameter,
    ('a_abs','sigma_prime_abs','sigma_second_abs','sigma_third','p_abs'))
h, v, X, L, V, d, tau = map(amplitude,
    ('J','v','X','Z','V_gauge','d','tau'))
need(ell == F(1,10**6) and min(rho,a,sigma1,sigma2,sigma3,p,h,v,X,L,V,d,tau)>=0,
     'nonnegative inherited magnitude caps')
# These are magnitudes of the eight ACTUAL profiles, not their x derivatives.
# Time is t=ell*x; G itself already contains exactly one ell from its adapter.
beta = a*h + sigma1*v/3
alpha = a*h*h/4 + sigma1*h*v/6 + sigma2*v*v/6
pi = p*h*h/4 + sigma2*h*v/2 + sigma3*v*v/2
g0 = (F(4,3)*sigma2*v*v + 4*sigma1*h*v + 10*a*h*h
      + 8*a*h*tau + 4*a*(4*X*X + 12*L*L + 64*L*d))
profiles = [ell*rho*g0,
            8*ell*rho*a*V*V,
            16*ell*rho*a*d*d,
            2*ell*rho*(alpha + 2*beta*X),
            ell*rho*pi,
            2*ell*rho*beta*V,
            2*ell*rho*(alpha + 2*beta*L),
            4*ell*rho*beta*d]
rows = [sum(abs(entry)*cap for entry,cap in zip(row, profiles)) for row in B]
norm2 = sum(abs(S[i][j])*profiles[i]*profiles[j] for i in range(8) for j in range(8))
need(min(profiles) >= 0 and norm2 >= 0, 'ordered source bounds')
cap = math.isqrt(norm2.numerator//norm2.denominator)
if cap*cap*norm2.denominator <= norm2.numerator:
    cap += 1
need(norm2 < F(cap)**2 and cap <= 487122, 'readable strict source H cap')
need(max(rows) < 1508337, 'readable strict modal cap')
need(all(rows[i] == 0 for i in (3,4,5,15,16)), 'exact Dirichlet source zeros')
need(rows[21] == rows[22] and rows[21] < F(4933,10**24), 'TT bound')
need(rows[12] < 94279, 'w0 source cap')
x = F(9,8)
need(0 <= x <= 3, 'point covered by uniform shell enclosure')
# A zero-collar derivative estimate is also valid when the inherited flat
# first-field collar proves G=0 up to x=1. Its readable H bound is worse here:
# (x-1)*58,807,538 > 487,122, so it is not substituted for this value bound.
need((x-1)*58807538 > cap, 'uniform source value cap is tighter than readable collar derivative cap')
print(json.dumps({
    'status':'PASS',
    'point':[str(x),'0'],
    'source':'actual registered affine full23 shell G',
    'method':'uniform inherited actual-field amplitude enclosure; no point-sign claim',
    'uniform_domain':'0<=x<=3,r=0',
    'input_sha256':{shell_path.name:shell_hash,full_path.name:full_hash},
    'source_profile_abs_caps':list(map(str,profiles)),
    'source_H_norm_square_upper':str(norm2),
    'source_H_norm_strict_integer_upper':str(cap),
    'source_max_modal_strict_integer_upper':'1508337',
    'rows':[{'label':label,'abs_upper':str(value),
             'interval':[str(-value),str(value)]} for label,value in zip(labels,rows)],
    'source_gram_entries_checked':64,
    'dirichlet_rows_exactly_zero':5,
    'scope':'conditional on authenticated inherited first-field/background bounds and source factorization; not independent PDE evolution, pointwise actual source evaluation, force, or physical observable'
}, sort_keys=True))
