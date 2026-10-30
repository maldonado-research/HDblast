#!/usr/bin/env python3
"""Independent rational endpoint/norm audit of actual time coefficient boxes."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
OP_PIN = 'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def iv(z):
    a, b = map(Q, z)
    need(a <= b, 'ordered rational interval')
    return a, b


def add(a, b):
    return a[0]+b[0], a[1]+b[1]


def neg(a):
    return -a[1], -a[0]


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    z = [x*y for x in a for y in b]
    return min(z), max(z)


def scale(a, b):
    return mul(a, (b, b))


def inside(a, b):
    return b[0] <= a[0] <= a[1] <= b[1]


def total(values):
    out = Q(0), Q(0)
    for z in values:
        out = add(out, z)
    return out


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--report', type=Path, required=True)
args = ap.parse_args()
raw = args.report.read_bytes()
d = json.loads(raw)
path = args.repo_root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py'
op_bytes = path.read_bytes()
need(hashlib.sha256(op_bytes).hexdigest() == OP_PIN, 'INPUT_HASH:primitive operator')
spec = importlib.util.spec_from_file_location('time_audit_primitive', path)
op = importlib.util.module_from_spec(spec)
spec.loader.exec_module(op)
m = op.matrices()
ell = Q(1, 10**6)
weights = list(map(Q, d['Gram']))
blocks = d['orthogonal_blocks']
need(d['status'] == 'PASS' and d['labels'] == op.LABELS and weights == op.gram(), 'same complete registered space')
need(sorted(i for block in blocks for i in block) == list(range(23)), 'orthogonal block partition')
audits = []
for domain in ('global', 'selected_shell_point'):
    b = d[domain]
    rho, a, p, tau, s = [iv(b[k]) for k in ('rho_interval', 'a_interval', 'p_interval', 'tanh_interval', 'sech_interval')]
    x = iv(b['x_interval'])
    tlo, thi = ell*x[0], ell*x[1]
    need(0 <= tlo <= thi < 1, 'small nonnegative physical time')
    need(tau[0] <= tlo-tlo**3/3 and tau[1] >= thi, 'analytic tanh bracket covered')
    need(s[0] <= 1-thi**2/2 and s[1] >= 1, 'analytic sech bracket covered')
    s2, tau2 = mul(s, s), mul(tau, tau)
    coeff = {
        'Delta': scale(mul(s2, tau), -2*ell**3),
        'Ca': scale(mul(mul(mul(rho, a), s), tau), -2*ell**3),
        'Ch': scale(mul(s, sub(s2, tau2)), 2*ell**3),
        'Sah': scale(mul(mul(rho, a), s2), ell**3),
        'O': scale(mul(mul(rho, a), s2), -3*ell**3),
        'Sh': scale(mul(tau, s2), 2*ell**3),
        'W_L': scale(mul(mul(rho, p), s2), 3*ell**3),
        'Ephi': scale(mul(mul(mul(rho, p), s), tau), -ell**3),
    }
    recorded = {key: iv(z) for key, z in b['D_x_primitive_coefficients'].items()}
    need(set(recorded) == set(coeff), 'all eight primitive coefficients')
    for key, z in coeff.items():
        need(inside(z, recorded[key]), 'recorded outward analytic coefficient:'+key)
    matrix = [[iv(z) for z in row] for row in b['matrix_intervals']['D_x']]
    for i in range(23):
        for j in range(23):
            plain = total(scale(recorded[key], m[key][i][j]) for key in recorded)
            need(inside(plain, matrix[i][j]), 'recorded outward complete matrix sum')
    squares = []
    for block, recorded_square, cap in zip(blocks, b['D_x_block_weighted_Frobenius_squared_upper'], b['D_x_block_H_operator_caps']):
        square = sum(weights[i]/weights[j]*max(abs(v) for v in matrix[i][j])**2 for i in block for j in block)
        need(square == Q(recorded_square), 'independent exact weighted Frobenius square')
        need(square <= Q(cap)**2, 'independent outward block cap')
        squares.append(str(square))
    cap = Q(b['H_operator_caps']['D_x'])
    need(cap == max(map(Q, b['D_x_block_H_operator_caps'])), 'orthogonal sum operator cap maximum')
    for key in ('B_x', 'M_x'):
        need(all(iv(z) == (0, 0) for row in b['matrix_intervals'][key] for z in row), 'exact zero matrix '+key)
        need(Q(b['H_operator_caps'][key]) == 0, 'exact zero cap '+key)
    ax = [[iv(z) for z in row] for row in b['matrix_intervals']['A_x']]
    expected_ax = (-3*ell**2, -3*ell**2*s[0]**2)
    need(all(ax[i][j] == (expected_ax if i == j else (0, 0)) for i in range(23) for j in range(23)), 'exact analytic scalar A_x boxes')
    need(Q(b['H_operator_caps']['A_x']) == 3*ell**2, 'sharp dimension-free A_x operator cap')
    audits.append({'domain': domain, 'weighted_Frobenius_squares': squares, 'matrix_entry_containments': 529,
                   'coefficient_containments': 8, 'A_x_operator_cap': str(3*ell**2), 'D_x_operator_cap': str(cap)})
global_caps = d['global']['H_operator_caps']
comm = Q(global_caps['A_x'])+Q(3, 4)*Q(global_caps['D_x'])
need(comm == Q(d['coefficient_commutator_X_to_H_upper']) and comm < Q('3.000177e-12'), 'coefficient commutator with zero-wedge value factor')
need(Q(global_caps['D_x']) < Q('2.354e-16'), 'readable complete D_x norm')
need(path.read_bytes() == op_bytes, 'primitive readback unchanged')
print(json.dumps({'status': 'PASS_INDEPENDENT_TIME_COEFFICIENT_REPORT_AUDIT',
    'coefficient_report_sha256': hashlib.sha256(raw).hexdigest(), 'primitive_sha256': OP_PIN,
    'audits': audits, 'weighted_Frobenius_blocks_checked': 10,
    'coefficient_commutator_X_to_H_upper': str(comm),
    'scope': 'independent rational interval containment and exact block-norm recomputation; underlying actual background certificate remains inherited; no actual F_x, C2 error or physical residual established'}, sort_keys=True))
