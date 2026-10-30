#!/usr/bin/env python3
"""Continuous actual full23 scaled-time coefficient derivative ledger.

Differentiation is analytic before interval substitution. No interval endpoint,
frozen matrix, source magnitude cap, or numerical stencil is differentiated.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
DRIVER_PIN = 'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73'
PRIMITIVE_PIN = 'a008be0d89e2ac44940c273cfd2cbc5ba0d09de8fc172d5d81ef29121ab60189'
PROFILE_PIN = '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return raw


def sqrt_up(q):
    need(q >= 0, 'positive norm square')
    if not q:
        return Q(0)
    g = 1 << 112
    n = math.isqrt(q.numerator*g*g//q.denominator)
    if Q(n, g)**2 < q:
        n += 1
    out = Q(n, g)
    need(out*out >= q, 'outward rational square root')
    return out


def encode(z):
    if isinstance(z, Q):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (list, tuple)):
        return [encode(v) for v in z]
    return z


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
args = ap.parse_args()
root = args.repo_root.resolve()
base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
paths = {'coefficient_driver': base/'compute/full23_response.py',
         'primitive_operator': root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py',
         'profiles': base/'compute/response-normal/RESPONSE_PROFILES.json'}
pins = {'coefficient_driver': DRIVER_PIN, 'primitive_operator': PRIMITIVE_PIN, 'profiles': PROFILE_PIN}
raw = {k: pin(v, pins[k]) for k, v in paths.items()}
profiles = json.loads(raw['profiles'])
spec = importlib.util.spec_from_file_location('timejet_authenticated_response', paths['coefficient_driver'])
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
response = mod.Response(root)  # authenticates seven complete scientific roots
ctx, I = response.context, response.I
ell, weights = response.ell, list(response.weights)
need(ell == Q(1, 10**6), 'registered scaled time t=ell*x')
need(response.labels == profiles['labels'] and weights == list(map(Q, profiles['gram'])) and len(weights) == 23,
     'same complete full23 space and Gram')
need(min(weights) > 0 and sorted(i for b in response.blocks for i in b) == list(range(23)), 'positive orthogonal block partition')
m = response.primitive


def norm_blocks(matrix):
    for i in range(23):
        for j in range(23):
            if response.block_of[i] != response.block_of[j]:
                need(matrix[i][j].lo == matrix[i][j].hi == 0, 'time derivative degree-block closure')
    squares = [sum(weights[i]/weights[j]*max(abs(matrix[i][j].lo), abs(matrix[i][j].hi))**2
                   for i in b for j in b) for b in response.blocks]
    caps = [sqrt_up(z) for z in squares]
    need(all(c*c >= z for c, z in zip(caps, squares)), 'block weighted Frobenius bound')
    return squares, caps


def enclose(x, r):
    need(0 <= x.lo <= x.hi <= Q(3, 2) and 0 <= r.lo <= r.hi <= Q(7, 4), 'certified coefficient rectangle')
    bg = ctx.actual_background(r)
    rho, a, p = [bg[k] for k in ('rho', 'a', 'p')]
    need(rho.lo > 0 and p.hi < 0, 'actual background signs')
    t = ell*x
    # For t>=0, tanh(t)<=t and tanh(t)>=t-t^3/3; sech(t)<=1
    # and sech(t)>=1-t^2/2. All lower endpoints are positive here.
    tau = I(t.lo-t.lo**3/3, t.hi)
    sech = I(1-t.hi**2/2, 1)
    need(0 <= tau.lo <= tau.hi < 1 and 0 < sech.lo <= sech.hi <= 1, 'analytic hyperbolic ranges')
    s2 = sech*sech
    coefficients = {
        'Delta': -2*ell**3*s2*tau,
        'Ca': -2*ell**3*rho*a*sech*tau,
        'Ch': 2*ell**3*sech*(s2-tau*tau),
        'Sah': ell**3*rho*a*s2,
        'O': -3*ell**3*rho*a*s2,
        'Sh': 2*ell**3*tau*s2,
        'W_L': 3*ell**3*rho*p*s2,
        'Ephi': -ell**3*rho*p*sech*tau,
    }
    dx = response.matrix_sum(*[(value, m[key]) for key, value in coefficients.items()])
    # Keep this sharp analytic interval as exact Fractions instead of passing
    # its lower endpoint through the inherited dyadic interval constructor.
    # The latter is valid but would round -3ell^2 slightly farther outward.
    ax_scalar = [-3*ell**2, -3*ell**2*sech.lo**2]
    ax = [[ax_scalar[:] if i == j else [Q(0), Q(0)] for j in range(23)] for i in range(23)]
    zero = response.matrix_sum((I(0), m['I']))
    squares, caps = norm_blocks(dx)
    need(all(ax[i][j] == [0, 0] for i in range(23) for j in range(23) if i != j), 'scalar-identity A_x')
    # A_x is a shared scalar times I, so its operator bound is the scalar
    # magnitude, with no dimension-dependent Frobenius factor.
    acap = max(map(abs, ax_scalar))
    need(acap <= 3*ell**2, 'sharp global A_x cap')
    need(all(v.lo == v.hi == 0 for row in zero for v in row), 'exact B_x and M_x zero')
    return {
        'x_interval': [x.lo, x.hi], 'r_interval': [r.lo, r.hi],
        'rho_interval': [rho.lo, rho.hi], 'a_interval': [a.lo, a.hi], 'p_interval': [p.lo, p.hi],
        'tanh_interval': [tau.lo, tau.hi], 'sech_interval': [sech.lo, sech.hi],
        'D_x_primitive_coefficients': {key: [value.lo, value.hi] for key, value in coefficients.items()},
        'matrix_intervals': {'A_x': ax, **{key: [[[v.lo, v.hi] for v in row] for row in val]
                             for key, val in (('B_x', zero), ('D_x', dx), ('M_x', zero))}},
        'D_x_block_weighted_Frobenius_squared_upper': squares,
        'D_x_block_H_operator_caps': caps,
        'H_operator_caps': {'A_x': acap, 'B_x': Q(0), 'D_x': max(caps), 'M_x': Q(0)},
    }


global_box = enclose(I(0, Q(3, 2)), I(0, Q(7, 4)))
point_box = enclose(I(Q(9, 8)), I(0))
need(global_box['H_operator_caps']['A_x'] == Q(3, 10**12), 'whole-rectangle A_x maximum at x=0')
commutator = global_box['H_operator_caps']['A_x'] + Q(3, 4)*global_box['H_operator_caps']['D_x']
readable_Dx = Q('2.354e-16')
readable_commutator = Q('3.000177e-12')
need(0 < global_box['H_operator_caps']['D_x'] < readable_Dx and
     0 < commutator < readable_commutator, 'strict readable whole-domain derivative ledger')
ctx.reauthenticate()
need(mod.authenticate_local() == response.local_before, 'local dependency readback')
for key, path in paths.items():
    pin(path, pins[key])
print(json.dumps(encode({'status': 'PASS', 'coordinate_derivative': 'partial_x at fixed inward scaled r; t=ell*x, partial_x=ell*partial_t',
    'derivation': 'analytic product rule on the actual homogeneous coefficient assembly before interval substitution; radial background fixed in x',
    'input_sha256': pins, 'labels': list(response.labels), 'Gram': weights,
    'block_names': list(mod.BLOCK_NAMES), 'orthogonal_blocks': response.blocks,
    'domain': 'D+:0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x; bounded using containing rectangle',
    'global': global_box, 'selected_shell_point': point_box,
    'coefficient_commutator_X_to_H_upper': commutator,
    'readable_caps': {'A_x_inclusive_upper': Q('3e-12'), 'B_x_exact': Q(0),
                      'D_x_strict_upper': readable_Dx, 'M_x_exact': Q(0),
                      'coefficient_commutator_strict_upper': readable_commutator},
    'commutator_meaning': '||(A_x e_x+B_x e_r+D_x e)||_H <= coefficient_commutator_X_to_H_upper*X(e), using inherited zero-wedge value bound ||e||_H<=(3/4)X(e)',
    'source_derivative_supplied': False,
    'scope': 'actual full23 homogeneous time-coefficient derivatives only; no complete bulk-source derivative, second-initial-jet error, quantitative C2 theorem, Einstein/constraint residual, nonlinear family or physical observable'},), sort_keys=True))
