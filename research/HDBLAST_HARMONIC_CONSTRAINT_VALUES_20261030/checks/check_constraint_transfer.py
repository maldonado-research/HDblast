#!/usr/bin/env python3
"""Bound the linear harmonic covector of the registered candidate error.

The nonlinear quadratic gauge source cancels only under the same-fixed-source
and formal subsidiary completion premises stated in the report. This script
does not infer those premises from the norm certificate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
DRIVER_PIN = 'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73'
BOUND_PIN = '321a9fc2917f1d6680dab825a087c03858ebf39de7e3fde22507df1dfc4b957f'
PROFILE_PIN = '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'
FORMULA_PIN = '8ac9f891973b2f592583eab49ea9b68c2b33b5606bfde928d45c7e580f87ed0e'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return raw


def sqrt_up(z):
    need(z >= 0, 'nonnegative norm square')
    g = 1 << 112
    n = math.isqrt(z.numerator*g*g//z.denominator)
    if F(n, g)**2 < z:
        n += 1
    out = F(n, g)
    need(out*out >= z, 'outward rational square root')
    return out


def encode(z):
    if isinstance(z, F):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (list, tuple)):
        return [encode(v) for v in z]
    return z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--bulk-capsule', type=Path)
    args = ap.parse_args()
    root = args.repo_root.resolve()
    base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
    bulk = args.bulk_capsule or root/'research/HDBLAST_ACTUAL_BULK_UPDATE_20261030'
    driver = base/'compute/full23_response.py'
    profile_path = base/'compute/response-normal/RESPONSE_PROFILES.json'
    bound_path = bulk/'reference/BULK_BOUND.json'
    formula_path = root/'research/HDBLAST_COMMON_GAUGE_FORCING_20261016/evolution-review/EVOLUTION_REVIEW.md'
    pin(driver, DRIVER_PIN)
    profile = json.loads(pin(profile_path, PROFILE_PIN))
    bound = json.loads(pin(bound_path, BOUND_PIN))
    pin(formula_path, FORMULA_PIN)
    spec = importlib.util.spec_from_file_location('constraint_authenticated_response', driver)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    response = mod.Response(root)
    ctx, I = response.context, response.I
    labels, weights = list(response.labels), list(response.weights)
    need(labels == profile['labels'] == bound['labels'] and weights == list(map(F, profile['gram'])), 'registered full23 metric norm')
    need(len(labels) == 23 and min(weights) > 0 and bound['status'] == 'PASS', 'successful inherited full23 certificate')
    ell = response.ell
    need(ell == F(1, 10**6), 'registered coordinate scale')
    E = F(bound['one_update_readable_error_strict_upper'])
    need(E == 3541*10**22 and F(bound['one_update_retained_target_error_X_cap']) < E, 'same readable strict global C1 error')
    idx = {z: i for i, z in enumerate(labels)}
    output = [f'{ch}:{l}' for ch in ('R', 'T') for l in (0, 2, 4)] + ['A:2', 'A:4']
    oid = {z: i for i, z in enumerate(output)}
    S = {0: F(1), 2: F(1, 12), 4: F(19, 720)}
    outweights = [S[int(z.split(':')[1])]*(int(z.split(':')[1])*(int(z.split(':')[1])+2) if z.startswith('A:') else 1) for z in output]
    Rx, Rr = [[[F(0)]*23 for _ in output] for _ in range(2)]
    def put(m, out, inc, degree, value):
        if f'{inc}:{degree}' in idx:
            m[oid[f'{out}:{degree}']][idx[f'{inc}:{degree}']] += value
    for degree in (0, 2, 4):
        for ch, c in [('K', -F(1, 2)), ('N', -F(1, 2)), ('J', F(3, 2))]:
            put(Rr, 'R', ch, degree, c)
        put(Rx, 'R', 'L', degree, -1)
        put(Rr, 'T', 'L', degree, -1)
        for ch in ('K', 'N', 'J'):
            put(Rx, 'T', ch, degree, -F(3 if ch == 'J' else 1, 2))
        if degree:
            put(Rr, 'A', 'v', degree, -1)
            put(Rx, 'A', 'u', degree, -1)

    # Whitened combined derivative Gram. The cross R/T entries cancel only
    # after retaining the metric's correct mixed-component weight and signs.
    joint = [[sum((Rx[i][k]*Rx[j][k]+Rr[i][k]*Rr[j][k])/weights[k] for k in range(23)) for j in range(8)] for i in range(8)]
    for i in range(8):
        for j in range(8):
            want = (F(1) if output[i].startswith('A:') else F(7, 4))/outweights[i] if i == j else 0
            need(joint[i][j] == want, 'exact combined derivative Gram')
    derivative_cap = sqrt_up(F(7, 2))
    PN = [not z.startswith(('L:', 'v:')) for z in labels]
    shell_rows = [i for i, z in enumerate(output) if not z.startswith('R:')]
    for matrix, project, norm2 in ((Rx, PN, F(5, 4)), (Rr, [not z for z in PN], F(1, 2))):
        for i in shell_rows:
            for j in shell_rows:
                z = sum(matrix[i][k]*matrix[j][k]/weights[k] for k in range(23) if project[k])
                need(z == 0 if i != j else z*outweights[i] <= norm2, 'shell derivative operator bound')
    shell_derivative_cap = sqrt_up(F(5, 4)) + sqrt_up(F(1, 2))

    # Check the exact normal-gauge cancellation by independent symbolic
    # coefficients. Interval entry boxes cannot prove a shared-symbol identity.
    zero = [F(0)]*23
    for degree in (0, 2, 4):
        rr = Rr[oid[f'R:{degree}']]
        va = [F(0)]*23
        for ch, c in [('K', 4), ('N', 1), ('J', -3)]:
            va[idx[f'{ch}:{degree}']] = F(c)
        for name in ('By_a', 'By_p', 'By_sigma1', 'By_sigma2'):
            result = [sum(rr[k]*response.primitive[name][k][j] for k in range(23)) for j in range(23)]
            need(result == (va if name == 'By_a' else zero), 'exact shared-symbol shell normal-gauge cancellation')

    def background(x, r):
        state = ctx.actual_background(r)
        need(state['rho'].lo > 0, 'positive actual warp factor')
        t = ell*x
        ht = I(t.lo-t.lo**3/3, t.hi)
        sech = I(1-t.hi**2/2, 1)
        need(sech.lo > 0, 'positive exact sech enclosure')
        return state, ell*state['rho']*state['a'], ell*ht, ell*sech

    def value_matrix(A, T, Z):
        mat = [[I(0) for _ in labels] for _ in output]
        def pv(out, inc, degree, value):
            if f'{inc}:{degree}' in idx:
                i, j = oid[f'{out}:{degree}'], idx[f'{inc}:{degree}']
                mat[i][j] = mat[i][j] + value
        for degree in (0, 2, 4):
            lam = degree*(degree+2)
            for ch, c in [('K', 4), ('N', 1), ('J', -3)]:
                pv('R', ch, degree, c*A)
            pv('R', 'L', degree, -3*T)
            pv('R', 'v', degree, -lam*Z)
            pv('T', 'L', degree, 5*A)
            for ch in ('N', 'J'):
                pv('T', ch, degree, -3*T)
            pv('T', 'u', degree, -lam*Z)
            if degree:
                pv('A', 'v', degree, 5*A)
                pv('A', 'u', degree, -4*T)
                for ch, c in [('K', -F(1, 2)), ('N', F(1, 2)), ('J', -F(1, 2)), ('H', 2-F(2, 3)*lam)]:
                    pv('A', ch, degree, c*Z)
        return mat

    def block_caps(mat, shell=False):
        caps = {}
        for degree in (0, 2, 4):
            rows = [i for i, z in enumerate(output) if z.endswith(':'+str(degree)) and (not shell or not z.startswith('R:'))]
            columns = [j for j, z in enumerate(labels) if z.endswith(':'+str(degree)) and (not shell or PN[j])]
            need(all(mat[i][j].lo == mat[i][j].hi == 0 for i in rows for j in range(23) if j not in columns), 'zero-order degree/projector closure')
            square = sum(outweights[i]*max(abs(mat[i][j].lo), abs(mat[i][j].hi))**2/weights[j] for i in rows for j in columns)
            caps[str(degree)] = sqrt_up(square)
        return caps

    full_state, A, T, Z = background(I(0, F(3, 2)), I(0, F(7, 4)))
    full_matrix = value_matrix(A, T, Z)
    full_caps = block_caps(full_matrix)
    full_zero = max(full_caps.values())
    full_scaled_constant = derivative_cap + F(3, 4)*full_zero
    full_scaled_cap = full_scaled_constant*E
    full_frame_cap = full_scaled_cap/(ell*full_state['rho'].lo)
    point_state, Ap, Tp, Zp = background(I(F(9, 8)), I(0))
    shell_matrix = value_matrix(I(0), Tp, Zp)
    # The discarded a terms multiply only zero PD traces or the exactly
    # canceled radial constraint. This is not an actual-a midpoint replacement.
    shell_caps = block_caps(shell_matrix, True)
    shell_zero = max(shell_caps.values())
    shell_scaled_constant = shell_derivative_cap + F(1, 8)*shell_zero
    shell_scaled_cap = shell_scaled_constant*E
    shell_frame_cap = shell_scaled_cap/(ell*point_state['rho'].lo)
    readable = {'global_scaled': F('6.625e25'), 'global_frame': F('5.100e29'),
                'shell_scaled': F('6.463e25'), 'shell_frame': F('4.975e29')}
    for name, value in (('global_scaled', full_scaled_cap), ('global_frame', full_frame_cap),
                        ('shell_scaled', shell_scaled_cap), ('shell_frame', shell_frame_cap)):
        need(0 < value < readable[name], 'strict readable constraint-value upper bound')
    ctx.reauthenticate()
    need(mod.authenticate_local() == response.local_before, 'local dependency readback')
    for p, h in ((driver, DRIVER_PIN), (profile_path, PROFILE_PIN), (bound_path, BOUND_PIN), (formula_path, FORMULA_PIN)):
        pin(p, h)
    result = {'status': 'PASS', 'input_sha256': {'coefficient_driver': DRIVER_PIN, 'global_error_report': BOUND_PIN, 'profiles': PROFILE_PIN, 'radial_formula_note': FORMULA_PIN},
        'constraint': 'mathfrak_C^(2)[W1]=DC[W1-U], conditional on same fixed actual first field/Q_C and exact formal subsidiary completion mathfrak_C^(2)[U]=0',
        'norm': 'positive L2(S3) comparison norm of radially parallel frame covector components, normalized unit-sphere average; not Lorentz invariant',
        'scaled_constraint': 'ell*rho*DC', 'input_labels': labels, 'input_Gram': weights, 'output_labels': output, 'output_Gram': outweights,
        'derivative_matrices': {'Rx': Rx, 'Rr': Rr}, 'combined_derivative_Gram': joint,
        'global': {'domain': 'D+:0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x',
            'containing_rectangle': ['0','3/2','0','7/4'], 'rho_interval': [full_state['rho'].lo, full_state['rho'].hi],
            'ell_rho_a_interval': [A.lo,A.hi], 'ell_tanh_interval': [T.lo,T.hi], 'ell_sech_interval': [Z.lo,Z.hi],
            'value_matrix_intervals': [[[v.lo,v.hi] for v in row] for row in full_matrix],
            'zero_order_block_H_operator_caps': full_caps, 'joint_derivative_constant_upper': derivative_cap,
            'value_error_factor': F(3,4), 'scaled_constraint_operator_constant_upper': full_scaled_constant,
            'scaled_constraint_H_cap': full_scaled_cap, 'frame_constraint_H_cap': full_frame_cap},
        'shell_point': {'point': ['9/8','0'], 'rho_interval': [point_state['rho'].lo,point_state['rho'].hi],
            'exact_radial_component_error': ['0','0'], 'tangential_value_matrix_intervals': [[[v.lo,v.hi] for v in row] for row in shell_matrix],
            'tangential_zero_order_block_H_operator_caps': shell_caps,
            'derivative_constant_upper': shell_derivative_cap, 'value_error_factor': F(1,8),
            'scaled_constraint_operator_constant_upper': shell_scaled_constant,
            'scaled_constraint_H_cap': shell_scaled_cap, 'frame_constraint_H_cap': shell_frame_cap},
        'same_global_error_readable_strict_upper': E,
        'readable_strict_constraint_H_caps': readable,
        'support_premise': 'same target/candidate vanish on inherited x+r<=1; global value error <=(3/4)E, shell value error <=(1/8)E, separately from magnitude-only source caps',
        'shell_boundary_premise': 'same fixed actual M and G; error has PN e_r=M PN e and PD e=PD e_x=0; PD e_r remains generic and is retained',
        'scope': 'conditional harmonic covector VALUE upper bound; no new direct Q_C evaluation, lower bound, relative smallness, constraint derivative, Hamiltonian/momentum residual, full Einstein residual, quantitative C2 certificate, nonlinear family or physical observable'}
    print(json.dumps(encode(result), sort_keys=True))


if __name__ == '__main__':
    main()
