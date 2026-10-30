#!/usr/bin/env python3
"""Audit the actual initial radial curvature error at the retained endpoint.

The actual initial jet is computed by an autonomous ODE recurrence at an
authenticated actual state enclosure, independently of the paired Taylor-jet
implementation.  This is a coordinate second-initial-derivative calculation,
not a spacetime C2 estimate or an Einstein residual.
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

PINS = {
    'producer': '2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1',
    'initial_report': 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad',
    'bump_report': 'add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd',
    'old_smoothing_checker': 'c830e0d09f1e522c4945e589fb9ca5912fb6bc8572c5d653144db0fede566a71',
}


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    data = path.read_bytes()
    need(hashlib.sha256(data).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return data


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0) for i in range(n)])


def pscale(p, c):
    return trim([c*x for x in p])


def pder(p):
    return trim([i*x for i, x in enumerate(p) if i] or [Q(0)])


def pmul(a, b):
    out = [Q(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def shifted(p, c):
    return [sum((p[k]*math.comb(k, j)*c**(k-j) for k in range(j, len(p))), Q(0)) for j in range(len(p))]


def integral_01(p):
    return sum((x/Q(i+1) for i, x in enumerate(p)), Q(0))


def encode(z):
    if isinstance(z, Q):
        return str(z)
    if isinstance(z, dict):
        return {k: encode(v) for k, v in z.items()}
    if isinstance(z, (list, tuple)):
        return [encode(v) for v in z]
    return z


def algebra_checks(polys):
    groups, controls = [], []
    def checked(name, ok):
        need(ok, name)
        groups.append(name)
    checked('bump_base_polynomial', polys[0] == [Q(1)])
    for n in range(6):
        difference = padd(polys[n], pscale(pder(polys[n]), -1))
        factor = [Q(0), Q(0), Q(1)] if n % 2 == 0 else [Q(0), Q(-4), Q(1)]
        want = pmul(factor, difference)
        if n % 2:
            want = padd(want, pscale(polys[n], -2))
        checked('bump_parity_recurrence_'+str(n+1), trim(polys[n+1]) == want)
    for n in range(7):
        monomial = [Q(0)]*(2*n)+[Q(1)]
        for sign in (-1, 1):
            # s>=4. Nonnegative coefficients after s=4+t prove both
            # s^(2n)+P_n(s)>=0 and s^(2n)-P_n(s)>=0 continuously.
            coefficients = shifted(padd(monomial, pscale(polys[n], sign)), Q(4))
            checked('bump_continuous_polynomial_envelope_'+str(n)+'_'+str(sign), min(coefficients) >= 0)
    kernel = [Q(35, 32)*x for x in (1, 0, -3, 0, 3, 0, -1)]
    half_mass = integral_01(kernel)
    half_first = integral_01([Q(0)]+kernel)
    half_second = integral_01([Q(0), Q(0)]+kernel)
    checked('positive_even_kernel_half_mass', half_mass == Q(1, 2))
    checked('positive_even_kernel_half_first_moment', half_first == Q(35, 256))
    checked('positive_even_kernel_half_second_moment', half_second == Q(1, 18))
    checked('kernel_full_moments', 2*half_first == Q(35, 128) and 2*half_second == Q(1, 9))
    # f(r)=r_+^2/2 is C1; piecewise f'''=0 but f'' jumps. Convolving at
    # zero gives half of the unit curvature, not either one-sided value.
    checked('C1_piecewise_quadratic_endpoint_half_curvature', half_mass == Q(1, 2))
    need(half_mass != 1, 'omitted exterior half mass must be rejected')
    controls.append('omitted_exterior_half_mass')
    need(half_mass > 0, 'zero piecewise third derivative cannot suppress a second-derivative jump')
    controls.append('piecewise_M3_alone_at_a_curvature_jump')
    checked('quadratic_extension_value_change_factor', half_second/2 == Q(1, 36))
    checked('quadratic_extension_first_change_factor', half_first == Q(35, 256))
    checked('quadratic_extension_second_change_factor', half_mass == Q(1, 2))
    need(half_first != Q(35, 128), 'full instead of half first moment must be rejected')
    controls.append('full_instead_of_half_first_moment')
    S = list(map(Q, (0,0,0,0,0,126,-420,540,-315,70)))
    deriv = S
    for j in range(5):
        at0, at1 = deriv[0], sum(deriv)
        checked('cutoff_join_order_'+str(j), (at0, at1) == ((Q(0), Q(1)) if j == 0 else (Q(0), Q(0))))
        deriv = pder(deriv)
    checked('cutoff_fifth_derivative_can_jump', deriv[0] == sum(deriv) == 15120)
    return {'groups': groups, 'group_count': len(groups), 'wrong_formula_controls': controls,
            'control_count': len(controls)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--robin-capsule', type=Path)
    args = ap.parse_args()
    root = args.repo_root.resolve()
    base = root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
    robin = args.robin_capsule or root/'research/HDBLAST_ROBIN_RECONSTRUCTION_20261030'
    paths = {
        'producer': base/'full23_initial.py',
        'initial_report': base/'results/INITIAL_REPORT.json',
        'bump_report': root/'research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016/numeric/inputs/BUMP_DERIVATIVE_BOUNDS.json',
        'old_smoothing_checker': robin/'checks/check_robin_resummed.py',
    }
    raw = {k: pin(path, PINS[k]) for k, path in paths.items()}
    initial, bump = [json.loads(raw[k]) for k in ('initial_report', 'bump_report')]
    need(initial['status'] == 'PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER', 'successful continuous initial inputs')
    need(initial['derivative_convention']['second_derivatives'][1] == 'partial_r^3 U', 'one-sided nominal third radial derivative caps')
    need(bump['status'] == 'PASS' and bump['maximum_derivative_order'] >= 6, 'fifth and sixth packet polynomials available')
    polys = [list(map(Q, row)) for row in bump['polynomial_coefficients_ascending']]
    controls = algebra_checks(polys)
    spec = importlib.util.spec_from_file_location('endpoint_authenticated_initial', paths['producer'])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    ctx = mod.Context(root)  # Authenticating all seven roots precedes scientific imports.
    need(ctx.data['bump'] == bump, 'same authenticated packet derivative table')
    I, J, ell = ctx.I, ctx.J, ctx.ell
    need(ell == Q(1, 10**6), 'registered coordinate scaling')
    R, delta, order = Q(7, 4), Q(1, 2**30), 5
    need(R-delta > Q(66,64), 'local convolution interval lies in one analytic nominal branch')
    state_values = ctx.actual_background(I(R))
    state = {k: J([v]+[0]*order) for k, v in state_values.items()}
    # ODE Taylor recurrence around the actual state at R. No derivative of
    # a state interval or of an unknown Taylor remainder is taken.
    for n in range(order):
        rho, eta, p, a, alpha = [state[k] for k in ('rho','eta','p','a','alpha')]
        V1 = ctx.init.poly(ctx.init.UP[1], eta)
        rhs = {'rho': -rho*rho*a, 'eta': -rho*p,
               'p': -rho*(V1-4*a*p), 'a': 1/rho+rho*p*p/3,
               'alpha': rho*(Q(5,2)*a-V1/p)*alpha}
        for key in state:
            state[key].c[n+1] = rhs[key].c[n]/(n+1)
    s = Q(1)/(Q(1,4)-(R-Q(3,2))**2)
    v = 3-2*R
    exp_packet = ctx.init.expbox(I(s-4))
    packet = J([ctx.init.poly(polys[n], I(s))*exp_packet*(v if n%2 else 1)
                /(ell**n*math.factorial(n)) for n in range(order+1)])
    exp_E = ctx.init.expbox(I(ell*R))
    E = J([(-1)**n*exp_E/math.factorial(n) for n in range(order+1)])
    actual_u, _ = ctx._registered(state, E, packet)
    nominal_u, _ = ctx.nominal_profiles(I(R), order=6)
    need(all(z.n >= 2 for z in actual_u) and len(actual_u) == 23, 'actual fifth-input jet supplies full23 second initial derivative')
    M3 = list(map(Q, initial['nominal_second_derivative_component_caps'][1]))
    need(len(M3) == 23 and min(M3) >= 0, 'one-sided nominal third derivative caps')
    rows, lower_square = [], Q(0)
    for i, label in enumerate(ctx.op.LABELS):
        actual = 2*ell**2*actual_u[i].c[2]
        nominal = 2*ell**2*nominal_u[i].c[2]
        remainder = M3[i]*35*delta/256
        smooth = nominal/2+I(-remainder, remainder)
        fixed_error = actual-smooth
        limit_error = actual-nominal/2
        minimum = Q(0) if fixed_error.lo <= 0 <= fixed_error.hi else min(abs(fixed_error.lo),abs(fixed_error.hi))
        lower_square += ctx.weights[i]*minimum**2
        rows.append({'label': label, 'input_Gram': str(ctx.weights[i]),
            'actual_u_rr': actual.row(), 'nominal_u_rr_left': nominal.row(),
            'nominal_u_rrr_piecewise_H_component_cap': str(M3[i]),
            'mollifier_endpoint_half_first_moment_error': str(remainder),
            'affine_smoothed_u_rr': smooth.row(),
            'actual_minus_affine_smoothed_u_rr': fixed_error.row(),
            'zero_width_limit_error': limit_error.row()})
    need(rows[12]['label'] == 'w:0', 'registered scalar mean')
    scalar = [Q(v) for v in rows[12]['actual_minus_affine_smoothed_u_rr']]
    limit = [Q(v) for v in rows[12]['zero_width_limit_error']]
    lower, upper = Q('-1.6246e33'), Q('-1.6245e33')
    need(lower < scalar[0] <= scalar[1] < upper < 0, 'strict scalar mean actual endpoint curvature error')
    need(lower < limit[0] <= limit[1] < upper < 0, 'strict nonzero limiting endpoint curvature error')
    need(lower_square > Q('1.6245e33')**2, 'strict complete positive H lower bound')
    ctx.reauthenticate()
    for key, path in paths.items():
        pin(path, PINS[key])
    print(json.dumps({'status': 'PASS', 'endpoint': ['0', '7/4'], 'delta_u': str(delta),
        'definition': 'actual U0_rr minus k_delta convolved nominal initial displacement_rr with affine tangent continuation beyondR',
        'actual_method': 'autonomous background ODE recurrence through order5 around actual value enclosure atR, actual bump polynomial jet through5, exact exponential enclosure; independent of paired Taylor-difference jets',
        'normalization': 'ordinary scaled radial second derivative partial_r^2 U0 = 2 ell^2 times physical-z normalized Taylor coefficient c2; not partial_y^2',
        'kernel_endpoint_identity': 'u_s_rr(R)=nominal_u_rr(R-)/2 + remainder, |remainder_i|<=M3_i*35*delta_u/256',
        'zero_width_limit': 'actual_u_rr(R) - nominal_u_rr(R-)/2; actual and nominal functions fixed whiledelta tends0+',
        'affine_endpoint_scalar_mean_strict_error_bounds': ['-1.6246e33','-1.6245e33'],
        'complete_H_endpoint_error_strict_lower': '1.6245e33',
        'complete_H_endpoint_error_squared_lower': str(lower_square),
        'actual_background_value_intervals': {k: val.row() for k,val in state_values.items()},
        'actual_background_normalized_order5_jets': {k: [z.row() for z in val.c] for k,val in state.items()},
        'rows': rows, 'algebra_checks': controls, 'input_sha256': PINS,
        'scope': 'conditional coordinate initial-curvature error lower bound for unchanged affine-continued candidate data, including fixeddelta and nonzero zero-width limit; no actual spacetime C2 upper certificate, differentiated bulk-source bound, Einstein residual, nonlinear family or Big Bang evidence'}, sort_keys=True))


if __name__ == '__main__':
    main()
