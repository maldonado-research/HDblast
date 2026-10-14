#!/usr/bin/env python3
"""Independent exact material junction algebra; no inherited science imports.

Checks are grouped identities, not independent physical results. Uses no
assert, and deliberately fails under normal and optimized execution when
--inject-failure is requested. Does not solve a hyperbolic interface PDE.
"""
import argparse
import json
from pathlib import Path
import sympy as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--inject-failure', action='store_true')
    args = parser.parse_args()
    checks, controls = [], []
    eps = sp.Symbol('eps')

    def zero(z):
        if isinstance(z, sp.MatrixBase):
            return all(zero(v) for v in z)
        if isinstance(z, (tuple, list)):
            return all(zero(v) for v in z)
        return sp.expand(sp.sympify(z)) == 0

    def check(name, residual):
        passed = zero(residual)
        checks.append({'name': name, 'passed': passed})
        if not passed:
            raise RuntimeError('positive check failed: ' + name)

    def control(name, residual):
        detected = not zero(residual)
        controls.append({'name': name, 'detected': detected})
        if not detected:
            raise RuntimeError('wrong formula escaped: ' + name)

    def symmetric(prefix):
        symbols = {(i, j): sp.Symbol(f'{prefix}{i}{j}') for i in range(4) for j in range(i, 4)}
        return sp.Matrix(4, 4, lambda i, j: symbols[tuple(sorted((i, j)))])

    def derivative0(q, order):
        if isinstance(q, sp.MatrixBase):
            return q.applyfunc(lambda z: sp.diff(z, eps, order).subs(eps, 0))
        return sp.diff(q, eps, order).subs(eps, 0)

    # This non-diagonal rational metric has exactly one negative eigenvalue:
    # its first 2x2 block has negative determinant; the last is positive.
    h0 = sp.Matrix([[-2, sp.Rational(1, 3), 0, 0],
                    [sp.Rational(1, 3), 3, 0, 0],
                    [0, 0, 2, sp.Rational(1, 5)],
                    [0, 0, sp.Rational(1, 5), 4]])
    inv0 = h0.inv()
    if not (h0[:2, :2].det() < 0 and h0[2:, 2:].det() > 0 and h0[2, 2] > 0):
        raise RuntimeError('test metric is not Lorentzian')
    h1, h2, K2 = symmetric('h1_'), symmetric('h2_'), symmetric('k2_')
    h = h0 + eps * h1 + eps**2 * h2 / 2
    s0, s1, s2, s3, psi1, psi2 = sp.symbols('sigma0 sigma1 sigma2 sigma3 psi1 psi2')
    dv = eps * psi1 + eps**2 * psi2 / 2
    tension = s0 + s1 * dv + s2 * dv**2 / 2 + s3 * dv**3 / 6
    tension_prime = s1 + s2 * dv + s3 * dv**2 / 2
    check('tension_second_derivative_chain_rule', derivative0(tension, 2) - s1 * psi2 - s2 * psi1**2)
    check('tension_prime_second_derivative_chain_rule', derivative0(tension_prime, 2) - s2 * psi2 - s3 * psi1**2)

    # Exact invertible trace reversal for every symmetric tensor S on h0.
    S = symmetric('sumK_')
    trS = sp.trace(inv0 * S)
    israel = S - trS * h0 + s0 * h0
    trace_israel = sp.trace(inv0 * israel)
    reversed_israel = israel - trace_israel * h0 / 3
    R = S - s0 * h0 / 3
    check('full_ten_component_Lorentz_trace_reversal', reversed_israel - R)
    check('invertible_trace_reversal_inverse', R - sp.trace(inv0 * R) * h0 - israel)
    control('wrong_intrinsic_trace_dimension', israel - trace_israel * h0 / 4 - R)
    control('double_common_shell_tension', S - trS * h0 + 2 * s0 * h0 - israel)

    # Sum of two outward curvatures, with correct first-order matching.
    Ksum = s0 * h0 / 3 + eps * (s0 * h1 / 3 + s1 * psi1 * h0 / 3) + eps**2 * K2 / 2
    residual2 = derivative0(Ksum - tension * h / 3, 2)
    linear2 = K2 - s0 * h2 / 3 - s1 * psi2 * h0 / 3
    metric_source = 2 * s1 * psi1 * h1 / 3 + s2 * psi1**2 * h0 / 3
    check('complete_second_order_ten_metric_junction_components', residual2 - (linear2 - metric_source))
    control('omit_metric_product_factor_two', residual2 - (linear2 - s1 * psi1 * h1 / 3 - s2 * psi1**2 * h0 / 3))
    p2 = sp.Symbol('sum_p2')
    psum = -s1 - eps * s2 * psi1 + eps**2 * p2 / 2
    scalar2 = derivative0(psum + tension_prime, 2)
    check('complete_second_order_scalar_junction', scalar2 - p2 - s2 * psi2 - s3 * psi1**2)
    control('wrong_scalar_third_derivative', scalar2 - p2 - s2 * psi2 - s2 * psi1**2)
    control('wrong_scalar_source_sign', scalar2 - p2 - s2 * psi2 + s3 * psi1**2)

    # Z2 reduction follows from exact copy equality, with one shell tension.
    check('exact_Z2_metric_division', residual2 / 2 - (K2 / 2 - s0 * h2 / 6 - s1 * psi2 * h0 / 6 - s1 * psi1 * h1 / 3 - s2 * psi1**2 * h0 / 6))
    check('exact_Z2_scalar_division', scalar2 / 2 - p2 / 2 - s2 * psi2 / 2 - s3 * psi1**2 / 2)
    control('single_copy_uses_two_copy_tension', residual2 / 2 - (K2 / 2 - s0 * h2 / 3 - s1 * psi2 * h0 / 3 - 2 * s1 * psi1 * h1 / 3 - s2 * psi1**2 * h0 / 3))

    # Actual registered tension; affine detuning cannot affect sigma'''.
    v, detuning, cstar = sp.symbols('varphi detuning cstar')
    z = v - 1
    registered = 2 * (sp.Rational(1, 3) + z**2 + z**3 / 3) + detuning * (1 + cstar + cstar * z)
    check('registered_tension_second_derivative', sp.diff(registered, v, 2) - 4 * v)
    check('registered_tension_third_derivative', sp.diff(registered, v, 3) - 4)
    check('registered_scalar_source_minus_four', -sp.diff(registered, v, 3) * psi1**2 + 4 * psi1**2)

    # Derive pullback jets by direct graph evaluation, including bending.
    ell, zeta, eta = sp.symbols('ell zeta eta')
    f0, f01, f02, f1, f11, f2 = sp.symbols('phi0 phi0_r phi0_rr phi1 phi1_r phi2')
    field = f0 + f01 * ell + f02 * ell**2 / 2 + eps * (f1 + f11 * ell) + eps**2 * f2 / 2
    pullback = field.subs(ell, eps * zeta + eps**2 * eta / 2)
    check('material_graph_first_scalar_derivative', derivative0(pullback, 1) - f1 - f01 * zeta)
    check('material_graph_second_scalar_derivative', derivative0(pullback, 2) - f2 - 2 * f11 * zeta - f01 * eta - f02 * zeta**2)
    control('discard_graph_bending_first_source', derivative0(pullback, 1).subs(f1, 0))
    control('omit_material_pullback_mixing_factor_two', derivative0(pullback, 2) - f2 - f11 * zeta - f01 * eta - f02 * zeta**2)

    # Independent inverse-metric expansion and varying-metric trace jets.
    inv = inv0 - eps * inv0 * h1 * inv0 + eps**2 * (inv0 * h1 * inv0 * h1 * inv0 - inv0 * h2 * inv0 / 2)
    check('inverse_metric_through_second_order', [derivative0(h * inv, j) - (sp.eye(4) if j == 0 else sp.zeros(4)) for j in range(3)])
    D0, D1, D2 = symmetric('deltaK0_'), symmetric('deltaK1_'), symmetric('deltaK2_')
    D = D0 + eps * D1 + eps**2 * D2 / 2
    trD2 = derivative0(sp.trace(inv * D), 2)
    trace_formula = sp.trace(inv0 * D2) - 2 * sp.trace(inv0 * h1 * inv0 * D1) + sp.trace((2 * inv0 * h1 * inv0 * h1 * inv0 - inv0 * h2 * inv0) * D0)
    check('normal_force_trace_includes_all_metric_variations', trD2 - trace_formula)
    control('omit_inverse_metric_cross_term', trD2 - (trace_formula + 2 * sp.trace(inv0 * h1 * inv0 * D1)))

    # Exact off-shell Gauss difference, using arbitrary two-sided tensors.
    plus, minus = symmetric('Kplus_'), symmetric('Kminus_')
    pp, pm, up = sp.symbols('pplus pminus sigma_prime')
    delta = plus - minus
    trplus, trminus = sp.trace(inv0 * plus), sp.trace(inv0 * minus)
    normplus = sp.trace(inv0 * plus * inv0 * plus)
    normminus = sp.trace(inv0 * minus * inv0 * minus)
    twice_bulk_difference = trplus**2 - trminus**2 - normplus + normminus - pp**2 + pm**2
    I = plus + minus - (trplus + trminus) * h0 + s0 * h0
    P = pp + pm + up
    force = s0 * sp.trace(inv0 * delta) + up * (pp - pm)
    check('exact_off_shell_normal_force_identity', force - twice_bulk_difference - sp.trace(inv0 * I * inv0 * delta) - P * (pp - pm))
    control('wrong_outward_normal_sum_sign', force - twice_bulk_difference - sp.trace(inv0 * I * inv0 * delta) + P * (pp - pm))

    # The Codazzi identity reduces to the shell-tension chain rule.
    ya = sp.symbols('y0:4')
    intrinsic_phi = sum((i + 1) * ya[i]**2 for i in range(4))
    sigmaphi = s0 + s1 * intrinsic_phi + s2 * intrinsic_phi**2 / 2
    check('tangential_Codazzi_tension_product_rule', [sp.diff(sigmaphi, q) - (s1 + s2 * intrinsic_phi) * sp.diff(intrinsic_phi, q) for q in ya])
    # Pure conformal material h1 gives a pure-trace metric source on h0.
    conformal = sp.Symbol('conformal')
    source_conf = 2 * s1 * psi1 * conformal * h0 / 3 + s2 * psi1**2 * h0 / 3
    check('conformal_first_metric_source_has_zero_tracefree_part', source_conf - sp.trace(inv0 * source_conf) * h0 / 4)

    # A nonzero ambient-time shell point: X=(3,5,0,0,0), s=16.
    # A normalized tangent frame independently proves conformal pullback.
    eta5 = sp.diag(-1, 1, 1, 1, 1)
    X = sp.Matrix([3, 5, 0, 0, 0])
    tangents = sp.Matrix([[sp.Rational(5, 4), 0, 0, 0],
                          [sp.Rational(3, 4), 0, 0, 0],
                          [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    F, Omega = sp.symbols('masterF Omega', nonzero=True)
    flat_induced = tangents.T * eta5 * tangents
    covX = eta5 * X
    dg = 2 * F * (16 * eta5 - 3 * covX * covX.T)
    check('even_master_shell_metric_pullback_is_conformal',
          [X.T * eta5 * tangents, tangents.T * dg * tangents - 32 * F * flat_induced])

    r = sp.Symbol('r')
    rho, t2 = sp.Function('rho')(r), sp.Function('tensor2')(r)
    warped_h2 = rho**2 * t2
    flux = sp.diff(warped_h2, r) / 2 - sp.diff(rho, r) * warped_h2 / rho
    check('tensor_Robin_flux_becomes_rescaled_Neumann', flux - rho**2 * sp.diff(t2, r) / 2)
    control('bare_metric_Neumann_omits_warp_rescaling', flux - sp.diff(warped_h2, r) / 2)

    if args.inject_failure:
        check('deliberate_failure', 1)
    report = {'status': 'PASS', 'positive_checks': len(checks), 'wrong_formula_controls': len(controls),
              'scope': 'Exact material junction algebra and Gauss/Codazzi consistency; no shell-contact PDE solve',
              'checks': checks, 'controls': controls, 'registered_sigma_third_derivative': 4,
              'registered_two_copy_scalar_source': '-4 psi1^2',
              'sympy_version': sp.__version__}
    payload = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(payload)
    print(payload, end='')


if __name__ == '__main__':
    main()
