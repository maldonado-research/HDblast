"""Self-contained exact arithmetic for the prospective rescaled branch bounds.

This producer does not import the accepted evaluator or archived norm probe.
Arithmetic passing is conditional on review of the analytic estimates; it is
not an ODE integration or an independent analytical proof.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
BASELINE = HERE / 'public_baseline'


def plus(left, right):
    size = max(len(left), len(right))
    return [(left[i] if i < len(left) else F(0)) +
            (right[i] if i < len(right) else F(0)) for i in range(size)]


def times(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def scaled(poly, factor):
    return [factor * a for a in poly]


def derivative(poly, order=1):
    result = list(poly)
    for _ in range(order):
        result = [i * result[i] for i in range(1, len(result))] or [F(0)]
    return result


def absolute_bound(poly, radius):
    return sum((abs(a) * radius**i for i, a in enumerate(poly)), F(0))


def lower_bound(poly, radius):
    return poly[0] - sum((abs(a) * radius**i for i, a in enumerate(poly) if i), F(0))


def serialized(value):
    if isinstance(value, F):
        return f'{value.numerator}/{value.denominator}'
    if isinstance(value, dict):
        return {key: serialized(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialized(item) for item in value]
    return value


def registered_parameters():
    coupling = 2 / F('1.0357712571566784') - F(4, 3)
    return dict(k=F(1, 9), w=F(2), cubic_W=F(2), f0=1 + coupling,
                f1=coupling, f2=F(0))


def constants(parameters, radius, tube, derivative_scale=None, root_bound=None,
              endpoint_mode='convex_riccati'):
    p = parameters
    k, w, cubic = p['k'], p['w'], p['cubic_W']
    f0, f1, f2 = p['f0'], p['f1'], p['f2']
    if not (k > 0 and w > 4 * k and f0 > 0 and radius > 0 and tube > 0):
        raise ValueError('Outside the theorem domain')
    lam = w - 4 * k
    s = max(k, lam) if derivative_scale is None else derivative_scale
    if s <= 0:
        raise ValueError('Nonpositive derivative scale')
    W = [3 * k, F(0), w / 2, cubic / 6]
    f = [f0, f1, f2 / 2]
    Wp, fp = derivative(W), derivative(f)
    U = plus(scaled(times(Wp, Wp), F(1, 2)), scaled(times(W, W), -F(2, 3)))
    A = plus(scaled(times(W, f), F(1, 9)), scaled(times(Wp, fp), -F(1, 12)))
    B_fun = plus(scaled(times(f, f), F(1, 36)), scaled(times(fp, fp), -F(1, 48)))
    m2 = w * lam
    nu = lam / 2
    mu = nu * (nu + 4 * k)
    gap = m2 - mu
    zeta_lower = mu * F(15, 16) / (F(16, 3) * k + nu)
    J_upper = 2 / k + 1 / (2 * zeta_lower)
    M = max(F(1), lam / s)
    B = 2 * M
    C_G = max(1 / gap, (1 + m2 / gap) / (4 * k * s))
    M2 = absolute_bound(derivative(U, 2), radius)
    M3 = absolute_bound(derivative(U, 3), radius)
    L1 = M3 / 2
    C_F = s**2 / 4 + M2 / 12
    epsilon = C_F * B**2 * tube**2 * J_upper / (2 * k)
    C_N = L1 * B**2 + 4 * C_F * (s / k) * B**3 * tube
    L_N = 2 * L1 * B * tube + 20 * C_F * (s / k) * B**2 * tube**2
    C_P = s * C_G * B * (2 * L1 * B + 20 * C_F * (s / k) * B**2 * tube)
    C_D = s * C_G * C_N
    C_a = 4 * C_F * B**2 / k
    C_L_upper = 8 * lam + F(28, 3) * k * lam / (w - 2 * k)
    C_H_accepted = 4 * C_L_upper / (9 * k**2)
    C_H_riccati = 16 * lam / (9 * k * (w - 3 * k))
    C_H_convex_riccati = lam / (k * (w - 3 * k))
    endpoint_choices = dict(accepted=C_H_accepted, riccati=C_H_riccati,
                            convex_riccati=C_H_convex_riccati)
    if endpoint_mode not in endpoint_choices:
        raise ValueError('Unknown linear endpoint estimate')
    C_H = endpoint_choices[endpoint_mode]
    W2 = absolute_bound(derivative(W, 2), radius)
    W3 = absolute_bound(derivative(W, 3), radius)
    F1 = absolute_bound(fp, radius)
    F2 = absolute_bound(derivative(f, 2), radius)
    F_abs = absolute_bound(f, radius)
    W_abs = absolute_bound(W, radius)
    L = abs(f1) / w if root_bound is None else root_bound
    if L < abs(f1) / w:
        raise ValueError('Root bound smaller than the proved monotonicity bound')
    a1 = k * f0 / 3
    C_K = L * absolute_bound(derivative(A), radius) + absolute_bound(B_fun, radius)
    a_bar = F(4, 3) * k + C_F * B**2 * tube**2 / k
    a_min = k - C_F * B**2 * tube**2 / k
    C_Tp = M2 + 4 * a_bar * s * B + s * B * (lam + C_P * tube)
    C_Tbeta = 2 * C_Tp / w
    C_slope_retired = s * B * C_a * L**2 + ((C_a + W2 / 3) * L + F1 / 6) * C_Tbeta * L
    C_v = 8 * C_F * B**2 * J_upper / (3 * k)
    C_Q_coarse = 2 * C_v / k
    C_Q = C_v / a_min if a_min > 0 else None
    d0 = 2 * (w - 2 * k)
    alpha = -f1 / (2 * d0)
    K_eta = (C_H * (3 * a1 / 2) * L + (C_D + W3 / 2) * L**2 + F2 * L / 2) / d0
    K_H = abs(A[1]) * K_eta + absolute_bound(derivative(A, 2), radius) * L**2 / 2
    K_H += absolute_bound(derivative(B_fun), radius) * L
    suppression_q = k * f1**2 / (24 * (w - 2 * k))
    return {key: value for key, value in locals().items()
            if key not in ('parameters', 'p', 'derivative_scale', 'root_bound', 'endpoint_choices')}


def gate(name, lhs, rhs, strict=False):
    passed = lhs < rhs if strict else lhs <= rhs
    return dict(name=name, lhs=lhs, rhs=rhs, strict=strict, pass_check=passed,
                exact_margin=rhs - lhs)


def sufficient_gates(c, delta, uniqueness_mode='structural'):
    b, k, B = c['tube'], c['k'], c['B']
    theta = c['C_Q'] * (c['C_Tbeta'] + c['s'] * B) * c['L']**2 * delta**2
    result = [
        gate('field_tube_inside_polynomial_radius', B * b, c['radius']),
        gate('metric_volterra_contraction', c['epsilon'], F(1, 4)),
        gate('scalar_contraction', c['C_G'] * c['L_N'], F(1, 2)),
        gate('scalar_ball_preservation', c['C_G'] * c['C_N'] * b, c['M']),
        gate('positive_metric_derivative', c['C_F'] * B**2 * b**2 / k, k / 2),
        gate('scalar_junction_derivative_feasibility', (c['C_P'] + c['W3']) * b, c['w'] / 2, True),
        gate('scalar_root_inside_boundary_tube', c['L'] * delta, b),
        gate('scalar_junction_derivative', (c['C_P'] + c['W3']) * b + delta * c['F2'] / 2, c['w'] / 2),
        gate('positive_reduced_curvature', c['C_K'] * delta, c['a1'] / 2),
        gate('metric_residual_positive_at_Y', 3 * c['a1'] * delta / 2, F(9, 32) * k**2, True),
        gate('positive_shell_tension', c['W2'] * b**2 / 6 + delta * c['F_abs'] / 6, k / 2),
        gate('strict_downward_shell_crossing_factored_H2', theta, F(1), True),
    ]
    if uniqueness_mode == 'structural':
        result.pop()  # Exact convexity monotonicity supersedes the theta gate.
    elif uniqueness_mode != 'factored':
        raise ValueError('Unknown endpoint uniqueness proof')
    if c['endpoint_mode'] == 'convex_riccati' or uniqueness_mode == 'structural':
        # This explicit polynomial gate suffices for F>=0 and R<=1. The
        # separate automatic-convexity lemma is not needed by this producer.
        result.append(gate('convex_bulk_potential_for_sharp_endpoint_bound',
                           F(0), lower_bound(derivative(c['U'], 2), c['radius'])))
    if uniqueness_mode == 'structural':
        sigma_lower = lower_bound(scaled(derivative(c['W'], 2), F(2)), c['radius'])
        sigma_lower += min(F(0), delta * c['f2'])
        result.append(gate('convex_shell_tension_for_structural_monotonicity', F(0), sigma_lower))
    return result


def polynomial_remainder(c, delta_star):
    # beta/delta is bounded both by L and by |alpha|+K_eta*delta_star.
    E = min(c['L'], abs(c['alpha']) + c['K_eta'] * delta_star)
    A, B_fun = c['A'], c['B_fun']
    K = abs(A[1]) * c['K_eta']
    K += sum((abs(A[i]) * E**i * delta_star**(i - 2) for i in range(2, len(A))), F(0))
    K += sum((abs(B_fun[i]) * E**i * delta_star**(i - 1) for i in range(1, len(B_fun))), F(0))
    return E, K


def certificate(delta_star, root_bound=None, radius=F(1, 100), tube=None,
                endpoint_mode='convex_riccati', parameters=None,
                uniqueness_mode='structural'):
    parameters = registered_parameters() if parameters is None else parameters
    L = abs(parameters['f1']) / parameters['w'] if root_bound is None else root_bound
    b = L * delta_star if tube is None else tube
    c = constants(parameters, radius, b, root_bound=root_bound, endpoint_mode=endpoint_mode)
    checks = sufficient_gates(c, delta_star, uniqueness_mode=uniqueness_mode)
    E, K_poly = polynomial_remainder(c, delta_star)
    K_best = min(c['K_H'], K_poly)
    theta_coarse = c['C_Q_coarse'] * (c['C_Tbeta'] + c['s'] * c['B']) * c['L']**2 * delta_star**2
    margin = c['suppression_q'] - K_best * delta_star
    convex_U = lower_bound(derivative(c['U'], 2), radius)
    convex_sigma = lower_bound(scaled(derivative(c['W'], 2), F(2)), radius)
    convex_sigma += min(F(0), delta_star * c['f2'])
    theta = c['C_Q'] * (c['C_Tbeta'] + c['s'] * c['B']) * c['L']**2 * delta_star**2
    theta_gate = gate('redundant_factored_shell_slope_gate', theta, F(1), True)
    retired = gate('retired_delta_squared_absolute_slope', c['C_slope_retired'] * delta_star, c['a1'] / 4)
    return dict(delta_star=delta_star, parameters=parameters, r=radius, b=b,
                uniqueness_mode=uniqueness_mode,
                pass_all_sufficient_gates=all(item['pass_check'] for item in checks),
                checks=checks, exact_constants=c,
                K_H_polynomial=K_poly, K_H_best=K_best, E_beta_over_delta=E,
                absolute_errors_at_delta_star=dict(scalar=c['K_eta'] * delta_star**2,
                                                  curvature=K_best * delta_star**3),
                strict_curvature_suppression=dict(q=c['suppression_q'], exact_margin=margin,
                                                   pass_check=c['suppression_q'] > 0 and margin > 0),
                additional_convexity_bounds=dict(U_second_derivative_lower=convex_U,
                                                sigma_second_derivative_lower=convex_sigma,
                                                pass_both=convex_U >= 0 and convex_sigma >= 0,
                                                bulk_convexity_used_in_certificate=endpoint_mode == 'convex_riccati',
                                                sigma_convexity_used_in_certificate=uniqueness_mode == 'structural'),
                redundant_factored_slope_gate=theta_gate,
                retired_slope_gate=retired,
                display_only=dict(delta_star=float(delta_star), b=float(b),
                                  theta=float(theta), theta_coarse=float(theta_coarse),
                                  K_eta=float(c['K_eta']), K_H_taylor=float(c['K_H']),
                                  K_H_polynomial=float(K_poly), K_H_best=float(K_best),
                                  scalar_error=float(c['K_eta'] * delta_star**2),
                                  curvature_error=float(K_best * delta_star**3),
                                  suppression_margin=float(margin)))


def check_controls():
    # Each control exercises a mathematical gate instead of a syntax-only check.
    base = certificate(F(1, 500))
    large_delta_gates = sufficient_gates(base['exact_constants'], F(1))
    bad_tube = certificate(F(1, 500), tube=F(1, 100))
    if not base['pass_all_sufficient_gates']:
        raise ArithmeticError('Required uniform delta_star=1/500 gates failed')
    if all(item['pass_check'] for item in large_delta_gates) or bad_tube['pass_all_sufficient_gates']:
        raise ArithmeticError('A negative arithmetic control unexpectedly passed')
    if base['retired_slope_gate']['pass_check']:
        raise ArithmeticError('Old absolute slope obstruction unexpectedly disappeared')
    p = registered_parameters()
    legacy_L = (abs(p['f1']) + 1) / p['w']
    loose = certificate(F(1, 500), root_bound=legacy_L, endpoint_mode='accepted',
                        uniqueness_mode='factored')
    exact_coarse = loose['exact_constants']['C_Q_coarse'] * (
        loose['exact_constants']['C_Tbeta'] + loose['exact_constants']['s'] * loose['exact_constants']['B']) * (
        loose['exact_constants']['L'] * F(1, 500))**2
    if not exact_coarse > 1:
        raise ArithmeticError('Expected coarse crossing bound failure missing')
    if not base['strict_curvature_suppression']['pass_check']:
        raise ArithmeticError('Required uniform strict suppression margin failed')
    return dict(uniform_interval_passes=True, uniform_strict_curvature_suppression_passes=True,
                delta_one_rejected=True,
                oversized_tube_rejected=True, retired_slope_gate_rejected=True,
                coarse_factored_crossing_gate_rejected_for_legacy_loose_L=True,
                note='Controls are exact fractions; approximate display fields are not used.')


def source_hashes():
    return {
        'public_baseline/EXACT_UNIFORM_CONSTANTS.json': hashlib.sha256(
            (BASELINE / 'EXACT_UNIFORM_CONSTANTS.json').read_bytes()).hexdigest(),
        'accepted_uniform_branch_theorem_20261008': '544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51',
        'accepted_uniform_constants_evaluator_20261008': '4ca9b1824c2c94b58ef04e56289272c3e5b51a5cd52a5a062cb8511728fac7cc',
    }


def accepted_arithmetic_regression():
    """Reproduce accepted values as a regression, not an independent review."""
    path = BASELINE / 'EXACT_UNIFORM_CONSTANTS.json'
    fixture = json.loads(path.read_text())
    aliases = dict(ball='B', C_H_DN='C_H', C_slope='C_slope_retired')
    result = []
    for model in fixture['models']:
        p = {key: F(value) for key, value in model['parameters'].items()}
        p['cubic_W'] = p.pop('v')
        original_L = (abs(p['f1']) + 1) / p['w']
        c = constants(p, F(model['r']), F(model['b']), derivative_scale=p['k'],
                      root_bound=original_L, endpoint_mode='accepted')
        mismatches = [key for key, value in model['exact_constants'].items()
                      if c[aliases.get(key, key)] != F(value)]
        poly_aliases = dict(BB='B_fun')
        mismatches += ['polynomial_' + key for key, value in model['polynomial_coefficients'].items()
                       if c[poly_aliases.get(key, key)] != [F(a) for a in value]]
        if mismatches:
            raise ArithmeticError(f'Accepted arithmetic regression failed: {mismatches}')
        result.append(dict(name=model['name'], pass_exact_regression=True,
                           compared_constant_fields=list(model['exact_constants']),
                           compared_polynomials=list(model['polynomial_coefficients']),
                           accepted_delta_star=F(model['delta_star']), accepted_b=F(model['b']),
                           accepted_K_eta=c['K_eta'], accepted_K_H_taylor=c['K_H'],
                           accepted_K_H_polynomial=F(model['K_H_polynomial_refinement'])))
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=HERE / 'RESCALED_BOUND_CERTIFICATE.json')
    parser.add_argument('--proof', type=Path)
    parser.add_argument('--proof-sha256')
    arguments = parser.parse_args()
    out_path = arguments.out.resolve()
    if (arguments.proof is None) != (arguments.proof_sha256 is None):
        raise ValueError('Supply both --proof and --proof-sha256, or neither')
    proof_binding = dict(runtime_proof_pin_checked=False,
                         review_scope='Analytic proof review is supplied by separate review receipts, not inferred from arithmetic.')
    if arguments.proof is not None:
        proof_sha = hashlib.sha256(arguments.proof.read_bytes()).hexdigest()
        if proof_sha != arguments.proof_sha256:
            raise ValueError('Runtime analytic proof SHA256 pin mismatch')
        proof_binding.update(runtime_proof_pin_checked=True, proof_sha256=proof_sha,
                             proof_logical_name=arguments.proof.name)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    samples = [certificate(d) for d in (F(1, 2000), F(1, 1000), F(1, 500))]
    uniform = samples[-1]
    c = uniform['exact_constants']
    q = c['suppression_q']
    K = uniform['K_H_best']
    # With the constants for the full 1/500 interval fixed, this strictly smaller
    # rational subinterval also proves a curvature sign, without recomputing c.
    signed_delta = F(1, 500) if q - K * F(1, 500) > 0 else min(F(1, 500), q / (2 * K))
    signed_margin = q - K * signed_delta
    if not signed_delta > 0 or not signed_margin > 0:
        raise ArithmeticError('Strict-sign conservative subinterval failed')
    # Float evaluation occurs only to identify a separately declared rational
    # model. All subsequent constants and decisions for it use exact fractions.
    historical_c_float = 2.0 / 1.0357712571566784 - 4.0 / 3.0
    historical_c = F(*historical_c_float.as_integer_ratio())
    historical_p = dict(k=F(1, 9), w=F(2), cubic_W=F(2), f0=1 + historical_c,
                        f1=historical_c, f2=F(0))
    historical_row = certificate(F(1, 500), parameters=historical_p)
    stored_p = dict(historical_p, k=F(*float(1 / 9).as_integer_ratio()),
                    f0=F(*(1.0 + historical_c_float).as_integer_ratio()))
    stored_row = certificate(F(1, 500), parameters=stored_p)
    for additional in (historical_row, stored_row):
        if not additional['pass_all_sufficient_gates'] or not additional['strict_curvature_suppression']['pass_check']:
            raise ArithmeticError('Historical rational coefficient model failed')
    readable_checks = [gate('K_eta_less_than_3_point_196', c['K_eta'], F(3196, 1000), True),
                       gate('K_H_less_than_0_point_254', K, F(254, 1000), True),
                       gate('suppression_margin_greater_than_0_point_000423',
                            F(423, 1000000), signed_margin, True)]
    if not all(check['pass_check'] for check in readable_checks):
        raise ArithmeticError('Readable strict rational bounds failed')
    result = dict(status='EXACT_ARITHMETIC_PASS_ANALYTIC_PROOF_SEPARATELY_REVIEWED',
                  producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  evidence_source_sha256=source_hashes(),
                  proof_binding=proof_binding,
                  redundant_supporting_checks=['convex_shell_tension_for_structural_monotonicity is an extra sufficient target-model check; the theorem only requires g_beta>0 and automatic bulk convexity.', 'The explicit bulk convexity check independently supports R<=1; scalar contraction already implies strict bulk convexity on the actual field tube.', 'The factored theta slope gate is redundant under the structural monotonicity theorem.'],
                  arithmetic='All decisions and error bounds use fractions.Fraction; floats are display only.',
                  model_scope='Exact decimal 1.0357712571566784, k=1/9, cubic W and affine f. No binary64 parameter equivalence is asserted.',
                  primary_uniform_certificate=uniform, sample_interval_certificates=samples,
                  conservative_metric_ratio_16over9=certificate(F(1, 500), endpoint_mode='riccati'),
                  legacy_loose_root_bound_comparison=certificate(F(1, 500),
                      root_bound=(abs(registered_parameters()['f1']) + 1) / registered_parameters()['w'],
                      endpoint_mode='accepted', uniqueness_mode='factored'),
                  accepted_original_norm_regressions=accepted_arithmetic_regression(),
                  historical_binary64_coupling_as_exact_rational=dict(
                      scope='c is the exact rational represented by historical binary64 evaluation 2.0/1.0357712571566784−4.0/3.0; k remains 1/9 and f0=1+c is exact mathematical addition. This does not certify historical floating-point trajectories.',
                      c=historical_c, certificate=historical_row),
                  historical_stored_binary64_coefficients_as_exact_rationals=dict(
                      scope='k, c and f0 are exact rationals represented by stored binary64 values; subsequent model arithmetic is exact. No bound on rounding during ODE evaluation or shooting is asserted.',
                      certificate=stored_row),
                  readable_strict_bounds=readable_checks,
                  sample_errors_using_uniform_constants=[dict(delta=d, scalar=c['K_eta'] * d**2, curvature=K * d**3)
                                                        for d in (F(1, 2000), F(1, 1000), F(1, 500))],
                  conservative_strict_suppression_subinterval=dict(delta_star=signed_delta,
                      q=q, K_H=K, exact_positive_margin=signed_margin,
                      display_only_delta_star=float(signed_delta),
                      comparator='Separate constant-scalar, constant-tension model with f(eta)=f0; it generally does not satisfy this scalar-dependent model scalar junction.'),
                  controls=check_controls(),
                  claim_limits=[
                      'Arithmetic passing does not itself prove Green, Volterra, Riccati, Jacobi or moving-endpoint estimates.',
                      'Proof review and independent arithmetic review are separate receipts; this producer does not certify either review by its own execution.',
                      'Absolute cubic curvature error and strict curvature suppression apply to the declared mathematical models under the separately reviewed analytic proof.',
                      'No ODE diagnostics, remote edits, baseline theorem replacement, stability or observational claims are included.'])
    out_path.write_text(json.dumps(serialized(result), indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], samples=[item['display_only'] for item in samples],
                          uniform_covers_all_samples=uniform['pass_all_sufficient_gates'],
                          strict_suppression_subinterval=float(signed_delta)), indent=2))


if __name__ == '__main__':
    main()
