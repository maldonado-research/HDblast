"""Exact-rational error and leading-response enclosures; no waveform solver.

The analytic hypotheses and residual enclosures are mathematical inputs.
This program encloses arithmetic consequences, not the PDE theorem itself.
"""
from fractions import Fraction as F
from math import isqrt


def sqrt_bounds(value, bits=96):
    value = F(value)
    if value < 0 or not isinstance(bits, int) or bits < 8:
        raise ValueError('Nonnegative radicand and at least eight bits required')
    scale = 1 << bits
    integer = isqrt(value.numerator * scale * scale // value.denominator)
    lower = F(integer, scale)
    upper = lower if lower * lower == value else F(integer + 1, scale)
    return lower, upper


def exp_upper(value):
    """exp(x) <= (1-x/n)^(-n), with rational 0<=x/n<=1/2."""
    value = F(value)
    if not 0 <= value <= 1024:
        raise ValueError('Exponent outside the documented arithmetic budget')
    count = max(1, (2 * value.numerator + value.denominator - 1) // value.denominator)
    return (1 - value / count) ** (-count)


def arctan_bounds(inverse, terms=32):
    x = F(1, inverse)
    partial = sum(((-1) ** k * x ** (2 * k + 1) / (2 * k + 1)
                   for k in range(terms)), F(0))
    next_term = (-1) ** terms * x ** (2 * terms + 1) / (2 * terms + 1)
    return min(partial, partial + next_term), max(partial, partial + next_term)


def pi_bounds():
    # Machin identity: pi = 16 arctan(1/5) - 4 arctan(1/239).
    a, b = arctan_bounds(5), arctan_bounds(239)
    return 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]


def response_certificate(data):
    """Uniform-residual version of the shrinking-domain energy/trace bound.

    Exact rational scalar inputs are required; see the manuscript for their
    norms, sign conventions and independent analytic obligations.
    """
    if any(not isinstance(v, (str, int, F)) or isinstance(v, bool) for v in data.values()):
        raise ValueError('Use exact integers or rational strings, not floats')
    d = {key: F(value) for key, value in data.items()}
    required = {'T', 'L', 'lambda', 'nu', 'mu_lower', 'W_abs', 'beta_abs',
                'c_abs', 'F_l2_sup', 'G_abs_sup', 'X0_sq_upper',
                'incoming_sq_upper', 'delta_c', 'trial_q_sq_upper',
                'trial_D_sq_lower', 'trial_D_sq_upper',
                'rho_sq_lower', 'rho_sq_upper', 'p_sq_lower', 'p_sq_upper'}
    if set(d) != required:
        raise ValueError('Unexpected or missing certificate fields')
    if any(v < 0 for v in d.values()):
        raise ValueError('All scalar bounds must be nonnegative')
    if not (0 < d['T'] < d['L'] and d['lambda'] > 0 and d['nu'] > 0
            and d['mu_lower'] > 0 and 0 < d['rho_sq_lower'] <= d['rho_sq_upper']
            and 0 < d['p_sq_lower'] <= d['p_sq_upper']
            and d['trial_D_sq_lower'] <= d['trial_D_sq_upper']):
        raise ValueError('Positive inertia, denominator bounds and strict domain margin required')
    sl = sqrt_bounds(d['lambda'])[0]
    smn = sqrt_bounds(d['mu_lower'] * d['nu'])[0]
    if sl == 0 or smn == 0:
        raise ValueError('Increase square-root precision for very small inputs')
    a = max((d['lambda'] + d['W_abs']) / sl,
            (d['nu'] + d['beta_abs']) / smn)
    h_sq = d['F_l2_sup'] ** 2 + d['G_abs_sup'] ** 2 / d['mu_lower']
    h = sqrt_bounds(h_sq)[1]
    x0 = sqrt_bounds(d['X0_sq_upper'])[1]
    exponential = exp_upper(a * d['T'] / 2)
    B = exponential * (x0 + d['T'] * h)
    incoming = sqrt_bounds(d['incoming_sq_upper'])[1]
    trace_factor_sq = d['T'] * (1 / d['mu_lower'] + d['c_abs'] ** 2 / d['nu'])
    delta = (incoming + d['T'] * (d['F_l2_sup'] + d['W_abs'] * B / sl)
             + sqrt_bounds(trace_factor_sq)[1] * B
             + d['delta_c'] * sqrt_bounds(d['trial_q_sq_upper'])[1])
    norm_lower = max(F(0), sqrt_bounds(d['trial_D_sq_lower'])[0] - delta)
    norm_upper = sqrt_bounds(d['trial_D_sq_upper'])[1] + delta
    pil, piu = pi_bounds()
    prefactor_lower = 3 * pil ** 2 / (d['rho_sq_upper'] * d['p_sq_upper'])
    prefactor_upper = 3 * piu ** 2 / (d['rho_sq_lower'] * d['p_sq_lower'])
    return {'a_upper': a, 'forcing_squared_upper': h_sq,
            'exp_upper': exponential, 'energy_norm_upper': B,
            'trace_error_upper': delta, 'true_D_norm_lower': norm_lower,
            'true_D_norm_upper': norm_upper,
            'C_lower': prefactor_lower * norm_lower ** 2,
            'C_upper': prefactor_upper * norm_upper ** 2,
            'positive_lower_bound': norm_lower > 0}


def comparator_inputs(ell=F(1, 10**10)):
    """Analytic Dirichlet comparator and the separately proved collar box."""
    ell = F(ell)
    T, L, M = 3 * ell, 8 * ell, F(10**6)
    if not 0 < L <= F(1, 1250000000):
        raise ValueError('Comparator exceeds the validated collar depth')
    f = M * (1 + T) * sqrt_bounds(L)[1]
    return {'T': T, 'L': L, 'lambda': F(1), 'nu': F(1),
            'mu_lower': F(1, 1000), 'W_abs': M, 'beta_abs': F(1000),
            'c_abs': F(250), 'F_l2_sup': f, 'G_abs_sup': 16 / ell + 1,
            'X0_sq_upper': F(0), 'incoming_sq_upper': F(0),
            'delta_c': F(0), 'trial_q_sq_upper': F(0),
            'trial_D_sq_lower': 2 / ell,
            'trial_D_sq_upper': ell * (16 / ell + 1) ** 2,
            'rho_sq_lower': F(16641), 'rho_sq_upper': F(17161),
            'p_sq_lower': F(9, 625000000), 'p_sq_upper': F(49, 2500000000)}


def rational_json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: rational_json(item) for key, item in value.items()}
    return value
