"""Constructive interval tail budget for the frozen causal pulse experiment.

Importing this module performs no source, mode, root or response evaluation.
Call only after the encompassing registration is publicly frozen. This is a
directed-rounding interval implementation of analytic total-variation bounds,
not an assertion that every part of a numerical response run is certified.
"""

from fractions import Fraction
from functools import lru_cache
from math import inf, nextafter


INTERVAL_DECIMAL_PRECISION = 60
ROOT_BISECTIONS = 256
AMPLITUDE = Fraction(1, 10000)
OBSERVATIONS = tuple(Fraction(v) for v in ("-5.5", "-4.5", "-4", "-3.5", "-2.5", "-1.5"))
ROOT_POLYNOMIALS = {
    "positive_B": (2, -14, 21, -6),
    "signed_uB": (4, -32, 64, -36, 3),
}
ROOT_BRACKETS = {
    "positive_B": ((Fraction(3, 2), Fraction(13, 8)), (Fraction(5), Fraction(81, 16))),
    "signed_uB": ((Fraction(7, 4), Fraction(15, 8)), (Fraction(5), Fraction(6))),
}


def _polynomial(coefficients, value):
    answer = Fraction(0)
    for coefficient in coefficients:
        answer = answer * value + coefficient
    return answer


@lru_cache(maxsize=None)
def _root_bracket(source, index):
    left, right = ROOT_BRACKETS[source][index]
    coefficients = ROOT_POLYNOMIALS[source]
    left_value = _polynomial(coefficients, left)
    if left_value * _polynomial(coefficients, right) >= 0:
        raise RuntimeError("The registered critical-root bracket has no sign change")
    for _ in range(ROOT_BISECTIONS):
        midpoint = (left + right) / 2
        midpoint_value = _polynomial(coefficients, midpoint)
        if midpoint_value == 0:
            return midpoint, midpoint
        if left_value * midpoint_value < 0:
            right = midpoint
        else:
            left, left_value = midpoint, midpoint_value
    return left, right


def _rational(iv, value):
    return iv.mpf(value.numerator) / iv.mpf(value.denominator)


def _root_interval(iv, source, index):
    left, right = _root_bracket(source, index)
    # Rational division is outward rounded before the endpoint hull is made.
    return iv.mpf([_rational(iv, left).a, _rational(iv, right).b])


def _outward_pair(value):
    return [nextafter(float(value.a), -inf), nextafter(float(value.b), inf)]


def _primitive(iv, source, q):
    if source == "positive_B":
        polynomial = 4*q**4 - 12*q**3 + 6*q**2
        return iv.exp(1-q) * polynomial
    polynomial = 4*q**4 - 12*q**3 + 2*q**2
    return iv.sqrt(1-1/q) * iv.exp(1-q) * polynomial


def _budget_intervals(source, eta):
    if source not in ROOT_POLYNOMIALS:
        raise ValueError("Source must be positive_B or signed_uB")
    time = Fraction(str(eta))
    if time not in OBSERVATIONS:
        raise ValueError("Observation time is outside this frozen six-point registration")
    from mpmath.ctx_iv import MPIntervalContext
    iv = MPIntervalContext()
    iv.dps = INTERVAL_DECIMAL_PRECISION
    if time == Fraction(-11, 2):
        return iv, iv.mpf(0), iv.mpf(0), iv.mpf(0), {}
    q_minimum = _root_interval(iv, source, 0)
    q_maximum = _root_interval(iv, source, 1)
    minimum = _primitive(iv, source, q_minimum)
    maximum = _primitive(iv, source, q_maximum)
    separation = maximum - minimum
    exp_third = iv.exp(-iv.mpf(1)/3)
    if time in (Fraction(-9, 2), Fraction(-4)):
        normalized_A = 2*separation
    elif time == Fraction(-7, 2):
        normalized_A = 2*separation
        if source == "positive_B":
            normalized_A += -4 + iv.mpf(832)/81*exp_third
        else:
            normalized_A += iv.mpf(992)/81*exp_third
    else:
        normalized_A = 4*separation - (4 if source == "positive_B" else 0)
    if time == Fraction(-4):
        local_absolute = iv.mpf(1 if source == "positive_B" else 0)
    elif time in (Fraction(-9, 2), Fraction(-7, 2)):
        local_absolute = exp_third if source == "positive_B" else exp_third/2
    else:
        local_absolute = iv.mpf(0)
    amplitude = _rational(iv, AMPLITUDE)
    details = {
        "q_minimum_interval": _outward_pair(q_minimum),
        "q_maximum_interval": _outward_pair(q_maximum),
        "normalized_second_derivative_minimum_interval": _outward_pair(minimum),
        "normalized_second_derivative_maximum_interval": _outward_pair(maximum),
        "normalized_third_derivative_L1_interval": _outward_pair(4*separation-(4 if source == "positive_B" else 0)),
    }
    return iv, normalized_A*amplitude, normalized_A, local_absolute*amplitude, details


def derivative_tail_budget(source, eta):
    """Return an upper bound for |s''(eta)| + integral_{-6}^eta |s'''|.

    A_upper includes epsilon=1e-4. normalized_A_upper removes that amplitude.
    Every returned binary64 upper endpoint is advanced toward +infinity.
    """
    _, actual_A, normalized_A, source_absolute, details = _budget_intervals(source, eta)
    if Fraction(str(eta)) == Fraction(-11, 2):
        actual_upper = normalized_upper = source_upper = 0.0
    else:
        actual_upper = _outward_pair(actual_A)[1]
        normalized_upper = _outward_pair(normalized_A)[1]
        source_upper = _outward_pair(source_absolute)[1] if source_absolute != 0 else 0.0
    return {
        "source": source, "eta": str(eta), "A_upper": actual_upper,
        "normalized_A_upper": normalized_upper, "source_abs_upper": source_upper,
        "method": "analytic total variation; exact rational critical-root bisection; directed-rounding interval primitive evaluation",
        "interval_dps": INTERVAL_DECIMAL_PRECISION, "root_bisections": ROOT_BISECTIONS,
        "amplitude": "1/10000", "details": details,
    }


def combined_tail_bound(source, eta, cutoff):
    """Enclose the complete analytic tail expression, including pi and a.

    This covers the removed momentum interval only. It supplies no estimate
    for finite-cutoff quadrature, mode integration, or response cancellation.
    """
    cutoff = Fraction(str(cutoff))
    if cutoff <= 0:
        raise ValueError("The momentum cutoff must be positive")
    iv, actual_A, _, source_absolute, _ = _budget_intervals(source, eta)
    if Fraction(str(eta)) == Fraction(-11, 2):
        return 0.0
    a = -1/_rational(iv, Fraction(str(eta)))
    M_squared = 2*a*a
    bound = (3*source_absolute*M_squared/2 + actual_A/4)/(16*iv.pi**2*a*a*_rational(iv, cutoff)**2)
    return _outward_pair(bound)[1]


def normalized_combined_tail_bound(source, eta, cutoff):
    """Bound the removed tail of y=a^2*delta_Q/epsilon directly in intervals."""
    cutoff = Fraction(str(cutoff))
    if cutoff <= 0:
        raise ValueError("The momentum cutoff must be positive")
    iv, _, normalized_A, source_absolute, _ = _budget_intervals(source, eta)
    if Fraction(str(eta)) == Fraction(-11, 2):
        return 0.0
    a = -1/_rational(iv, Fraction(str(eta)))
    M_squared = 2*a*a
    normalized_source_absolute = source_absolute/_rational(iv, AMPLITUDE)
    bound = (3*normalized_source_absolute*M_squared/2 + normalized_A/4)/(16*iv.pi**2*_rational(iv, cutoff)**2)
    return _outward_pair(bound)[1]
