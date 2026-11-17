"""Exact finite inequality ledgers supporting the effective observer proofs.

Only bounded rational arithmetic and affine exponents a*P+b are evaluated,
where P=10^500207 is a SYMBOL with P>=10^6. Neither P nor 10^P is built.
These checks certify the displayed algebra and inequality directions, not
the continuum metric, observer, support, curvature or Jacobi estimates.
No ancestor source producer is imported or executed. This module has no CLI.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import factorial, prod
import re


P_INDEX = 500207
ELL = F(1, 1000000)
INTEGER_DEFAULTS = {
    "amplitude_power_index": 500209,
    "growth_power_index": 500212,
    "M_prefactor_power": 7007411,
    "embedding_power": 7,
    "metric_gain_power": 11,
    "Taylor_order": 11,
}
INTEGER_BOUNDS = {
    "amplitude_power_index": (500200, 500220),
    "growth_power_index": (500200, 500220),
    "M_prefactor_power": (0, 10000000),
    "embedding_power": (0, 20),
    "metric_gain_power": (0, 30),
    "Taylor_order": (1, 14),
}
RATIONAL_DEFAULTS = {
    "P_floor": "1000000",
    "growth_rounding_coefficient": "1",
    "M_decay_coefficient": "10",
    "coframe_ceiling": "1936",
    "delta_cap": "1/100000000000000000000",
    "inverse_operator_ceiling": "2",
    "inverse_difference_factor": "2",
    "inverse_derivative_ceiling": "5",
    "inverse_difference_derivative_factor": "8",
    "connection_difference_factor": "6",
    "connection_derivative_factor": "48",
    "ambient_connection_ceiling": "2",
    "ambient_connection_derivative_ceiling": "15",
    "Q_factor": "25/16",
    "T_factor": "13/8",
    "qstar_factor": "6145/4096",
    "stopped_velocity_squared_ceiling": "2",
    "Dini_forcing_factor": "24",
    "velocity_error_factor": "48",
    "drift_error_factor": "24",
    "past_time_slope": "7/8",
    "phase_speed_ceiling": "16",
    "phase_Lipschitz_ceiling": "256",
    "mixed_curvature_error_factor": "144",
    "covariant_curvature_error_factor": "304",
    "actual_metric_tidal_error_factor": "512",
    "frame_growth_rate": "4",
    "frame_norm_ceiling": "2",
    "Jacobi_E_error_factor": "32768",
    "Jacobi_operator_bootstrap": "2",
    "Jacobi_S_error_factor": "65536",
    "singular_value_floor": "3/4",
    "cut_metric_floor": "1/2",
    "relative_density_error_factor": "4",
    "relative_area_error_factor": "262144",
    "Taylor_kernel_denominator": "39916800",
}
MODE_DEFAULTS = {
    "source_norm_mode": "full_arrays_distinct_operator_norms",
    "causal_shadow_mode": "shadow_zero_and_upper_bound_only",
    "duration_scope": "certified_Q_only",
    "homogeneous_zero_mode": "weak_at_zero",
    "Jacobi_sign_mode": "hyperbolic_plus",
}
MODES = {
    "source_norm_mode": {"full_arrays_distinct_operator_norms", "identify_Frobenius_with_operator"},
    "causal_shadow_mode": {"shadow_zero_and_upper_bound_only", "nonzero_arrival_from_upper_bound"},
    "duration_scope": {"certified_Q_only", "import_unverified_long_duration"},
    "homogeneous_zero_mode": {"weak_at_zero", "strict_at_zero"},
    "Jacobi_sign_mode": {"hyperbolic_plus", "oscillatory_minus"},
}
DEFAULTS = {**INTEGER_DEFAULTS, **RATIONAL_DEFAULTS, **MODE_DEFAULTS}
FACTORS = tuple(DEFAULTS)
_WRONG = {
    "amplitude_power_index": 500206, "growth_power_index": 500215,
    "M_prefactor_power": 7007410, "embedding_power": 6,
    "metric_gain_power": 10, "Taylor_order": 12,
    "P_floor": "1", "growth_rounding_coefficient": "1/4",
    "M_decay_coefficient": "100", "coframe_ceiling": "1935", "delta_cap": "1/2",
    "inverse_operator_ceiling": "1", "inverse_difference_factor": "1",
    "inverse_derivative_ceiling": "4", "inverse_difference_derivative_factor": "7",
    "connection_difference_factor": "3", "connection_derivative_factor": "24",
    "ambient_connection_ceiling": "1/100000", "ambient_connection_derivative_ceiling": "13",
    "Q_factor": "13/8", "T_factor": "25/16", "qstar_factor": "6145/2048",
    "stopped_velocity_squared_ceiling": "3/2", "Dini_forcing_factor": "12",
    "velocity_error_factor": "24", "drift_error_factor": "12", "past_time_slope": "1",
    "phase_speed_ceiling": "14", "phase_Lipschitz_ceiling": "220",
    "mixed_curvature_error_factor": "120", "covariant_curvature_error_factor": "200",
    "actual_metric_tidal_error_factor": "256", "frame_growth_rate": "3",
    "frame_norm_ceiling": "1", "Jacobi_E_error_factor": "16384",
    "Jacobi_operator_bootstrap": "1", "Jacobi_S_error_factor": "32768",
    "singular_value_floor": "1", "cut_metric_floor": "1",
    "relative_density_error_factor": "3", "relative_area_error_factor": "131072",
    "Taylor_kernel_denominator": "479001600",
}
MUTATIONS = tuple([(name, _WRONG[name]) for name in (*INTEGER_DEFAULTS, *RATIONAL_DEFAULTS)]
                  + [(name, next(value for value in sorted(MODES[name]) if value != default))
                     for name, default in MODE_DEFAULTS.items()])
_RATIONAL_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def rational_text(value):
    value = F(value)
    return str(value.numerator) if value.denominator == 1 else str(value.numerator) + "/" + str(value.denominator)


def rational(value):
    if type(value) is not str or len(value) > 80 or not _RATIONAL_RE.fullmatch(value):
        raise ValueError("bounded canonical rational string required")
    p, separator, q = value.partition("/")
    numerator, denominator = int(p), int(q) if separator else 1
    if abs(numerator) > 10 ** 30 or denominator > 10 ** 30:
        raise ValueError("rational exceeds finite parser bound")
    result = F(numerator, denominator)
    if rational_text(result) != value:
        raise ValueError("noncanonical rational")
    return result


def ten_small(exponent):
    if type(exponent) is not int or not -40 <= exponent <= 40:
        raise ValueError("only small decimal powers may be materialized")
    return F(10 ** exponent) if exponent >= 0 else F(1, 10 ** (-exponent))


def decimal_lower_power(value):
    """Return k with 10^k<=value, using at most forty tiny powers."""
    if value <= 0:
        raise ValueError("positive comparison threshold required")
    for magnitude in range(41):
        if ten_small(-magnitude) <= value:
            return -magnitude
    raise ValueError("threshold exceeds the bounded rational resolution")


@dataclass(frozen=True)
class Affine:
    """An exponent a*P+b; P is never evaluated."""
    a: F = F(0)
    b: F = F(0)

    def __add__(self, other):
        return Affine(self.a + other.a, self.b + other.b)

    def __sub__(self, other):
        return Affine(self.a - other.a, self.b - other.b)

    def at(self, floor):
        return self.a * floor + self.b

    def report(self):
        return {"P_coefficient": rational_text(self.a), "constant": rational_text(self.b)}


def uniform_less(left, right, floor):
    gap = right - left
    return gap.a >= 0 and gap.at(floor) > 0


def exact(value):
    if isinstance(value, F):
        return rational_text(value)
    if isinstance(value, Affine):
        return value.report()
    if isinstance(value, dict):
        return {key: exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [exact(item) for item in value]
    return value


def polynomial_add(left, right):
    result = dict(left)
    for powers, coefficient in right.items():
        result[powers] = result.get(powers, F(0)) + coefficient
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def polynomial_scale(polynomial, factor):
    return {powers: coefficient * factor for powers, coefficient in polynomial.items() if coefficient * factor}


def polynomial_multiply(left, right):
    result = {}
    for powers_a, coefficient_a in left.items():
        for powers_b, coefficient_b in right.items():
            powers = tuple(a + b for a, b in zip(powers_a, powers_b))
            result[powers] = result.get(powers, F(0)) + coefficient_a * coefficient_b
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def evaluate(claims, check):
    """Name-first check(name,passed) callback; no source producer execution."""
    records = []

    def need(name, condition, witness=None):
        passed = bool(condition)
        full = "effective_oracles/" + name
        check(full, passed)
        row = {"name": name, "status": "pass" if passed else "fail"}
        if witness is not None:
            row["witness"] = exact(witness)
        records.append(row)
        if not passed:
            raise ValueError("effective oracle rejects " + full)

    def upper(name, used, ceiling, strict=False):
        need(name, used < ceiling if strict else used <= ceiling,
             {"used": used, "relation": "<" if strict else "<=", "ceiling": ceiling})

    def exponent_upper(name, left, right, floor):
        gap = right - left
        need(name, uniform_less(left, right, floor),
             {"left_exponent": left, "right_exponent": right,
              "right_minus_left": gap, "P_lower": floor,
              "gap_at_lower": gap.at(floor), "gap_slope_nonnegative": gap.a >= 0})

    need("schema", type(claims) is dict and claims.get("schema") == "hdblast-effective-oracle-claims-v1" and set(claims) == {"schema", *FACTORS})
    need("exact_engine", F(1, 3) * 3 == 1 and F(2, 4) == F(1, 2))
    for name, (low, high) in INTEGER_BOUNDS.items():
        need("type/" + name, type(claims[name]) is int and low <= claims[name] <= high)
    values = {}
    for name in RATIONAL_DEFAULTS:
        try:
            values[name] = rational(claims[name])
        except (ValueError, TypeError):
            need("parser/" + name, False)
    for name in MODE_DEFAULTS:
        need("type/" + name, type(claims[name]) is str and claims[name] in MODES[name])
    f = values.__getitem__
    need("interval/positive_declared_constants", all(value > 0 for value in values.values()))

    # Full-array and operator norms are different even on the identity.
    reported_identity_F_squared = F(3) if claims["source_norm_mode"] == "full_arrays_distinct_operator_norms" else F(1)
    need("norms/full_array_operator_distinction", reported_identity_F_squared == 3,
         {"identity3_Frobenius_squared": 3, "identity3_operator_squared": 1,
          "claimed_Frobenius_squared": reported_identity_F_squared})
    p_floor, delta = f("P_floor"), f("delta_cap")
    need("exponents/P_positive_lower_witness", 0 < p_floor <= 10 ** 30 and P_INDEX > 30)
    need("interval/nonnegative_perturbation_cap", 0 < delta <= 1)
    amplitude_coefficient = ten_small(claims["amplitude_power_index"] - P_INDEX)
    growth_coefficient = 2 * f("T_factor") * ten_small(claims["growth_power_index"] - 6 - P_INDEX)
    # n!>=2^(n-1) for n>=2 gives e<=2+(1/2)/(1-1/2)=3<10.
    Euler_upper = 2 + F(1, 2) / (1 - F(1, 2))
    upper("series/Euler_decimal_base_rounding", Euler_upper, F(10), strict=True)
    upper("exponents/growth_rounding", growth_coefficient, f("growth_rounding_coefficient"), strict=True)
    # 2*10^(7003210+4200)<10^7007411; compare decimal indices, not huge integers.
    need("exponents/M_prefactor_rounding", claims["M_prefactor_power"] - (7003210 + 4200) >= 1,
         {"leading_factor": 2, "source_decimal_power": 7007410,
          "rounded_decimal_power": claims["M_prefactor_power"]})
    need("exponents/positive_decay_coefficient", f("M_decay_coefficient") > 0)
    M_log_upper = Affine(-amplitude_coefficient + f("growth_rounding_coefficient"), F(claims["M_prefactor_power"]))
    M_log_target = Affine(-f("M_decay_coefficient"), F(0))
    exponent_upper("exponents/M_uniform_decay", M_log_upper, M_log_target, p_floor)

    coframe_used = F(16 * sum(3 ** order for order in range(5)))
    upper("metric/C4_coframe_count", coframe_used, f("coframe_ceiling"))
    upper("source/boundary_inclusive_full_tensor_H3", F(10 ** 7), ten_small(claims["embedding_power"]))
    metric_conversion = f("coframe_ceiling") * ten_small(claims["embedding_power"])
    upper("metric/C4_decimal_rounding", metric_conversion, ten_small(claims["metric_gain_power"]), strict=True)
    delta_log = Affine(-f("M_decay_coefficient"), F(claims["metric_gain_power"]))
    cap_decimal = decimal_lower_power(delta)
    exponent_upper("exponents/metric_delta_below_cap", delta_log, Affine(F(0), F(cap_decimal)), p_floor)
    scalar_log = Affine(-f("M_decay_coefficient"), F(claims["embedding_power"]))
    exponent_upper("exponents/scalar_value_upper", scalar_log, Affine(F(0), F(-20)), p_floor)
    need("scalar/C4_array_count", sum(2 ** order for order in range(5)) == 31)

    Q, T, qstar = f("Q_factor") * ELL, f("T_factor") * ELL, f("qstar_factor") * ELL
    need("domain/positive_horizons", Q > 0 and T > 0 and qstar > 0)
    upper("domain/Q_below_two_ell", Q, 2 * ELL, strict=True)
    upper("domain/Q_below_half", Q, F(1, 2), strict=True)
    upper("domain/qstar_covered", qstar, Q, strict=True)
    need("domain/old_slab_buffer", T - Q >= ELL / 16,
         {"T_minus_Q": T - Q, "required_minimum": ELL / 16})
    need("domain/qstar_gap_identity", F(25, 16) - F(6145, 4096) == F(255, 4096))

    # All homogeneous perturbation statements are weak at delta=0.
    upper("metric/Neumann_quarter_cap", delta, F(1, 4), strict=True)
    inverse_ratio = 1 / (1 - delta)
    upper("metric/inverse_operator", inverse_ratio, f("inverse_operator_ceiling"), strict=True)
    upper("metric/inverse_difference_coefficient", inverse_ratio, f("inverse_difference_factor"))
    G, K = f("inverse_operator_ceiling"), f("inverse_difference_factor")
    upper("metric/inverse_derivative", G ** 2 * (1 + delta), f("inverse_derivative_ceiling"), strict=True)
    DK_used = K * F(1, 2) * G + G ** 2 + K
    upper("metric/inverse_difference_derivative", DK_used, f("inverse_difference_derivative_factor"))
    Gamma_used = F(3, 2) * (G + K)
    upper("metric/connection_difference", Gamma_used, f("connection_difference_factor"))
    DGamma_used = F(3, 2) * (f("inverse_difference_derivative_factor") + f("inverse_derivative_ceiling") + 8 * K + G)
    upper("metric/connection_derivative", DGamma_used, f("connection_derivative_factor"))
    upper("metric/background_first_derivative", 8 * T, F(1, 2), strict=True)
    upper("metric/background_connection", 12 * T, F(1), strict=True)
    upper("metric/actual_metric_operator", F(3, 2) + delta, F(2), strict=True)
    upper("metric/actual_metric_full_Frobenius", 2 * (F(3, 2) + delta), F(4), strict=True)
    upper("metric/actual_connection", 12 * T + f("connection_difference_factor") * delta, f("ambient_connection_ceiling"), strict=True)
    upper("metric/actual_connection_derivative", F(27, 2) + f("connection_derivative_factor") * delta, f("ambient_connection_derivative_ceiling"), strict=True)

    stop_squared = F(5, 4) ** 2 + F(1, 4) ** 2
    upper("phase/stopped_full_velocity_squared", stop_squared, f("stopped_velocity_squared_ceiling"), strict=True)
    upper("phase/background_Dini_linear_rounding", 3 * T, F(1), strict=True)
    Dini_used = 2 * f("connection_difference_factor") * f("stopped_velocity_squared_ceiling")
    upper("phase/two_acceleration_output_blocks", Dini_used, f("Dini_forcing_factor"))
    # exp(q)<=1/(1-q); integrate exp(s) on [0,q] using the common Q bound.
    velocity_used = f("Dini_forcing_factor") / (1 - Q)
    upper("phase/Dini_velocity_comparison", velocity_used, f("velocity_error_factor"))
    upper("phase/integrated_drift", f("velocity_error_factor") / 2, f("drift_error_factor"))
    upper("phase/strict_velocity_stop_improvement", f("velocity_error_factor") * delta * Q, F(1, 8), strict=True)
    upper("phase/velocity_integer_witness", 16 * f("velocity_error_factor"), 1 / (delta * ELL), strict=True)
    upper("phase/drift_integer_witness", 128 * f("drift_error_factor"), 1 / (delta * ELL), strict=True)
    upper("phase/drift_before_inner_slab", f("drift_error_factor") * delta * Q ** 2, ELL / 32, strict=True)
    need("phase/past_time_slope", 1 - f("velocity_error_factor") * delta * Q > f("past_time_slope"),
         {"derived_slope_floor": 1 - f("velocity_error_factor") * delta * Q, "claimed_floor": f("past_time_slope")})
    need("phase/closed_slab_retention", T - Q - f("drift_error_factor") * delta * Q ** 2 > ELL / 32)

    A, DA = f("ambient_connection_ceiling"), f("ambient_connection_derivative_ceiling")
    speed_used = F(3, 2) + F(1, 2) + 2 * A * 3 + F(1, 2) ** 2 * F(3, 2)
    upper("continuation/global_phase_speed", speed_used, f("phase_speed_ceiling"), strict=True)
    Ht, Homega, Ha = 3 * DA, 6 * DA + 4 * A + 6 * A, 4 * A
    derivative_columns = (2 * Ht, 2 * Homega + F(1, 4), 2 * Ha + 1, 2 * Ha + F(3, 2) + 1)
    upper("continuation/all_phase_input_output_slots", max(derivative_columns), f("phase_Lipschitz_ceiling"), strict=True)
    # For omega'=W, W'=-|W|^2 omega-A^S with omega.A^S=0:
    # c1'=2c2 and c2'=s-s*r=-s*(r-1), as exact polynomial coefficients.
    one, r, s = {(0, 0): F(1)}, {(1, 0): F(1)}, {(0, 1): F(1)}
    direct_constraint_second = polynomial_add(s, polynomial_scale(polynomial_multiply(s, r), -1))
    expected_constraint_second = polynomial_scale(polynomial_multiply(s, polynomial_add(r, polynomial_scale(one, -1))), -1)
    need("sphere/past_phase_constraint_rows", direct_constraint_second == expected_constraint_second,
         {"c1prime": "2*c2", "c2prime_direct": "s-s*r", "c2prime_constraint": "-s*(r-1)",
          "expanded_coefficient_residual_is_zero": polynomial_add(direct_constraint_second, polynomial_scale(expected_constraint_second, -1)) == {}})

    Gamma, DGamma = f("connection_difference_factor"), f("connection_derivative_factor")
    mixed_used = 2 * DGamma + 4 * Gamma + 2 * Gamma ** 2 * F(1, 4)
    upper("curvature/mixed_full_array", mixed_used, f("mixed_curvature_error_factor"))
    background_R_F = F(2 * 2 * 4)
    upper("curvature/lower_with_actual_metric", 2 * f("mixed_curvature_error_factor") + background_R_F, f("covariant_curvature_error_factor"))
    W_difference = F(2 * (4 + 4))
    upper("curvature/actual_metric_W_telescoping", W_difference, F(32))
    upper("curvature/compare_to_actual_metric_W", f("covariant_curvature_error_factor") + 32, f("actual_metric_tidal_error_factor"))

    upper("frame/connection_velocity_generator", 2 * A, f("frame_growth_rate"))
    frame_exponent = f("frame_growth_rate") * Q
    upper("frame/exponent_below_half", frame_exponent, F(1, 2), strict=True)
    frame_series_upper = 1 / (1 - frame_exponent)
    upper("frame/exp_geometric_majorant", frame_series_upper, f("frame_norm_ceiling"), strict=True)
    upper("frame/all_four_vector_caps", f("stopped_velocity_squared_ceiling"), f("frame_norm_ceiling") ** 2, strict=True)
    E_used = 3 * f("frame_norm_ceiling") ** 4 * f("actual_metric_tidal_error_factor")
    upper("Jacobi/full_three_by_three_error", E_used, f("Jacobi_E_error_factor"))

    # W(g)(E_i,V,V,E_j)=delta_ij in a unit timelike orthonormal screen.
    # The background cosh(q)I has S''(0)=I and cannot solve S''=-S.
    sign = 1 if claims["Jacobi_sign_mode"] == "hyperbolic_plus" else -1
    basis4 = tuple(tuple(F(int(i == j)) for i in range(4)) for j in range(4))
    V4, screens = tuple(-x for x in basis4[0]), basis4[1:]
    dot4 = lambda left, right: -left[0] * right[0] + sum(left[i] * right[i] for i in range(1, 4))
    background_tidal = tuple(tuple(dot4(left, V4) * dot4(V4, right) - dot4(left, right) * dot4(V4, V4) for right in screens) for left in screens)
    cosh_second_jet = factorial(2) * F(1, factorial(2))
    residual_tidal = tuple(tuple(cosh_second_jet * int(i == j) - sign * background_tidal[i][j] for j in range(3)) for i in range(3))
    need("Jacobi/background_curvature_and_cosh_sign", all(value == 0 for row in residual_tidal for value in row),
         {"S0": "Id", "Sprime0": "0", "Ssecond0": "Id", "equation_sign": sign,
          "universal_unit_tidal_matrix": background_tidal, "initial_equation_residual": residual_tidal})
    # Positive factorial-series coefficients yield exp(x)<=1/(1-x).
    # For cosh, coefficients after the constant are <=1/2; at x<=1/2,
    # cosh(x)-1 <=x^2/[2(1-x^2)] <=x^2.
    upper("series/cosh_quadratic_majorant", 1 / (2 * (1 - F(1, 4))), F(1))
    upper("series/sinh_linear_kernel", F(5, 4), F(2))
    S_error_used = 2 * f("Jacobi_E_error_factor") * f("Jacobi_operator_bootstrap") / 2
    upper("Jacobi/Volterra_kernel_and_operator_bootstrap", S_error_used, f("Jacobi_S_error_factor"))
    d_uniform = f("Jacobi_S_error_factor") * delta / 4
    upper("Jacobi/error_below_quarter", d_uniform, F(1, 4), strict=True)
    upper("Jacobi/bootstrap_strict_improvement", F(5, 4) + d_uniform, F(3, 2), strict=True)
    upper("Jacobi/closed_operator_bootstrap", F(3, 2), f("Jacobi_operator_bootstrap"), strict=True)
    need("cut/minimum_singular_value", 1 - d_uniform > f("singular_value_floor"),
         {"derived_floor": 1 - d_uniform, "claimed_floor": f("singular_value_floor")})
    need("cut/full_positive_metric_gap", f("singular_value_floor") ** 2 > f("cut_metric_floor"),
         {"squared_singular_floor": f("singular_value_floor") ** 2, "cut_floor": f("cut_metric_floor")})

    upper("area/upper_cubic_density_deviation", 3 + 3 * F(1, 4) + F(1, 4) ** 2, f("relative_density_error_factor"))
    upper("area/lower_cubic_density_deviation", F(3), f("relative_density_error_factor"))
    upper("area/integrated_relative_error", f("relative_density_error_factor") * f("Jacobi_S_error_factor"), f("relative_area_error_factor"))
    relative_area_gap = f("relative_area_error_factor") * delta * Q ** 2
    upper("area/positive_total_area", relative_area_gap, F(1), strict=True)
    # A_static=2*pi^2*rho_b^3, A_bg=A_static*cosh(q)^3. No pi/rho value
    # or actual integral is evaluated; only this conditional scale cancels.
    static_area_multiplier = (1 + relative_area_gap) * (1 + Q ** 2) ** 3
    upper("area/static_background_volume_upper", static_area_multiplier, F(2), strict=True)
    # Positive area weights have norm one after normalization. This normalized
    # exact two-point witness also admits zero response under the same bound.
    normalized_mean = (2 * F(-1) + 3 * F(1)) / 5
    upper("scalar/positive_weight_mean_value_bound", abs(normalized_mean), F(1))

    # Restart at the exact causal shadow, supplied separately by support and
    # ODE uniqueness. A support upper bound is not an arrival assertion.
    front = (ELL / 2) / F(20, 19)
    need("causal/front_from_support_speed", front == F(19, 40) * ELL)
    elapsed_after_front = max(Q - front, F(0))
    delayed_elapsed = elapsed_after_front + f("drift_error_factor") * delta * elapsed_after_front ** 2
    # q-front is affine. Its maximum on [0,front] is its endpoint maximum.
    pre_shadow_affine_maximum = max(-front, F(0))
    need("causal/exact_zero_restart_is_separate_input", front > 0 and pre_shadow_affine_maximum == 0,
         {"analytic_dependencies": ["exact pre-shadow support zero", "background phase and ODE uniqueness through the shadow"],
          "front": front, "maximum_of_q_minus_front_on_pre_shadow_interval": pre_shadow_affine_maximum,
          "restart": "d=(q-front)+", "delayed_elapsed": "D=d+drift*delta*d^2"})
    need("causal/restart_retained_in_old_slab", delayed_elapsed < T - front)
    order = claims["Taylor_order"]
    need("Taylor/pointwise_time_jet_reserve", 14 - order >= 3,
         {"mixed_order": 14, "time_derivative_order": order, "remaining_Sobolev_order": 14 - order, "required_boundary_inclusive_embedding": 3})
    exact_kernel = F(1, order * factorial(order - 1))
    claimed_kernel = 1 / f("Taylor_kernel_denominator")
    upper("Taylor/integrated_kernel_coefficient", exact_kernel, claimed_kernel)
    # f(z)=z_+^(23/2)/C is C11 with zero jets through order11 at zero;
    # f^(11)(z)=sqrt(z)<=1 on [0,1], but no bounded twelfth derivative.
    # At z=1/4 both function values are exact rationals.
    z, sqrt_z = F(1, 4), F(1, 2)
    C11 = prod(F(23, 2) - j for j in range(11))
    f_example = z ** 11 * sqrt_z / C11
    need("Taylor/C11_zero_extension_derivative_certificate", all(F(23, 2) - j > 0 for j in range(12))
         and prod(F(23, 2) - j for j in range(12)) / C11 == F(1, 2)
         and F(23, 2) - 12 == -F(1, 2),
         {"positive_derivative_powers_through_order_eleven": [F(23, 2) - j for j in range(12)],
          "twelfth_derivative_coefficient": F(1, 2), "twelfth_derivative_power": -F(1, 2)})
    need("Taylor/C11_zero_extension_refutes_order_twelve", f_example > z ** 12 / factorial(12) and sqrt_z <= 1 and f_example <= z ** 11 / factorial(11),
         {"z": z, "eleventh_derivative": sqrt_z, "function_value": f_example,
          "false_twelfth_bound": z ** 12 / factorial(12), "valid_eleventh_bound": z ** 11 / factorial(11)})
    need("Taylor/zero_before_shadow", pre_shadow_affine_maximum ** order == 0)
    scalar_Taylor_multiplier = delayed_elapsed ** order * claimed_kernel
    upper("Taylor/envelope_factor_below_one", scalar_Taylor_multiplier, F(1), strict=True)

    zero_homogeneous_ok = (F(0) <= F(0)) if claims["homogeneous_zero_mode"] == "weak_at_zero" else (F(0) < F(0))
    need("zero/homogeneous_bounds_include_zero", zero_homogeneous_ok,
         {"endpoint_delta": 0, "correct_relation": "0<=0", "strict_relation_is_false": True})
    response = F(0)
    causal_ok = response <= 1 if claims["causal_shadow_mode"] == "shadow_zero_and_upper_bound_only" else response > 0
    need("causal/upper_bounds_do_not_imply_nonzero_arrival", causal_ok,
         {"zero_response_witness": response, "positive_upper_bound": 1,
          "actual_shadow_zero_is_a_separate_analytic_input": True})
    requested_horizon = Q if claims["duration_scope"] == "certified_Q_only" else 2 * Q
    upper("scope/no_unverified_long_duration", requested_horizon, Q)

    return {"schema": "hdblast-effective-oracle-report-v1", "status": "pass",
            "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(), "checks": records,
            "source_bound_model": {"state_fields": 16, "first_jet_directions": 5, "copy_count": 2,
                                   "intrinsic_dimension": 4, "cut_dimension": 3,
                                   "individual_member_mixed_order": 14,
                                   "parameter_field_mixed_order": 10,
                                   "C2_parameter_into_top_order_claimed": False,
                                   "numerical_parameter_derivative_bounds": False,
                                   "full_Frobenius_arrays": True, "all_angular_modes": True,
                                   "zero_perturbation_included": True,
                                   "actual_metric_or_curvature_arrays_evaluated": False},
            "exponent_DAG": exact({"P": {"operation": "symbolic_pow10", "decimal_index": P_INDEX, "lower_bound": p_floor},
                                   "amplitude_power_relative_to_P": amplitude_coefficient,
                                   "growth_exponent_relative_to_P": growth_coefficient,
                                   "M_log_upper": M_log_upper, "M_log_target": M_log_target,
                                   "metric_delta_log_upper": delta_log, "scalar_value_log_upper": scalar_log,
                                   "huge_integer_or_huge_power_materialized": False}),
            "finite_interval": exact({"ell": ELL, "Q": Q, "T": T, "qstar": qstar,
                                      "delta_range": [F(0), delta], "Q_minus_qstar": Q - qstar,
                                      "T_minus_Q": T - Q, "causal_front": front,
                                      "restarted_elapsed_at_Q": elapsed_after_front,
                                      "delayed_elapsed_at_Q": delayed_elapsed}),
            "computed_budgets": exact({"coframe": coframe_used, "DK_coefficient": DK_used,
                                       "D_connection_coefficient": DGamma_used,
                                       "global_phase_speed": speed_used,
                                       "global_phase_derivative_columns": derivative_columns,
                                       "mixed_curvature": mixed_used, "W_metric_error": W_difference,
                                       "frame_exponential_upper": frame_series_upper,
                                       "Jacobi_E": E_used, "Jacobi_S_error": S_error_used,
                                       "relative_area_error": relative_area_gap,
                                       "static_area_multiplier": static_area_multiplier,
                                       "Taylor_kernel": exact_kernel,
                                       "scalar_Taylor_multiplier_at_Q": scalar_Taylor_multiplier}),
            "scope": "Exact rational inequality directions, source norm/exponent DAG arithmetic, zero-endpoint checks and conditional scalar/area/Taylor bounds supporting the written effective-domain proof. P=10^500207 and its exponential are never materialized. These finite checks do not prove the continuum embedding, metric/source estimates, observer continuation, support theorem, Jacobi comparison or exact pre-shadow zero; those are analytic inputs. Causal shadow and upper response envelopes supply no arrival, nonzero signal, long-duration validity or experimental calibration."}


def controls(claims):
    """Fresh claim evaluation only; fresh capsule processes belong to wrappers."""
    records = []
    for field, value in MUTATIONS:
        bad = copy.deepcopy(claims)
        bad[field] = value
        changed = [key for key in claims if claims[key] != bad[key] or type(claims[key]) is not type(bad[key])]
        if changed != [field]:
            raise ValueError("fault is not single-field: " + field)
        raw, failures = canonical(bad), []
        try:
            evaluate(json.loads(raw), lambda name, passed: failures.append(name) if not passed else None)
        except ValueError as exc:
            if len(failures) != 1 or any(token in failures[0] for token in ("/schema", "/exact_engine", "/parser/", "/type/")):
                raise ValueError("fault missed its mathematical/dependency check: " + field) from exc
            records.append({"fault": field, "wrong_value": value, "changed_fields": changed,
                            "self_rehashed_sha256": hashlib.sha256(raw).hexdigest(),
                            "failed_check": failures[0], "diagnostic": str(exc)})
        else:
            raise ValueError("surviving effective fault: " + field)
    return {"schema": "hdblast-effective-oracle-controls-v1", "status": "pass", "controls": records,
            "scope": "Fresh single-field claim evaluations in one process. No fresh-subprocess or closed-capsule authentication claim is made here."}
