"""Exact finite arithmetic for the conditional nonlinear Picard limit.

Only bounded *relative* decimal powers are constructed.  Infinite Volterra,
factorial-tail, PDE-limit, and causal-support arguments remain proof premises.
The report distinguishes their symbolic numerical interfaces from finite checks.
No assertion is used, so optimized Python performs the same validation.
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import re
import sys


INTEGER_FACTORS = (
    "component_count", "spatial_dimension", "first_jet_directions", "copy_count",
    "recurrence_start", "Cstar_exponent", "coefficient_total_order",
    "embedding_exponent", "R_exponent", "admitted_ball_divisor", "star_factor",
    "p_exponent", "ell_exponent", "FV_group_count", "FP_cross_count",
    "principal_difference_count", "principal_difference_sign",
    "intermediate_Z1_multiplier", "Ldelta_multiplier", "pair_squared_factor",
    "energy_lower_denominator", "energy_upper_multiplier", "sqrtE_to_X1_multiplier",
    "X1_to_sqrtE_multiplier", "energy_production_multiplier", "root_force_multiplier",
    "energy_growth_multiplier", "volterra_multiplier", "stability_prefactor",
    "stability_rate_multiplier", "time_product_exponent", "Nstar_exponent",
    "factorial_e_ceiling", "geometric_sum_multiplier", "growth_threshold_multiplier",
    "Y_initial_exponent", "energy_unknown_exponent", "energy_predecessor_exponent",
    "normal_Y_exponent", "normal_unknown_exponent", "normal_predecessor_exponent",
    "C_J_exponent", "C_lift_exponent", "C_seed_exponent", "data_slope_exponent",
    "initial_jet_guard_exponent", "seed_amplitude_exponent", "D0_safe_exponent",
    "picard_C0_exponent", "picard_H_exponent", "smooth_approximation_factor",
    "iterated_amplitude_power_exponent", "growth_power_exponent",
    "nested_seed_test_exponent", "majorant_combined_linear_exponent",
    "nested_gap_coefficient", "majorant_size_exponent", "nested_gap_test_exponent",
    "causal_F_multiplier", "causal_forcing_multiplier", "causal_production_multiplier",
)
RATIONAL_FACTORS = (
    "point_state_radius", "p_comparison", "T_factor", "source_scalar_ceiling",
    "GT_coefficient", "betaT_coefficient", "stabilityT_coefficient",
    "first_difference_R_factor", "pair_difference_R_factor", "general_tail_R_factor",
    "twice_z_ratio", "factorial_e_ratio", "growth_ratio", "log2_lower",
    "binary_growth_power", "geometric_ratio_ceiling", "geometric_tail_R_factor",
    "factorial_binary_power", "final_binary_power", "A_ceiling", "B_upper",
    "radical_ceiling", "radial_speed", "outer_support_start", "outer_support_end",
    "outer_wall", "outer_wall_margin", "inner_gap", "inner_nonarrival_time",
)
BOOLEAN_FACTORS = (
    "global_copy_norms", "full_tensor_derivative_slots", "G_symmetric",
    "first_difference_has_PDE", "zero_initial_difference", "requires_small_betaT",
    "extra_amplitude_restriction", "varying_energy_metrics_identified", "source_C_zero",
    "finite_propagation_is_text_premise",
)
FACTORS = ("review_status", *INTEGER_FACTORS, *RATIONAL_FACTORS, *BOOLEAN_FACTORS)

# Every control changes one flat claims field and is self-rehashed before replay.
MUTATIONS = (
    ("component_count", True), ("spatial_dimension", 3), ("first_jet_directions", 4),
    ("copy_count", 1), ("global_copy_norms", False),
    ("full_tensor_derivative_slots", False), ("G_symmetric", False),
    ("first_difference_has_PDE", True), ("zero_initial_difference", False),
    ("requires_small_betaT", True), ("extra_amplitude_restriction", True),
    ("varying_energy_metrics_identified", True), ("recurrence_start", 0),
    ("Cstar_exponent", 499999), ("coefficient_total_order", 15),
    ("embedding_exponent", 6), ("R_exponent", -11), ("admitted_ball_divisor", 2),
    ("star_factor", 1), ("p_exponent", -6), ("point_state_radius", "1/10000"),
    ("p_comparison", "1/1000000"), ("ell_exponent", -5), ("T_factor", "3/2"),
    ("FV_group_count", 2), ("FP_cross_count", 1),
    ("principal_difference_count", 1), ("principal_difference_sign", -1),
    ("source_scalar_ceiling", "1"), ("intermediate_Z1_multiplier", 1),
    ("Ldelta_multiplier", 1), ("pair_squared_factor", 1),
    ("energy_lower_denominator", 2), ("energy_upper_multiplier", 0),
    ("sqrtE_to_X1_multiplier", 0), ("X1_to_sqrtE_multiplier", 2),
    ("energy_production_multiplier", 20), ("root_force_multiplier", 0),
    ("energy_growth_multiplier", 10), ("volterra_multiplier", 11),
    ("stability_prefactor", 2), ("stability_rate_multiplier", 27),
    ("time_product_exponent", 499993), ("GT_coefficient", "25"),
    ("betaT_coefficient", "19"), ("stabilityT_coefficient", "45"),
    ("first_difference_R_factor", "1/4"), ("pair_difference_R_factor", "1/4"),
    ("general_tail_R_factor", "1/4"), ("Nstar_exponent", 499995),
    ("twice_z_ratio", "39/10000"), ("factorial_e_ceiling", 2),
    ("factorial_e_ratio", "117/20000"), ("log2_lower", "1"),
    ("binary_growth_power", "1/4"), ("growth_threshold_multiplier", 3),
    ("growth_ratio", "13/1250"), ("geometric_ratio_ceiling", "1/4"),
    ("geometric_sum_multiplier", 1), ("geometric_tail_R_factor", "1/2"),
    ("factorial_binary_power", "2"), ("final_binary_power", "3/4"),
    ("Y_initial_exponent", 29), ("energy_unknown_exponent", 500100),
    ("energy_predecessor_exponent", 500099), ("normal_Y_exponent", 109),
    ("normal_unknown_exponent", 500207), ("normal_predecessor_exponent", 500207),
    ("C_J_exponent", 6502999), ("C_lift_exponent", 254), ("C_seed_exponent", 6503254),
    ("data_slope_exponent", 4199), ("initial_jet_guard_exponent", -6503019),
    ("seed_amplitude_exponent", -6507484), ("D0_safe_exponent", 7003208),
    ("picard_C0_exponent", 7003208), ("picard_H_exponent", 500211),
    ("smooth_approximation_factor", 1), ("iterated_amplitude_power_exponent", 500206),
    ("growth_power_exponent", 500206), ("nested_seed_test_exponent", 6),
    ("majorant_combined_linear_exponent", 7007409), ("nested_gap_coefficient", 98),
    ("majorant_size_exponent", -12), ("nested_gap_test_exponent", 4),
    ("A_ceiling", "1/20"), ("B_upper", "1"), ("radical_ceiling", "2"),
    ("radial_speed", "1"), ("outer_support_start", "2"),
    ("outer_support_end", "4"), ("outer_wall", "5"),
    ("outer_wall_margin", "2"), ("inner_gap", "1/4"),
    ("inner_nonarrival_time", "1/2"), ("source_C_zero", False),
    ("causal_F_multiplier", 1), ("causal_forcing_multiplier", 7),
    ("causal_production_multiplier", 28), ("finite_propagation_is_text_premise", False),
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def pow10(exponent):
    """Bounded ordinary relative powers only, never a huge or nested power."""
    if type(exponent) is not int or abs(exponent) > 128:
        raise ValueError("large decimal expansion forbidden")
    return F(10**exponent) if exponent >= 0 else F(1, 10**(-exponent))


def rational(text):
    """No scientific notation, enormous integers, nonfinite values, or booleans."""
    if (type(text) is not str or len(text) > 80
            or re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/[1-9][0-9]*)?", text) is None):
        raise ValueError("bounded canonical rational string required")
    value = F(text)
    if str(value) != text:
        raise ValueError("rational string must be reduced and canonical")
    return value


def positive_sum_below(terms, exponent, strict=True):
    """Remove a common power, rounding very small positive terms upward.

    Relative powers below -64 are replaced by -64, an upper bound.  No
    cancellation occurs and no absolute large power is constructed.
    """
    if type(exponent) is not int:
        raise ValueError("integer output exponent required")
    total = F(0)
    for coefficient, power in terms:
        coefficient = F(coefficient)
        if coefficient < 0 or type(power) is not int:
            raise ValueError("nonnegative exact monomial required")
        if coefficient == 0:
            continue
        relative = power-exponent
        if relative > 128:
            return False
        total += coefficient*pow10(max(relative, -64))
    return total < 1 if strict else total <= 1


def evaluate(claims, check):
    records = []

    def need(name, condition):
        boolean = bool(condition)
        check("convergence_arithmetic/"+name, boolean)
        records.append({"name": name, "status": "pass" if boolean else "fail"})
        if not boolean:
            raise ValueError("convergence arithmetic rejects "+name)

    need("schema_and_exact_keys", type(claims) is dict
         and claims.get("schema") == "hdblast-convergence-arithmetic-claims-v1"
         and set(claims) == {"schema", *FACTORS})
    need("exact_bounded_integer_types", all(type(claims[k]) is int
         and abs(claims[k]) <= 20_000_000 for k in INTEGER_FACTORS))
    need("exact_boolean_types", all(type(claims[k]) is bool for k in BOOLEAN_FACTORS))
    need("exact_rational_string_types", all(type(claims[k]) is str for k in RATIONAL_FACTORS))
    q = {k: rational(claims[k]) for k in RATIONAL_FACTORS}
    need("canonical_rational_values", all(str(q[k]) == claims[k] for k in q))
    need("review_status", type(claims["review_status"]) is str
         and claims["review_status"] == "accepted_conditional_analytic_inputs")
    need("full_bundle_and_five_jet_slots", claims["component_count"] == 16
         and claims["spatial_dimension"] == 4 and claims["first_jet_directions"] == 5
         and claims["first_jet_directions"] == 1+claims["spatial_dimension"])
    need("global_one_and_identical_two_copy_convention", claims["copy_count"] == 2
         and claims["global_copy_norms"] and claims["full_tensor_derivative_slots"])
    need("source_symmetry_and_zero", claims["G_symmetric"] and claims["source_C_zero"])
    need("correct_iteration_start_and_seed_scope", claims["recurrence_start"] == 1
         and not claims["first_difference_has_PDE"] and claims["zero_initial_difference"])
    need("no_unjustified_common_energy_metric", not claims["varying_energy_metrics_identified"])
    need("no_small_kernel_or_amplitude_hypothesis", not claims["requires_small_betaT"]
         and not claims["extra_amplitude_restriction"])

    # These are inherited analytic interfaces, not newly derived estimates.
    cs = claims["Cstar_exponent"]; radius = claims["R_exponent"]
    embed = claims["embedding_exponent"]; ball = claims["admitted_ball_divisor"]
    need("retained_coefficient_and_embedding_interfaces", cs == 500000
         and claims["coefficient_total_order"] == 16 and embed == 7)
    need("retained_admitted_ball", radius == -12 and ball == 4)
    need("same_actual_slab", claims["ell_exponent"] == -6 and q["T_factor"] == F(13,8))
    need("same_state_patch", q["point_state_radius"] == F(1,1000)
         and claims["star_factor"] == 2)
    pexp = claims["p_exponent"]
    need("pointwise_jet_exponent", pexp == embed+radius)
    p = pow10(pexp)
    need("joint_segment_in_state_patch", claims["star_factor"]*p < q["point_state_radius"])
    need("small_pointwise_jet_comparison", p <= q["p_comparison"] <= F(1,100))
    need("complete_source_derivative_groups", claims["FV_group_count"] == 3
         and claims["FP_cross_count"] == 2)
    need("both_positive_principal_differences", claims["principal_difference_count"] == 2
         and claims["principal_difference_sign"] == 1)
    fv = 1+p+p*p
    fp = 1+claims["FP_cross_count"]*p
    state = fv+claims["principal_difference_count"]*p
    ceiling = q["source_scalar_ceiling"]
    need("joint_mean_value_state_coefficient", state == 1+3*p+p*p and state <= ceiling)
    need("joint_mean_value_full_jet_coefficient", fp == 1+2*p and fp <= ceiling)
    need("eleven_tenths_source_ceiling", ceiling == F(11,10))
    zfactor = claims["intermediate_Z1_multiplier"]
    need("source_Cauchy_Schwarz_without_roots", zfactor > 0 and 2*ceiling**2 <= zfactor**2)
    ld = claims["Ldelta_multiplier"]
    need("recorded_Ldelta", ld == 4 and ld >= zfactor)

    ed = claims["energy_lower_denominator"]
    eupper = claims["energy_upper_multiplier"]
    pair = claims["pair_squared_factor"]
    need("retained_energy_and_pair_comparisons", ed == 4 and eupper == 1 and pair == 2)
    upper = claims["sqrtE_to_X1_multiplier"]; lower = claims["X1_to_sqrtE_multiplier"]
    need("sqrtE_bounded_by_pair", upper > 0 and upper**2 >= eupper)
    need("pair_bounded_by_three_sqrtE", lower == 3 and lower**2 >= pair*ed)
    prod = claims["energy_production_multiplier"]
    need("retained_energy_production", prod == 21)
    rootforce = claims["root_force_multiplier"]
    growth = claims["energy_growth_multiplier"]
    need("regularized_root_force_rounding", rootforce == 1 and 2*rootforce**2 >= 1)
    need("regularized_root_growth_rounding", growth == 16 and 2*growth >= prod)
    beta = claims["volterra_multiplier"]
    need("varying_metric_volterra_conversion", beta == lower*rootforce*ld == 12)
    pref = claims["stability_prefactor"]; rate = claims["stability_rate_multiplier"]
    need("nonzero_initial_pair_stability_prefactor", pref == lower*upper == 3)
    need("nonzero_initial_pair_stability_rate", rate == growth+beta == 28)

    tp = claims["time_product_exponent"]; tf = q["T_factor"]
    need("full_slab_symbolic_power", tp == cs+claims["ell_exponent"])
    g = q["GT_coefficient"]; z = q["betaT_coefficient"]; total = q["stabilityT_coefficient"]
    need("full_slab_energy_growth", g == growth*tf == 26)
    need("full_slab_volterra_size", z == beta*tf == F(39,2))
    need("full_slab_stability_exponent", total == g+z == rate*tf == F(91,2))
    need("kernel_need_not_be_small", tp >= 0 and z > 1)
    start = q["first_difference_R_factor"]
    need("seed_difference_from_admission_only", start == F(2,ball))
    need("independent_pair_tail_ceiling", q["pair_difference_R_factor"] == F(2,ball))
    need("general_factorial_tail_prefactor", q["general_tail_R_factor"] == start)
    nstar = claims["Nstar_exponent"]
    need("recorded_sufficient_index", nstar == 499997)
    relative = pow10(tp-nstar)
    need("geometric_tail_threshold", q["twice_z_ratio"] == 2*z*relative == F(39,1000)
         and 0 < q["twice_z_ratio"] < 1)
    eceil = claims["factorial_e_ceiling"]
    need("recorded_Euler_number_upper_premise", eceil == 3)
    need("factorial_binary_threshold", q["factorial_e_ratio"] == eceil*z*relative == F(117,2000)
         and 0 < q["factorial_e_ratio"] < F(1,2))
    loglower = q["log2_lower"]; binarygrowth = q["binary_growth_power"]
    need("recorded_log_two_lower_premise", loglower == F(1,2))
    need("recorded_binary_growth_allocation", binarygrowth == F(1,2))
    gmult = claims["growth_threshold_multiplier"]
    need("exponential_to_binary_threshold_factor", gmult*loglower*binarygrowth == 1)
    need("exponential_to_binary_numeric_threshold", q["growth_ratio"] == gmult*g*relative == F(104,1000)
         and 0 < q["growth_ratio"] < 1)
    gratio = q["geometric_ratio_ceiling"]
    need("tail_successive_ratio_ceiling", gratio == F(1,2)
         and z*relative < gratio < 1)
    gsum = claims["geometric_sum_multiplier"]
    need("geometric_series_multiplier", gsum == 1/(1-gratio))
    need("geometric_tail_prefactor", q["geometric_tail_R_factor"] == start*gsum == 1)
    need("factorial_binary_decay_power", q["factorial_binary_power"] == 1)
    need("final_uniform_binary_tail_power", q["final_binary_power"]
         == q["factorial_binary_power"]-binarygrowth == F(1,2))

    # Rejoin the six accepted outputs and the SAME original amplitude.
    six = tuple(claims[k] for k in ("Y_initial_exponent", "energy_unknown_exponent",
        "energy_predecessor_exponent", "normal_Y_exponent", "normal_unknown_exponent",
        "normal_predecessor_exponent"))
    need("six_inherited_frozen_outputs", six == (30,500101,500100,110,500208,500208))
    cy, ku, kv, nn, lu, lv = six
    cj = claims["C_J_exponent"]; lift = claims["C_lift_exponent"]
    seed = claims["C_seed_exponent"]
    need("inherited_initial_jet_and_lift", cj == 6503000 and lift == 255)
    need("joined_compatible_seed_exponent", seed == cj+lift)
    slope = claims["data_slope_exponent"]; guard = claims["initial_jet_guard_exponent"]
    seedamp = claims["seed_amplitude_exponent"]
    need("actual_parameter_data_slope", slope == 4200)
    need("inherited_initial_jet_guard", guard == -6503020)
    need("inherited_compatible_seed_interval", seedamp == -6507485)
    d0 = claims["D0_safe_exponent"]
    need("D0_majorant_from_six_outputs", positive_sum_below(
        [(1,nn+cy+cj),(1,lu+cj),(1,lv+cj)],d0))
    c0 = claims["picard_C0_exponent"]; hh = claims["picard_H_exponent"]
    need("C0_covers_initial_and_seed", c0 == 7003210 and d0 < c0 and seed < c0)
    need("Hu_majorant", positive_sum_below([(1,nn+ku),(1,lu)],hh))
    need("Hv_majorant", positive_sum_below([(1,nn+kv),(1,lv)],hh))
    need("retained_common_H", hh == 500212)
    approx = claims["smooth_approximation_factor"]
    need("smooth_supported_approximation_factor", approx == 2)
    ap = claims["iterated_amplitude_power_exponent"]
    need("unchanged_original_amplitude", ap == 500209)
    seedtest = claims["nested_seed_test_exponent"]
    need("bounded_nested_seed_test", 0 <= seedtest <= 8 and ap >= seedtest)
    testpower = pow10(seedtest)
    need("iterated_interval_inside_seed_interval", testpower > -seedamp)
    need("approximate_data_inside_initial_jet_guard", positive_sum_below(
        [(approx,slope-int(testpower))],guard))
    gp = claims["growth_power_exponent"]
    need("nested_growth_exponent", gp == 500207)
    need("full_slab_two_H_growth", positive_sum_below([(2*tf,hh+claims["ell_exponent"])],gp))
    need("Euler_ten_premise_sufficient", eceil < 10)
    linear = claims["majorant_combined_linear_exponent"]
    need("linear_admission_exponent_sum", linear == c0+slope)
    gap = claims["nested_gap_coefficient"]
    need("nested_power_gap_exact", gap == pow10(ap-gp)-1 == 99)
    size = claims["majorant_size_exponent"]; gaptest = claims["nested_gap_test_exponent"]
    need("bounded_nested_gap_test", 0 <= gaptest <= 8 and gp >= gaptest)
    need("nested_gap_dominates_linear_remainder", gap*pow10(gaptest) > linear-size)
    need("admitted_quarter_ball_margin", positive_sum_below([(ball*approx,size)],radius))

    # Finite translation of the separate causal-support proof (no PDE proof here).
    aa = q["A_ceiling"]; bb = q["B_upper"]; radical = q["radical_ceiling"]
    speed = q["radial_speed"]
    need("retained_normalized_speed_coefficients", aa == F(1,19) and bb == F(20,19))
    need("exact_speed_discriminant", radical > 0 and radical**2 == aa**2+4*bb)
    need("positive_characteristic_speed", speed == (aa+radical)/2 == F(20,19))
    need("retained_support_and_wall_locations", q["outer_support_start"] == F(5,2)
         and q["outer_wall"] == 6 and q["inner_gap"] == F(1,2))
    need("outer_support_at_full_slab", q["outer_support_end"]
         == q["outer_support_start"]+speed*tf == F(80,19))
    need("strict_outer_wall_margin", q["outer_wall_margin"]
         == q["outer_wall"]-q["outer_support_end"] == F(34,19) and q["outer_wall_margin"] > 0)
    need("inner_gap_nonarrival_lower_bound_only", q["inner_nonarrival_time"]
         == q["inner_gap"]/speed == F(19,40))
    cforce = claims["causal_F_multiplier"]
    need("nonlinear_zero_state_source_point_bound", cforce > 0
         and cforce**2 >= 1+(1+p)**2)
    cpair = claims["causal_forcing_multiplier"]
    need("source_pairing_energy_conversion", cpair >= cforce*ed)
    need("local_energy_production_rounding", claims["causal_production_multiplier"] >= prod+cpair)
    need("causal_support_argument_is_not_arithmetic", claims["finite_propagation_is_text_premise"])

    # These checks exercise the expansion guard without ever expanding a power.
    rejected = []
    for forbidden in (cs, nstar, ap, True):
        try:
            pow10(forbidden)
        except ValueError:
            rejected.append(forbidden)
    need("huge_absolute_and_boolean_power_requests_rejected", len(rejected) == 4)
    need("positive_tiny_term_rounding_is_upward", positive_sum_below([(1,-1_000_000)],0)
         and not positive_sum_below([(1,0),(1,-1_000_000)],0,strict=False))

    return {
        "schema": "hdblast-convergence-arithmetic-report-v1", "status": "pass",
        "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(),
        "named_check_count": len(records), "checks": records,
        "exact_rationals": {"p": str(p), "F_V_coefficient": str(fv),
            "F_Q_coefficient": str(fp), "full_source_state_coefficient": str(state),
            "twice_z_over_Nstar": str(q["twice_z_ratio"]),
            "three_z_over_Nstar": str(q["factorial_e_ratio"]),
            "four_g_over_Nstar": str(q["growth_ratio"]),
            "radial_speed": str(speed), "outer_wall_margin_in_ell": str(q["outer_wall_margin"])},
        "symbolic_interfaces": {
            "L_delta": {"coefficient": ld, "decimal_exponent": cs},
            "beta": {"coefficient": beta, "decimal_exponent": cs},
            "G_en": {"coefficient": growth, "decimal_exponent": cs},
            "S_low": {"prefactor": pref, "exp_argument": {"coefficient": str(total), "decimal_exponent": tp}},
            "individual_difference": "(R/2) exp(G_en s) (beta s)^n/n!, n>=0",
            "finite_uniform_tail": "min(R/2,(R/2) exp((G_en+beta)T) (beta T)^N/N!), M>N>=0",
            "N_star": {"decimal_exponent": nstar},
            "large_index_uniform_tail": "R 2^(-N/2), M>N>=N_star",
            "same_original_amplitude": {"base": 10, "outer_exponent_sign": -1,
                "inner_power_base": 10, "inner_power_exponent": ap},
            "causal_speed": str(speed),
            "inner_nonarrival_lower_bound_in_ell": str(q["inner_nonarrival_time"]),
        },
        "analytic_premises_not_proved_by_this_replay": [
            "Authenticated coefficient arrays, tensor embedding, admitted smooth iterates and homogeneous boundary rows.",
            "Joint mean-value source identity, energy identity and legitimate regularized-root integration.",
            "Volterra induction for every n, infinite factorial-tail summation, factorial integral bound, e<3 and log(2)>=1/2.",
            "Banach completion, finite-Sobolev PDE identification, top regularity and nonlinear uniqueness.",
            "Symmetric normal-flux finite propagation with all cross slots and zero background source.",
        ],
        "scope": "Exact finite rational translation only; symbolic large powers; no new amplitude restriction or analytic theorem.",
    }


def controls(claims):
    if len({field for field, _ in MUTATIONS}) != len(MUTATIONS):
        raise ValueError("mutation fields must be distinct and flat")
    results = []
    for field, value in MUTATIONS:
        changed = copy.deepcopy(claims)
        changed[field] = value
        payload = canonical(changed)
        fresh = json.loads(payload)
        try:
            evaluate(fresh, lambda name, ok: None)
        except (ValueError, TypeError, KeyError) as error:
            results.append({"field": field, "value": value, "status": "rejected",
                "self_rehashed_sha256": hashlib.sha256(payload).hexdigest(), "diagnostic": str(error)})
        else:
            raise ValueError("arithmetic mutation survived: "+field)
    return {"schema": "hdblast-convergence-arithmetic-controls-v1", "status": "pass",
        "control_count": len(results), "controls": results}


def _pairs(items):
    result = {}
    for key, value in items:
        if key in result:
            raise ValueError("duplicate JSON key: "+key)
        result[key] = value
    return result


def _invalid_constant(value):
    raise ValueError("nonfinite JSON token forbidden: "+value)
