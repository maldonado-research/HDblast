"""Exact finite ledger for conditional frozen energy and normal recovery.

This checks finite rational margins, word counts, scale allocations and the
conditional smooth-Picard majorant. It does not prove the analytic product,
embedding, frozen-existence or convergence premises. Nested decimal powers
are compared symbolically, never constructed as large integers or decimals.
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial
from pathlib import Path
import sys


FACTORS = (
    "review_status", "component_count", "full_tensor_derivative_slots",
    "source_C_zero", "unknown_ball_required", "inverse_derivatives_required",
    "Cstar_exponent", "coefficient_total_order", "predecessor_radius_exponent",
    "embedding_exponent", "ell", "T_factor", "retained_copy_count",
    "point_state_radius", "metric_patch_radius", "inverse_difference_ceiling",
    "cosh_ceiling", "A_ceiling", "B_lower", "B_upper", "Bzz_lower",
    "Bzz_reciprocal_ceiling", "energy_lower_denominator", "mixed_order",
    "max_original_tensor_rank", "rotation_generator_count", "rotation_bridge_factor",
    "inverse_angular_bridge", "rotation_decimal_exponent", "word_count_sqrt_exponent",
    "Y_initial_exponent", "coefficient_low_leaf_exponent", "coefficient_high_leaf_exponent",
    "coefficient_low_infinity_exponent", "coefficient_high_L2_exponent",
    "source_leaf_exponent", "source_covariant_exponent", "source_per_word_exponent",
    "source_aggregate_exponent", "commutator_subset_exponent", "unknown_low_infinity_exponent",
    "energy_commutator_per_word_exponent", "energy_commutator_aggregate_exponent",
    "normal_commutator_per_word_exponent", "product_exponent", "energy_derivative_coefficient",
    "energy_divergence_coefficient", "energy_mixed_production_coefficient",
    "energy_production_coefficient", "energy_aggregate_production_coefficient",
    "energy_unknown_exponent", "energy_predecessor_exponent", "normal_base_squared_cost",
    "normal_base_factor", "normal_leading_cost", "normal_stage_count",
    "normal_row_sqrt_exponent", "normal_step_exponent", "normal_lower_step_exponent",
    "normal_geometric_sum_exponent", "normal_W_Y_exponent", "normal_W_lower_exponent",
    "mixed_norm_bridge_factor", "mixed_norm_bridge_exponent", "normal_Y_exponent",
    "normal_unknown_exponent", "normal_predecessor_exponent", "C_J_exponent",
    "C_lift_exponent", "C_seed_exponent", "data_d_slope_exponent", "initial_jet_guard_exponent",
    "seed_amplitude_exponent", "seed_data_exponent", "seed_I14_exponent", "seed_X14_exponent",
    "D0_safe_exponent", "picard_C0_exponent", "picard_H_exponent", "smooth_approximation_factor",
    "iterated_amplitude_power_exponent", "growth_power_exponent", "nested_seed_comparison_exponent",
    "majorant_combined_linear_exponent", "nested_gap_coefficient", "majorant_size_exponent",
    "nested_gap_test_exponent",
)
RATIONAL_FACTORS = (
    "ell", "T_factor", "point_state_radius", "metric_patch_radius",
    "inverse_difference_ceiling", "cosh_ceiling", "A_ceiling", "B_lower",
    "B_upper", "Bzz_lower", "Bzz_reciprocal_ceiling",
)
BOOLEAN_FACTORS = (
    "full_tensor_derivative_slots", "source_C_zero", "unknown_ball_required",
    "inverse_derivatives_required",
)


MUTATIONS = (
    ("Cstar_exponent", 499999),
    ("coefficient_total_order", 15),
    ("predecessor_radius_exponent", -11),
    ("embedding_exponent", 6),
    ("point_state_radius", "1/100000"),
    ("metric_patch_radius", "1/39"),
    ("inverse_difference_ceiling", "1/40"),
    ("cosh_ceiling", "9/8"),
    ("A_ceiling", "1/20"),
    ("B_lower", "3/4"),
    ("B_upper", "21/20"),
    ("Bzz_lower", "24/25"),
    ("Bzz_reciprocal_ceiling", "1"),
    ("energy_lower_denominator", 2),
    ("rotation_bridge_factor", 39),
    ("rotation_decimal_exponent", 22),
    ("coefficient_low_leaf_exponent", 8),
    ("coefficient_high_leaf_exponent", 11),
    ("coefficient_low_infinity_exponent", 39),
    ("coefficient_high_L2_exponent", 34),
    ("source_leaf_exponent", 13),
    ("source_covariant_exponent", 13),
    ("source_per_word_exponent", 37),
    ("source_aggregate_exponent", 43),
    ("commutator_subset_exponent", 4),
    ("unknown_low_infinity_exponent", 29),
    ("energy_commutator_per_word_exponent", 69),
    ("energy_commutator_aggregate_exponent", 76),
    ("normal_commutator_per_word_exponent", 67),
    ("product_exponent", 500076),
    ("word_count_sqrt_exponent", 5),
    ("energy_derivative_coefficient", 1),
    ("energy_divergence_coefficient", 3),
    ("energy_mixed_production_coefficient", 11),
    ("energy_production_coefficient", 20),
    ("energy_aggregate_production_coefficient", 10),
    ("Y_initial_exponent", 21),
    ("energy_unknown_exponent", 500100),
    ("energy_predecessor_exponent", 500099),
    ("normal_base_squared_cost", 111),
    ("normal_base_factor", 10),
    ("normal_leading_cost", 10),
    ("normal_stage_count", 12),
    ("normal_row_sqrt_exponent", 4),
    ("normal_step_exponent", 7),
    ("normal_lower_step_exponent", 6),
    ("normal_geometric_sum_exponent", 95),
    ("normal_W_Y_exponent", 105),
    ("normal_W_lower_exponent", 103),
    ("mixed_norm_bridge_factor", 400),
    ("mixed_norm_bridge_exponent", 3),
    ("normal_Y_exponent", 109),
    ("normal_unknown_exponent", 500207),
    ("normal_predecessor_exponent", 500207),
    ("C_J_exponent", 6502999),
    ("C_lift_exponent", 254),
    ("C_seed_exponent", 6503254),
    ("data_d_slope_exponent", 4199),
    ("initial_jet_guard_exponent", -6503019),
    ("seed_amplitude_exponent", -6507484),
    ("seed_data_exponent", -6503286),
    ("seed_I14_exponent", -286),
    ("seed_X14_exponent", -31),
    ("D0_safe_exponent", 7003208),
    ("picard_C0_exponent", 7003208),
    ("picard_H_exponent", 500211),
    ("smooth_approximation_factor", 1),
    ("iterated_amplitude_power_exponent", 500206),
    ("growth_power_exponent", 500206),
    ("majorant_combined_linear_exponent", 7007409),
    ("nested_gap_coefficient", 98),
    ("majorant_size_exponent", -12),
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def pow10(exponent):
    """Only bounded ordinary powers; nested powers never enter this helper."""
    if type(exponent) is not int or abs(exponent) > 1000:
        raise ValueError("large decimal expansion forbidden")
    return F(10**exponent) if exponent >= 0 else F(1, 10**(-exponent))


def positive_sum_below(terms, exponent, strict=True):
    """A conservative rational upper sum after removing a common power.

    A relative exponent below -100 is rounded UP to -100. Thus a passing
    comparison proves the actual positive-monomial comparison, while huge
    exponents are never expanded. There is no cancellation in this helper.
    """
    total = F(0)
    if type(exponent) is not int:
        raise ValueError("integer output exponent required")
    for coefficient, power in terms:
        coefficient = F(coefficient)
        if coefficient < 0 or type(power) is not int:
            raise ValueError("nonnegative exact monomial required")
        relative = power-exponent
        if relative > 1000:
            return False
        total += coefficient*pow10(max(relative, -100))
    return total < 1 if strict else total <= 1


def evaluate(claims, check):
    records = []

    def need(name, condition):
        boolean = bool(condition)
        check("energy_arithmetic/"+name, boolean)
        records.append({"name": name, "status": "pass" if boolean else "fail"})
        if not boolean:
            raise ValueError("energy arithmetic rejects "+name)

    need("schema", isinstance(claims,dict)
         and claims.get("schema") == "hdblast-variable-energy-arithmetic-claims-v1"
         and set(claims) == {"schema",*FACTORS})
    need("exact_integer_ledger_types", all(type(claims[field]) is int for field in FACTORS
         if field not in {*RATIONAL_FACTORS,*BOOLEAN_FACTORS,"review_status"}))
    need("rational_string_ledger_types", all(type(claims[field]) is str for field in RATIONAL_FACTORS))
    need("boolean_ledger_types", all(type(claims[field]) is bool for field in BOOLEAN_FACTORS))
    need("review_status", claims["review_status"] == "accepted_conditional_analytic_inputs")
    need("scalar_principal_component_count", claims["component_count"] == 16)
    need("full_covariant_slot_convention", claims["full_tensor_derivative_slots"] is True)
    need("zero_state_source_extraction", claims["source_C_zero"] is True)
    need("no_unknown_ball_hypothesis", claims["unknown_ball_required"] is False)
    need("interior_divide_after_differentiation", claims["inverse_derivatives_required"] is False)

    # Authenticated source interfaces are distinguished from derived ceilings.
    cs = claims["Cstar_exponent"]; order = claims["coefficient_total_order"]
    radius = claims["predecessor_radius_exponent"]; embed = claims["embedding_exponent"]
    need("retained_full_array_coefficient_bound", cs == 500000)
    need("retained_total_coefficient_order", order == 16)
    need("retained_predecessor_radius", radius == -12)
    need("retained_full_slot_embedding", embed == 7)
    ell = F(claims["ell"]); tfactor = F(claims["T_factor"])
    need("same_prescribed_slab", ell == F(1,10**6) and tfactor == F(13,8))
    need("slab_within_half", 0 < tfactor*ell < F(1,2))
    copies = claims["retained_copy_count"]
    need("retained_copy_convention", copies == 2)
    # pi^2<10 is sufficient and charges the entire retained two-copy volume.
    need("two_copy_volume_below_one", copies*6*ell*2*10 < 1)
    need("same_point_state_patch", F(claims["point_state_radius"]) == F(1,1000))
    need("predecessor_in_actual_point_patch", positive_sum_below(
        [(2,embed+radius)], -3))
    need("low_field_factor_less_than_one", positive_sum_below([(1,embed+radius)],0))

    # Same-point inverse ratios: rho^{-2} cancels BEFORE taking extrema.
    eps = F(claims["metric_patch_radius"])
    diff = F(claims["inverse_difference_ceiling"])
    need("same_all_state_metric_patch", eps == F(1,40))
    need("inverse_neumann_difference", diff >= eps/(1-eps))
    ilower = 1-diff; iupper = 1+diff
    need("inverse_signed_diagonal_margins", ilower == F(38,39) and iupper == F(40,39))
    cosh = F(claims["cosh_ceiling"])
    # cosh(t)<=sum_k(t^2/2)^k<=1/(1-1/8)=8/7 for |t|<=1/2.
    need("factorial_cosh_geometric_majorant", all(factorial(2*k) >= 2**k for k in range(19)))
    need("cosh_margin", cosh >= 1/(1-F(1,8)))
    A = F(claims["A_ceiling"]); Blo = F(claims["B_lower"])
    Bhi = F(claims["B_upper"]); Bzzlo = F(claims["Bzz_lower"])
    recip = F(claims["Bzz_reciprocal_ceiling"])
    need("normalized_A_ratio", A >= 2*diff/ilower)
    need("fixed_angular_B_lower", Blo <= ilower/iupper/cosh**2)
    need("normalized_B_upper", Bhi >= iupper/ilower)
    need("radial_Bzz_lower", Bzzlo <= ilower/iupper)
    need("zeroth_Bzz_reciprocal_only", recip >= 1/Bzzlo and recip < 2)
    need("coarse_principal_margins", A < 1 and Blo > F(1,2) and Bhi < 2)
    Eden = claims["energy_lower_denominator"]
    need("energy_coercive_comparison", F(1,Eden) <= min(F(1,2),Blo/2))
    need("energy_upper_comparison", max(F(1,2),Bhi/2) <= 1)

    # Complete natural SO(4) lifts: six ordered generators, no cutoff.
    m = claims["mixed_order"]; r = claims["max_original_tensor_rank"]
    rotations = claims["rotation_generator_count"]
    bridge = claims["rotation_bridge_factor"]
    need("actual_orders_and_original_rank", m == 14 and r == 2 and rotations == 6)
    need("all_rotated_derivative_slots", r+m <= 16)
    need("natural_rotation_upper_induction", 1+6*16**2 < bridge**2 and bridge >= 40)
    need("inverse_angular_bridge_constant_one", claims["inverse_angular_bridge"] == 1)
    drot = bridge**m; dexp = claims["rotation_decimal_exponent"]
    need("rotation_exact_integer", drot == 26843545600000000000000)
    need("rotation_decimal_rounding", drot < 10**dexp)
    energy_order = m-1; normal_order = m-2
    need("normal_and_energy_word_orders", energy_order == 13 and normal_order == 12)
    nwords = sum((m-k)*rotations**k for k in range(m))
    need("energy_exact_ordered_word_count", nwords == 18807399380)
    sqrtwords = claims["word_count_sqrt_exponent"]
    need("energy_word_square_sum_rounding", nwords < 10**(2*sqrtwords))
    need("energy_time_length_pairs", sum(range(1,m+1)) == 105)
    need("rotation_energy_aggregate_rounding", 315*bridge**26 < 10**46)
    cy = claims["Y_initial_exponent"]
    need("Y_initial_and_word_rounding", 315*bridge**26 < 10**(2*cy))

    # Ordered covariant chain leaves, followed by low-rank output rotations.
    lowleaf = claims["coefficient_low_leaf_exponent"]
    highleaf = claims["coefficient_high_leaf_exponent"]
    need("coefficient_ordered_leaf_counts", all(factorial(n+1) >= 1 for n in range(14)))
    need("coefficient_low_Sobolev_leaf_rounding", all(
        (n+1)*factorial(n+1)**2 < 10**(2*lowleaf) for n in range(11)))
    need("coefficient_high_Sobolev_leaf_rounding", all(
        (n+1)*factorial(n+1)**2 < 10**(2*highleaf) for n in range(14)))
    lowcoef = claims["coefficient_low_infinity_exponent"]
    highcoef = claims["coefficient_high_L2_exponent"]
    need("coefficient_low_infinity_rounding", positive_sum_below(
        [(1,embed+dexp+lowleaf),(1,embed+dexp+lowleaf+radius)],lowcoef))
    need("coefficient_high_L2_rounding", dexp+highleaf <= highcoef)
    need("low_coefficient_embedding_degree", 7+3 <= 10 and r+7+3 <= 16)
    need("source_mean_value_certificate_order", energy_order+1 <= order)
    sourceleaves = 3*factorial(energy_order+3)//factorial(3)
    need("source_exact_three_group_leaf_count", sourceleaves == 10461394944000)
    sleaf = claims["source_leaf_exponent"]
    need("source_leaf_rounding", sourceleaves < 10**sleaf)
    scov = claims["source_covariant_exponent"]
    need("source_covariant_hierarchy_rounding", 14*sourceleaves**2 < 10**(2*scov))
    need("source_single_high_factor", energy_order+2 <= 15 and 2*8 > energy_order+2)
    need("source_low_field_embedding_rank", r+7+3 <= 16)
    need("normal_source_individual_field_order", normal_order+1 <= 13)
    sper = claims["source_per_word_exponent"]
    sagg = claims["source_aggregate_exponent"]
    need("source_per_word_rotation_rounding", scov+dexp <= sper)
    need("source_word_square_sum_rounding", sper+sqrtwords <= sagg)
    subset = 2*(2**energy_order-1)
    subsetexp = claims["commutator_subset_exponent"]
    need("actual_two_principal_subset_count", subset == 16382)
    need("commutator_subset_rounding", subset < 10**subsetexp)
    need("principal_coefficient_hit_unknown_order", all(
        n-q+2 <= 14 for n in range(energy_order+1) for q in range(1,n+1)))
    need("normal_coefficient_hit_unknown_order", all(
        n-q+2 <= 13 for n in range(normal_order+1) for q in range(1,n+1)))
    need("energy_high_coefficient_low_unknown", all(
        n-q+2 <= 7 for n in range(energy_order+1) for q in range(8,n+1)))
    need("normal_high_coefficient_low_unknown", all(
        n-q+2 <= 6 for n in range(normal_order+1) for q in range(8,n+1)))
    uinf = claims["unknown_low_infinity_exponent"]
    need("unknown_low_infinity_rounding", embed+dexp <= uinf)
    commper = claims["energy_commutator_per_word_exponent"]
    need("energy_commutator_per_word_rounding", max(lowcoef+dexp,highcoef+uinf)+subsetexp <= commper)
    # Independent unsaturated integers: these also check that adding the
    # two high/low allocations has ample room in the displayed 10^70 row.
    exactlow = pow10(embed)*drot*4*factorial(11)*(1+pow10(radius))*drot
    exacthigh = drot*4*factorial(14)*pow10(embed)*drot
    need("coefficient_H10_exact_four_factorial", 11*factorial(11)**2 <= (4*factorial(11))**2)
    need("coefficient_H13_exact_four_factorial", 14*factorial(14)**2 <= (4*factorial(14))**2)
    need("unsaturated_commutator_two_allocation_sum", subset*(exactlow+exacthigh) < 10**commper)
    commagg = claims["energy_commutator_aggregate_exponent"]
    need("energy_commutator_ball_square_sum_rounding", positive_sum_below(
        [(1,commper+sqrtwords),(1,commper+sqrtwords+radius)],commagg))
    commnormal = claims["normal_commutator_per_word_exponent"]
    need("normal_commutator_per_word_rounding", positive_sum_below(
        [(subset,lowcoef+dexp),(subset,highcoef+uinf),
         (subset,highcoef+uinf+radius)],commnormal))
    need("unsaturated_normal_commutator_two_allocation_sum",
         subset*(exactlow+exacthigh*(1+pow10(radius))) < 10**commnormal)
    product = claims["product_exponent"]
    need("common_product_interface", product >= cs+max(sagg,commagg,sper,commnormal))
    need("product_decimal_reserve", product-(cs+max(sagg,commagg,sper,commnormal)) >= 23)

    # Energy production from the full interior equation and zero flux rows.
    chain = claims["energy_derivative_coefficient"]
    div = claims["energy_divergence_coefficient"]
    need("energy_coefficient_first_chain", F(chain) >= 1+pow10(embed+radius))
    need("four_direction_divergence_contraction", div >= 2*chain)
    mixedproduction = claims["energy_mixed_production_coefficient"]
    production = claims["energy_production_coefficient"]
    need("energy_mixed_production_sqrt_rounding", (2*div)**2*2 < mixedproduction**2)
    need("energy_total_production_rounding", production >= 1+div+2*chain+mixedproduction)
    aggregate = claims["energy_aggregate_production_coefficient"]
    need("energy_square_root_production", F(production,2) <= aggregate)
    ku = claims["energy_unknown_exponent"]; kv = claims["energy_predecessor_exponent"]
    need("energy_unknown_final_rounding", positive_sum_below(
        [(1,product),(aggregate,cs+cy)],ku))
    need("energy_predecessor_final_rounding", product <= kv)

    # Interior triangular recovery: all coefficient differentiations occur
    # in the unisolated equation. Only the zeroth reciprocal is used.
    basecost = claims["normal_base_squared_cost"]
    base = claims["normal_base_factor"]
    need("normal_base_four_energy_groups", 2+4+2+104 <= basecost)
    need("normal_base_square_root_rounding", basecost < base**2)
    leading = claims["normal_leading_cost"]
    need("normal_interior_Az_and_double_BzA", 1+1+1+4+4 <= leading)
    need("angular_B_Hilbert_contraction", (2**2)*3 < 4**2)
    need("mixed_B_symmetry_factor", 2*2 <= 4)
    stages = claims["normal_stage_count"]
    need("all_thirteen_normal_stages", stages == m-1 == 13)
    totalrows = sum(comb(m+2-k,2)*rotations**k for k in range(m+1))
    rows = [sum((m+1-j-k)*rotations**k for k in range(m-j+1))
            for j in range(2,m+1)]
    need("normal_exact_full_row_count", totalrows == 135413275530)
    need("normal_exact_stage_row_counts", rows == [3134566561,522427758,
        87071291,14511880,2418645,403106,67183,11196,1865,310,51,8,1])
    rowsqrt = claims["normal_row_sqrt_exponent"]
    need("all_normal_row_square_sums", totalrows < 10**(2*rowsqrt)
         and all(n < 10**(2*rowsqrt) for n in rows))
    step = claims["normal_step_exponent"]
    lowstep = claims["normal_lower_step_exponent"]
    # Use 2 rather than the sharper 20/19, preserving the proof's ledger.
    need("normal_step_triangle_rounding", 1+2*leading*10**rowsqrt < 10**step)
    need("normal_lower_step_rounding", 2*10**rowsqrt < 10**lowstep)
    gsum = claims["normal_geometric_sum_exponent"]
    need("normal_geometric_series_last_power", gsum >= step*(stages-1))
    need("normal_geometric_series_factor_two", sum(F(1,10**(step*j))
         for j in range(stages)) < 2)
    wy = claims["normal_W_Y_exponent"]; wl = claims["normal_W_lower_exponent"]
    need("normal_thirteen_stage_Y_rounding", positive_sum_below([(base,step*stages)],wy))
    need("normal_thirteen_stage_lower_rounding", positive_sum_below([(2,lowstep+gsum)],wl))
    fullbridge = claims["mixed_norm_bridge_factor"]
    fullbexp = claims["mixed_norm_bridge_exponent"]
    need("full_spatial_derivative_slot_placements", 15*2**14 < fullbridge**2)
    need("literal_time_sum_bridge", 15*2**7 <= fullbridge)
    need("full_mixed_norm_decimal_rounding", fullbridge < 10**fullbexp)
    ny = claims["normal_Y_exponent"]; nu = claims["normal_unknown_exponent"]
    nv = claims["normal_predecessor_exponent"]
    need("normal_Y_final_rounding", wy+fullbexp <= ny)
    need("normal_unknown_final_rounding", wl+fullbexp+product <= nu)
    need("normal_predecessor_final_rounding", wl+fullbexp+product <= nv)

    # Accepted initial/seed interfaces. Their analytic theorems are inputs.
    cj = claims["C_J_exponent"]; lift = claims["C_lift_exponent"]
    seed = claims["C_seed_exponent"]; slope = claims["data_d_slope_exponent"]
    guard = claims["initial_jet_guard_exponent"]
    aseed = claims["seed_amplitude_exponent"]
    need("retained_completed_initial_jet_constant", cj == 6503000)
    need("retained_jet_table_lift_constant", lift == 255)
    need("pair_to_seed_constant", seed == cj+lift)
    need("retained_original_parameter_slope", slope == 4200)
    need("retained_initial_jet_guard", guard == -6503020)
    need("retained_seed_parameter_interval", aseed == -(cj+lift+slope+30))
    sd = claims["seed_data_exponent"]; si = claims["seed_I14_exponent"]
    sx = claims["seed_X14_exponent"]
    need("seed_interval_data_size", slope+aseed <= sd)
    need("seed_interval_initial_guard", sd < guard)
    need("seed_interval_initial_jets", cj+sd <= si)
    need("seed_interval_X14_size", lift+si <= sx)
    need("seed_interval_point_state", embed+sx <= -23)

    # Fundamental-theorem-of-calculus join, then scalar Volterra comparison.
    d0 = claims["D0_safe_exponent"]
    need("D0_three_initial_terms", positive_sum_below(
        [(1,ny+cy+cj),(1,nu+cj),(1,nv+cj)],d0))
    need("seed_below_D0_safe", seed < d0)
    c0 = claims["picard_C0_exponent"]; H = claims["picard_H_exponent"]
    need("safe_C0_strict_reserve", d0 < c0)
    need("Volterra_unknown_kernel", positive_sum_below([(1,ny+ku),(1,nu)],H))
    need("Volterra_predecessor_kernel", positive_sum_below([(1,ny+kv),(1,nv)],H))
    need("Volterra_supersolution_factor_two", 2*pow10(0) >= 1+1)

    # a_iter=10^(-10^P), exp(2HT)<10^(10^G): compare P and G
    # symbolically. No call to pow10(P), pow10(G), or integer 10**P occurs.
    approximation = claims["smooth_approximation_factor"]
    need("supported_smooth_approximation_allowance", approximation == 2)
    P = claims["iterated_amplitude_power_exponent"]
    G = claims["growth_power_exponent"]
    need("selected_nested_amplitude", P == 500209)
    comparison = claims["nested_seed_comparison_exponent"]
    need("nested_amplitude_below_seed_interval", P >= comparison == 7
         and 10**comparison > -aseed)
    need("smooth_approximation_initial_jet_guard", positive_sum_below(
        [(approximation,slope-10**comparison)],guard))
    growthordinary = H-6
    need("growth_power_rounding", positive_sum_below([(2*tfactor,growthordinary)],G))
    need("e_less_than_ten_finite_geometric_proof", F(1)+1+sum(F(1,2**(n-1))
         for n in range(2,20))+F(1,2**18) < 10)
    combined = claims["majorant_combined_linear_exponent"]
    need("majorant_linear_exponent_join", c0+slope <= combined)
    gap = claims["nested_gap_coefficient"]
    need("nested_decimal_power_difference", 0 <= P-G <= 100
         and gap == 10**(P-G)-1)
    small = claims["majorant_size_exponent"]
    need("selected_majorant_reporting_scale", small == -15)
    tinytest = claims["nested_gap_test_exponent"]
    need("nested_gap_finite_lower_comparison", G >= tinytest == 6
         and gap*10**tinytest > combined-small)
    need("majorant_strictly_inside_quarter_ball", positive_sum_below(
        [(4*approximation,small)],radius))

    return {"schema": "hdblast-variable-energy-arithmetic-v1", "status": "pass",
        "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(),
        "checks": records, "named_check_count": len(records),
        "exact": {"rotation_integer": str(drot), "energy_word_count": nwords,
            "source_leaf_count": sourceleaves, "normal_row_count": totalrows,
            "normal_stage_rows": rows, "six_output_exponents": [cy,ku,kv,ny,nu,nv],
            "temporal_seed_exponents": [cj,lift,seed],
            "Volterra_safe_exponents": [c0,H],
            "nested_amplitude": {"base": 10, "outer_sign": -1, "power_exponent": P},
            "growth_power_exponent": G, "nested_gap_coefficient": gap,
            "smooth_data_factor": approximation, "majorant_size_exponent": small},
        "scope": "Finite arithmetic for accepted conditional coefficient/product, frozen energy, "
            "interior normal recovery and compatible smooth-Picard admission proofs. Continuous "
            "analytic premises, finite-Sobolev convergence, top-norm continuity and evolved "
            "parameter/observer assertions are not proved by this runner."}


def controls(claims):
    records = []
    if len({field for field,value in MUTATIONS}) != len(MUTATIONS):
        raise ValueError("duplicate arithmetic mutation field")
    for field,value in MUTATIONS:
        bad = copy.deepcopy(claims); bad[field] = value
        payload = canonical(bad)
        try:
            evaluate(json.loads(payload), lambda name,ok: None)
        except ValueError as exc:
            records.append({"fault": field,
                "self_rehashed_sha256": hashlib.sha256(payload).hexdigest(),
                "rejected_by": "fresh exact variable-energy arithmetic",
                "diagnostic": str(exc)})
        else:
            raise ValueError("surviving variable-energy arithmetic fault: "+field)
    return {"schema": "hdblast-variable-energy-arithmetic-controls-v1",
            "status": "pass", "control_count": len(records), "controls": records}
