"""Exact finite ledger for nonlinear initial jets and their compatible seed.

The graph-domain, spectral trace, covariant product, support, and symmetry
statements are analytic premises proved in the companion texts. This module
checks their finite integer/rational coefficient and exponent allocations.
It does not certify a nonlinear evolved solution or a Picard sequence.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from math import comb, factorial, prod
from pathlib import Path
import sys


MUTATIONS = (
    ("q_ordinary_exponent", 17013),
    ("connection_exponent", 50),
    ("first_jet_substitution_exponent", 150120),
    ("a_reciprocal_zeroth_ceiling", 17000),
    ("a_nonconstant_aggregate_exponent", 18000),
    ("a_reciprocal_aggregate_exponent", 288277),
    ("a_reciprocal_ordinary_exponent", 288300),
    ("normalized_ordinary_exponent", 450200),
    ("normalized_aggregate_exponent", 450225),
    ("connection_aggregate_exponent", 60),
    ("affine_state_aggregate_ceiling", 15),
    ("OmegaU_aggregate_exponent", 72),
    ("horizontal_slot_ceiling", 39),
    ("horizontal_step_exponent", 70),
    ("after_horizontal_exponent", 451516),
    ("final_array_conversion_exponent", 19),
    ("converted_covariant_exponent", 451546),
    ("covariant_rounding_exponent", 451546),
    ("all_state_coefficient_exponent", 451000),
    ("term_count_exponent", 12),
    ("low_field_order", 6),
    ("embedding_rank_ceiling", 11),
    ("row_exponent", 500012),
    ("recursion_step_exponent", 500200),
    ("recursion_steps", 12),
    ("sharp_jet_map_exponent", 6502612),
    ("C_J_exponent", 6502612),
    ("jet_data_guard_exponent", -6503000),
    ("a_jets_exponent", -6507219),
    ("graph_exponent", 23),
    ("Gaussian_exponent", 48),
    ("bump_exponent", 52),
    ("bump_integral_lower", "1/5"),
    ("cutoff_jet_exponent", 149),
    ("cutoff_multiplier_exponent", 154),
    ("cutoff_plateau_left", "5/8"),
    ("cutoff_plateau_right", "19/8"),
    ("lift_exponent", 251),
    ("seed_from_data_exponent", 6503000),
    ("a_seed_exponent", -6507484),
    ("seed_data_exponent", -6503286),
    ("seed_jet_exponent", -286),
    ("seed_size_exponent", -31),
    ("seed_point_exponent", -24),
    ("Gamma_averaging_weight", "1/17"),
    ("Gamma_projection_norm", "1/2"),
    ("seed_union_norm_factor", "1/2"),
    ("cutoff_after_heat", False),
    ("projection_after_lift", False),
    ("same_original_parameter", False),
    ("retained_nonzero_first_velocity", False),
    ("smooth_approximants_recompute_nonlinear_jets", False),
    ("initial_curve_zero_at_parameter_zero", False),
    ("input_jet_table_Gamma_invariant", False),
    ("Gamma_natural_full_slot_action", False),
    ("Gamma_preserves_boundary_split", False),
    ("initial_pair_fixed_zero_collars", False),
)
FACTORS = (
    'A_spectral_shift', 'C_J_exponent', 'DD_first_frequency_index',
    'D_components', 'Frobenius_star_factor', 'Gamma_averaging_weight',
    'Gamma_natural_full_slot_action', 'Gamma_order', 'Gamma_preserves_boundary_split',
    'Gamma_projection_norm', 'Gamma_reflection_count', 'Gaussian_exponent',
    'Gaussian_integer', 'I14_is_global_union_norm', 'K_sum_factor',
    'ND_frequency_offset', 'N_components', 'OmegaU_aggregate_exponent',
    'a_jets_exponent', 'a_nonconstant_aggregate_exponent', 'a_reciprocal_aggregate_exponent',
    'a_reciprocal_ordinary_exponent', 'a_reciprocal_zeroth_ceiling', 'a_seed_exponent',
    'affine_state_aggregate_ceiling', 'after_horizontal_exponent', 'all_state_coefficient_exponent',
    'all_state_horizontal_conversion', 'base_dimension', 'bounded_linear_lift',
    'bump_derivative_integer', 'bump_exponent', 'bump_integral_lower',
    'coefficient_chain_branches', 'coefficient_state_patch_exponent', 'coefficient_variables',
    'connection_aggregate_exponent', 'connection_exponent', 'converted_covariant_exponent',
    'covariant_rounding_exponent', 'cutoff_after_heat', 'cutoff_jet_exponent',
    'cutoff_left_shift', 'cutoff_multiplier_exponent', 'cutoff_plateau_left',
    'cutoff_plateau_right', 'cutoff_right_shift', 'cutoff_scale_factor',
    'cutoff_support_left', 'cutoff_support_right', 'cylinder_length_in_ell',
    'ell', 'embedding_exponent', 'embedding_order',
    'embedding_rank_ceiling', 'exact_zero_state_source', 'external_derivative_order',
    'final_array_conversion_exponent', 'first_jet_substitution_exponent', 'full_angular_spectrum',
    'graph_bridge_integer', 'graph_exponent', 'graph_word_stage_factor',
    'heat_polynomial_coefficient_integer', 'hierarchy_size_at_jet_guard_exponent', 'horizontal_slot_ceiling',
    'horizontal_step_exponent', 'horizontal_steps', 'initial_curve_zero_at_parameter_zero',
    'initial_field_factor_count', 'initial_first_derivative_exponent', 'initial_pair_fixed_zero_collars',
    'input_jet_table_Gamma_invariant', 'jet_data_guard_exponent', 'jet_support_left',
    'jet_support_right', 'lift_exponent', 'low_field_order',
    'max_order', 'maximum_field_factors', 'maximum_intrinsic_order',
    'normalized_aggregate_exponent', 'normalized_groups', 'normalized_ordinary_exponent',
    'old_initial_amplitude_exponent', 'old_initial_size_exponent', 'ordinary_coefficient_order',
    'original_rank_ceiling', 'orthonormal_Frobenius_components', 'partial_sum_guard_exponent',
    'projection_after_lift', 'q_aggregate_exponent', 'q_ordinary_exponent',
    'raw_lower_exponent', 'recursion_step_exponent', 'recursion_steps',
    'retained_copy_count', 'retained_nonzero_first_velocity', 'review_status',
    'row_exponent', 'same_original_parameter', 'seed_data_exponent',
    'seed_from_data_exponent', 'seed_jet_exponent', 'seed_point_exponent',
    'seed_size_exponent', 'seed_union_norm_factor', 'sharp_jet_map_exponent',
    'smooth_approximants_recompute_nonlinear_jets', 'state_dimension', 'term_count_exponent',
    'term_count_integer', 'time_derivative_count', 'zero_collar_graph_domain',
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def pow10(exponent):
    if type(exponent) is not int or abs(exponent) > 1000:
        raise ValueError("large decimal expansion forbidden")
    return F(10**exponent) if exponent >= 0 else F(1, 10**(-exponent))


def positive_sum_below(terms, exponent, strict=True):
    """Conservatively round relative exponents below -1000 upward.

    This proves passing comparisons without expanding any million-decade
    power. Each exact coefficient is nonnegative and the rounded sum is an
    upper bound. Excessive positive relative powers reject immediately.
    """
    if type(exponent) is not int:
        raise ValueError("integer target exponent required")
    total = F(0)
    for coefficient, power in terms:
        coefficient = F(coefficient)
        if coefficient < 0 or type(power) is not int:
            raise ValueError("nonnegative exact monomial required")
        relative = power - exponent
        if relative > 1000:
            return False
        total += coefficient * pow10(max(relative, -1000))
    return total < 1 if strict else total <= 1


def evaluate(claims, check):
    records = []

    def need(name, condition):
        ok = bool(condition)
        check("initialization_arithmetic/" + name, ok)
        records.append({"name": name, "status": "pass" if ok else "fail"})
        if not ok:
            raise ValueError("initialization arithmetic rejects " + name)

    need("schema", claims["schema"] == "hdblast-temporal-seed-arithmetic-claims-v1"
         and set(claims) == {"schema", *FACTORS})
    need("review_status", claims["review_status"] == "accepted_conditional_analytic_inputs")
    m = claims["max_order"]; p = claims["ordinary_coefficient_order"]
    need("actual_orders", m == 14 and p == 16)
    state = claims["state_dimension"]; base = claims["base_dimension"]
    variables = claims["coefficient_variables"]
    need("actual_dimensions", state == 16 and base == 5 and variables == state + base)
    ell = F(claims["ell"])
    need("actual_cylinder_scale", ell == F(1, 10**6) and claims["cylinder_length_in_ell"] == 6)
    need("actual_source_guards", claims["same_original_parameter"] is True
         and claims["retained_nonzero_first_velocity"] is True
         and claims["exact_zero_state_source"] is True
         and claims["all_state_horizontal_conversion"] is True
         and claims["initial_pair_fixed_zero_collars"] is True)

    # All-state ordinary coefficients: use the sharper principal q branch.
    qagg = claims["q_aggregate_exponent"]; qord = claims["q_ordinary_exponent"]
    raw = claims["raw_lower_exponent"]; conn = claims["connection_exponent"]
    need("retained_analytic_coefficient_inputs", qagg == 17000 and raw == 150000)
    need("q_Taylor_to_ordinary", positive_sum_below([(factorial(p), qagg)], qord))
    connection_integer = 2*10**4*factorial(17)*100**17
    need("whole_connection_Cauchy_bound", positive_sum_below([(connection_integer, 0)], conn))
    substitution = claims["first_jet_substitution_exponent"]
    need("full_first_jet_substitution", positive_sum_below(
        [(100*5**p, raw+2*conn), (100*5**p, qord+2*conn)], substitution))
    a0 = claims["a_reciprocal_zeroth_ceiling"]
    need("all_state_lapse_reciprocal_input", a0 >= 18000 and a0 < 10**5)
    count = comb(variables+p, p)
    need("twenty_one_variable_multiindices", count == 12875774670 and count < 10**12)
    an = claims["a_nonconstant_aggregate_exponent"]
    need("nonconstant_a_Taylor_aggregate", positive_sum_below([(count, qord)], an))
    ar = claims["a_reciprocal_aggregate_exponent"]
    need("reciprocal_finite_geometric_sum", an+5 >= 0 and positive_sum_below(
        [(2**p, 5+p*(an+5))], ar))
    ari = claims["a_reciprocal_ordinary_exponent"]
    need("reciprocal_restore_ordinary_derivatives", positive_sum_below([(factorial(p), ar)], ari))
    ordinary = claims["normalized_ordinary_exponent"]
    need("five_normalized_groups", positive_sum_below(
        [(100*2**p, ari+substitution), (100*2**p, ari+qord)], ordinary))

    # Horizontal conversion acts on the state and every tensor derivative slot.
    aggregate = claims["normalized_aggregate_exponent"]
    ca = claims["connection_aggregate_exponent"]
    need("ordinary_to_coefficient_aggregate", positive_sum_below([(count, ordinary)], aggregate))
    need("connection_Taylor_aggregate", positive_sum_below([(count, conn)], ca))
    affine = F(claims["affine_state_aggregate_ceiling"])
    need("affine_state_coordinate_jet", affine >= state + F(state, 1000))
    omega = claims["OmegaU_aggregate_exponent"]
    need("state_transport_coefficient", positive_sum_below([(affine, ca)], omega))
    slots = claims["horizontal_slot_ceiling"]
    need("every_input_and_derivative_slot", slots >= 8+2*p)
    step = claims["horizontal_step_exponent"]
    need("full_horizontal_array_operation", positive_sum_below(
        [(base*p, 0), (base*p, omega), (base*slots, ca)], step))
    need("full_vertical_array_operation", positive_sum_below([(state*p, 0)], step))
    steps = claims["horizontal_steps"]
    need("horizontal_derivative_order", steps == p)
    after = claims["after_horizontal_exponent"]
    need("all_ordered_conversion_words", aggregate+steps*step <= after)
    frame = claims["final_array_conversion_exponent"]
    need("full_array_orthonormal_conversion", positive_sum_below([(2**64, 0)], frame))
    converted = claims["converted_covariant_exponent"]
    need("conversion_final_exponent", after+frame <= converted)
    rounded = claims["covariant_rounding_exponent"]
    Cstar = claims["all_state_coefficient_exponent"]
    need("covariant_coefficient_rounding", converted < rounded <= Cstar)

    # Low/high products and thirteen triangular rows. The source zero identity
    # ensures every row term has an actual field factor (an analytic input).
    groups = claims["normalized_groups"]; initial_fields = claims["initial_field_factor_count"]
    branches = claims["coefficient_chain_branches"]; derivatives = claims["external_derivative_order"]
    need("finite_row_structure", groups == 5 and initial_fields == 2 and branches == 2
         and derivatives == m-2)
    leaves = groups*prod(k+initial_fields+branches for k in range(derivatives))
    need("ordered_leaf_count", leaves == 5*factorial(15)//factorial(3)
         == int(claims["term_count_integer"]))
    need("term_count_rounding", positive_sum_below([(leaves, 0)], claims["term_count_exponent"]))
    intrinsic = claims["maximum_intrinsic_order"]
    low = claims["low_field_order"]; embed_order = claims["embedding_order"]
    need("highest_field_factor_order", intrinsic == m and claims["maximum_field_factors"] == 14)
    need("at_most_one_unembedded_high_factor", 2*(low+1) > intrinsic and low >= 7)
    need("spatial_embedding_reserve", embed_order == 3 and low+embed_order <= m)
    rank = claims["original_rank_ceiling"]; rank_ceiling = claims["embedding_rank_ceiling"]
    need("full_tensor_embedding_rank", rank == 2 and rank+low+embed_order <= rank_ceiling == 16)
    embed = claims["embedding_exponent"]; partial = claims["partial_sum_guard_exponent"]
    patch = claims["coefficient_state_patch_exponent"]
    need("retained_embedding_and_patch", embed == 7 and patch == -3)
    need("small_factors_in_Linfinity", embed+partial < 0)
    need("coefficient_state_admission", positive_sum_below(
        [(F(claims["Frobenius_star_factor"]), embed+partial)], patch))
    row = claims["row_exponent"]
    need("square_summed_output_row", positive_sum_below(
        [(13*leaves*leaves, 2*Cstar)], 2*row))
    growth = claims["recursion_step_exponent"]
    need("one_plus_row_growth", positive_sum_below([(1, 0), (1, row)], growth))
    nsteps = claims["recursion_steps"]
    need("number_of_triangular_rows", nsteps == m-1 == 13)
    sharp = claims["sharp_jet_map_exponent"]; CJ = claims["C_J_exponent"]
    need("thirteen_row_exact_exponent", nsteps*growth <= sharp)
    need("jet_map_rounding", sharp <= CJ)
    guard = claims["jet_data_guard_exponent"]
    guarded = claims["hierarchy_size_at_jet_guard_exponent"]
    need("noncircular_hierarchy_guard", CJ+guard <= guarded < partial)
    need("relative_open_pair_ball", positive_sum_below([(2, CJ+guard)], partial))
    first = claims["initial_first_derivative_exponent"]
    olda = claims["old_initial_amplitude_exponent"]
    need("same_completed_initial_curve", first == 4200 and olda == -6402
         and claims["old_initial_size_exponent"] == -2200
         and claims["initial_curve_zero_at_parameter_zero"] is True)
    ajets = claims["a_jets_exponent"]
    need("jet_parameter_interval", first+ajets <= guard and ajets < olda)

    # Auxiliary A uses the full scalar spherical tower in the quaternion frame.
    need("selfadjoint_component_boundary_split", claims["N_components"] == 12
         and claims["D_components"] == 4 and claims["N_components"]+claims["D_components"] == state
         and F(claims["ND_frequency_offset"]) == F(1, 2)
         and claims["DD_first_frequency_index"] == 1
         and claims["A_spectral_shift"] == 1)
    need("graph_domain_scope", claims["zero_collar_graph_domain"] is True
         and claims["orthonormal_Frobenius_components"] is True
         and claims["full_angular_spectrum"] is True)
    stage = claims["graph_word_stage_factor"]
    need("all_covariant_derivative_slot_stages", all(
        1+2*(r+j) <= stage*(j+1) for r in range(rank+1) for j in range(m)))
    need("word_to_graph_binomials", all(1 <= comb(n,j) <= 2**n
         for n in range(m+1) for j in range(n+1)))
    need("graph_sum_and_binomial_factors", m+1 < 4**2 and 2**(m//2)*4 == 512)
    CB = 512*stage**m*factorial(m)
    need("graph_bridge_exact_integer", CB == int(claims["graph_bridge_integer"]))
    graph = claims["graph_exponent"]
    need("two_direction_graph_rounding", positive_sum_below([(CB, 0)], graph))

    # Gaussian multiplier bound, finite counterterms, and exact spectral powers.
    polynomial = factorial(29)//factorial(15)
    need("Hermite_coefficient_product", max(prod(j+i+2 for i in range(k))
         for j in range(m+1) for k in range(m+1)) == polynomial
         == int(claims["heat_polynomial_coefficient_integer"]))
    need("Hermite_degree_ceiling", max(j+k for j in range(m+1) for k in range(m+1)) == 28)
    CH = 2*factorial(28)*polynomial
    need("Gaussian_exact_integer", CH == int(claims["Gaussian_integer"]))
    gaussian = claims["Gaussian_exponent"]
    need("Gaussian_rounding", positive_sum_below([(CH, 0)], gaussian))
    need("all_spectral_weight_balances", all((m-k)+(k-j)+2*r == m-(j-2*r)
         for k in range(m+1) for j in range(m+1) for r in range(j//2+1)))
    need("finite_heat_counterterm_cancellation", all(sum(
        (F((-1)**r, factorial(r)*factorial(h-r)) for r in range(h+1)), F(0))
        == (1 if h == 0 else 0) for h in range(m//2+1)))
    # An exact upper tail proves e<3 and hence all finite counterterm sums.
    e_upper = sum((F(1, factorial(r)) for r in range(9)), F(0)) + F(10, 9*factorial(9))
    Kfactor = claims["K_sum_factor"]; times = claims["time_derivative_count"]
    need("counterterm_norm_sum", e_upper < Kfactor == 3)
    need("all_time_derivatives_counted", times == m+1 == 15)

    # Flat bump and complete H14 multiplier, not just a zeroth flat weight.
    exp_half_upper = sum((F(1, 2**r*factorial(r)) for r in range(17)), F(0))
    exp_half_upper += F(1, 2**17*factorial(17))*F(36, 35)
    integral = F(claims["bump_integral_lower"])
    need("bump_normalization_lower_bound", exp_half_upper < 2 and integral <= F(1, 6))
    need("Bell_composition_coefficient_bound", all(sum(
        (F(comb(n-1,r-1), factorial(r)) for r in range(1,n+1)), F(0)) <= 2**n
        for n in range(1,m)))
    need("bump_exponential_prefactor", 3**4 < 100)
    BQ = 100*factorial(13)*2**13*26**26
    need("bump_derivative_majorant", max(100*factorial(n)*2**n*(2*n)**(2*n)
         for n in range(1,m)) == BQ == int(claims["bump_derivative_integer"]))
    need("bump_derivative_rounding", positive_sum_below([(BQ, 0)], claims["bump_exponent"]))
    scale = claims["cutoff_scale_factor"]
    need("actual_cutoff_affine_arguments", scale == 8
         and claims["cutoff_left_shift"] == 2 and claims["cutoff_right_shift"] == 22)
    support_l = F(claims["cutoff_support_left"]); support_r = F(claims["cutoff_support_right"])
    plateau_l = F(claims["cutoff_plateau_left"]); plateau_r = F(claims["cutoff_plateau_right"])
    need("cutoff_support_and_plateau", support_l == F(2,scale) and support_r == F(22,scale)
         and plateau_l == F(3,scale) and plateau_r == F(21,scale)
         and 0 < support_l < plateau_l < plateau_r < support_r < 6)
    need("original_jet_support_inside_plateau", plateau_l < F(claims["jet_support_left"])
         == F(1,2) < F(claims["jet_support_right"]) == F(5,2) < plateau_r)
    need("smooth_support_enlargement_inside_plateau", plateau_l < F(7,16) < F(1,2)
         and F(5,2) < F(41,16) < plateau_r)
    need("disjoint_cutoff_transitions", plateau_l < plateau_r and scale/ell > 1)
    Mchi = (1/integral)*BQ*(scale/ell)**m
    Cchi = 2**m*4*Mchi
    need("complete_cutoff_jet_bound", positive_sum_below([(Mchi, 0)], claims["cutoff_jet_exponent"]))
    need("complete_H14_cutoff_multiplier", positive_sum_below([(Cchi, 0)], claims["cutoff_multiplier_exponent"]))

    # These guard checks encode the analytic construction being replayed; the
    # norm-one isometries and density theorem themselves are analytic premises.
    need("Gamma_average_identity", claims["Gamma_reflection_count"] == 4
         and claims["Gamma_order"] == 2**claims["Gamma_reflection_count"] == 16
         and F(claims["Gamma_averaging_weight"])*claims["Gamma_order"] == 1
         and claims["input_jet_table_Gamma_invariant"] is True
         and claims["Gamma_natural_full_slot_action"] is True
         and claims["Gamma_preserves_boundary_split"] is True)
    need("Gamma_projection_norm_charge", F(claims["Gamma_projection_norm"]) >= 1)
    need("after_lift_cutoff_and_projection", claims["cutoff_after_heat"] is True
         and claims["projection_after_lift"] is True)
    need("compatible_smooth_approximation_guard", claims["smooth_approximants_recompute_nonlinear_jets"] is True)
    need("retained_union_input_norm", claims["retained_copy_count"] == 2
         and claims["I14_is_global_union_norm"] is True and F(claims["seed_union_norm_factor"]) >= 1)
    lift = claims["lift_exponent"]
    exact_lift = times*Kfactor*CB**2*CH*Cchi*F(claims["Gamma_projection_norm"])*F(claims["seed_union_norm_factor"])
    need("complete_lift_exact_product", positive_sum_below([(exact_lift, 0)], lift))
    need("lift_from_displayed_decimal_ceilings", positive_sum_below(
        [(times*Kfactor*F(claims["Gamma_projection_norm"])*F(claims["seed_union_norm_factor"]),
          2*graph+gaussian+claims["cutoff_multiplier_exponent"])], lift))
    joined = claims["seed_from_data_exponent"]
    need("distinct_jet_and_data_input_spaces", CJ+lift <= joined and joined > lift)
    aseed = claims["a_seed_exponent"]; dseed = claims["seed_data_exponent"]
    iseed = claims["seed_jet_exponent"]; seed = claims["seed_size_exponent"]
    point = claims["seed_point_exponent"]
    need("selected_seed_interval_join", first+aseed <= dseed and aseed < ajets < olda)
    need("seed_inside_jet_data_guard", dseed < guard)
    need("seed_initial_hierarchy_join", CJ+dseed <= iseed)
    need("seed_lift_smallness_join", lift+iseed <= seed <= -30)
    need("direct_seed_from_data_join", joined+dseed <= seed)
    need("seed_point_state_guard", embed+seed <= point < patch)
    need("relative_open_C2_admission", positive_sum_below([(1,dseed)], guard)
         and claims["bounded_linear_lift"] is True)

    return {"schema": "hdblast-temporal-seed-arithmetic-v1", "status": "pass",
        "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(),
        "checks": records, "named_check_count": len(records),
        "exact": {"multiindex_count": count, "connection_integer": str(connection_integer),
                  "ordered_leaf_count": str(leaves), "sharp_jet_map_exponent": sharp,
                  "graph_integer": str(CB), "Gaussian_integer": str(CH),
                  "bump_integer": str(BQ), "cutoff_jet_integer": str(Mchi.numerator),
                  "cutoff_multiplier_integer": str(Cchi.numerator),
                  "complete_lift_integer": str(exact_lift.numerator),
                  "lift_exponent": lift, "seed_from_data_exponent": joined,
                  "seed_admission_exponents": [aseed,dseed,iseed,seed,point]},
        "scope": "Exact finite allocations conditional on the accepted full-state coefficient, "
                 "initial-jet and compatible-seed analytic arguments. Spectral domains, strong "
                 "traces, covariant products, symmetry isometries and density are not inferred "
                 "from arithmetic. No evolved solution or Picard-sequence patch is certified."}


def controls(claims):
    records = []
    for field, value in MUTATIONS:
        bad = copy.deepcopy(claims); bad[field] = value
        payload = canonical(bad)
        try:
            evaluate(json.loads(payload), lambda name, ok: None)
        except ValueError as exc:
            records.append({"fault": field, "self_rehashed_sha256": hashlib.sha256(payload).hexdigest(),
                            "rejected_by": "fresh finite initialization ledger", "diagnostic": str(exc)})
        else:
            raise ValueError("surviving initialization arithmetic fault: " + field)
    return {"schema": "hdblast-temporal-seed-arithmetic-controls-v1", "status": "pass",
            "control_count": len(records), "controls": records}
