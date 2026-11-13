"""Exact finite ledger for the uniform initial Cauchy map.

Continuous graded products, Nemytskii estimates, support and existence are
written analytic arguments. This replay checks their explicit finite
coefficient, branch-count, scale and parameter-chain allocations.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
from math import factorial, prod
from pathlib import Path
import sys


MUTATIONS = (
    ("initial_amplitude_exponent", -6300),
    ("initial_K_exponent", 6300),
    ("potential_radius_exponent", -6405),
    ("potential_second_derivative_factor", 4),
    ("combined_correction_exponent", 29),
    ("five_field_state_exponent", -6202),
    ("five_field_first_exponent", 200),
    ("five_field_second_exponent", 6431),
    ("rho_scalar_jet_exponent", 27),
    ("rho_oneform_jet_exponent", 26),
    ("rho_multiplier_exponent", 32),
    ("halfdet_first_factor", 3),
    ("halfdet_second_factor", 15),
    ("state_monomial_ceiling", 31),
    ("first_monomial_ceiling", 767),
    ("second_monomial_ceiling", 10751),
    ("coefficient_contraction_ceiling", 8),
    ("graded_state_exponent", 37),
    ("graded_first_exponent", 38),
    ("graded_second_exponent", 139),
    ("common_graded_state_first_exponent", 38),
    ("physical_multiplier_exponent", 200),
    ("raw_state_exponent", -5936),
    ("raw_first_exponent", 466),
    ("raw_second_exponent", 6698),
    ("raw_to_product_exponent", 3215),
    ("flow_total_order", 16),
    ("flow_Hilbert_jet_exponent", 126),
    ("root_jet_exponent", 14),
    ("Bell_exponent", 22),
    ("point_chart_jet_exponent", 419),
    ("functional_branch_exponent", 16),
    ("functional_D1_exponent", 436),
    ("functional_D2_exponent", 443),
    ("functional_D3_exponent", 450),
    ("raw_chart_admission_exponent", -3227),
    ("chart_linear_exponent", 3660),
    ("chart_quadratic_exponent", 6883),
    ("copy_linear_factor", 1),
    ("copy_quadratic_factor", 2),
    ("completed_Hessian_exponent", 10361),
    ("direct_initial_size_exponent", -2282),
    ("completed_first_derivative_exponent", 4128),
    ("Taylor_initial_size_exponent", -1202),
)


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def pow10(exponent):
    if type(exponent) is not int or abs(exponent) > 1000:
        raise ValueError("large decimal expansion forbidden")
    return F(10**exponent) if exponent >= 0 else F(1, 10**(-exponent))


def positive_sum_below(terms, exponent, strict=True):
    """Use a rational upper sum, rounding relative exponents below -100 up.

    Thus no integer/denominator with thousands of digits is constructed.
    A passing comparison is a conservative proof of the actual comparison.
    """
    total = F(0)
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
        check("cauchy_arithmetic/" + name, boolean)
        records.append({"name": name, "status": "pass" if boolean else "fail"})
        if not boolean:
            raise ValueError("Cauchy arithmetic rejects " + name)

    need("schema", claims["schema"] == "hdblast-uniform-cauchy-arithmetic-claims-v1")
    need("review_status", claims["review_status"] == "accepted_conditional_analytic_inputs")

    # Exact retained scales and same selected initial interval.
    ell = F(claims["ell"]); kappa = F(claims["weight_kappa"])
    need("selected_weight_scales", ell == F(1,10**6) and kappa == ell/100)
    aexp = claims["initial_amplitude_exponent"]
    Kexp = claims["initial_K_exponent"]
    rexp = claims["potential_radius_exponent"]
    need("retained_initial_parameter_scale", aexp == -6402 and Kexp == 6400)
    need("retained_potential_radius", rexp == Kexp+2*aexp == -6404)
    uprime = F(claims["potential_first_derivative_bound"])
    usecond = F(claims["potential_second_derivative_factor"])
    need("retained_first_parameter_bound", uprime == F(101,5000))
    need("retained_second_parameter_bound", usecond >= 8)

    # Rebuild the all-point polynomial-vector-field coefficient reserve.
    U = [F(v) for v in claims["potential_coefficients"]]
    expected_U = [F(-2,27),F(0),F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)]
    need("actual_potential", U == expected_U)
    Uprime_l1 = sum((abs(i*U[i]) for i in range(1,len(U))), F(0))
    need("potential_derivative_l1", Uprime_l1 == 12)
    vector_norms = [4*kappa**2/ell, 4*kappa**2,
                    4*kappa**2*F(4,3), 4*kappa**2*(Uprime_l1+4), 4*kappa**2]
    vector_ceiling = F(claims["vectorfield_coefficient_ceiling"])
    degree = claims["vectorfield_degree_ceiling"]
    need("actual_polynomial_vectorfield_ceiling", max(vector_norms) <= vector_ceiling == F(1,10**8))
    need("actual_polynomial_vectorfield_degree", degree == 9)
    need("polynomial_D_degree_increment", degree-1 == claims["D_degree_increment"] == 8)
    need("primitive_weight_jets_through18", all(vector_ceiling*(10+8*i) < 1 for i in range(18)))
    need("primitive_root_chart_weight_coefficients", F(6000006,10**16) < 1
         and F(8,10**10) < 1 and F(4,10**16) < 1)

    # Five-field canonical/lapse envelope; no output Cauchy Hessian is used.
    CH = 8*17**14*factorial(14); CM = 8*17**13*factorial(13)
    correction_exp = claims["combined_correction_exponent"]
    need("accepted_correction_sum", CH+CM < 10**correction_exp)
    Bexp = claims["seed_exponent"]
    need("accepted_seed_bound", Bexp == 200)
    delta0 = claims["five_field_state_exponent"]
    delta1 = claims["five_field_first_exponent"]
    delta2 = claims["five_field_second_exponent"]
    need("five_field_state_from_seed_and_curve", positive_sum_below(
        [(2,Bexp+aexp),(1,correction_exp+Kexp+2*aexp)], delta0))
    need("five_field_first_from_seed_and_curve", positive_sum_below(
        [(2,Bexp),(uprime,correction_exp)], delta1))
    need("five_field_second_from_curve", positive_sum_below(
        [(usecond,correction_exp+Kexp)], delta2))
    Qexp = claims["graded_algebra_exponent"]
    need("same_graded_algebra", Qexp == 100)
    guard = F(claims["graded_ball_denominator"])
    need("whole_state_in_composition_ball", guard >= 1024
         and positive_sum_below([(guard,Qexp+delta0)],0))

    # Rho/background coefficient tensors, including derivative slots.
    need("rho_scalar_covariant_stage", all(k+5*j <= 5*(j+1)
         for k in range(4) for j in range(14)))
    need("rho_oneform_covariant_stage", all(5*j+6 <= 6*(j+1) for j in range(14)))
    rho_scalar = 131**3*5**14*factorial(14)
    rho_oneform = 131**2*6**14*factorial(14)
    need("rho_scalar_exact_integer", rho_scalar == 1196195873059687500000000000)
    need("rho_oneform_exact_integer", rho_oneform == 117238012873166550151987200)
    rs = claims["rho_scalar_jet_exponent"]; rv = claims["rho_oneform_jet_exponent"]
    need("rho_scalar_jet_rounding", rho_scalar < 10**rs)
    need("rho_oneform_jet_rounding", rho_oneform < 10**rv)
    multiplier = claims["rho_multiplier_exponent"]
    need("rho_multiplier_Leibniz_sum", (2**15-1)*10**max(rs,rv) < 10**multiplier)

    # Nonlinear inverse determinant, not the scalar or zero-state formula.
    need("positive_halfdet_branch", F(4,1024) < F(1,2))
    need("halfdet_first_direction", F(1,2)*2*4 <= claims["halfdet_first_factor"])
    need("halfdet_second_direction", 2*(F(1,4)*4*4+F(1,2)*8)
         <= claims["halfdet_second_factor"])
    state_max = claims["state_monomial_ceiling"]
    first_max = claims["first_monomial_ceiling"]
    second_max = claims["second_monomial_ceiling"]
    need("finite_state_product_rule", max(2**(n-1) for n in range(1,7)) <= state_max)
    need("finite_first_product_rule", max(n*4*2**(n-1) for n in range(1,7)) <= first_max)
    need("finite_second_product_rule", max(n*16*2**(n-1)+n*(n-1)*16*F(2)**(n-2)
         for n in range(1,7)) <= second_max)
    count = claims["monomial_count_ceiling"]; coefficient = claims["coefficient_contraction_ceiling"]
    need("accepted_finite_row_reserves", count == 32 and coefficient >= 4*2*2)
    need("free_trace_identity_insertion", F(4,3)*2*2 < coefficient)
    totals = (count*coefficient*state_max, count*coefficient*first_max, count*coefficient*second_max)
    need("finite_state_total", totals[0] == 16384)
    need("finite_first_total", totals[1] == 393216)
    need("finite_second_total", totals[2] == 5505024)
    gs = claims["graded_state_exponent"]; gf = claims["graded_first_exponent"]
    g2 = claims["graded_second_exponent"]
    need("graded_state_rounding", positive_sum_below([(totals[0],multiplier)],gs))
    need("graded_first_rounding", positive_sum_below([(totals[1],multiplier)],gf))
    need("graded_second_rounding", positive_sum_below([(totals[2],multiplier+Qexp)],g2))
    common = claims["common_graded_state_first_exponent"]
    need("common_graded_rounding", common >= max(gs,gf))

    # Strip exactly one already differentiated flat multiplier.
    physical = claims["physical_multiplier_exponent"]
    need("flat_multiplier_reconstruction", 24*8 == 192 and 24**24 < 10**34 and physical >= 192+34)
    raw_linear = claims["raw_linear_map_exponent"]; raw_quadratic = claims["raw_quadratic_map_exponent"]
    need("raw_linear_conversion", raw_linear >= physical+common)
    need("raw_quadratic_conversion", raw_quadratic >= physical+g2)
    raw0 = claims["raw_state_exponent"]; raw1 = claims["raw_first_exponent"]
    raw2 = claims["raw_second_exponent"]
    need("direct_raw_state", raw_linear+delta0 <= raw0)
    need("raw_first_chain", raw_linear+delta1 <= raw1)
    need("raw_second_chain", positive_sum_below(
        [(1,raw_linear+delta2),(1,raw_quadratic+2*delta1)],raw2))

    # The full initial boundary chart at a nonzero state.
    Iexp = claims["raw_frame_exponent"]; Lexp = claims["normalization_jet_exponent"]
    mult_exp = claims["linear_multiplier_exponent"]; Aexp = claims["raw_to_product_exponent"]
    need("retained_generic_linear_inputs", Iexp == 3000 and Lexp == 210)
    need("linear_multiplier_factor", 2**28*15 < 10**(2*mult_exp))
    need("raw_to_product_join", Iexp+Lexp+mult_exp <= Aexp)
    flow_order = claims["flow_total_order"]
    need("flow_derivative_budget", flow_order == 18
         and max(m+k for m,k in [(14,1),(14,2),(14,3),(13,4)]) <= flow_order)
    fs = claims["flow_star_jet_exponent"]; fe = claims["flow_Hilbert_jet_exponent"]
    need("retained_flow_jet_input", fs == 126)
    need("flow_Hilbert_norm_conversion", positive_sum_below([(3,fs)],fe))
    root_radius = F(claims["root_state_radius"])
    need("proper_y_root_branch", 0 < root_radius <= F(1,1000)
         and 4*(1-root_radius) >= 1)
    need("root_high_derivative_guard", 1/(1-18*root_radius) < 2)
    root_exp = claims["root_jet_exponent"]
    need("root_high_derivative_factorial", factorial(17) < 10**root_exp)
    Bell_exp = claims["Bell_exponent"]
    need("Bell_partition_ceiling", 18**18 < 10**Bell_exp)
    point = claims["point_chart_jet_exponent"]
    need("whole_state_point_chart_jets", fe+Bell_exp+18*root_exp <= point)
    branch = claims["functional_branch_exponent"]
    need("functional_Faa_di_Bruno_branch_count", all(
        prod(k+j+2 for j in range(m)) <= factorial(18)
        for m,k in [(13,1),(13,2),(13,3),(14,1),(14,2),(14,3),(13,4)]))
    need("functional_derivative_square_sum", 4*factorial(18) < 10**branch)
    need("all_field_low_high_derivative_split", all(2*(m-2)>m for m in (13,14)))
    need("embedded_field_derivative_rank", 2+11 <= 16)
    embedding = claims["cylinder_embedding_exponent"]
    need("retained_cylinder_embedding", embedding == 7)
    functional = [claims["functional_D1_exponent"], claims["functional_D2_exponent"],
                  claims["functional_D3_exponent"]]
    for k,value in enumerate(functional,1):
        need("functional_D"+str(k)+"_rounding", point+branch+embedding*(k-1) <= value)
    need("functional_C2_continuity_margin", point+branch+3*embedding <= claims["functional_D4_H13_exponent"])
    base = claims["product_base_exponent"]; admission = claims["raw_chart_admission_exponent"]
    need("uniform_raw_chart_ball", Aexp+admission <= base)
    need("direct_curve_chart_admission", raw0 <= admission)
    point_exp = claims["unrooted_point_exponent"]
    need("point_state_embedding", embedding+base <= point_exp)
    need("small_base_for_all_products", embedding+base < 0)
    need("root_domain_contains_point_states", positive_sum_below(
        [(1,point_exp)],-3,strict=False) and root_radius == F(1,1000))
    flow_L = F(claims["flow_generator_Lipschitz"])
    need("flow_initial_Gronwall_factor", flow_L == 400 and
         1/(1-flow_L*8*ell) < F(101,100))
    need("reconstructed_initial_inner_patch", positive_sum_below([(F(101,100),point_exp)],-3))
    need("initial_velocity_needs_only_H13", Aexp+raw0 <= claims["product_velocity_base_exponent"] == 0)

    # Full parameter chain; Q0<=1 and all four velocity Hessian terms.
    cauchy_value = claims["chart_value_exponent"]
    need("chart_value_join", Aexp+functional[0] <= cauchy_value)
    prelinear = claims["chart_prelinear_exponent"]; prequadratic = claims["chart_prequadratic_exponent"]
    need("nonzero_velocity_linear_terms", positive_sum_below([(1,functional[0]),(1,functional[1])],prelinear))
    need("nonzero_velocity_quadratic_terms", positive_sum_below([(1,functional[1]),(1,functional[2])],prequadratic))
    linear = claims["chart_linear_exponent"]; quadratic = claims["chart_quadratic_exponent"]
    need("chart_linear_coefficient", Aexp+prelinear <= linear)
    need("chart_quadratic_coefficient", 2*Aexp+prequadratic <= quadratic)
    copy_count = claims["retained_copy_count"]
    copy_linear = F(claims["copy_linear_factor"])
    copy_quadratic = F(claims["copy_quadratic_factor"])
    need("retained_identical_two_copy_scope", copy_count == 2)
    need("copy_linear_Hilbert_reserve", copy_linear > 0 and copy_linear**2 >= copy_count)
    need("copy_quadratic_chain_reserve", copy_quadratic >= copy_linear**2)
    completed = claims["completed_Hessian_exponent"]
    need("completed_uniform_Hessian", positive_sum_below(
        [(copy_linear,linear+raw2),(copy_quadratic,quadratic+2*raw1)],completed))
    need("direct_initial_size", positive_sum_below(
        [(copy_linear,cauchy_value+raw0)],claims["direct_initial_size_exponent"]))
    need("whole_interval_first_derivative", positive_sum_below(
        [(copy_linear,linear+raw1)],claims["completed_first_derivative_exponent"]))
    zero_first = claims["retained_zero_first_derivative_exponent"]
    need("retained_first_zero_jet", zero_first == 5200)
    need("initial_Taylor_size", positive_sum_below(
        [(1,aexp+zero_first),(F(1,2),2*aexp+completed)],claims["Taylor_initial_size_exponent"]))

    return {"schema": "hdblast-uniform-cauchy-arithmetic-v1", "status": "pass",
        "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(),
        "checks": records, "named_check_count": len(records),
        "exact": {"rho_scalar_integer": str(rho_scalar), "rho_oneform_integer": str(rho_oneform),
                  "correction_metric_integer": str(CH), "correction_momentum_integer": str(CM),
                  "row_totals": list(totals), "raw_exponents": [raw0,raw1,raw2],
                  "chart_coefficient_exponents": [cauchy_value,linear,quadratic],
                  "completed_Hessian_exponent": completed,
                  "initial_size_exponents": [claims["direct_initial_size_exponent"],
                                             claims["Taylor_initial_size_exponent"]]},
        "scope": "Finite rational/integer coefficient, branch and scale ledger for the accepted "
                 "conditional initial Cauchy-map proofs. No continuous analytic theorem or evolved "
                 "spacetime/observer tube is validated by this runner."}


def controls(claims):
    records = []
    for field, value in MUTATIONS:
        bad = copy.deepcopy(claims); bad[field] = value
        payload = canonical(bad)
        try:
            evaluate(json.loads(payload), lambda name, ok: None)
        except ValueError as exc:
            records.append({"fault": field, "self_rehashed_sha256": hashlib.sha256(payload).hexdigest(),
                "rejected_by": "fresh exact Cauchy coefficient and scale ledger", "diagnostic": str(exc)})
        else:
            raise ValueError("surviving Cauchy arithmetic fault: " + field)
    return {"schema": "hdblast-uniform-cauchy-arithmetic-controls-v1", "status": "pass", "controls": records}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims", type=Path, default=Path(__file__).with_name("arithmetic_claims.json"))
    parser.add_argument("--report", type=Path); parser.add_argument("--controls", type=Path)
    args = parser.parse_args(argv)
    try:
        claims = json.loads(args.claims.read_text())
        report = evaluate(claims, lambda name, ok: None)
        payload = canonical(report)
        if args.report: args.report.write_bytes(payload)
        if args.controls: args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(payload)
        return 0
    except (ValueError, KeyError, TypeError, OSError) as exc:
        sys.stderr.write(str(exc)+"\n")
        return 1
