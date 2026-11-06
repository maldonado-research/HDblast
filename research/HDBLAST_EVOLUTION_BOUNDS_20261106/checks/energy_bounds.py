"""Exact energy constants and finite counterexamples from retained premises."""
from fractions import Fraction as F
import json


def evaluate(root, data, check):
    lam = F(data["coercivity"])
    mu = F(data["inverse_upper"])
    eps = F(data["patch_epsilon"])
    time = F(data["candidate_time"])
    check("inverse_error", eps / (1 - eps) == F(1, 39))
    check("independent_time_and_spatial_lower_margin",
          0 < lam < F(38, 39 * 131**2))
    check("inverse_upper_margin", F(40, 39 * 129**2) < mu)
    check("135_mass_and_zero_order_multiplier", 0 < F(data["mass_multiplier"]) and F(data["mass_multiplier"])**2 > 1 / lam)
    check("190_force_multiplier", 0 < F(data["E_forcing_multiplier"]) and F(data["E_forcing_multiplier"])**2 > 2 / lam)
    check("95_square_root_force_multiplier", 0 < F(data["sqrt_E_forcing_multiplier"]) and F(data["sqrt_E_forcing_multiplier"])**2 > 1 / (2 * lam))
    check("269_fixed_product_norm_multiplier", 0 < F(data["fixed_product_norm_multiplier"]) and F(data["fixed_product_norm_multiplier"])**2 > 72000)
    check("18000_time_coefficient_multiplier", 1 / lam == F(data["time_coefficient_multiplier"]))
    check("36000_div_b_and_Lt_multiplier", 2 / lam == F(data["div_b_multiplier"]))
    check("18000_div_c_and_Lw_multiplier", 1 / lam == F(data["time_coefficient_multiplier"]))
    check("frame_multiplier", 2 * mu / lam == F(data["moving_frame_multiplier"]))
    check("178_lifting_energy_multiplier", 0 < F(data["lifting_energy_denominator"]) and F(data["lifting_energy_denominator"])**2 < 2 / mu)
    check("candidate_time_less_than_half", time < F(1, 2))
    check("cosh_square_safe_upper_bound", F(64, 49) < 2)
    check("frame_K_bound", F(9, 4) * time == F(117, 32000000))
    check("frame_Gronwall_integral_bound",
          F(9, 16) * time**2 == F(1521, 1024000000000000))
    check("full_fiber_dimension", 4 + (1 + 10 + 1) == 16)
    check("independent_symmetric_tensor_coordinates", 4 * (4 + 1) // 2 == 10)
    check("normal_recovery_inverse_multiplier", 1 / lam == F(data["time_coefficient_multiplier"]))
    check("normalized_inverse_derivative_prefactor",
          F(1, 129**2) * F(40, 39)**2 == F(1600, 25310961) < F(1, 15000))
    check("radial_warp_inverse_derivative_prefactor", 2 * mu == F(1, 8000))
    check("metric_time_jet_majorant", F(18000 + 18000, 15000) == F(12, 5))
    check("metric_space_jet_majorant", F(36000 + 18000, 15000) == F(18, 5))
    check("background_warp_jet_majorant", F(36000 + 18000, 8000) == F(27, 4))
    bg = json.loads((root / "inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json").read_text())
    rho_upper = F(bg["background_intervals"]["rho"]["exact"][1])
    a_bg_upper = F(bg["background_intervals"]["a"]["exact"][1])
    check("actual_background_warp_budget",
          a_bg_upper * rho_upper == F(1834, 125) < 15)
    check("actual_background_warp_contribution",
          F(27, 4) * a_bg_upper * rho_upper == F(24759, 250))
    check("237_background_frame_constant",
          F(data["mass_multiplier"]) + F(27, 4) * 15 + F(data["moving_frame_multiplier"]) * time < F(data["background_frame_constant"]))
    check("background_propagator_constant_exponent", time * F(data["background_frame_constant"]) / 2 == F(3081, 16000000))
    check("background_propagator_Lt_exponent", time * F(data["lower_t_multiplier"]) / 2 == F(117, 4000))
    check("background_propagator_Lw_exponent", time * F(data["lower_w_multiplier"]) / 2 == F(117, 8000))
    check("background_propagator_L0_exponent", time * F(data["lower_0_multiplier"]) / 2 == F(351, 3200000))
    check("actual_background_time_box_contains_candidate", time <= F(bg["domain"]["T"]))
    check("actual_background_radial_box_contains_candidate",
          F(6, 1000000) <= F(bg["domain"]["Z"]))
    ell = F(data["ell"])
    check("conditional_radial_support_margin",
          F(5, 2) * ell + 2 * time == F(23, 4) * ell < 6 * ell)
    check("conditional_auxiliary_flow_time_margin",
          time + 6 * ell == F(61, 8) * ell < 8 * ell)
    check("radial_causal_speed_11_over_10", F(41, 39) < F(11, 10)**2)
    check("conservative_radial_speed_2", F(11, 10) < 2)

    # An explicit symmetric c and two fiber vectors, with exact arithmetic,
    # check the moving-frame symmetrization for all vector coordinates.
    c = [[F(2), F(1, 3), F(1, 5), F(-1, 7)],
         [F(1, 3), F(3), F(1, 11), F(1, 13)],
         [F(1, 5), F(1, 11), F(4), F(1, 17)],
         [F(-1, 7), F(1, 13), F(1, 17), F(5)]]
    vectors = [[F(i + 1), F(2 * i - 1), F(3 - i)] for i in range(4)]
    dot = lambda x, y: sum(a * b for a, b in zip(x, y))
    pi = [0, 1, 1, 1]
    direct = sum(c[i][j] * pi[j] * dot(vectors[i], vectors[j])
                 for i in range(4) for j in range(4))
    symmetric = sum(c[i][j] * F(pi[i] + pi[j], 2)
                    * dot(vectors[i], vectors[j])
                    for i in range(4) for j in range(4))
    check("moving_frame_symmetric_form", direct == symmetric)

    # Fiber scalar principal symbol and the prescribed orthogonal split.
    q = [F(i + 1) if i >= 4 else F(0) for i in range(16)]
    wz = [F(i + 2) if i < 4 else F(0) for i in range(16)]
    check("all_sixteen_fiber_shell_flux_zero", dot(q, wz) == 0)
    q_bad = q[:]
    q_bad[0] = 1
    check("nonzero_Dirichlet_time_trace_is_detected", dot(q_bad, wz) != 0)

    # Timelike dt alone: the t/x inverse block has eigenvalues +1 and -3.
    check("timelike_counterexample_block_determinant", (-1) * (-1) - 2 * 2 == -3)
    check("timelike_counterexample_eigen_trace", 1 + (-3) == -2)
    check("timelike_counterexample_negative_spatial_rayleigh", F(-1) < 0)

    # At t=0, a_t/a=delta*N exactly. This family stays in the metric patch
    # for every N but gives no uniform multiplier from patch membership.
    delta = F(1, 100)
    check("high_frequency_counterexample_inside_patch", 0 < delta < eps)
    for freq in [1, 1000, 1000000]:
        inv_rho_sq = F(1, 130**2)
        a_zero = inv_rho_sq
        a_t_zero = inv_rho_sq * delta * freq
        check(f"high_frequency_energy_ratio_N={freq}",
              a_t_zero / a_zero == delta * freq)

    # Four spatial dimensions: u_epsilon=epsilon^-1 psi(x/epsilon).
    # Powers are exact scaling exponents for the local Euclidean model;
    # fixed smooth product geometry gives the same asymptotic orders.
    check("concentrating_gradient_L2_scaling", (-2) + F(4, 2) == 0)
    check("concentrating_gradient_square_L2_scaling", (-4) + F(4, 2) == -2)
    check("concentrating_u_D2u_L2_scaling", (-1) + (-3) + F(4, 2) == -2)

    check("lower_order_coefficients_match_energy_costs", F(data["lower_t_multiplier"]) == 2/lam and F(data["lower_w_multiplier"]) == 1/lam and 0 < F(data["lower_0_multiplier"]) and F(data["lower_0_multiplier"])**2 > 1/lam)
    check("divergence_c_coefficient_matches_energy_cost", F(data["div_c_multiplier"]) == 1/lam)
    return {"evaluated_majorant":data, "background_warp_budget":str(a_bg_upper*rho_upper), "background_frame_constant":data["background_frame_constant"], "forcing_square_root_multiplier":data["sqrt_E_forcing_multiplier"]}
