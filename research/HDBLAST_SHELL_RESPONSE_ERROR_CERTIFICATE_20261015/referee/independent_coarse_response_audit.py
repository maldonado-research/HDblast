#!/usr/bin/env python3
"""Independent coarse rational audit of the specified analytic comparator.

This program does not import the primary certificate implementation. Its
inequalities intentionally discard accuracy and use only integer/rational
supersolutions. Analytic energy/trace estimates and retained background
existence remain mathematical inputs, not machine-certified conclusions.
"""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json


def require(condition, name, results):
    results[name] = bool(condition)
    if not condition:
        raise ValueError(name)


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json"
    data = json.loads(source.read_text())
    bounds = {key: F(value) for key, value in
              data["collar"]["simple_bounds_exact"].items()}
    results = {}
    require(all(data["collar"]["gates"].values()), "retained_collar_gates", results)
    require(bounds["mu_lower"] >= F(1, 1000), "mu_floor", results)
    require(bounds["absolute_W_upper"] <= 10**6, "potential_ceiling", results)
    require(bounds["absolute_boundary_beta_upper"] <= 1000, "beta_ceiling", results)
    require(bounds["absolute_C_b_over_two_upper"] <= 250, "gamma_ceiling", results)
    require(bounds["rho_b_squared_lower"] >= 16641 and
            bounds["rho_b_squared_upper"] <= 17161, "rho_square_box", results)
    require(bounds["p_b_squared_lower"] >= F(9, 625000000) and
            bounds["p_b_squared_upper"] <= F(49, 2500000000), "p_square_box", results)

    ell = F(1, 10**10)
    time, depth = 3 * ell, 8 * ell
    require(0 < time < depth <= F(data["collar"]["inward_depth_exact"]),
            "strict_causal_triangle_inside_collar", results)

    # f=M(1+T)sqrt(L), g=16/ell+1, h=sqrt(f^2+g^2/mu).
    f_bound, g_bound, h_bound = F(31), F(161000000000), F(5100000000000)
    require(10**12 * (1 + time)**2 * depth < f_bound**2,
            "bulk_forcing_bound_31", results)
    require(16 / ell + 1 < g_bound, "boundary_forcing_bound", results)
    require(f_bound**2 + 1000 * g_bound**2 < h_bound**2,
            "combined_forcing_bound", results)

    # sqrt(1/1000)>1/32. Thus a<=max(1+M,32(1+B)).
    a_bound = F(1000001)
    require(F(1, 32)**2 < F(1, 1000) and 32 * 1001 < a_bound,
            "energy_growth_bound", results)
    exponent = a_bound * time / 2
    # exp(x)<=1/(1-x)<=1001/1000 for 0<=x<=1/1001.
    exp_bound = F(1001, 1000)
    require(0 <= exponent <= F(1, 1001) and
            1 / (1 - exponent) <= exp_bound,
            "rational_exponential_majorant", results)
    energy_bound = F(1532)
    require(exp_bound * time * h_bound < energy_bound,
            "uniform_energy_norm_bound_1532", results)

    trace_factor_bound = F(11, 2500)
    require(time * (1000 + 250**2) < trace_factor_bound**2,
            "boundary_trace_factor", results)
    delta_bound = time * (f_bound + 10**6 * energy_bound) + trace_factor_bound * energy_bound
    require(delta_bound < 8, "normal_trace_error_less_than_8", results)

    # Comparator: norm(D)^2>=2/ell (the independent note proves >=4/ell).
    require(F(140000)**2 < 2 / ell, "comparator_norm_above_140000", results)
    require(ell * (16 / ell + 1)**2 < F(1600001)**2,
            "comparator_norm_below_1600001", results)
    require(139000 < 140000 - 8 and 1600001 + 8 < 1601000,
            "true_trace_norm_coarse_interval", results)

    # The elementary Archimedean bounds 3<pi<22/7 suffice here.
    require(F(27) / (17161 * F(49, 2500000000)) > 80000,
            "prefactor_above_80000", results)
    require(3 * F(22, 7)**2 / (16641 * F(9, 625000000)) < 125000,
            "prefactor_below_125000", results)
    require(80000 * 139000**2 > 10**15,
            "leading_response_above_1e15", results)
    require(125000 * 1601000**2 < 4 * 10**17,
            "leading_response_below_4e17", results)

    output = {
        "status": "PASS_INDEPENDENT_COARSE_RATIONAL_CONSEQUENCES",
        "input_path": str(source.relative_to(root)),
        "input_sha256": sha256(source.read_bytes()).hexdigest(),
        "checks": results,
        "ell_exact": str(ell),
        "observation_time_exact": str(time),
        "causal_initial_depth_exact": str(depth),
        "physical_master_trace_error_upper": "8",
        "leading_response_interval_exact": [str(10**15), str(4 * 10**17)],
        "scope": "Conditional on retained actual-background and exact-family hypotheses, the analytic energy/characteristic trace argument, and the specified short-bump normalization. Not a waveform, physical-units prediction, nonlinear amplitude radius, or formal PDE proof.",
    }
    destination = Path(__file__).with_name("INDEPENDENT_COARSE_RESPONSE_AUDIT.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "checks": len(results),
                      "response_interval_exact": output["leading_response_interval_exact"]}))


if __name__ == "__main__":
    main()
