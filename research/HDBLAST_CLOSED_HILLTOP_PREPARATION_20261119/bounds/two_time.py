"""State-independent, two-time localization for physical n=2 hilltop modes.

Example (from the preparation study directory):
  python bounds/two_time.py --input prior/hilltop/HILLTOP.json --output bounds

Uses NumPy, SciPy and mpmath. This is a calculation on the exactly stationary
closed-FRW hilltop background, not a simulation of nonlinear quantum gravity.
"""
import argparse
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq


N = 9
V9 = 32 * math.pi**4 / 945


def prolate_lambda(c, order=96):
    """Nyström approximation to the top eigenvalue of the EXACT stated kernel."""
    if c == 0:
        return 0.0
    nodes, weights = leggauss(order)
    difference = nodes[:, None] - nodes[None, :]
    # np.sinc(z) = sin(pi*z)/(pi*z), including its removable singularity.
    kernel = (c / math.pi) * np.sinc(c * difference / math.pi)
    matrix = np.sqrt(weights[:, None] * weights[None, :]) * kernel
    return float(eigh(matrix, subset_by_index=[order - 1, order - 1],
                      check_finite=False, driver="evr")[0][0])


def prolate_threshold(probability, order=96):
    target = (2 * probability - 1)**2
    return brentq(lambda c: prolate_lambda(c, order) - target,
                  .001, 10.0, xtol=5e-14)


def transfer_mp(time, h, mu_squared, n=2):
    """Independent closed-form A,B from the Pöschl-Teller equation."""
    y = mp.mpf(h) * time
    beta = mp.sqrt(mp.mpf(mu_squared) / mp.mpf(h)**2 + mp.mpf("2.25"))
    nu = mp.mpf(n) + mp.mpf("0.5")
    z = -mp.sinh(y)**2
    prefactor = mp.cosh(y)**(-n - 2)
    a = prefactor * mp.hyp2f1((beta - nu) / 2,
                             (-beta - nu) / 2, mp.mpf("0.5"), z)
    b = mp.sinh(y) / h * prefactor * mp.hyp2f1(
        (1 + beta - nu) / 2, (1 - beta - nu) / 2, mp.mpf("1.5"), z)
    return a, b


def transfer_direct(h, mu_squared, end):
    def rhs(time, state):
        friction = 3 * h * math.tanh(h * time)
        restoring = mu_squared - 8 * h*h / math.cosh(h*time)**2
        a, ap, b, bp = state
        return [ap, restoring*a-friction*ap, bp, restoring*b-friction*bp]
    solution = solve_ivp(rhs, (0, end), [1, 0, 0, 1], method="DOP853",
                         rtol=2e-12, atol=2e-14, max_step=.025,
                         dense_output=True)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution.sol


def probability_record(c, probability):
    overlap_upper = min(1.0, math.sqrt(2*c/math.pi))
    rigorous_excess = overlap_upper / 2
    eigenvalue = prolate_lambda(c)
    sharp_excess = math.sqrt(max(0.0, eigenvalue)) / 2
    log_ball_overlap = math.log(V9) + (N/2)*(math.log(c)-math.log(2*math.pi))
    ball_overlap = min(1.0, math.exp(log_ball_overlap))
    if overlap_upper >= math.sqrt(probability):
        conditional_upper = 1.0
        conditional_excess = probability
    else:
        conditional_excess = overlap_upper**2*(2*probability-1) + 2*overlap_upper*math.sqrt(
            (1-overlap_upper**2)*probability*(1-probability))
        conditional_upper = 1-probability+conditional_excess
    return {
        "c_equal_width": c,
        "rigorous_one_mode_upper_on_min_marginal": .5 + rigorous_excess,
        "rigorous_one_mode_excess_over_half": rigorous_excess,
        "numerical_sharp_one_mode_upper_on_min_marginal": .5 + sharp_excess,
        "numerical_sharp_one_mode_excess_over_half": sharp_excess,
        "rigorous_ball_upper_on_min_marginal": .5 + ball_overlap/2,
        "rigorous_ball_excess_over_half": ball_overlap/2,
        "requested_both_marginal_probability": probability,
        "initial_marginal_probability_minimum": probability,
        "rigorous_upper_on_final_marginal_given_initial_minimum": conditional_upper,
        "rigorous_final_excess_over_one_minus_initial_minimum": conditional_excess,
    }


def calculate(source, tolerances, probabilities):
    constants = {}
    validation = {"prolate_thresholds": {}, "branches": {}}
    for p in probabilities:
        sharp = prolate_threshold(p, 96)
        higher = prolate_threshold(p, 160)
        constants[f"{p:g}"] = {
            "sharp_one_mode_c": sharp,
            "HS_one_mode_necessary_c": math.pi*(2*p-1)**2/2,
            "HS_nine_mode_ball_necessary_c": 2*math.pi*((2*p-1)/V9)**(2/N),
        }
        validation["prolate_thresholds"][f"{p:g}"] = {
            "order_96_c": sharp, "order_160_c": higher,
            "absolute_difference": abs(sharp-higher)}
    records = {}
    for name, branch in source["branches"].items():
        mu_squared = branch["negative_mass_squared_over_Htop_squared"]
        growth = branch["limiting_linear_growth_rate"]
        h = (mu_squared-growth**2)/(3*growth)
        if mu_squared <= 8*h*h:
            raise ValueError("Critical-time monotonicity requires n=2 tachyonic at bounce")
        end = max(c["primary"]["x_match"] for c in branch["detunings"].values())
        direct = transfer_direct(h, mu_squared, end)
        samples = sorted(set([0, .001, .2, 1., 3., end]))
        differences = []
        for time in samples:
            with mp.workdps(65):
                analytic_a, analytic_b = map(float, transfer_mp(time, h, mu_squared))
            numeric = direct(time)
            differences.append({"x": time,
                                "A_relative_error": abs(numeric[0]-analytic_a)/max(1.,abs(analytic_a)),
                                "B_relative_error": abs(numeric[2]-analytic_b)/max(1e-300,abs(analytic_b)) if time else 0.})
        # Check W = sech^3(hx) at high precision; cancellation makes a double
        # precision determinant an inappropriate long-time validation.
        wronskians = []
        for time in [.2, 1., 3.]:
            with mp.workdps(100):
                a, b = transfer_mp(time, h, mu_squared)
                ap = mp.diff(lambda x: transfer_mp(x, h, mu_squared)[0], mp.mpf(time))
                bp = mp.diff(lambda x: transfer_mp(x, h, mu_squared)[1], mp.mpf(time))
                expected = mp.cosh(mp.mpf(h)*time)**(-3)
                wronskians.append({"x": time,
                    "relative_error": float(abs((a*bp-ap*b)/expected-1))})
        validation["branches"][name] = {"transfer_comparison": differences,
                                       "analytic_wronskian": wronskians}

        def b_value(time):
            # Direct integration is independent of the exact hypergeometric
            # transfer used to report horizon coefficients.
            return float(direct(time)[2])

        detunings = {}
        for detuning, case in branch["detunings"].items():
            archive_sigma = case["archive_displacement_scale"]
            htop_mpl = 2*math.pi*archive_sigma
            kappa = htop_mpl**2*h**3/(2*math.pi**2)
            time = case["primary"]["x_match"]
            with mp.workdps(65):
                a, b = map(float, transfer_mp(time, h, mu_squared))
            c_scale = abs(kappa*b)
            tolerance_records = {}
            for epsilon in tolerances:
                eps_records = {
                    "one_mode_RMS_critical_time": brentq(
                        lambda t: b_value(t)-2*epsilon**2/kappa, 0, end,
                        xtol=1e-14),
                    "nine_mode_vector_RMS_critical_time": brentq(
                        lambda t: b_value(t)-2*epsilon**2/(N*kappa), 0, end,
                        xtol=1e-14),
                    "at_reference_time": probability_record(epsilon**2/c_scale, max(probabilities)),
                    "RMS_initial_width_required_for_final_epsilon": c_scale/(2*epsilon),
                    "probability_targets": {},
                }
                for p in probabilities:
                    constant = constants[f"{p:g}"]
                    sharp_c = constant["sharp_one_mode_c"]
                    hs_c = constant["HS_one_mode_necessary_c"]
                    ball_c = constant["HS_nine_mode_ball_necessary_c"]
                    eps_records["probability_targets"][f"{p:g}"] = {
                        "sharp_one_mode_last_feasible_time": brentq(
                            lambda t: b_value(t)-epsilon**2/(kappa*sharp_c), 0, end, xtol=1e-14),
                        "HS_one_mode_latest_necessary_time": brentq(
                            lambda t: b_value(t)-epsilon**2/(kappa*hs_c), 0, end, xtol=1e-14),
                        "HS_nine_mode_ball_latest_necessary_time": brentq(
                            lambda t: b_value(t)-epsilon**2/(kappa*ball_c), 0, end, xtol=1e-14),
                        "sharp_one_mode_initial_halfwidth_required_at_reference_time": sharp_c*c_scale/epsilon,
                    }
                tolerance_records[f"{epsilon:g}"] = eps_records
            detunings[detuning] = {
                "Htop_over_Mpl": htop_mpl, "kappa": kappa,
                "reference_homogeneous_departure_time_x_match": time,
                "A_at_reference_time": a, "B_at_reference_time": b,
                "absolute_two_time_commutator_coefficient": c_scale,
                "tolerances": tolerance_records,
            }
        records[name] = {"h_star": h, "a_b": 1/h,
                         "mu_squared": mu_squared, "n": 2,
                         "detunings": detunings}
    max_transfer_error = max(row[key]
        for branch in validation["branches"].values()
        for row in branch["transfer_comparison"]
        for key in ["A_relative_error", "B_relative_error"])
    max_wronskian_error = max(row["relative_error"]
        for branch in validation["branches"].values()
        for row in branch["analytic_wronskian"])
    max_prolate_error = max(row["absolute_difference"]
        for row in validation["prolate_thresholds"].values())
    validation.update(max_transfer_relative_error=max_transfer_error,
                      max_wronskian_relative_error=max_wronskian_error,
                      max_prolate_threshold_absolute_difference=max_prolate_error)
    if max_transfer_error > 2e-9 or max_wronskian_error > 1e-55 or max_prolate_error > 1e-10:
        raise AssertionError("Independent numerical validation failed")
    return {"probability_constants": constants, "branches": records,
            "validation": validation}


def report(result, input_path):
    first_branch = next(iter(result["branches"].values()))
    first_case = next(iter(first_branch["detunings"].values()))
    tolerances = list(first_case["tolerances"])
    tolerance_header = "| branch | detuning | " + " | ".join(
        "epsilon=" + value for value in tolerances) + " |"
    tolerance_separator = "|---|---:|" + "---:|"*len(tolerances)
    lines = ["# Two-time state-independent localization on the stationary hilltop", "",
             f"Input: `{input_path}`. All times are x = Htop t. n=2 means eigenvalue 8 on the unit S³.", "",
             "Equal interval halfwidths are epsilon; probabilities concern one initial density operator, evolved without intermediate measurement.", "",
             "| branch | detuning | kappa | homogeneous x_match (reference) | absolute kappa B |",
             "|---|---:|---:|---:|---:|"]
    for name, branch in result["branches"].items():
        for detuning, case in branch["detunings"].items():
            lines.append(f"| {name} | {detuning} | {case['kappa']:.8e} | {case['reference_homogeneous_departure_time_x_match']:.8f} | {case['absolute_two_time_commutator_coefficient']:.8e} |")
    for p in result["probability_constants"]:
        lines += ["", f"## Last time admitting both one-mode marginal probabilities >= {p}", "",
                  "This uses the sharp interval/prolate bound, evaluated numerically. The exact operator identity and rigorous simpler bounds are in THEORY.md.", "",
                  tolerance_header, tolerance_separator]
        for name, branch in result["branches"].items():
            for detuning, case in branch["detunings"].items():
                values = [row["probability_targets"][p]["sharp_one_mode_last_feasible_time"]
                          for row in case["tolerances"].values()]
                lines.append(f"| {name} | {detuning} | " + " | ".join(f"{x:.9g}" for x in values) + " |")
    lines += ["", "## At the homogeneous departure reference time", "",
              "The rigorous HS upper bound on min(p0,pT) equals 1/2 plus the values below. The exact sharp prolate ceiling is never larger. This is not a bound on one specified marginal alone.", "",
              tolerance_header, tolerance_separator]
    for name, branch in result["branches"].items():
        for detuning, case in branch["detunings"].items():
            values = [row["at_reference_time"]["rigorous_one_mode_excess_over_half"]
                      for row in case["tolerances"].values()]
            lines.append(f"| {name} | {detuning} | " + " | ".join(f"{x:.6e}" for x in values) + " |")
    initial_minimum = first_case["tolerances"][tolerances[0]]["at_reference_time"]["initial_marginal_probability_minimum"]
    largest_final_excess = max(row["at_reference_time"]["rigorous_final_excess_over_one_minus_initial_minimum"]
        for branch in result["branches"].values()
        for case in branch["detunings"].values()
        for row in case["tolerances"].values())
    lines += ["", f"If the initial unmeasured marginal probability is already at least {initial_minimum:g}, the stronger Landau-Pollak angle relation bounds the final unmeasured marginal by {1-initial_minimum:.6g} plus at most {largest_final_excess:.6e} across these reference cases. The JSON preserves the small corrections separately."]
    validation = result["validation"]
    lines += ["", "## Validation", "",
              f"Direct DOP853 transfer vs independent hypergeometric solution: maximum relative discrepancy {validation['max_transfer_relative_error']:.3e}.",
              f"100-digit analytic Wronskian check: maximum relative discrepancy {validation['max_wronskian_relative_error']:.3e}.",
              f"Prolate 96- vs 160-node critical-c difference: maximum {validation['max_prolate_threshold_absolute_difference']:.3e}.", "",
              "The prolate quadrature convergence is numerical evidence, not a rigorous interval eigenvalue certificate. The Hilbert-Schmidt probability inequality itself is rigorous.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tolerances", nargs="+", type=float, default=[.01, .001, .0001])
    parser.add_argument("--probabilities", nargs="+", type=float, default=[.95, .99])
    args = parser.parse_args()
    if any(not 0 < e < 1 for e in args.tolerances):
        parser.error("tolerances must lie in (0,1)")
    if any(not .5 < p < 1 for p in args.probabilities):
        parser.error("probabilities must lie in (.5,1)")
    result = calculate(json.loads(args.input.read_text()), args.tolerances, args.probabilities)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "BOUNDS.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    (args.output / "BOUNDS.md").write_text(report(result, args.input))
    print(json.dumps({key: value for key, value in result["validation"].items()
                      if key.startswith("max_")}, indent=2))


if __name__ == "__main__":
    main()
