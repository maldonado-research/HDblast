"""Check the hilltop-preparation capsule and optionally regenerate its numerical results."""
import argparse
import csv
import hashlib
import json
import math
import os
import runpy
from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
BASE = "b0e8280d7bac82e105b85ee33df0ae8745228a5f"
SOURCE_MANIFEST_SHA256 = "8affe8e0819d5e7fbbba8f47dff5a18999ff117695494cfc754989d9f17fd024"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def authenticate():
    inventory = HERE / "SHA256SUMS"
    listed = set()
    for line in inventory.read_text().splitlines():
        checksum, filename = line.split("  ", 1)
        path = HERE / filename
        require(not Path(filename).is_absolute() and path.resolve() == path and path.is_relative_to(HERE) and path != inventory
                and path not in listed, f"Invalid inventory entry {filename}")
        require(digest(path) == checksum, f"Invalid study checksum for {filename}")
        listed.add(path)
    present = {path for path in HERE.rglob("*") if path.is_file() and path != inventory
               and "__pycache__" not in path.relative_to(HERE).parts}
    require(listed == present, f"Study inventory differs: {present ^ listed}")
    manifest = HERE / "SOURCES.json"
    require(digest(manifest) == SOURCE_MANIFEST_SHA256, "Historical source manifest changed")
    pins = read(manifest)
    require(pins["base"] == BASE, "Historical revision mismatch")
    for source in pins["sources"]:
        path = HERE / source["bundled"]
        require(path in listed and path.stat().st_size == source["bytes"] and digest(path) == source["sha256"],
                f"Historical input differs from {source['path']}")
    print(f"Verified {len(listed)} study files, including {len(pins['sources'])} historical inputs.", flush=True)


def number_tolerance(path):
    # Keep a zero absolute tolerance on physical values, including widths below 1e-90.
    # Only diagnostic differences and derivatives reconstructed by cancellation have absolute floors.
    key = path[-1]
    if path[0] == "BOUNDS" and "validation" in path:
        if "analytic_wronskian" in path or key == "max_wronskian_relative_error":
            return 1e-6, 1e-55 if key != "x" else 0
        if "error" in key or "difference" in key:
            return 1e-6, 2e-9 if "error" in key else 1e-10
    if key in {"stationary_theta", "top"} and path[0] in {"NONLINEAR", "HARMONICS"}:
        return 1e-6, 1e-14
    if path[0] == "NONLINEAR":
        floors = {"basis_orthonormality_error": 2e-14, "removed_force_at_top": 5e-14,
                  "max_abs_energy_continuity_residual": 1e-12, "energy_residual_over_largest_component": 1e-8,
                  "stop_x_difference": 1e-9, "stop_relative_mean_difference": 1e-8, "odd_to_even": 1e-7,
                  "odd_mean": 1e-16, "max_relative_retarded_integral_error": 1e-8}
        return 1e-6, 1e-13 if "mode_coefficients" in path else floors.get(key, 0)
    if path[0] == "VELOCITY":
        floors = {"max_scaled_Hamiltonian_constraint": 2e-9, "max_Friedmann_error": 2e-12,
                  "max_relative_symplectic_drift": 2e-9, "max_relative_canonical_transfer_difference": 2e-9,
                  "horizon_refinement_difference": 1e-9, "max_horizon_change_under_center_bracket": 1e-9}
        return 1e-6, floors.get(key, 0)
    if path[0] == "ROLLING":
        floors = {"constraint_max_scaled_hamiltonian": 2e-8, "independent_trace_max_scaled_momentum": 2e-9,
                  "max_friedmann_error": 2e-10, "trace_vs_momentum_fractional_difference": 2e-9,
                  "symplectic_max_relative_drift": 2e-8, "metric_source_fractional_change_delta_theta": 1e-9,
                  "metric_source_fractional_change_delta_theta_prime": 1e-9}
        return 1e-6, floors.get(key, 0)
    if path[0] == "PRECISION":
        floors = {"match_q_relative_error": 2e-7, "match_loga_absolute_error": 2e-9,
                  "center_refinement_difference": 1e-30}
        return 1e-6, floors.get(key, 0)
    if path[0] == "MANIFOLD":
        floors = {"relative_match_p_error": 5e-8, "match_h_error": 2e-10,
                  "center_refinement_difference": 1e-11, "direct_ivp_relative_difference": 2e-6,
                  "quadratic_residual": 5e-15, "max_constraint": 2e-9,
                  "finite_difference_quadratic": 1e-8, "Uthird_right": 1e-7, "Uthird_left": 1e-7}
        return 1e-6, floors.get(key, 0)
    if ("error" in key or "Wronskian_error" in key) and path[0] in {"MINIMAX", "HARMONICS"}:
        return 1e-6, 1e-8
    return 1e-6, 0


def compare(actual, expected, path):
    label = "/".join(path)
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), f"Object fields differ in {label}")
        for key, value in expected.items():
            compare(actual[key], value, (*path, key))
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), f"List size differs in {label}")
        for index, (new, saved) in enumerate(zip(actual, expected, strict=True)):
            compare(new, saved, (*path, str(index)))
    elif path[0] == "PRECISION" and isinstance(expected, str):
        require(isinstance(actual, str), f"Decimal string required for {label}")
        new, old = Decimal(actual), Decimal(expected)
        # The first correction measures the difference from a double-precision BVP starting guess.
        # Later corrections compare successive finite-time stable boundary value problems.
        floor = Decimal("5e-14") if "finite_T_changes" in path and path[-1] == "0" else Decimal("1e-44")
        require(new.is_finite() and old.is_finite() and abs(new-old) <= max(floor, abs(old)*Decimal("0.000001")),
                f"{label}: computed {new}, archived {old}")
    elif isinstance(expected, float):
        rel, floor = number_tolerance(path)
        require(type(actual) in (int, float) and math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual, expected, rel_tol=rel, abs_tol=floor),
                f"{label}: computed {actual!r}, archived {expected!r}")
    else:
        require(type(actual) is type(expected) and actual == expected, f"Value differs in {label}")


def check_transfer_csv(new_path, saved_path):
    rows = {}
    with new_path.open(newline="") as new_stream, saved_path.open(newline="") as saved_stream:
        fresh, saved = csv.DictReader(new_stream), csv.DictReader(saved_stream)
        require(fresh.fieldnames == saved.fieldnames, "Transfer columns differ")
        for index, (new, old) in enumerate(zip(fresh, saved, strict=True), 1):
            for field in fresh.fieldnames:
                if field in {"branch", "n", "degeneracy", "k_squared"}:
                    require(new[field] == old[field], f"Transfer field {field} differs at row {index}")
                else:
                    compare(float(new[field]), float(old[field]), ("TRANSFER", str(index), field))
            n = int(new["n"])
            key = (new["branch"], float(new["x"]), n)
            require(key not in rows and int(new["degeneracy"]) == (n + 1)**2
                    and int(new["k_squared"]) == n*(n + 2), f"Transfer harmonic mismatch at row {index}")
            rows[key] = new
    require(bool(rows), "Empty transfer table")
    return rows


def near(value, reference, label, relative=2e-8):
    require(math.isfinite(value) and math.isfinite(reference)
            and math.isclose(value, reference, rel_tol=relative), f"{label}: {value} versus {reference}")


def cross_checks(all_results, transfer):
    historical = read(HERE / "prior/hilltop/HILLTOP.json")["branches"]
    for name, branch in historical.items():
        modes = all_results["HARMONICS"]["branches"][name]
        bounds = all_results["BOUNDS"]["branches"][name]
        stable = all_results["MANIFOLD"]["branches"][name]
        near(modes["mu_squared_over_Htop_squared"], bounds["mu_squared"], f"Hilltop mass in {name}")
        beta = bounds["mu_squared"] / bounds["h_star"]**2
        levels = [n for n in range(1 + math.ceil(math.sqrt(beta))) if n*(n + 3) < beta]
        degeneracy = sum((n+1)*(n+2)*(2*n+3)//6 for n in levels)
        require((max(levels), degeneracy) == {"phi_m=0.5": (12, 3185), "phi_m=-0.5": (10, 1716)}[name],
                f"Euclidean scalar spectrum in {name}")

        nonlinear = all_results["NONLINEAR"]["branches"][name]
        near(nonlinear["mu_squared"], bounds["mu_squared"], "Nonlinear potential curvature")
        for coefficient, side in zip(nonlinear["sided_quadratic_force_coefficients"], ("left", "right"), strict=True):
            near(coefficient, -1.5*stable[f"Uthird_{side}"], f"Nonlinear source on {name}, {side}")
        require(all(value == 0 for value in nonlinear["zero_control"]["at_stop"]["mode_coefficients"]),
                "Stationary nonlinear control")
        precision = all_results["PRECISION"][name]
        sample = next(row for row in stable["samples"] if row["v"] == precision["v"])
        near(float(abs(Decimal(precision["fine"]["signed_separation_decimal"]))), sample["conditional_width"],
             f"Resolved off-center launch width in {name}")
        require(abs(precision["center_refinement_difference"]) < abs(sample["conditional_width"])*1e-5,
                f"Finite-time decimal shoot refinement in {name}")
        rest_seed = stable["critical_zero_v_distance"]
        near(rest_seed, math.exp(branch["detunings"]["0.001"]["primary"]["log_abs_bounce_seed"]),
             f"Critical zero-velocity seed in {name}")
        for sample in stable["samples"]:
            near(sample["conditional_width"], rest_seed*sample["conditional_width_ratio"], "Conditional width")
        for launch in all_results["VELOCITY"]["branches"][name]:
            sample = next(row for row in stable["samples"] if row["v"] == launch["velocity"])
            near(launch["calibrated_distance_to_center"], sample["conditional_width"], "Compensated launch width")
            near(launch["center_offset"], sample["center"], "Compensated launch center")
            near(launch["kappa_in_q_P"]*launch["B_in_q_P"], launch["kappa_in_q_qprime"]*launch["B_in_q_qprime"],
                 "Agreement of canonical charts at the RMS horizon")
            near(launch["canonical_rms_floor"], all_results["VELOCITY"]["rms_goal"], "Compensated RMS horizon")
            require(launch["horizon_x"] < launch["prior_paired_match_x"], "Compensated horizon precedes handoff")
        for event in all_results["EVENTS"][name]:
            before_turn = event["first_event"] == "equality"
            require(before_turn == (event["offset_multiplier"] < 1)
                    and before_turn == (event["ratio_at_maximum"] > 1), f"Event bracket in {name}")

        for detuning, history in branch["detunings"].items():
            case = bounds["detunings"][detuning]
            optimum = all_results["MINIMAX"]["branches"][name][detuning]
            rolling_case = all_results["ROLLING"]["branches"][name][detuning]
            hcase = modes["detunings"][detuning]
            time = history["primary"]["x_match"]
            near(hcase["homogeneous_linear_handoff_x"], time, "Finite-band reference time")
            near(optimum["time_required"], time, "Minimum RMS reference time")
            near(case["reference_homogeneous_departure_time_x_match"], time, "Probability reference time")
            row = transfer[(name, time, 2)]
            near(float(row["q_from_q"]), case["A_at_reference_time"], "Independent displacement transfer")
            B = float(row["q_from_v"])
            for B_test in (case["B_at_reference_time"], optimum["n2_B"],
                           rolling_case["bounce_quantization"]["samples"]["at_prior_match"]["B_stationary"]):
                near(B, B_test, f"Independent velocity transfer in {name} at {detuning}")
            kappa = case["kappa"]
            near(kappa*abs(B), case["absolute_two_time_commutator_coefficient"], "Two-time commutator")
            near(kappa, optimum["kappa"], "Field normalization")
            near(kappa, rolling_case["bounce_quantization"]["kappa_per_real_harmonic"], "Rolling launch normalization")
            for label, band in (("n2", [2]), ("n2_to_8", range(2, 9))):
                sum_B = sum((n+1)**2 * float(transfer[(name, time, n)]["q_from_v"]) for n in band)
                near(math.sqrt(kappa*sum_B/2), optimum["band_results"][label]["minimax_rms_at_required_time"],
                     f"All-time band RMS from the harmonic CSV for {name} at {detuning}")
            horizons = optimum["band_results"]["n2"]["rms_horizons"]
            for epsilon, row in case["tolerances"].items():
                near(row["nine_mode_vector_RMS_critical_time"], horizons[epsilon], "RMS horizon from two solvers")
            near(horizons["0.01"], rolling_case["bounce_quantization"]["rolling_horizon_x"], "Rolling RMS horizon")
            near(math.exp(history["primary"]["log_abs_bounce_seed"]), rolling_case["abs_seed"], "Rolling seed")


def replay(jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    prior = HERE / "prior"
    with tempfile.TemporaryDirectory(prefix="hdblast-preparation-") as directory:
        fresh = Path(directory)

        def run(script, *arguments):
            print(f"Recomputing {script}", flush=True)
            subprocess.run([sys.executable, "-I", "-B", str(HERE/script), *map(str, arguments)],
                           cwd=HERE, env=env, check=True)

        products = [("harmonics", "HARMONICS"), ("minimax", "MINIMAX"), ("bounds", "BOUNDS"),
                    ("rolling", "ROLLING"), ("manifold", "MANIFOLD"), ("manifold", "EVENTS"),
                    ("manifold", "PRECISION"), ("nonlinear", "NONLINEAR"), ("rolling", "VELOCITY")]
        with ThreadPoolExecutor(max_workers=jobs) as workers:
            manifold = workers.submit(run, "manifold/manifold.py", "--prior", prior, "--output", fresh/"MANIFOLD.json")
            tasks = [
                workers.submit(run, "harmonics/harmonics.py", "--prior", prior, "--output-dir", fresh),
                workers.submit(run, "minimax/minimax.py", "--prior", prior, "--output", fresh/"MINIMAX.json"),
                workers.submit(run, "bounds/two_time.py", "--input", "prior/hilltop/HILLTOP.json", "--output", fresh),
                workers.submit(run, "rolling/rolling.py", "--prior", prior, "--output", fresh/"ROLLING.json"),
                workers.submit(run, "nonlinear/nonlinear.py", "--prior", prior, "--output", fresh/"NONLINEAR.json"),
            ]
            manifold.result()
            tasks.append(workers.submit(run, "manifold/events.py", "--prior", prior,
                                        "--manifold", fresh/"MANIFOLD.json", "--output", fresh/"EVENTS.json"))
            tasks.append(workers.submit(run, "manifold/precision.py", "--prior", prior,
                                        "--manifold", fresh/"MANIFOLD.json", "--output", fresh/"PRECISION.json"))
            tasks.append(workers.submit(run, "rolling/velocity.py", "--prior", prior,
                                        "--manifold", fresh/"MANIFOLD.json", "--output", fresh/"VELOCITY.json"))
            for task in tasks:
                task.result()

        result = {}
        for section, name in products:
            data = read(fresh/(name + ".json"))
            compare(data, read(HERE/section/(name + ".json")), (name,))
            result[name] = data
        rows = check_transfer_csv(fresh/"TRANSFER.csv", HERE/"harmonics/TRANSFER.csv")
        cross_checks(result, rows)
        # The generator prints rounded floating-point validation residuals. Check the publication table
        # against the archived numerical values after the fresh values have passed their replay tolerances.
        formatter = runpy.run_path(str(HERE/"bounds/two_time.py"))["report"]
        saved_bounds = read(HERE/"bounds/BOUNDS.json")
        require(formatter(saved_bounds, Path("prior/hilltop/HILLTOP.json")) == (HERE/"bounds/BOUNDS.md").read_text(),
                "Published probability table disagrees with its archived numerical data")
        print(f"Verified {len(result)} numerical JSON results, {len(rows)} transfer rows, and cross-study checks.",
              flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Run numerical calculations and compare all result fields")
    parser.add_argument("--jobs", default=2, type=int, help="Maximum simultaneous numerical processes")
    args = parser.parse_args()
    require(args.jobs >= 1, "--jobs must be positive")
    authenticate()
    if args.verify:
        replay(args.jobs)
