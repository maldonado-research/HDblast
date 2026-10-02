#!/usr/bin/env python3
"""Post-run archive audit; no source evaluation or mode evolution is performed."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import traceback

import numpy as np


PUBLIC_FREEZE = "a01d17e015f070be2ad2018745e877148fe2f4e4"
REGISTRATION_SHA256 = "f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6"
AUDIT_ARITHMETIC_ALLOWANCE = 128 * np.finfo(float).eps
PANEL_MOMENT_ALLOWANCE = 1e-10


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def close(left, right, tolerance, label):
    require(math.isfinite(left) and math.isfinite(right)
            and abs(left-right) <= tolerance, f"{label}: {left!r} versus {right!r}")


def indexed_rows(rows):
    result = {(row["source"], row["eta"]): row for row in rows}
    require(len(result) == len(rows), "Duplicate source/observation rows")
    return result


def inspect(checkpoint, primary_path, modes_path, checks_path):
    registration_path = checkpoint / "FULL_REGISTRATION.json"
    require(sha(registration_path) == REGISTRATION_SHA256,
            "Full registration does not match the public freeze anchor")
    registration = read_json(registration_path)
    for name, digest in registration["files"].items():
        path = (checkpoint / name).resolve()
        require(path.is_relative_to(checkpoint), "Frozen input escaped checkpoint")
        require(sha(path) == digest, f"Frozen payload hash mismatch: {name}")
    experiment = read_json(checkpoint / "EXPERIMENT.json")
    primary = read_json(primary_path)
    modes = read_json(modes_path)
    checks = read_json(checks_path)
    gates = experiment["gates"]
    epsilon = experiment["source"]["epsilon"]
    expected_grid = {(source, eta) for source in experiment["source"]["ids"]
                     for eta in experiment["observations"]}
    primary_rows = indexed_rows(primary["rows"])
    mode_rows = indexed_rows(modes["rows"])
    require(set(primary_rows) == set(mode_rows) == expected_grid, "Observation grid mismatch")
    require(checks["status"] == "passed", "Original numerical validator did not pass")
    require(modes["status"] == "passed_internal_gates", "Mode producer did not pass")
    require(modes["freeze_commit"] == primary["provenance"]["public_freeze_commit"]
            == PUBLIC_FREEZE, "Public freeze identifier mismatch")
    require(modes["epsilon"] == epsilon, "Mode source amplitude mismatch")
    require(modes["runtime"]["numpy"] == np.__version__ == "2.2.6", "Pinned NumPy version mismatch")
    require(platform.python_version() == "3.12.14", "Pinned Python version mismatch")
    evidence = checks["evidence"]
    require(evidence["modes_sha256"] == sha(modes_path), "Validator modes evidence mismatch")
    require(evidence["primary_sha256"] == sha(primary_path), "Validator primary evidence mismatch")
    require(evidence["registration_sha256"] == REGISTRATION_SHA256,
            "Validator registration evidence mismatch")
    require(evidence["validator_sha256"] == sha(checkpoint / "code/validate_causal.py"),
            "Validator source evidence mismatch")
    controls = checks["wrong_formula_controls"]
    require(len(controls) == 17 and len({c["name"] for c in controls}) == 17
            and all(c["rejected"] is True for c in controls), "Mutation status/count mismatch")
    require(len(checks["comparisons"]) == 36, "Validator comparison count mismatch")
    manifest_path = checkpoint / "independent/MANIFEST.json"
    manifest = read_json(manifest_path)
    require(modes["provenance"]["manifest_sha256"] == sha(manifest_path), "Mode manifest mismatch")
    require(modes["provenance"]["files"] == manifest["files"], "Mode file-pin list mismatch")
    require(modes["resources"]["elapsed_seconds"] <= 900
            and modes["resources"]["peak_rss_kib"] <= 128*1024, "Recorded mode resource gate")
    required_runs = {(source, setting) for source in experiment["source"]["ids"]
                     for setting in ("coarse", "fine")}
    require({(r["source"], r["setting"]) for r in modes["runs"]} == required_runs
            and len(modes["runs"]) == 4, "Mode run coverage mismatch")

    run_audits = []
    reconstructed = {}
    peak_integral_difference = 0.0
    peak_amplitude_wronskian = 0.0
    peak_physical_wronskian = 0.0
    peak_moment_error = 0.0
    panel_moment_checks = 0
    for run in modes["runs"]:
        source, setting = run["source"], run["setting"]
        fine = setting == "fine"
        width = 0.25 if fine else 0.5
        count = 16384 if fine else 8192
        n_panels = count // 16
        require(run["momentum_nodes"] == count and run["momentum_panel_width"] == width,
                "Reported momentum settings mismatch")
        require(run["time_steps"] == (1152 if fine else 576)
                and run["dt"] == 1/(256 if fine else 128), "Reported time settings mismatch")
        require(run["time_gauss_nodes"] == 8
                and run["local_time_momentum_pairs"] == run["time_steps"]*8*count,
                "Reported local forcing budget mismatch")
        for label in ("wronskian_max_scaled", "physical_wronskian_max_scaled"):
            require(0 <= run[label] <= gates["mode_wronskian_over_epsilon"],
                    "Reported global Wronskian gate failed")
        archive_path = (modes_path.parent / run["archive"]["path"]).resolve()
        require(archive_path.is_relative_to(modes_path.parent.resolve()), "Archive path escaped modes directory")
        require(sha(archive_path) == run["archive"]["sha256"], "Archive SHA-256 mismatch")
        observations = []
        with np.load(archive_path, allow_pickle=False) as data:
            expected_names = {"k", "momentum_weights"}
            expected_names |= {f"{name}_{i}" for name in ("u", "w", "subtraction", "combined")
                               for i in range(len(experiment["observations"]))}
            require(set(data.files) == expected_names, "Unexpected archive members")
            k, weights = data["k"], data["momentum_weights"]
            require(k.shape == weights.shape == (count,) and k.dtype == weights.dtype == np.float64,
                    "Momentum archive shape/dtype mismatch")
            require(np.all(np.isfinite(k)) and np.all(np.isfinite(weights))
                    and np.all(k > 0) and np.all(k < 256) and np.all(np.diff(k) > 0)
                    and np.all(weights > 0), "Invalid archived momentum grid")
            left = np.arange(n_panels)[:, None] * width
            x = (k.reshape(n_panels, 16)-left) * (2/width)-1
            scaled_weights = weights.reshape(n_panels, 16) * (2/width)
            require(np.all(np.abs(x) < 1), "Momentum node outside its panel")
            powers = np.ones_like(x)
            for degree in range(32):
                expected_moment = 2/(degree+1) if degree % 2 == 0 else 0.0
                measured = np.sum(scaled_weights*powers, axis=1)
                error = float(np.max(np.abs(measured-expected_moment)))
                peak_moment_error = max(peak_moment_error, error)
                require(error <= PANEL_MOMENT_ALLOWANCE, "GL16 panel polynomial moment failed")
                powers *= x
                panel_moment_checks += n_panels
            require([row["eta"] for row in run["rows"]] == experiment["observations"],
                    "Raw observation ordering mismatch")
            for i, row in enumerate(run["rows"]):
                eta = row["eta"]
                require(row["source"] == source, "Raw observation source mismatch")
                canonical = mode_rows[(source, eta)]
                primary_row = primary_rows[(source, eta)]
                a = -1/eta
                close(row["a"], a, 0.0, "Raw scale factor")
                close(canonical["a"], a, 0.0, "Summary scale factor")
                close(row["s_over_epsilon"], canonical["s_over_epsilon"], 0.0, "Raw/summary source")
                close(row["s_over_epsilon"], primary_row["s_over_epsilon"],
                      AUDIT_ARITHMETIC_ALLOWANCE, "Independent source records")
                u, w = data[f"u_{i}"], data[f"w_{i}"]
                subtraction, combined = data[f"subtraction_{i}"], data[f"combined_{i}"]
                require(u.dtype == w.dtype == np.complex128
                        and subtraction.dtype == combined.dtype == np.float64,
                        "Mode archive dtype mismatch")
                require(all(array.shape == k.shape and np.all(np.isfinite(array))
                            for array in (u, w, subtraction, combined)), "Invalid archived mode array")
                s0 = epsilon*row["s_over_epsilon"]
                mass2 = 2*a*a
                expected_subtraction = s0*k*k / (k*k+mass2)**1.5
                require(np.max(np.abs(subtraction-expected_subtraction)) <=
                        AUDIT_ARITHMETIC_ALLOWANCE*epsilon, "Archived local subtraction mismatch")
                mode_term = 4*k*u.real
                require(np.array_equal(mode_term+subtraction, combined),
                        "Archived pointwise combined integrand mismatch")
                wronskian = float(np.max(np.abs(2*u.real-w.imag/k))/epsilon)
                phase = np.cos(k*eta)-1j*np.sin(k*eta)
                v = phase/np.sqrt(2*k)
                dv = v*u
                vp = -1j*k*v
                dvp = v*(w-1j*k*u)
                physical = float(np.max(np.abs(dv*np.conj(vp)+v*np.conj(dvp)
                                               -dvp*np.conj(v)-vp*np.conj(dv)))/epsilon)
                require(wronskian <= gates["mode_wronskian_over_epsilon"]
                        and physical <= gates["mode_wronskian_over_epsilon"],
                        "Reconstructed observation Wronskian gate failed")
                close(wronskian, row["wronskian_scaled"], AUDIT_ARITHMETIC_ALLOWANCE,
                      "Archived amplitude Wronskian record")
                close(physical, row["physical_wronskian_scaled"], AUDIT_ARITHMETIC_ALLOWANCE,
                      "Archived physical Wronskian record")
                require(wronskian <= run["wronskian_max_scaled"]+AUDIT_ARITHMETIC_ALLOWANCE
                        and physical <= run["physical_wronskian_max_scaled"]+AUDIT_ARITHMETIC_ALLOWANCE,
                        "Observation exceeds reported all-step Wronskian maximum")
                peak_amplitude_wronskian = max(peak_amplitude_wronskian, wronskian)
                peak_physical_wronskian = max(peak_physical_wronskian, physical)
                if eta < -5:
                    require(all(np.count_nonzero(array) == 0 for array in (u, w, subtraction, combined)),
                            "Nonzero pre-pulse mode snapshot")
                high = k >= 128
                high_mode = float(np.max(np.abs(mode_term[high])))
                high_subtraction = float(np.max(np.abs(subtraction[high])))
                high_combined = float(np.max(np.abs(combined[high])))
                scaled_high = float(np.max(np.abs(combined[high])*k[high]**3)/epsilon)
                budget = primary_row["tail_derivative_budget"]
                coefficient_bound = (1.5*budget["source_abs_upper"]/epsilon*mass2
                                     + budget["normalized_A_upper"]/4)
                require(scaled_high <= coefficient_bound+1e-8,
                        "Archived high-k integrand exceeds analytic coefficient diagnostic")
                finite_k = []
                for point, summary in zip(row["finite_k"], canonical["finite_k"]):
                    cutoff = point["K"]
                    require(cutoff == summary["K"], "Raw/summary cutoff mismatch")
                    # Different summation grouping from the producer: directly
                    # sum all individual weighted nodes below the cutoff.
                    y = math.fsum(float(z) for z in (weights[k < cutoff]*combined[k < cutoff]))
                    y /= 8*math.pi**2*epsilon
                    delta = abs(y-point["y"])
                    peak_integral_difference = max(peak_integral_difference, delta)
                    close(y, point["y"], AUDIT_ARITHMETIC_ALLOWANCE, "Archived finite-K reconstruction")
                    close(point["y"], summary["y" if fine else "y_coarse"], 0.0,
                          "Raw/summary response record")
                    reconstructed[(source, eta, setting, cutoff)] = y
                    finite_k.append({"K": cutoff, "reconstructed_y": y,
                                     "producer_difference": delta})
                observations.append({"eta": eta, "wronskian_scaled": wronskian,
                                     "physical_wronskian_scaled": physical,
                                     "high_k_max_abs_mode": high_mode,
                                     "high_k_max_abs_subtraction": high_subtraction,
                                     "high_k_max_abs_combined": high_combined,
                                     "high_k_cancellation_ratio": max(high_mode, high_subtraction)/high_combined
                                         if high_combined else None,
                                     "high_k_max_abs_k_cubed_integrand_over_epsilon": scaled_high,
                                     "analytic_coefficient_bound": coefficient_bound,
                                     "finite_k": finite_k})
        run_audits.append({"source": source, "setting": setting,
                           "archive": {"name": archive_path.name, "sha256": sha(archive_path)},
                           "momentum_panels": n_panels, "momentum_nodes": count,
                           "observations": observations})

    comparisons = []
    maximum_cross_route = 0.0
    maximum_refinement = 0.0
    maximum_bound_ratio = 0.0
    reported_comparisons = {(r["source"], r["eta"], r["K"]): r for r in checks["comparisons"]}
    require(len(reported_comparisons) == 36, "Duplicate validator comparison")
    for (source, eta), row in primary_rows.items():
        for primary_point, mode_point in zip(row["finite_k"], mode_rows[(source, eta)]["finite_k"]):
            cutoff = primary_point["K"]
            fine = reconstructed[(source, eta, "fine", cutoff)]
            coarse = reconstructed[(source, eta, "coarse", cutoff)]
            refinement = abs(fine-coarse)
            cross = abs(fine-primary_point["y"])
            require(refinement <= gates["mode_refinement_absolute"], "Reconstructed refinement gate")
            require(cross <= gates["mode_cross_route_absolute"]+primary_point["quadrature_error_estimate"],
                    "Reconstructed finite-K cross-route gate")
            close(mode_point["refinement_difference"], abs(mode_point["y"]-mode_point["y_coarse"]),
                  0.0, "Recorded refinement difference")
            close(mode_point["estimated_numerical_error"], 2*mode_point["refinement_difference"],
                  0.0, "Recorded empirical error estimate")
            bound = primary_point["analytic_tail_bound"]
            budget = row["tail_derivative_budget"]
            coefficient = (1.5*budget["source_abs_upper"]/epsilon*(2/eta**2)
                           +budget["normalized_A_upper"]/4)
            reconstructed_bound = coefficient/(16*math.pi**2*cutoff**2)
            close(bound, reconstructed_bound, max(1e-20, abs(bound)*1e-12),
                  "Recorded analytic omitted-tail formula")
            difference_from_continuum = abs(fine-row["continuum"]["y"])
            allowance = (bound+gates["mode_cross_route_absolute"]
                         +primary_point["quadrature_error_estimate"]
                         +row["continuum"]["quadrature_error_estimate"]
                         +gates["roundoff_allowance"])
            require(difference_from_continuum <= allowance, "Reconstructed modes-to-continuum gate")
            reported = reported_comparisons[(source, eta, cutoff)]
            close(reported["independent_mode_difference"], abs(mode_point["y"]-primary_point["y"]),
                  0.0, "Validator cross-route comparison record")
            close(reported["mode_refinement_difference"], mode_point["refinement_difference"],
                  0.0, "Validator refinement comparison record")
            close(reported["analytic_tail_bound"], bound, 0.0, "Validator tail-bound record")
            maximum_cross_route = max(maximum_cross_route, cross)
            maximum_refinement = max(maximum_refinement, refinement)
            if bound:
                maximum_bound_ratio = max(maximum_bound_ratio, difference_from_continuum/bound)
            comparisons.append({"source": source, "eta": eta, "K": cutoff,
                                "reconstructed_refinement_difference": refinement,
                                "reconstructed_primary_difference": cross,
                                "reconstructed_continuum_difference": difference_from_continuum,
                                "analytic_tail_bound": bound})
    return {
        "status": "passed", "classification": "post_run_archive_audit",
        "public_freeze_commit": PUBLIC_FREEZE,
        "provenance": {"registration_sha256": REGISTRATION_SHA256,
                       "primary_sha256": sha(primary_path), "modes_sha256": sha(modes_path),
                       "checks_sha256": sha(checks_path), "audit_script_sha256": sha(__file__)},
        "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                    "python_optimization": sys.flags.optimize},
        "counts": {"frozen_payload_pins": len(registration["files"]), "archives": 4,
                   "mode_observation_snapshots": 24, "panel_moment_checks": panel_moment_checks,
                   "reconstructed_finite_K_responses": 72, "fine_cross_route_comparisons": 36,
                   "bound_mutation_status_records": 17},
        "maximums": {"archive_integral_reconstruction_difference": peak_integral_difference,
                     "panel_moment_error": peak_moment_error,
                     "observation_amplitude_wronskian_over_epsilon": peak_amplitude_wronskian,
                     "observation_physical_wronskian_over_epsilon": peak_physical_wronskian,
                     "reconstructed_primary_difference": maximum_cross_route,
                     "reconstructed_refinement_difference": maximum_refinement,
                     "continuum_difference_over_tail_bound": maximum_bound_ratio},
        "post_run_audit_tolerances": {"archive_arithmetic_absolute": AUDIT_ARITHMETIC_ALLOWANCE,
                                      "GL_panel_moment_absolute": PANEL_MOMENT_ALLOWANCE,
                                      "high_k_coefficient_diagnostic_absolute": 1e-8},
        "run_audits": run_audits, "comparisons": comparisons,
        "limits": [
            "This is a post-run audit; its arithmetic checks are not new prospectively registered scientific gates.",
            "No source values, forcing histories, or mode solutions were newly evaluated or evolved.",
            "Archived snapshots independently verify observation Wronskians; intermediate-step maxima are recorded producer evidence, not complete archived trajectories.",
            "The 17 mutation statuses are hash-bound validator evidence; this audit does not rerun those mutations.",
            "Finite-step and finite-momentum-quadrature accuracy remain empirical refinement/cross-route estimates, not interval-certified errors.",
            "Only the omitted UV tail uses the analytically enclosed derivative budget; its bound does not certify the numerical integrals.",
            "Scope is fixed-geometry linear scalar variance: no stress response, stability, particle yield, or heating claim."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--primary", required=True, type=Path)
    parser.add_argument("--modes", required=True, type=Path)
    parser.add_argument("--checks", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), "Refusing to overwrite audit output")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = {"status": "started", "classification": "post_run_archive_audit"}
    try:
        report = inspect(args.checkpoint.resolve(), args.primary.resolve(),
                         args.modes.resolve(), args.checks.resolve())
    except Exception as exc:
        report.update({"status": "failed", "exception": str(exc), "traceback": traceback.format_exc()})
        raise
    finally:
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"status": report["status"], "output": str(args.output),
                      "counts": report["counts"], "maximums": report["maximums"]}))


if __name__ == "__main__":
    main()
