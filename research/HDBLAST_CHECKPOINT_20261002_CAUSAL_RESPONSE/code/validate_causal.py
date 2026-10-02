#!/usr/bin/env python3
"""Validate the frozen primary and independent mode evidence, including mutations.

No producer functions are imported. Explicit exceptions remain active under -O.
Reported quadrature/refinement uncertainties are estimates, not certificates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import traceback

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def close(actual, expected, tolerance, label):
    require(math.isfinite(actual) and math.isfinite(expected), label+": nonfinite")
    require(abs(actual-expected) <= tolerance, label+": mismatch")


def require_negative(value, margin, label):
    require(math.isfinite(value) and value < -margin, label+": failed negative margin")


def reject_mutation(name, callable_check, records):
    try:
        callable_check()
    except RuntimeError:
        records.append({"name": name, "rejected": True})
    else:
        raise RuntimeError("Wrong-formula mutation survived: "+name)


def check_registration(pin):
    path = ROOT/"FULL_REGISTRATION.json"
    require(sha(path) == pin, "Registration hash pin mismatch")
    registration = json.loads(path.read_text())
    for name, expected in registration["files"].items():
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT), "Registration path escape")
        require(sha(path) == expected, "Registered file changed: "+name)
    for required in ("code/causal_primary.py", "code/validate_causal.py", "EXPERIMENT.json", "theory/tail_bounds.py"):
        require(required in registration["files"], "Required input absent from registration: "+required)
    return registration


def index_rows(document, experiment, label):
    required = {(s, t) for s in experiment["source"]["ids"] for t in experiment["observations"]}
    rows = {(r["source"], r["eta"]): r for r in document["rows"]}
    require(len(rows) == len(document["rows"]) and set(rows) == required, label+": row grid changed")
    for row in rows.values():
        require([v["K"] for v in row["finite_k"]] == experiment["cutoffs"], label+": cutoff grid changed")
    return rows


def check_mode_provenance(modes, mode_path, registration, experiment):
    require(modes["status"] == "passed_internal_gates", "Independent producer did not pass its gates")
    require(modes["route"] == "independent_forced_canonical_modes", "Unexpected independent route")
    require(modes["epsilon"] == experiment["source"]["epsilon"], "Independent source amplitude changed")
    require(modes["normalization"] == "y=a^2 deltaQ/epsilon", "Independent normalization changed")
    require(modes["runtime"]["numpy"] == "2.2.6", "Independent NumPy version changed")
    manifest_path = ROOT/"independent/MANIFEST.json"
    require("independent/MANIFEST.json" in registration["files"], "Independent manifest is not in full registration")
    require("independent/forced_modes.py" in registration["files"], "Independent producer is not in full registration")
    require(modes["provenance"]["manifest_sha256"] == sha(manifest_path), "Independent manifest evidence pin mismatch")
    manifest = json.loads(manifest_path.read_text())
    require(modes["provenance"]["files"] == manifest["files"], "Independent input file list changed")
    for entry in manifest["files"]:
        path = (manifest_path.parent/entry["path"]).resolve()
        require(path.is_relative_to(ROOT), "Independent input escaped frozen package")
        require(sha(path) == entry["sha256"], "Independent input hash mismatch")
    resources = modes["resources"]
    require(0 <= resources["elapsed_seconds"] <= 900, "Independent wall budget")
    require(0 <= resources["peak_rss_kib"] <= 128*1024, "Independent peak memory budget")
    runs = {(run["source"], run["setting"]): run for run in modes["runs"]}
    required = {(source, setting) for source in experiment["source"]["ids"] for setting in ("coarse", "fine")}
    require(len(modes["runs"]) == len(runs) == 4 and set(runs) == required, "Independent run coverage changed")
    directory = mode_path.resolve().parent
    summary = index_rows(modes, experiment, "independent summary")
    for (source, setting), run in runs.items():
        fine = setting == "fine"
        steps, nodes = (1152, 16384) if fine else (576, 8192)
        require(run["time_steps"] == steps and run["momentum_nodes"] == nodes, "Independent integration grid changed")
        require(run["dt"] == 1/(256 if fine else 128), "Independent time step changed")
        require(run["momentum_panel_width"] == (0.25 if fine else 0.5), "Independent momentum panels changed")
        require(run["time_gauss_nodes"] == 8 and run["local_time_momentum_pairs"] == steps*8*nodes,
                "Independent source evaluation budget changed")
        for name in ("wronskian_max_scaled", "physical_wronskian_max_scaled"):
            require(0 <= run[name] <= experiment["gates"]["mode_wronskian_over_epsilon"], "Global independent Wronskian gate")
        archive = Path(run["archive"]["path"])
        require(not archive.is_absolute(), "Independent archive must have a relative path")
        archive = (directory/archive).resolve()
        require(archive.is_relative_to(directory), "Independent archive escaped its output directory")
        require(sha(archive) == run["archive"]["sha256"], "Independent raw archive hash mismatch")
        require([row["eta"] for row in run["rows"]] == experiment["observations"], "Independent raw run observation grid changed")
        for row in run["rows"]:
            require(row["source"] == source, "Independent raw run source mismatch")
            main = summary[(source, row["eta"])]
            require([point["K"] for point in row["finite_k"]] == experiment["cutoffs"], "Independent raw run cutoff grid changed")
            close(row["a"], main["a"], 1e-15, "Independent raw/summary scale factor")
            close(row["s_over_epsilon"], main["s_over_epsilon"], 1e-15, "Independent raw/summary source")
            for point, summarized in zip(row["finite_k"], main["finite_k"]):
                close(point["y"], summarized["y" if fine else "y_coarse"], 0.0, "Independent raw/summary response")
            for name, peak in (("wronskian_scaled", "wronskian_max_scaled"),
                               ("physical_wronskian_scaled", "physical_wronskian_max_scaled")):
                require(0 <= row[name] <= run[peak], "Independent observation exceeds recorded global Wronskian maximum")
                if fine:
                    close(row[name], main[name], 0.0, "Independent raw/summary Wronskian")
    refinement_max = max(point["refinement_difference"] for row in summary.values() for point in row["finite_k"])
    close(modes["maximum_refinement_difference"], refinement_max, 0.0, "Global independent refinement maximum")


def check(primary, modes, experiment):
    gates = experiment["gates"]
    rows = index_rows(primary, experiment, "primary")
    mode_rows = index_rows(modes, experiment, "modes")
    comparisons = []
    mutations = []
    for key, row in rows.items():
        source_id, eta = key
        a = -1.0/eta
        close(row["a"], a, 1e-15, "Scale factor")
        u = eta+4.0
        f0 = 0.0 if abs(u) >= 1 else math.exp(-u*u/(1-u*u))
        if source_id == "signed_uB":
            f0 *= u
        close(row["s_over_epsilon"], f0, 1e-14, "Registered source")
        continuum = row["continuum"]
        finitepart = row["finite_part"]
        for quantity in (continuum, finitepart, *row["finite_k"]):
            error = quantity["quadrature_error_estimate"]
            require(0 <= error <= gates["primary_normalized_quadrature_estimate_max"], "Primary quadrature uncertainty gate")
        close(continuum["y"], finitepart["y"], gates["primary_memory_forms_absolute_agreement"], "Two memory representations")
        if eta < -5:
            require(continuum["y"] == finitepart["y"] == 0, "Strict retarded causality")
        if eta > -3:
            require_negative(continuum["y"], gates["positive_postpulse_negative_response_margin"] if source_id == "positive_B" else 0,
                             "Post-pulse memory sign")
        norm_A = row["tail_derivative_budget"]["normalized_A_upper"]
        require(0 <= norm_A <= (127 if source_id == "positive_B" else 121), "Analytic derivative-budget fallback ceiling")
        for point, mode in zip(row["finite_k"], mode_rows[key]["finite_k"]):
            K = point["K"]
            expected_bound = (1.5*abs(f0)*(2/eta**2)+norm_A/4)/(16*math.pi**2*K*K)
            tail = point["analytic_tail_bound"]
            close(tail, expected_bound, max(1e-20, abs(expected_bound)*2e-14), "Actual combined analytic tail")
            if K == max(experiment["cutoffs"]):
                require(tail <= gates["largest_cutoff_analytic_tail_bound_max"], "Largest-cutoff tail precision gate")
            budget = tail+point["quadrature_error_estimate"]+continuum["quadrature_error_estimate"]+gates["roundoff_allowance"]
            difference = abs(point["y"]-continuum["y"])
            require(difference <= budget, "Finite-K to continuum mismatch beyond distinct tail/numerical budgets")
            close(mode["refinement_difference"], abs(mode["y"]-mode["y_coarse"]), 1e-16, "Mode refinement accounting")
            require(mode["refinement_difference"] <= gates["mode_refinement_absolute"], "Mode refinement gate")
            close(mode["estimated_numerical_error"], 2*mode["refinement_difference"], 1e-16, "Mode numerical estimate accounting")
            cross = abs(mode["y"]-point["y"])
            if eta < -5:
                require(point["y"] == mode["y"] == mode["y_coarse"] == tail == 0.0,
                        "Finite-regulator and forced-mode strict pre-pulse zero")
            require(cross <= gates["mode_cross_route_absolute"]+point["quadrature_error_estimate"], "Independent forced-mode finite-K gate")
            require(abs(mode["y"]-continuum["y"]) <= tail+gates["mode_cross_route_absolute"]+
                    point["quadrature_error_estimate"]+continuum["quadrature_error_estimate"]+gates["roundoff_allowance"],
                    "Independent modes to removed-cutoff memory gate")
            comparisons.append({"source": source_id, "eta": eta, "K": K,
                                "memory_to_finite_k_difference": difference, "analytic_tail_bound": tail,
                                "primary_quadrature_estimate": point["quadrature_error_estimate"],
                                "independent_mode_difference": cross,
                                "mode_refinement_difference": mode["refinement_difference"]})
        for key_w in ("wronskian_scaled", "physical_wronskian_scaled"):
            require(0 <= mode_rows[key][key_w] <= gates["mode_wronskian_over_epsilon"], "Independent Wronskian gate")
    controls = primary["controls"]
    for error in (controls["advanced_quadrature_error_estimate"],
                  controls["conformal_input_mutation_at_minus2p5"]["quadrature_error_estimate"]):
        require(0 <= error <= gates["primary_normalized_quadrature_estimate_max"], "Wrong-formula diagnostic quadrature uncertainty gate")
    # The exact stationary expression is independently evaluated here.
    with mp.workdps(70):
        reference = float(-(2*mp.euler+mp.log(2))/(16*mp.pi**2))
    analytic = controls["analytic_stationary"]
    for name in ("memory_Qx", "common_action_Qx"):
        close(analytic[name], reference, gates["analytic_stationary_absolute"], "Fixed-r stationary susceptibility")
    for name in ("omit_plus_one_Qx", "omit_euler_gamma_Qx", "moving_reference_Qx"):
        require(abs(analytic[name]-reference) > 0.001, "Stationary mutation lacks declared finite margin")
        reject_mutation(name, lambda name=name: close(analytic[name], reference, 1e-14, name), mutations)
    # Separated-time normalization is exact, hence independent of any local term.
    separation = 0.1
    kernel = -1/(8*math.pi**2*separation)
    for name, wrong in (("reverse_kubo_sign", -kernel), ("omit_wick_factor_two", kernel/2), ("replace_2k_with_k", 2*kernel)):
        reject_mutation(name, lambda wrong=wrong: close(wrong, kernel, 1e-14, "Exact separated kernel"), mutations)
    post = rows[("positive_B", -2.5)]["continuum"]["y"]
    margin = gates["positive_postpulse_negative_response_margin"]
    for name, wrong in (("flipped_postpulse_sign", -post), ("instantaneous_only", 0.0), ("memory_window_excludes_pulse", 0.0)):
        reject_mutation(name, lambda wrong=wrong: require_negative(wrong, margin, "Post-pulse sign and memory"), mutations)
    advanced = controls["advanced_response_y_at_minus5p5"]
    require_negative(advanced, 0.004, "Advanced-control sensitivity")
    for name, wrong in (("advanced_kernel", advanced), ("time_symmetric_kernel", advanced/2)):
        reject_mutation(name, lambda wrong=wrong: close(wrong, 0.0, 1e-14, "Strict causality"), mutations)
    wrong_input = controls["conformal_input_mutation_at_minus2p5"]["y"]
    require(abs(wrong_input-post) > 0.004, "Missing source conformal factor lacks sensitivity")
    reject_mutation("omit_input_a_squared", lambda: close(wrong_input, post, 1e-10, "Input conformal factor"), mutations)
    wrong_output = post/(-2.5)**2
    require(abs(wrong_output-post) > 0.002, "Missing output conformal factor lacks sensitivity")
    reject_mutation("omit_output_a_inverse_squared", lambda: close(wrong_output, post, 1e-10, "Output conformal factor"), mutations)
    require([d["separation"] for d in controls["state_diagnostics"]] == [0.1, 0.4], "State diagnostic grid changed")
    for item in controls["state_diagnostics"]:
        delta = item["separation"]
        with mp.workdps(60):
            def occupation(k):
                z = (k-1)/mp.mpf(".25")
                return mp.mpf(".01")*mp.exp(1-1/(1-z*z)) if abs(z) < 1 else mp.mpf(0)
            state_reference = float(-mp.quad(lambda k: occupation(k)*mp.sin(2*k*delta), [.75, 1, 1.25])/(2*mp.pi**2))
            beta_reference = float(mp.quad(lambda k: k*occupation(k)*mp.cos(2*k*delta), [.75, 1, 1.25])/(2*mp.pi**2))
        state = item["a2_extra_occupied_kernel"]
        beta = item["a2_initial_beta_variance"]
        for error in (item["occupied_quadrature_error_estimate"], item["beta_quadrature_error_estimate"]):
            require(0 <= error <= gates["primary_normalized_quadrature_estimate_max"], "State diagnostic quadrature uncertainty gate")
        close(state, state_reference, 1e-12, "Independent occupation integral")
        close(beta, beta_reference, 1e-12, "Independent initial-beta integral")
        require_negative(state, gates["state_extra_kernel_margin"], "Changed occupied state")
        require(beta > gates["initial_beta_margin"], "Initial-state boundary term")
        reject_mutation("silently_return_vacuum_for_occupied_state_"+str(delta),
                        lambda: close(0.0, state_reference, 1e-12, "Occupation response"), mutations)
        reject_mutation("omit_initial_beta_boundary_"+str(delta),
                        lambda: close(0.0, beta_reference, 1e-12, "Initial-state boundary"), mutations)
    require(len(mutations) == 17, "Wrong-formula control count changed")
    return {"status": "passed", "comparisons": comparisons, "wrong_formula_controls": mutations,
            "analytic_stationary_Qx": reference,
            "uncertainty_scope": "Only the UV-tail component is analytically enclosed. QUADPACK and independent integration refinement remain numerical uncertainty estimates.",
            "claim_boundary": experiment["scope"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--modes", type=Path, required=True)
    parser.add_argument("--registration-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    require(not output.exists(), "Validation output already exists")
    require(not output.is_relative_to(ROOT), "Validation output must be outside the frozen checkpoint")
    registration = check_registration(args.registration_sha256)
    evidence = {"registration_sha256": args.registration_sha256,
                "primary_sha256": sha(args.primary), "modes_sha256": sha(args.modes),
                "validator_sha256": sha(__file__), "python_optimization": sys.flags.optimize}
    try:
        primary = json.loads(args.primary.read_text())
        modes = json.loads(args.modes.read_text())
        require(primary["provenance"]["registration_sha256"] == args.registration_sha256, "Primary provenance registration mismatch")
        require(primary["provenance"]["producer_sha256"] == sha(ROOT/"code/causal_primary.py"), "Primary producer pin mismatch")
        require(primary["provenance"]["experiment_sha256"] == sha(ROOT/"EXPERIMENT.json"), "Primary experiment pin mismatch")
        require(primary["provenance"]["public_freeze_commit"] == modes["freeze_commit"], "Independent freeze commit mismatch")
        experiment = json.loads((ROOT/"EXPERIMENT.json").read_text())
        for package, version in (("numpy", "2.2.6"), ("scipy", "1.15.3"), ("mpmath", "1.3.0")):
            require(primary["provenance"][package] == version, "Primary numerical dependency version changed")
        require(0 <= primary["elapsed_seconds"] <= experiment["execution"]["wall_time_budget_seconds"], "Primary wall budget")
        check_mode_provenance(modes, args.modes, registration, experiment)
        result = check(primary, modes, experiment)
        result["evidence"] = evidence
    except BaseException as exc:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps({"status": "failed", "evidence": evidence, "exception": repr(exc),
                                      "traceback": traceback.format_exc()}, indent=2)+"\n")
        raise
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
    print("Causal validation passed;", len(result["comparisons"]), "route comparisons and", len(result["wrong_formula_controls"]), "mutations")


if __name__ == "__main__":
    main()
