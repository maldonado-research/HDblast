"""Recognize one preserved scientific FAIL without changing its status or exit.

This check is a prerequisite to the separately registered K=384 follow-up. It
never certifies the original matrix as PASS. Any extra failure is an error.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

PROTOCOL_SHA256 = "37fda5b30332d0079ec44db708e266788fdf5859303707bb709750fdcf59d898"
PRIMARY_SHA256 = "32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a"
OBS = ("rho", "p", "Q")


def read(path):
    def reject(value):
        raise ValueError("Nonfinite JSON constant: " + value)
    data = json.loads(path.read_text(), parse_constant=reject)
    def finite(item):
        if isinstance(item, float) and not math.isfinite(item):
            raise ValueError("Nonfinite JSON number in " + str(path))
        if isinstance(item, dict):
            for value in item.values():
                finite(value)
        elif isinstance(item, list):
            for value in item:
                finite(value)
    finite(data)
    return data


def require(condition, message):
    if not condition:
        raise ValueError(message)


def gate(record, should_pass=True, atol=2e-5, rtol=.005):
    difference, scale, limit = (record[key] for key in ("absolute_difference", "signal_scale", "limit"))
    require(difference >= 0 and scale >= 0, "Negative difference/scale")
    require(math.isclose(limit, atol + rtol * scale, rel_tol=1e-12, abs_tol=1e-15), "Changed gate tolerance")
    require(record["passed"] is (difference <= limit), "Stored observable gate is inconsistent with its values")
    if should_pass is not None:
        require(record["passed"] is should_pass, "Unexpected observable gate result")


def exchange_gate(record, must_pass):
    for prefix, atol, rtol in (("ledger", 1e-5, .005), ("local", 2e-4, .01)):
        maximum, scale, limit = (record[prefix + suffix] for suffix in ("_residual_max", "_scale", "_limit"))
        require(maximum >= 0 and scale >= 0, "Negative exchange residual/scale")
        require(math.isclose(limit, atol + rtol * scale, rel_tol=1e-12, abs_tol=1e-15), "Changed exchange tolerance")
        require(record[prefix + "_passed"] is (maximum <= limit), "Stored exchange gate is inconsistent")
        if must_pass:
            require(record[prefix + "_passed"] is True, "Additional fine-grid exchange failure")


def bare_work_gate(run_root, record):
    path = run_root / record["selected"] / "arrays.npz"
    with np.load(path, allow_pickle=False) as arrays:
        energy = arrays["a"]**4 * arrays["K192_rho"]
        ledger = arrays["K192_ode_ren_ledger"]
    require(energy.shape == ledger.shape == (801,), "Bare-work curve shape changed")
    require(np.isfinite(energy).all() and np.isfinite(ledger).all(), "Nonfinite bare-work curve")
    maximum = float(np.max(abs(energy - energy[0] - ledger)))
    scale = float(max(np.max(abs(energy - energy[0])), np.max(abs(ledger))))
    limit = 1e-5 + .005 * scale
    stored = record["selected_diagnostics"]
    require(math.isclose(stored["ode_bare_work_residual_max"], maximum, rel_tol=1e-11, abs_tol=1e-15), "Bare-work maximum differs from saved curves")
    require(math.isclose(stored["ode_bare_work_limit"], limit, rel_tol=1e-12, abs_tol=1e-15), "Changed bare-work tolerance")
    require(stored["ode_bare_work_passed"] is True and maximum <= limit, "Additional bare-work failure")


def observable_gates(records, expected_failure=None):
    require(set(records) == set(OBS), "Empty or altered observable comparison")
    for name in OBS:
        gate(records[name], should_pass=(name != expected_failure))


def check(run_root, exit_code):
    require(exit_code == 1, "Original matrix must retain actual scientific-failure exit 1")
    require(not (run_root / "failure.json").exists(), "Original matrix had an execution exception")
    summary_path = run_root / "summary.json"
    summary = read(summary_path)
    require(summary["status"] == "FAIL", "Original scientific status must remain FAIL")
    require(summary["provenance"]["protocol_sha256"] == PROTOCOL_SHA256, "Original protocol changed")
    require(summary["provenance"]["source_sha256"] == PRIMARY_SHA256, "Original primary implementation changed")
    require(set(summary["amplitudes"]) == {"0.0", "0.2"}, "Original matrix incomplete")
    checked = []
    for amplitude, expected_status in (("0.0", "PASS"), ("0.2", "FAIL")):
        record = summary["amplitudes"][amplitude]
        require(record["status"] == expected_status, "Unexpected amplitude verdict")
        require(record["selected"].endswith("_quadrature_K192"), "Expected registered K192 contingency did not execute")
        observable_gates(record["final_cutoff_gate"], expected_failure="p" if amplitude == "0.2" else None)
        observable_gates(record["solver_refinement"])
        observable_gates(record["quadrature_refinement"])
        gate(record["time_refinement"], atol=1e-5, rtol=.005)
        require(record["exchange_fine"]["time_nodes"] == 801 and record["exchange_coarse"]["time_nodes"] == 401, "Exchange resolutions changed")
        exchange_gate(record["exchange_fine"], must_pass=True)
        exchange_gate(record["exchange_coarse"], must_pass=False)  # Registered reporting only.
        selected = record["selected_diagnostics"]
        require(selected["exchange_fine"] == record["exchange_fine"] and selected["exchange_coarse"] == record["exchange_coarse"], "Selected exchange reports disagree")
        bare_work_gate(run_root, record)
        gate(selected["future_spectral_energy"])
        gate(selected["trace"]["direct_derivatives"], atol=2e-6, rtol=5e-6)
        sampled = selected["trace"]["sampled_derivatives"]
        require(len(sampled) == 2 and {entry["time_nodes"] for entry in sampled} == {401, 801}, "Trace diagnostics incomplete")
        for entry in sampled:
            gate(entry, should_pass=True if entry["time_nodes"] == 801 else None, atol=2e-4, rtol=.01)
        contingency = record["contingency"]
        initial = record["initial_cutoff_K48_to_K96"]
        require(set(initial) == set(OBS), "Initial cutoff comparison empty or altered")
        for entry in initial.values():
            gate(entry, should_pass=None)
        require(any(entry["passed"] is False for entry in initial.values()), "Stored K48-to-K96 gates did not trigger extension")
        require(contingency["trigger"] == "K48-to-K96 cutoff gate failed", "Contingency trigger changed")
        observable_gates(contingency["continuity_at_K96"])
        require(contingency["cutoff_K96_to_K192"] == record["final_cutoff_gate"], "Cutoff gate disagreement")
        observable_gates(contingency["quadrature_at_K192"])
        observable_gates(contingency["solver_step_refinement_at_K192"])
        require(contingency["all_extended_wronskians_passed"] is True, "Additional extended Wronskian failure")
        checked.append(dict(amplitude=float(amplitude), scientific_status=expected_status,
                            final_cutoff_gate=record["final_cutoff_gate"],
                            reported_coarse_exchange=record["exchange_coarse"],
                            reported_coarse_trace=next(entry for entry in sampled if entry["time_nodes"] == 401)))
    diagnostics = sorted(run_root.glob("*/diagnostics.json"))
    expected_labels = {f"A{amplitude}_{label}_K{cutoff}" for amplitude in ("0", "0.2")
                       for label, cutoff in (("primary", 96), ("tight", 96), ("quadrature", 96),
                                             ("tight", 192), ("quadrature", 192), ("tail_step", 192))}
    expected_labels.add("A0_tight_K96_static")
    require({path.parent.name for path in diagnostics} == expected_labels, "Registered runs missing or unexpected")
    for path in diagnostics:
        diagnostic = read(path)
        expected_limit = 2e-9 if "_primary_" in path.parent.name else 2e-10
        require(diagnostic["wronskian_limit"] == expected_limit, "Changed Wronskian limit")
        require(0 <= diagnostic["wronskian_max"] <= expected_limit and diagnostic["wronskian_passed"] is True,
                "Additional Wronskian failure: " + path.parent.name)
    static = summary["static_negative_control"]
    require(static["passed"] is True and static["absolute_limit"] == 2e-5, "Original static control failed or gate changed")
    require(set(static["maxima"]) == set(OBS) and all(0 <= value <= 2e-5 for value in static["maxima"].values()), "Static floor outside original gate")
    return dict(validation_status="EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED", scientific_status="FAIL",
                original_command_exit_code=exit_code, original_summary_sha256=hashlib.sha256(summary_path.read_bytes()).hexdigest(),
                protocol_sha256=PROTOCOL_SHA256, primary_source_sha256=PRIMARY_SHA256,
                failed_acceptance_gate="A=0.2 final K96-to-K192 pressure cutoff change",
                all_other_acceptance_gates_passed=True, diagnostic_runs_checked=len(diagnostics), amplitudes=checked,
                limitation="This recognizes the original failure as prerequisite evidence only; a separate prospective follow-up is required and may fail.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--matrix-exit-code", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = check(args.run_root, args.matrix_exit_code)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
