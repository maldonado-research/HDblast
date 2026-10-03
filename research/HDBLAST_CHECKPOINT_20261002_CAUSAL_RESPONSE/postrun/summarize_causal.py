#!/usr/bin/env python3
"""Post-run presentation of frozen raw data; no source or response evaluation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


FREEZE = "a01d17e015f070be2ad2018745e877148fe2f4e4"
REGISTRATION_SHA = "f59abda221a6a0296e839bc3022298bdf81fac913feb3877ac8f1485f4d71da6"


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--modes", type=Path, required=True)
    parser.add_argument("--checks", type=Path, required=True)
    parser.add_argument("--checks-optimized", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.checks_optimized is None:
        args.checks_optimized = args.checks.parent/"CHECKS_OPTIMIZED.json"
    require(not args.output.exists(), "Refusing to overwrite post-run summary output")
    require(sha(args.checkpoint/"FULL_REGISTRATION.json") == REGISTRATION_SHA, "Registration digest changed")
    registration = load(args.checkpoint/"FULL_REGISTRATION.json")
    for name, digest in registration["files"].items():
        require(sha(args.checkpoint/name) == digest, "Frozen input changed: "+name)
    primary, modes, checks, optimized = map(load, (args.primary, args.modes, args.checks, args.checks_optimized))
    require(primary["provenance"]["public_freeze_commit"] == modes["freeze_commit"] == FREEZE, "Freeze mismatch")
    require(primary["provenance"]["registration_sha256"] == REGISTRATION_SHA, "Primary registration mismatch")
    for label, report, optimization in (("normal", checks, 0), ("optimized", optimized, 1)):
        require(report["status"] == "passed", label+" validator failed")
        require(report["evidence"]["primary_sha256"] == sha(args.primary), label+" primary raw digest mismatch")
        require(report["evidence"]["modes_sha256"] == sha(args.modes), label+" modes raw digest mismatch")
        require(report["evidence"]["python_optimization"] == optimization, "Optimization provenance mismatch")
        require(len(report["comparisons"]) == 36 and len(report["wrong_formula_controls"]) == 17, "Gate coverage changed")
        require(all(c["rejected"] for c in report["wrong_formula_controls"]), "Unrejected wrong formula")
    require({k:v for k,v in checks.items() if k != "evidence"} == {k:v for k,v in optimized.items() if k != "evidence"},
            "Normal and optimized validation contents differ")
    comparisons = checks["comparisons"]
    mode_index = {(r["source"], r["eta"]): r for r in modes["rows"]}
    rows = []
    for row in primary["rows"]:
        mode = mode_index[(row["source"], row["eta"])]
        points = []
        for finite, independent in zip(row["finite_k"], mode["finite_k"]):
            require(finite["K"] == independent["K"], "Cutoff mismatch")
            points.append({"K": finite["K"], "primary_y": finite["y"], "mode_y": independent["y"],
                           "primary_to_continuum_difference": abs(finite["y"]-row["continuum"]["y"]),
                           "mode_to_continuum_difference": abs(independent["y"]-row["continuum"]["y"]),
                           "mode_to_primary_difference": abs(independent["y"]-finite["y"]),
                           "mode_refinement_difference": independent["refinement_difference"],
                           "analytic_tail_bound": finite["analytic_tail_bound"],
                           "quadrature_error_estimate": finite["quadrature_error_estimate"]})
        rows.append({"source": row["source"], "eta": row["eta"], "s_over_epsilon": row["s_over_epsilon"],
                     "continuum_y": row["continuum"]["y"],
                     "continuum_quadrature_error_estimate": row["continuum"]["quadrature_error_estimate"],
                     "finite_part_y": row["finite_part"]["y"], "finite_k": points})
    points = [p for r in rows for p in r["finite_k"]]
    last_points = [r["finite_k"][-1] for r in rows]
    nonzero_tail = [p["analytic_tail_bound"] for p in last_points if p["analytic_tail_bound"] > 0]
    metrics = {
        "registered_source_observation_points": len(rows), "finite_k_comparisons": len(comparisons),
        "rejected_wrong_formula_controls": len(checks["wrong_formula_controls"]),
        "maximum_primary_finite_k_to_continuum_difference": max(p["primary_to_continuum_difference"] for p in points),
        "maximum_primary_K256_to_continuum_difference": max(p["primary_to_continuum_difference"] for p in last_points),
        "maximum_mode_to_primary_finite_k_difference": max(p["mode_to_primary_difference"] for p in points),
        "maximum_mode_refinement_difference": max(p["mode_refinement_difference"] for p in points),
        "maximum_memory_form_difference": max(abs(r["continuum"]["y"]-r["finite_part"]["y"]) for r in primary["rows"]),
        "maximum_primary_response_quadrature_estimate": max(q["quadrature_error_estimate"] for r in primary["rows"] for q in (r["continuum"],r["finite_part"],*r["finite_k"])),
        "maximum_primary_finite_k_quadrature_estimate": max(p["quadrature_error_estimate"] for p in points),
        "minimum_nonzero_K256_analytic_tail_bound": min(nonzero_tail),
        "maximum_K256_analytic_tail_bound": max(nonzero_tail),
        "maximum_linear_wronskian_residual_over_epsilon": max(r["wronskian_max_scaled"] for r in modes["runs"]),
        "maximum_reconstructed_wronskian_residual_over_epsilon": max(r["physical_wronskian_max_scaled"] for r in modes["runs"]),
        "analytic_stationary_Qx": checks["analytic_stationary_Qx"],
        "primary_elapsed_seconds": primary["elapsed_seconds"],
        "mode_elapsed_seconds": modes["resources"]["elapsed_seconds"],
        "mode_peak_rss_kib": modes["resources"]["peak_rss_kib"],
        "all_registered_prepulse_values_exactly_zero": all(r["continuum_y"] == 0 and all(p["primary_y"] == p["mode_y"] == 0 for p in r["finite_k"]) for r in rows if r["eta"] < -5),
        "all_registered_postpulse_values_strictly_negative": all(r["continuum_y"] < 0 for r in rows if r["eta"] > -3),
    }
    hashes = {"FULL_REGISTRATION.json": REGISTRATION_SHA,
              "outputs/primary/results.json": sha(args.primary), "outputs/independent/results.json": sha(args.modes),
              "outputs/CHECKS.json": sha(args.checks), "outputs/CHECKS_OPTIMIZED.json": sha(args.checks_optimized),
              "postrun/summarize_causal.py": sha(__file__)}
    for run in modes["runs"]:
        archive = args.modes.parent/run["archive"]["path"]
        require(sha(archive) == run["archive"]["sha256"], "Mode archive digest changed")
        hashes["outputs/independent/"+archive.name] = sha(archive)
    for relative in ("started.json", "EXECUTION.json", "partial_results.json"):
        hashes["outputs/primary/"+relative] = sha(args.primary.parent/relative)
    summary = {"status": "registered_calibration_passed", "public_freeze_commit": FREEZE,
               "normalization": "y=a(eta)^2 delta_Q/epsilon; epsilon=1e-4; a=-1/eta",
               "postrun_only": True, "new_source_response_or_mode_evaluations": 0,
               "metrics": metrics, "registered_rows": rows, "raw_and_presentation_sha256": hashes,
               "frozen_source_sha256": registration["files"],
               "controls": primary["controls"], "uncertainty_scope": checks["uncertainty_scope"],
               "claim_boundary": checks["claim_boundary"]}
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output/"SUMMARY.json").write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)+"\n")
    m = metrics
    text = ["# Registered fixed-geometry causal-response calibration", "",
            "The bounded calibration passed all 36 finite-cutoff comparisons and all 17 wrong-formula controls, both normally and under Python `-O`. The independent forced modes agree with the exact finite-cutoff expression at the twelve registered source/observation points. Both pulses leave a negative response at the two registered observations after the source has vanished.", "",
            "These are responses of the scalar variance on a prescribed de Sitter geometry. The test establishes one matched retarded susceptibility; it does not establish stress response, shell stability, relaxation, particle production, radiation transfer, thermalization, or heating.", "",
            "## Frozen experiment and normalization", "",
            "H=1, r=2, a(eta)=-1/eta, epsilon=1e-4. With u=eta+4 and B=exp(1-1/(1-u²)) on |u|<1, the two sources are s=epsilon B and s=epsilon uB. Both vanish outside -5<eta<-3. The incoming BD state and geometry are fixed. The observation grid is {-5.5,-4.5,-4,-3.5,-2.5,-1.5}; cutoffs are K={64,128,256}.", "",
            "Every response and response discrepancy below uses **y=a² deltaQ/epsilon**. To recover the physical response, deltaQ=epsilon y/a². The small changes between finite K and removed cutoff are regulator errors, not changes in geometry or source amplitude.", "",
            "## Numerical evidence", "",
            "| Quantity | Recorded value | Interpretation |", "|---|---:|---|",
            f"| Maximum modes / finite-K expression difference | {m['maximum_mode_to_primary_finite_k_difference']:.9e} | Cross-route numerical agreement |",
            f"| Maximum coarse / fine mode difference | {m['maximum_mode_refinement_difference']:.9e} | Empirical refinement estimate |",
            f"| Maximum two-memory-form difference | {m['maximum_memory_form_difference']:.9e} | Internal representation check |",
            f"| Maximum finite-K / removed-cutoff difference, all K | {m['maximum_primary_finite_k_to_continuum_difference']:.9e} | Actual regulator discrepancy |",
            f"| Maximum K=256 / removed-cutoff difference | {m['maximum_primary_K256_to_continuum_difference']:.9e} | Actual regulator discrepancy |",
            f"| K=256 nonzero analytic tail bounds | {m['minimum_nonzero_K256_analytic_tail_bound']:.9e} to {m['maximum_K256_analytic_tail_bound']:.9e} | Constructive UV remainder enclosure |",
            f"| Maximum primary quadrature estimate | {m['maximum_primary_response_quadrature_estimate']:.9e} | QUADPACK estimate, not a proof |",
            f"| Maximum linear Wronskian residual / epsilon | {m['maximum_linear_wronskian_residual_over_epsilon']:.9e} | All independent mode steps/nodes |",
            f"| Maximum reconstructed Wronskian residual / epsilon | {m['maximum_reconstructed_wronskian_residual_over_epsilon']:.9e} | Observation-time canonical modes |", "",
            "The analytic tail bound encloses only the omitted combined momentum integral. Quadrature estimates, mode refinement differences, and cross-route differences remain numerical evidence; their small observed values do not certify all displayed digits or make the entire response a rigorous interval. Before support the response and tail are exactly zero. Post-pulse regulator differences are much smaller than their conservative analytic bounds.", "",
            "## All registered observations", "",
            "| Source | eta | s/epsilon | Removed-cutoff y | abs(y256-y) | Analytic K=256 tail |", "|---|---:|---:|---:|---:|---:|"]
    for row in rows:
        p = row["finite_k"][-1]
        text.append(f"| {row['source']} | {row['eta']:g} | {row['s_over_epsilon']:.9g} | {row['continuum_y']:.12g} | {p['primary_to_continuum_difference']:.5e} | {p['analytic_tail_bound']:.5e} |")
    text += ["", "The zero source at eta=-2.5 and -1.5 coexists with a negative variance response for both histories. For the positive pulse this is the signed retarded-memory test. For the signed pulse, pairing u and -u makes the post-pulse convolution positive before the common overall minus sign. The finite-part local contact vanishes at these post-pulse observations; it cannot account for their nonzero response.", "",
             "## References and controls", "",
             f"The independent analytic stationary identity gives Q_x=-(2 gamma_E+ln2)/(16 pi²)={m['analytic_stationary_Qx']:.15g}. This is a separate past-infinite-source identity, not a substitution of an instantaneous equilibrium response for the pulse histories. Normalization, matching contacts, fixed reference, strict causality, both conformal factors, occupation-state sensitivity, and the initial beta boundary term were tested. The frozen symbolic suite supplies the additional exact mass-law and coordinate-rescaling identities.", "",
             "The deliberately occupied state and initial-beta diagnostics are separate controls; they do not alter the BD state of the registered primary experiment. Every numerical wrong-formula control was rejected under normal and optimized Python execution.", "",
             "## Figures and reproducibility", "",
             "`causal_memory.svg/.png/.pdf` shows only the six registered times per source. Source markers are raw recorded values; shaded support comes from registration. Response lines guide the eye and are not additional evaluated histories. `causal_agreement.svg/.png/.pdf` compares the registered K=256 results with the removed-cutoff memory response and its analytic UV bound. Exact pre-pulse zeros are omitted on the logarithmic axes. `causal_cutoff.svg/.png/.pdf` reduces the existing points to maxima at each of the three registered cutoffs.", "",
             "The figures and this report are explicitly post-run presentation. No additional source, response, mode, or denser time-grid evaluations were performed. `SUMMARY.json` contains full-precision values, source hashes, raw archive hashes, and the registration digest. The renderer has separate presentation provenance and does not modify any frozen scientific input.", "",
             f"Public prospective freeze: `{FREEZE}`. Full registration SHA256: `{REGISTRATION_SHA}`. Primary runtime: {m['primary_elapsed_seconds']:.6f} s. Independent mode runtime: {m['mode_elapsed_seconds']:.6f} s; peak RSS: {m['mode_peak_rss_kib']} KiB. These timings are execution metadata, not portable performance guarantees.", "",
             "## Exact raw and presentation hashes", "", "| Artifact | SHA256 |", "|---|---|"]
    text += [f"| {name} | `{digest}` |" for name,digest in hashes.items()]
    text += ["", "All frozen source hashes are preserved in `SUMMARY.json` under `frozen_source_sha256` and in `FULL_REGISTRATION.json`.", ""]
    (args.output/"NUMERICAL_RESULTS.md").write_text("\n".join(text))
    print(json.dumps(metrics, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
