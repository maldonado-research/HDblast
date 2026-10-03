#!/usr/bin/env python3
"""Read-only exact-rational audit of saved result JSON; never import a producer."""
import argparse
import hashlib
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def pin(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def display(value):
    with localcontext() as ctx:
        ctx.prec = 20
        return format(Decimal(value.numerator) / Decimal(value.denominator), ".12E")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("checkpoint", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    names = [
        "reports/REGISTERED_RESULT.md",
        "reports/REGISTERED_RESULT_SUMMARY.json",
        "outputs/registered_repaired_run/primary/diagnostic.json",
        "outputs/registered_repaired_run/independent/diagnostic.json",
        "outputs/registered_repaired_run/validation/VALIDATION.json",
        "outputs/registered_repaired_run/validation_optimized/VALIDATION.json",
    ]
    before = {name: pin(args.checkpoint / name) for name in names}
    loaded = {name: json.loads((args.checkpoint / name).read_text())
              for name in names if name.endswith(".json")}
    summary = loaded[names[1]]
    primary = loaded[names[2]]
    independent = loaded[names[3]]
    validators = [loaded[names[4]], loaded[names[5]]]
    canonical = {}
    for record in primary["records"]:
        controls = [c for c in record["controls"] if c["quadrature_order"] == 32]
        if len(controls) != 1:
            raise ValueError("Canonical primary order must occur once")
        for row in controls[0]["levels"]["100"]:
            key = (record["source"], record["setting"], row["K"])
            if key in canonical:
                raise ValueError("Duplicate canonical case")
            canonical[key] = row
    other = {(r["source"], r["setting"], r["K"]): r["precisions"]["100"]["fields"]
             for r in independent["rows"]}
    expected = {(s, setting, k) for s in ("positive_B", "signed_uB")
                for setting in ("coarse", "fine") for k in (64, 128, 256)}
    if set(canonical) != expected or set(other) != expected:
        raise ValueError("The twelve-case universe must be complete")
    errors, cases, witnesses = [], [], []
    tolerance = Fraction("1e-12")
    q = Fraction
    maxima = {k: max(abs(q(r[k])) for r in canonical.values())
              for k in summary["maxima_exact_rational"]}
    if any(maxima[k] != q(v) for k, v in summary["maxima_exact_rational"].items()):
        errors.append("summary_maximum_mismatch")
    serial_summary = {(r["source"], r["setting"], r["K"]): r
                      for r in summary["canonical_primary_rows"]}
    if set(serial_summary) != expected:
        errors.append("summary_case_universe")
    for key in sorted(expected):
        r, d = canonical[key], other[key]
        closures = {
            "primary_sampled_accounting": q(r["D_S"]) - q(r["D_cont"]) - q(r["E_Q"]),
            "primary_matched_components": q(r["D_cont"]) - q(r["E_flow"]) - q(r["E_operator"]) - q(r["E_reconstruction"]),
            "primary_raw_components": q(r["D_cont_A"]) - q(r["E_flow"]) - q(r["E_momentum"]) - q(r["E_operator"]) - q(r["E_reconstruction"]),
            "primary_contact_matching": q(r["I_ab"]) - q(r["I_A"]) - q(r["E_momentum"]),
            "independent_sampled_accounting": q(d["D_S"]) - q(d["D_cont"]) - q(d["E_Q"]),
            "independent_components": q(d["D_cont"]) - q(d["E_flow"]) - q(d["E_operator"]) - q(d["E_evolution_ledger"]),
        }
        if any(abs(v) > tolerance for v in closures.values()):
            errors.append("serialized_accounting:" + repr(key))
        witness = all(abs(q(t["D_S"])) > q("2e-6")
                      and abs(q(t["D_cont"])) <= q("0.1") * abs(q(t["D_S"]))
                      and abs(q(t["E_Q"])) >= q("0.9") * abs(q(t["D_S"]))
                      for t in (r, d))
        if witness:
            witnesses.append({"source": key[0], "setting": key[1], "K": key[2]})
        for field, value in serial_summary.get(key, {}).items():
            if field not in ("source", "setting", "K") and q(value) != q(r[field]):
                errors.append("summary_scalar_mismatch:" + repr((key, field)))
        triangle_p = q(r["E_flow_triangle_bound"])
        triangle_i = q(d["triangle_bound"])
        if triangle_p < abs(q(r["E_flow"])) or triangle_i < abs(q(d["E_flow"])):
            errors.append("triangle_projection_mismatch:" + repr(key))
        cases.append({"source": key[0], "setting": key[1], "K": key[2],
                      "same_case_witness": witness,
                      "closure_max_exact_rational": str(max(abs(v) for v in closures.values())),
                      "matched_integral_gap_exact_rational": str(q(r["I_ab"]) - q(d["I_ab"])),
                      "flow_signed_to_computed_triangle_fraction": {
                          "primary": str(abs(q(r["E_flow"])) / triangle_p) if triangle_p else None,
                          "independent": str(abs(q(d["E_flow"])) / triangle_i) if triangle_i else None},
                      "computed_flow_triangle": {"primary": str(triangle_p), "independent": str(triangle_i)}})
    witnessed = {(r["source"], r["setting"], r["K"]) for r in witnesses}
    for validator in validators:
        if validator["classification"] != "LEDGER_ERROR_DEMONSTRATED" or validator["failures"] or validator["fatal_failures"]:
            errors.append("authoritative_validator_outcome")
        if witnessed != {(r["source"], r["setting"], r["K"]) for r in validator["attribution_cases"]}:
            errors.append("authoritative_witness_mismatch")
        if validator["old_metric_status"] != "FAIL":
            errors.append("historical_failure_changed")
        if q(validator["cross_route_max_scientific_gap_exact_rational"]) != q(summary["cross_route_max_scientific_gap_exact_rational"]):
            errors.append("validator_maximum_mismatch")
    after = {name: pin(args.checkpoint / name) for name in names}
    if before != after:
        errors.append("sealed_input_changed_during_review")
    receipt = {
        "schema_version": 1,
        "status": "SAVED_RESULT_AUDIT_PASS" if not errors else "SAVED_RESULT_AUDIT_FAIL",
        "scope": "Exact-rational review of already serialized outputs only. No producer import, source callback, checkpoint array decode, numerical rerun, or new acceptance gate.",
        "reviewed_files": before,
        "reviewed_files_unchanged_after_read": before == after,
        "classification": summary["classification"],
        "old_metric_status": summary["old_metric_status"],
        "fresh_zip_replay_status_in_sealed_report": summary["fresh_zip_replay"],
        "case_count": len(cases), "same_case_witnesses": witnesses,
        "canonical_primary_maxima": {k: {"exact_rational": str(v), "decimal": display(v)} for k, v in maxima.items()},
        "canonical_independent_max_abs_D_cont": display(max(abs(q(r["D_cont"])) for r in other.values())),
        "max_computed_flow_triangle": {route: display(max(q(c["computed_flow_triangle"][route]) for c in cases)) for route in ("primary", "independent")},
        "min_flow_signed_to_computed_triangle_fraction": {route: display(min(q(c["flow_signed_to_computed_triangle_fraction"][route]) for c in cases if c["flow_signed_to_computed_triangle_fraction"][route] is not None)) for route in ("primary", "independent")},
        "errors": errors, "cases": cases,
        "limits": ["Computed triangle expressions concern retained mode differences at fixed anchors; they are not certified total-error or incoming-state bounds.", "Canonical signed residuals and cross-route agreement do not establish all-source/all-phase accuracy, full-history contact equality, continuum momentum convergence, or coupled gravitational dynamics.", "The original full-global-prefix Simpson ledger is compared by the two registered readers; this audit does not decode or replay that history."]}
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: receipt[k] for k in ("status", "case_count", "max_computed_flow_triangle", "min_flow_signed_to_computed_triangle_fraction", "canonical_independent_max_abs_D_cont", "errors")}, indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
