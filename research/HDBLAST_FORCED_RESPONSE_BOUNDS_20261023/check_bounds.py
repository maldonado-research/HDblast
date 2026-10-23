#!/usr/bin/env python3
"""Authenticate published inputs and check conditional response-bound arithmetic."""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from math import factorial
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
NUMERIC_REPORT_SHA256 = "beff99f6d7c2792cdcfc1540c8949036ca454ae36418e5ea5a80f51ec453b969"


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, "DUPLICATE_JSON_KEY:" + key)
        result[key] = value
    return result


def load_json(raw):
    return json.loads(raw, object_pairs_hook=unique_object)


def authenticate(root=ROOT):
    pins = load_json((root / "INPUT_PINS.json").read_bytes())
    need(pins["version"] == 1, "INPUT_PIN_VERSION")
    needed = {
        "inputs/ACTUAL_SOURCES_NORMAL.json",
        "inputs/FULL_SECOND_RESPONSE_CONSTRUCTION.md",
        "inputs/REFINED_RESPONSE_NORMAL.json",
        "inputs/angular_generic.py",
        "inputs/angular_project.py",
        "inputs/EXACT_GEOMETRIC_PROJECTION.md",
        "inputs/UPSTREAM_MANIFEST.json",
        "inputs/FORCED_INCREMENT_BOUND.md",
    }
    need(set(pins["files"]) == needed, "INPUT_PIN_FILE_SET")
    for name, row in pins["files"].items():
        path = root / name
        need(path.is_file() and not path.is_symlink(), "INPUT_NOT_REGULAR:" + name)
        raw = path.read_bytes()
        need(len(raw) == row["bytes"] and digest(raw) == row["sha256"],
             "INPUT_AUTHENTICATION:" + name)
    upstream = load_json((root / "inputs/UPSTREAM_MANIFEST.json").read_bytes())
    need(upstream["capsule"] == "HDBLAST_GEOMETRIC_BULK_SOURCES_20261016",
         "UPSTREAM_CAPSULE")
    need(digest((root / "inputs/UPSTREAM_MANIFEST.json").read_bytes()) ==
         "3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2",
         "UPSTREAM_MANIFEST_EXTERNAL_PIN")
    indexed = {row["path"]: row for row in upstream["files"]}
    need(len(indexed) == len(upstream["files"]), "UPSTREAM_DUPLICATE_PATH")
    for name, row in pins["files"].items():
        original = row.get("upstream_manifest_path")
        if original:
            need(original in indexed, "UPSTREAM_MEMBER_MISSING:" + name)
            expected = indexed[original]
            need(row["bytes"] == expected["bytes"] and row["sha256"] == expected["sha256"],
                 "UPSTREAM_MEMBER_PIN:" + name)
    return pins


def validate_conventions(conventions):
    expected = {
        "principal_operator": "partial_x^2-partial_r^2",
        "retarded_K_factor": "1/2",
        "metric_forcing_factor": "-2 ell^2 rho^2 QEred",
        "derivative_weight": "1",
        "value_weight": "3",
        "mean_H_dual_norm": "1",
        "lower_transfer": "1+q",
        "cutoff_ratio": "h/4",
        "full_coupled_Robin_retained": True,
    }
    need(conventions == expected, "REGISTERED_CONVENTION_MISMATCH")


def sphere_moment(exponent):
    if any(power % 2 for power in exponent):
        return Q(0)
    halves = [power // 2 for power in exponent]
    numerator = 1
    for power in halves:
        numerator *= factorial(2 * power) // factorial(power)
    return Q(numerator, 4 ** sum(halves) * factorial(1 + sum(halves)))


def basis_means(root):
    # The reference module is loaded only after its complete bytes authenticate.
    # This imports the retained basis definitions; it does not run a projection
    # producer or a source evaluator. Moments below are computed independently.
    path = root / "inputs/angular_generic.py"
    spec = importlib.util.spec_from_file_location("retained_angular_basis", path)
    need(spec is not None and spec.loader is not None, "BASIS_LOADER")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    need(module.SCALAR_NAMES == ["1", "Z", "F"], "SCALAR_BASIS_LABELS")
    means = [sum((Q(value) * sphere_moment(exponent)
                  for exponent, value in tensor[()].items()), Q(0))
             for tensor in module.SCALARS]
    need(means == [Q(1), Q(0), Q(0)], "SCALAR_BASIS_MEANS")
    tree = ast.parse((root / "inputs/angular_project.py").read_text())
    assignments = [node for node in tree.body if isinstance(node, ast.Assign)
                   and any(isinstance(t, ast.Name) and t.id == "FIELD_LABELS"
                           for t in node.targets)]
    need(len(assignments) == 1, "PACKING_LABEL_ASSIGNMENT")
    fields = ast.literal_eval(assignments[0].value)
    need(fields["yy"] == ("K:0", "K:2", "K:4"), "YY_PACKING_CONVENTION")
    need(sum(len(labels) for labels in fields.values()) == 23, "FULL_23_SLOT_COUNT")
    return means


def numerical_crosscheck(root, conventions, report_path=None):
    saved = [root/"compute"/name for name in
             ("ACTUAL_CONE_RESPONSE_NORMAL.json", "ACTUAL_CONE_RESPONSE_OPTIMIZED.json")]
    for path in saved:
        need(path.is_file() and not path.is_symlink() and digest(path.read_bytes()) == NUMERIC_REPORT_SHA256,
             "ACCEPTED_NUMERICAL_REPORT_PIN")
    need(saved[0].read_bytes() == saved[1].read_bytes(), "ACCEPTED_NUMERICAL_MODE_EQUALITY")
    chosen = saved[0] if report_path is None else Path(report_path)
    raw = chosen.read_bytes()
    need(digest(raw) == NUMERIC_REPORT_SHA256, "NUMERICAL_CROSSCHECK_REPORT_PIN")
    report = load_json(raw)
    need(report["status"] == "CERTIFIED_CAUSAL_RESTRICTED_ACTUAL_FORCED_MEAN_INTERVAL" and
         report["strict_positive"] is True, "ACCEPTED_NUMERICAL_CERTIFICATE_STATUS")
    need(report["source_authentication"] == {
        "manifest_sha256": "3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2",
        "files_authenticated": 169}, "NUMERICAL_SOURCE_AUTHENTICATION")
    lower, upper = map(Q, report["actual_W_yy_probability_mean_interval"])
    need(0 < lower < upper, "POSITIVE_NUMERICAL_INTERVAL")
    value_weight = Q(conventions["value_weight"])
    need(value_weight == 3, "POINT_TO_GLOBAL_NORM_FACTOR_THREE")
    global_lower = lower/value_weight
    need(global_lower > 8800000, "STRONGER_GLOBAL_NORM_LOWER")
    rectangle = report["continuous_rectangle_corollary"]
    rect_lower, rect_upper = map(Q, rectangle["actual_W_yy_probability_mean_interval"])
    need(rectangle["strict_positive"] and rect_lower > 17000000 and rect_upper < 65000000000,
         "REGISTERED_RECTANGLE_ROUNDING")
    return {"accepted_numeric_report_sha256": NUMERIC_REPORT_SHA256,
            "point_mean_interval": [str(lower), str(upper)],
            "point_to_global_X_factor": str(value_weight),
            "global_X_nonstrict_lower": str(global_lower),
            "global_X_strict_lower": "8800000",
            "logical_transfer": "abs(mean W_yy at apex)<=3 ||W||_X",
            "rectangle_mean_interval": [str(rect_lower), str(rect_upper)],
            "scope": "Arithmetic bridge from the separately certified coordinate mean; this checker does not recalculate its source cells or couplings."}


def calculate(root=ROOT, conventions=None, numeric_report=None):
    pins = authenticate(root)
    actual = load_json((root / "inputs/ACTUAL_SOURCES_NORMAL.json").read_bytes())
    refined = load_json((root / "inputs/REFINED_RESPONSE_NORMAL.json").read_bytes())
    conventions = load_json((root / "CONVENTIONS.json").read_bytes()) if conventions is None else conventions
    validate_conventions(conventions)
    need(actual["normalization"] ==
         "Second parameter derivatives along the affine actual-first-pair path; metric_forcing=-2 ell^2 rho^2 QEred; scalar_forcing=ell^2 rho^2 Qphi. QC is the actual-metric-lowered quadratic gauge Hessian.",
         "ACTUAL_FORCING_NORMALIZATION")
    bulk = actual["bulk_witness"]
    x0, x1 = map(Q, bulk["x"])
    r0, r1 = map(Q, bulk["r"])
    h = x1 - x0
    need(h == r1-r0 == Q(1, 4096), "WITNESS_SQUARE")
    need(0 < x0 < x1 < 3 and 0 < r0 < r1 and x1+r1 < 8 and r0 > x1,
         "WITNESS_INTERIOR_EXTERIOR_LOCATION")
    cminus, cplus = map(Q, bulk["sources"]["metric_forcing"]["yy"][0])
    c = Q(16000000000)
    q = Q(refined["retained_join"]["q"])
    FH = Q(actual["global_cover"]["summary"]["positive_bundle_bulk_forcing_H_bound"]["F_H_upper"])
    GH = Q(refined["retained_join"]["G_H"])
    need(q == Q(2109, 1000000) and FH == 17387513156573 and GH == 3851400,
         "REGISTERED_OPERATOR_AND_SOURCE_CAPS")
    need(0 < q < 1 and 0 < c < cminus <= cplus, "SOURCE_SLACK_AND_CONTRACTION")
    means = basis_means(root)
    lower = cminus*h/(4*(1+q))
    strict = c*h/(4*(1+q))
    need(strict == Q(976562500000, 1002109) > 974000 and lower > strict,
         "LOWER_BOUND_ARITHMETIC")
    Scap = 3*FH + GH
    upper = Scap/(1-q)
    need(Scap == 52162543321119 and upper == Q(52162543321119000000, 997891)
         and upper < Q(5228, 1000)*10**13, "UPPER_BOUND_ARITHMETIC")
    tails = [q**n*upper for n in range(1, 5)]
    need(tails[3] == Q(1031965259713929363466110159, 997891000000000000000000)
         and tails[3] < 1035, "FOUR_TERM_TAIL")
    numerical = numerical_crosscheck(root, conventions, numeric_report)
    need(Q(numerical["global_X_nonstrict_lower"]) > lower, "NUMERICAL_NORM_BOUND_STRENGTHENS_WEAK_BOUND")
    return {
        "status": "PASS_EXACT_CONDITIONAL_RESPONSE_BOUNDS",
        "scope": "Authenticated published premises, basis convention, weak/upper/tail arithmetic and point-to-global norm bridge from the separately certified numerical report; no upstream scientific proof replay.",
        "input_files_authenticated": len(pins["files"]),
        "input_pins_sha256": digest((root/"INPUT_PINS.json").read_bytes()),
        "proof_sha256": digest((root/"CONDITIONAL_RESPONSE_BOUNDS.md").read_bytes()),
        "producer_sha256": digest((root/"check_bounds.py").read_bytes()),
        "conventions_sha256": digest((root/"CONVENTIONS.json").read_bytes()),
        "basis_probability_means": [str(value) for value in means],
        "yy_constant_slot": "K:0",
        "numerical_mean_to_global_norm": numerical,
        "exact": {
            "witness_x": [str(x0), str(x1)], "witness_r": [str(r0), str(r1)],
            "h": str(h), "minimum_r_minus_x": str(r0-x1),
            "cminus": str(cminus), "cplus": str(cplus),
            "c_rounded": str(c), "positive_source_slack": str(cminus-c),
            "q": str(q), "F_H_upper": str(FH), "full_G_H_upper": str(GH),
            "cutoff_limiting_mass": str(h*h),
            "cutoff_limiting_summed_first_variation": str(4*h),
            "S_nonstrict_lower": str(cminus*h/4),
            "S_strict_lower": str(c*h/4),
            "W_nonstrict_lower": str(lower), "W_strict_lower": str(strict),
            "S_X_upper": str(Scap), "W_X_upper": str(upper),
            "Neumann_tail_X_upper": {str(n): str(tails[n-1]) for n in range(1, 5)},
        },
        "analytic_claims_not_reproved": ["actual source enclosure", "causal inverse and coupled Robin construction", "weak test-function argument", "full all-point source projection"],
        "claims_not_made": ["full evaluated waveform", "complete numerical residual on the full domain", "physical force", "separate physical constraint solution", "nonlinear continuation", "cosmological-origin evidence", "novelty"],
    }


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         "USE_PYTHON_I_S_B")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--numeric-report", type=Path)
    args = parser.parse_args()
    raw = json.dumps(calculate(numeric_report=args.numeric_report), indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        need(not output.is_relative_to(ROOT), "OUTPUT_INSIDE_SOURCE")
        with output.open("x") as handle:
            handle.write(raw)
    else:
        sys.stdout.write(raw)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit("RESPONSE_BOUND_CHECK_FAILED:" + str(error))
