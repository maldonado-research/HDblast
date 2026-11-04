#!/usr/bin/env python3
"""Portable exact-rational intrinsic preparation-class certificate.

Only Python's standard library is used. No oracle, evolved-field evaluation,
quadrature or bump evaluation is performed. Analytic ODE and mean-resolvent
premises are documented in PROOF.md; this script verifies their finite rational
consequences and authenticated retained coefficient boxes.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import argparse
import hashlib
import json
import sys


DESIGN_SHA = "18b6defb16a5f79548666ceb539cc6152c6cfeaf36b68f4fc26c0cd31081900a"
PINS = {
    "static": ("research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json",
               "3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22"),
    "base": ("research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json",
             "1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3"),
    "direction": ("research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/ORIGINAL_REPORT.json",
                  "78d63abd8caa6c183bd1d7ee7ba897be393bde14bcca1a37ee25073ea8d8b880"),
}
FOLLOW_SHA = "d9fd349165af7f7fd495b4979f1edaae49aaa75a3cadc0aaebbf8d1e07bdbb4c"
SOURCE_FROZEN_REPORT_SHA = "5b1b693b96087db54f419d232941587925eec1b68d5c41d9550526a25bbaa605"
SCIENCE_SHA = "e940b780cb5a020371575f42deef7e080a69c0e605807d44f152c30893ea69ce"
SCIENCE_EXCLUSIONS = ["authenticated_inputs", "frozen_design", "provenance", "scope"]
COMPONENT = Path(__file__).resolve().parent


class Rejected(Exception):
    pass


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise Rejected("invalid arguments")


class Checks:
    def __init__(self):
        self.entries = []

    def need(self, condition, name):
        if not condition:
            raise Rejected(name)
        self.entries.append(name)


def read_pinned(path, sha, checks, name):
    try:
        payload = path.read_bytes()
    except OSError:
        raise Rejected("source unreadable: " + name) from None
    checks.need(hashlib.sha256(payload).hexdigest() == sha, name + " byte pin")
    return json.loads(payload), {"sha256": sha, "bytes": len(payload)}


def box(node):
    return tuple(map(F, node["exact"]))


def rec(value):
    if isinstance(value, tuple):
        return {"exact": [str(x) for x in value], "display_only": [float(x) for x in value]}
    return {"exact": str(value), "display_only": float(value)}


def subbox(node, lower, upper, checks, label):
    lo, hi = box(node)
    checks.need(lower <= lo <= hi <= upper, label)


def runtime_need(ok, message):
    if not ok:
        raise Rejected(message)


def validate_output(path, repo_root):
    runtime_need(path.is_absolute(), "absolute output path")
    # Reject symlinks before resolving, including a symlinked ancestor.
    runtime_need(not any(p.is_symlink() for p in [path, *path.parents]),
                 "output path contains symlink")
    runtime_need(path == path.resolve(), "normalized output path")
    runtime_need(path.parent.is_dir(), "output parent exists")
    runtime_need(not path.exists(), "output absent")
    resolved_repo = repo_root.resolve()
    runtime_need(not path.is_relative_to(COMPONENT) and not path.is_relative_to(resolved_repo),
                 "output outside code and repository")


def run(args):
    # Runtime path guards are separate from the unchanged 38 scientific gates.
    runtime_need(args.repo_root.is_absolute(), "absolute repository path")
    validate_output(args.output, args.repo_root)
    c = Checks()
    c.need(args.output.is_absolute(), "absolute output path")
    c.need(args.output == args.output.resolve(), "normalized output path")
    c.need(args.output.parent.is_dir(), "output parent exists")
    c.need(not args.output.exists() and not args.output.is_symlink(), "output absent")
    design, design_auth = read_pinned(args.design, DESIGN_SHA, c, "frozen design")
    sources, auth = {}, {}
    for key, (relative, sha) in PINS.items():
        sources[key], auth[key] = read_pinned(args.repo_root / relative, sha, c, key)
    follow, auth["follow_on"] = read_pinned(args.follow_on, FOLLOW_SHA, c, "follow-on")
    design_auth["path"] = "FROZEN_DESIGN.json"
    for key, (relative, unused_sha) in PINS.items():
        auth[key]["path"] = relative
    auth["follow_on"]["path"] = "inputs/FOLLOW_ON_REFERENCE.json"
    st, old, direction = (sources[key] for key in ("static", "base", "direction"))

    # Constants and final gates were frozen before this script's first run.
    ell, target, delta, q = F(1, 10**6), F(5, 4), F(1, 4), F(1201, 3200000)
    rho_lo, rho_hi = F(129), F(131)
    p_lo, p_hi = F(3, 25000), F(1, 7000)
    a_lo, a_hi, py_hi, ay_hi = F(111, 1000), F(14, 125), F(1, 4000), F(1, 10000)
    Fprime_hi = F(1, 1000)
    class_radius, N_per_R_guard, N_class_guard = F(2 * 10**33), F(1, 2000), F(10**30)
    mu, exp_guard, degree = F(20001, 20000), F(4, 3), 6
    c.need(F(design["class"]["radius"]) == class_radius, "frozen class radius")
    gates = design["fixed_strict_gates"]
    c.need(F(gates["N_per_R_upper"]) == N_per_R_guard and
           F(gates["N_class_upper"]) == N_class_guard, "frozen norm gates")
    c.need(target - 1 == delta and q < 1, "target and retained contraction")
    c.need(F(st["domain"]["t"]) == 0 and F(st["domain"]["ell"]) == ell,
           "static initial slice and scale")
    c.need(F(st["domain"]["z_min"]) <= ell and
           F(st["domain"]["z_max"]) >= ell * target,
           "continuous annulus covers causal initial segment")
    c.need(st["domain"]["coverage"] == "Every point of the original packet annulus; no samples.",
           "retained continuous-domain declaration")
    c.need(st["normalization"] == "coefficient j = D_z^j(value)/j! or D_y^j(value)/j!, as named",
           "first y-jet entries equal actual first derivatives")
    subbox(st["y_jets"]["rho"][0], rho_lo, rho_hi, c, "raw static rho cap")
    subbox(st["y_jets"]["p"][0], -p_hi, -p_lo, c, "raw static nonzero negative p cap")
    subbox(st["y_jets"]["a"][0], a_lo, a_hi, c, "raw static positive a cap")
    subbox(st["y_jets"]["p"][1], -py_hi, py_hi, c, "raw static p_y cap")
    subbox(st["y_jets"]["a"][1], -ay_hi, ay_hi, c, "raw static a_y cap")
    subbox(old["first"]["rho"], rho_lo, rho_hi, c, "raw shell rho_b cap")
    subbox(old["first"]["Wprime"], -Fprime_hi, -F(1, 1200), c, "raw shell Fprime cap")
    c.need(F(old["q"]["exact"]) == q == F(direction["q"]["exact"]),
           "published coupled q agrees")
    c.need(F(direction["ell"]["exact"]) == ell and
           F(direction["target"]["exact"]) == target, "published scale and target agree")

    # The physical D-profile reconstruction gives z_r=3bz-h and
    # h_r=3ell^2 z+(-b+g)h-(2/3)gD, with b,g positive.
    b = ell * rho_hi * a_hi
    g = ell * rho_hi * p_hi**2 / (3 * a_lo)
    matrix_row_1 = 3 * b + 1
    matrix_row_2 = 3 * ell**2 + b + g
    c.need(max(matrix_row_1, matrix_row_2) <= mu, "matrix infinity norm cap")
    x = mu * delta
    # The sum through degree six plus a geometric upper bound on its tail:
    # after the first omitted term, ratios are <= x/(degree+2).
    lower_exp = sum((x**j / factorial(j) for j in range(degree + 1)), F())
    tail_first = x**(degree + 1) / factorial(degree + 1)
    tail_ratio = x / (degree + 2)
    upper_exp = lower_exp + tail_first / (1 - tail_ratio)
    c.need(0 <= tail_ratio < 1, "exponential Taylor tail ratio")
    c.need(upper_exp < exp_guard, "strict rational exponential cap")

    # All caps here are per unit R=max(sup|D|,sup|D_r|).
    H = delta * exp_guard * F(2, 3) * g
    Z = H
    Hr = 3 * ell**2 * Z + (b + g) * H + F(2, 3) * g
    k = p_hi / (2 * a_lo)
    kr = ell * rho_hi * (py_hi / (2 * a_lo) + p_hi * ay_hi / (2 * a_lo**2))
    S = k * (H + F(2, 3))
    Sr = kr * (H + F(2, 3)) + k * (Hr + F(2, 3))
    V = Hr + b * H
    D0, Dr, Dv = 2 * H + S, 2 * Hr + Sr, 2 * V
    N0, N1 = (D0 + target * Dv) / 3, Dr + Dv
    N_per_R = max(N0, N1)
    c.need(N_per_R <= N_per_R_guard, "fixed N/R guard")
    c.need(class_radius * N_per_R_guard == N_class_guard, "class radius gives fixed N guard")

    # Exact zero wedge gives |m(x)| <= (x-1) N/(1-q). Kernel integration
    # gives the complete observable differences, with every mean channel in q.
    Cv = (rho_hi * ell)**2 * Fprime_hi * delta**3 / (12 * (1 - q))
    Ct = rho_hi * ell * Fprime_hi * delta**2 / (4 * (1 - q))
    old_bounds = follow["general_homogeneous_operator_bound"]
    c.need(Cv == F(old_bounds["volume_slope_absolute_per_initial_N"]), "follow-on exact volume transfer")
    c.need(Ct == F(old_bounds["expansion_slope_absolute_per_initial_N"]), "follow-on exact expansion transfer")
    base_v, base_t = box(old["proper_v2"]), box(old["proper_Theta2"])
    headline_base_v = tuple(map(F, gates["base_v2_outward"]))
    headline_base_t = tuple(map(F, gates["base_Theta2_outward"]))
    c.need(headline_base_v[0] <= base_v[0] <= base_v[1] <= headline_base_v[1] < 0,
           "frozen outward negative base volume")
    c.need(headline_base_t[0] <= base_t[0] <= base_t[1] <= headline_base_t[1] < 0,
           "frozen outward negative base expansion")
    threshold_v, threshold_t = -base_v[1] / Cv, -base_t[1] / Ct
    threshold = min(threshold_v, threshold_t)
    c.need(threshold == F(old_bounds["sufficient_N_for_both_base_coefficients_negative_strictly_below"]),
           "exact old robustness threshold reconstructed")
    c.need(N_class_guard < threshold, "class N strictly below old exact robustness threshold")
    transfer_v, transfer_t = Cv * N_class_guard, Ct * N_class_guard
    new_v = (base_v[0] - transfer_v, base_v[1] + transfer_v)
    new_t = (base_t[0] - transfer_t, base_t[1] + transfer_t)
    out_v = tuple(map(F, gates["new_v2_outward"]))
    out_t = tuple(map(F, gates["new_Theta2_outward"]))
    c.need(out_v[0] <= new_v[0] <= new_v[1] <= out_v[1] < 0,
           "frozen outward whole-class volume sign gate")
    c.need(out_t[0] <= new_t[0] <= new_t[1] <= out_t[1] < 0,
           "frozen outward whole-class expansion sign gate")

    result = {
        "status": "PASS_CONDITIONAL_PHYSICAL_INITIAL_PROFILE_CLASS_ROBUSTNESS",
        "scope": "Portable exact rational postprocessing; retained analytic premises; no new oracle, PDE solve or bump query. Local transfer reference is a recomputed cross-check, not an analytic premise.",
        "provenance": {"source_frozen_report_sha256": SOURCE_FROZEN_REPORT_SHA,
                       "scientific_projection_sha256": SCIENCE_SHA,
                       "scientific_projection_excluded_top_level_keys": SCIENCE_EXCLUSIONS,
                       "source_frozen_design_sha256": DESIGN_SHA,
                       "local_transfer_reference_sha256": FOLLOW_SHA,
                       "source_paths": "Published source paths are repository-relative; local input labels are component-relative.",
                       "local_reference_role": "Non-premise cross-check of Cv, Ct and exact threshold; producer recomputes them from pinned old published data."},
        "frozen_design": design_auth,
        "authenticated_inputs": auth,
        "class": {"profile": design["physical_profile"], "norm": design["class"]["norm"],
                  "radius": rec(class_radius), "causal_initial_segment": ["1", "5/4"],
                  "complete_increment_required": True},
        "constants": {"ell": rec(ell), "T": rec(target), "delta": rec(delta), "q": rec(q),
                      "rho": rec((rho_lo, rho_hi)), "minus_p": rec((p_lo, p_hi)),
                      "a": rec((a_lo, a_hi)), "abs_p_y": rec(py_hi), "abs_a_y": rec(ay_hi),
                      "abs_Fprime": rec(Fprime_hi)},
        "ode_caps": {"bmax": rec(b), "gmax": rec(g), "matrix_row_1": rec(matrix_row_1),
                     "matrix_row_2": rec(matrix_row_2), "matrix_inf_norm_guard": rec(mu),
                     "Z_per_R": rec(Z), "H_per_R": rec(H), "Hr_per_R": rec(Hr)},
        "exponential_certificate": {"argument": rec(x), "Taylor_degree": degree,
                                    "partial_sum": rec(lower_exp), "first_tail_term": rec(tail_first),
                                    "tail_ratio_upper": rec(tail_ratio), "exponential_upper": rec(upper_exp),
                                    "strict_upper_guard": rec(exp_guard)},
        "conditioned_component_caps_per_R": {"k": rec(k), "abs_k_r": rec(kr),
                                              "h": rec(H), "h_r": rec(Hr), "s": rec(S),
                                              "s_r": rec(Sr), "V_L": rec(V)},
        "Gram_norm_caps_per_R": {"D0": rec(D0), "Dr": rec(Dr), "Dv": rec(Dv),
                                  "N0": rec(N0), "N1": rec(N1), "N": rec(N_per_R),
                                  "fixed_N_per_R_guard": rec(N_per_R_guard)},
        "class_N": {"actual_derived_cap": rec(class_radius * N_per_R), "fixed_guard": rec(N_class_guard)},
        "observable_operator_caps_per_N": {"volume": rec(Cv), "expansion": rec(Ct)},
        "old_exact_N_threshold": {"volume": rec(threshold_v), "expansion": rec(threshold_t),
                                  "both": rec(threshold)},
        "absolute_observable_difference_caps": {"volume": rec(transfer_v), "expansion": rec(transfer_t)},
        "base_exact_intervals": {"v2": rec(base_v), "Theta2": rec(base_t)},
        "whole_class_exact_intervals": {"v2": rec(new_v), "Theta2": rec(new_t)},
        "frozen_outward_intervals": {"v2": rec(out_v), "Theta2": rec(out_t)},
        "gates": {"N_per_R_at_most_1_over_2000": True, "N_at_most_1e30": True,
                  "N_strictly_below_old_threshold": True, "all_class_v2_strict_negative": True,
                  "all_class_Theta2_strict_negative": True},
        "verification": {"finite_checks": len(c.entries), "check_names": c.entries,
                         "oracle_calls": 0, "new_PDE_solves": 0, "new_bump_queries": 0},
        "limits": design["limits"],
    }
    scientific = {key: value for key, value in result.items() if key not in SCIENCE_EXCLUSIONS}
    scientific_bytes = json.dumps(scientific, sort_keys=True, separators=(",", ":")).encode()
    runtime_need(hashlib.sha256(scientific_bytes).hexdigest() == SCIENCE_SHA,
                 "frozen scientific projection mismatch")
    payload = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode()
    with args.output.open("xb") as stream:
        stream.write(payload)
    print(json.dumps({"status": result["status"], "finite_checks": len(c.entries),
                      "output_sha256": hashlib.sha256(payload).hexdigest()}, sort_keys=True))


def main():
    parser = Parser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--design", type=Path, default=Path(__file__).resolve().parent / "FROZEN_DESIGN.json")
    parser.add_argument("--follow-on", type=Path,
                        default=COMPONENT / "inputs/FOLLOW_ON_REFERENCE.json")
    run(parser.parse_args())


if __name__ == "__main__":
    try:
        main()
    except Rejected as exc:
        print(json.dumps({"status": "REJECTED_INCONCLUSIVE", "reason": str(exc)}, sort_keys=True))
        sys.exit(2)
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError):
        print(json.dumps({"status": "REJECTED_INCONCLUSIVE", "reason": "input or filesystem error"}, sort_keys=True))
        sys.exit(2)
