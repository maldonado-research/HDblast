#!/usr/bin/env python3
"""Audit manuscript displays with Fraction; no producer imports or oracle calls."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re

HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--witness", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    tex = HERE / "INTRINSIC_PREPARATION_ROBUSTNESS.tex"
    source = tex.read_text()
    report = json.loads(args.report.read_text())
    witness = json.loads(args.witness.read_text())
    checks = []

    def need(condition, name):
        if not condition:
            raise ValueError(name)
        checks.append(name)

    def scalar(path):
        out = report
        for part in path.split("."):
            out = out[part]
        return F(out["exact"])

    def interval(path):
        out = report
        for part in path.split("."):
            out = out[part]
        return tuple(map(F, out["exact"]))

    def printed_fraction(numerator, denominator, expected, name):
        expression = re.compile(r"\\frac\{\s*" + numerator + r"\s*\}\s*\{\s*" + denominator + r"\s*\}")
        need(expression.search(source) is not None or numerator + "/" + denominator in source,
             name + " present in source")
        need(F(int(numerator), int(denominator)) == expected, name + " agrees with exact report")

    def sci(s):
        s = s.replace(" ", "")
        coefficient, power = s.split(r"\times10^")
        return F(coefficient) * F(10) ** int(power.strip("{}"))

    for quantity in ("v_2", r"\Theta_2"):
        match = re.search(re.escape(quantity) + r"\\in\[([^\]]+)\]", source)
        need(match is not None, quantity + " displayed interval found")
        printed = tuple(map(sci, match.group(1).split(",")))
        key = "v2" if quantity == "v_2" else "Theta2"
        exact = interval("whole_class_exact_intervals." + key)
        frozen = interval("frozen_outward_intervals." + key)
        need(printed == frozen, key + " printed interval equals frozen outward interval")
        need(printed[0] <= exact[0] <= exact[1] <= printed[1] < 0,
             key + " exact whole-class interval contained and strictly negative")

    printed_fraction("798916710779559466175931203077",
                     "1861847711577000000000000000000000",
                     scalar("Gram_norm_caps_per_R.N"), "c_N")
    printed_fraction("17161", "767711760000000000",
                     scalar("observable_operator_caps_per_N.volume"), "C_v")
    printed_fraction("131", "63975980000",
                     scalar("observable_operator_caps_per_N.expansion"), "C_Theta")
    printed_fraction("131", "73426500000000",
                     scalar("ode_caps.H_per_R"), "H")
    printed_fraction("20001", "20000",
                     scalar("ode_caps.matrix_inf_norm_guard"), "matrix norm guard")
    printed_fraction("1201", "3200000", scalar("constants.q"), "q")
    need(scalar("ode_caps.bmax") == F(917, 62500000) and
         "b_+=917/62500000" in source, "printed b cap")
    need(scalar("ode_caps.gmax") == F(131, 16317000000000) and
         "g_+=131/16317000000000" in source, "printed g cap")
    need(scalar("conditioned_component_caps_per_R.k") == F(1, 1554)
         and r"\frac1{1554}" in source, "printed k cap")
    need(scalar("Gram_norm_caps_per_R.N") < F("0.000429098849392") < F(1, 2000),
         "printed c_N strict upper guards")
    need(scalar("class.radius") == 2 * 10**33 and
         scalar("class_N.fixed_guard") == 10**30 and
         scalar("class.radius") / 2000 == scalar("class_N.fixed_guard"),
         "printed class radius gives printed N guard")
    need(F("1.1864e30") < scalar("old_exact_N_threshold.both") and
         10**30 < F("1.1864e30"), "printed sufficient strict threshold")
    need(scalar("absolute_observable_difference_caps.volume") < F("2.23535e16"),
         "printed absolute volume increment cap")
    need(scalar("absolute_observable_difference_caps.expansion") < F("2.04765e21"),
         "printed absolute expansion increment cap")
    for key, display in [("v2", ("-5.4520e16", "-2.6521e16")),
                         ("Theta2", ("-9.0208e21", "-4.4458e21"))]:
        lo, hi = interval("base_exact_intervals." + key)
        need(F(display[0]) <= lo <= hi <= F(display[1]) < 0,
             "printed base " + key + " outward enclosure")
    caps = {
        "rho": ("129", "131"), "minus_p": ("3/25000", "1/7000"),
        "a": ("111/1000", "14/125")
    }
    for key, value in caps.items():
        need(interval("constants." + key) == tuple(map(F, value)),
             "printed continuous " + key + " caps")
    for key, value in [("abs_p_y", "1/4000"), ("abs_a_y", "1/10000"),
                       ("abs_Fprime", "1/1000"), ("T", "5/4"),
                       ("delta", "1/4"), ("ell", "1/1000000")]:
        need(scalar("constants." + key) == F(value), "printed " + key)
    d_lo, d_hi = map(F, witness["direction_D"]["exact"])
    dr_lo, dr_hi = map(F, witness["direction_D_r"]["exact"])
    need(F("3.7621e12") < d_lo <= d_hi, "printed old-witness D strict lower bound")
    need(dr_lo <= dr_hi < F("-6.870e13"), "printed old-witness derivative strict upper bound")
    need(F("3.7621e12") * 10**21 == F("3.7621e33") > 2 * 10**33,
         "printed scaled witness exceeds class radius")
    need(F("6.870e13") * 10**21 == F("6.870e34"), "printed scaled witness derivative")
    need(report["verification"]["finite_checks"] == 38 and
         "verifies 38 finite checks" in source, "printed finite-check count")
    need(all(report["verification"][x] == 0 for x in
             ["oracle_calls", "new_PDE_solves", "new_bump_queries"]),
         "new class arithmetic uses zero science oracle calls")
    result = {
        "status": "PASS_DISPLAYED_BOUNDS",
        "method": "stdlib Fraction; no producer imports; no oracle calls",
        "inputs": {name: {"filename": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                   for name, path in [("source", tex), ("report", args.report),
                                      ("witness", args.witness)]},
        "finite_checks": len(checks), "checks": checks
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "finite_checks": len(checks)}, sort_keys=True))


if __name__ == "__main__":
    main()
