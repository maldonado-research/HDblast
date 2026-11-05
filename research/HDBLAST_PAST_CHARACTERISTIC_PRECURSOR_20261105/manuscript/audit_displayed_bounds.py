#!/usr/bin/env python3
"""Exact display audit for the local past characteristic precursor manuscript."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json
import re

HERE = Path(__file__).resolve().parent
STEM = "PAST_CHARACTERISTIC_PRECURSOR"
REPORT_SHA = "2e7f708932ebfbfa5dad30d589fad1db35f1a71544a1f4b4ce42adc31a006913"
PROOF_SHA = "469d2fcc3872aeea001271d01dc03c62e71bc84e604eeb9ab3f7cb66d28596ac"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--optimized", required=True, type=Path)
    ap.add_argument("--proof", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    tex = HERE / (STEM + ".tex")
    source = tex.read_text()
    report = json.loads(args.report.read_text())
    checks = []

    def need(ok, name):
        if not ok:
            raise ValueError(name)
        checks.append(name)

    def val(path):
        obj = report
        for part in path.split("."):
            obj = obj[part]
        return F(obj)

    def interval(path):
        obj = report
        for part in path.split("."):
            obj = obj[part]
        return tuple(map(F, obj))

    def printed_fraction(n, d, expected, name):
        pattern = r"\\frac\{\s*" + str(n) + r"\s*\}\{\s*" + str(d) + r"\s*\}"
        need(re.search(pattern, source) is not None or str(n)+"/"+str(d) in source,
             name + " printed")
        need(F(n, d) == expected, name + " exact agreement")

    need(sha256(args.report.read_bytes()).hexdigest() == REPORT_SHA,
         "frozen primary canonical report pin")
    need(args.report.read_bytes() == args.optimized.read_bytes(),
         "normal and optimized canonical reports byte identical")
    need(sha256(args.proof.read_bytes()).hexdigest() == PROOF_SHA,
         "frozen mathematical proof pin")
    need(report["status"] == "PASS_ACTUAL_REGISTERED_PAST_MASTER_PRECONTACT",
         "primary actual-packet report passes")
    ell, Z, mu = F(1, 10**6), F(8, 10**6), F(103, 10**9)
    need(Z == val("static_collar_length") and 131*Z < 9,
         "wide collar length and proper endpoint margin")
    outer = {"rho": (F(129), F(131)),
             "eta": (-F(9, 10**5), -F(8, 10**5)),
             "p": (-F(7, 50000), -F(3, 25000)),
             "a": (F(111, 1000), F(14, 125))}
    for key, (lo, hi) in outer.items():
        row = report["first_exit"][key]
        initial = tuple(map(F, row["shell"]))
        rhs = tuple(map(F, row["independent_monomial_rhs"]))
        flow = (initial[0]+min(F(0), Z*rhs[0]),
                initial[1]+max(F(0), Z*rhs[1]))
        need(flow == tuple(map(F, row["flow_tube"])),
             key + " exact first-exit interval image")
        need(lo < flow[0] <= flow[1] < hi and
             flow[0]-lo == F(row["strict_lower_margin"]) and
             hi-flow[1] == F(row["strict_upper_margin"]),
             key + " two strict printed-box margins")
    need(interval("derived_h_interval") == (F(7, 4), F(7, 3)),
         "polynomial quotient guards")
    need(interval("derived_U2_interval") == (F(31, 10), F(25, 8)),
         "second potential derivative guards")
    need(interval("derived_Q_interval") == (F(1, 3), F(6)),
         "static potential coefficient guards")
    printed_fraction(1933333317, 5000000000,
                     interval("independent_Q_majorants")[0], "lower static potential majorant")
    printed_fraction(2680267, 450000, interval("independent_Q_majorants")[1],
                     "upper static potential majorant")
    need(F(129)**2/3-F(35, 4)*ell**2 > 5500 and
         F(131)**2*6+F(35, 4)*ell**2 < 103000,
         "fresh continuous negative-time potential outward bounds")
    need(interval("fresh_negative_time_W_interval") == (F(5500), F(103000)) and
         "5500<\\mathcal W<103000" in source,
         "negative-time potential interval printed")
    need(interval("negative_time_interval") == (-ell, F(0)),
         "new negative-time coefficient rectangle")
    printed_fraction(103, 1000000000, val("fresh_scaled_potential_upper"),
                     "scaled potential guard")
    need(ell**2*103000 == mu, "scaling of actual master potential")
    beta, B = mu/2, (1+ell)/(1-mu/2)
    printed_fraction(103, 2000000000,
                     val("uniform_estimate.volterra_contraction_upper"),
                     "uniform Cauchy Volterra contraction")
    printed_fraction(2000002000, 1999999897,
                     val("uniform_estimate.C0_and_relative_L2_upper"),
                     "Cauchy C0 and relative L2 factor")
    need(beta < 1 and B < F("1.000002"), "printed uniform Cauchy norm guard")
    delta = ell/2+mu*B/2
    printed_fraction(2206000103, 3999999794000000,
                     val("uniform_estimate.kappa"), "principal relative norm error")
    need(delta == val("uniform_estimate.kappa") < F(1, 10**6),
         "uniform relative error is below one millionth")
    need(1-delta == val("uniform_estimate.branch_norm_lower_factor") and
         1+delta == val("uniform_estimate.branch_norm_upper_factor"),
         "two principal norm factors reconstructed")
    lower = (1-delta)**2/((1-delta)**2+(1+delta)**2)
    upper = (1+delta)**2/((1-delta)**2+(1+delta)**2)
    need(lower == val("uniform_estimate.principal_share_lower") and
         upper == val("uniform_estimate.principal_share_upper"),
         "normalized principal quadratic shares reconstructed")
    need(F("0.499999") < lower <= upper < F("0.500001") and
         interval("uniform_estimate.safe_strict_share_interval") ==
         (F("0.499999"), F("0.500001")),
         "printed strict near-half interval is outward")
    need("0.499999<\\frac{E_\\pm(x)}{E_+(x)+E_-(x)}<0.500001" in source,
         "both-share quantitative display present")
    need(report["uniform_past_scaled_time"]["interval"] == "-1<x<=0" and
         interval("uniform_past_scaled_time.quarter_support") ==
         (F(3, 4), F(9, 4)), "open precontact and quarter support binding")

    beta_half = mu/8
    B_half = (1+ell/2)/(1-beta_half)
    epsilon = beta_half*B_half
    need(B_half == val("constructive_past_diamond.B_half"),
         "half-time Cauchy factor")
    printed_fraction(206000103, 15999999794000000,
                     val("constructive_past_diamond.C0_free_error_upper"),
                     "two-face C0 free-wave error")
    need(epsilon == val("constructive_past_diamond.C0_free_error_upper"),
         "half-time free-wave error reconstructed")
    mid = F(1, 2)-ell/4-epsilon
    need(mid == val("constructive_past_diamond.two_face_midpoint_lower") >
         F("0.499999"), "actual two-face midpoint lower display")
    exponent = F(4, 3)
    expupper = 1+exponent+exponent**2/2+(exponent**3/6)/(1-exponent/4)
    printed_fraction(103, 27, val("constructive_past_diamond.exp_4_over_3_upper"),
                     "positive Taylor-tail exponential guard")
    need(expupper == F(103, 27) < 4, "analytic bump lower-bound Taylor arithmetic")
    corner = -ell*(1-2*ell)/16+epsilon
    need(corner == val("constructive_past_diamond.common_apex_upper") <
         -F(1, 25000000) == -F("4e-8"), "negative common apex upper display")
    printed_fraction(103, 16000000000,
                     val("constructive_past_diamond.goursat_contraction_upper"),
                     "local two-face Goursat contraction")
    need(mu/16 == val("constructive_past_diamond.goursat_contraction_upper") < 1,
         "null integration rectangle contraction")
    need("2f_v=J_-" in source and "2g_u=J_+" in source,
         "exact face tangent and principal factors")

    polygon = re.search(r"\\fill\[blue!10\]\s*(.*?)--cycle;", source)
    need(polygon is not None, "half-diamond figure polygon present")
    points = [tuple(map(F, xy)) for xy in
              re.findall(r"\(([0-9.-]+),([0-9.-]+)\)", polygon.group(1))]
    need(points == [(F(1), F(0)), (F(2), F(0)), (F(3, 2), -F(1, 2))],
         "actual plotted base and interior apex")
    for r, x in points:
        u, v = x+r, x-r
        need(u >= 1 and v >= -2 and u+v <= 0 and r > 0,
             "plotted vertex is inside actual contact-free null domain " + str((r, x)))
    area = abs(sum((a[0]*b[1]-a[1]*b[0] for a, b in
                    zip(points, points[1:]+points[:1])), F(0)))/2
    need(area == F(1, 4), "physical coordinate half-diamond area")
    for r, x, u, v in [(F(5, 4), -F(1, 4), F(1), -F(3, 2)),
                        (F(7, 4), -F(1, 4), F(3, 2), -F(2))]:
        need(x+r == u and x-r == v, "face midpoint and null labels " + str((r, x)))
    need(report["activity"]["exact_background_ODE_image_recomputations"] == 4 and
         all(report["activity"][key] == 0 for key in
             ["background_oracle_calls", "bump_queries", "legacy_checker_invocations",
              "numerical_PDE_solves", "project_imports", "source_oracle_calls"]),
         "reported exact ODE arithmetic and zero scientific oracle/PDE activity")
    need(report["selected_local_inputs_authenticated"] == 5,
         "five selected local input copies")
    result = {
        "status": "PASS_DISPLAYED_PAST_CHARACTERISTIC_PRECURSOR_BOUNDS",
        "method": "stdlib Fraction; printed constant, report, first-exit image and diagram audit; no producer imports",
        "finite_checks": len(checks), "checks": checks,
        "audit_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "inputs": {name: {"filename": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                   for name, path in [("source", tex), ("primary_report", args.report),
                                      ("optimized_report", args.optimized), ("primary_proof", args.proof)]},
        "scope": "Finite displayed bounds only; analytic branch, reconstruction and smooth Volterra arguments are reviewed separately."
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "finite_checks": len(checks)}, sort_keys=True))


if __name__ == "__main__":
    main()
