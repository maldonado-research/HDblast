#!/usr/bin/env python3
"""Audit manuscript displays against frozen source/theory reports, with Fraction."""
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, isqrt
from pathlib import Path
import argparse
import json
import re

HERE = Path(__file__).resolve().parent
STEM = "INCOMING_PREPARATION_LEDGER"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", required=True, type=Path)
    ap.add_argument("--theory", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    tex = HERE / (STEM + ".tex")
    text = tex.read_text()
    src = json.loads(args.sources.read_text())
    theory = json.loads(args.theory.read_text())
    checks = []

    def need(condition, name):
        if not condition:
            raise ValueError(name)
        checks.append(name)

    def value(record, path):
        for part in path.split("."):
            record = record[part]
        return F(record)

    def printed_fraction(n, d, expected, name):
        pattern = r"\\frac\{\s*" + str(n) + r"\s*\}\{\s*" + str(d) + r"\s*\}"
        need(re.search(pattern, text) is not None or str(n) + "/" + str(d) in text,
             name + " printed")
        need(F(n, d) == expected, name + " exact agreement")

    def area(poly):
        return abs(sum((a[0]*b[1]-a[1]*b[0] for a, b in
                        zip(poly, poly[1:]+poly[:1])), F(0)))/2

    def clip(poly, threshold, keep_above):
        out = []
        for a, b in zip(poly, poly[1:]+poly[:1]):
            va, vb = sum(a)-threshold, sum(b)-threshold
            inside_a = va >= 0 if keep_above else va <= 0
            inside_b = vb >= 0 if keep_above else vb <= 0
            if inside_a:
                out.append(a)
            if inside_a != inside_b:
                weight = -va/(vb-va)
                out.append(tuple(a[j]+weight*(b[j]-a[j]) for j in range(2)))
        return out

    need(src["status"] == "PASS_EXACT_ZERO_ORACLE_SOURCE_LEDGER",
         "canonical source postprocessing passes")
    need(theory["status"] == "PASS_CONDITIONAL_INCOMING_PREPARATION_IDENTITIES",
         "canonical conditional theory passes")
    need(theory["operator_caps_remain_conditional_future_assertions"] is True,
         "continuous caps retain conditional future domain")
    need(src["mean_gram"] == [1, 2, 1, 3, 1], "actual five-channel Gram")
    ell, T, delta, k, alpha = F(1, 10**6), F(5, 4), F(1, 4), F(1, 2), F(1, 2)
    need(ell**2/4 == F(1, 4000000000000) and
         ell**2/8 == F(1, 8000000000000), "principal imbalance and half-weight guards")
    need("4000000000000" in text and "8000000000000" in text,
         "principal rational guards printed")
    need(T == value(src, "exact_geometry.T") and T*T/2 == F(25, 32),
         "full future cone area")
    need((T*T-1)/2 == value(src, "exact_geometry.active_area") == F(9, 32),
         "exact zero wedge leaves area nine thirty-seconds")
    need(T-1 == value(src, "exact_geometry.active_boundary_length") == delta,
         "active shell segment length")
    fills = re.findall(r"\\fill\[(?!gray)[^\]]+\]\s*(.*?)--cycle;", text)
    need(len(fills) == 5, "five active polygons embedded in figure")
    rows = {tuple(row["old_cell_index"]): row
            for row in src["retained_bulk_cells_intersecting_cone"]}
    cell_order = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 0)]
    figure_areas = []
    for vertices, index in zip(fills, cell_order):
        poly = [tuple(map(F, xy)) for xy in
                re.findall(r"\(([0-9.]+),([0-9.]+)\)", vertices)]
        i, j = index
        square = [(F(i, 2), F(j, 2)), (F(i+1, 2), F(j, 2)),
                  (F(i+1, 2), F(j+1, 2)), (F(i, 2), F(j+1, 2))]
        exact = clip(clip(square, T, False), F(1), True)
        figure_areas.append(area(poly))
        need(area(poly) == area(exact) == F(rows[index]["active_area"]),
             "figure and independently clipped cell " + str(index))
    need(sum(figure_areas, F(0)) == F(9, 32), "five displayed polygon areas sum")
    need(F(rows[(0, 0)]["active_area"]) == 0, "lower left old cell inactive")

    a = b = F(1, 10000)
    d, mu = F(1, 1000000), F(1, 2000)
    cv = a+(b+k*a)/2+d/2+F(1, 2)
    cp = k*b+(b+k*a)/2+k*d/2
    cu = d/2+k*d/2+F(1, 2)
    bulk = max(2*cv/(1-k), 2*cp/(1-k), 2*cu)
    gamma = alpha+(1-k)*mu
    boundary = 2/alpha*(gamma**2/k+(1-k)*mu**2+alpha*mu)
    need(bulk == value(theory, "forced_energy.homogeneous_bulk_rate") ==
         F("2.000702"), "homogeneous bulk rate reconstructed")
    need(boundary == value(theory, "forced_energy.homogeneous_boundary_rate") ==
         F("2.003001"), "nonsymmetric Robin rate reconstructed")
    need("2.000702" in text and "2.003001" in text, "homogeneous rates printed")
    C, CG = F(301, 100), 1+(F(1, 2)+mu)**2
    need(max(bulk+1, boundary+1) < C, "forced growth rate guard")
    printed_fraction(301, 100, value(theory, "forced_energy.forced_growth_guard"),
                     "forced growth")
    printed_fraction(5002001, 4000000,
                     value(theory, "forced_energy.boundary_square_factor"),
                     "actual Robin source square factor")
    need(CG == value(src, "energy_coefficients_conditional_on_accompanying_forced_proof.Cg"),
         "source and theory boundary energy factors agree")
    need(F(1, 2) == value(theory, "forced_energy.bulk_square_factor") ==
         value(src, "energy_coefficients_conditional_on_accompanying_forced_proof.Cf"),
         "sharp bulk source square factor")
    argument = C*T
    n = int(theory["forced_energy"]["Taylor_order"])
    expupper = sum((argument**j/F(factorial(j)) for j in range(n+1)), F(0))
    expupper += argument**(n+1)/F(factorial(n+1))/(1-argument/(n+2))
    need(argument == F(301, 80) == value(theory, "forced_energy.exponential_argument"),
         "growth exponential argument")
    need(expupper == value(theory, "forced_energy.exponential_rational_upper") < 44,
         "rational Taylor-tail exponential bound")
    need(2/alpha*44 == 176 and "|w(x,0)|^2\\le176" in text,
         "printed shell trace square factor")

    FH, GH = (value(src, "registered_caps." + key) for key in ["F_H", "G_H"])
    need(str(FH.numerator) in text and str(GH.numerator) in text,
         "retained global actual-source caps printed")
    exact_bulk = sum((F(row["active_area"])*F(row["chosen_square_cap"])
                      for row in rows.values()), F(0))
    exact_bd = GH**2*delta
    budget = src["budgets"]["wedge_zero_retained_cell_caps"]
    need(exact_bulk == F(budget["bulk_L2_squared_upper"]),
         "all five retained cell source squares integrated exactly")
    need(exact_bd == F(budget["boundary_L2_squared_upper"]) == F(3708320490000),
         "active actual Robin source integral")
    exact_B = exact_bulk/2+CG*exact_bd
    need(exact_B == F(budget["energy_source_ledger_upper"]),
         "continuous energy source ledger recomputed")
    labels = [("bulk L2 budget", 8526418129275428520152526, exact_bulk),
              ("boundary L2 budget", 3708320490000, exact_bd),
              ("energy source ledger", 4263209064642351515776089, exact_B)]
    for label, integer, exact in labels:
        need(str(integer) in text, label + " outward integer printed")
        need(integer-1 < exact <= integer, label + " safe integer ceiling")
    need(F(9, 32)*FH**2 < F("8.50291e25") and
         exact_bulk < F("8.52642e24"), "printed scientific source comparisons")
    trace_squared = 176*exact_B
    need(trace_squared == F(budget["conditional_zero_initial_mean_trace_squared_upper"]),
         "auxiliary zero-state trace squared bound")
    root = isqrt(trace_squared.numerator//trace_squared.denominator)
    if root*root < trace_squared:
        root += 1
    need(root == 27392057158547 and str(root) in text,
         "printed minimal outward integer auxiliary trace")
    need(root**2 >= 176*F(4263209064642351515776089),
         "printed trace safely follows printed energy ledger ceiling")

    bM, cM = F(917, 62500000), F(131, 10500000000)
    Hw, Hrw, Vw = cM/3, (3*ell**2+bM)*(cM/3)+cM, 3*ell**2*(cM/3)+cM
    for key, v in [("cM", cM), ("Hw", Hw), ("Hrw", Hrw), ("Vw", Vw)]:
        need(v == value(theory, "inverse_history." + key), "inverse initial map " + key)
    coefficient = max((1+2*Hw+2*T*Vw)/3, 1+2*Hrw+2*Vw)
    Ng = F(20000001, 20000000)
    need(coefficient == value(theory, "inverse_history.N_exact_majorant") < Ng,
         "complete initial norm inverse guard")
    printed_fraction(20000001, 20000000, value(theory, "inverse_history.N_guard"),
                     "inverse initial norm guard")
    q = F(1201, 3200000)
    printed_fraction(1201, 3200000, value(theory, "inverse_history.q_retained"),
                     "retained coupled resolvent")
    kap = q*Ng/(1-q)
    printed_fraction(24020001201, 63975980000000,
                     value(theory, "inverse_history.kappa_exact"), "history inverse contraction")
    printed_fraction(47, 125000, value(theory, "inverse_history.kappa_guard"),
                     "history inverse decimal guard")
    need(kap == value(theory, "inverse_history.kappa_exact") < F(47, 125000) < 1,
         "restricted history inverse is a strict contraction")
    need(all(n == 0 for n in src["scientific_activity"].values()) and
         all(theory["verification"][key] == 0 for key in
             ["new_PDE_solves", "new_background_or_evolved_field_queries", "new_bump_queries"]),
         "reports record zero new scientific oracles and PDE solves")
    result = {
        "status": "PASS_DISPLAYED_INCOMING_PREPARATION_BOUNDS",
        "method": "stdlib Fraction; independent clipping and display arithmetic; no producer imports",
        "audit_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "finite_checks": len(checks), "checks": checks,
        "inputs": {name: {"filename": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                   for name, path in [("source", tex), ("sources_report", args.sources),
                                      ("theory_report", args.theory)]},
        "scope": "Printed finite constants only; retained analytic and physical premises are not reproved."
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "finite_checks": len(checks)}, sort_keys=True))


if __name__ == "__main__":
    main()
