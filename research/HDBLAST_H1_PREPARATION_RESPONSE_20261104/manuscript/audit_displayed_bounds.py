#!/usr/bin/env python3
"""Fraction audit of the manuscript's exact constants and outward displays."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
from math import factorial
import argparse
import json
import re

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--design", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    tex = HERE / "H1_PREPARATION_RESPONSE.tex"
    source = tex.read_text()
    report, base, design = (json.loads(p.read_text()) for p in
                            [args.report, args.base, args.design])
    checks = []

    def need(condition, name):
        if not condition:
            raise ValueError(name)
        checks.append(name)

    def scalar(path):
        out = report
        for part in path.split("."):
            out = out[part]
        return F(out)

    def interval(path):
        out = report
        for part in path.split("."):
            out = out[part]
        return tuple(map(F, out))

    def printed_fraction(numerator, denominator, expected, name):
        pattern = re.compile(r"\\frac\{\s*" + str(numerator) +
                             r"\s*\}\s*\{\s*" + str(denominator) + r"\s*\}")
        need(pattern.search(source) is not None or
             str(numerator) + "/" + str(denominator) in source,
             name + " present in source")
        need(F(numerator, denominator) == expected, name + " exact agreement")

    def sci(text):
        coefficient, power = text.replace(" ", "").split(r"\times10^")
        return F(coefficient) * F(10) ** int(power.strip("{}"))

    need(report["status"] == "PASS_CONDITIONAL_H1_PREPARATION_RESPONSE_BOUNDS",
         "conditional report passes")
    need(report["continuous_operator_caps_are_retained_conditional_assertions"] is True,
         "continuous coefficient caps explicitly remain premises")
    need(report["mean_Gram"] == [1, 2, 1, 3, 1], "five-channel Gram")
    for quantity, key in [("v_2", "v2"), (r"\Theta_2", "Theta2")]:
        match = re.search(re.escape(quantity) + r"\\in\[([^\]]+)\]", source)
        need(match is not None, key + " printed interval found")
        printed = tuple(map(sci, match.group(1).split(",")))
        exact = interval("whole_class_exact_intervals." + key)
        need(printed == interval("whole_class_outward_intervals." + key),
             key + " printed interval equals certified outward display")
        need(printed[0] <= exact[0] <= exact[1] <= printed[1] < 0,
             key + " exact interval contained and strictly negative")

    printed_fraction(17161, 19200000000000000000,
                     scalar("observable_difference_caps_per_R.v2"), "C_v")
    printed_fraction(131, 2400000000000,
                     scalar("observable_difference_caps_per_R.Theta2"), "C_Theta")
    printed_fraction(2003001, 1000000,
                     scalar("energy.boundary_growth_exact_upper"), "Robin growth")
    printed_fraction(1000351, 500000,
                     scalar("energy.bulk_growth_exact_upper"), "bulk growth")
    printed_fraction(13, 1250000, scalar("energy.boundary_value_square_guard"),
                     "boundary trace square")
    k = alpha = F(1, 2)
    a = b = F(1, 10000)
    d = F(1, 1000000)
    mu = F(1, 2000)
    gamma = alpha + (1-k)*mu
    cbdry = 2/alpha*(gamma**2/k+(1-k)*mu**2+alpha*mu)
    cv = a+(b+k*a)/2+d/2+F(1, 2)
    cp = k*b+(b+k*a)/2+k*d/2
    cu = d/2+k*d/2+F(1, 2)
    cbulk = max(2*cv/(1-k), 2*cp/(1-k), 2*cu)
    C, T, delta, ell = F(201, 100), F(5, 4), F(1, 4), F(1, 10**6)
    need(cbdry == scalar("energy.boundary_growth_exact_upper") < C,
         "nonsymmetric Robin Young rate reconstructed")
    need(cbulk == scalar("energy.bulk_growth_exact_upper") < C,
         "bulk rate reconstructed")
    need(scalar("energy.growth_guard") == C and
         scalar("energy.tilt_k") == k and scalar("energy.boundary_alpha") == alpha,
         "printed energy constants")
    x = C*T
    expupper = sum((x**j/F(factorial(j)) for j in range(25)), F(0))
    expupper += (x**25/F(factorial(25)))/(1-x/26)
    need(expupper == scalar("energy.exponential_rational_upper") < 13,
         "printed positive Taylor exponential certificate")
    need(scalar("energy.exponential_Taylor_order") == 24,
         "printed Taylor truncation order")
    H = F(131, 73426500000000)
    bmax, gmax = F(917, 62500000), F(131, 16317000000000)
    Hr = 3*ell**2*H+(bmax+gmax)*H+F(2, 3)*gmax
    V = Hr+bmax*H
    h2, hr2, vel2 = F(21, 20)*H, F(21, 20)*Hr, F(21, 20)*V
    beta = F(1, 1554)
    betar = ell*131*(F(1, 4000)/(2*F(111, 1000))+
                      F(1, 7000)*F(1, 10000)/(2*F(111, 1000)**2))
    w2 = beta*(h2+F(2, 3))
    wr2 = betar*(h2+F(2, 3))+beta*(hr2+F(2, 3))
    caps = {"h": h2, "h_r": hr2, "V_L": vel2, "w": w2, "w_r": wr2,
            "Gram_U0": 2*h2+w2, "Gram_U0_r": 2*hr2+wr2, "Gram_U0_x": 2*vel2}
    for name, value in caps.items():
        need(value == scalar("initial_L2_caps_per_R." + name),
             "initial L2 map " + name + " reconstructed")
    e0 = ((2*h2+w2)**2+(2*hr2+wr2)**2+(2*vel2)**2)/2
    need(e0 == scalar("energy.initial_energy_per_R_squared") < F(1, 5000000),
         "complete initial energy gate")
    need(scalar("initial_L2_caps_per_R.initial_tilt_cross_term") == 0 and
         scalar("initial_L2_caps_per_R.initial_boundary_value") == 0,
         "initial tilt and boundary terms exactly zero")
    need(F(17, 4) < F(21, 10)**2, "printed Sobolev rational guard")
    need(2*13/(alpha*5000000) == F(13, 1250000) < F(1, 300**2),
         "printed scalar trace factor")
    Cv = (131*ell)**2*F(1, 1000)*delta**2/(4*300)
    Ct = 131*ell*F(1, 1000)*delta/(2*300)
    need(Cv == scalar("observable_difference_caps_per_R.v2") and
         Ct == scalar("observable_difference_caps_per_R.Theta2"),
         "proper-time transfer factors and kernel integrals")
    radius = F(2*10**31)
    need(radius == scalar("class.radius_closed") ==
         F(design["parameters"]["closed_H1_radius"]), "predeclared closed class radius")
    for key, cap, printed in [("v2", Cv*radius, F("1.78761e16")),
                              ("Theta2", Ct*radius, F("1.09167e21"))]:
        need(cap == scalar("whole_class_absolute_difference_caps." + key) < printed,
             "printed absolute " + key + " increment cap")
        bkey = "proper_v2" if key == "v2" else "proper_Theta2"
        lo, hi = map(F, base[bkey]["exact"])
        need((lo-cap, hi+cap) == interval("whole_class_exact_intervals." + key),
             key + " complete base-plus-difference interval")
    for bkey, lo, hi in [("proper_v2", "-5.4520e16", "-2.6521e16"),
                         ("proper_Theta2", "-9.0208e21", "-4.4458e21")]:
        blo, bhi = map(F, base[bkey]["exact"])
        need(F(lo) <= blo <= bhi <= F(hi) < 0, bkey + " printed base enclosure")
    need(report["verification"]["finite_checks"] == 24 and
         "The 24 finite arithmetic gates" in source, "printed energy gate count")
    need(all(report["verification"][name] == 0 for name in
             ["new_PDE_solves", "new_background_or_evolved_field_queries", "new_bump_queries"]),
         "zero new science oracle/PDE/bump calls")
    result = {
        "status": "PASS_DISPLAYED_H1_RESPONSE_BOUNDS",
        "method": "stdlib Fraction; no producer imports or new science oracle calls",
        "inputs": {name: {"filename": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}
                   for name, path in [("source", tex), ("report", args.report),
                                      ("base", args.base), ("design", args.design)]},
        "finite_checks": len(checks), "checks": checks
    }
    args.output.write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "finite_checks": len(checks)}, sort_keys=True))


if __name__ == "__main__":
    main()
