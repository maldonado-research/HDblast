#!/usr/bin/env python3
"""Exact rational spot checks supplementary to the displayed all-point algebra."""
import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)


def run():
    rows = []
    for n in (1, 2, 4):
        rho = Q(17+n, 3)
        ir = 1/rho
        a = Q(2+n, 19)
        p = -Q(3+n, 23)
        h = Q(5+n, 11)
        u = Q(7+n, 13)
        v = -Q(4+n, 17)
        vy = Q(9+n, 29)
        E = Q(13+n, 31)
        U2 = Q(8+n, 37)
        d = 3*ir/p
        py = p*(h-4*a)
        ay = -ir**2-p**2/3
        dy = d*(3*a-h)
        uy = -p*v/3-2*a*u
        uyy = -py*v/3-p*vy/3-2*ay*u-2*a*uy
        ut = E*u
        uyt = E*(uy+ir*u)
        vt = E*(v-d*u)
        vyt = E*(vy-d*uy-dy*u+ir*(v-d*u))
        Fmom = ir*(18*u*ut-v*vt)
        Gmom = ir*(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)
        Fmom_y = ir*(18*uy*ut+18*u*uyt-vy*vt-v*vyt)-a*Fmom
        Khat = E*ir/3*(18*u*u-v*v+d*u*v)
        Pihat = E*ir/p*(d*u*vy+(d*(h-ir)+10*p)*u*v
                           +(30*a-6*ir)*u*u+(2*ir-a)*v*v)
        Aoriginal = (-24*uy*uy-288*a*u*uy-24*u*uyy
                     -(192*ay+384*a*a+32*p*p)*u*u
                     -ir*ir*(336*u*u+12*ut*ut+2*vt*vt)
                     -2*vy*vy-16*p*u*vy-2*U2*v*v)
        Areduced = (-2*vy*vy-8*p*u*vy+8*p*(h+2*a)*u*v
                    +(16*p*p-192*ir*ir)*u*u
                    -(2*U2+Q(8,3)*p*p)*v*v
                    -ir*ir*(12*ut*ut+2*vt*vt))
        checks = {"K_equals_Fmom_over_3": Khat == Fmom/3,
                  "Pi_equals_Gmom_minus_Fmom_y_over_p": Pihat == (Gmom-Fmom_y)/p,
                  "Hamiltonian_A_reduction": Aoriginal == Areduced}
        for key, value in checks.items():
            need(value, key)
        rows.append({"rational_jet_label": n, "checks": checks})
    return {"status": "PASS_EXACT_RATIONAL_FORMULA_SPOT_CHECKS",
            "scope": "Supplementary generic-jet checks; all-point proof is algebra in INITIAL_FORMULA_REVIEW.md.",
            "numeric_capsule_imports": False, "rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output
    need(not output.exists(), "refuse overwrite")
    result = run()
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])
