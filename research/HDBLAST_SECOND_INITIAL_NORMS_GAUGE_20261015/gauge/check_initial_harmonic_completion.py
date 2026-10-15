#!/usr/bin/env python3
"""Exact point-jet checks of the explicit INITIAL relative-harmonic completion.

Independent implementation: form all five-dimensional Christoffel symbols
from metric coefficient jets, contract their difference with the perturbed
inverse metric, and compare both Taylor orders. No field evolution is tested.
"""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path

ZERO = (Q(0), Q(0), Q(0))


def jet(*x):
    return tuple(Q(v) for v in x) + (Q(0),) * (3 - len(x))


def add(*xs):
    return tuple(sum(x[k] for x in xs) for k in range(3))


def scale(c, x):
    return tuple(c * v for v in x)


def mul(x, y):
    return tuple(sum(x[j] * y[k-j] for j in range(k+1)) for k in range(3))


def matmul(a, b):
    return [[add(*(mul(a[i][k], b[k][j]) for k in range(5)))
             for j in range(5)] for i in range(5)]


def inverse(g):
    # Coefficient recursion for a matrix series with diagonal constant term.
    inv = [[[Q(0) for _ in range(5)] for _ in range(5)] for _ in range(3)]
    for i in range(5):
        inv[0][i][i] = 1 / g[i][i][0]
    for n in (1, 2):
        for i in range(5):
            for j in range(5):
                inv[n][i][j] = -inv[0][i][i] * sum(
                    g[i][k][r] * inv[n-r][k][j]
                    for r in range(1, n+1) for k in range(5))
    return [[tuple(inv[n][i][j] for n in range(3))
             for j in range(5)] for i in range(5)]


def christoffel(ginv, dg):
    return [[[scale(Q(1, 2), add(*(
        mul(ginv[a][d], add(dg[b][d][c], dg[c][b][d], scale(-1, dg[d][b][c])))
        for d in range(5)))) for c in range(5)] for b in range(5)] for a in range(5)]


MUTATIONS = (
    "omit_lapse_quadratic", "halve_radial_U_Uy", "omit_radial_fy",
    "halve_angular_quadratic", "omit_mixed_metric_cross",
    "omit_second_metric_lapse_value",
)


def residual(values, mutation=None, longitudinal=False):
    r, a, U, Uy, T, f, fy, q, DU, Df = values
    r2 = r*r
    b1y = -6*r2*(Uy+4*a*U)
    b2y = r2*(-74*U*Uy-Q(3, 2)*fy-194*a*U*U-3*a*f)
    b2ang = [22*U*du-Q(1, 2)*df for du, df in zip(DU, Df)]
    ndot2 = 4*r2*q-24*r*U*T
    if longitudinal:
        b1y = Q(0)
    if mutation == "omit_lapse_quadratic":
        ndot2 = 4*r2*q
    elif mutation == "halve_radial_U_Uy":
        b2y += 37*r2*U*Uy
    elif mutation == "omit_radial_fy":
        b2y += Q(3, 2)*r2*fy
    elif mutation == "halve_angular_quadratic":
        b2ang = [11*U*du-Q(1, 2)*df for du, df in zip(DU, Df)]

    g = [[ZERO for _ in range(5)] for _ in range(5)]
    dg = [[[ZERO for _ in range(5)] for _ in range(5)] for _ in range(5)]
    g[0][0] = jet(-r2, -2*r2*U, -r2*U*U)
    g[1][1] = jet(1, -4*U)
    for j in range(2, 5):
        g[j][j] = jet(r2, 2*r2*U, r2*f/2)

    # Spatial derivatives of the actual initial lapse/spatial metric.
    dg[1][0][0] = jet(-2*a*r2, -2*r2*(2*a*U+Uy), -2*r2*(a*U*U+U*Uy))
    dg[1][1][1] = jet(0, -4*Uy)
    for j in range(2, 5):
        dg[1][j][j] = jet(2*a*r2, r2*(4*a*U+2*Uy), r2*(a*f+fy/2))
    for k, (du, df) in enumerate(zip(DU, Df), 2):
        dg[k][0][0] = jet(0, -2*r2*du, -2*r2*U*du)
        dg[k][1][1] = jet(0, -4*du)
        for j in range(2, 5):
            dg[k][j][j] = jet(0, 2*r2*du, r2*df/2)

    # Physical K=+1/2 L_n gamma: gamma_t=2 N K at zero initial shift.
    dg[0][1][1] = jet(0, -4*T, r*q-4*U*T)
    for j in range(2, 5):
        dg[0][j][j] = jet(0, 2*r2*T, r2*(r*q+2*U*T))
    # g_tt,t=-2 N N_t and g_tI,t=gamma_IJ beta_t^J.
    dg[0][0][0] = jet(0, -2*r2*T, -r*ndot2-2*r2*U*T)
    cross = -4*U*b1y if mutation != "omit_mixed_metric_cross" else 0
    dg[0][0][1] = dg[0][1][0] = jet(0, b1y, b2y/2+cross)
    for j, b2 in enumerate(b2ang, 2):
        dg[0][0][j] = dg[0][j][0] = jet(0, 0, r2*b2/2)
    if mutation == "omit_second_metric_lapse_value":
        g[0][0] = jet(*g[0][0][:2])
        for d in range(1, 5):
            dg[d][0][0] = jet(*dg[d][0][0][:2])

    ginv = inverse(g)
    ident = matmul(g, ginv)
    for i in range(5):
        for j in range(5):
            if ident[i][j] != jet(int(i == j)):
                raise RuntimeError("inverse series check failed")
    Gamma = christoffel(ginv, dg)
    # The reference connection is fixed, not contracted with g0 in C.
    return [add(*(mul(ginv[b][c], add(Gamma[A][b][c],
                    jet(-Gamma[A][b][c][0])))
                  for b in range(5) for c in range(5))) for A in range(5)]


def cases():
    for j in (1, 2, 3):
        yield (Q(129+j), Q(2+j, 7), Q(3+j, 11), Q(-2-j, 13),
               Q(5+j, 17), Q(-7-j, 19), Q(9+j, 23), Q(4+j, 29),
               tuple(Q(j+2*k+1, 31+k) for k in range(3)),
               tuple(Q(-j-k-2, 37+k) for k in range(3)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", required=True)
    ap.add_argument("--run-controls", action="store_true")
    ap.add_argument("--mutation", choices=MUTATIONS)
    args = ap.parse_args()
    witnesses = []
    for idx, data in enumerate(cases(), 1):
        C = residual(data, args.mutation)
        if any(x != ZERO for x in C):
            raise SystemExit("FAIL: nonzero initial relative-harmonic residual under "
                             + str(args.mutation) + " at point-jet " + str(idx))
        Clong = residual(data, longitudinal=True)
        r, a, U, Uy = data[:4]
        target = -6*(Uy+4*a*U)
        if Clong[1][1] != target or any(Clong[i][1] for i in (0, 2, 3, 4)):
            raise SystemExit("FAIL: longitudinal first gauge defect")
        witnesses.append({"point_jet": idx, "C0_C1_C2_all_zero": True,
                          "longitudinal_C1_y": str(target)})
    controls = []
    if args.run_controls:
        for name in MUTATIONS:
            bad = [i for i, data in enumerate(cases(), 1)
                   if any(x != ZERO for x in residual(data, name))]
            if not bad:
                raise SystemExit("FAIL: wrong-formula control survived: " + name)
            controls.append({"mutation": name, "nonzero_residual_point_jets": bad})
    report = {
        "scope": "Initial gauge algebra only; no gauge wave or boundary evolution solved",
        "factorial_convention": "series epsilon^2 coefficients equal one half second variation",
        "independent_method": "direct five-dimensional inverse metric and Christoffel contraction",
        "positive_checks": witnesses, "wrong_formula_controls": controls,
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(args.output).write_text(json.dumps(report, indent=2) + "\n")
    print("PASS: three exact point-jet checks; " + str(len(controls)) + " controls")


if __name__ == "__main__":
    main()
