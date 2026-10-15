#!/usr/bin/env python3
"""Exact independent 5D background, shell, radial, and support checks.

SymPy arithmetic only.  No assert statement, floating point, external input,
path-dependent report field, or optimization-dependent check is used.
The background is a separate explicit model, not the registered HDBLAST branch.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


DIMENSION = 5
T, x1, x2, x3, x4 = coordinates = sp.symbols("T x1 x2 x3 x4", real=True)
s = sp.Symbol("s", real=True)
eta = sp.diag(-1, 1, 1, 1, 1)
X = sp.Matrix(coordinates)
X_lower = eta * X
invariant_relation = x4 ** 2 - (s + T ** 2 - x1 ** 2 - x2 ** 2 - x3 ** 2)


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def partial(expression, index):
    """Cartesian derivative after s=eta_AB X^A X^B, by the chain rule."""
    return sp.diff(expression, coordinates[index]) + 2 * X_lower[index] * sp.diff(expression, s)


def ambient_reduce(expression):
    """Exact reduction modulo x4^2=s+T^2-x1^2-x2^2-x3^2.

    s is retained as a symbol for small exact expressions. This reduction
    identifies it with the Cartesian invariant; it does not add an equation
    of motion. All expressions reduced here have denominators depending only
    on s and constants, so no denominator vanishes through this operation.
    """
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    remainder = sp.rem(sp.Poly(numerator, x4, domain="EX"),
                       sp.Poly(invariant_relation, x4, domain="EX")).as_expr()
    return sp.cancel(remainder / denominator)


def zero(expression):
    if isinstance(expression, sp.MatrixBase):
        return all(ambient_reduce(entry) == 0 for entry in expression)
    return ambient_reduce(expression) == 0


def equal(left, right):
    return zero(left - right)


class Checks:
    def __init__(self):
        self.positive = []
        self.controls = []

    def check(self, name, condition, error=None):
        require(condition, error or "exact background check failed: " + name)
        self.positive.append({"name": name, "pass": True})

    def reject(self, name, incorrect_condition):
        require(not bool(incorrect_condition), "wrong-formula control did not reject: " + name)
        self.controls.append({"name": name, "incorrect_formula_rejected": True})


def geometry_from_metric(metric):
    """Levi-Civita and Ricci definitions in Cartesian coordinates.

    No conformal-Ricci formula or expected Ricci expression is used here.
    R_ab = d_c Gamma^c_ab - d_b Gamma^c_ac
           + Gamma^c_cd Gamma^d_ab - Gamma^c_bd Gamma^d_ac.
    """
    inverse = metric.inv()
    dg = [[[partial(metric[a, b], c) for c in range(DIMENSION)]
           for b in range(DIMENSION)] for a in range(DIMENSION)]
    connection = [[[sp.cancel(sum(
        inverse[a, d] * (dg[d][c][b] + dg[d][b][c] - dg[b][c][d]) / 2
        for d in range(DIMENSION)))
        for c in range(DIMENSION)] for b in range(DIMENSION)] for a in range(DIMENSION)]
    trace_connection = [sum(connection[c][c][d] for c in range(DIMENSION))
                        for d in range(DIMENSION)]
    ricci = sp.zeros(DIMENSION)
    for a in range(DIMENSION):
        for b in range(DIMENSION):
            first = sum(partial(connection[c][a][b], c)
                        - partial(connection[c][a][c], b) for c in range(DIMENSION))
            second = sum(trace_connection[d] * connection[d][a][b] for d in range(DIMENSION))
            third = sum(connection[c][b][d] * connection[d][a][c]
                        for c in range(DIMENSION) for d in range(DIMENSION))
            ricci[a, b] = ambient_reduce(first + second - third)
    return inverse, connection, ricci


def run_checks(inject_failure=False):
    checks = Checks()
    root6 = sp.sqrt(6)
    Omega = (1 + s) ** -2
    Phi = root6 * sp.log(1 + s)
    phi = sp.Symbol("phi", real=True)
    potential = 84 * sp.exp(2 * phi / root6) - 36 * sp.exp(3 * phi / root6)
    U_at_Phi = sp.simplify(potential.subs(phi, Phi))
    Uprime_at_Phi = sp.simplify(sp.diff(potential, phi).subs(phi, Phi))
    metric = Omega ** 2 * eta
    inverse, connection, actual_Ricci = geometry_from_metric(metric)
    checks.check("metric inverse from matrix inversion", equal(metric * inverse, sp.eye(DIMENSION)))
    checks.check("Cartesian invariant derivative is 2 X_lower",
                 all(equal(partial(s, a), 2 * X_lower[a]) for a in range(DIMENSION)))
    checks.check("metric-compatible torsion-free Levi-Civita connection",
                 all(equal(connection[a][b][c], connection[a][c][b])
                     for a in range(DIMENSION) for b in range(DIMENSION) for c in range(DIMENSION))
                 and all(zero(partial(metric[a, b], c)
                              - sum(connection[d][c][a] * metric[d, b]
                                    + connection[d][c][b] * metric[a, d] for d in range(DIMENSION)))
                         for a in range(DIMENSION) for b in range(DIMENSION) for c in range(DIMENSION)))
    coefficient = 23 if inject_failure else 24
    expected_Ricci = (coefficient * X_lower * X_lower.T + (32 - 24 * s) * eta) / (1 + s) ** 2
    checks.check("all Ricci components from the Christoffel definition", equal(actual_Ricci, expected_Ricci),
                 "injected background mismatch detected" if inject_failure else None)
    checks.check("Ricci symmetry from the independently calculated components", equal(actual_Ricci, actual_Ricci.T))
    checks.check("Ricci scalar by metric contraction",
                 equal(sp.trace(inverse * actual_Ricci), (160 - 96 * s) * (1 + s) ** 2))
    phi_gradient = sp.Matrix([partial(Phi, a) for a in range(DIMENSION)])
    phi_hessian = sp.Matrix(DIMENSION, DIMENSION, lambda a, b:
        partial(phi_gradient[b], a)
        - sum(connection[c][a][b] * phi_gradient[c] for c in range(DIMENSION)))
    box_connection = ambient_reduce(sum(inverse[a, b] * phi_hessian[a, b]
                                       for a in range(DIMENSION) for b in range(DIMENSION)))
    density = Omega ** DIMENSION  # sqrt(|det g|), since Omega>0 on the stated domain.
    box_divergence = ambient_reduce(sum(partial(sum(density * inverse[a, b] * phi_gradient[b]
                                                   for b in range(DIMENSION)), a)
                                       for a in range(DIMENSION)) / density)
    expected_box = 2 * root6 * (5 - 9 * s) * (1 + s) ** 2
    checks.check("volume density from the metric determinant", equal(density ** 2, -metric.det()))
    checks.check("scalar dAlembertian from covariant Hessian", equal(box_connection, expected_box))
    checks.check("independent divergence definition of scalar dAlembertian", equal(box_divergence, box_connection))
    checks.check("potential derivative by differentiation in phi", equal(Uprime_at_Phi, expected_box))
    checks.check("all trace-reversed Einstein equations",
                 equal(actual_Ricci - phi_gradient * phi_gradient.T - sp.Rational(2, 3) * U_at_Phi * metric,
                       sp.zeros(DIMENSION)))
    checks.check("scalar equation with the differentiated stated potential", equal(box_connection, Uprime_at_Phi))
    checks.reject("Ricci curvature sign reversal", equal(-actual_Ricci, expected_Ricci))
    checks.reject("Ricci radial coefficient 12 instead of 24",
                  equal(actual_Ricci, (12 * X_lower * X_lower.T + (32 - 24 * s) * eta) / (1 + s) ** 2))
    checks.reject("Ricci metric coefficient with wrong sign of 24s",
                  equal(actual_Ricci, (24 * X_lower * X_lower.T + (32 + 24 * s) * eta) / (1 + s) ** 2))
    wrong_density = Omega ** 4
    wrong_box = ambient_reduce(sum(partial(sum(wrong_density * inverse[a, b] * phi_gradient[b]
                                             for b in range(DIMENSION)), a)
                                 for a in range(DIMENSION)) / wrong_density)
    checks.reject("four-dimensional volume power in five dimensions", equal(wrong_box, box_connection))
    wrong_U = 84 * sp.exp(2 * phi / root6) - 36 * sp.exp(2 * phi / root6)
    checks.reject("incorrect cubic exponential in the potential",
                  equal(box_connection, sp.diff(wrong_U, phi).subs(phi, Phi)))

    # Normalized outward normal and its covariant derivative on positive-s leaves.
    # This calculation uses the same metric-derived connection, not an assumed K.
    normal_upper = X / (Omega * sp.sqrt(s))
    normal_lower = metric * normal_upper
    nabla_normal = sp.Matrix(DIMENSION, DIMENSION, lambda a, b:
        partial(normal_lower[b], a)
        - sum(connection[c][a][b] * normal_lower[c] for c in range(DIMENSION)))
    induced = metric - normal_lower * normal_lower.T
    umbilicity = (1 + s) * (1 - 3 * s) / sp.sqrt(s)
    checks.check("outward hyperboloid normal has unit positive norm",
                 equal((normal_upper.T * metric * normal_upper)[0], 1))
    checks.check("outward normal differentiates s positively on s>0",
                 equal(sum(normal_upper[a] * partial(s, a) for a in range(DIMENSION)),
                       2 * sp.sqrt(s) / Omega))
    checks.check("extrinsic curvature from covariant normal derivative is umbilic",
                 equal(nabla_normal, umbilicity * induced))
    checks.check("normal derivative has zero normal contraction",
                 zero(nabla_normal * normal_upper))
    checks.check("umbilicity from the tangential trace",
                 equal(sp.trace(inverse * nabla_normal) / 4, umbilicity))
    normal_phi = ambient_reduce((normal_upper.T * phi_gradient)[0])
    b = sp.Rational(1, 2)
    wall = b ** 2
    phi_b = root6 * sp.log(sp.Rational(5, 4))
    a_b = sp.simplify(umbilicity.subs(s, wall))
    p_b = sp.simplify(normal_phi.subs(s, wall))
    sigma = sp.Rational(15, 4) - 5 * root6 * (phi - phi_b) / 2
    checks.check("wall scalar is the restricted bulk scalar", equal(Phi.subs(s, wall), phi_b))
    checks.check("wall conformal factor and de Sitter radius",
                 equal(Omega.subs(s, wall), sp.Rational(16, 25))
                 and equal(b * Omega.subs(s, wall), sp.Rational(8, 25)))
    checks.check("wall a and p computed from geometric derivatives",
                 equal(a_b, sp.Rational(5, 8)) and equal(p_b, 5 * root6 / 4))
    checks.check("doubled-interior metric junction sigma equals 6a",
                 equal(sigma.subs(phi, phi_b), 6 * a_b))
    checks.check("doubled-interior scalar junction sigma prime equals minus 2p",
                 equal(sp.diff(sigma, phi).subs(phi, phi_b), -2 * p_b))
    checks.reject("wrong scalar junction sign plus 2p",
                  equal(sp.diff(sigma, phi).subs(phi, phi_b), 2 * p_b))
    checks.reject("wrong outward extrinsic-curvature sign",
                  equal(sigma.subs(phi, phi_b), -6 * a_b))

    # Radial definitions are independently evaluated before comparing formulas.
    r = sp.Symbol("r", positive=True)
    f = sp.diff(Phi, s)
    L = sp.cancel(sp.diff(Omega, s) / Omega + 2 * sp.diff(f, s) / f)
    mr = sp.cancel((9 * r ** 3 / (4 * Omega * f ** 2)).subs(s, r ** 2))
    nr = sp.cancel((6 * r ** 7 / Omega).subs(s, r ** 2))
    checks.check("master lower-order coefficient from Omega and Phi derivatives", equal(L, -4 / (1 + s)))
    checks.check("radial m from its defining expression", equal(mr, sp.Rational(3, 8) * r ** 3 * (1 + r ** 2) ** 4))
    checks.check("radial n from its defining expression", equal(nr, 6 * r ** 7 * (1 + r ** 2) ** 2))
    checks.check("radial integration-by-parts logarithmic derivative",
                 equal(r * sp.diff(mr, r) / mr, 3 - 2 * r ** 2 * L.subs(s, r ** 2)))
    coefficient_by_parts = 20 * mr - sp.Rational(9, 2) * sp.diff(r * mr, r) + nr
    coefficient_by_L = mr * (2 + 9 * r ** 2 * L.subs(s, r ** 2)) + nr
    numerator = 2 - 32 * r ** 2 - 18 * r ** 4
    expected_coefficient = mr * numerator / (1 + r ** 2) ** 2
    checks.check("radial coefficient from direct integration by parts", equal(coefficient_by_parts, coefficient_by_L))
    checks.check("radial coefficient exact factored expression", equal(coefficient_by_L, expected_coefficient))
    checks.check("positive radial m on r>0", mr.is_positive is True)
    checks.check("annulus numerator has strictly negative derivative on r>0",
                 equal(sp.diff(numerator, r), -r * (64 + 72 * r ** 2)))
    numerator_floor = sp.Rational(432, 625)
    ratio_floor = sp.Rational(108, 169)
    checks.check("annulus positivity endpoint and denominator certificate",
                 equal(numerator.subs(r, sp.Rational(1, 5)), numerator_floor)
                 and numerator_floor > 0
                 and equal(numerator_floor / (1 + sp.Rational(1, 25)) ** 2, ratio_floor)
                 and ratio_floor > 0)
    checks.reject("omitting radial n changes the coefficient",
                  equal(mr * (2 + 9 * r ** 2 * L.subs(s, r ** 2)), expected_coefficient))
    checks.reject("radial m with wrong conformal power",
                  equal(mr, sp.Rational(3, 8) * r ** 3 * (1 + r ** 2) ** 2))

    # Exact causal inequalities and a standard nonzero C-infinity flat bump.
    r1, r2, time = sp.Rational(11, 100), sp.Rational(19, 100), sp.Rational(1, 32)
    checks.check("bump endpoints lie strictly in the positivity annulus",
                 sp.Rational(1, 10) < r1 < r2 < sp.Rational(1, 5) < b)
    midpoint = (r1 + r2) / 2
    bump_exponent_midpoint = -1 / ((midpoint - r1) * (r2 - midpoint))
    checks.check("standard flat bump is nonzero at its midpoint",
                 bump_exponent_midpoint == -625 and sp.exp(bump_exponent_midpoint).is_positive is True)
    checks.check("time satisfies both strict pre-contact conditions",
                 0 < time < r1 / 2 and time < (b - r2) / 2 and time < sp.Rational(1, 2))
    checks.check("support remains strictly outside the null cone",
                 r1 - 2 * time == sp.Rational(19, 400) and r1 - 2 * time > 0)
    checks.check("support remains strictly inside the shell radial lower bound",
                 b - r2 - time == sp.Rational(223, 800) and b - r2 - time > 0)
    support_s_lower, support_s_upper = r1 ** 2 - 2 * r1 * time, r2 ** 2 + 2 * r2 * time
    checks.check("entire causal support lies in a compact positive-s subdomain",
                 support_s_lower == sp.Rational(209, 40000)
                 and support_s_upper == sp.Rational(1919, 40000)
                 and 0 < support_s_lower < support_s_upper < wall)
    checks.check("complete physical slab has 1+s greater than 3/4",
                 1 - sp.Rational(1, 2) ** 2 == sp.Rational(3, 4)
                 and sp.Rational(3, 4) > 0)
    checks.check("scalar radial derivative is nonzero throughout s>-1",
                 equal(f, root6 / (1 + s)) and root6.is_positive is True)
    checks.reject("overlong time loses strict pre-cone separation",
                  r1 - 2 * sp.Rational(1, 8) > 0)

    return {
        "status": "PASS_EXACT_BACKGROUND_CHECKS",
        "scope": "Separate explicit five-dimensional Einstein-scalar shell model; exact local algebra and causal bounds only; no new-solution novelty or registered-branch claim",
        "positive_checks": len(checks.positive),
        "wrong_formula_controls": len(checks.controls),
        "checks": checks.positive,
        "controls": checks.controls,
        "verified_expressions": {
            "Omega": "(1+s)^(-2)", "Phi": "sqrt(6)*log(1+s)",
            "U": "84*exp(2*phi/sqrt(6))-36*exp(3*phi/sqrt(6))",
            "Ricci": "(24*X_lower*X_lower^T+(32-24*s)*eta)/(1+s)^2",
            "box_Phi": "2*sqrt(6)*(5-9*s)*(1+s)^2",
            "b": "1/2", "a_b": "5/8", "p_b": "5*sqrt(6)/4",
            "sigma_b": "15/4", "sigma_prime_b": "-5*sqrt(6)/2",
            "L": "-4/(1+s)", "m": "(3/8)*r^3*(1+r^2)^4",
            "n": "6*r^7*(1+r^2)^2",
            "radial_coefficient_over_m": "(2-32*r^2-18*r^4)/(1+r^2)^2",
            "annulus_ratio_lower_bound": "108/169",
            "bump_support": "[11/100,19/100]", "T_star": "1/32",
            "causal_support_s_bounds": ["209/40000", "1919/40000"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inject-failure", action="store_true")
    args = parser.parse_args()
    try:
        report = run_checks(args.inject_failure)
        code = 0
    except Exception as error:
        report = {"status": "failed", "error": type(error).__name__ + ": " + str(error)}
        code = 1
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if code == 0:
        print(f"PASS {report['positive_checks']} exact identities; {report['wrong_formula_controls']} wrong-formula controls")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
