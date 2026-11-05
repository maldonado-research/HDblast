#!/usr/bin/env python3
"""Exact checks of selected identities in ODD_REMAINDER.md; no dependencies.

These finite algebra checks do not replace the proof and do not estimate a
physical observable. In particular all constants below are synthetic inputs.
"""

import argparse
import json
from fractions import Fraction as F
from pathlib import Path


CHECK_COUNT = 0


def check(condition):
    """A live finite guard, also when invoked with python -O."""
    global CHECK_COUNT
    if not condition:
        raise RuntimeError("exact arithmetic guard failed")
    CHECK_COUNT += 1


def triangular_monomial(power, h=F(1)):
    """h^-1 int_0^h (h-s)s^power ds, integer power."""
    return h ** (power + 1) / ((power + 1) * (power + 2))


def mul(p, q):
    """Bivariate polynomials in e and y, with exact rational coefficients."""
    out = {}
    for (i, j), v in p.items():
        for (k, l), w in q.items():
            key = i + k, j + l
            out[key] = out.get(key, F(0)) + v * w
    return {key: value for key, value in out.items() if value}


def add(p, q):
    out = p.copy()
    for key, value in q.items():
        out[key] = out.get(key, F(0)) + value
    return {key: value for key, value in out.items() if value}


def sub(p, q):
    return add(p, {key: -value for key, value in q.items()})


def evaluate(p, e, y=F(0)):
    return sum((value * e**i * y**j for (i, j), value in p.items()), F(0))


def mean_s3(p):
    # Only these Y moments occur below. Reflection/quarter-turn symmetry gives
    # the odd moments, and <Y^2>=1/6 is the supplied round-S3 identity.
    moments = {0: F(1), 1: F(0), 2: F(1, 6), 3: F(0)}
    out = {}
    for (i, j), value in p.items():
        check(j in moments)
        key = i, 0
        out[key] = out.get(key, F(0)) + value * moments[j]
    return {key: value for key, value in out.items() if value}


def equivariant(p):
    # Joint transform e -> -e, Y -> -Y leaves an equivariant polynomial fixed.
    return all((i + j) % 2 == 0 for i, j in p)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="absolute destination of the canonical JSON report")
    args = parser.parse_args()
    if not args.output.is_absolute():
        parser.error("--output must be an absolute path")
    groups = []
    for alpha in (F(0), F(1, 4), F(1, 2), F(3, 4), F(1)):
        coeff = F(1) / (alpha + 1) - F(1) / (alpha + 2)
        check(coeff == F(1) / ((alpha + 1) * (alpha + 2)))
        # Differentiating |e|^(alpha+2)/denominator twice gives |e|^alpha.
        check(coeff * (alpha + 2) * (alpha + 1) == 1)
    groups.append("Holder integral coefficients: alpha=0,1/4,1/2,3/4,1")

    # Exact all-sign evaluation of the sharp anchored example. Set h=t^q
    # so h^(alpha+1) is rational when alpha=p/q. Skip alpha=0 as a C2 model.
    S, L, t = F(-7, 3), F(5, 2), F(2, 5)
    for alpha in (F(1, 4), F(1, 2), F(3, 4), F(1)):
        p, q = alpha.numerator, alpha.denominator
        h = t**q
        err = L * t ** (p + q) / ((alpha + 1) * (alpha + 2))
        for sign in (-1, 1):
            e = sign * h
            value = e * S + sign * L * t ** (p + 2 * q) / ((alpha + 1) * (alpha + 2))
            check(value / e - S == err)
            second_derivative = sign * L * t**p
            check(second_derivative * sign == L * t**p)
    groups.append("both parameter signs: normalized remainder is even")

    # Exact boundary-layer integral at h=1. On [0,delta], M''=L*s/delta;
    # on [delta,1], M''=L. Integrate the two polynomial pieces independently.
    for k in (2, 3, 10, 100):
        delta = F(1, k)
        first = L / delta * (delta**2 / 2 - delta**3 / 3)
        second = L * (F(1, 2) - delta + delta**2 / 2)
        value = first + second
        check(value == L * (F(1, 2) - delta / 2 + delta**2 / 6))
        check(0 < L / 2 - value < L * delta / 2)
    groups.append("alpha=0: C2 boundary layers approach, but do not attain, L/2")

    # Check the integrated raw residual coefficients, including both area terms.
    check(3 * triangular_monomial(1) == F(1, 2))
    check(F(1, 2) * triangular_monomial(2) == F(1, 24))
    for k in (1, 2, 3):
        check(triangular_monomial(k) == F(1, (k + 1) * (k + 2)))
    groups.append("raw residual coefficients: E2*J1/2 and E2*J2/24")

    # Algebraic angular example with a genuinely varying even denominator:
    # A=e*S*Y+C*e^2/2, j=1+b*e*Y+g*e^2. It obeys joint equivariance.
    # Adding g lets the checker exercise the full residual quotient identity.
    S, C, b, g = F(7), F(6), F(3, 5), F(2, 3)
    A = {(1, 1): S, (2, 0): C / 2}
    j = {(0, 0): F(1), (1, 1): b, (2, 0): g}
    Y, Y2 = {(0, 1): F(1)}, {(0, 2): F(1)}
    n = {key: value * 6 for key, value in mean_s3(mul(Y, mul(A, j))).items()}
    d = {key: value * 6 for key, value in mean_s3(mul(Y2, j)).items()}
    H = mul(sub(A, {(1, 1): S}), j)
    q = {key: value * 6 for key, value in mean_s3(mul(Y, H)).items()}
    check(equivariant(A) and equivariant(j) and equivariant(H))
    check(n == {(1, 0): S, (3, 0): S * g + C * b / 2})
    check(d == {(0, 0): F(1), (2, 0): g})
    check(q == {(3, 0): C * b / 2})
    check(q == sub(n, mul({(1, 0): S}, d)))
    for e in (F(-1, 8), F(1, 8)):
        nv, dv, qv = evaluate(n, e), evaluate(d, e), evaluate(q, e)
        error = nv / (e * dv) - S
        check(error == qv / (e * dv) == C * b * e**2 / (2 * (1 + g * e**2)))
        check(1 - abs(b * e) + g * e**2 > 0)  # sufficient j>0 for |Y|<=1
    # If the nonlinear scalar residual is absent, varying area cancels exactly.
    A_linear = {(1, 1): S}
    n_linear = {key: value * 6 for key, value in
                mean_s3(mul(Y, mul(A_linear, j))).items()}
    check(n_linear == mul({(1, 0): S}, d))
    groups.append("angular/area example: N odd, D even; coupling survives projection")
    groups.append("quotient residual preserves exact A=e*S*Y cancellation")
    report = {
        "schema_version": 1,
        "status": "pass",
        "finite_exact_guard_count": CHECK_COUNT,
        "checks": groups,
        "method": "Python standard-library Fraction arithmetic with live guards",
        "scope": "Finite algebra checks only; no PDE estimates or physical constants; not a formal proof",
    }
    canonical = json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n"
    args.output.write_text(canonical, encoding="utf-8")
    print(canonical, end="")


if __name__ == "__main__":
    main()
