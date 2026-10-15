#!/usr/bin/env python3
"""Exact arithmetic for the independent uniform waveform certificate.

No archived producer is imported. No floating value participates in a decision.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / "inputs" / "previous"
PINS = {
    "SHELL_RESPONSE_ERROR_CERTIFICATE.md": "5f01be536c77123e9f667f4cd5514b52067641b2b389ba4e6eb655cb24654022",
    "sol/RESIDUAL_TO_PHYSICAL_RESPONSE.md": "b8c77c810c70aa286227c5dabe010a53942af84067eb03961f0186419a5ea96a",
    "archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md": "0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c",
    "archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
}


def rat(value: F | int) -> str:
    value = F(value)
    return f"{value.numerator}/{value.denominator}"


def positive_exp_upper(x: F, n: int = 16) -> F:
    """Positive Taylor sum plus rigorous geometric tail, for x>=0."""
    if not isinstance(x, F) or x < 0 or n < 0:
        raise ValueError("nonnegative rational and nonnegative order required")
    term = F(1)
    total = term
    for k in range(1, n + 1):
        term = term * x / k
        total += term
    next_term = term * x / (n + 1)
    ratio = x / (n + 2)
    if ratio >= 1:
        raise ValueError("geometric tail ratio must be less than one")
    return total + next_term / (1 - ratio)


def run(inject_error: str | None = None) -> dict:
    checks = []

    def check(name: str, predicate: bool) -> None:
        if predicate is not True:
            raise RuntimeError(f"failed check: {name}")
        checks.append(name)

    sources = {}
    for relative, expected in PINS.items():
        path = FROZEN / relative
        data = path.read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        check(f"pinned source {relative}", actual == expected)
        sources[relative] = {"sha256": actual, "bytes": len(data)}

    ell = F(1, 10**10)
    T = 3 * ell
    L = 8 * ell
    m = F(1, 1000)
    M, B, Gamma = F(10**6), F(1000), F(250)
    check("positive inertia lower bound", m > 0)
    check("causal initial line margin", L - T == 5 * ell > 0)
    check("initial support has positive inner margin", ell > 0)
    check("initial support has positive outer margin", L - 2 * ell == 6 * ell > 0)
    check("strict evolution margin", 7 * ell / 2 - T == ell / 2 > 0)
    check("strict nonlinear correction margin", L - (3 * ell + F(5, 4) * 7 * ell / 2) == 5 * ell / 8 > 0)

    a = 1 + T
    r = M * T**2 / 2
    s = B * T**2 / (2 * m)
    det = (1 - r) * (1 - s) - T * r / m
    if inject_error == "drop-determinant-coupling":
        det = (1 - r) * (1 - s)
    check("r equals displayed exact value", r == F(9, 200000000000000))
    check("s equals displayed exact value", s == r)
    check("positive matrix diagonal", 1 - r > 0 and 1 - s > 0)
    check("positive coupled determinant", det > 0)
    A = a * (1 - s + T / m) / det
    V = a / (m * det)
    check("first majorant equation", (1 - r) * A - T * V == a)
    check("second majorant equation", (1 - s) * V - (r / m) * A == a / m)
    check("positive inverse entries", all(x > 0 for x in ((1 - s) / det, T / det, r / (m * det), (1 - r) / det)))
    check("displayed A upper bound", A < F(100000031, 100000000))
    check("displayed V upper bound", V < F(1000000001, 1000000))

    eta = ell * (M * T * A + (1 + Gamma * T) * V)
    eta_display = F(10000004, 10**14)
    check("strict global uniform output bound", eta < eta_display)
    V_tail = ell / m**2 + ell / m + (r / m) * A + s * V
    eta_tail = ell * (M * T * A + V_tail + Gamma * T * V)
    eta_tail_display = F(4752, 10**17)
    check("strict tail uniform output bound", eta_tail < eta_tail_display)
    check("tail bound tighter than global bound", 0 < eta_tail < eta)
    zeta = F(10000005, 10**14)
    check("support-aware observation norm bound", eta**2 + eta_tail**2 < zeta**2)

    # Exact boundary jet at t=1: mu=3,beta=2,gamma=5,
    # v=t^2, q=t^3/3, incoming J=mu*v'+v+beta*q.
    q, v, vp, mu, beta, gamma = F(1, 3), F(1), F(2), F(3), F(2), F(5)
    J = mu * vp + v + beta * q
    if inject_error == "flip-boundary-damping":
        J = mu * vp - v + beta * q
    wz = J - v
    D = wz - gamma * q
    if inject_error == "use-full-C":
        D = wz - 2 * gamma * q
    check("boundary jet correct sign identity", -wz + mu * vp + beta * q == 0)
    check("wrong damping sign rejected", mu * vp - v + beta * q != J)
    check("wrong beta sign rejected", mu * vp + v - beta * q != J)
    check("wrong incoming reflection rejected", J + v != wz)
    check("wrong observation gamma sign rejected", wz + gamma * q != D)
    check("wrong gamma equals full C rejected", wz - 2 * gamma * q != D)

    wrong_det = (1 - r) * (1 - s)
    wrong_A = a * (1 - s + T / m) / wrong_det
    wrong_V = a / (m * wrong_det)
    check("uncoupled determinant fails first supersolution", (1 - r) * wrong_A - T * wrong_V < a)
    check("uncoupled determinant fails second supersolution", (1 - s) * wrong_V - (r / m) * wrong_A < a / m)

    # At x=3/2,Q=1,Q'=0, so dropping the prescribed initial velocity
    # would incorrectly replace positive ell*exp(-3ell/2) by zero.
    velocity_lower = ell * (1 - 3 * ell / 2)
    if inject_error == "omit-velocity-input":
        velocity_lower = F(0)
    check("dropped velocity input rejected at peak", velocity_lower > 0)

    # For f(s)=s(1-s) on [0,1], at t=2,mu=1 the derivative forcing
    # equals -integral exp(-(2-s))*f(s) ds.  e^2<8 proves its
    # magnitude > (1/8)*integral f=1/48, despite f(2)=0.
    exp2_upper = positive_exp_upper(F(2))
    check("rigorous exponential diagnostic", exp2_upper < 8)
    memory_magnitude_lower = F(1, 6) / exp2_upper
    if inject_error == "delete-boundary-memory":
        memory_magnitude_lower = F(0)
    check("deleted boundary memory rejected", memory_magnitude_lower > F(1, 48))

    own = {}
    for name in ("UNIFORM_WAVEFORM_CERTIFICATE.md", "check_uniform_waveform.py"):
        data = (HERE / name).read_bytes()
        own[name] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}

    return {
        "schema": "hdblast-sol-uniform-waveform-v1",
        "status": "PASS",
        "scope": "Exact-rational arithmetic conditional on the written characteristic/kernel proof and retained analytic hypotheses; no formal PDE verification.",
        "decision_arithmetic": "fractions.Fraction; no floating-point inputs or calculations",
        "python_optimization": sys.flags.optimize,
        "sources": sources,
        "owned_inputs": own,
        "parameters": {k: rat(x) for k, x in {"ell": ell, "T": T, "L": L, "m": m, "M": M, "B": B, "Gamma": Gamma}.items()},
        "exact": {k: rat(x) for k, x in {"a": a, "r": r, "s": s, "determinant": det, "A_upper": A, "V_upper": V, "eta": eta, "eta_display_upper": eta_display, "V_tail_upper": V_tail, "eta_tail": eta_tail, "eta_tail_display_upper": eta_tail_display, "support_aware_L2_error_upper": zeta, "velocity_input_lower_at_peak": velocity_lower, "diagnostic_exp2_upper": exp2_upper, "diagnostic_memory_magnitude_lower": memory_magnitude_lower}.items()},
        "checks": checks,
        "check_count": len(checks),
        "unproved_by_program": ["retained actual-background theorem", "coupled reconstruction and physical gauge bridge", "nonlinear exact-data/evolution hypotheses", "analytic characteristic and positive-kernel estimates", "physical observational scale or effective nonlinear epsilon"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="UNIFORM_WAVEFORM_CHECKS.json")
    parser.add_argument("--inject-error", choices=("drop-determinant-coupling", "flip-boundary-damping", "use-full-C", "omit-velocity-input", "delete-boundary-memory"))
    args = parser.parse_args()
    filename = Path(args.output)
    if filename.name != args.output:
        raise ValueError("output must be a filename within the owned Sol directory")
    result = run(args.inject_error)
    (HERE / filename).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"PASS: {result['check_count']} exact checks; output {filename.name}")


if __name__ == "__main__":
    main()
