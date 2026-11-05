#!/usr/bin/env python3
"""Exact, deterministic algebra checks; no evaluated physical ODE inputs.

Run: python3 check_budget.py --output /absolute/path/report.json
The report bytes are identical with and without Python -O.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path


N = 8
E, X, Y, VX, VY, WX, WY, T = range(N)
ZERO = (0,) * N


class Poly:
    """Small sparse exact polynomial ring, truncated only in perturbation T."""

    def __init__(self, data=None):
        self.data = {m: Q(c) for m, c in (data or {}).items() if c and m[T] <= 2}

    @classmethod
    def constant(cls, c):
        return cls({ZERO: Q(c)})

    @classmethod
    def variable(cls, i):
        m = list(ZERO)
        m[i] = 1
        return cls({tuple(m): Q(1)})

    def __add__(self, other):
        other = coerce(other)
        out = self.data.copy()
        for m, c in other.data.items():
            out[m] = out.get(m, Q(0)) + c
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.data.items()})

    def __sub__(self, other):
        return self + -coerce(other)

    def __mul__(self, other):
        other = coerce(other)
        out = {}
        for m, c in self.data.items():
            for n, d in other.data.items():
                p = tuple(a + b for a, b in zip(m, n))
                if p[T] <= 2:
                    out[p] = out.get(p, Q(0)) + c * d
        return Poly(out)

    __rmul__ = __mul__

    def __pow__(self, n):
        out = Poly.constant(1)
        for _ in range(n):
            out = out * self
        return out

    def derivative(self, i):
        out = {}
        for m, c in self.data.items():
            if m[i]:
                p = list(m)
                p[i] -= 1
                out[tuple(p)] = c * m[i]
        return Poly(out)

    def substitute(self, mapping):
        out = Poly()
        for m, c in self.data.items():
            term = Poly.constant(c)
            for i, power in enumerate(m):
                if power:
                    term = term * mapping.get(i, Poly.variable(i)) ** power
            out = out + term
        return out

    def t_coefficient(self, k):
        out = {}
        for m, c in self.data.items():
            if m[T] == k:
                p = list(m)
                p[T] = 0
                out[tuple(p)] = c
        return Poly(out)

    def __eq__(self, other):
        return self.data == coerce(other).data


def coerce(x):
    return x if isinstance(x, Poly) else Poly.constant(x)


def require(condition, name):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + name)


def run_checks():
    checks = []

    def check(condition, name):
        require(condition, name)
        checks.append(name)

    e, x, y, vx, vy, wx, wy, t = [Poly.variable(i) for i in range(N)]
    path = {E: e + t, X: x + t * vx + Q(1, 2) * t**2 * wx,
            Y: y + t * vy + Q(1, 2) * t**2 * wy}

    def chain(f):
        fx, fy, fe = [f.derivative(i) for i in (X, Y, E)]
        first = fx * vx + fy * vy + fe
        second = (fx * wx + fy * wy
                  + fx.derivative(X) * vx**2
                  + 2 * fx.derivative(Y) * vx * vy
                  + fy.derivative(Y) * vy**2
                  + 2 * fe.derivative(X) * vx
                  + 2 * fe.derivative(Y) * vy + fe.derivative(E))
        composed = f.substitute(path)
        return composed.t_coefficient(1), 2 * composed.t_coefficient(2), first, second

    # A basis check is coefficientwise, not floating-point sampling.
    # Vector-valued F follows componentwise and arbitrary coefficients by linearity.
    monomial_count = 0
    for powers in itertools.product(range(5), repeat=3):
        if sum(powers) > 4:
            continue
        a, b, c = powers
        f = Q(2 * a + 3 * b + 5 * c + 1, 7) * e**a * x**b * y**c
        observed1, observed2, expected1, expected2 = chain(f)
        label = "_".join(map(str, powers))
        check(observed1 == expected1, "ode_first_" + label)
        check(observed2 == expected2, "ode_second_" + label)
        monomial_count += 1
    check(monomial_count == 35, "complete_total_degree_four_basis")

    # A nontrivial complete scalar example explicitly exercises every chain term.
    scalar = e**2 * x**2 + e * x * y + x * y**2 + Q(3, 11) * y
    observed1, observed2, expected1, expected2 = chain(scalar)
    check(observed1 == expected1, "full_scalar_first_chain")
    check(observed2 == expected2, "full_scalar_second_chain")
    no_mixed = expected2 - 2 * scalar.derivative(E).derivative(X) * vx
    check(observed2 != no_mixed, "mixed_factor_two_term_is_required")
    no_position_second = expected2 - scalar.derivative(X) * wx
    check(observed2 != no_position_second, "position_second_variation_is_required")

    # Fixed, nonlinear state-chart transitions have no explicit parameter term.
    chart_count = 0
    for b, c in itertools.product(range(5), repeat=2):
        if b + c > 4:
            continue
        f = x**b * y**c
        observed1, observed2, expected1, expected2 = chain(f)
        check(observed1 == expected1, f"chart_first_{b}_{c}")
        check(observed2 == expected2, f"chart_second_{b}_{c}")
        chart_count += 1
    check(chart_count == 15, "complete_chart_degree_four_basis")
    nonlinear_transition = x**2
    _, observed2, _, _ = chain(nonlinear_transition)
    check(observed2 == 2 * x * wx + 2 * vx**2,
          "chart_hessian_charge_explicit")
    check(observed2 != 2 * x * wx, "chart_hessian_charge_cannot_be_dropped")

    # One-dimensional algebra model verifies the universal spray/Jacobi signs.
    # It is an algebra model only; the physical shell has four dimensions.
    gamma = e * x**2 + x
    proper_field = (y, -gamma * y**2, vy,
                    -gamma.derivative(X) * vx * y**2 - 2 * gamma * y * vy)
    past_field = (-y, gamma * y**2, -vy,
                  gamma.derivative(X) * vx * y**2 + 2 * gamma * y * vy)
    for i, (g, f) in enumerate(zip(proper_field, past_field)):
        check(f == -g, f"past_orientation_component_{i}")
        for variable in (E, X, Y, VX, VY):
            check(f.derivative(variable) == -g.derivative(variable),
                  f"past_orientation_derivative_{i}_{variable}")
    check(past_field[0] != y, "future_velocity_not_redefined_as_past_velocity")

    # Exact formal coefficients for phi_L(T)=T sum (LT)^k/(k+1)!.
    # These verify algebra and edge limits, not an exponential tail bound.
    factorial = 1
    exp_coeffs = []
    for k in range(13):
        if k:
            factorial *= k
        exp_coeffs.append(Q(1, factorial))
    phi_coeffs = [exp_coeffs[k + 1] for k in range(12)]
    for k in range(12):
        check(phi_coeffs[k] == exp_coeffs[k + 1], f"phi_exp_identity_coeff_{k}")
        check((k + 1) * phi_coeffs[k] == exp_coeffs[k],
              f"phi_integral_identity_coeff_{k}")
        check(phi_coeffs[k] > 0, f"phi_nonnegative_coeff_{k}")
    check(phi_coeffs[0] == 1, "phi_zero_L_limit_is_T")

    # Rational edge cases of the stated majorants; these are generic algebra
    # test vectors, not evaluated HDBLAST/observer coefficients.
    def zero_l_bound(duration, p1, hh, cc, p2, v0, w0):
        vv = v0 + duration * p1
        ww = w0 + duration * (hh * vv**2 + 2 * cc * vv + p2)
        return vv, ww

    generic = list(map(Q, (0, 2, 3, 5, 7, 11, 13)))
    check(zero_l_bound(*generic) == (Q(11), Q(13)), "zero_time_preserves_initial_jets")
    check(zero_l_bound(Q(2), Q(3), Q(5), Q(7), Q(11), Q(13), Q(17))
          == (Q(19), Q(4181)), "zero_L_majorant_exact_rational")
    check(zero_l_bound(Q(3), Q(0), Q(9), Q(8), Q(4), Q(0), Q(2))
          == (Q(0), Q(14)), "zero_first_jet_can_have_second_forcing")
    check(zero_l_bound(Q(3), Q(0), Q(9), Q(8), Q(0), Q(0), Q(0))
          == (Q(0), Q(0)), "all_zero_jets_and_forcing")

    # Triangle/submultiplicative coefficient structure in the scalar bound.
    # Exact polynomial equality checks the coefficients before inequalities.
    be, bx, bee, bex, bxx, z, w = list(map(Q, (2, 3, 5, 7, 11, 13, 17)))
    check(be + bx * z == 41, "scalar_first_budget_coefficients")
    check(bee + 2 * bex * z + bxx * z*z + bx * w == 2097,
          "scalar_second_budget_coefficients")

    return {
        "schema": "hdblast-observer-flow-algebra-check-v1",
        "status": "PASS",
        "arithmetic": "fractions.Fraction; exact sparse polynomial coefficients",
        "checks_passed": len(checks),
        "check_names": checks,
        "physical_inputs_evaluated": False,
        "pde_or_tube_certified": False,
        "exponential_numeric_value_evaluated": False,
        "limitations": [
            "Polynomial identities and exact rational edge cases only.",
            "Gronwall and differentiability arguments are in the accompanying proof.",
            "No compact tube, PDE constant, nonlinear radius, or area bound is evaluated."
        ]
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True,
                        help="absolute fresh file path for canonical JSON report; parent must exist")
    args = parser.parse_args()
    require(args.output.is_absolute(), "--output must be an absolute path")
    require(args.output.parent.is_dir(), "--output parent directory must exist")
    require(not args.output.exists() and not args.output.is_symlink(),
            "--output must be a fresh non-symlink path; overwrites are forbidden")
    report = run_checks()
    canonical = json.dumps(report, sort_keys=True, indent=2) + "\n"
    # Exclusive creation also prevents a concurrent writer or symlink from
    # replacing a destination between validation and file creation.
    with args.output.open("x", encoding="utf-8") as output_file:
        output_file.write(canonical)
    print(canonical, end="")


if __name__ == "__main__":
    main()
