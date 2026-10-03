#!/usr/bin/env python3
"""Exact arithmetic checks of analytic bound constants; no source evaluation."""

from fractions import Fraction as F
from math import factorial


def check_constants():
    rho = F(5, 8)
    radius = F(1, 8)
    half_width = F(1, 128)
    denominator_bound = 1 - rho**2
    lapse_bound = F(8, 27)
    bump_bound = 1 + rho**2
    log_derivative_bound = 2 * rho / denominator_bound**2
    second_derivative_bound = (
        log_derivative_bound**2
        + 2 / denominator_bound**2
        + 8 * rho**2 / denominator_bound**3
    )
    positive_bound = bump_bound * (
        4 * lapse_bound**2
        + 2 * lapse_bound * log_derivative_bound
        + second_derivative_bound
    )
    signed_bound = bump_bound * (
        4 * lapse_bound**2 * rho
        + 2 * lapse_bound * (1 + rho * log_derivative_bound)
        + 2 * log_derivative_bound
        + rho * second_derivative_bound
    )
    assert denominator_bound == F(39, 64)
    assert bump_bound == F(89, 64)
    assert log_derivative_bound == F(5120, 1521)
    assert second_derivative_bound == F(70623232, 2313441)
    assert positive_bound == F(951819044, 20820969) < 64
    assert signed_bound == F(3227905133, 83283876) < 64
    assert lapse_bound * 64 == F(512, 27) < 32
    e_upper = F(5, 2) + F(1, 6) / (1 - F(1, 4))
    assert e_upper == F(49, 18) < F(11, 4)
    assert F(11, 4) ** 8 == F(214358881, 65536) < 4096

    centers = [F(-9, 2) + F(2 * j + 1, 128) for j in range(64)]
    assert centers[0] - half_width == F(-9, 2)
    assert centers[-1] + half_width == F(-7, 2)
    assert all(b - a == 2 * half_width for a, b in zip(centers, centers[1:]))
    assert all(abs(c + 4) + radius <= rho for c in centers)
    assert all(abs(c) - radius >= F(27, 8) for c in centers)
    assert 64 * 2 * half_width == 1

    q = half_width / radius
    assert q == F(1, 16)
    for degree in (20, 24, 28, 32):
        remainder = 64 * q ** (degree + 1) / (1 - q)
        exact = F(1, 15 * 2 ** (4 * degree - 6))
        assert remainder == exact
        print(f"degree={degree} uniform_tail={remainder} L1_tail={remainder / (degree + 2)}")

    for cutoff, numerator, exponent in ((48, 1133, 36), (56, 1101, 45), (64, 3514, 55)):
        tail = F(2 * 4 ** (cutoff + 1), factorial(cutoff + 1)) / (1 - F(4, cutoff + 2))
        assert tail < F(numerator, 10**exponent)
        print(f"moment_cutoff={cutoff} tail_upper={numerator}/10^{exponent}")
    print("EXACT_ANALYTIC_BOUND_CONSTANT_CHECKS_PASS")


if __name__ == "__main__":
    check_constants()
