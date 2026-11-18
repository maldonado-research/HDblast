"""Enclose shell curvature slopes and weighted curvature contrasts in model units."""
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path

PINS = {
    "SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "FUTURE_ENCLOSURE.json": "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8",
    "PAST_TRANSFER.json": "0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564",
}


def read(directory, name):
    blob = (directory / name).read_bytes()
    if hashlib.sha256(blob).hexdigest() != PINS[name]:
        raise ValueError(name)
    return blob


def box(low, high=None):
    low = F(low)
    return (low, low if high is None else F(high))


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def mul(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def horner(coefficients, x):
    value = box(0)
    for coefficient in reversed(coefficients):
        value = add(box(coefficient), mul(value, x))
    return value


def square(a):
    low, high = a
    return 0 if low <= 0 <= high else min(low*low, high*high), max(low*low, high*high)


def intersect(a, b):
    low, high = max(a[0], b[0]), min(a[1], b[1])
    if low > high:
        raise ValueError("disjoint enclosures of the same quantity")
    return low, high


def calculation(directory):
    shell = json.loads(read(directory, "SHELL_BOX.json"))
    at_shell = lambda key: box(*shell["boxes"][key]["exact"])
    d, c = (F(shell["parameters_exact"][name]) for name in ("delta", "c"))
    horizon = F(25,16)
    meta = json.loads(read(directory, "FUTURE_ENCLOSURE.json"))
    past = json.loads(read(directory, "PAST_TRANSFER.json"))
    ell = F(meta["state"]["ell"])
    error = F(past["difference"]["safe_strict_global_upper"])
    grid = meta["cell_data"]["time_grid_denominator"]
    denominator = meta["cell_data"]["value_denominator"]
    if (d, c, ell, grid, denominator, error) != (
            F(1,1000), F(96714321425692,161839258930731),
            F(1,10**6), 4096, 2**112, F(1,40)):
        raise ValueError("registered background, clock or data grid")

    # F(eta) = d*(4/3 W*(1+c+c eta)-c*W') + d^2*((1+c+c eta)^2/3-c^2/4).
    b = 1 + c
    polynomial = [
        d*(4*b/9) + d*d*(b*b/3 - c*c/4),
        d*(4*c/9 - 2*c) + d*d*2*b*c/3,
        d*(4*b/3 - c) + d*d*c*c/3,
        d*(4*c/3 + 4*b/9),
        d*4*c/9,
    ]
    derivative = [k*v for k, v in enumerate(polynomial) if k]
    second = [k*v for k, v in enumerate(derivative) if k]
    eta, h2 = at_shell("eta_b"), at_shell("H_squared")
    r_background = intersect(horner(polynomial, eta), mul(box(12), h2))
    slope = intersect(horner(derivative, eta), mul(box(2),mul(at_shell("p_b"),at_shell("B_b"))))
    curvature2 = horner(second, eta)
    if not (slope[1] < 0 and curvature2[0] > 0):
        raise ValueError("sign of effective curvature derivatives")
    ratio = mul(mul(box(4), h2), box(1-(ell*horizon)**2, 1))
    # The explicit angular factor is twice 1/3. The cut weight w=3-4u kills Y^2.
    prefactor = mul(box(F(2,3)), h2)
    rows = list(csv.DictReader(io.StringIO(read(directory, "FUTURE_CELLS.csv").decode("ascii"))))

    def window(left, right):
        linear, projected, variance, length = box(0), box(0), box(0), F(0)
        count, edge = 0, F(0)
        for row in rows:
            a, b = F(int(row["left_grid"]), grid), F(int(row["right_grid"]), grid)
            if a != edge or b <= a:
                raise ValueError("cell coverage")
            edge = b
            lo, hi = max(a,left), min(b,right)
            if hi <= lo:
                continue
            s = (F(int(row["ell_S_lo"]),denominator)/ell-error,
                 F(int(row["ell_S_hi"]),denominator)/ell+error)
            if b <= 1:
                s = intersect(s, box(0))
            if s[0] > s[1]:
                raise ValueError("reversed scalar interval")
            width = hi-lo
            scl = box(width)
            linear = add(linear, mul(scl,mul(slope,s)))
            power = square(s)
            spatial_factor = mul(prefactor,box(1-(ell*hi)**2,1))
            projected = add(projected,mul(scl,mul(spatial_factor,power)))
            variance = add(variance,mul(box(width/6),power))
            length += width
            count += 1
        if length != right-left or edge != 3:
            raise ValueError("window coverage")
        scale = box(1/length)
        return {"interval_x": [left,right], "cells": count,
                "mean_scalar_curvature_quadrupole": mul(linear,scale),
                "mean_quadratic_projected_gauss_residual": mul(projected,scale),
                "mean_leading_scalar_variance": mul(variance,scale)}

    return {"source_sha256": PINS, "F_background": r_background,
            "F_prime": slope, "F_second": curvature2,
            "weighted_residual_to_variance_coefficient": ratio,
            "windows": [window(F(5,4),F(11,8)), window(F(49,32),F(25,16)),
                        window(F(0),horizon)]}


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent / "inputs"
    print(json.dumps(calculation(directory), indent=2, default=str))
