"""Two rational lower bounds on the null scalar flux from the registered shell trace."""

import argparse
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path


PINS = {
    "FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "FUTURE_ENCLOSURE.json": "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8",
    "PAST_TRANSFER.json": "0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564",
}
A, B = F(49, 32), F(99, 64)
SPAN = B - A
ELL = F(1, 10**6)
GRID = 4096
K = 1 + F(1, 10**11)
C = 1 - F(1, 10**11)
TRANSFER = F(3, 10**6)


def require(ok, text):
    if not ok:
        raise ValueError(text)


def read_inputs(directory):
    verified = {}
    for name, digest in PINS.items():
        raw = (directory / name).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == digest, name)
        verified[name] = raw
    future = json.loads(verified["FUTURE_ENCLOSURE.json"])
    past = json.loads(verified["PAST_TRANSFER.json"])
    grid = future["cell_data"]
    require(grid["sha256"] == PINS["FUTURE_CELLS.csv"], "cell association")
    den = grid["value_denominator"]
    require(grid["time_grid_denominator"] == GRID and den == 2**112, "source grid")
    require(F(past["ell"]) == F(future["state"]["ell"]) == ELL, "scaled clock")
    gain_left, gain_right = map(F, future["physical_gain_interval"])
    error = F(past["difference"]["safe_strict_global_upper"])
    require(gain_left < gain_right < 0 and error == F(1, 40), "gain and past transfer")
    require(1 / (1 - (ELL*B)**2) < K and (1 - (ELL*B)**2)**2 > C, "cosh and dilution floors")
    require(ELL*error / (-gain_right*C) < TRANSFER, "transfer to conjugated units")
    cells = []
    for row in csv.DictReader(io.StringIO(verified["FUTURE_CELLS.csv"].decode("ascii"))):
        x0, x1 = (F(int(row[key]), GRID) for key in ("left_grid", "right_grid"))
        if x1 <= A or x0 >= B:
            continue
        require(x0 == (cells[-1][1] if cells else A) and x1 == x0+F(1, GRID), "cell coverage")
        low, high = (F(int(row[key]), den) for key in ("H_lo", "H_hi"))
        require(low <= high, "ordered cell")
        cells.append((x0, x1, low - TRANSFER, high + TRANSFER))
    require(len(cells) == 64 and cells[-1][1] == B, "complete optical interval")
    return cells, -gain_right


def apply_geometry(low, high):
    # The product of the three positive geometric factors lies in [C^3, 1].
    return min(low, C**3 * low), max(high, C**3 * high)


def moment_bound(cells, gain):
    lo, hi = F(0), F(0)
    for x0, x1, low, high in cells:
        left, right = (x0 - A)/SPAN, (x1 - A)/SPAN
        low, high = apply_geometry(low, high)
        # Primitive for 2y^3-4y^5, integrated with its sign on each cell.
        p = (right**4 - left**4)/2
        n = -F(2, 3)*(right**6 - left**6)
        if right**2 <= F(1, 2) or left**2 >= F(1, 2):
            weights = (p + n,)
        else:
            weights = (p, n)
        for w in weights:
            lo += min(w*low, w*high)
            hi += max(w*low, w*high)
    require(lo > 0, "sign of weighted moment")
    norm = F(1, 4) + F(1, 5) - F(1, 6) - F(1, 7)
    require(norm == F(59, 420), "polynomial norm")
    scale = (gain*ELL*SPAN**2)**2/(3*K**3)
    return (lo, hi), scale * lo**2 / norm


def point_bound(cells, gain):
    points = []
    for first, second in zip(cells, cells[1:]):
        x = first[1]
        low = max(first[2], second[2])
        high = min(first[3], second[3])
        require(low <= high, "closed cell intersection")
        low, high = apply_geometry(low, high)
        points.append((x, low*(x-A)**2, high*(x-A)**2))
    total = F(0)
    for (left, lo1, hi1), (right, lo2, hi2) in zip(points, points[1:]):
        gap = max(F(0), lo2 - hi1, lo1 - hi2)
        w0, w1 = (left-A)*(B-left), (right-A)*(B-right)
        # The reciprocal of (x-A)(B-x) is convex on this interior segment.
        reciprocal_integral_upper = (right-left)*(1/w0 + 1/w1)/2
        total += gap**2/reciprocal_integral_upper
    return total * (gain*ELL)**2 / (3*SPAN*K**3)


def compute(directory):
    cells, gain = read_inputs(directory)
    moment, main = moment_bound(cells, gain)
    point = point_bound(cells, gain)
    cut = F(4, 5)
    cap = (F(19, 40)*(1-cut**5) - F(7, 12)*(1-cut**3) + F(3, 8)*(1-cut))/(1-cut)
    full = F(2, 5) - F(2, 15)
    require(cap == F(10283, 18750) and full == F(4, 15), "sky moments")
    require(main > F(857864, 10**30) and point > F(84783, 10**29), "axis checks")
    return {
        "clock": [str(A), str(B)],
        "cells": len(cells),
        "gain_abs_lower": str(gain),
        "moment_interval": [str(v) for v in moment],
        "axis_polynomial_lower": str(main),
        "axis_interior_pairs_lower": str(point),
        "interior_pairs": len(cells) - 2,
        "cap_mean_quadrupole_squared": str(cap),
        "cap_lower": str(main * cap),
        "sky_mean_quadrupole_squared": str(full),
        "sky_lower": str(main * full),
        "input_sha256": PINS,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="?",
                        default=Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118" / "inputs",
                        help="joint-curvature edition inputs directory")
    print(json.dumps(compute(parser.parse_args().inputs), indent=2))
