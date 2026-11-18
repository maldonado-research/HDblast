"""Interval certificate for two finite exposures of reciprocal displacement releases."""
import csv
import hashlib
import io
import json
import sys
from fractions import Fraction as F
from pathlib import Path

PINS = {
    "FUTURE_CELLS.csv": "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7",
    "FUTURE_ENCLOSURE.json": "d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8",
    "PAST_TRANSFER.json": "0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564",
}
ELL = F(1, 10**6)
END = F(25, 16)
H = F(1, 32)
LEFT = F(5, 4)
WIDTH = F(1, 8)
FLOOR = F(1) - F(1, 10**11)  # For sech(ell*x)^(3/2) on the observed horizon.
RADIUS = F(3, 10**6)
FACTOR = F(1053, 1000)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def shutter_primitive(x):
    """Antiderivative of normalized triangular shutter on [LEFT,LEFT+WIDTH]."""
    t = x - LEFT
    return (2 / WIDTH**2) * (
        max(t, 0)**2 - 2*max(t - WIDTH/2, 0)**2 + max(t - WIDTH, 0)**2
    )


def calculate(directory):
    def read(name):
        content = (directory / name).read_bytes()
        require(hashlib.sha256(content).hexdigest() == PINS[name], name)
        return content

    meta = json.loads(read("FUTURE_ENCLOSURE.json"))
    transfer = json.loads(read("PAST_TRANSFER.json"))
    grid = meta["cell_data"]
    require(grid["sha256"] == PINS["FUTURE_CELLS.csv"] and
            grid["time_grid_denominator"] == 4096 and grid["value_denominator"] == 2**112,
            "continuous cell grid")
    require(F(meta["state"]["ell"]) == F(transfer["ell"]) == ELL, "clock scale")
    physical_error = F(transfer["difference"]["safe_strict_global_upper"])
    gain_left, gain_right = map(F, meta["physical_gain_interval"])
    require(physical_error == F(1, 40) and gain_left < gain_right < 0, "past transfer")
    require(1 - 3*(ELL*END)**2/4 > FLOOR and
            ELL*physical_error / (-gain_right*FLOOR) < RADIUS, "physical corrections")

    # The release may occur at any s in [0, 1-H]. Numbers bound this entire range.
    gamma = 1/FLOOR                          # y=cosh^(3/2)(ell*x)*a / cosh^(3/2)(ell*s)
    residual = F(35, 4) * ELL**4 * END**2   # Omega^2 minus conjugated potential
    slope = F(3, 2) * ELL**2 * (1-H)        # y'(s) for a(s)=1, a'(s)=0
    early_error = slope*H + gamma*residual*H**2/2
    high_error = 3*H*slope + 6*early_error + H**2*gamma*residual
    low_lower = F(19, 20) - early_error
    high_lower = F(19, 20) - high_error
    low_upper = gamma
    high_upper = gamma*(16*H**2/WIDTH**2 + H**2*residual)
    require(low_lower > 0 and high_lower > 0, "positive early response")
    require(max(low_upper/low_lower, high_upper/high_lower) < FACTOR, "all masses")

    total = F(0)
    numerator_low = numerator_high = box_low = F(0)
    count = 0
    edge = F(0)
    rows = csv.DictReader(io.StringIO(read("FUTURE_CELLS.csv").decode("ascii")))
    for row in rows:
        lft = F(int(row["left_grid"]), 4096)
        rgt = F(int(row["right_grid"]), 4096)
        require(lft == edge < rgt, "continuous coverage")
        edge = rgt
        start, end = max(lft, LEFT), min(rgt, LEFT+WIDTH)
        if start >= end:
            continue
        lo = F(int(row["H_lo"]), 2**112) - RADIUS
        hi = F(int(row["H_hi"]), 2**112) + RADIUS
        require(0 < lo <= hi, "positive late signal on each cell")
        weight = shutter_primitive(end) - shutter_primitive(start)
        require(weight > 0, "positive shutter")
        numerator_low += weight*lo
        numerator_high += weight*hi
        box_low += (end-start)*lo*FLOOR/WIDTH
        total += weight
        count += 1
    require(edge == 3 and total == 1 and count == 512, "complete late exposure")
    # Numerator is the conjugated trace. A normalized physical shutter divides
    # by the mean of cosh^(3/2), between 1 and 1/FLOOR.
    late = [numerator_low*FLOOR, numerator_high]
    discrepancy = late[0]/(1+FACTOR)

    return {
        "input_sha256": PINS,
        "early_duration": H, "release_range": [F(0), 1-H], "late_duration": WIDTH,
        "late_start": LEFT, "constant": FACTOR,
        "conjugated_potential_shift_upper": residual,
        "early_flat_error_upper": early_error, "early_scaled_flat_error_upper": high_error,
        "ratio_upper_by_frequency_case": [low_upper/low_lower, high_upper/high_lower],
        "late_target_normalized_exposure": late, "late_cells": count,
        "point_and_box_equal_error_floor": box_low/2,
        "two_exposure_equal_error_floor": discrepancy,
        "two_exposure_relative_error_floor": 1/(1+FACTOR),
    }


if __name__ == "__main__":
    if len(sys.argv) > 2:
        raise SystemExit("usage: reciprocal.py [capsule-input-directory]")
    source = Path(sys.argv[1]) if len(sys.argv) == 2 else Path(__file__).resolve().parent / "inputs"
    print(json.dumps(calculate(source), indent=2, default=str))
