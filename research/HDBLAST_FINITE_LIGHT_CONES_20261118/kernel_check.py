"""Independent positive-kernel check of the registered optical axis coefficient."""
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PARENT = ROOT.parent / "HDBLAST_JOINT_CURVATURE_20261118"
CSV_SHA = "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7"
METRIC_SHA = "d7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42"
BEAM_SHA = "40fd88fdf92806c78b9d84dbfc7f256e9a9ec396dba24e5b5b76e87021db0652"
START, STOP = F(49, 32), F(99, 64)
SPAN = STOP - START


def authenticated(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"input digest mismatch: {path.name}")
    return data


def primitive(t):
    # Exchanging the optical-weight and acceleration integrals gives a positive
    # kernel K(t)=SPAN/2 + (START-t)/3 for t <= START.
    # For START <= t <= STOP, y=(t-START)/SPAN and K=SPAN*(1/2-y/3-y**4/6).
    # These primitives agree at START. The second kernel decreases to zero.
    offset = t - START
    if t <= START:
        return SPAN * offset / 2 - offset**2 / 6
    y = offset / SPAN
    return SPAN**2 * (y / 2 - y**2 / 6 - y**5 / 30)


def compute():
    csv_data = authenticated(PARENT / "inputs/FUTURE_CELLS.csv", CSV_SHA)
    metric = json.loads(authenticated(PARENT / "METRIC.json", METRIC_SHA))
    beam = json.loads(authenticated(ROOT / "BEAM.json", BEAM_SHA))
    if (F(beam["vertex_x"]), F(beam["endpoint_background_x"])) != (START, STOP):
        raise ValueError("clock registration")

    low = high = F(0)
    edge = F(1)
    count = 0
    radius = F(3, 10**6)
    for row in csv.DictReader(io.StringIO(csv_data.decode("ascii"))):
        left = max(F(1), F(int(row["left_grid"]), 4096))
        right = min(STOP, F(int(row["right_grid"]), 4096))
        if left >= right:
            continue
        mass = primitive(right) - primitive(left)
        if left != edge or mass <= 0:
            raise ValueError("continuous cell partition or kernel mass")
        lower = F(int(row["H_lo"]), 2**112) - radius
        upper = F(int(row["H_hi"]), 2**112) + radius
        if upper < lower:
            raise ValueError("acceleration enclosure")
        low += mass * lower
        high += mass * upper
        count += 1
        edge = right
    if edge != STOP:
        raise ValueError("incomplete interval")

    ell, horizon = F(1, 10**6), F(25, 16)
    z_bound = F(metric["z_bound"])
    m_error = F(metric["absolute_mean_error"]) + (ell * horizon)**2 * z_bound
    # The exact normalized affine weight integrates to 1/3. Its ratio to the
    # uniform-clock quadratic weight lies in [cosh_max**-3, cosh_max**3].
    cosh_max = 1 + F(1, 10**11)
    if 1 / (1 - (ell * STOP)**2) >= cosh_max:
        raise ValueError("affine Jacobian range")
    error = (5 * m_error + 2 * (cosh_max**3 - 1) * z_bound) / 3
    axis = (low - error, high + error)
    outer_low, outer_high = map(F, beam["refined_axis_coefficient"])
    if not outer_low < axis[0] <= axis[1] < outer_high:
        raise ValueError("independent result outside the optical enclosure")
    return {
        "csv_sha256": CSV_SHA,
        "metric_sha256": METRIC_SHA,
        "beam_sha256": BEAM_SHA,
        "cells": count,
        "flat_axis": [str(low), str(high)],
        "physical_time_and_metric_error": str(error),
        "axis_coefficient": list(map(str, axis)),
        "axis_display_only": list(map(float, axis)),
    }


if __name__ == "__main__":
    print(json.dumps(compute(), indent=2))
