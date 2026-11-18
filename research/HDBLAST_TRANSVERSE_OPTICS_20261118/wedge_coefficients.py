"""Enclose the scalar terms for an inward wedge near a tangent shell ray."""
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
    "SHELL_BOX.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
}
START, STOP = F(49, 32), F(99, 64)
SPAN = STOP - START
ELL = F(1, 10**6)
K = 1 + F(1, 10**11)
C = 1 - F(1, 10**11)
TRANSFER = F(3, 10**6)


def require(check, message):
    if not check:
        raise ValueError(message)


def times(interval, positive_lower, positive_upper):
    lo, hi = interval
    return min(positive_lower*lo, positive_upper*lo), max(positive_lower*hi, positive_upper*hi)


def calculate(directory):
    blobs = {}
    for name, expected in PINS.items():
        raw = (directory / name).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == expected, name)
        blobs[name] = raw
    future = json.loads(blobs["FUTURE_ENCLOSURE.json"])
    previous = json.loads(blobs["PAST_TRANSFER.json"])
    shell = json.loads(blobs["SHELL_BOX.json"])
    rho = tuple(map(F, shell["boxes"]["rho_b"]["exact"]))
    p_out = tuple(map(F, shell["boxes"]["p_b"]["exact"]))
    p_in = -p_out[1], -p_out[0]
    gain = tuple(map(F, future["physical_gain_interval"]))
    grid = future["cell_data"]
    require(grid["sha256"] == PINS["FUTURE_CELLS.csv"], "cell register")
    require(grid["time_grid_denominator"] == 4096 and grid["value_denominator"] == 2**112, "grid")
    require(F(previous["ell"]) == F(future["state"]["ell"]) == ELL, "clock")
    require(F(previous["difference"]["safe_strict_global_upper"]) == F(1, 40), "past transfer")
    require(0 < rho[0] <= rho[1] and 0 < p_in[0] <= p_in[1] and gain[0] < gain[1] < 0,
            "outward sign")
    require(1 / (1 - (ELL*STOP)**2) < K and (1 - (ELL*STOP)**2)**2 > C, "background geometry")
    require(ELL*F(1,40)/(-gain[1]*C) < TRANSFER, "transfer normalization")

    left_edge = START
    moment = [F(0), F(0)]
    clock_error = F(0)
    flux_moment_lower = F(0)
    count = 0
    for row in csv.DictReader(io.StringIO(blobs["FUTURE_CELLS.csv"].decode("ascii"))):
        xlo, xhi = (F(int(row[key]), 4096) for key in ("left_grid", "right_grid"))
        if xhi <= START or xlo >= STOP:
            continue
        require(left_edge == xlo and xhi-xlo == F(1, 4096), "domain tiling")
        left_edge = xhi
        count += 1
        lo, hi = (F(int(row[key]), 2**112) for key in ("H_lo", "H_hi"))
        require(lo <= hi, "ordered input")
        # Unit-interval integrand uses y^2 u g (dt/dy).
        physical = times((lo - TRANSFER, hi + TRANSFER), C**3/K, K)
        left, right = (xlo-START)/SPAN, (xhi-START)/SPAN
        require(right <= F(1,2) or left >= F(1,2), "weight sign on cell")
        mass = (right**3-left**3)/3
        weight = (right**4-left**4)/2 - mass
        signed = sorted(value*weight for value in physical)
        moment[0] += signed[0]
        moment[1] += signed[1]
        clock_error += 2*(K-1)*max(map(abs, physical))*mass
        # Repeat the scalar flux moment in the common gain convention.
        fluence_trace = times((lo-TRANSFER, hi+TRANSFER), C**3, F(1))
        even = (right**4-left**4)/2
        odd = -F(2,3)*(right**6-left**6)
        weights = (even+odd,) if right**2 <= F(1,2) or left**2 >= F(1,2) else (even, odd)
        for coefficient in weights:
            flux_moment_lower += min(coefficient*u for u in fluence_trace)
    require(left_edge == STOP and count == 64, "complete exposure")

    lo, hi = moment[0]-clock_error, moment[1]+clock_error
    require(lo < hi < 0, "negative normalized moment")
    # L/(rho ell SPAN) belongs to [1,K], and p_in*G=-3/rho^3.
    mixed = times((-hi,-lo), ELL**2*SPAN**3/rho[1]**2, K*ELL**2*SPAN**3/rho[0]**2)
    distance = rho[0]*ELL*SPAN, rho[1]*ELL*SPAN*K
    passive = (p_in[0]*distance[0])**2/54, (p_in[1]*distance[1])**2/54
    require(passive[0] > 0 and mixed[0] > 0, "strict positive coefficients")
    require(flux_moment_lower > 0, "positive flux moment")
    norm = F(59,420)
    cross_ratio = F(9,2)*K**3*norm*lo**2/flux_moment_lower**2
    abstract_residual = 1-cross_ratio
    flux_lower = (gain[1]*ELL*SPAN**2)**2*flux_moment_lower**2/(3*K**3*norm)
    residual_lower = abstract_residual*flux_lower
    require(cross_ratio < F(589961,10**6) and residual_lower > F(351758,10**30),
            "quadratic separation")

    return {
        "clock_interval": [str(START), str(STOP)],
        "cells": count,
        "p_in_interval": list(map(str, p_in)),
        "L_interval": list(map(str, distance)),
        "normalized_moment_interval": list(map(str, (lo,hi))),
        "clock_moment_error": str(clock_error),
        "P_interval": list(map(str, passive)),
        "H_interval": list(map(str, mixed)),
        "H_squared_over_4PF_upper": str(cross_ratio),
        "residual_fraction_lower": str(abstract_residual),
        "residual_axis_lower": str(residual_lower),
        "input_sha256": PINS,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", type=Path, nargs="?",
                        default=Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118" / "inputs",
                        help="joint-curvature edition inputs directory")
    print(json.dumps(calculate(parser.parse_args().inputs), indent=2))
