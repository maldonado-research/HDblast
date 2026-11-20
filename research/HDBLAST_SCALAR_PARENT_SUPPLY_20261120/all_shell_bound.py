"""Bound the neutral conformal scalar population on the archived precursor wall.

For j=n+1, the oscillator frequency in conformal time is
omega_j = sqrt(j**2 + c**2), where c = a*M, M=2*R(psi) and mL=1.
An instantaneous-vacuum initial state obeys
    asinh(abs(beta_j)) <= total_variation(log(omega_j))/2.
A real scalar on the round three-sphere has j**2 independent modes.

The primary HDblast mass law is constant at M=1.9 for psi<=0.755,
equals 2 at psi=0.77, and rises monotonically to M=2.1 at psi>=0.785.
The checks below prove that c decreases to one minimum and then increases
on the archived Hermite interpolation. Endpoints therefore bound the total
variation for every j without differentiating a high-frequency solution.

Use a positional source path to override the adjacent HDblast archive.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from math import ceil, isqrt
from pathlib import Path


ARCHIVE_SHA256 = "5dbed28e158584403a4cd238d392aab15c00480522cf489d1ee0a36f1d550f90"
DEFAULT_SOURCE = (Path(__file__).resolve().parent.parent
                  / "HDBLAST_NEUTRAL_PARENT_DECAY_20261120" / "prior" / "CHARGED.json")
GATE_START = Fraction("0.755")
GATE_MIDPOINT = Fraction("0.77")
INITIAL_MASS = Fraction("1.9")
INITIAL_C = Fraction("3.8")
MINIMUM_C_LOWER = Fraction("3.5")
LAST_SUMMED_J = 1000
BAND_LAST_J = 129


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def hermite_controls(left, right, key, left_slope, right_slope):
    """Return cubic value controls and quadratic proper-time derivative controls."""
    duration = right["t"] - left["t"]
    require(duration > 0, "Unordered precursor trace")
    start, end = left[key], right[key]
    controls = [start, start + duration*left_slope/3,
                end - duration*right_slope/3, end]
    slopes = [left_slope, 3*(end-start)/duration-left_slope-right_slope, right_slope]
    return controls, slopes


def prefix_controls(controls, fraction):
    """de Casteljau coefficients of the prefix [0,fraction] of one segment."""
    require(0 <= fraction <= 1, "Invalid prefix fraction")
    row = controls
    prefix = [row[0]]
    while len(row) > 1:
        row = [(1-fraction)*left + fraction*right for left, right in zip(row, row[1:])]
        prefix.append(row[0])
    return prefix


def quadratic_halves(controls):
    left, mid, right = controls
    center = (left + 2*mid + right)/4
    return (left, (left+mid)/2, center), (center, (mid+right)/2, right)


def geometry_certificate(rows, turnaround, at_turnaround):
    """Certify positive radius, a single minimum in c, and its endpoint bounds.

    A polynomial lies inside the convex hull of its Bernstein coefficients.
    A constant sign for all derivative coefficients, on either half when
    necessary, proves monotonicity. A strictly increasing sequence of
    derivative coefficients with opposite endpoint signs proves a unique
    minimum on that segment. The mass gate is constant until after it.
    """
    require(rows[0]["t"] == 0 and rows[0]["a"]*INITIAL_MASS == INITIAL_C
            and rows[0]["psi"] < GATE_START, "Unexpected initial conditions")
    radius_upper = Fraction(10) if at_turnaround else Fraction("9.6")
    increasing = False
    bounce_index = None
    for index, (left, right) in enumerate(zip(rows, rows[1:])):
        radius, radius_slope = hermite_controls(
            left, right, "a", left["a"]*left["H"], right["a"]*right["H"])
        psi, psi_slope = hermite_controls(left, right, "psi", left["pi"], right["pi"])
        contains_turnaround = left["t"] < turnaround <= right["t"]
        if at_turnaround and contains_turnaround:
            fraction = (turnaround-left["t"])/(right["t"]-left["t"])
            radius, radius_slope, psi, psi_slope = (
                prefix_controls(values, fraction) for values in
                (radius, radius_slope, psi, psi_slope))
        require(INITIAL_MASS*min(radius) > MINIMUM_C_LOWER and max(radius) < radius_upper,
                "Radius left its certified interval")
        require(max(psi) <= GATE_START or min(psi_slope) >= 0,
                "Cannot prove the gate is constant or increasing")
        if increasing:
            require(all(min(part) >= 0 for part in quadratic_halves(radius_slope)),
                    "Unexpected oscillation after bounce")
        else:
            require(max(psi) < GATE_START, "Mass gate active before the bounce")
            if radius_slope[0] < 0 < radius_slope[-1]:
                require(radius_slope[0] < radius_slope[1] < radius_slope[2],
                        "Unresolved bounce")
                bounce_index = index
                increasing = True
            else:
                require(all(max(part) <= 0 for part in quadratic_halves(radius_slope)),
                        "Unexpected oscillation during contraction")
        contains_opening = left["psi"] < GATE_MIDPOINT <= right["psi"]
        finished = contains_turnaround if at_turnaround else contains_opening
        if finished:
            require(increasing, "No bounce before requested endpoint")
            # At the first crossing the continuous gate gives exactly M=2.
            mass_upper = Fraction("2.1") if at_turnaround else Fraction(2)
            return radius_upper*mass_upper, bounce_index, index
    raise ArithmeticError("Requested endpoint absent")


def ceiling_decimal(value, digits):
    """Format with outward rounding using integers only."""
    scale = 10**digits
    whole, fraction = divmod(ceil(value*scale), scale)
    return f"{whole}.{fraction:0{digits}d}"


def population_certificate(final_c_upper):
    """Enclose the finite and infinite parts of the state-count upper bound."""
    finite = Fraction(0)
    finite_above_band = Fraction(0)
    sqrt_scale = 10**10
    for j in range(1, LAST_SUMMED_J+1):
        # b_j=¼ log U bounds the squeezing through this endpoint. If v=sqrt(U),
        # sinh²(b_j)=(v-1)²/(4v), increasing for v>=1.
        ratio = ((j*j + INITIAL_C**2)*(j*j + final_c_upper**2)
                 / (j*j + MINIMUM_C_LOWER**2)**2)
        require(ratio >= 1, "Invalid logarithmic envelope")
        root_floor = isqrt(ratio.numerator*sqrt_scale**2 // ratio.denominator)
        root_upper = Fraction(root_floor+1, sqrt_scale)
        shell = j*j*(root_upper-1)**2/(4*root_upper)
        shell_upper = Fraction(ceil(shell*10**9), 10**9)
        finite += shell_upper
        if j > BAND_LAST_J:
            finite_above_band += shell_upper

    # log(1+x)<=x gives b_j<=D/j². For 0<=x<1, sinh(x)<=x/(1-x).
    # The integral test gives sum(j^-2 for j>J)<=1/J.
    variation = (INITIAL_C**2+final_c_upper**2-2*MINIMUM_C_LOWER**2)/4
    first_tail_bound = variation/(LAST_SUMMED_J+1)**2
    require(0 < first_tail_bound < 1, "Invalid tail envelope")
    tail = (variation/(1-first_tail_bound))**2/LAST_SUMMED_J
    return {
        "finite_sum_upper_through_j1000": ceiling_decimal(finite, 8),
        "tail_upper_above_j1000": ceiling_decimal(tail, 8),
        "unresolved_upper_above_j129": ceiling_decimal(finite_above_band+tail, 8),
        "all_modes_upper": ceiling_decimal(finite+tail, 8),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    blob = args.source.read_bytes()
    require(hashlib.sha256(blob).hexdigest() == ARCHIVE_SHA256, "Unexpected source SHA256")
    data = json.loads(blob, parse_float=Fraction, parse_int=Fraction)
    run = data["runs"]["neutral_precursor"]
    rows = [dict(zip(data["trace_columns"], row, strict=True)) for row in run["trace"]]
    results = {"source_sha256": ARCHIVE_SHA256, "trace_rows": len(rows)}
    for name, at_turnaround in (("opening", False), ("turnaround", True)):
        final_c, bounce, end = geometry_certificate(rows, run["maximum_tau"], at_turnaround)
        certificate = population_certificate(final_c)
        require(Fraction(certificate["all_modes_upper"]) < (699 if at_turnaround else 525),
                "Population ceiling exceeded")
        results[name] = {"bounce_segment_index": bounce, "last_segment_index": end,
                         "final_a_times_mass_upper": str(final_c), **certificate}
    output = json.dumps(results, indent=2) + "\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
