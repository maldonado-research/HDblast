"""Enclose the linear area distances of two finite cones at fixed affine length."""

import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

PARENT = Path(__file__).resolve().parent.parent / "HDBLAST_JOINT_CURVATURE_20261118"
MANIFEST_SHA = "f44f963be6065b11753e69dcb3f67d443d1213b0859db68ad9d5e72539dd10dd"
ELL = F(1, 10**6)
START, STOP = F(49, 32), F(99, 64)
SPAN = STOP - START
TRANSFER_ERROR = F(3, 10**6)
COSH_MAX = 1 + F(1, 10**11)


def verify(condition, message):
    if not condition:
        raise ValueError(message)


def verify_inputs():
    manifest = (PARENT / "SHA256SUMS").read_bytes()
    verify(hashlib.sha256(manifest).hexdigest() == MANIFEST_SHA, "antecedent manifest")
    pins = {}
    for line in manifest.decode().splitlines():
        checksum, name = line.split("  ", 1)
        pins[name] = checksum
    for name in ("inputs/FUTURE_CELLS.csv", "inputs/SHELL_BOX.json", "METRIC.json"):
        verify(hashlib.sha256((PARENT / name).read_bytes()).hexdigest() == pins[name], name)


def weighted_cell(start, dt, position, velocity, acceleration):
    """Integrate (x-START)^2 times a quadratic trajectory, divided by SPAN^3."""
    offset = start - START
    weight_coefficients = (offset**2, 2 * offset, F(1))
    trajectory_coefficients = (position, velocity, acceleration / 2)
    total = F(0)
    for i, weight in enumerate(weight_coefficients):
        for j, coefficient in enumerate(trajectory_coefficients):
            power = i + j + 1
            total += weight * coefficient * dt**power / power
    return total / SPAN**3


def integrate_cells():
    position, velocity, integrated = [[F(0), F(0)] for _ in range(3)]
    end = F(1)
    with (PARENT / "inputs/FUTURE_CELLS.csv").open() as stream:
        for row in csv.DictReader(stream):
            left = max(F(1), F(int(row["left_grid"]), 4096))
            right = min(STOP, F(int(row["right_grid"]), 4096))
            if right <= left:
                continue
            verify(left == end, "grid coverage")
            dt = right - left
            acceleration = (F(int(row["H_lo"]), 2**112) - TRANSFER_ERROR,
                            F(int(row["H_hi"]), 2**112) + TRANSFER_ERROR)
            verify(acceleration[0] <= acceleration[1], "ordered cell")
            for side, acc in enumerate(acceleration):
                if right > START:
                    verify(left >= START, "exposure aligned to grid")
                    integrated[side] += weighted_cell(left, dt, position[side], velocity[side], acc)
                position[side] += dt * velocity[side] + dt**2 * acc / 2
                velocity[side] += dt * acc
            end = right
    verify(end == STOP, "endpoint")
    return position, integrated


def compute():
    verify_inputs()
    metric = json.loads((PARENT / "METRIC.json").read_text())
    low, high = map(F, metric["uniform_null_abs_intervals"][1])
    rough = (low - 2 * high / 3, high - 2 * low / 3)
    verify(rough[0] > F("0.05161") and rough[1] < F("0.10710"), "uniform enclosure")

    position, integrated = integrate_cells()
    z_bound = F(metric["z_bound"])
    # The parent's uniform bound for M-z0 plus the loss from sech(q)^2.
    m_error = F(metric["absolute_mean_error"]) + (ELL * F(25, 16))**2 * z_bound
    verify(1 / (1 - (ELL * STOP)**2) < COSH_MAX, "affine-clock Jacobian")
    error = F(5, 3) * m_error + F(2, 3) * (COSH_MAX**3 - 1) * z_bound
    axis = (position[0] - 2 * integrated[1] - error,
            position[1] - 2 * integrated[0] + error)
    verify(rough[0] < axis[0] <= axis[1] < rough[1], "refined coefficient")

    angular = F(18, 25)
    cap = tuple(angular * coefficient for coefficient in axis)
    # L/rho = integral cosh(q) dq / cosh(q0). The integrand ranges from 1 to COSH_MAX.
    length_low, length_high = ELL * SPAN, COSH_MAX * ELL * SPAN
    box = json.loads((PARENT / "inputs/SHELL_BOX.json").read_text())["boxes"]
    b_low, b_high = map(F, box["B_b"]["exact"])
    rho_low, rho_high = map(F, box["rho_b"]["exact"])
    relative = (cap[0] * ELL * b_low / rho_high * length_low**2,
                cap[1] * ELL * b_high / rho_low * length_high**2)
    return {
        "parent_manifest_sha256": MANIFEST_SHA,
        "vertex_x": str(START),
        "endpoint_background_x": str(STOP),
        "rough_axis_coefficient": list(map(str, rough)),
        "refined_axis_coefficient": list(map(str, axis)),
        "solid_angle_fraction": "1/10",
        "cap_angular_factor": str(angular),
        "refined_cap_coefficient": list(map(str, cap)),
        "fractional_distance_first_variation_positive_cap": list(map(str, relative)),
        "display_only": {
            "axis": list(map(float, axis)),
            "cap": list(map(float, cap)),
            "fractional_distance": list(map(float, relative)),
        },
    }


if __name__ == "__main__":
    print(json.dumps(compute(), indent=2))
