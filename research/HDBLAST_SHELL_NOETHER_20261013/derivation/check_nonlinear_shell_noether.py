#!/usr/bin/env python3
"""Exact algebra checks for the nonlinear shell Noether calculation.

This checks finite-dimensional tensor/algebra identities, not nonlinear PDE
existence.  SymPy is the only nonstandard dependency.  No assertion is used as
a correctness guard, so running Python with -O does not remove a check.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import sympy as sp


HALF = sp.Rational(1, 2)


def symmetric_matrix(prefix: str, size: int) -> sp.Matrix:
    entries = {}
    for a in range(size):
        for b in range(a, size):
            entries[a, b] = sp.Symbol(f"{prefix}{a}{b}")
    return sp.Matrix(size, size, lambda a, b: entries[min(a, b), max(a, b)])


def reduced(expression):
    if isinstance(expression, sp.MatrixBase):
        return expression.applyfunc(sp.expand)
    return sp.expand(expression)


def is_zero(expression) -> bool:
    result = reduced(expression)
    if isinstance(result, sp.MatrixBase):
        return all(value == 0 for value in result)
    return result == 0


def check_zero(name: str, expression, results: dict) -> None:
    if not is_zero(expression):
        raise RuntimeError(f"{name}: exact identity failed: {reduced(expression)}")
    components = len(expression) if isinstance(expression, sp.MatrixBase) else 1
    results[name] = {"status": "pass", "exact_residual": "0", "components": components}


def check_detected(name: str, expression, results: dict) -> None:
    result = reduced(expression)
    if is_zero(result):
        raise RuntimeError(f"{name}: the deliberately incorrect control was not detected")
    if isinstance(result, sp.MatrixBase):
        first = next(value for value in result if value != 0)
    else:
        first = result
    results[name] = {"status": "detected", "one_nonzero_exact_residual": str(first)}


def gaussian_ghy_tensor_check(results: dict, controls: dict, inject_wrong_corner_sign: bool = False) -> None:
    """Build theta and delta K from covariant derivatives at a Gaussian point.

    Background first jets have d_r h_ab=2K_ab and vanishing tangential metric
    jets.  All metric-variation components and their first jets are arbitrary.
    This does not freeze the varied lapse, shift, normal, or induced metric.
    """
    size = 5
    inverse = sp.diag(1, -1, 1, 1, 1)
    induced_inverse = inverse[1:, 1:]
    extrinsic = symmetric_matrix("k", 4)
    variation = symmetric_matrix("q", 5)
    jets = [symmetric_matrix(f"d{a}q", size) for a in range(size)]
    connection = [[[sp.S.Zero for _ in range(size)] for _ in range(size)] for _ in range(size)]
    for a in range(1, size):
        for b in range(1, size):
            connection[0][a][b] = -extrinsic[a - 1, b - 1]
            connection[a][0][b] = inverse[a, a] * extrinsic[a - 1, b - 1]
            connection[a][b][0] = connection[a][0][b]

    def nabla_q(a: int, b: int, c: int):
        return jets[a][b, c] - sum(
            connection[d][a][b] * variation[d, c]
            + connection[d][a][c] * variation[b, d]
            for d in range(size)
        )

    def delta_connection(a: int, b: int, c: int):
        return HALF * inverse[a, a] * (
            nabla_q(b, c, a) + nabla_q(c, b, a) - nabla_q(a, b, c)
        )

    theta_normal = HALF * sum(
        inverse[b, b] * (nabla_q(b, 0, b) - nabla_q(0, b, b))
        for b in range(size)
    )
    # n_A is the normalized conormal: delta n_r=q_rr/2, delta n_a=0.
    delta_extrinsic = sp.Matrix(
        4, 4,
        lambda a, b: HALF * extrinsic[a, b] * variation[0, 0]
        - delta_connection(0, a + 1, b + 1),
    )
    q_induced = variation[1:, 1:]
    raised_q_induced = induced_inverse * q_induced * induced_inverse
    trace_k = sp.trace(induced_inverse * extrinsic)
    trace_q = sp.trace(induced_inverse * q_induced)
    contraction_kq = sum(extrinsic[a, b] * raised_q_induced[a, b] for a in range(4) for b in range(4))
    delta_trace_k = sp.trace(induced_inverse * delta_extrinsic) - contraction_kq
    delta_ghy_density = delta_trace_k + HALF * trace_k * trace_q
    canonical_density = -HALF * contraction_kq + HALF * trace_k * trace_q
    divergence_shift = sum(inverse[a, a] * jets[a][0, a] for a in range(1, size))
    corner_divergence = (HALF if inject_wrong_corner_sign else -HALF) * divergence_shift
    check_zero(
        "unrestricted_ghy_from_connection_and_normal_variations",
        theta_normal + delta_ghy_density - canonical_density - corner_divergence,
        results,
    )
    check_detected(
        "wrong_sign_ghy_corner",
        theta_normal + delta_ghy_density - canonical_density + corner_divergence,
        controls,
    )
    frozen_normal_density = delta_ghy_density - HALF * trace_k * variation[0, 0]
    check_detected(
        "omit_unit_normal_variation",
        theta_normal + frozen_normal_density - canonical_density - corner_divergence,
        controls,
    )


def junction_and_charge_checks(results: dict, controls: dict) -> None:
    h = sp.diag(-1, 1, 1, 1)
    kp = symmetric_matrix("kp", 4)
    km = symmetric_matrix("km", 4)
    sigma, sigma_prime, pp, pm = sp.symbols("sigma sigma_prime pp pm")
    trp = sp.trace(h * kp)
    trm = sp.trace(h * km)
    israel = kp + km - (trp + trm) * h + sigma * h
    scalar = pp + pm + sigma_prime
    tri = sp.trace(h * israel)
    check_zero("four_dimensional_israel_trace", trp + trm - (4 * sigma - tri) / 3, results)
    check_zero("four_dimensional_israel_trace_reversal", kp + km - israel - (sigma - tri) * h / 3, results)

    def squared_norm(matrix):
        return sp.trace(h * matrix * h * matrix)

    # Common intrinsic R, tangential scalar gradient, and U cancel exactly.
    twice_delta_e_nn = trp**2 - trm**2 - squared_norm(kp) + squared_norm(km) - pp**2 + pm**2
    delta_k = kp - km
    contraction_i_delta_k = sp.trace(h * israel * h * delta_k)
    force = sigma * (trp - trm) + sigma_prime * (pp - pm)
    source = twice_delta_e_nn + contraction_i_delta_k + scalar * (pp - pm)
    check_zero("full_nonlinear_normal_force_factorization", force - source, results)
    check_detected("wrong_scalar_junction_sign_in_force", force - (
        twice_delta_e_nn + contraction_i_delta_k
        + (pp + pm - sigma_prime) * (pp - pm)
    ), controls)
    check_detected("omit_nonlinear_israel_force_product", force - (twice_delta_e_nn + scalar * (pp - pm)), controls)

    # A = tangential projection of n·nabla xi is independent on each side.
    xi = sp.Matrix(sp.symbols("xi0:4"))
    ap = sp.Matrix(sp.symbols("ap0:4"))
    am = sp.Matrix(sp.symbols("am0:4"))
    q_vector = -HALF * (ap + h * kp * xi + am + h * km * xi)
    beta_vector = -HALF * (ap - h * kp * xi + am - h * km * xi)
    ell_vector = (trp + trm - sigma) * xi
    completed = q_vector + ell_vector - beta_vector
    expected = -h * israel * xi
    check_zero("komar_ghy_single_tension_completed_charge", completed - expected, results)
    check_detected("reverse_ghy_corner_sign_in_charge", q_vector + ell_vector + beta_vector - expected, controls)
    check_detected("two_tensions_instead_of_one", q_vector + (trp + trm - 2 * sigma) * xi - beta_vector - expected, controls)
    check_detected("omit_one_outward_gravitational_copy", -HALF * (ap + h * kp * xi) + (trp - sigma) * xi + HALF * (ap - h * kp * xi) - expected, controls)

    # At a shell-normal frame point, Codazzi itself supplies C_i,a; product
    # differentiation of sigma(phi) must cancel the stress exchange term.
    dsum_k = sp.Matrix(sp.symbols("codazzi0:4"))
    dphi = sp.Matrix(sp.symbols("dphi0:4"))
    sum_e_na = dsum_k - (pp + pm) * dphi
    divergence_israel = dsum_k + sigma_prime * dphi
    check_zero("codazzi_full_tangential_residual", sum_e_na - (divergence_israel - scalar * dphi), results)
    check_detected("omit_scalar_force_in_tangential_identity", sum_e_na - divergence_israel, controls)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="also save the exact JSON report to this path")
    parser.add_argument(
        "--inject-wrong-corner-sign", action="store_true",
        help="deliberately reverse the tested GHY corner sign; must fail with exit status 1",
    )
    options = parser.parse_args()
    identities = {}
    controls = {}
    try:
        gaussian_ghy_tensor_check(identities, controls, options.inject_wrong_corner_sign)
        junction_and_charge_checks(identities, controls)
    except RuntimeError as error:
        print(f"check failed: {error}", file=sys.stderr)
        return 1
    report = {
        "schema_version": 1,
        "arithmetic": "exact symbolic rational polynomial arithmetic",
        "dimension": {"bulk": 5, "shell": 4, "closed_cut": 3},
        "identities": identities,
        "wrong_controls": controls,
        "summary": {
            "identity_count": len(identities),
            "wrong_control_count": len(controls),
            "all_identities_zero": True,
            "all_wrong_controls_detected": True,
        },
        "scope": "Finite algebra checks only; no PDE existence or full symbolic verification of the differential-form proof.",
    }
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if options.output is not None:
        options.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
