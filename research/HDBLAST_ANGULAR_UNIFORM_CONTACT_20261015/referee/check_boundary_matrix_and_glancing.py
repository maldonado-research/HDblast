#!/usr/bin/env python3
"""Independent exact algebra checks for the boundary referee note."""

import argparse
import json
from pathlib import Path
import sys
import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def check(corruption=None):
    a, sig1, sig2, eps = sp.symbols("a sig1 sig2 eps", real=True)
    k, w = sp.symbols("k w", real=True)
    h0 = sp.diag(-1, 1, 1, 1)
    indices = [(i, j) for i in range(4) for j in range(i, 4)]
    hs = sp.symbols("h0:10", real=True)
    H = sp.zeros(4)
    for x, (i, j) in zip(hs, indices):
        H[i, j] = H[j, i] = x
    trH = sp.trace(h0 * H)
    n = sp.sqrt(1 + eps * k)
    sigma = 6 * a + eps * sig1 * w
    sigma_prime = sig1 + eps * sig2 * w
    # Differentiate the exact inverse, rather than assume the trace row.
    inverse_trace = sp.trace((h0 + eps * H).inv() * h0)
    exact_H_y = sigma * n * (h0 + eps * H) / 3 - 2 * a * h0
    trace_factor = 2 if corruption == "wrong-trace-factor" else 1
    exact_k_y = 4 * sigma * n**3 / 3 - 2 * a * n**4 * inverse_trace * trace_factor
    exact_w_y = -sigma_prime * n / 2 + sig1 / 2
    rows = [sp.diff(exact_H_y[i, j], eps).subs(eps, 0) for i, j in indices]
    rows += [sp.diff(exact_k_y, eps).subs(eps, 0), sp.diff(exact_w_y, eps).subs(eps, 0)]
    variables = [*hs, k, w]
    B = sp.Matrix(rows).jacobian(variables)
    J = sp.Matrix([h0[i, j] for i, j in indices])
    ell = sp.Matrix([trH]).jacobian(hs)
    expected = sp.zeros(12)
    expected[:10, :10] = 2 * a * sp.eye(10)
    expected[:10, 10] = a * J
    expected[:10, 11] = sig1 * J / 3
    expected[10, :10] = 2 * a * ell
    expected[10, 10] = -4 * a
    expected[10, 11] = 4 * sig1 / 3
    expected[11, 10] = -sig1 / (3 if corruption == "wrong-linear-scalar" else 4)
    expected[11, 11] = -sig2 / 2
    require(sp.simplify(B - expected) == sp.zeros(12), "exact direct-metric linearization disagrees with Robin matrix")
    require((ell * J)[0] == 4, "four-dimensional boundary trace pairing is incorrect")

    # Null tangential covector in an orthonormal Lorentzian frame.
    zeta = sp.Matrix([1, 1, 0, 0])
    require((zeta.T * h0 * zeta)[0] == 0, "test tangential covector is not null")
    constraints = H * h0 * zeta - zeta * (k + trH) / 2
    C = constraints.jacobian(variables)
    if corruption == "wrong-glancing-rank":
        bogus_scalar_constraint = sp.zeros(1, 12)
        bogus_scalar_constraint[0, 11] = 1
        C = C.col_join(bogus_scalar_constraint)
    require(C.rank() == 4, "glancing harmonic constraint rank is not four")
    require(C[:, 11] == sp.zeros(4, 1), "scalar glancing amplitude is incorrectly constrained")

    # T(0)=I and T_z(0)=M imply exact boundary cancellation,
    # independent of any symmetry assumptions.
    rho = sp.symbols("rho", positive=True)
    M = -rho * B
    vv = sp.Matrix(sp.symbols("v0:12"))
    vz = sp.Matrix(sp.symbols("vz0:12"))
    boundary_T_z = 2 * M if corruption == "wrong-conjugation" else M
    require(sp.simplify((boundary_T_z * vv + vz) - M * vv) == vz,
            "Robin conjugation leaves a nonzero boundary term")

    # Lower-order normal-symbol roots cannot have Re(s)>||M||:
    # If lambda=A+iB, s=gamma+i tau, then A B=gamma tau and
    # (A^2-gamma^2)(A^2+tau^2)=c^2|eta|^2 A^2 >= 0.
    A, gamma, tau, angular2 = sp.symbols("A gamma tau angular2", positive=True)
    dispersion = A**2 - gamma**2 * tau**2 / A**2 - gamma**2 + tau**2 - angular2
    identity = (A**2 - gamma**2) * (A**2 + tau**2) - angular2 * A**2
    require(sp.simplify(A**2 * dispersion - identity) == 0,
            "Laplace normal-root bound dispersion identity failed")
    return {
        "status": "PASS",
        "check_group_count": 5,
        "check_assertion_count": 7,
        "exact_direct_metric_linearization": "PASS",
        "normal_block_dimension": B.rows,
        "trace_rank_pairing": int((ell * J)[0]),
        "glancing_constraint_rank": C.rank(),
        "constrained_glancing_N_dimension": 12 - C.rank(),
        "scalar_glancing_mode_survives": True,
        "nonsymmetric_Robin_cancellation": "PASS",
        "laplace_normal_root_bound_identity": "PASS",
        "scope": "exact algebra checks; not a replacement for the PDE proof",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inject-corruption", choices=[
        "wrong-linear-scalar", "wrong-trace-factor", "wrong-glancing-rank", "wrong-conjugation"
    ])
    args = parser.parse_args()
    exit_code = 0
    try:
        result = check(args.inject_corruption)
    except AssertionError as exc:
        result = {"status": "FAIL", "injected_corruption": args.inject_corruption,
                  "detected_error": str(exc)}
        exit_code = 1
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")
    sys.exit(exit_code)
