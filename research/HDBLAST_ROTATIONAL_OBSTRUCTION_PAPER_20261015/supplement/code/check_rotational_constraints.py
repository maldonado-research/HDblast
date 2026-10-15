#!/usr/bin/env python3
"""Independent exact controls for the rotational shell constraint proof.

No scientific code from the public capsule is imported. Explicit exceptions,
not assert statements, enforce checks under normal and optimized Python.
JSON uses only local basenames; no network or repository mutation is needed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


checks: list[dict[str, object]] = []
controls: list[dict[str, object]] = []


def zero(name: str, expression: sp.Expr) -> None:
    residual = sp.simplify(sp.expand(expression))
    if residual != 0:
        raise RuntimeError(f"{name}: nonzero residual {residual}")
    checks.append({"name": name, "passed": True})


def nonzero(name: str, expression: sp.Expr) -> None:
    residual = sp.simplify(sp.expand(expression))
    if residual == 0:
        raise RuntimeError(f"{name}: wrong formula was not detected")
    controls.append({"name": name, "detected": True,
                     "exact_residual": str(residual)})


def nonzero_witness(name: str, expression: sp.Expr) -> None:
    residual = sp.simplify(sp.expand(expression))
    if residual == 0:
        raise RuntimeError(f"{name}: expected nonzero witness vanished")
    checks.append({"name": name, "passed": True,
                   "exact_nonzero_witness": str(residual)})


def sphere_integral(poly: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    """Integrate polynomial in four unit-sphere coordinates on S^3."""
    result = sp.S.Zero
    for exponents, coefficient in sp.Poly(poly, *variables).terms():
        if any(e % 2 for e in exponents):
            continue
        half = tuple(e // 2 for e in exponents)
        moment = 2 * sp.prod(sp.gamma(sp.Rational(1, 2) + a) for a in half)
        moment /= sp.gamma(2 + sum(half))
        result += coefficient * moment
    return sp.simplify(result)


def run() -> dict[str, object]:
    # The boost calculation is independent of the Cauchy reconstruction.
    eta, eta_j, shell_mixed = sp.symbols("eta eta_J K_tau_J", real=True)
    c, s = sp.cosh(eta), sp.sinh(eta)
    lorentz = sp.diag(-1, 1)
    tau = sp.Matrix([1, 0])
    normal = sp.Matrix([0, 1])
    u = c * tau - s * normal
    spatial_normal = c * normal - s * tau
    zero("tilt_u_unit", (u.T * lorentz * u)[0] + 1)
    zero("tilt_s_unit", (spatial_normal.T * lorentz * spatial_normal)[0] - 1)
    zero("tilt_u_s_orthogonal", (u.T * lorentz * spatial_normal)[0])
    derivative_u = sp.diff(u, eta) * eta_j
    # <nabla_J tau,n>=-K(tau,J), <nabla_J n,tau>=+K(tau,J).
    frame_mixed = -(c**2 - s**2) * shell_mixed
    coefficient_mixed = (derivative_u.T * lorentz * spatial_normal)[0]
    zero("tilted_cut_curvature", frame_mixed + coefficient_mixed
         + shell_mixed + eta_j)
    nonzero("wrong_boost_angle_sign", frame_mixed + coefficient_mixed
            + shell_mixed - eta_j)

    # R/2 fixes pi=sqrt(q)(k^ij-k q^ij)/2. All entries are computed.
    xs = sp.symbols("x1:5", real=True)
    x = sp.Matrix(xs)
    r2 = sum(v**2 for v in xs)
    Om, q, p, q2, p2, f = sp.symbols("Omega q p q2 p2 f", nonzero=True)
    identity = sp.eye(4)
    dq = 2 * q * (r2 * identity - 3 * x * x.T)
    dk = p / Om * (r2 * identity - 3 * x * x.T)
    background_inverse = identity / Om**2
    dk_trace = sp.trace(background_inverse * dk)
    zero("ADM_extrinsic_trace", dk_trace - r2 * p / Om**3)
    dpi = Om**4 / 2 * (background_inverse * dk * background_inverse
                       - dk_trace * background_inverse)
    dpi_expected = -3 * p / (2 * Om) * x * x.T
    for a in range(4):
        for b in range(4):
            zero(f"ADM_momentum_{a}{b}", dpi[a, b] - dpi_expected[a, b])
    contraction = sum(dpi[a, b] * dq[a, b] for a in range(4) for b in range(4))
    zero("ADM_gravity_pairing", contraction - 6 * r2**2 * p * q / Om)
    antisym = (contraction.subs(q, q2)
               - contraction.subs({p: p2}))
    zero("ADM_gravity_physical_orientation", antisym
         + 6 * r2**2 / Om * (q * p2 - q2 * p))
    nonzero("wrong_gravity_symplectic_sign", antisym
            - 6 * r2**2 / Om * (q * p2 - q2 * p))

    qD, pD, q2D, p2D = sp.symbols("Dq Dp Dq2 Dp2")
    chi1 = -3 * (qD + 4 * q) / (2 * Om**2 * f)
    chi2 = -3 * (q2D + 4 * q2) / (2 * Om**2 * f)
    chidot1 = -3 * (pD + 5 * p) / (2 * Om**2 * f)
    chidot2 = -3 * (p2D + 5 * p2) / (2 * Om**2 * f)
    scalar_can = Om**3 * (chidot1 * chi2 - chidot2 * chi1)
    scalar_phys = 9 / (4 * Om * f**2) * (
        (qD + 4 * q) * (p2D + 5 * p2)
        - (q2D + 4 * q2) * (pD + 5 * p))
    zero("ADM_scalar_physical_orientation", scalar_can + scalar_phys)
    wrong_chidot = -3 * (pD + 4 * p) / (2 * Om**2 * f)
    nonzero("missing_D_plus_5", Om**3 * (wrong_chidot - chidot1) * chi2)

    # Direct operator commutator, with enough polynomial dependence to fail.
    T = sp.symbols("T", real=True)
    test = T**3 * xs[0]**2 + T * xs[1] + xs[2]**4
    def euler(expr: sp.Expr) -> sp.Expr:
        return T * sp.diff(expr, T) + sum(v * sp.diff(expr, v) for v in xs)
    zero("ambient_E_commutator", sp.diff(euler(test), T)
         - euler(sp.diff(test, T)) - sp.diff(test, T))
    nonzero("missing_commutator_unit", sp.diff(euler(test), T)
            - euler(sp.diff(test, T)))

    # Finite-dimensional skew rotation checks moment map sign and factor.
    a, b, pa, pb = sp.symbols("a b pa pb")
    gen = sp.Matrix([[0, -1], [1, 0]])
    pos = sp.Matrix([a, b])
    mom = sp.Matrix([pa, pb])
    moment = (mom.T * gen * pos)[0]
    omega = (mom.T * gen * pos)[0] - ((gen * mom).T * pos)[0]
    zero("moment_map_Hessian_factor", omega - 2 * moment)
    nonzero("lost_Hessian_factor_two", omega - moment)

    # Exact S^3 harmonic and angular constants, independently integrated.
    n = sp.symbols("n1:5", real=True)
    H = n[0]**2 - n[1]**2
    def rotation(expr: sp.Expr) -> sp.Expr:
        return n[0] * sp.diff(expr, n[1]) - n[1] * sp.diff(expr, n[0])
    JH = rotation(H)
    J2H = rotation(JH)
    zero("degree_two_ambient_harmonic", sum(sp.diff(H, v, 2) for v in n))
    zero("degree_two_JH", JH + 4 * n[0] * n[1])
    zero("degree_two_J2H", J2H + 4 * H)
    volume = sphere_integral(sp.S.One, n)
    HH = sphere_integral(H**2, n)
    JHJH = sphere_integral(JH**2, n)
    HJH = sphere_integral(H * JH, n)
    HJ2H = sphere_integral(H * J2H, n)
    zero("sphere_volume", volume - 2 * sp.pi**2)
    zero("sphere_HH", HH - sp.pi**2 / 3)
    zero("sphere_JHJH", JHJH - 4 * sp.pi**2 / 3)
    zero("sphere_HJH", HJH)
    zero("sphere_angular_integration_by_parts", HJ2H + JHJH)
    w, Krad = sp.symbols("omega K_rad", nonzero=True)
    physical_one = w * Krad * (HJ2H - JHJH)
    physical_double = 2 * physical_one
    Hessian_double = 4 * w * Krad * JHJH
    zero("one_interior_rotational_candidate", physical_one
         + 2 * w * Krad * JHJH)
    zero("two_interior_constraint_orientation", Hessian_double + physical_double)
    zero("explicit_l2_Hessian", Hessian_double - 16 * sp.pi**2 * w * Krad / 3)
    nonzero_witness("rotational_second_order_obstruction", Hessian_double)

    # Radial IBP reduces the 9 Q DQ cross term. No background is frozen.
    radial_r, L = sp.symbols("r L", positive=True)
    rm_log_derivative = 4 - 2 * radial_r**2 * L
    radial_coefficient = 20 - sp.Rational(9, 2) * rm_log_derivative
    zero("radial_diagonal_coefficient", radial_coefficient - (2 + 9 * radial_r**2 * L))
    wrong_diagonal = 16 - 4 * rm_log_derivative
    nonzero("radial_missing_velocity_unit", wrong_diagonal - radial_coefficient)

    # Nonlinear normal shell force identity from full Gauss constraint.
    h = sp.diag(-1, 1, 1, 1)
    symbols_plus = sp.symbols("kp0:10")
    symbols_minus = sp.symbols("km0:10")
    def symmetric(entries: tuple[sp.Symbol, ...]) -> sp.Matrix:
        out = sp.zeros(4)
        counter = 0
        for i in range(4):
            for j in range(i, 4):
                out[i, j] = out[j, i] = entries[counter]
                counter += 1
        return out
    Kp, Km = symmetric(symbols_plus), symmetric(symbols_minus)
    trp, trm = sp.trace(h * Kp), sp.trace(h * Km)
    sigma, sigmap, pp, pm = sp.symbols("sigma sigma_prime p_plus p_minus")
    Israel = Kp + Km - (trp + trm) * h + sigma * h
    scalar_residual = pp + pm + sigmap
    delta_K = Kp - Km
    double_E_difference = (trp**2 - trm**2
                           - sp.trace(h * Kp * h * Kp)
                           + sp.trace(h * Km * h * Km)
                           - pp**2 + pm**2)
    mixed = sp.trace(h * Israel * h * delta_K)
    force = sigma * (trp - trm) + sigmap * (pp - pm)
    zero("nonlinear_force_full_Gauss_identity", force
         - double_E_difference - mixed - scalar_residual * (pp - pm))
    nonzero("missing_scalar_force_residual", force - double_E_difference - mixed)
    nonzero("double_counted_shell_tension", 2 * sigma * (trp - trm)
            + sigmap * (pp - pm) - double_E_difference - mixed
            - scalar_residual * (pp - pm))
    # Tangential tension projection is exactly zero, independently of sigma.
    tau = sp.Matrix([1, 0, 0, 0])
    J = sp.Matrix([0, 1, 0, 0])
    zero("single_tension_mixed_projection", sigma * (tau.T * h * J)[0])
    zero("Israel_mixed_projection", (tau.T * Israel * J)[0]
         - (tau.T * (Kp + Km) * J)[0])

    # Closed-cut controls: constant background area density kills arbitrary w.
    theta = sp.symbols("theta", real=True)
    area0 = sp.symbols("area0", positive=True)
    w_eta = sp.sin(theta) + 3 * sp.cos(2 * theta)
    zero("arbitrary_second_order_cut_angle", sp.integrate(
        area0 * sp.diff(w_eta, theta), (theta, 0, 2 * sp.pi)))
    amp, tilt = sp.symbols("area_amplitude tilt_amplitude")
    eps = sp.symbols("epsilon")
    moving_area = area0 + eps * amp * sp.cos(theta)
    moving_eta = eps * tilt * sp.sin(theta) + eps**2 * w_eta / 2
    corner_moment = sp.integrate(moving_area * sp.diff(moving_eta, theta),
                                 (theta, 0, 2 * sp.pi))
    zero("corner_background_first_derivative", sp.diff(corner_moment, eps).subs(eps, 0))
    zero("corner_Hessian_with_second_order_correction",
         sp.diff(corner_moment, eps, 2).subs(eps, 0) - 2 * sp.pi * amp * tilt)
    zero("compact_first_order_corner_Hessian", sp.diff(corner_moment, eps, 2)
         .subs({eps: 0, amp: 0, tilt: 0}))
    nonzero("omitted_tilted_cut_angular_term", sp.integrate(
        (1 + sp.cos(theta)) * sp.diff(sp.sin(theta), theta),
        (theta, 0, 2 * sp.pi)))

    source = Path(__file__).resolve()
    note = source.with_name("ROTATIONAL_SECOND_ORDER_SHELL_OBSTRUCTION.md")
    return {
        "schema": "independent-rotational-shell-constraints-v1",
        "scientific_status": "exact_algebra_controls_passed",
        "scope": "spatial constraint Hessian; PDE and differentiability assumptions remain analytic",
        "sympy_version": sp.__version__,
        "checks": checks,
        "wrong_controls": controls,
        "counts": {"checks": len(checks), "wrong_controls": len(controls)},
        "exact_results": {
            "cut_normal_identity": "k_i(s_i,J) = -K_i(tau,J) - J eta_i",
            "full_spatial_constraint": "C_J = M_J + integral_S sqrt(c) J(eta_plus+eta_minus)",
            "Hessian_orientation": "D2 C_J(v,v) = -Omega_phys(v,Jv) for compact initial data",
            "sphere_HH": str(HH), "sphere_JHJH": str(JHJH),
            "two_interior_l2_Hessian": str(Hessian_double),
        },
        "files": [
            {"name": source.name, "sha256": hashlib.sha256(source.read_bytes()).hexdigest()},
            {"name": note.name, "sha256": hashlib.sha256(note.read_bytes()).hexdigest()},
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--inject-failure", action="store_true")
    args = parser.parse_args()
    try:
        result = run()
        if args.inject_failure:
            zero("explicit_injected_failure", sp.S.One)
    except Exception as error:
        failure = {"schema": "independent-rotational-shell-constraints-v1",
                   "scientific_status": "failed", "error": str(error),
                   }
        rendered = json.dumps(failure, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.write_text(rendered)
        else:
            sys.stdout.write(rendered)
        return 1
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered)
    else:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
