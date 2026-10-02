#!/usr/bin/env python3
"""Exact symbolic stationary-model checks; no numerical sources or ODE solve.

Requires SymPy. Every test remains active with Python optimization enabled.
An explicit, unused output file is required so replay preserves prior evidence.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import sympy as s


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite existing evidence: {args.output}")

    checks: list[dict] = []
    controls: list[dict] = []

    def equal(name, actual, expected=0):
        residual = s.simplify(s.expand(actual - expected))
        entry = {"name": name, "pass": residual == 0, "residual": s.sstr(residual)}
        checks.append(entry)
        if residual != 0:
            raise RuntimeError(f"Exact identity failed: {name}: {residual}")

    def wrong(name, actual, correct=0):
        defect = s.simplify(s.expand(actual - correct))
        entry = {"name": name, "detected": defect != 0, "defect": s.sstr(defect)}
        controls.append(entry)
        if defect == 0:
            raise RuntimeError(f"Wrong-formula control escaped detection: {name}")

    z, x, r, L, kap, gamma = s.symbols("z x r L kappa5 gamma", positive=True)
    phi, phi0, b, delta, c = s.symbols("phi phi0 b delta c", real=True)
    q, w, eta = s.symbols("q w eta", real=True)
    C = 1 / (16 * s.pi**2)
    W = s.Function("W")(z, x, r)
    Q = 2 * s.diff(W, x)
    rho = W - z * s.diff(W, z) / 2
    mass_factor = 1 + b * (phi - phi0) / 2
    xr = r * mass_factor**2
    xp = b * r * mass_factor
    xpp = b**2 * r / 2
    j = xp * Q / 2

    # Normalization: the action has -(grad phi)^2/(2 kappa5^2).
    gradPhi = s.Symbol("gradPhi", real=True)
    equal("canonical_bulk_scalar_no_sqrt2", -(kap * gradPhi)**2 / (2 * kap**2), -gradPhi**2 / 2)
    Rho, Cur = s.symbols("rho j", real=True)
    equal("metric_source_dimension", L * (gamma * L**3) * (Rho / L**4), gamma * Rho)
    equal("scalar_source_dimension", L * (gamma * L**3) * (Cur / L**4), gamma * Cur)
    equal("mass_current_dimension", (xp / L**2) * (Q / L**2) / 2, j / L**4)
    equal("mass_law_first_derivative", s.diff(xr, phi), xp)
    equal("mass_law_second_derivative", s.diff(xr, phi, 2), xpp)
    equal("mass_law_reference_mass", xr.subs(phi, phi0), r)
    equal("mass_law_reference_slope", xp.subs(phi, phi0), b * r)
    equal("quadratic_original_interaction", xr, b**2 * r * (phi - phi0 + 2 / b)**2 / 4)
    # For b != 0 and fixed z>0, insert only the known leading zero-mode term.
    distance = s.Symbol("phi_minus_phi_chi", nonzero=True)
    mass_near_zero = r * b**2 * distance**2 / 4
    slope_near_zero = r * b**2 * distance / 2
    Q_zero_mode = 3 * z**2 / (8 * s.pi**2 * mass_near_zero)
    equal("quadratic_mass_zero_current_leading_coefficient",
          distance * slope_near_zero * Q_zero_mode / 2, 3 * z**2 / (8 * s.pi**2))

    equal("common_action_pairing", s.diff(rho, x), Q / 2 - z * s.diff(Q, z) / 4)
    equal("mass_law_current_pairing", xp * s.diff(rho, x), j - z * s.diff(j, z) / 2)
    equal("current_phi_full_chain", s.diff(j, phi) + xp * s.diff(j, x), (xpp * Q + xp**2 * s.diff(Q, x)) / 2)
    equal("density_z_from_action", s.diff(rho, z), s.diff(W, z) / 2 - z * s.diff(W, z, 2) / 2)
    equal("current_z_from_action", s.diff(j, z), xp * s.diff(Q, z) / 2)

    u = x / z
    v = r / z
    Pfun = s.Function("Psi")
    Avec = Pfun(u) - s.log(v)
    Qc = C * z * ((u - 2) * Avec - u + v + s.Rational(4, 3))
    # Direct derivative of Psi(x/z), expressed without any numerical special function.
    pux = z * s.diff(Pfun(u), x)
    Qx = C * (Avec + (u - 2) * pux - 1)
    Qz = C * (-2 * Avec - u * (u - 2) * pux + u - s.Rational(2, 3))
    equal("closed_current_x_derivative", s.diff(Qc, x), Qx)
    equal("closed_current_z_derivative_fixed_reference", s.diff(Qc, z), Qz)
    poly = (x - r)**2 / 2 - 2 * z * (x - r) + s.Rational(29, 15) * z**2
    Du = u**2 / 2 - 2 * u + s.Rational(29, 15)
    rhoc = C * z**2 / 4 * (u * (u - 2) * Avec - u**2 + s.Rational(4, 3) * u
                            + 2 * u * v - 2 * v - v**2 / 2 - Du)
    equal("closed_source_trace_check", 4 * rhoc, x * Qc - C * poly)
    equal("closed_source_pairing_check", s.diff(rhoc, x), Qc / 2 - z * Qz / 4)
    equal("closed_density_z_derivative", s.diff(rhoc, z), (x * Qz - C * (-2 * (x - r) + s.Rational(58, 15) * z)) / 4)
    equal("exact_reference_current", Qc.subs({x: 2 * z, r: 2 * z}, simultaneous=True), z / (12 * s.pi**2))
    equal("exact_reference_density", rhoc.subs({x: 2 * z, r: 2 * z}, simultaneous=True), 11 * z**2 / (960 * s.pi**2))
    equal("exact_reference_scalar_current", (xp * Qc / 2).subs(phi, phi0).subs({x: 2 * z, r: 2 * z}, simultaneous=True), b * z**2 / (12 * s.pi**2))

    # Differentiate the inherited determinant normalization through its Hessian.
    zz, xx = s.symbols("zz xx", positive=True)
    qq = xx / zz - s.Rational(9, 4)
    B1fun = s.Function("B1")
    B00 = s.Symbol("B00")
    B0 = B00 + qq / 8 + qq**2 / 4
    P = B1fun(qq) - B0 * s.log(zz)
    Pq = zz * s.diff(P, xx)
    Pqq = zz**2 * s.diff(P, xx, 2)
    B0q = s.Rational(1, 8) + qq / 2
    Det = -C * zz**2 * P + C * (r * xx - r**2 / 4 - 2 * r * zz
          - (xx**2 / 2 - 2 * xx * zz + s.Rational(29, 15) * zz**2) * s.log(r)) / 2
    uu = xx / zz
    equal("determinant_Qx_hessian", 2 * s.diff(Det, xx, 2), -2 * C * Pqq - C * s.log(r))
    equal("determinant_Qz_hessian", 2 * s.diff(Det, xx, zz), -2 * C * (Pq - uu * Pqq - B0q) + 2 * C * s.log(r))
    Wzz_expected = -C * (2 * P - 2 * uu * Pq - 3 * B0 + uu**2 * Pqq + 2 * uu * B0q) - C * s.Rational(29, 15) * s.log(r)
    equal("determinant_Wzz_hessian", s.diff(Det, zz, 2), Wzz_expected)

    B = 1 - phi + phi**3 / 3
    Bphi = s.diff(B, phi)
    U = Bphi**2 / 2 - 2 * B**2 / 3
    sig = 2 * B + delta * (1 + c * phi)
    S = sig + gamma * Rho
    T = s.diff(sig, phi) + gamma * Cur
    E1 = q - S / 6
    E2 = w + T / 2
    Dconstraint = z + w**2 / 12 - U / 6 - S**2 / 36
    tau = delta * (1 + c * phi) + gamma * Rho
    Dstable = z + (w - Bphi) * (w + Bphi) / 12 - B * tau / 9 - tau**2 / 36
    equal("stable_factored_metric_residual", Dconstraint, Dstable)
    equal("shifted_superpotential", B.subs(phi, 1 + eta), s.Rational(1, 3) + eta**2 + eta**3 / 3)
    equal("shifted_superpotential_derivative", Bphi.subs(phi, 1 + eta), 2 * eta + eta**2)
    equal("squared_unsquared_relation", q**2 - S**2 / 36, (q + S / 6) * E1)
    equal("endpoint_constraint_from_both_junctions", Dconstraint.subs(w, -T / 2), z - S**2 / 36 + T**2 / 48 - U / 6)

    qa, wa, phia, za, rhoa, ja = s.symbols("q_a w_a phi_a z_a rho_a j_a")
    rz, rp, jz, jp = s.symbols("rho_z rho_phi j_z j_phi")
    sigp, sigpp, Up = s.symbols("sigma_phi sigma_phiphi U_phi")
    chainrho = rz * za + rp * phia
    chainj = jz * za + jp * phia
    Sa = sigp * phia + gamma * chainrho
    M1 = qa - Sa / 6
    M2 = wa + sigpp * phia / 2 + gamma * chainj / 2
    equal("full_metric_chain_rule", M1, qa - sigp * phia / 6 - gamma * chainrho / 6)
    equal("full_scalar_chain_rule", M2, wa + sigpp * phia / 2 + gamma * chainj / 2)
    Ss = s.Symbol("S")
    equal("squared_residual_offshell_jacobian", 2 * q * qa - Ss * Sa / 18,
          (q + Ss / 6) * M1 + (q - Ss / 6) * (qa + Sa / 6))
    equal("squared_residual_onshell_row_scale", (2 * q * qa - Ss * Sa / 18).subs(Ss, 6 * q), 2 * q * M1)
    E2vac = w + sigp / 2
    M1y = M1.subs({qa: -z - w**2 / 3, za: -2 * z * q, phia: w})
    equal("moving_endpoint_metric_column", M1y, -z - w * E2vac / 3 - gamma * (-2 * z * q * rz + w * rp) / 6)
    equal("moving_endpoint_scalar_column", M2.subs({wa: Up - 4 * q * w, za: -2 * z * q, phia: w}),
          Up - 4 * q * w + sigpp * w / 2 + gamma * (-2 * z * q * jz + w * jp) / 2)
    equal("sourced_endpoint_metric_pairing", M1y.subs({sigp: -2 * w - gamma * Cur, rp: Cur - z * jz / 2}),
          -z + gamma * z * q * rz / 3 + gamma * z * w * jz / 12)

    k, vv, f, fp, RR, II = s.symbols("k v f fp R I")
    weighted = 2 * z * k - 2 * w * fp / 3 + (Up - 4 * q * w) * f / 3 + w * fp / 3 + 4 * q * (vv + w * f / 3)
    # Differentiated constraint is w*f'=12*q*v+12*z*k+U_phi*f.
    # A dummy w*f' keeps this a polynomial identity, including w=0.
    wfp = s.Symbol("w_times_fprime")
    weighted_polynomial = s.expand(weighted).subs(w * fp, wfp)
    equal("bulk_weighted_variation_identity", weighted_polynomial.subs(wfp, 12 * q * vv + 12 * z * k + Up * f), -2 * z * k)
    equal("bulk_weighted_integral_derivative", (-2 * RR**4 * z * k).subs(z, 1 / RR**2), -2 * RR * (RR * k))
    equal("offshell_ell_metric_derivative", (-2 * II / RR**4 - w * f / 3) - sigp * f / 6,
          -2 * II / RR**4 - E2vac * f / 3)

    # Differentiate the full sourced endpoint constraint, not a frozen source.
    zg, pg, gg = s.symbols("dz dphi dgamma")
    Ts, Us = s.symbols("T U")
    dS = sigp * pg + gamma * (rz * zg + rp * pg) + Rho * gg
    dT = sigpp * pg + gamma * (jz * zg + jp * pg) + Cur * gg
    dg = zg - Ss * dS / 18 + Ts * dT / 24 - Up * pg / 6
    expected_dg = (1 - gamma * (Ss * rz / 18 - Ts * jz / 24)) * zg \
        - (Ss * (sigp + gamma * rp) / 18 - Ts * (sigpp + gamma * jp) / 24 + Up / 6) * pg \
        - (Ss * Rho / 18 - Ts * Cur / 24) * gg
    equal("full_endpoint_constraint_differential", dg, expected_dg)

    V = (x**2 * s.log(x / r) - s.Rational(3, 2) * x**2 + 2 * r * x - r**2 / 2) / (64 * s.pi**2)
    F = (x * s.log(x / r) - x + r) / (192 * s.pi**2)
    for order in range(3):
        equal(f"fixed_reference_V_matching_order_{order}", s.diff(V, x, order).subs(x, r))
    for order in range(2):
        equal(f"fixed_reference_F_matching_order_{order}", s.diff(F, x, order).subs(x, r))
    alpha = s.Function("alpha")(x)
    Wlocal = V - 12 * z * F + alpha * z**2
    equal("local_density_curvature_variation", Wlocal - z * s.diff(Wlocal, z) / 2, V - 6 * z * F)
    equal("curvature_squared_scalar_variation", 2 * s.diff(alpha * z**2, x), 2 * z**2 * s.diff(alpha, x))
    kT = U / 6 - w**2 / 12
    kN = U / 6 + w**2 / 4
    equal("bulk_Ricci_trace_from_sectional_curvatures", 12 * kT + 8 * kN, w**2 + s.Rational(10, 3) * U)
    Ee = s.Symbol("E", positive=True)
    equal("N1_gravity_power_counting", (Ee / gamma**(-s.Rational(1, 3)))**3, gamma * Ee**3)

    # Wrong-formula controls are symbolic nonidentities, not numerical probes.
    wrong("omit_metric_radius_variation", W, rho)
    wrong("treat_current_as_density_phi_derivative", j, xp * s.diff(rho, x))
    wrong("track_reference_with_mass_in_scalar_variation", Q + 2 * s.diff(W, r), Q)
    wrong("track_reference_with_curvature_in_metric_variation", rho - r * s.diff(W, r) / 2, rho)
    wrong("omit_quadratic_second_mass_derivative", xp**2 * s.diff(Q, x) / 2, s.diff(j, phi) + xp * s.diff(j, x))
    wrong("substitute_exponential_mass_hessian_at_reference", b**2 * r, xpp)
    wrong("freeze_quantum_density_in_jacobian", qa - sigp * phia / 6, M1)
    wrong("omit_quantum_scalar_jacobian", wa + sigpp * phia / 2, M2)
    wrong("scale_only_density_source", w + sigp / 2, w + sigp / 2 + gamma * Cur / 2)
    wrong("reverse_current_sign", w + sigp / 2 - gamma * Cur / 2, w + sigp / 2 + gamma * Cur / 2)
    wrong("omit_moving_endpoint_metric_column", 0, M1y)
    wrong("omit_offshell_squared_jacobian_term", (q + Ss / 6) * M1, 2 * q * qa - Ss * Sa / 18)
    wrong("accept_negative_squared_root", (q - Ss / 6).subs(q, -Ss / 6), 0)
    wrong("set_vacuum_scalar_residual_zero_at_sourced_root", -2 * II / RR**4,
          -2 * II / RR**4 + gamma * Cur * f / 6)
    wrong("freeze_sources_in_endpoint_constraint", zg - Ss * (sigp * pg + Rho * gg) / 18
          + Ts * (sigpp * pg + Cur * gg) / 24 - Up * pg / 6, dg)
    wrong("double_count_local_potential_density", 2 * V - 6 * z * F, V - 6 * z * F)
    wrong("drop_curvature_squared_scalar_source", 0, 2 * z**2 * s.diff(alpha, x))
    wrong("omit_induced_Einstein_scalar_current", 2 * s.diff(V, x), 2 * s.diff(V - 12 * z * F, x))
    wrong("omit_determinant_log_in_second_mass_derivative", -2 * C * (Pqq + s.log(zz) / 2) - C * s.log(r), -2 * C * Pqq - C * s.log(r))
    wrong("wrong_bulk_scalar_sqrt2_normalization", -(s.sqrt(2) * kap * gradPhi)**2 / (2 * kap**2), -gradPhi**2 / 2)

    report = {
        "status": "PASS",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "python_optimization": sys.flags.optimize,
        "sympy": s.__version__,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": checks,
        "negative_controls": controls,
        "scope": "Exact symbolic identities and deliberate nonidentities only; no numerical source evaluation, radial integration, root search, or existence certification.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(f"PASS: {len(checks)} exact identities; {len(controls)} wrong-formula controls detected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
