#!/usr/bin/env python3
"""Exact algebra for the pre-contact dynamical continuation.

Imports no checkpoint scientific code. The checks concern ADM conversion,
initial-time kinetic forcing, an SO(4) tensor-sector witness, and causal
support geometry. They do NOT test a numerical Einstein-shell evolution.
All validation uses explicit exceptions, so python -O changes no checks.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path)
    ap.add_argument("--inject-failure", action="store_true")
    args = ap.parse_args()
    checks = []
    controls = []

    def zero(z):
        if isinstance(z, s.MatrixBase):
            return all(s.simplify(v) == 0 for v in z)
        if isinstance(z, (list, tuple)):
            return all(zero(v) for v in z)
        return s.simplify(z) == 0

    def check(name, residual):
        ok = zero(residual)
        checks.append({"name": name, "passed": ok})
        if not ok:
            raise RuntimeError(f"positive check failed: {name}")

    def control(name, wrong_residual):
        detected = not zero(wrong_residual)
        controls.append({"name": name, "detected": detected})
        if not detected:
            raise RuntimeError(f"wrong formula escaped: {name}")

    # Independently derive the time-only spacetime Ricci tensor. This checks
    # the universal kinetic terms of ADM, without an inherited ADM engine.
    # Here gamma_i = 1 + 2 eps k_i tau + eps^2 a_i tau^2/4; thus
    # d_eps^2 d_tau^2 gamma_i|0 = a_i.
    t, eps = s.symbols("tau eps")
    k = s.symbols("k1:5")
    a = s.symbols("a1:5")
    q = [1 + 2 * eps * k[i] * t + eps**2 * a[i] * t**2 / 4 for i in range(4)]
    metric = [-s.Integer(1)] + q
    inv = [1 / z for z in metric]

    def d(z, mu):
        return s.diff(z, t) if mu == 0 else s.Integer(0)

    def g(i, j):
        return metric[i] if i == j else s.Integer(0)

    def christ(i, j, h):
        return inv[i] * (d(g(i, h), j) + d(g(i, j), h) - d(g(j, h), i)) / 2

    ricci = []
    for i in range(1, 5):
        z = sum(d(christ(h, i, i), h) - d(christ(h, i, h), i)
                + sum(christ(h, h, l) * christ(l, i, i)
                      - christ(h, i, l) * christ(l, i, h) for l in range(5))
                for h in range(5))
        ricci.append(s.simplify(z.subs(t, 0)))
    predicted_ricci = [eps**2 * (a[i] / 4 + sum(k) * k[i] - 2 * k[i]**2)
                       for i in range(4)]
    check("direct_5d_ricci_kinetic_identity", [ricci[i] - predicted_ricci[i] for i in range(4)])
    accel = [8 * k[i]**2 - 4 * sum(k) * k[i] for i in range(4)]
    check("metric_second_variation_acceleration", [ricci[i].subs(a[i], accel[i]) for i in range(4)])
    control("wrong_second_variation_factor", [ricci[i].subs(a[i], accel[i] / 2) for i in range(4)])
    control("wrong_trace_K_sign", [ricci[i].subs(a[i], 8 * k[i]**2 + 4 * sum(k) * k[i]) for i in range(4)])

    u, ac = s.symbols("u scalar_acceleration")
    ph = eps * u * t + eps**2 * ac * t**2 / 4
    # The scalar d'Alembertian is evaluated directly using the Christoffels.
    box = -s.diff(ph, t, 2) - sum(inv[i] * christ(0, i, i) * s.diff(ph, t) for i in range(1, 5))
    box0 = s.simplify(box.subs(t, 0))
    check("direct_scalar_kinetic_identity", box0 + eps**2 * (ac / 2 + sum(k) * u))
    check("scalar_second_variation_acceleration", box0.subs(ac, -2 * sum(k) * u))
    control("omit_scalar_Ku_term", box0.subs(ac, 0))

    # Densitized ADM reconstruction, independently converted to K and u.
    x = s.Matrix(s.symbols("x1:5"))
    rad = (x.T * x)[0]
    om, p, f, kr, kp = s.symbols("Omega p Phi_prime K K_prime", nonzero=True)
    gamma = om**2 * s.eye(4)
    pi = -3 * p * x * x.T / (2 * om)
    J = om**4
    trpi = s.trace(gamma * pi)
    K1 = 2 / J * (gamma * pi * gamma - gamma * trpi / 3)
    Kexpected = p / om * (rad * s.eye(4) - 3 * x * x.T)
    check("density_to_extrinsic_curvature", K1 - Kexpected)
    trK = s.trace(gamma.inv() * K1)
    check("quadrupole_trace_K", trK - rad * p / om**3)
    forcing = 8 * K1 * gamma.inv() * K1 - 4 * trK * K1
    target = 4 * rad * p**2 / om**4 * (rad * s.eye(4) + 9 * x * x.T)
    check("closed_quadrupole_metric_forcing", forcing - target)
    control("wrong_density_Omega_power", K1 - p / om**2 * (rad * s.eye(4) - 3 * x * x.T))
    H = x[0]**2 - x[1]**2
    u1 = -3 * (2 * rad * kp + 7 * kr) * H / (2 * om**3 * f)
    scalarforcing = -2 * trK.subs(p, kr * H) * u1
    scalar_target = 3 * rad * kr * (2 * rad * kp + 7 * kr) * H**2 / (om**6 * f)
    check("closed_quadrupole_scalar_forcing", scalarforcing - scalar_target)

    # An explicit transverse, trace-free tensor harmonic witness on S^3.
    lam = [1, -1, 0, 0]
    sectional = {(i, j): s.Rational(lam[i] * lam[j])
                 + s.Rational(lam[i]**2 + lam[j]**2, 2) - s.Rational(1, 3)
                 for i in range(4) for j in range(i + 1, 4)}
    def c(i, j):
        return sectional[tuple(sorted((i, j)))] if i != j else s.Integer(0)
    TT = s.Matrix(4, 4, lambda i, j: sum(c(i, h) * x[h]**2 for h in range(4))
                  if i == j else -c(i, j) * x[i] * x[j])
    check("weyl_tensor_tangency_and_trace", [TT * x, s.trace(TT)])
    proj = s.eye(4) - x * x.T / rad
    div = s.Matrix([sum(proj[j, h] * s.diff(TT[i, j], x[h])
                         for j in range(4) for h in range(4)) for i in range(4)])
    check("weyl_tensor_sphere_divergence", proj * div)
    gradH = 2 * s.diag(*lam) * x - 2 * H * x / rad
    contraction = s.expand(s.cancel((gradH.T * TT * gradH)[0]))
    predicted = 4 * (-s.Rational(4, 3) * x[0]**2 * x[1]**2
                     + s.Rational(1, 6) * (x[0]**2 + x[1]**2) * (x[2]**2 + x[3]**2))
    check("tensor_stress_projection_polynomial", contraction - predicted)
    # Exact uniform S^3 moment of a degree-four polynomial.
    def sphere_average(poly):
        val = s.Integer(0)
        for powers, coeff in s.Poly(s.expand(poly), *x).terms():
            if any(j % 2 for j in powers):
                continue
            n = sum(powers) // 2
            den = s.prod(4 + 2 * j for j in range(n))
            num = s.prod(s.factorial2(j - 1) if j else 1 for j in powers)
            val += coeff * num / den
        return s.simplify(val)
    ttmean = sphere_average(contraction)
    check("nonzero_tensor_projection_exact_sphere_moment", ttmean + s.Rational(1, 9))
    control("discard_all_tensor_harmonics", ttmean)
    H2 = x[0]**2 + x[1]**2 - rad / 2
    H4 = H**2 - rad * (x[0]**2 + x[1]**2) / 2 + rad**2 / 12
    def tangential_tf(M):
        M = proj * M * proj
        return M - s.trace(M) * proj / 3
    # Hessians of homogeneous ambient harmonics restrict to sphere Hessians;
    # the pure-trace normal term disappears in their trace-free parts.
    scalar2 = tangential_tf(s.hessian(H2, x))
    scalar4 = tangential_tf(s.hessian(H4, x)) / rad
    stress = 4 / rad * (s.diag(*lam) * x - H * x / rad) * (s.diag(*lam) * x - H * x / rad).T
    stress_tf = stress - s.trace(stress) * proj / 3
    check("exact_scalar_stress_tensor_decomposition",
          stress_tf - scalar2 / 2 - scalar4 / 3 + 4 * TT / (3 * rad))
    ttnorm = sphere_average(s.trace(TT * TT))
    check("tensor_witness_exact_sphere_norm", ttnorm - s.Rational(1, 12))
    check("scalar_quadrupole_square_decomposition",
          H**2 - ((H**2 - rad * (x[0]**2 + x[1]**2) / 2 + rad**2 / 12)
                  + rad * (x[0]**2 + x[1]**2 - rad / 2) / 2 + rad**2 / 6))

    b, upper, lower, T = s.symbols("b upper lower T", positive=True)
    Tc = (b**2 - upper**2) / (2 * upper)
    check("outer_support_shell_gap", b**2 + T**2 - (upper + T)**2 - 2 * upper * (Tc - T))
    check("inner_support_positive_s_gap", (lower - T)**2 - T**2 - lower * (lower - 2 * T))
    control("static_shell_contact_time", (b**2 + T**2 - (upper + T)**2).subs(T, b - upper))

    if args.inject_failure:
        check("deliberate_failure", s.Integer(1))
    report = {
        "status": "PASS",
        "scope": "Exact initial-time ADM kinetic algebra, angular tensor witness, and causal support geometry; no shell-contact evolution test",
        "positive_checks": len(checks),
        "wrong_formula_controls": len(controls),
        "checks": checks,
        "controls": controls,
        "tensor_projection_normalized_sphere_mean": str(ttmean),
        "tensor_witness_normalized_sphere_norm_squared": str(ttnorm),
        "scalar_stress_TT_projection_coefficient": "-4/3",
        "sympy_version": s.__version__,
    }
    payload = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()
