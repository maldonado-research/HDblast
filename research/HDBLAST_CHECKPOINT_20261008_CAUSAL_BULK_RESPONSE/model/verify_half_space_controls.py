#!/usr/bin/env python3
"""Manufactured symbolic controls; reads no HDBLAST input or observational data.

Exact SymPy/rational identities check the displayed algebra. mpmath integral
values are explicitly diagnostic, not interval certificates or physical runs.
Failure checks use explicit exceptions and remain active under python -O.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def run() -> dict:
    checks = []

    def zero(name, expression):
        result = sp.simplify(expression)
        require(result == 0, f"{name}: nonzero residual {result}")
        checks.append({"name": name, "kind": "EXACT_SYMBOLIC", "residual": "0"})

    p, c, s, M, m2, g, y = sp.symbols("p c s M m2 g y", positive=True)
    q, ph, qd, phd, source = sp.symbols("q ph qd phd source", real=True)
    om, Q = sp.symbols("omega Q", positive=True)

    # The outward normal is -dy, so the bulk variation contributes +Phi_y.
    boundary_potential = -c * ph**2 / 2 + g * q * ph
    zero("boundary_variation_sign", c * ph - g * q + sp.diff(boundary_potential, ph))
    zero("brane_force_from_same_action", sp.diff(boundary_potential, q) - g * ph)
    zero("retarded_boundary_elimination", (-s - c) * g * q / (c + s) + g * q)

    # Square completion after the exact half-line trace inequality.
    left = (c + M) * ph**2 - 2 * g * q * ph + m2 * q**2
    right = (c + M) * (ph - g * q / (c + M))**2 + (m2 - g**2 / (c + M)) * q**2
    zero("boundary_trace_square_completion", left - right)

    up = sp.sqrt(2 / sp.pi) * (p * sp.cos(p * y) + c * sp.sin(p * y)) / sp.sqrt(p**2 + c**2)
    zero("Robin_generalized_eigenfunction", sp.diff(up, y).subs(y, 0) - c * up.subs(y, 0))
    zero("Robin_mode_bulk_eigenvalue", sp.diff(up, y, 2) + p**2 * up)
    zero("Robin_boundary_spectral_weight", up.subs(y, 0)**2 - 2 * p**2 / (sp.pi * (p**2 + c**2)))
    integrand = p**2 / ((p**2 + c**2) * (p**2 + s**2))
    partial_fractions = (s**2 / (p**2 + s**2) - c**2 / (p**2 + c**2)) / (s**2 - c**2)
    zero("Stieltjes_partial_fractions", integrand - partial_fractions)
    zero("integrated_Stieltjes_identity", (s - c) / (s**2 - c**2) - 1 / (c + s))
    zero("coincident_s_equals_c_limit", sp.limit((s - c) / (s**2 - c**2), s, c) - 1 / (2 * c))
    zero("mass_squared_measure_Jacobian", up.subs(y, 0)**2 / (2 * p) - p / (sp.pi * (p**2 + c**2)))

    kernel_cut = 1 / (c - sp.I * p)
    re_kernel, im_kernel = sp.re(kernel_cut), sp.im(kernel_cut)
    zero("retarded_cut_real_part", re_kernel - c / (c**2 + p**2))
    zero("retarded_cut_imaginary_part", im_kernel - p / (c**2 + p**2))
    gamma = g**2 * im_kernel
    zero("inverse_absorption_affine_restriction", p / gamma - (c**2 + p**2) / g**2)
    zero("response_real_imaginary_restriction", g**2 * re_kernel / gamma - c / p)
    zero("absorption_peak_derivative", sp.diff(gamma, p).subs(p, c))
    zero("absorption_peak_value", gamma.subs(p, c) - g**2 / (2 * c))
    zero("spectral_continuum_response_identity", sp.im(1 / (m2 - p**2 - M**2 - g**2 * kernel_cut)) - gamma / sp.Abs(m2 - p**2 - M**2 - g**2 * kernel_cut)**2)

    # Bulk flux, Robin energy and interaction energy all matter.
    eb_dot = -phd * (c * ph - g * q)
    er_dot = c * ph * phd
    eq_dot = source * qd + g * ph * qd
    eint_dot = -g * (qd * ph + q * phd)
    zero("complete_signed_energy_ledger", eb_dot + er_dot + eq_dot + eint_dot - source * qd)
    amp = g * Q / (c - sp.I * p)
    flux = -sp.re((-sp.I * om * amp) * sp.conjugate(sp.I * p * amp)) / 2
    self_force_work = sp.re((g**2 * kernel_cut * Q) * sp.conjugate(-sp.I * om * Q)) / 2
    zero("outgoing_bulk_energy_flux", flux - om * gamma * Q**2 / 2)
    zero("flux_equals_minus_self_force_work", flux + self_force_work)

    # Every manufactured inequality below is exact rational arithmetic.
    def validate_stability(mass2, bulk_mass, robin, coupling2):
        if not (mass2 > 0 and bulk_mass > 0 and robin > 0 and 0 < coupling2 < mass2 * (robin + bulk_mass)):
            raise ValueError("strict_stable_domain_required")

    validate_stability(sp.Rational(9), sp.Rational(2), sp.Rational(3), sp.Rational(1))
    require(sp.Rational(9) - 4 - sp.Rational(1, 3) > 0, "no-bound manufactured example")
    validate_stability(sp.Rational(13, 4), sp.Rational(2), sp.Rational(3), sp.Rational(1))
    zb, sb = sp.Rational(3), sp.Rational(1)
    zero("manufactured_bound_mode_mass", sp.Rational(13, 4) - zb - 1 / (3 + sb))
    zero("manufactured_bound_mode_residue", 1 / (1 + 1 / (2 * sb * (3 + sb)**2)) - sp.Rational(32, 33))
    require(sp.Rational(32, 33) > 0, "bound residue positivity")
    # This exact unstable example has M=2,c=3,m²=4,g²=54,z=-5,s=3.
    zero("manufactured_tachyon_when_stability_fails", sp.Rational(4) - (-5) - sp.Rational(54, 6))
    require(sp.Rational(54) > 4 * (3 + 2), "unstable example lies outside contract")
    negatives = []
    for name, values in [
        ("reject_tachyonic_domain", (4, 2, 3, 54)),
        ("reject_marginal_massless_boundary", (sp.Rational(1, 5), 2, 3, 1)),
        ("reject_decoupled_branch_cut_claim", (9, 2, 3, 0)),
    ]:
        try:
            validate_stability(*map(sp.Rational, values))
        except ValueError:
            negatives.append({"name": name, "status": "EXPECTED_REJECTION"})
        else:
            raise RuntimeError(f"negative control was accepted: {name}")

    wrong_branch_im = sp.im(1 / (sp.Rational(3) + sp.I * sp.Rational(2)))
    require(wrong_branch_im < 0, "wrong branch must violate positive-frequency passivity")
    negatives.append({"name": "wrong_outgoing_branch_has_negative_spectral_density", "status": "EXPECTED_REJECTION", "imaginary_part": str(wrong_branch_im)})
    for pv in [sp.Rational(1, 10), sp.Rational(1), sp.Rational(3), sp.Rational(10)]:
        weight = pv / (sp.pi * (9 + pv**2))
        require(weight > 0, f"spectral weight at manufactured p={pv}")
    checks.append({"name": "manufactured_positive_spectral_weights", "kind": "EXACT_RATIONAL_SIGN_WITH_PI_POSITIVE", "status": "PASS"})

    # Finite positive spectral Gram example: diagnostic of the same positivity
    # mechanism, not a discretization of actual quantum/observational data.
    times = [0, sp.pi / 2, sp.pi, 3 * sp.pi / 2]
    frequencies = [1, 2, 3]
    weights = [sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(3, 4)]
    cov = sp.Matrix([[sum(w * sp.cos(f * (ti - tj)) for f, w in zip(frequencies, weights)) for tj in times] for ti in times])
    eigen = cov.eigenvals()
    require(all(ev >= 0 for ev in eigen), "manufactured covariance must be positive semidefinite")
    checks.append({"name": "manufactured_spectral_covariance_Gram_positivity", "kind": "EXACT_SYMBOLIC", "eigenvalues_with_multiplicity": {str(k): int(v) for k, v in eigen.items()}})

    # Independent quadrature checks of (11), kept explicitly non-certified.
    mp.mp.dps = 80
    diagnostics = []
    for real, imag in [(1, 0), (3, 0), (2, 1), (1, 4), (7, -2)]:
        sv = mp.mpc(real, imag)
        cv = mp.mpf(3)
        actual = 2 / mp.pi * mp.quad(lambda pv: pv**2 / ((pv**2 + cv**2) * (pv**2 + sv**2)), [0, 1, mp.inf])
        expected = 1 / (cv + sv)
        error = abs(actual - expected)
        require(error < mp.mpf("1e-65"), "manufactured spectral quadrature diagnostic")
        diagnostics.append({"c": "3", "s": [str(real), str(imag)], "absolute_discrepancy": mp.nstr(error, 8)})
    beta, omega = mp.mpf("1.5"), mp.mpf("4")
    coth = mp.coth(beta * omega / 2)
    require(coth > 1 and mp.coth(-beta * omega / 2) == -coth, "FDT sign and thermal factor")

    return {
        "schema_version": 1,
        "status": "PASS_MANUFACTURED_CONTROLS_ONLY",
        "scope": "Exact symbolic identities and rational sign checks plus explicitly non-certified numerical diagnostics for a separate flat half-space scalar model.",
        "checks": checks,
        "negative_controls": negatives,
        "quadrature_diagnostics": {"kind": "MPMATH_NON_CERTIFIED_DIAGNOSTIC", "decimal_precision": 80, "samples": diagnostics},
        "versions": {"sympy": sp.__version__, "mpmath": mp.__version__},
        "control_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "physical_source_callbacks": 0,
        "physical_target_evaluations": 0,
        "retained_array_reads": 0,
        "retained_array_decodes": 0,
        "observational_likelihood_evaluations": 0,
        "project_statuses": {"metric_calibration": "FAIL", "full_continuous_certificate": "UNRESOLVED", "higher_dimensional_cause": "NOT_ESTABLISHED", "external_novelty": "NOT_ASSESSED"},
        "limitations": ["Symbolic controls do not replace independent review of functional-analysis hypotheses.", "Numerical diagnostics are not interval enclosures.", "No real physical or observational input is opened.", "No finite raw equal-time bath noise variance is asserted.", "No physical source is derived from a prescribed pulse."]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "exact_checks": len(result["checks"]), "negative_controls": len(result["negative_controls"]), "output": str(args.output)}))


if __name__ == "__main__":
    main()
