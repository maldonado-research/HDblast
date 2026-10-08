#!/usr/bin/env python3
"""Exact manufactured controls for the additive incoming-scattering theorem.

No network, physical source, physical target, retained array, or observed data
is opened. SymPy identities and rational hypothesis tests remain active under
python -O. The packet decay proof is analytic; no sampled trajectory is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def run():
    exact = []
    negative = []

    def zero(name, expr):
        result = sp.simplify(expr)
        require(result == 0, f"{name}: {result}")
        exact.append({"name": name, "residual": "0", "kind": "EXACT_SYMBOLIC"})

    p, c, g, omega, sb = sp.symbols("p c g omega sb", positive=True)
    D = sp.symbols("D", real=True)
    F = (sp.I * p - c) * D + g**2
    R = ((sp.I * p + c) * D - g**2) / F
    Q = 2 * sp.I * p * g / F
    zero("boundary_equation", sp.I * p * (R - 1) - c * (1 + R) + g * Q)
    zero("brane_equation", D * Q - g * (1 + R))
    zero("boundary_trace", 1 + R - 2 * sp.I * p * D / F)
    zero("denominator_positive_sum_of_squares", F * sp.conjugate(F) - ((g**2 - c * D)**2 + p**2 * D**2))
    zero("reflection_phase_identity", R + sp.conjugate(F) / F)
    zero("unit_reflection", R * sp.conjugate(R) - 1)
    zero("D_zero_reflection_finite", R.subs(D, 0) + 1)
    zero("D_zero_brane_amplitude_finite", Q.subs(D, 0) - 2 * sp.I * p / g)
    zero("D_zero_trace", (1 + R).subs(D, 0))
    zero("D_zero_derivative_boundary_condition", (sp.I * p * (R - 1) + g * Q).subs(D, 0))

    R0 = (sp.I * p + c) / (sp.I * p - c)
    G = 1 / (c - sp.I * p)
    chi = 1 / (D - g**2 * G)
    zero("free_Robin_trace_phase", 1 + R0 - 2 * sp.I * p / (sp.I * p - c))
    zero("susceptibility_matches_incident_brane_motion", Q - chi * g * (1 + R0))
    zero("coherent_reflection_plus_radiation", R - R0 - g * G * Q)
    zero("decoupled_limit_away_from_D_zero_reflection", sp.limit(R, g, 0) - R0)
    zero("decoupled_limit_away_from_D_zero_brane", sp.limit(Q, g, 0))
    zero("incoming_plus_outgoing_flux", omega * p * (R * sp.conjugate(R) - 1) / 2)
    zero("stationary_work_into_brane", -omega * g * sp.im((1 + R) * sp.conjugate(Q)) / 2)

    # Radiation alone is nonzero. Its interference with the direct Robin
    # reflection is essential for the exact flux identity.
    rrad = g * G * Q
    rad_flux = omega * p * rrad * sp.conjugate(rrad) / 2
    interference_flux = omega * p * sp.re(R0 * sp.conjugate(rrad))
    zero("radiation_interference_cancellation", rad_flux + interference_flux)

    # Exact continuum/bound orthogonality uses the full bulk+brane Hilbert norm.
    Dbound = g**2 / (c + sb) - sb**2 - p**2
    pairing = Q + g / (c + sb) * (1 / (sb + sp.I * p) + R / (sb - sp.I * p))
    zero("continuum_bound_orthogonality", pairing.subs(D, Dbound))

    # Same signed energy ledger in real or complex Cartesian components.
    q, ph, qd, phd = sp.symbols("q ph qd phd", real=True)
    eb_dot = -phd * (c * ph - g * q)
    er_dot = c * ph * phd
    eq_dot = g * ph * qd
    eint_dot = -g * (qd * ph + q * phd)
    zero("closed_system_energy_ledger", eb_dot + er_dot + eq_dot + eint_dot)

    # Local nonstationary-phase integration-by-parts identity.
    t, massgap = sp.symbols("t massgap", positive=True)
    wp = sp.sqrt(massgap**2 + p**2)
    amp = sp.Function("f")(p)
    phase = sp.exp(-sp.I * t * wp)
    Lf = sp.diff(amp / sp.diff(wp, p), p)
    zero("packet_integration_by_parts_identity", sp.diff(amp * phase / sp.diff(wp, p), p) - Lf * phase + sp.I * t * amp * phase)
    zero("positive_group_velocity_formula", sp.diff(wp, p) - p / wp)
    zero("complex_packet_Plancherel_energy_weight", (wp**2 + p**2 + massgap**2) / 2 - wp**2)

    # Two exact manufactured points, one at the bare-brane frequency and one
    # with a coupled bound mode, use no physical source or data values.
    point = {p: sp.Rational(2), c: sp.Rational(3), g: sp.Rational(1), D: sp.Rational(1)}
    rv, qv = sp.simplify(R.subs(point)), sp.simplify(Q.subs(point))
    require(sp.simplify(rv * sp.conjugate(rv)) == 1, "manufactured reflected flux")
    exact.append({"name": "manufactured_scattering_values", "kind": "EXACT_RATIONAL_COMPLEX", "R": str(rv), "Q": str(qv)})
    zero("manufactured_bound_norm", (1 + g**2 / (2 * sb * (c + sb)**2)).subs({g: 1, c: 3, sb: 1}) - sp.Rational(33, 32))

    def domain(pv, cv, gv, m2v, Mv):
        if not (pv > 0 and cv > 0 and Mv > 0 and m2v > 0 and 0 < gv**2 < m2v * (cv + Mv)):
            raise ValueError("outside_strict_incident_contract")

    domain(sp.Rational(2), sp.Rational(3), sp.Rational(1), sp.Rational(9), sp.Rational(2))
    for name, values in [
        ("reject_threshold_packet", (0, 3, 1, 9, 2)),
        ("reject_decoupled_degenerate_formula", (sp.sqrt(5), 3, 0, 9, 2)),
        ("reject_unstable_scalar_action", (2, 3, sp.sqrt(54), 4, 2)),
    ]:
        try:
            domain(*values)
        except ValueError:
            negative.append({"name": name, "status": "EXPECTED_REJECTION"})
        else:
            raise RuntimeError(f"negative domain accepted: {name}")

    wrong_boundary_residual = sp.simplify((sp.I * p * (R - 1) - c * (1 + R) - g * Q).subs(point))
    require(wrong_boundary_residual != 0, "wrong coupling sign should fail matching")
    negative.append({"name": "wrong_boundary_coupling_sign", "status": "EXPECTED_REJECTION", "nonzero_residual": str(wrong_boundary_residual)})
    wrong_q = chi * g  # Missing the Robin incoming trace 1+R0.
    missing_phase_residual = sp.simplify((wrong_q - Q).subs(point))
    require(missing_phase_residual != 0, "missing incoming Robin normalization must fail")
    negative.append({"name": "missing_Robin_incoming_phase", "status": "EXPECTED_REJECTION", "nonzero_residual": str(missing_phase_residual)})
    rad_value = sp.simplify(rad_flux.subs(point).subs(omega, 3))
    require(rad_value > 0, "outgoing radiation-only flux must be positive")
    negative.append({"name": "omit_reflection_radiation_interference", "status": "EXPECTED_REJECTION", "spurious_unbalanced_flux": str(rad_value)})

    # An explicitly added bound component invalidates the local-decay premise.
    # M=2,c=3,g=1,m²=13/4 has z_b=3, norm_b²=33/32.
    bound_coefficient = sp.Rational(1, 7)
    persistent_energy = sp.Rational(3) * sp.Rational(33, 32) * bound_coefficient**2
    require(persistent_energy > 0, "bound contamination persists")
    negative.append({"name": "nonzero_initial_bound_projection", "status": "EXPECTED_REJECTION_OF_UNCONDITIONAL_DECAY", "persistent_complex_fiber_energy": str(persistent_energy)})
    # Endpoint terms cannot be dropped for a constant amplitude on [1,2].
    boundary_term_at_t0 = sp.sqrt(5) / 2 - sp.sqrt(2)
    require(boundary_term_at_t0 != 0, "constant compact-interval amplitude has nonzero IBP boundary term")
    negative.append({"name": "nonvanishing_packet_endpoint_amplitude", "status": "EXPECTED_REJECTION_OF_RAPID_DECAY_PROOF", "boundary_term_at_t0": str(boundary_term_at_t0)})

    return {
        "schema_version": 1,
        "status": "PASS_EXACT_MANUFACTURED_INCIDENT_CONTROLS",
        "scope": "Algebra and hypothesis controls for the separate autonomous flat scalar half-space model; no sampled physical trajectories.",
        "exact_controls": exact,
        "negative_controls": negative,
        "versions": {"sympy": sp.__version__},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "physical_source_callbacks": 0,
        "physical_target_evaluations": 0,
        "retained_array_reads": 0,
        "retained_array_decodes": 0,
        "observational_likelihood_evaluations": 0,
        "network_calls": 0,
        "project_statuses": {"metric_calibration": "FAIL", "full_continuous_certificate": "UNRESOLVED", "higher_dimensional_cause": "NOT_ESTABLISHED", "external_novelty": "NOT_ASSESSED"},
        "limitations": ["Analytic packet-decay proof requires smooth compact support away from p=0 and zero bound projection.", "Fixed-k energies are Fourier-fiber energies; finite global brane energy needs the stated smooth k packet.", "No physical origin, thermalization, or finite-data model selection is established."]
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "exact_controls": len(result["exact_controls"]), "negative_controls": len(result["negative_controls"]), "output": str(args.output)}))


if __name__ == "__main__":
    main()
