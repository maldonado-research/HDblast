#!/usr/bin/env python3
"""Independent exact-rational L-infinity trace certificate, no L2 inference.

The analytic characteristic and retained-background statements remain inputs.
No primary waveform or response implementation is imported.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json

BACKGROUND_SHA256 = "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8"


def compute(background):
    raw = background.read_bytes()
    if sha256(raw).hexdigest() != BACKGROUND_SHA256:
        raise ValueError("Actual-background source pin mismatch")
    source = json.loads(raw)
    bounds = {k: F(v) for k, v in source["collar"]["simple_bounds_exact"].items()}
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise ValueError(name)

    ell = F(1, 10**10)
    time, depth = 3 * ell, 8 * ell
    a, potential, beta, gamma = F(1000), F(10**6), F(1000), F(250)
    check("all_retained_collar_gates", all(source["collar"]["gates"].values()))
    check("positive_inverse_inertia_bound", bounds["mu_lower"] > 0 and
          1 / bounds["mu_lower"] <= a)
    check("actual_potential_bound", bounds["absolute_W_upper"] <= potential)
    check("actual_boundary_potential_bound", bounds["absolute_boundary_beta_upper"] <= beta)
    check("actual_physical_trace_shift", bounds["absolute_C_b_over_two_upper"] <= gamma)
    check("strict_causal_domain", 0 < time < depth <= F(source["collar"]["inward_depth_exact"]))
    bulk_source = potential * (1 + time)
    boundary_source = 16 / ell + 1
    theta = a * time**2 * (time * potential + beta) + potential * time**2 / 2
    check("strict_absorption_gate", 0 <= theta < 1)
    displacement = time**2 * (a * boundary_source + (a * time + F(1, 2)) * bulk_source) / (1 - theta)
    velocity = a * time * (boundary_source + time * bulk_source +
                           (time * potential + beta) * displacement)
    incoming = time * (bulk_source + potential * displacement)
    error = ell * (incoming + velocity + gamma * displacement)
    check("displacement_less_than_15e_minus6", displacement < F(15, 10**6))
    check("boundary_velocity_less_than_49000", velocity < 49000)
    check("incoming_characteristic_less_than_1e_minus3", incoming < F(1, 1000))
    check("scaled_uniform_trace_error_less_than_5e_minus6", error < F(1, 200000))
    check("positive_computed_error", error > 0)
    # Independent signed-forcing refinement, using integration by parts of
    # the known q_ell' forcing rather than its supremum.
    r = potential * time**2 / 2
    s = a * beta * time**2 / 2
    determinant = (1 - r) * (1 - s) - a * r * time
    check("sharp_bulk_absorption_gate", 0 <= r < 1)
    check("sharp_boundary_absorption_gate", 0 <= s < 1)
    check("sharp_positive_inverse_gate", determinant > 0)
    trial_upper = 1 + time
    sharp_displacement = trial_upper * (r * (1 - s) + a * time * (1 + r)) / determinant
    sharp_velocity = a * trial_upper / determinant
    sharp_master = trial_upper + sharp_displacement
    check("sharp_master_formula_identity", sharp_master ==
          trial_upper * (1 - s + a * time) / determinant)
    check("sharp_linear_inequality_one", (1 - r) * sharp_displacement -
          time * sharp_velocity == r * trial_upper)
    check("sharp_linear_inequality_two", -a * r * sharp_displacement +
          (1 - s) * sharp_velocity == a * trial_upper * (1 + r))
    sharp_error = ell * (potential * time * sharp_master +
                        (1 + gamma * time) * sharp_velocity)
    check("sharp_scaled_error_strict_bound", sharp_error < F(10000004, 10**14))
    check("sharp_bound_improves_generic_bound", 0 < sharp_error < error)
    return {
        "status": "PASS_INDEPENDENT_UNIFORM_CHARACTERISTIC_TRACE_BOUND",
        "background_sha256": BACKGROUND_SHA256,
        "finite_check_groups": len(checks), "checks": checks,
        "exact_inputs": {k: str(v) for k, v in {
            "ell": ell, "time": time, "initial_depth": depth,
            "inverse_mu_upper": a, "potential_abs_upper": potential,
            "boundary_beta_abs_upper": beta, "physical_gamma_abs_upper": gamma,
            "bulk_source_abs_upper": bulk_source, "boundary_source_abs_upper": boundary_source}.items()},
        "exact_outputs": {k: str(v) for k, v in {
            "absorption_theta": theta, "bulk_displacement_upper": displacement,
            "boundary_velocity_upper": velocity, "incoming_characteristic_upper": incoming,
            "scaled_uniform_trace_error_upper": error,
            "generic_strict_uniform_error_upper": F(1, 200000),
            "sharp_determinant": determinant,
            "sharp_displacement_upper": sharp_displacement,
            "sharp_boundary_velocity_upper": sharp_velocity,
            "sharp_master_upper": sharp_master,
            "sharp_scaled_uniform_trace_error_upper": sharp_error,
            "sharp_strict_uniform_error_upper": F(10000004, 10**14)}.items()},
        "scope": "Exact-rational consequences of an independent pointwise characteristic proof and the pinned conditional actual-background enclosure. No L2-to-pointwise inference, numerical waveform solver, finite-epsilon nonlinear response bound, or formal PDE certification."
    }


def main():
    parser = ArgumentParser()
    parser.add_argument("--background", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compute(args.background)
    encoded = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
