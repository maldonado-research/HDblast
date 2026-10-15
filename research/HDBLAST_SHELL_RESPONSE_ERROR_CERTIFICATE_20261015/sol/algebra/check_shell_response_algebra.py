#!/usr/bin/env python3
"""Exact symbolic audit of the shell-response algebra.

Only the frozen PHYSICAL_SHELL_RESPONSE.md is used as a mathematical input.
No numerical background, numerical solution, or proof certification is claimed.
Run: python check_shell_response_algebra.py
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as s


HERE = Path(__file__).resolve().parent
POSITIVE: list[dict[str, str]] = []
NEGATIVE: list[dict[str, str]] = []


def normal(expr):
    """Canonical exact hyperbolic/rational simplification, without floats."""
    return s.factor(s.simplify(s.trigsimp(s.expand(expr))))


def check(name, expr):
    result = normal(expr)
    if result != 0:
        raise AssertionError(f"{name}: expected zero, got {result}")
    POSITIVE.append({"name": name, "residual": "0"})


def reject(name, wrong_minus_correct, witness):
    """Reject a proposed identity using an exact nonzero symbolic witness."""
    defect = normal(wrong_minus_correct)
    value = normal(defect.subs(witness))
    if defect == 0 or value == 0 or value.has(s.Float):
        raise AssertionError(f"{name}: negative control was not rejected")
    NEGATIVE.append(
        {"name": name, "defect": str(defect), "exact_witness_defect": str(value)}
    )


t, z = s.symbols("t z", real=True)
rho, rho_b, p, p_b, a, p_y = s.symbols(
    "rho rho_b p p_b a p_y", nonzero=True, real=True
)
C, C_z, C_b = s.symbols("C C_z C_b", real=True)
mu, mu_hat = s.symbols("mu mu_hat", positive=True)
beta_hat, W, W_hat = s.symbols("beta_hat W W_hat", real=True)
q, v, acc, r = s.symbols("q v acc r", real=True)
jp, jm = s.symbols("jp jm", real=True)
q_hat, acc_hat, r_hat = s.symbols("q_hat acc_hat r_hat", real=True)
w0, wt, wz, wtt, wzz, wtz = s.symbols(
    "w wt wz wtt wzz wtz", real=True
)
f_hat, b_hat = s.symbols("f_hat b_hat", real=True)

H = 3 * s.tanh(t)
D = 8 / s.cosh(t) ** 2 - 4
k = s.Rational(35, 4) / s.cosh(t) ** 2 - s.Rational(25, 4)
m = s.Rational(2, 3) * rho**2 * p**2
beta = mu * k + C_b / 2
gamma, gamma_hat, d = s.symbols("gamma gamma_hat d", real=True)

# Equation (2), dy=-rho dz: alpha_z/alpha=-C/2.
alpha_log_z = s.Rational(3, 2) * (-rho**2 * a) / rho + (-rho * p_y) / p
check("alpha_log_derivative", alpha_log_z + rho * (3 * a + 2 * p_y / p) / 2)

# Equation (7) multiplied by -rho^2 gives equation (10).
# P_yy=(V_zz+rho*a*V_z)/rho^2 and P_y=-V_z/rho.
minus_rho2_bulk = -(
    wzz + rho * a * wz + rho * (2 * a + 2 * p_y / p) * wz
    - wtt - H * wt - 8 / s.cosh(t) ** 2 * w0 + 4 * w0 - m * w0
)
target_bulk = wtt - wzz + H * wt - rho * (3 * a + 2 * p_y / p) * wz + (D + m) * w0
check("y_to_z_master_operator", minus_rho2_bulk - target_bulk)

# Full differential conjugation, keeping all independent derivative jets.
f = s.cosh(t) ** (-s.Rational(3, 2))
f_log_t = normal(s.diff(f, t) / f)
f_tt_over_f = normal(s.diff(f, t, 2) / f)
alpha_log = -C / 2
alpha_zz_over_alpha = -C_z / 2 + C**2 / 4
conjugated = (
    wtt - wzz + (2 * f_log_t + H) * wt
    + (-2 * alpha_log - C) * wz
    + (f_tt_over_f + H * f_log_t - alpha_zz_over_alpha - C * alpha_log + D + m) * w0
)
potential = k + C_z / 2 + C**2 / 4 + m
check("conjugated_bulk_operator", conjugated - (wtt - wzz + potential * w0))
check("time_drift_removed", 2 * f_log_t + H)
check("space_drift_removed", -2 * alpha_log - C)
check("time_potential_k", f_tt_over_f + H * f_log_t + D - k)
check("initial_displacement_factor", f.subs(t, 0) - 1)
check("initial_velocity_factor", s.diff(f, t).subs(t, 0))

# At the shell alpha_b=1. Here r denotes w_z(t,0), acc=q''.
R = r - C_b * q / 2
B = -r + mu * acc + beta * q
transformed_original_boundary = acc + (f_tt_over_f + H * f_log_t + D) * q - R / mu
check("original_to_conjugated_boundary", transformed_original_boundary - B / mu)
check("physical_boundary_identity_with_residual", R - mu * (acc + k * q) + B)
check("physical_boundary_identity_on_solution", (R - mu * (acc + k * q)).subs(r, mu * acc + beta * q))

# Characteristic equations and boundary ODE, equation (11a).
E = wtt - wzz + W * w0
wf = s.Function("w")(t, z)
jpf = s.diff(wf, t) + s.diff(wf, z)
jmf = s.diff(wf, t) - s.diff(wf, z)
Ef = s.diff(wf, t, 2) - s.diff(wf, z, 2) + W * wf
check("plus_characteristic_residual", s.diff(jpf, t) - s.diff(jpf, z) + W * wf - Ef)
check("minus_characteristic_residual", s.diff(jmf, t) + s.diff(jmf, z) + W * wf - Ef)
check("boundary_characteristic_ode_residual", (mu * acc + v + beta * q - jp).subs(jp, v + r) - B)
check("boundary_reflection_identity", (jm - 2 * v + jp).subs({jp: v + r, jm: v - r}))
check("time_reconstruction", (jp + jm).subs({jp: wt + wz, jm: wt - wz}) / 2 - wt)
check("space_reconstruction", (jp - jm).subs({jp: wt + wz, jm: wt - wz}) / 2 - wz)
check("physical_trace_from_incoming_characteristic", (jp - v - C_b * q / 2 - R).subs(jp, v + r))

# Candidate equations L_hat(w_hat)=f_hat, B_hat(q_hat)=b_hat.
# The residuals in the *true* equations include coefficient discrepancies.
candidate_acc = (b_hat + r_hat - beta_hat * q_hat) / mu_hat
G = b_hat + (mu - mu_hat) * acc_hat + (beta - beta_hat) * q_hat
true_B_hat = -r_hat + mu * acc_hat + beta * q_hat
check(
    "true_boundary_residual_includes_coefficients",
    (true_B_hat - G).subs(b_hat, -r_hat + mu_hat * acc_hat + beta_hat * q_hat),
)
check(
    "candidate_boundary_acceleration",
    (-r_hat + mu_hat * acc_hat + beta_hat * q_hat - b_hat).subs(acc_hat, candidate_acc),
)
F = f_hat + (W - W_hat) * w0
candidate_second_difference = f_hat - W_hat * w0
check("true_bulk_residual_includes_coefficients", candidate_second_difference + W * w0 - F)

# Error orientation is exact minus candidate. Then L(e)=-F and B(e)=-G.
# e_r is eliminated using B(e)=-G to check both resulting identities.
e_q, e_v, e_acc, e_r, e_jp = s.symbols("e_q e_v e_acc e_r e_jp", real=True)
e_r_from_boundary = mu * e_acc + beta * e_q + G
check(
    "error_boundary_ode_sign",
    (mu * e_acc + e_v + beta * e_q - e_jp + G).subs(e_jp, e_v + e_r).subs(e_r, e_r_from_boundary),
)
check(
    "error_physical_trace_residual_sign",
    (e_r - C_b * e_q / 2 - mu * (e_acc + k * e_q) - G).subs(e_r, e_r_from_boundary),
)

# A coherent normalized residual uses the candidate's own mu_hat, gamma_hat.
# Its raw positive true-boundary residual is G, and signed error forcing g=-G.
R_hat = r_hat - gamma_hat * q_hat
check("candidate_trace_coefficient_conversion", (r_hat - gamma * q_hat) - (R_hat - (gamma - gamma_hat) * q_hat))
normalized_d = acc_hat + k * q_hat - R_hat / mu_hat
G_normalized = mu * d + (mu - mu_hat) * R_hat / mu_hat + (gamma - gamma_hat) * q_hat
G_true_coherent = -r_hat + mu * acc_hat + (mu * k + gamma) * q_hat
check("normalized_boundary_residual_identity", G_normalized.subs(d, normalized_d) - G_true_coherent)
g_normalized = -mu * d - (mu - mu_hat) * R_hat / mu_hat - (gamma - gamma_hat) * q_hat
check("normalized_signed_error_forcing", g_normalized.subs(d, normalized_d) + G_true_coherent)

# An artificial W=beta=gamma=k=0 wave family tests a proposed norm control.
# This is not a substitute physical background or a source-model solution.
omega = s.symbols("omega", positive=True)
u = (s.cos(omega * z) / (mu * omega**2) - s.sin(omega * z) / omega) * s.cos(omega * t)
uq = u.subs(z, 0)
ur = s.diff(u, z).subs(z, 0)
check("artificial_high_frequency_bulk_residual", s.diff(u, t, 2) - s.diff(u, z, 2))
check("artificial_high_frequency_boundary_residual", -ur + mu * s.diff(uq, t, 2))
check("artificial_high_frequency_trace", ur + s.cos(omega * t))
amplitude_bound = 1 / (mu * omega**2) + 1 / omega
check("artificial_high_frequency_value_bound_tends_to_zero", s.limit(amplitude_bound, omega, s.oo))
check("artificial_high_frequency_trace_at_zero", ur.subs(t, 0) + 1)
n = s.symbols("n", positive=True, integer=True)
check("artificial_high_frequency_trace_l2_squared", s.integrate(ur.subs(omega, n)**2, (t, 0, 2 * s.pi)) - s.pi)

# Observable: equation (20), equation (9), alpha_b=1, and V_z=f R.
S = 3 * f * R / (rho_b**3 * p_b)
density_from_source = s.pi**2 * rho_b**4 * s.cosh(t) ** 3 * S**2 / 3
density = 3 * s.pi**2 * R**2 / (rho_b**2 * p_b**2)
check("weighted_observable_cosh_cancellation", density_from_source - density)

# Angular normalization uses exact S^3 second/fourth coordinate moments.
angular_norm = 2 * s.pi**2 * (s.Rational(1, 8) + s.Rational(1, 8) - 2 * s.Rational(1, 24))
check("quadrupole_norm", angular_norm - s.pi**2 / 3)

# The coefficient of eps^2 uses only first scalar variation and zeroth measure.
eps, phi1, phi2, measure0, measure1, measure2 = s.symbols("eps phi1 phi2 measure0 measure1 measure2")
expanded_integrand = (eps * phi1 + eps**2 * phi2 / 2) ** 2 * (measure0 + eps * measure1 + eps**2 * measure2 / 2)
check("quadratic_variation_independent_of_second_completion", s.expand(expanded_integrand).coeff(eps, 2) - phi1**2 * measure0)

# Deliberate wrong formulas must each have an exact nonzero witness.
wrong_f = s.cosh(t) ** s.Rational(3, 2)
reject("wrong_positive_cosh_conjugation", 2 * s.diff(wrong_f, t) / wrong_f + H, {t: s.log(2)})
wrong_S = 3 * s.cosh(t) ** s.Rational(3, 2) * R / (rho_b**3 * p_b)
wrong_density = s.pi**2 * rho_b**4 * s.cosh(t) ** 3 * wrong_S**2 / 3
reject("wrong_cosh_observable_factor", wrong_density - density, {t: s.acosh(2), r: 1, q: 0, rho_b: 1, p_b: 1})
reject("wrong_value_trace_w_instead_of_normal_trace", q - R, {q: 0, r: 1})
reject("wrong_uncorrected_wz_trace", r - R, {C_b: 2, q: 1})
reject("wrong_missing_boundary_coefficient_residual", G - b_hat, {mu: 1, mu_hat: 1, C_b: 0, beta_hat: 0, t: 0, q_hat: 1})
reject("wrong_missing_boundary_mass_residual", G - (b_hat + (beta - beta_hat) * q_hat), {mu: 2, mu_hat: 1, acc_hat: 1})
reject("wrong_acceleration_sign_in_trace", mu * (-acc + k * q) - mu * (acc + k * q), {mu: 1, acc: 1})
reject("wrong_acceleration_sign_in_boundary", (-r - mu * acc + beta * q) - B, {mu: 1, acc: 1})
reject("wrong_error_boundary_residual_sign", (e_r - C_b * e_q / 2 - mu * (e_acc + k * e_q) + G).subs(e_r, e_r_from_boundary), {b_hat: 1, q_hat: 0, acc_hat: 0})

report = {
    "scope": "Exact symbolic identities and deliberately false controls only; no numerical physical background or solution.",
    "sympy_version": s.__version__,
    "positive_checks_passed": len(POSITIVE),
    "negative_controls_rejected": len(NEGATIVE),
    "positive_checks": POSITIVE,
    "negative_controls": NEGATIVE,
}
out = HERE / "symbolic_check_results.json"
out.write_text(json.dumps(report, indent=2) + "\n")
print(f"PASS: {len(POSITIVE)} exact identities; {len(NEGATIVE)} deliberately wrong controls rejected.")
print(f"Results: {out}")
