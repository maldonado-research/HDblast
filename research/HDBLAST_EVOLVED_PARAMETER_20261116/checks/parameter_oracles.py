#!/usr/bin/env python3
"""Bounded exact witnesses for the actual-parameter proof interfaces.

The source-bound model is the normalized 16-field scalar-principal system in
four spatial directions with all five fixed-connection first-jet slots.  The
symbolic coefficient witness below has TWO fields, nonlinear state-dependent
C,E,G,A,B and all FIVE jet slots.  It checks coefficient identities, not the
evaluation of the actual 16-field coefficient arrays or a continuum theorem.
The finite temporal-row and integrating-factor examples have the same limited
scope.  Analytic estimates and existence are established only by the proofs.
No ancestor producer is imported or executed.
"""
from __future__ import annotations

import argparse
import copy
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import re
import sys

import sympy as S
import mpmath


_RATIONAL_DEFAULTS = {
    "F_V_C_factor": "1", "F_V_E_factor": "1", "F_V_G_factor": "1",
    "F_Q_E_factor": "1", "F_Q_G_left_factor": "1", "F_Q_G_right_factor": "1",
    "F_VV_C_factor": "1", "F_VV_E_factor": "1", "F_VV_G_factor": "1",
    "F_VQ_E_factor": "1", "F_VQ_G_left_factor": "1", "F_VQ_G_right_factor": "1",
    "F_QV_E_factor": "1", "F_QV_G_left_factor": "1", "F_QV_G_right_factor": "1",
    "F_QQ_G_factor": "2",
    "A_V_first_factor": "1", "B_V_first_factor": "1",
    "A_VV_second_factor": "1", "B_VV_second_factor": "1",
    "A_V_cross_factor": "2", "B_V_cross_factor": "2",
    "linear_z_factor": "1", "quadratic_z_extra_factor": "0",
    "residual_A_difference_factor": "1", "residual_B_difference_factor": "1",
    "linearized_A_difference_factor": "1", "linearized_B_difference_factor": "1",
    "linearized_FV_difference_factor": "1", "linearized_FQ_difference_factor": "1",
    "linearized_AV_difference_factor": "1", "linearized_BV_difference_factor": "1",
    "linearized_AV_base_difference_factor": "1", "linearized_BV_base_difference_factor": "1",
    "linearized_segment_FVV_factor": "1", "linearized_segment_FVQ_factor": "1",
    "linearized_segment_FQV_factor": "1", "linearized_segment_FQQ_factor": "1",
    "inverse_first_factor": "-1", "inverse_cross_factor": "-1",
    "inverse_lapse_Hessian_factor": "-1", "inverse_square_factor": "2",
    "normal_tau_multiplier": "4", "normal_tau_offset": "2",
    "gronwall_exponent_factor": "-1", "gronwall_damping_factor": "-1",
}
_INTEGER_DEFAULTS = {
    "source_fields": 16, "spatial_directions": 4, "firstjet_directions": 5,
    "angular_generators": 6, "Gamma_generators": 4,
    "norm_bundle_copies": 2, "auxiliary_fields": 2,
    "initial_grade": 14, "temporal_last_recurrence": 12,
    "temporal_input_offset": 1, "temporal_output_offset": 2,
    "parameter_derivative_order": 2, "coefficient_derivative_reserve": 16,
    "initial_mean_value_slots": 1, "initial_Hessian_continuity_slots": 1,
    "dirichlet_last_trace": 13, "neumann_last_trace": 12,
    "energy_source_loss": 1, "normal_source_loss": 2,
}
_BOOLEAN_DEFAULTS = {
    "scalar_principal": True, "top_L2_boundary_trace": False,
    "initial_explicit_base_time": True,
    "independent_high_temporal_jets": False, "new_numeric_estimate_constant": False,
    "new_small_time_assumption": False,
}
_MODE_DEFAULTS = {
    "parameter_binding": "initial_only_fixed_background",
    "canonical_binding": "half_action_fixed_density",
    "jet_mode": "all_five",
    "quadratic_slot_mode": "combined_symmetric",
    "mixed_norm_mode": "all_time_all_rotations_normal",
    "initial_boundary_mode": "fixed_primal_chart",
    "pair_temporal_mode": "PDE_rows_from_the_pair",
    "residual_second_derivative_state": "a",
    "linearized_compared_field_state": "u_a",
    "initial_symmetry_support": "same_Gamma_same_A0",
}
_MODES = {
    "parameter_binding": {"initial_only_fixed_background", "vary_background"},
    "canonical_binding": {"half_action_fixed_density", "doubled_action"},
    "jet_mode": {"all_five", "spatial_four"},
    "quadratic_slot_mode": {"combined_symmetric", "direction_diagonal"},
    "mixed_norm_mode": {"all_time_all_rotations_normal", "spatial_only"},
    "initial_boundary_mode": {"fixed_primal_chart", "physical_N_without_chain_rule"},
    "pair_temporal_mode": {"PDE_rows_from_the_pair", "unrelated_temporal_jets"},
    "residual_second_derivative_state": {"a", "b"},
    "linearized_compared_field_state": {"u_a", "u_b"},
    "initial_symmetry_support": {"same_Gamma_same_A0", "arbitrary_support"},
}
_DEFAULTS = {**_RATIONAL_DEFAULTS, **_INTEGER_DEFAULTS, **_BOOLEAN_DEFAULTS, **_MODE_DEFAULTS}
FACTORS = tuple(_DEFAULTS)
MUTATIONS = tuple(
    [(name, "0" if value != "0" else "1") for name, value in _RATIONAL_DEFAULTS.items()]
    + [(name, value + 1) for name, value in _INTEGER_DEFAULTS.items()]
    + [(name, not value) for name, value in _BOOLEAN_DEFAULTS.items()]
    + [(name, next(value for value in sorted(_MODES[name]) if value != default))
       for name, default in _MODE_DEFAULTS.items()]
)
_RATIONAL_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _bounded_rational(value):
    """Parse reduced canonical rationals only; never sympify arbitrary input."""
    if type(value) is not str or len(value) > 24 or not _RATIONAL_RE.fullmatch(value):
        raise ValueError("expected a bounded canonical rational string")
    numerator, separator, denominator = value.partition("/")
    p, q = int(numerator), int(denominator) if separator else 1
    if abs(p) > 1000000 or q > 1000000:
        raise ValueError("rational exceeds the finite parser bound")
    result = S.Rational(p, q)
    if str(result) != value:
        raise ValueError("rational is not in reduced canonical form")
    return result


@lru_cache(maxsize=4096)
def _scalar_zero(expression):
    numerator = S.cancel(expression).as_numer_denom()[0]
    return S.Poly(S.expand(numerator), *sorted(numerator.free_symbols, key=str)).is_zero if numerator.free_symbols else numerator == 0


def zero(expression):
    if isinstance(expression, S.MatrixBase):
        return all(_scalar_zero(entry) for entry in expression)
    return _scalar_zero(S.sympify(expression))


def _sum_vectors(vectors, n=2):
    return sum(vectors, S.zeros(n, 1))


@lru_cache(maxsize=1)
def _model():
    """An independent nonlinear two-field/five-jet symbolic witness."""
    n, directions = 2, 5
    v = S.Matrix(S.symbols("v0:2"))
    w = S.Matrix(S.symbols("w0:2"))
    z = S.Matrix(S.symbols("z0:2"))
    Q = S.Matrix(directions, n, lambda mu, a: S.Symbol("Q_%d_%d" % (mu, a)))
    p = S.Matrix(directions, n, lambda mu, a: S.Symbol("p_%d_%d" % (mu, a)))
    r = S.Matrix(directions, n, lambda mu, a: S.Symbol("r_%d_%d" % (mu, a)))

    def dv(expression, direction):
        return expression.applyfunc(lambda e: sum(S.diff(e, v[a]) * direction[a] for a in range(n))) if isinstance(expression, S.MatrixBase) else sum(S.diff(expression, v[a]) * direction[a] for a in range(n))

    def dq(expression, direction):
        return expression.applyfunc(lambda e: sum(S.diff(e, Q[mu, a]) * direction[mu, a] for mu in range(directions) for a in range(n)))

    def projection(q, output, kind):
        return sum(S.Rational((mu + 1) * (a + 2) + output + kind + 1, 11 + output + 2 * kind + mu + a) * q[mu, a]
                   for mu in range(directions) for a in range(n))

    C = S.Matrix([(a + 1) * v[0] + (a + 2) * v[1] + (a + 1) * v[0] ** 2 * v[1] + (a + 2) * v[0] * v[1] ** 2 + v[0] ** 3 / (a + 2) for a in range(n)])

    def E(q):
        return S.Matrix([(1 + v[0] + (a + 1) * v[1] ** 2) * projection(q, a, 0)
                         + (v[0] * v[1] + (a + 2) * v[0] ** 2) * projection(q, a, 1) for a in range(n)])

    def G(q1, q2):
        return S.Matrix([(1 + v[0] ** 2 + (a + 1) * v[1]) * projection(q1, a, 0) * projection(q2, a, 0)
                         + (v[0] * v[1] + (a + 1) * v[1] ** 2) *
                         (projection(q1, a, 0) * projection(q2, a, 1) + projection(q1, a, 1) * projection(q2, a, 0)) / 2
                         + (1 + (a + 1) * v[0] + v[1] ** 2) * projection(q1, a, 1) * projection(q2, a, 1) for a in range(n)])

    def G_diagonal(q1, q2):
        # Meaningful wrong model: retain only equal base-direction pairs.
        out = S.zeros(n, 1)
        for mu in range(directions):
            q1mu, q2mu = S.zeros(directions, n), S.zeros(directions, n)
            for a in range(n):
                q1mu[mu, a], q2mu[mu, a] = q1[mu, a], q2[mu, a]
            out += G(q1mu, q2mu)
        return out

    F = C + E(Q) + G(Q, Q)
    Fv, Fq = dv(F, w), dq(F, p)
    expected = {
        "F_V": Fv, "F_Q": Fq, "F_VV": dv(Fv, w),
        "F_VQ": dq(Fv, p), "F_QV": dv(Fq, w), "F_QQ": dq(Fq, p),
    }

    A = [S.Rational(i + 1, 7) + v[0] ** 2 / (i + 2) + v[0] * v[1] / (i + 3) + v[1] ** 3 / (i + 4) for i in range(4)]
    B = S.Matrix(4, 4, lambda i, j: S.Rational(i + j + 2, 11) + v[0] * v[1] / (i + j + 2)
                 + v[0] ** 2 / (i + j + 3) + v[1] ** 2 / (i + j + 5) + v[0] * v[1] ** 2 / (i + j + 7))
    mixed = {key: [S.Matrix(S.symbols("%s_%d_0:2" % (key, i))) for i in range(4)] for key in ("Vts", "Wts", "Zts")}
    hessian = {key: [[None] * 4 for _ in range(4)] for key in ("VH", "WH", "ZH")}
    for key, table in hessian.items():
        for i in range(4):
            for j in range(i, 4):
                table[i][j] = table[j][i] = S.Matrix(S.symbols("%s_%d_%d_0:2" % (key, i, j)))

    def Aprincipal(coefs, jets):
        return _sum_vectors(coefs[i] * jets[i] for i in range(4))

    def Bprincipal(coefs, jets):
        return _sum_vectors(coefs[i, j] * jets[i][j] for i in range(4) for j in range(4))

    Av, Bv = [dv(a, w) for a in A], B.applyfunc(lambda b: dv(b, w))
    Avv, Bvv = [dv(a, w) for a in Av], Bv.applyfunc(lambda b: dv(b, w))
    Avz, Bvz = [dv(a, z) for a in A], B.applyfunc(lambda b: dv(b, z))
    principal = {
        "Av_base": Aprincipal(Av, mixed["Vts"]), "Bv_base": Bprincipal(Bv, hessian["VH"]),
        "Avv_base": Aprincipal(Avv, mixed["Vts"]), "Bvv_base": Bprincipal(Bvv, hessian["VH"]),
        "Av_first": Aprincipal(Av, mixed["Wts"]), "Bv_first": Bprincipal(Bv, hessian["WH"]),
        "Avz_base": Aprincipal(Avz, mixed["Vts"]), "Bvz_base": Bprincipal(Bvz, hessian["VH"]),
    }

    # Direct two-jet substitution, independent of the hand-assembled Frechet
    # formulas.  Twice the h^2 coefficient is its second h derivative at zero.
    h = S.Symbol("formal_h")
    substitutions = {v[a]: v[a] + h * w[a] + h ** 2 * z[a] / 2 for a in range(n)}
    substitutions.update({Q[mu, a]: Q[mu, a] + h * p[mu, a] + h ** 2 * r[mu, a] / 2 for mu in range(directions) for a in range(n)})
    Fh = F.subs(substitutions, simultaneous=True)
    Ah = [a.subs(substitutions, simultaneous=True) for a in A]
    Bh = B.subs(substitutions, simultaneous=True)
    Ah_product = Aprincipal(Ah, [mixed["Vts"][i] + h * mixed["Wts"][i] + h ** 2 * mixed["Zts"][i] / 2 for i in range(4)])
    Bh_product = Bprincipal(Bh, [[hessian["VH"][i][j] + h * hessian["WH"][i][j] + h ** 2 * hessian["ZH"][i][j] / 2 for j in range(4)] for i in range(4)])
    direct_first = (Fh + Ah_product + Bh_product).diff(h).subs(h, 0) - Aprincipal(A, mixed["Wts"]) - Bprincipal(B, hessian["WH"])
    direct_second = (Fh + Ah_product + Bh_product).diff(h, 2).subs(h, 0) - Aprincipal(A, mixed["Zts"]) - Bprincipal(B, hessian["ZH"])
    linear_z = dv(F, z) + dq(F, r) + principal["Avz_base"] + principal["Bvz_base"]
    mixed_time_angular = F.diff(Q[0, 0], Q[2, 1])

    return {"v": v, "w": w, "z": z, "Q": Q, "p": p, "r": r, "C": C,
            "E": E, "G": G, "G_diagonal": G_diagonal, "F": F, "dv": dv, "dq": dq,
            "expected": expected, "A": A, "B": B, "principal": principal,
            "direct_first": direct_first, "direct_second": direct_second,
            "linear_z": linear_z, "mixed_time_angular": mixed_time_angular}


def _source_formulas(model, f, claims):
    dv, C, E = model["dv"], model["C"], model["E"]
    G = model["G"] if claims["quadratic_slot_mode"] == "combined_symmetric" else model["G_diagonal"]
    Q, p, w = model["Q"], model["p"].copy(), model["w"]
    if claims["jet_mode"] == "spatial_four":
        for a in range(2):
            p[0, a] = 0
    return {
        "F_V": f("F_V_C_factor") * dv(C, w) + f("F_V_E_factor") * dv(E(Q), w) + f("F_V_G_factor") * dv(G(Q, Q), w),
        "F_Q": f("F_Q_E_factor") * E(p) + f("F_Q_G_left_factor") * G(p, Q) + f("F_Q_G_right_factor") * G(Q, p),
        "F_VV": f("F_VV_C_factor") * dv(dv(C, w), w) + f("F_VV_E_factor") * dv(dv(E(Q), w), w) + f("F_VV_G_factor") * dv(dv(G(Q, Q), w), w),
        "F_VQ": f("F_VQ_E_factor") * dv(E(p), w) + f("F_VQ_G_left_factor") * dv(G(p, Q), w) + f("F_VQ_G_right_factor") * dv(G(Q, p), w),
        "F_QV": f("F_QV_E_factor") * dv(E(p), w) + f("F_QV_G_left_factor") * dv(G(p, Q), w) + f("F_QV_G_right_factor") * dv(G(Q, p), w),
        "F_QQ": f("F_QQ_G_factor") * G(p, p),
    }


@lru_cache(maxsize=1)
def _residual_model():
    """Residual subtraction for two independent states, not assumed solutions."""
    a, b, Qa, Qb, Xa, Xb, Ha, Hb, TTa, TTb = S.symbols("a b Qa Qb Xa Xb Ha Hb TTa TTb")
    A = lambda u: 1 + u + u ** 2
    B = lambda u: 2 + u ** 2 + u ** 3
    F = lambda u, q: 1 + u ** 3 + (1 + u ** 2) * q + (1 + u) * q ** 2
    Ra = TTa - A(a) * Xa - B(a) * Ha - F(a, Qa)
    Rb = TTb - A(b) * Xb - B(b) * Hb - F(b, Qb)
    Pb_difference = TTb - TTa - A(b) * (Xb - Xa) - B(b) * (Hb - Ha)
    theta = S.Symbol("theta")
    u, q = S.symbols("interpolated_state interpolated_jet")
    interp = {u: a + theta * (b - a), q: Qa + theta * (Qb - Qa)}
    mean_F = S.integrate((S.diff(F(u, q), u) * (b - a) + S.diff(F(u, q), q) * (Qb - Qa)).subs(interp), (theta, 0, 1))
    mean_A = S.integrate(S.diff(A(u), u).subs(u, a + theta * (b - a)) * (b - a), (theta, 0, 1))
    mean_B = S.integrate(S.diff(B(u), u).subs(u, a + theta * (b - a)) * (b - a), (theta, 0, 1))
    return locals()


@lru_cache(maxsize=1)
def _inverse_model():
    # The two independent state directions test the bilinear Hessian, including
    # both numerator/inverse cross terms.  It applies to k=2b and k=c.
    h, l = S.symbols("inverse_h inverse_l")
    a0 = S.Symbol("positive_lapse_denominator", positive=True)
    ah, al, ahl, k0, kh, kl, khl = S.symbols("a_h a_l a_hl numerator numerator_h numerator_l numerator_hl")
    a = a0 + h * ah + l * al + h * l * ahl
    numerator = k0 + h * kh + l * kl + h * l * khl
    direct_first = S.diff(numerator / a, h).subs({h: 0, l: 0})
    direct_second = S.diff(numerator / a, h, l).subs({h: 0, l: 0})
    return locals()


@lru_cache(maxsize=2)
def _nonautonomous_initial_rows(include_base_time):
    # A finite nonautonomous scalar Cauchy row with nonlinear velocity source.
    # This verifies explicit base-time differentiation in a small model only.
    s, e = S.symbols("initial_time initial_parameter")
    initial = [S.Rational(1, 7) + e / 11 + e ** 2 / 13,
               S.Rational(2, 9) - e / 17 + e ** 2 / 19]
    for order in range(3):
        v = sum(initial[j] * s ** j / S.factorial(j) for j in range(order + 2))
        vs = S.diff(v, s)
        base = s if include_base_time else 0
        rhs = base * v + v ** 2 + (1 + base * v) * vs + (v + base) * vs ** 2
        initial.append(S.expand(S.diff(rhs, s, order).subs(s, 0)))
    return tuple(tuple(S.diff(jet, e, k).subs(e, 0) for jet in initial) for k in range(3))


def evaluate(claims, check):
    """Evaluate fresh claims; callback signature is check(name, passed)."""
    records = []

    def need(name, condition):
        passed = bool(condition)
        full_name = "parameter_oracles/" + name
        check(full_name, passed)
        records.append({"name": name, "status": "pass" if passed else "fail"})
        if not passed:
            raise ValueError("parameter oracle rejects " + full_name)

    need("schema", type(claims) is dict and claims.get("schema") == "hdblast-parameter-oracle-claims-v1" and set(claims) == {"schema", *FACTORS})
    need("exact_engine", S.__version__ == "1.14.0" and mpmath.__version__ == "1.3.0")
    parsed = {}
    for name in _RATIONAL_DEFAULTS:
        try:
            parsed[name] = _bounded_rational(claims[name])
        except (ValueError, TypeError):
            need("parser/" + name, False)
    for name in _INTEGER_DEFAULTS:
        need("type/" + name, type(claims[name]) is int and 0 <= claims[name] <= 64)
    for name in _BOOLEAN_DEFAULTS:
        need("type/" + name, type(claims[name]) is bool)
    for name in _MODE_DEFAULTS:
        need("type/" + name, type(claims[name]) is str and claims[name] in _MODES[name])
    f = parsed.__getitem__

    # These are source/interface guards, not auxiliary 16-array evaluations.
    for name, expected in (("source_fields", 16), ("spatial_directions", 4), ("firstjet_directions", 5),
                           ("angular_generators", 6), ("Gamma_generators", 4),
                           ("norm_bundle_copies", 2), ("auxiliary_fields", 2)):
        need("source_interface/" + name, claims[name] == expected)
    need("source_interface/scalar_principal", claims["scalar_principal"])
    need("source_interface/fixed_external_parameter", claims["parameter_binding"] == "initial_only_fixed_background")
    need("source_interface/current_canonical_binding", claims["canonical_binding"] == "half_action_fixed_density")
    need("source_interface/two_copies_full_mixed_norm", claims["mixed_norm_mode"] == "all_time_all_rotations_normal")

    model = _model()
    produced = _source_formulas(model, f, claims)
    names = {"F_V": "F_V_complete_state_derivative", "F_Q": "F_Q_complete_five_slot_derivative",
             "F_VV": "F_VV_state_Hessian", "F_VQ": "F_VQ_state_jet_derivative",
             "F_QV": "F_QV_jet_state_derivative", "F_QQ": "F_QQ_twice_symmetric_quadratic"}
    need("auxiliary/five_by_two_matrix_Q", model["Q"].shape == (5, 2))
    need("auxiliary/nonzero_time_angular_quadratic_coupling", not zero(model["mixed_time_angular"]))
    need("auxiliary/nonlinear_state_coefficient_Hessians", not zero(model["expected"]["F_VV"]))
    need("auxiliary/C_vanishes_at_zero_state", zero(model["C"].subs({entry: 0 for entry in model["v"]})))
    need("auxiliary/combined_slot_G_symmetry", zero(model["G"](model["Q"], model["p"]) - model["G"](model["p"], model["Q"])))
    for name, diagnostic in names.items():
        need(diagnostic, zero(produced[name] - model["expected"][name]))
    need("F_VQ_equals_F_QV", zero(produced["F_VQ"] - produced["F_QV"]))

    principal = model["principal"]
    first_principal = f("A_V_first_factor") * principal["Av_base"] + f("B_V_first_factor") * principal["Bv_base"]
    need("positive_A_V_first_principal", zero(f("A_V_first_factor") * principal["Av_base"] - principal["Av_base"]))
    need("positive_B_V_first_principal", zero(f("B_V_first_factor") * principal["Bv_base"] - principal["Bv_base"]))
    need("first_formal_h_coefficient", zero(produced["F_V"] + produced["F_Q"] + first_principal - model["direct_first"]))
    need("positive_A_VV_second_principal", zero((f("A_VV_second_factor") - 1) * principal["Avv_base"]))
    need("positive_B_VV_second_principal", zero((f("B_VV_second_factor") - 1) * principal["Bvv_base"]))
    need("factor_two_A_V_Dw_s", zero((f("A_V_cross_factor") - 2) * principal["Av_first"]))
    need("factor_two_B_V_Hsp_w", zero((f("B_V_cross_factor") - 2) * principal["Bv_first"]))
    need("same_linearized_operator_on_z", zero((f("linear_z_factor") - 1) * model["linear_z"]))
    need("no_linear_z_in_quadratic_source", zero(f("quadratic_z_extra_factor") * model["linear_z"]))
    quadratic = (produced["F_VV"] + produced["F_VQ"] + produced["F_QV"] + produced["F_QQ"]
                 + f("A_VV_second_factor") * principal["Avv_base"] + f("B_VV_second_factor") * principal["Bvv_base"]
                 + f("A_V_cross_factor") * principal["Av_first"] + f("B_V_cross_factor") * principal["Bv_first"]
                 + f("quadratic_z_extra_factor") * model["linear_z"])
    need("second_formal_twice_h2_coefficient", zero(f("linear_z_factor") * model["linear_z"] + quadratic - model["direct_second"]))

    inverse = _inverse_model()
    a0 = inverse["a0"]
    inverse_first = inverse["kh"] / a0 + f("inverse_first_factor") * inverse["ah"] * inverse["k0"] / a0 ** 2
    inverse_second = (inverse["khl"] / a0
                      + f("inverse_cross_factor") * (inverse["ah"] * inverse["kl"] + inverse["al"] * inverse["kh"]) / a0 ** 2
                      + f("inverse_lapse_Hessian_factor") * inverse["ahl"] * inverse["k0"] / a0 ** 2
                      + f("inverse_square_factor") * inverse["ah"] * inverse["al"] * inverse["k0"] / a0 ** 3)
    need("normalization/inverse_first_state_derivative", zero(inverse_first - inverse["direct_first"]))
    need("normalization/inverse_full_bilinear_state_Hessian", zero(inverse_second - inverse["direct_second"]))

    # Finite dependency graph for R_r: product H^(14-j), 0<=j<=r+1,
    # to H^(12-r).  Differentiating twice retains that same temporal graph.
    grade = claims["initial_grade"]
    last = claims["temporal_last_recurrence"]
    need("initial/retained_Y14_grade", grade == 14)
    need("initial/recurrence_through_J14", last + claims["temporal_output_offset"] == 14 and last == 12)
    need("initial/output_is_r_plus_two", claims["temporal_output_offset"] == 2)
    need("initial/input_ends_at_r_plus_one", claims["temporal_input_offset"] == 1)
    dependency_table = []
    for order in range(last + 1):
        inputs = list(range(order + claims["temporal_input_offset"] + 1))
        output = order + claims["temporal_output_offset"]
        target_sobolev = grade - output
        for spatial_order in range(13 - order):
            need("initial/intrinsic_budget_r%d_k%d" % (order, spatial_order), order + spatial_order + 2 <= 14)
        dependency_table.append({"r": order, "input_temporal_jets": inputs,
                                 "input_sobolev_orders": [grade - j for j in inputs],
                                 "output_jet": output, "output_sobolev_order": target_sobolev,
                                 "parameter_orders": [0, 1, 2]})
    need("initial/C2_parameter_hierarchy", claims["parameter_derivative_order"] == 2)
    need("initial/mean_value_extraction_slot", claims["initial_mean_value_slots"] == 1)
    need("initial/Hessian_continuity_slot", claims["initial_Hessian_continuity_slots"] == 1)
    need("initial/total_base_and_parameter_reserve", claims["coefficient_derivative_reserve"] == last + claims["parameter_derivative_order"] + claims["initial_mean_value_slots"] + claims["initial_Hessian_continuity_slots"] == 16)
    need("initial/no_independent_high_temporal_jets", claims["independent_high_temporal_jets"] is False)
    need("initial/fixed_chart_boundary_derivatives", claims["initial_boundary_mode"] == "fixed_primal_chart")
    need("initial/Dirichlet_trace_ceiling", claims["dirichlet_last_trace"] == grade - 1 == 13)
    need("initial/Neumann_trace_ceiling", claims["neumann_last_trace"] == grade - 2 == 12)
    need("initial/no_trace_of_top_L2_jet", claims["top_L2_boundary_trace"] is False)
    need("initial/same_symmetry_and_closed_support", claims["initial_symmetry_support"] == "same_Gamma_same_A0")
    need("initial/explicit_base_time_in_temporal_recursion", _nonautonomous_initial_rows(claims["initial_explicit_base_time"]) == _nonautonomous_initial_rows(True))

    # A finite, coupled constant-coefficient PDE Fourier row has time jets
    # J_(j+2)=M J_j+N J_(j+1).  This proves only the displayed dependency,
    # not compactness or a nonlinear estimate for the actual system.
    eta, free_high = S.symbols("pair_delta unrelated_high_jet")
    M = S.Matrix([[S.Rational(-7,5), S.Rational(2,7)], [S.Rational(3,11), S.Rational(-5,3)]])
    N = S.Matrix([[S.Rational(1,13), S.Rational(-2,17)], [S.Rational(4,19), S.Rational(3,23)]])
    Ja = [S.Matrix([S.Rational(1,3), S.Rational(2,5)]), S.Matrix([S.Rational(-3,7), S.Rational(4,11)])]
    Jb = [Ja[0] + eta * S.Matrix([S.Rational(1,7), S.Rational(2,13)]), Ja[1] + eta * S.Matrix([S.Rational(-1,17), S.Rational(3,19)])]
    for j in range(13):
        Ja.append(M * Ja[j] + N * Ja[j + 1])
        Jb.append(M * Jb[j] + N * Jb[j + 1])
    if claims["pair_temporal_mode"] == "unrelated_temporal_jets":
        Jb[2] += S.Matrix([free_high, 0])
    need("pair_compactness/actual_PDE_temporal_row_dependency", all(zero(Jb[j + 2] - M * Jb[j] - N * Jb[j + 1]) for j in range(13)))
    need("pair_compactness/pair_difference_controls_all_finite_rows", all(zero((Jb[j] - Ja[j]).subs(eta, 0)) for j in range(15)))
    need("pair_compactness/unrelated_jets_are_a_real_obstruction", not zero(S.Matrix([free_high, 0])))

    need("energy/full_time_rotation_source_order", claims["energy_source_loss"] == 1 and all(r - claims["energy_source_loss"] == expected for r, expected in ((13, 12), (12, 11))))
    need("energy/full_normal_time_rotation_source_order", claims["normal_source_loss"] == 2 and all(r - claims["normal_source_loss"] == expected for r, expected in ((13, 11), (12, 10))))
    residual = _residual_model()
    second_X = residual["Xa"] if claims["residual_second_derivative_state"] == "a" else residual["Xb"]
    second_H = residual["Ha"] if claims["residual_second_derivative_state"] == "a" else residual["Hb"]
    produced_difference = (residual["Pb_difference"] - (residual["F"](residual["b"], residual["Qb"]) - residual["F"](residual["a"], residual["Qa"]))
                           - f("residual_A_difference_factor") * (residual["A"](residual["b"]) - residual["A"](residual["a"])) * second_X
                           - f("residual_B_difference_factor") * (residual["B"](residual["b"]) - residual["B"](residual["a"])) * second_H)
    need("residual/P_b_difference_with_second_derivatives_of_a", zero(produced_difference - residual["Rb"] + residual["Ra"]))
    need("residual/joint_state_jet_mean_value_F", zero(residual["mean_F"] - residual["F"](residual["b"], residual["Qb"]) + residual["F"](residual["a"], residual["Qa"])))
    need("residual/mean_value_A_and_B", zero(residual["mean_A"] - residual["A"](residual["b"]) + residual["A"](residual["a"])) and zero(residual["mean_B"] - residual["B"](residual["b"]) + residual["B"](residual["a"])))

    # Exact difference of the linearized operators on the SECOND compared
    # field u=w_a.  The derivatives of d=v_b-v_a are separate terms.
    a, b, Qa, Qb, Xa, Xb, Ha, Hb = (residual[key] for key in ("a", "b", "Qa", "Qb", "Xa", "Xb", "Ha", "Hb"))
    u, Du, uX, uH, ubX, ubH, uTT = S.symbols("compared_u compared_Du compared_uX compared_uH other_uX other_uH compared_uTT")
    state, jet = S.symbols("lin_state lin_jet")
    AA, BB, FF = residual["A"], residual["B"], residual["F"]
    AV = S.diff(AA(state), state)
    BV = S.diff(BB(state), state)
    FV, FQ = S.diff(FF(state, jet), state), S.diff(FF(state, jet), jet)
    at = lambda expression, v, q: expression.subs({state: v, jet: q}, simultaneous=True)

    def linearized(v, q, base_X, base_H):
        return (uTT - AA(v) * uX - BB(v) * uH - at(FV, v, q) * u - at(FQ, v, q) * Du
                - AV.subs(state, v) * u * base_X - BV.subs(state, v) * u * base_H)

    selected_uX, selected_uH = (uX, uH) if claims["linearized_compared_field_state"] == "u_a" else (ubX, ubH)
    linearized_source = (
        f("linearized_A_difference_factor") * (AA(b) - AA(a)) * selected_uX
        + f("linearized_B_difference_factor") * (BB(b) - BB(a)) * selected_uH
        + f("linearized_FV_difference_factor") * (at(FV, b, Qb) - at(FV, a, Qa)) * u
        + f("linearized_FQ_difference_factor") * (at(FQ, b, Qb) - at(FQ, a, Qa)) * Du
        + f("linearized_AV_difference_factor") * (AV.subs(state, b) - AV.subs(state, a)) * u * Xa
        + f("linearized_BV_difference_factor") * (BV.subs(state, b) - BV.subs(state, a)) * u * Ha
        + f("linearized_AV_base_difference_factor") * AV.subs(state, b) * u * (Xb - Xa)
        + f("linearized_BV_base_difference_factor") * BV.subs(state, b) * u * (Hb - Ha))
    need("residual/linearized_second_field_and_base_difference_terms", zero(linearized_source - linearized(a, Qa, Xa, Ha) + linearized(b, Qb, Xb, Hb)))
    theta = residual["theta"]
    interpolate = {state: a + theta * (b - a), jet: Qa + theta * (Qb - Qa)}
    integrand = (
        f("linearized_segment_FVV_factor") * S.diff(FV, state) * (b - a) * u
        + f("linearized_segment_FVQ_factor") * S.diff(FV, jet) * u * (Qb - Qa)
        + f("linearized_segment_FQV_factor") * S.diff(FQ, state) * Du * (b - a)
        + f("linearized_segment_FQQ_factor") * S.diff(FQ, jet) * (Qb - Qa) * Du)
    integrated = S.integrate(integrand.subs(interpolate, simultaneous=True), (theta, 0, 1))
    need("residual/linearized_full_joint_segment_Hessian", zero(integrated - (at(FV, b, Qb) - at(FV, a, Qa)) * u - (at(FQ, b, Qb) - at(FQ, a, Qa)) * Du))

    # The finite normal-count weighting is independent of time-slab length.
    c = S.Symbol("qualitative_normal_recovery_constant", nonnegative=True)
    tau = 1 / (f("normal_tau_multiplier") * c + f("normal_tau_offset"))
    leading_weight = S.cancel(c * tau / (1 - tau))
    need("normal_recovery/weighted_leading_coefficient", zero(leading_weight - c / (4 * c + 1)))
    need("normal_recovery/strict_quarter_absorption", zero(S.Rational(1, 4) - leading_weight - 1 / (4 * (4 * c + 1))))
    weight = S.Symbol("normal_weight", positive=True)
    need("normal_recovery/finite_triangular_geometric_weight", zero(sum(weight ** j for j in range(1, 13)) - weight * (1 - weight ** 12) / (1 - weight)))

    # Exact integrating-factor residual identity: K is an AUXILIARY rational,
    # not an asserted PDE constant.  No K*T smallness is used.
    s = S.Symbol("s", real=True)
    K = S.Rational(3, 5)
    energy = 3 + 2 * s + s ** 2
    integrating_factor = S.exp(f("gronwall_exponent_factor") * K * s)
    chosen_residual = S.diff(energy, s) + f("gronwall_damping_factor") * K * energy
    need("gronwall/exact_integrating_factor_residual", S.simplify(S.diff(integrating_factor * energy, s) - S.exp(-K * s) * chosen_residual) == 0)
    need("analytic/no_invented_numeric_constant", claims["new_numeric_estimate_constant"] is False)
    need("analytic/no_new_small_time_assumption", claims["new_small_time_assumption"] is False)

    return {
        "schema": "hdblast-parameter-oracle-report-v1", "status": "pass",
        "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(), "checks": records,
        "source_bound_model": {"fields": 16, "spatial_directions": 4, "firstjet_directions": 5,
                               "angular_generators": 6, "Gamma_generators": 4,
                               "mixed_norm_bundle_copies": 2,
                               "scalar_principal": True, "actual_coefficient_arrays_evaluated": False},
        "finite_exact": {"auxiliary_fields": 2, "auxiliary_firstjet_matrix_shape": [5, 2],
                         "state_coefficients_nonlinear": True, "time_angular_quadratic_mix_nonzero": True,
                         "symbolic_coefficient_identity": True, "temporal_row_witness_fields": 2,
                         "nonautonomous_C2_initial_row_witness_last_jet": 4,
                         "temporal_row_witness_last_jet": 14, "rational_integrating_factor_witness_K": "3/5"},
        "initial_dependency_table": dependency_table,
        "scope": "Finite exact nonlinear auxiliary coefficient, principal-product, residual-subtraction, mean-value and integrating-factor identities, together with source-bound dimension/norm/initial-dependency guards. The two-field five-jet witness is not an evaluation of actual 16-field arrays. The finite PDE temporal-row example checks dependency only and does not prove compactness. No analytic estimate, solution existence, actual parameter differentiability, numerical PDE constant, amplitude change, shortened time slab, or observable response is established by this replay."
    }


def controls(claims):
    """Re-evaluate self-rehashed, single-field wrong claims in this process.

    The immutable symbolic auxiliary model is cached independently of claims.
    This is fresh CLAIM evaluation, not fresh subprocess or capsule/manifest
    authentication.  Those guarantees belong to the external replay wrapper.
    """
    records = []
    for field, value in MUTATIONS:
        bad = copy.deepcopy(claims)
        bad[field] = value
        changed = sorted(key for key in claims if claims[key] != bad[key] or type(claims[key]) is not type(bad[key]))
        if changed != [field]:
            raise ValueError("control is not an exact single-field change: " + field)
        raw = canonical(bad)
        failures = []

        def check(name, passed):
            if not passed:
                failures.append(name)

        try:
            evaluate(json.loads(raw), check)
        except ValueError as exc:
            if len(failures) != 1 or any(token in failures[0] for token in ("/schema", "/exact_engine", "/parser/", "/type/")):
                raise ValueError("wrong case did not reach its mathematical/interface diagnostic: " + field) from exc
            records.append({"fault": field, "changed_fields": changed,
                            "wrong_value": value, "self_rehashed_sha256": hashlib.sha256(raw).hexdigest(),
                            "failed_check": failures[0], "diagnostic": str(exc),
                            "rejected_by": "fresh claim evaluation: exact identity or explicit source/dependency guard"})
        else:
            raise ValueError("surviving parameter fault: " + field)
    return {"schema": "hdblast-parameter-oracle-controls-v1", "status": "pass", "controls": records}
