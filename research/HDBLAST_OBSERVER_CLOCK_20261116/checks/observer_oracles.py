#!/usr/bin/env python3
"""Small exact observer-interface witnesses, not actual coefficient evaluation.

The coordinate-probe clock and intrinsic normal-geodesic witnesses are
separate constructions. The source model has ambient dimension5, intrinsic dimension4 and cut
dimension3.  Clock and intrinsic-spray witnesses below have two coordinates;
moving metric/cut and tidal witnesses have three ambient coordinates and a
two-dimensional cut.  The five-dimensional Minkowski FW matrix is universal
frame algebra, not an evaluation of the actual five-dimensional curvature.
All continuum common-domain, positivity, degree and regularity results remain
written proofs. No ancestor producer is imported or executed.
"""
from __future__ import annotations

import copy
from functools import lru_cache
import hashlib
import json
import re

import sympy as S
import mpmath


RATIONAL_DEFAULTS = {
    "clock_origin_first_factor": "1",
    "clock_origin_second_cross_factor": "2",
    "clock_origin_second_time_factor": "1",
    "clock_origin_second_acceleration_factor": "1",
    "inverse_parameter_sign": "-1",
    "inverse_Hessian_sign": "-1",
    "inverse_Hessian_mixed_factor": "1",
    "inverse_Hessian_time_factor": "1",
    "past_velocity_sign": "-1",
    "ODE_first_state_factor": "1",
    "ODE_first_parameter_factor": "1",
    "ODE_second_linear_factor": "1",
    "ODE_second_state_Hessian_factor": "1",
    "ODE_second_mixed_factor": "2",
    "ODE_second_parameter_Hessian_factor": "1",
    "Jacobi_velocity_factor": "2",
    "Gauss_acceleration_factor": "1/6",
    "material_first_state_factor": "1",
    "material_second_mixed_factor": "2",
    "material_second_state_Hessian_factor": "1",
    "material_second_acceleration_factor": "1",
    "cut_first_frame_factor": "1",
    "cut_second_metric_factor": "1",
    "cut_second_mixed_factor": "2",
    "cut_second_frame_acceleration_factor": "1",
    "cut_second_frame_pair_factor": "2",
    "density_first_trace_factor": "1/2",
    "density_second_acceleration_trace_factor": "1/2",
    "density_second_trace_square_factor": "1/4",
    "density_second_square_trace_factor": "-1/2",
    "graph_time_sign": "-1",
    "tidal_first_velocity_factor": "1",
    "tidal_first_screen_factor": "1",
    "tidal_second_R1_slot_factor": "2",
    "tidal_second_slot_acceleration_factor": "1",
    "tidal_second_velocity_pair_factor": "2",
    "tidal_second_velocity_screen_factor": "2",
    "tidal_second_screen_pair_factor": "2",
    "sphere_return_factor": "1",
    "FW_velocity_normal_factor": "1/6",
    "FW_normal_velocity_factor": "1/6",
    "area_numerator_second_product_factor": "2",
    "area_mean_second_quotient_factor": "2",
}
MODE_DEFAULTS = {
    "source_model_mode": "ambient5_intrinsic4_cut3_full_norms",
    "q_velocity_mode": "rho_times_proper",
    "q_acceleration_mode": "rho_squared",
    "intrinsic_spray_mode": "intrinsic_LC",
    "graph_tangent_mode": "Gauss_orthogonal",
    "global_graph_gate": "Gauss_timelikeK_positivecut_degreeone",
    "cut_positivity_gate": "uniform_C0_from_initial_gap",
}
MODES = {
    "source_model_mode": {"ambient5_intrinsic4_cut3_full_norms", "radial_scalar_surrogate"},
    "q_velocity_mode": {"rho_times_proper", "proper_r"},
    "q_acceleration_mode": {"rho_squared", "rho"},
    "intrinsic_spray_mode": {"intrinsic_LC", "discard_intrinsic_connection"},
    "graph_tangent_mode": {"Gauss_orthogonal", "arbitrary_transverse"},
    "global_graph_gate": {"Gauss_timelikeK_positivecut_degreeone", "local_differential_only"},
    "cut_positivity_gate": {"uniform_C0_from_initial_gap", "positive_samples_only"},
}
DEFAULTS = {**RATIONAL_DEFAULTS, **MODE_DEFAULTS}
FACTORS = tuple(DEFAULTS)
MUTATIONS = tuple([(name, "0" if value != "0" else "1") for name, value in RATIONAL_DEFAULTS.items()]
                  + [(name, next(value for value in sorted(MODES[name]) if value != default))
                     for name, default in MODE_DEFAULTS.items()])
_RATIONAL_RE = re.compile(r"(?:0|-[1-9][0-9]*|[1-9][0-9]*)(?:/[1-9][0-9]*)?\Z")


def canonical(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def rational(value):
    if type(value) is not str or len(value) > 24 or not _RATIONAL_RE.fullmatch(value):
        raise ValueError("bounded canonical rational string required")
    p, separator, q = value.partition("/")
    numerator, denominator = int(p), int(q) if separator else 1
    if abs(numerator) > 1000000 or denominator > 1000000:
        raise ValueError("rational exceeds parser bound")
    answer = S.Rational(numerator, denominator)
    if str(answer) != value:
        raise ValueError("noncanonical rational")
    return answer


def zero(value):
    return all(S.cancel(entry) == 0 for entry in value) if isinstance(value, S.MatrixBase) else S.cancel(value) == 0


def addmat(items, rows, columns):
    return sum(items, S.zeros(rows, columns))


@lru_cache(maxsize=1)
def clock_model():
    e, s, t = S.symbols("clock_e clock_s clock_t")
    rate = 2 + s / 3 + s ** 2 / 5 + e * (1 + s ** 2) / 7 + e ** 2 * (2 + s) / 11
    origin = S.Rational(1, 5) + e / 13 + e ** 2 / 17
    offset = S.Rational(1, 19) + e / 23 + e ** 2 / 29
    primitive = S.integrate(rate, s)
    clock = offset + primitive.subs(s, t) - primitive.subs(s, origin)
    point = {e: 0, t: S.Rational(2, 3)}
    at = lambda value: S.cancel(value.subs(point))
    origin0, origin1, origin2 = (S.diff(origin, e, j).subs(e, 0) for j in range(3))
    end = {e: 0, s: origin0}
    first_integral = at(S.integrate(S.diff(rate, e), s).subs(s, t) - S.integrate(S.diff(rate, e), s).subs(s, origin0))
    second_integral = at(S.integrate(S.diff(rate, e, 2), s).subs(s, t) - S.integrate(S.diff(rate, e, 2), s).subs(s, origin0))
    first = [S.diff(offset, e).subs(e, 0) + first_integral, rate.subs(end) * origin1]
    second = [S.diff(offset, e, 2).subs(e, 0) + second_integral,
              S.diff(rate, e).subs(end) * origin1,
              S.diff(rate, s).subs(end) * origin1 ** 2,
              rate.subs(end) * origin2]
    b, ce, cee, cet, ctt = (at(expr) for expr in (S.diff(clock, t), S.diff(clock, e), S.diff(clock, e, 2), S.diff(clock, e, t), S.diff(clock, t, 2)))
    Te, Tr = -ce / b, 1 / b
    true_H = S.Matrix(2, 2, lambda i, j: -((cee if i == j == 0 else 0)
                       + cet * ((Te, Tr)[j] if i == 0 else 0)
                       + cet * ((Te, Tr)[i] if j == 0 else 0)
                       + ctt * (Te, Tr)[i] * (Te, Tr)[j]) / b)
    # Independent coefficient matching of f(e,T(e,r))-r at an exact base jet.
    h, k = S.symbols("inverse_e_delta inverse_tau_delta")
    formal_t = point[t] + h * Te + k * Tr + h ** 2 * true_H[0, 0] / 2 + h * k * true_H[0, 1] + k ** 2 * true_H[1, 1] / 2
    residual = clock.subs({e: h, t: formal_t}, simultaneous=True) - at(clock) - k
    true_inverse_residuals = [S.diff(residual, h, a, k, c).subs({h: 0, k: 0}) for a, c in ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2))]
    return locals()


@lru_cache(maxsize=1)
def intrinsic_model():
    e = S.Symbol("spray_e")
    x = S.Matrix(S.symbols("spray_t spray_x"))
    u = S.Matrix(S.symbols("spray_u0:2"))
    state = S.Matrix([*x, *u])
    eta = S.diag(-1, 1)
    phi = (x[0] + 2 * x[1]) / 13 + (x[0] ** 2 + x[1] ** 2 + x[0] * x[1]) / 5 + e * (x[0] + 2 * x[1]) / 7 + e ** 2 * (x[0] - x[1]) / 11
    dphi = [S.diff(phi, xi) for xi in x]
    Gamma = [[[int(a == b) * dphi[c] + int(a == c) * dphi[b] - eta[b, c] * eta[a, a] * dphi[a]
               for c in range(2)] for b in range(2)] for a in range(2)]
    gamma_contract = lambda left, right: S.Matrix([sum(Gamma[a][b][c] * left[b] * right[c] for b in range(2) for c in range(2)) for a in range(2)])
    spray = S.Matrix([*-u, *gamma_contract(u, u)])
    state0 = S.Matrix([0, 0, S.Rational(5, 3), S.Rational(4, 3)])
    V = S.Matrix([S.Rational(1, 7), S.Rational(2, 9), S.Rational(3, 11), -S.Rational(2, 13)])
    W = S.Matrix([S.Rational(2, 17), -S.Rational(1, 19), S.Rational(3, 23), S.Rational(4, 29)])
    point = {e: 0, **dict(zip(state, state0))}
    at = lambda value: value.subs(point)
    jac = spray.jacobian(state)
    groups = {
        "first_state": at(jac) * V, "first_parameter": at(spray.diff(e)),
        "second_linear": at(jac) * W,
        "second_Hessian": S.Matrix([sum(S.diff(row, state[i], state[j]).subs(point) * V[i] * V[j] for i in range(4) for j in range(4)) for row in spray]),
        "second_mixed": at(jac.diff(e)) * V, "second_parameter": at(spray.diff(e, 2)),
    }
    formal = {state[i]: state0[i] + e * V[i] + e ** 2 * W[i] / 2 for i in range(4)}
    varied = spray.subs(formal, simultaneous=True)
    first_direct, second_direct = (varied.diff(e, j).subs(e, 0) for j in (1, 2))
    P, Q = S.Matrix([S.Rational(2, 7), S.Rational(3, 11)]), S.Matrix([S.Rational(4, 13), -S.Rational(5, 17)])
    directional_Gamma = S.Matrix([sum(S.diff(Gamma[a][b][c], x[i]) * P[i] * u[b] * u[c] for i in range(2) for b in range(2) for c in range(2)) for a in range(2)])
    jacobi_state, jacobi_velocity = at(directional_Gamma), at(gamma_contract(u, Q))
    label = S.Symbol("label_delta")
    label_sub = {**{x[i]: state0[i] + label * P[i] for i in range(2)}, **{u[i]: state0[i + 2] + label * Q[i] for i in range(2)}, e: 0}
    jacobi_position_direct = spray[:2, 0].subs(label_sub, simultaneous=True).diff(label).subs(label, 0)
    jacobi_direct = gamma_contract(u, u).subs(label_sub, simultaneous=True).diff(label).subs(label, 0)
    # Exact coordinate proof of Gauss propagation from metric compatibility,
    # torsion freedom and a unit initial velocity variation.
    U0 = -state0[2:4, 0]
    P0 = S.Matrix([S.Rational(4, 3), S.Rational(5, 3)])
    phi_P = sum(dphi[i].subs(point) * P0[i] for i in range(2))
    Q0 = -phi_P * state0[2:4, 0] + P0 / 7
    metric_derivative = lambda direction: 2 * sum(dphi[i].subs(point) * direction[i] for i in range(2)) * eta
    Uprime = -at(gamma_contract(u, u))
    gauss_derivative = (U0.T * metric_derivative(U0) * P0)[0] + (Uprime.T * eta * P0)[0] + (U0.T * eta * (-Q0))[0]
    compatibility = [2 * dphi[i].subs(point) * eta[j, k] - sum(Gamma[a][i][j].subs(point) * eta[a, k] + Gamma[a][i][k].subs(point) * eta[j, a] for a in range(2)) for i in range(2) for j in range(2) for k in range(2)]
    return locals()


def seed(k):
    return S.Matrix([[S.Rational(k + 1, 11), S.Rational(k + 2, 13), S.Rational(k + 3, 17)],
                     [S.Rational(k + 2, 13), S.Rational(k + 4, 19), S.Rational(k + 5, 23)],
                     [S.Rational(k + 3, 17), S.Rational(k + 5, 23), S.Rational(k + 6, 29)]])


@lru_cache(maxsize=1)
def pullback_model():
    e = S.Symbol("pullback_e")
    H0 = S.Matrix([[-2, 0, 0], [0, 2, S.Rational(1, 5)], [0, S.Rational(1, 5), 3]])
    He, Hee = seed(1), seed(3)
    # Sparse matrices keep a fresh interpreter fast while every moving-point
    # direction and a time/angular mixed Hessian remain genuinely nonzero.
    Mi = [S.Matrix([[1, 1, 0], [1, 2, 0], [0, 0, 0]]) / 13,
          S.Matrix([[0, 0, 0], [0, 1, 1], [0, 1, 2]]) / 17,
          S.Matrix([[3, 0, 1], [0, 0, 0], [1, 0, 1]]) / 19]
    Ni = [S.Matrix([[0, 1, 0], [1, 1, 0], [0, 0, 0]]) / 7,
          S.Matrix([[0, 0, 0], [0, 2, 1], [0, 1, 0]]) / 11,
          S.Matrix([[1, 0, 2], [0, 0, 0], [2, 0, 1]]) / 13]
    Kij = [[S.zeros(3) for _ in range(3)] for _ in range(3)]
    Kij[0][2] = Kij[2][0] = S.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 1]]) / 17
    Kij[1][1] = S.diag(0, S.Rational(1, 19), 0)
    Y, Z = S.Matrix([1, -2, 3]) / 5, S.Matrix([2, 3, -4]) / 7
    X = e * Y + e ** 2 * Z / 2
    H = H0 + e * He + e ** 2 * Hee / 2 + addmat((X[i] * Mi[i] + e * X[i] * Ni[i] for i in range(3)), 3, 3) + addmat((X[i] * X[j] * Kij[i][j] / 2 for i in range(3) for j in range(3)), 3, 3)
    H1_state = addmat((Y[i] * Mi[i] for i in range(3)), 3, 3)
    H2_mixed = addmat((Y[i] * Ni[i] for i in range(3)), 3, 3)
    H2_Hessian = addmat((Y[i] * Y[j] * Kij[i][j] for i in range(3) for j in range(3)), 3, 3)
    H2_acceleration = addmat((Z[i] * Mi[i] for i in range(3)), 3, 3)
    H1, H2 = (H.diff(e, j).subs(e, 0) for j in (1, 2))
    P0 = S.Matrix([[0, 0], [1, 0], [0, 1]])
    P1 = S.Matrix([[1, 2], [3, -1], [2, 3]]) / 7
    P2 = S.Matrix([[2, -1], [3, 2], [-1, 4]]) / 11
    P = P0 + e * P1 + e ** 2 * P2 / 2
    k = P.T * H * P
    k0, k1, k2 = (k.diff(e, j).subs(e, 0) for j in range(3))
    cut_groups = {
        "first_metric": P0.T * H1 * P0,
        "first_frame": P1.T * H0 * P0 + P0.T * H0 * P1,
        "second_metric": P0.T * H2 * P0,
        "second_mixed": P1.T * H1 * P0 + P0.T * H1 * P1,
        "second_acceleration": P2.T * H0 * P0 + P0.T * H0 * P2,
        "second_pair": P1.T * H0 * P1,
    }
    # The reference density is fixed; choose kappa=k0 in this auxiliary model.
    # A determinant's first two derivatives depend only on this exact matrix
    # two-jet; retaining irrelevant powers through e^16 slows fresh replay.
    cut_twojet = k0 + e * k1 + e ** 2 * k2 / 2
    det_ratio = (k0.inv() * cut_twojet).det()
    det1, det2 = (S.diff(det_ratio, e, j).subs(e, 0) for j in (1, 2))
    # Direct sqrt(det) differentiation at det_ratio(0)=1, independent of trace
    # formulas.  Avoid expanding a large square-root expression.
    density1_direct, density2_direct = det1 / 2, det2 / 2 - det1 ** 2 / 4
    A = k0.inv() * k1
    traces = (S.trace(A), S.trace(k0.inv() * k2), S.trace(A) ** 2, S.trace(A * A))
    return locals()


@lru_cache(maxsize=1)
def tidal_model():
    e = S.Symbol("tidal_e")
    metric = S.diag(-1, 1, 1)
    L = S.Matrix([[0, S.Rational(1, 3), S.Rational(2, 5)], [S.Rational(1, 3), 0, S.Rational(1, 7)], [S.Rational(2, 5), -S.Rational(1, 7), 0]])
    B = S.eye(3) + e * L + e ** 2 * L ** 2 / 2
    slots = [B[:, 0], B[:, 1], B[:, 0], B[:, 2]]
    base, first, second = ([v.diff(e, order).subs(e, 0) for v in slots] for order in range(3))
    S0, S1, S2 = seed(2), seed(5), seed(9)
    curvature = S0 + e * S1 + e ** 2 * S2 / 2

    def contract(tensor, vectors):
        x, y, z, w = vectors
        dot = lambda left, mat, right: (left.T * mat * right)[0]
        return dot(x, metric, z) * dot(y, tensor, w) + dot(y, metric, w) * dot(x, tensor, z) - dot(x, metric, w) * dot(y, tensor, z) - dot(y, metric, z) * dot(x, tensor, w)

    replace = lambda i, values: [values if j == i else base[j] for j in range(4)]
    first_velocity = sum(contract(S0, replace(i, first[i])) for i in (0, 2))
    first_screen = sum(contract(S0, replace(i, first[i])) for i in (1, 3))
    second_R1_slot = sum(contract(S1, replace(i, first[i])) for i in range(4))
    second_slot_acceleration = sum(contract(S0, replace(i, second[i])) for i in range(4))

    def pair(i, j):
        values = list(base)
        values[i], values[j] = first[i], first[j]
        return contract(S0, values)

    second_velocity_pair = pair(0, 2)
    second_velocity_screen = sum(pair(i, j) for i, j in ((0, 1), (0, 3), (1, 2), (2, 3)))
    second_screen_pair = pair(1, 3)
    # In the direct contraction only the metric pairings' exact two-jets
    # contribute to derivatives of order <=2.  Keep that two-jet explicitly;
    # differentiating irrelevant degree10 expressions is needlessly costly.
    Gram = B.T * metric * B
    Gram_twojet = addmat((e ** j * Gram.diff(e, j).subs(e, 0) / S.factorial(j) for j in range(3)), 3, 3)
    dotS = lambda i, j: (slots[i].T * curvature * slots[j])[0]
    direct = (Gram_twojet[0, 0] * dotS(1, 3) + Gram_twojet[1, 2] * dotS(0, 2)
              - Gram_twojet[0, 2] * dotS(1, 2) - Gram_twojet[1, 0] * dotS(0, 3))
    first_direct, second_direct = (S.diff(direct, e, j).subs(e, 0) for j in (1, 2))
    R1_base, R2_base = contract(S1, base), contract(S2, base)
    return locals()


@lru_cache(maxsize=1)
def area_mean_model():
    # A genuine exact integral in a one-dimensional auxiliary domain; it is
    # not numerical quadrature or evaluation of the actual S3 area/mean.
    e, x = S.symbols("area_e area_x")
    density = 1 + x ** 2 / 7 + e * (1 + x) / 11 + e ** 2 * (2 + x ** 2) / 13
    scalar = 2 + x / 5 + e * (3 + x ** 2) / 17 + e ** 2 * (1 + x) / 19
    area = S.integrate(density, (x, -1, 1))
    numerator = S.integrate(scalar * density, (x, -1, 1))
    mean = numerator / area
    A0, A1, A2 = (S.diff(area, e, j).subs(e, 0) for j in range(3))
    N0, N1, N2 = (S.diff(numerator, e, j).subs(e, 0) for j in range(3))
    M0, M1, M2 = (S.diff(mean, e, j).subs(e, 0) for j in range(3))
    N2_pure = S.integrate((S.diff(scalar, e, 2) * density + scalar * S.diff(density, e, 2)).subs(e, 0), (x, -1, 1))
    N2_cross = S.integrate((S.diff(scalar, e) * S.diff(density, e)).subs(e, 0), (x, -1, 1))
    return locals()


def evaluate(claims, check):
    records = []

    def need(name, condition):
        passed = bool(condition)
        full = "observer_oracles/" + name
        check(full, passed)
        records.append({"name": name, "status": "pass" if passed else "fail"})
        if not passed:
            raise ValueError("observer oracle rejects " + full)

    need("schema", type(claims) is dict and claims.get("schema") == "hdblast-observer-oracle-claims-v1" and set(claims) == {"schema", *FACTORS})
    need("exact_engine", S.__version__ == "1.14.0" and mpmath.__version__ == "1.3.0")
    values = {}
    for name in RATIONAL_DEFAULTS:
        try:
            values[name] = rational(claims[name])
        except (ValueError, TypeError):
            need("parser/" + name, False)
    for name in MODE_DEFAULTS:
        need("type/" + name, type(claims[name]) is str and claims[name] in MODES[name])
    f = values.__getitem__
    need("source_contract/full_ambient_intrinsic_cut_model", claims["source_model_mode"] == "ambient5_intrinsic4_cut3_full_norms")

    c = clock_model()
    first = c["first"][0] - f("clock_origin_first_factor") * c["first"][1]
    second = (c["second"][0] - f("clock_origin_second_cross_factor") * c["second"][1]
              - f("clock_origin_second_time_factor") * c["second"][2]
              - f("clock_origin_second_acceleration_factor") * c["second"][3])
    need("clock/moving_origin_first", zero(first - c["ce"]))
    need("clock/moving_origin_full_second", zero(second - c["cee"]))
    b = c["b"]
    Te, Tr = f("inverse_parameter_sign") * c["ce"] / b, 1 / b
    need("clock/inverse_first_parameter", zero(b * Te + c["ce"]))
    need("clock/inverse_first_signed_clock", zero(b * Tr - 1))
    T = (Te, Tr)
    Hessian = S.Matrix(2, 2, lambda i, j: f("inverse_Hessian_sign") / b *
        ((c["cee"] if i == j == 0 else 0)
         + f("inverse_Hessian_mixed_factor") * c["cet"] * ((T[j] if i == 0 else 0) + (T[i] if j == 0 else 0))
         + f("inverse_Hessian_time_factor") * c["ctt"] * T[i] * T[j]))
    need("clock/inverse_full_joint_Hessian", zero(Hessian - c["true_H"]))
    need("clock/independent_implicit_twojet_coefficients", all(zero(v) for v in c["true_inverse_residuals"]))
    need("clock/nonzero_parameter_clock_cross_witness", c["cet"] != 0 and c["ctt"] != 0 and c["origin1"] != 0 and c["origin2"] != 0)
    need("clock/positive_auxiliary_rate", b > 0)

    intrinsic = intrinsic_model()
    eta = S.diag(-1, 1)
    future = intrinsic["state0"][2:4, 0]
    past = f("past_velocity_sign") * future
    need("proper_duration/past_velocity_is_minus_future", zero(past + future))
    need("proper_duration/r_has_unit_timelike_norm", zero((past.T * eta * past)[0] + 1))
    rho = S.Rational(3, 2)  # Auxiliary scale, not a new numerical observer bound.
    qvelocity = (rho if claims["q_velocity_mode"] == "rho_times_proper" else 1) * past
    need("proper_duration/q_velocity_retains_rho", zero(qvelocity - rho * past))
    need("proper_duration/q_is_not_unit_proper_duration", zero((qvelocity.T * eta * qvelocity)[0] + rho ** 2))
    sigma = S.Symbol("sigma")
    normal = S.Matrix([0, 0, 1])
    ambient_acceleration = sigma * normal / 6
    qacceleration = (rho ** 2 if claims["q_acceleration_mode"] == "rho_squared" else rho) * ambient_acceleration
    need("proper_duration/q_acceleration_has_rho_squared", zero(qacceleration - rho ** 2 * ambient_acceleration))

    spray0 = intrinsic["spray"].subs(intrinsic["point"])
    selected_spray = spray0 if claims["intrinsic_spray_mode"] == "intrinsic_LC" else S.Matrix([*spray0[:2, 0], 0, 0])
    need("intrinsic/geodesic_uses_intrinsic_connection", zero(selected_spray - spray0) and not zero(spray0[2:4, 0]) and zero(spray0[:2, 0] + future))
    g = intrinsic["groups"]
    first_ODE = f("ODE_first_state_factor") * g["first_state"] + f("ODE_first_parameter_factor") * g["first_parameter"]
    second_ODE = (f("ODE_second_linear_factor") * g["second_linear"] + f("ODE_second_state_Hessian_factor") * g["second_Hessian"]
                  + f("ODE_second_mixed_factor") * g["second_mixed"] + f("ODE_second_parameter_Hessian_factor") * g["second_parameter"])
    need("intrinsic/first_full_parameter_ODE", zero(first_ODE - intrinsic["first_direct"]))
    need("intrinsic/second_full_parameter_ODE", zero(second_ODE - intrinsic["second_direct"]))
    need("intrinsic/all_state_direction_components_present", all(v != 0 for v in intrinsic["V"]) and all(not zero(v) for v in g.values()))
    need("intrinsic/coordinate_Jacobi_velocity_pair", zero(intrinsic["jacobi_state"] + f("Jacobi_velocity_factor") * intrinsic["jacobi_velocity"] - intrinsic["jacobi_direct"]) and zero(intrinsic["jacobi_position_direct"] + intrinsic["Q"]))
    need("intrinsic/fixed_metric_compatibility", all(zero(v) for v in intrinsic["compatibility"]))
    need("intrinsic/Gauss_orthogonality_propagation", zero(intrinsic["gauss_derivative"]) and zero((intrinsic["U0"].T * eta * intrinsic["P0"])[0]))
    produced_acceleration = f("Gauss_acceleration_factor") * sigma * normal
    need("Gauss/normal_acceleration_plus_sigma_over_six", zero(produced_acceleration - ambient_acceleration))
    need("Gauss/past_and_future_have_same_normal_acceleration", zero((-1) ** 2 * ambient_acceleration - ambient_acceleration))

    p = pullback_model()
    H1 = p["He"] + f("material_first_state_factor") * p["H1_state"]
    H2 = (p["Hee"] + f("material_second_mixed_factor") * p["H2_mixed"]
          + f("material_second_state_Hessian_factor") * p["H2_Hessian"] + f("material_second_acceleration_factor") * p["H2_acceleration"])
    need("pullback/material_first_all_state_slots", zero(H1 - p["H1"]))
    need("pullback/material_full_second", zero(H2 - p["H2"]))
    need("pullback/nonzero_all_moving_point_directions", all(not zero(v) for v in p["Mi"]) and all(y != 0 for y in p["Y"]) and not zero(p["H2_Hessian"]) and not zero(p["Kij"][0][2]))
    kg = p["cut_groups"]
    produced_k1 = p["P0"].T * H1 * p["P0"] + f("cut_first_frame_factor") * kg["first_frame"]
    produced_k2 = (f("cut_second_metric_factor") * p["P0"].T * H2 * p["P0"]
                   + f("cut_second_mixed_factor") * kg["second_mixed"]
                   + f("cut_second_frame_acceleration_factor") * kg["second_acceleration"]
                   + f("cut_second_frame_pair_factor") * kg["second_pair"])
    need("pullback/cut_metric_first", zero(produced_k1 - p["k1"]))
    need("pullback/cut_metric_full_second", zero(produced_k2 - p["k2"]))
    trA, trK2, trA_squared, trAA = p["traces"]
    density1 = f("density_first_trace_factor") * trA
    density2 = (f("density_second_acceleration_trace_factor") * trK2
                + f("density_second_trace_square_factor") * trA_squared + f("density_second_square_trace_factor") * trAA)
    need("density/sqrt_det_first_trace", zero(density1 - p["density1_direct"]))
    need("density/sqrt_det_full_second_trace", zero(density2 - p["density2_direct"]))
    need("density/trace_square_differs_from_square_trace", not zero(trA_squared - trAA))
    need("density/positive_nondiagonal_auxiliary_cut", p["k0"][0, 0] > 0 and p["k0"].det() > 0 and p["k0"][0, 1] != 0)

    U = S.Matrix([-S.Rational(5, 3), S.Rational(4, 3)])
    tangent = S.Matrix([S.Rational(4, 3), -S.Rational(5, 3)])
    if claims["graph_tangent_mode"] == "arbitrary_transverse":
        tangent += U / 7
    angular_compensation = -U[1] / tangent[1]
    graph_time = U[0] + tangent[0] * angular_compensation
    Ktime = S.Matrix([1, 0])
    pairing = (Ktime.T * eta * U)[0]
    need("graph/Gauss_orthogonality_required", zero((U.T * eta * tangent)[0]))
    need("graph/noncrossing_time_identity", zero(graph_time * pairing - f("graph_time_sign")))
    need("graph/past_time_graph_strict_decrease", pairing > 0 and graph_time < 0)
    need("graph/global_injectivity_dependency_contract", claims["global_graph_gate"] == "Gauss_timelikeK_positivecut_degreeone")
    need("cut/uniform_positivity_dependency_contract", claims["cut_positivity_gate"] == "uniform_C0_from_initial_gap")

    tidal = tidal_model()
    t1 = tidal["R1_base"] + f("tidal_first_velocity_factor") * tidal["first_velocity"] + f("tidal_first_screen_factor") * tidal["first_screen"]
    t2 = (tidal["R2_base"] + f("tidal_second_R1_slot_factor") * tidal["second_R1_slot"]
          + f("tidal_second_slot_acceleration_factor") * tidal["second_slot_acceleration"]
          + f("tidal_second_velocity_pair_factor") * tidal["second_velocity_pair"]
          + f("tidal_second_velocity_screen_factor") * tidal["second_velocity_screen"]
          + f("tidal_second_screen_pair_factor") * tidal["second_screen_pair"])
    need("tidal/first_all_four_moving_slots", zero(t1 - tidal["first_direct"]))
    need("tidal/second_all_four_slots_and_six_pairs", zero(t2 - tidal["second_direct"]))
    need("tidal/nonzero_velocity_and_screen_pair_witnesses", all(not zero(tidal[key]) for key in ("first_velocity", "first_screen", "second_R1_slot", "second_slot_acceleration", "second_velocity_pair", "second_velocity_screen", "second_screen_pair")))
    need("tidal/frame_Gram_preserved_to_second_order", zero(tidal["L"].T * tidal["metric"] + tidal["metric"] * tidal["L"]) and zero(tidal["Gram_twojet"] - tidal["metric"]))

    omega = S.Matrix(S.symbols("sphere_omega0:3"))
    velocity = S.Matrix(S.symbols("sphere_U0:3"))
    fixed = S.Matrix([S.Rational(1, 3), S.Rational(2, 5), S.Rational(3, 7)])
    tangential_spray = omega.cross(fixed)
    omega_r = -velocity
    velocity_r = tangential_spray + f("sphere_return_factor") * velocity.dot(velocity) * omega
    c1, c2 = omega.dot(omega) - 1, omega.dot(velocity)
    c1r, c2r = 2 * omega.dot(omega_r), omega_r.dot(velocity) + omega.dot(velocity_r)
    need("sphere/constraint_first_row", zero(c1r + 2 * c2))
    need("sphere/constraint_return_row", zero(c2r - velocity.dot(velocity) * c1))

    minkowski = S.diag(-1, 1, 1, 1, 1)
    future5, normal5 = S.eye(5)[:, 0], S.eye(5)[:, 4]
    FW = (f("FW_normal_velocity_factor") * sigma * future5 * (normal5.T * minkowski)
          - f("FW_velocity_normal_factor") * sigma * normal5 * (future5.T * minkowski))
    need("FW/future_velocity_to_outward_normal", zero(FW * future5 - sigma * normal5 / 6))
    need("FW/outward_normal_to_future_velocity", zero(FW * normal5 - sigma * future5 / 6))
    need("FW/skew_adjoint_full_five_frame_Gram", zero(FW.T * minkowski + minkowski * FW))
    need("FW/all_three_screen_vectors_intrinsic_parallel", all(zero(FW * S.eye(5)[:, i]) for i in (1, 2, 3)))
    need("FW/past_parameter_retains_acceleration_and_normal_sign", zero((-FW) * (-future5) - sigma * normal5 / 6) and zero((-FW) * normal5 + sigma * future5 / 6))

    area = area_mean_model()
    N2 = area["N2_pure"] + f("area_numerator_second_product_factor") * area["N2_cross"]
    M1 = (area["N1"] - area["M0"] * area["A1"]) / area["A0"]
    M2 = (N2 - area["M0"] * area["A2"] - f("area_mean_second_quotient_factor") * M1 * area["A1"]) / area["A0"]
    need("area_mean/positive_auxiliary_total_area", area["A0"] > 0)
    need("area_mean/second_integrated_product", zero(N2 - area["N2"]))
    need("area_mean/first_positive_area_quotient", zero(M1 - area["M1"]))
    need("area_mean/second_positive_area_quotient", zero(M2 - area["M2"]))
    need("area_mean/nonzero_product_and_quotient_cross_witnesses", area["N2_cross"] != 0 and area["M1"] * area["A1"] != 0)

    return {"schema": "hdblast-observer-oracle-report-v1", "status": "pass",
            "claims_sha256": hashlib.sha256(canonical(claims)).hexdigest(), "checks": records,
            "source_bound_model": {"ambient_dimension": 5, "intrinsic_dimension": 4, "cut_dimension": 3,
                                   "all_angular_modes": True, "full_product_norms": True,
                                   "actual_metric_or_curvature_arrays_evaluated": False},
            "auxiliary_witnesses": {"clock_variables": ["parameter", "coordinate_time", "signed_probe_clock"],
                                    "intrinsic_coordinates": 2, "intrinsic_first_order_state_components": 4,
                                    "moving_metric_ambient_dimension": 3, "moving_cut_dimension": 2,
                                    "tidal_ambient_dimension": 3, "sphere_embedding_dimension": 3,
                                    "universal_Minkowski_FW_dimension": 5,
                                    "area_mean_integration_domain_dimension": 1,
                                    "auxiliary_clock_scale_rho": "3/2"},
            "scope": "Exact small auxiliary signed coordinate-probe clock/inverse two-jets, separate nonlinear intrinsic normal-geodesic spray parameter and coordinate-Jacobi identities, Gauss propagation/umbilic acceleration, moving metric/cut/density derivatives, graph algebra, all-slot tidal multilinearity, sphere constraint rows, universal Minkowski FW frame algebra and positive-area scalar-mean product/quotient identities on an auxiliary interval. The actual source has dimensions5/4/3; the small witnesses do not evaluate its coefficient arrays or S3 integrals. Uniform common short domains, cut positivity, global degree/injectivity, function-space regularity and physical observer conclusions are written proof obligations, represented here only by explicit dependency guards. No numerical proper-duration bound, empirical calibration or continuum theorem is proved by this replay."}


def controls(claims):
    """Fresh single-field claim evaluation in one process, not capsule auth."""
    records = []
    for field, value in MUTATIONS:
        bad = copy.deepcopy(claims)
        bad[field] = value
        changed = [key for key in claims if claims[key] != bad[key] or type(claims[key]) is not type(bad[key])]
        if changed != [field]:
            raise ValueError("control is not single-field: " + field)
        raw, failures = canonical(bad), []
        try:
            evaluate(json.loads(raw), lambda name, passed: failures.append(name) if not passed else None)
        except ValueError as exc:
            if len(failures) != 1 or any(token in failures[0] for token in ("/schema", "/exact_engine", "/parser/", "/type/")):
                raise ValueError("fault did not reach a mathematical/dependency diagnostic: " + field) from exc
            records.append({"fault": field, "wrong_value": value, "changed_fields": changed,
                            "self_rehashed_sha256": hashlib.sha256(raw).hexdigest(),
                            "failed_check": failures[0], "diagnostic": str(exc)})
        else:
            raise ValueError("surviving observer fault: " + field)
    return {"schema": "hdblast-observer-oracle-controls-v1", "status": "pass", "controls": records,
            "scope": "Fresh claim evaluations with cached immutable auxiliary models in one process. Fresh subprocess and closed-capsule authentication belong to the external wrapper."}


