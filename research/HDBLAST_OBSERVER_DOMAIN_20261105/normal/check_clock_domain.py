#!/usr/bin/env python3
"""Exact finite algebra checks for OBSERVER_CLOCK_DOMAIN.md; no physical inputs."""

import argparse
import itertools
import json
import os
from pathlib import Path
from fractions import Fraction as Q
import sys


class CheckFailure(Exception):
    pass


CHECKS = []


def check(condition, label):
    if not condition:
        raise CheckFailure(label)
    CHECKS.append(label)


class Poly:
    """Sparse commutative Laurent polynomials with rational coefficients."""

    def __init__(self, terms=None):
        self.terms = {k: Q(v) for k, v in (terms or {}).items() if v}

    @staticmethod
    def const(value):
        return Poly({(): Q(value)})

    @staticmethod
    def atom(name, exponent=1):
        return Poly({((name, exponent),): Q(1)})

    def __add__(self, other):
        other = aspoly(other)
        result = dict(self.terms)
        for key, value in other.terms.items():
            result[key] = result.get(key, Q(0)) + value
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-aspoly(other))

    def __rsub__(self, other):
        return aspoly(other) - self

    def __mul__(self, other):
        other = aspoly(other)
        result = {}
        for left, lv in self.terms.items():
            for right, rv in other.terms.items():
                powers = dict(left)
                for name, power in right:
                    powers[name] = powers.get(name, 0) + power
                key = tuple(sorted((name, power) for name, power in powers.items() if power))
                result[key] = result.get(key, Q(0)) + lv * rv
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Only nonnegative integer powers are supported")
        result = Poly.const(1)
        for unused in range(exponent):
            result = result * self
        return result

    def __eq__(self, other):
        return self.terms == aspoly(other).terms


def aspoly(value):
    return value if isinstance(value, Poly) else Poly.const(value)


def sign(perm):
    return -1 if sum(perm[i] > perm[j] for i in range(len(perm))
                    for j in range(i + 1, len(perm))) % 2 else 1


def determinant(matrix):
    size = len(matrix)
    result = Poly.const(0)
    for perm in itertools.permutations(range(size)):
        term = Poly.const(sign(perm))
        for row, col in enumerate(perm):
            term = term * matrix[row][col]
        result = result + term
    return result


def minor(matrix, row, col):
    return [[v for j, v in enumerate(line) if j != col]
            for i, line in enumerate(matrix) if i != row]


def gauss_contraction_identity():
    metric = [[Poly.atom("H%d%d" % (min(i, j), max(i, j)))
               for j in range(4)] for i in range(4)]
    u = [Poly.atom("U%d" % i) for i in range(4)]
    p = [[Poly.atom("P%d%d" % (i, j)) for j in range(3)] for i in range(4)]
    b = [Poly.atom("B%d" % i) for i in range(3)]
    v = [u[i] - sum((p[i][j] * b[j] for j in range(3)), Poly.const(0))
         for i in range(4)]

    def pairing(x, y):
        return sum((x[i] * metric[i][j] * y[j]
                    for i in range(4) for j in range(4)), Poly.const(0))

    orthogonal_terms = sum((b[j] * pairing([p[i][j] for i in range(4)], u)
                            for j in range(3)), Poly.const(0))
    check(pairing(v, u) - pairing(u, u) + orthogonal_terms == 0,
          "universal_gauss_fixed_angle_contraction_identity")
    rho = Poly.atom("rho")
    norm = Poly.atom("future_unit_norm")
    check((-rho) * (-rho) * norm == rho * rho * norm,
          "proper_past_clock_squared_norm_scaling")
    check(-rho + rho == 0, "proper_past_clock_velocity_sign")
    return {"metric_dimension": 4, "label_dimension": 3,
            "analytic_premises_not_machine_proved": [
                "geodesic preservation of normality", "Lorentzian opposite-orientation pairing positive",
                "positive cut metric", "future timelike coordinate vector"]}


def schur_identity():
    angular = [[Poly.atom("A%d%d" % (i, j)) for j in range(3)] for i in range(3)]
    b = [Poly.atom("b%d" % i) for i in range(3)]
    c = [Poly.atom("c%d" % i) for i in range(3)]
    d = Poly.atom("d")
    full = [[d] + c] + [[b[i]] + angular[i] for i in range(3)]
    adj = [[((-1) ** (i + j)) * determinant(minor(angular, j, i))
            for j in range(3)] for i in range(3)]
    rhs = d * determinant(angular)
    for i in range(3):
        for j in range(3):
            rhs = rhs - c[i] * adj[i][j] * b[j]
    check(determinant(full) == rhs, "universal_4x4_schur_identity")
    return {"variables": 16, "expanded_determinant_terms": len(determinant(full).terms)}


class Jet:
    """Two-variable degree-two jets whose coefficients are exact polynomials."""

    def __init__(self, terms=None):
        self.terms = {k: aspoly(v) for k, v in (terms or {}).items()
                      if aspoly(v) != 0 and sum(k) <= 2}

    @staticmethod
    def const(value):
        return Jet({(0, 0): value})

    def __add__(self, other):
        if not isinstance(other, Jet):
            other = Jet.const(other)
        result = dict(self.terms)
        for key, value in other.terms.items():
            result[key] = result.get(key, Poly.const(0)) + value
        return Jet(result)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, Jet):
            other = Jet.const(other)
        result = {}
        for left, lv in self.terms.items():
            for right, rv in other.terms.items():
                key = tuple(left[i] + right[i] for i in range(2))
                if sum(key) <= 2:
                    result[key] = result.get(key, Poly.const(0)) + lv * rv
        return Jet(result)

    __rmul__ = __mul__

    def coefficient(self, key):
        return self.terms.get(key, Poly.const(0))


def inverse_and_composition_jets():
    names = "a b c d e f g h k we wn wq wee wen wnn weq wnq wqq".split()
    atoms = {name: Poly.atom(name) for name in names}
    a, b, c, d, e, f, g, h, k = (atoms[name] for name in "a b c d e f g h k".split())
    ai = Poly.atom("a", -1)
    qe, qn = -b * ai, -c * ai
    qee = -(d + 2 * g * qe + k * qe * qe) * ai
    qen = -(e + g * qn + h * qe + k * qe * qn) * ai
    qnn = -(f + 2 * h * qn + k * qn * qn) * ai
    ey, ny = Jet({(1, 0): 1}), Jet({(0, 1): 1})
    dq = Jet({(1, 0): qe, (0, 1): qn, (2, 0): Q(1, 2) * qee,
              (1, 1): qen, (0, 2): Q(1, 2) * qnn})
    residual = (dq * a + ey * b + ny * c + ey * ey * (d * Q(1, 2))
                + ey * ny * e + ny * ny * (f * Q(1, 2))
                + ey * dq * g + ny * dq * h + dq * dq * (k * Q(1, 2)))
    for key in ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2)):
        check(residual.coefficient(key) == 0, "inverse_clock_jet_%d_%d" % key)
    w = atoms
    composed = (ey * w["we"] + ny * w["wn"] + dq * w["wq"]
                + ey * ey * (w["wee"] * Q(1, 2)) + ey * ny * w["wen"]
                + ny * ny * (w["wnn"] * Q(1, 2))
                + ey * dq * w["weq"] + ny * dq * w["wnq"]
                + dq * dq * (w["wqq"] * Q(1, 2)))
    expected = {
        (1, 0): w["we"] + w["wq"] * qe,
        (0, 1): w["wn"] + w["wq"] * qn,
        (2, 0): (w["wee"] + 2 * w["weq"] * qe + w["wqq"] * qe * qe
                 + w["wq"] * qee) * Q(1, 2),
        (1, 1): w["wen"] + w["weq"] * qn + w["wnq"] * qe
                 + w["wqq"] * qe * qn + w["wq"] * qen,
        (0, 2): (w["wnn"] + 2 * w["wnq"] * qn + w["wqq"] * qn * qn
                 + w["wq"] * qnn) * Q(1, 2),
    }
    for key, value in expected.items():
        check(composed.coefficient(key) == value, "angular_composition_jet_%d_%d" % key)
    qt, qtt = ai, -k * ai ** 3
    check(a * qt == 1, "inverse_coordinate_time_first_derivative")
    check(a * qtt + k * qt * qt == 0, "inverse_coordinate_time_second_derivative")
    return {"formal_variables": len(names), "inverse_jet_coefficients": 5,
            "angular_composition_coefficients": 5, "time_derivatives": 2}


def valid_guard(eta, sigma, alpha, beta, delta, margin):
    return (0 <= eta < 1 and sigma > 0 and alpha >= 0 and beta >= 0
            and delta >= 0 and margin > 0 and delta < margin
            and alpha * beta < (1 - eta) * sigma)


def guard_and_margins():
    ell = Q(1, 10 ** 6)
    tc, tq, qs = Q(13, 8) * ell, Q(25, 16) * ell, Q(6145, 4096) * ell
    slab, target = tc - tq, tq - qs
    check(slab == ell / 16, "exact_coordinate_slab_margin")
    check(target == Q(255, 4096) * ell, "exact_target_proper_clock_margin")
    check(slab - target == ell / 4096, "distinct_margins_differ_by_one_4096_ell")
    check(slab / tq == Q(1, 25), "derivative_only_slab_clock_threshold")
    check(target / tq == Q(255, 6400), "derivative_only_target_clock_threshold")
    eta, sigma, alpha, beta, delta = Q(1, 100), Q(7, 8), ell / 100, 1 / (10 * ell), ell / 100
    mu = 1 - eta
    c = mu - alpha * beta / sigma
    fixedtime_gap = sigma - alpha * beta / mu
    check(valid_guard(eta, sigma, alpha, beta, delta, slab), "synthetic_strict_guard_passes")
    check(delta < target, "synthetic_target_common_strip_passes")
    check(c > 0 and fixedtime_gap > 0, "graph_and_fixedtime_gaps_positive")
    check(c * sigma == fixedtime_gap * mu, "same_graph_and_fixedtime_guard")
    cases = {
        "clock_endpoint": (Q(1), sigma, alpha, beta, delta, slab),
        "angular_endpoint": (eta, Q(0), alpha, beta, delta, slab),
        "slab_endpoint": (eta, sigma, alpha, beta, slab, slab),
        "schur_endpoint": (eta, sigma, mu * sigma, Q(1), delta, slab),
        "negative_alpha": (eta, sigma, Q(-1), beta, delta, slab),
        "negative_delta": (eta, sigma, alpha, beta, Q(-1), slab),
    }
    for name, values in cases.items():
        check(not valid_guard(*values), "strict_guard_rejects_" + name)
    return {"ell": str(ell), "T_coord": str(tc), "T_q": str(tq),
            "q_star": str(qs), "slab_margin": str(slab), "target_margin": str(target),
            "derivative_only_eta_slab": "1/25", "derivative_only_eta_target": "255/6400",
            "synthetic_inputs_not_physical": {"eta_t": str(eta), "sigma": str(sigma),
                "alpha": str(alpha), "beta": str(beta), "delta_t": str(delta)},
            "synthetic_graph_gap": str(c), "synthetic_fixedtime_gap": str(fixedtime_gap),
            "strict_rejection_cases": sorted(cases)}


def counterexample():
    # q and theta are represented as exact rational coefficients of pi.
    q1, theta1, sin1 = Q(1, 4), Q(1, 2), Q(1)
    q2, theta2, sin2 = Q(1, 8), Q(1), Q(0)
    t1, t2 = -q1 * (1 - sin1 / 2), -q2 * (1 - sin2 / 2)
    w1, w2 = theta1 + 4 * q1, theta2 + 4 * q2
    check(q1 != q2 and theta1 != theta2, "counterexample_distinct_domain_points")
    check(t1 == t2 == Q(-1, 8), "counterexample_exact_time_collision_pi_units")
    check(w1 == w2 == Q(3, 2), "counterexample_exact_angle_collision_pi_units")
    check(1 - Q(1, 4) == Q(3, 4), "counterexample_uniform_cut_metric_lower_bound")
    check(16 - Q(3, 2) ** 2 == Q(55, 4), "counterexample_is_spacelike_not_proper_observer")
    # The two q-values are below 1 using the stated analytic fact 0<pi<4.
    check(q1 * 4 == 1 and q2 * 4 < 1, "counterexample_domain_using_pi_less_than_four")
    fixedtime_derivative = 1 + Q(4, 5) * Q(-4, 5) / (1 - Q(3, 5) / 2) ** 2
    check(fixedtime_derivative == Q(-15, 49), "counterexample_fixedtime_angular_jacobian_negative")
    return {"type": "smooth_map_not_normal_or_proper_geodesic",
            "q1_pi_units": str(q1), "q2_pi_units": str(q2),
            "theta1_pi_units": str(theta1), "theta2_pi_units": str(theta2),
            "collision_t_pi_units": str(t1), "collision_angle_pi_units": str(w1),
            "cut_metric_lower_bound": "3/4", "velocity_squared_lower_bound": "55/4",
            "negative_fixedtime_jacobian": str(fixedtime_derivative),
            "external_analytic_facts": ["0 < pi < 4", "sin(pi/2)=1", "sin(pi)=0"]}


def fresh_absolute_file(value):
    path = Path(value)
    if not path.is_absolute():
        raise CheckFailure("--output must be absolute")
    if path.name in ("", ".", "..") or ".." in path.parts:
        raise CheckFailure("--output must be a normalized file path")
    current = Path(path.anchor)
    for part in path.parts[1:-1]:
        current = current / part
        if current.is_symlink() or not current.is_dir():
            raise CheckFailure("--output parent chain must exist and contain no symlink")
    if path.is_symlink() or path.exists():
        raise CheckFailure("--output must name a fresh non-symlink file")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        destination = fresh_absolute_file(args.output)
        gauss = gauss_contraction_identity()
        schur = schur_identity()
        jets = inverse_and_composition_jets()
        guards = guard_and_margins()
        example = counterexample()
        report = {"schema": "hdblast.clock_domain.independent.v1", "status": "PASS",
                  "arithmetic": "exact fractions and formal Laurent polynomials",
                  "physical_amplitude_evaluated": False, "pde_verified": False,
                  "global_topological_theorem_machine_proved": False,
                  "checks_count": len(CHECKS), "checks": CHECKS,
                  "gauss_contraction_identity": gauss, "schur_identity": schur, "second_jets": jets,
                  "guards_and_margins": guards, "counterexample": example}
        raw = (json.dumps(report, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n").encode("ascii")
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(destination, flags, 0o600)
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(raw)
        sys.stdout.buffer.write(raw)
        return 0
    except (CheckFailure, OSError, ValueError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
