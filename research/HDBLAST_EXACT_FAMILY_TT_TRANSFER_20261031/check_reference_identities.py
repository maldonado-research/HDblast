#!/usr/bin/env python3
"""Bounded exact checks for the anchored-reference lemma.

Authenticates repository sources before importing SymPy. The checks certify
specified identities and reject named wrong formulas; they prove no analytic
existence, regularity, inverse-domain, or exact-family hypothesis.
"""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys


class SourceFailure(Exception):
    pass


class CheckFailure(Exception):
    pass


def authenticate(repo_root):
    receipt_path = Path(__file__).resolve().with_name("PUBLIC_SOURCE_RECEIPT.json")
    try:
        raw = receipt_path.read_bytes()
        receipt = json.loads(raw)
    except (OSError, ValueError) as exc:
        raise SourceFailure("Cannot read the sibling PUBLIC_SOURCE_RECEIPT.json") from exc
    if not isinstance(receipt, dict) or not isinstance(receipt.get("sources"), list):
        raise SourceFailure("Receipt must contain a sources list")
    if not receipt["sources"]:
        raise SourceFailure("Receipt sources list must not be empty")
    root = repo_root.resolve()
    if not root.is_dir():
        raise SourceFailure("Repository root is not a directory")
    checked = []
    seen = set()
    for entry in receipt["sources"]:
        if not isinstance(entry, dict):
            raise SourceFailure("Every source entry must be an object")
        name, expected = entry.get("path"), entry.get("sha256")
        if not isinstance(name, str) or not isinstance(expected, str):
            raise SourceFailure("Every source needs string path and sha256 fields")
        relative = PurePosixPath(name)
        if (relative.is_absolute() or not relative.parts
                or relative.parts[0] != "research" or ".." in relative.parts
                or "\\" in name or str(relative) != name):
            raise SourceFailure("Source path must be normalized and repository-relative: " + name)
        if name in seen:
            raise SourceFailure("Duplicate source path: " + name)
        if re.fullmatch(r"[0-9a-f]{64}", expected) is None:
            raise SourceFailure("Invalid SHA256 for source: " + name)
        seen.add(name)
        candidate = (root / Path(*relative.parts)).resolve()
        try:
            candidate.relative_to(root)
        except ValueError as exc:
            raise SourceFailure("Source escapes repository root: " + name) from exc
        try:
            observed = hashlib.sha256(candidate.read_bytes()).hexdigest()
        except OSError as exc:
            raise SourceFailure("Cannot read source: " + name) from exc
        if observed != expected:
            raise SourceFailure("SHA256 mismatch: " + name)
        checked.append({"path": name, "sha256": observed})
    return {
        "status": "pass",
        "receipt_sha256": hashlib.sha256(raw).hexdigest(),
        "sources": sorted(checked, key=lambda item: item["path"]),
        "source_count": len(checked),
    }


def run_checks(sp):
    checks = {}
    controls = {}

    def equal(label, actual, expected):
        difference = actual - expected
        entries = list(difference) if isinstance(difference, sp.MatrixBase) else [difference]
        if any(sp.simplify(entry) != 0 for entry in entries):
            raise CheckFailure(label)

    def reject(name, label, callback):
        try:
            callback()
        except CheckFailure as exc:
            if str(exc) != label:
                raise CheckFailure("Unexpected rejection for " + name) from exc
            controls[name] = {"status": "rejected", "failed_identity": label}
        else:
            raise CheckFailure("Wrong formula was accepted: " + name)

    e, perturbation = sp.symbols("e U", real=True)
    rho = sp.symbols("rho", positive=True)
    speed = 1 / (rho * (1 + e * perturbation))
    speed_jets = [sp.diff(speed, e, order).subs(e, 0) for order in range(3)]
    expected_jets = [1 / rho, -perturbation / rho, 2 * perturbation**2 / rho]

    def check_speed(candidate):
        for order in range(3):
            equal("reciprocal_lapse_normal_velocity_jets",
                  sp.diff(candidate, e, order).subs(e, 0), expected_jets[order])

    check_speed(speed)
    checks["reciprocal_lapse_normal_velocity_jets"] = {
        "status": "pass", "orders": [0, 1, 2], "jets": [str(v) for v in speed_jets]
    }
    reject("fixed_background_lapse_speed", "reciprocal_lapse_normal_velocity_jets",
           lambda: check_speed(1 / rho))

    # The independent symbols represent the two inverse-metric contractions.
    physical, target = sp.symbols("contract_Gamma_g contract_Gamma_g0")
    gauge = physical - target

    def check_tension(candidate):
        equal("identity_wave_map_tension_equals_minus_C", candidate, -gauge)

    check_tension(-physical + target)
    checks["identity_wave_map_tension_equals_minus_C"] = {
        "status": "pass", "identity": "tau(id)=-C; C=contract_Gamma_g-contract_Gamma_g0"
    }
    reject("wrong_target_connection_sign", "identity_wave_map_tension_equals_minus_C",
           lambda: check_tension(-physical - target))

    # Four tangential coordinates (time first), radial coordinate last.
    # H has positive spatial block and negative time Schur complement.
    H = sp.Matrix([[-4, 1, 0, 0], [1, 3, 1, 0],
                   [0, 1, 2, 0], [0, 0, 0, 5]])
    Hinv = H.inv()
    b = sp.Matrix(sp.symbols("b0:4"))
    c = sp.symbols("c", nonzero=True)
    schur = c - (b.T * Hinv * b)[0]
    covariant = H.row_join(b).col_join(b.T.row_join(sp.Matrix([[c]])))
    mixed = -Hinv * b / schur
    inverse = (Hinv + Hinv * b * b.T * Hinv / schur).row_join(mixed)
    inverse = inverse.col_join(mixed.T.row_join(sp.Matrix([[1 / schur]])))
    equal("mixed_metric_block_inverse", covariant * inverse, sp.eye(5))
    # Under invertible H and nonzero Schur complement, mixed=0 iff b=0.
    equal("mixed_metric_zero_equivalence", -schur * H * mixed, b)
    orthogonal_metric = sp.diag(1, 1, 1, 1, 9)
    orthogonal_metric[:4, :4] = H
    orthogonal_inverse = orthogonal_metric.inv()
    outward = -orthogonal_inverse[:, 4] / sp.sqrt(orthogonal_inverse[4, 4])
    observer = -orthogonal_inverse[:, 0] / sp.sqrt(-orthogonal_inverse[0, 0])
    equal("shell_normal_tangential_labels_zero", outward[:4, :], sp.zeros(4, 1))
    equal("shell_time_normal_is_tangent", observer[4], 0)
    equal("shell_spacelike_unit_normal", (outward.T * orthogonal_metric * outward)[0], 1)
    equal("shell_future_unit_time_normal", (observer.T * orthogonal_metric * observer)[0], -1)
    checks["mixed_metric_boundary_correspondence"] = {
        "status": "pass",
        "scope": "symbolic off-diagonal vector with fixed Lorentzian tangential block; nonzero Schur complement",
        "identities": ["g*g_inverse=I", "-Schur*H*g_inverse_tangent_radial=b",
                       "n(Phi^a)=0 at id iff g_Za=0", "U^Z=0 when g_Za=0"],
    }

    # Q is the positive square root of A=Q^2; Sylvester minors 2,3,4 prove Q>0.
    Q = sp.Matrix([[2, 1, 0], [1, 2, 1], [0, 1, 2]])
    minors = [Q[:i, :i].det() for i in range(1, 4)]
    if any(value <= 0 for value in minors):
        raise CheckFailure("positive_frame_fixture")
    A, B = Q * Q, Q.inv()
    D = sp.Matrix([[1, 1, 0], [0, 2, 1], [0, 0, 1]])
    R = sp.Rational(3, 2)
    physical_screen = R**2 * D.inv().T * A * D.inv()
    equal("screen_pullback_metric", D.T * physical_screen * D, R**2 * A)

    def check_frame(candidate):
        equal("normalized_screen_frame_isometry",
              candidate.T * physical_screen * candidate, sp.eye(3))

    check_frame(D * B / R)
    checks["normalized_screen_frame_isometry"] = {
        "status": "pass", "scope": "exact rational positive-definite fixture",
        "positive_square_root_principal_minors": [int(v) for v in minors],
        "identity": "J^T s J=I, J=(1/R)*D*B, B=(sbar/R^2)^(-1/2)",
    }
    reject("missing_reference_radius_frame_factor", "normalized_screen_frame_isometry",
           lambda: check_frame(D * B))

    x = sp.symbols("x", real=True)
    phi = x**2 / 2 + (1 + e + e**2) * x
    event = -(1 + e + e**2) + sp.sqrt((1 + e + e**2)**2 + 8)
    equal("explicit_inverse_event", phi.subs(x, event), 4)
    at_event = {e: 0, x: 2}
    jacobian = sp.diff(phi, x).subs(at_event)
    first = -sp.diff(phi, e).subs(at_event) / jacobian
    acceleration_terms = (
        sp.diff(phi, e, 2).subs(at_event),
        2 * sp.diff(phi, e, x).subs(at_event) * first,
        sp.diff(phi, x, 2).subs(at_event) * first**2,
    )
    equal("inverse_event_first_derivative", first, sp.diff(event, e).subs(e, 0))

    def check_inverse_acceleration(candidate):
        equal("inverse_event_second_derivative", candidate, sp.diff(event, e, 2).subs(e, 0))

    second = -sum(acceleration_terms) / jacobian
    check_inverse_acceleration(second)
    if any(term == 0 for term in acceleration_terms):
        raise CheckFailure("inverse_event_fixture_must_exercise_all_terms")
    checks["inverse_event_second_derivative"] = {
        "status": "pass", "scope": "explicit nonaffine one-dimensional varying map",
        "map": "Phi_e(x)=x^2/2+(1+e+e^2)*x", "target_event": 4,
        "event_at_zero": 2, "p_prime": str(first), "p_second": str(second),
        "acceleration_terms": [str(v) for v in acceleration_terms],
    }
    reject("missing_inverse_acceleration_cross_term", "inverse_event_second_derivative",
           lambda: check_inverse_acceleration(-(acceleration_terms[0] + acceleration_terms[2]) / jacobian))

    h = sp.symbols("h")
    gram = sp.diag(sp.Rational(1, 6), sp.Rational(1, 6))
    tensor_coefficients = sp.Matrix([h, -h])
    pairings = gram * tensor_coefficients
    recovered = 6 * pairings
    equal("orthogonal_template_projection", recovered, tensor_coefficients)

    def check_coefficient(candidate):
        equal("selected_coefficient_second_derivative", candidate, h)

    selected = (recovered[0] - recovered[1]) / 2
    check_coefficient(selected)
    checks["projection_and_second_derivative_normalization"] = {
        "status": "pass", "scope": "normalization identity using the sourced normalized L2 Gram matrix",
        "gram_diagonal": ["1/6", "1/6"], "gram_cross": "0",
        "identity": "Pi H''=h(Tplus-Tminus) implies mathfrak_h''=h",
    }
    reject("extra_half_in_second_derivative_coefficient", "selected_coefficient_second_derivative",
           lambda: check_coefficient(selected / 2))
    if len(controls) != 5 or any(row["status"] != "rejected" for row in controls.values()):
        raise CheckFailure("mutant_control_roster")
    return checks, controls


def emit(result, output):
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if output is not None:
        output.write_text(serialized, encoding="utf-8")
    else:
        sys.stdout.write(serialized)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        source_receipt = authenticate(args.repo_root)
    except SourceFailure as exc:
        emit({"status": "fail", "stage": "source_authentication", "error": str(exc),
              "sympy_imported": "sympy" in sys.modules}, args.output)
        return 2

    # Authentication intentionally precedes this first non-stdlib import.
    try:
        import sympy as sp
    except ImportError:
        emit({"status": "fail", "stage": "dependency", "error": "SymPy is required",
              "source_authentication": source_receipt}, args.output)
        return 2
    try:
        checks, controls = run_checks(sp)
    except CheckFailure as exc:
        emit({"status": "fail", "stage": "identity_checks", "error": str(exc),
              "source_authentication": source_receipt}, args.output)
        return 1
    emit({
        "status": "pass", "schema_version": 1,
        "scope": "Bounded exact identities and wrong-formula rejection only; no analytic existence or exact-family certification",
        "source_authentication": source_receipt, "sympy_version": sp.__version__,
        "identities": checks, "mutant_controls": controls,
        "assert_independent": True,
    }, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
