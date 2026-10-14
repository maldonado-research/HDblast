#!/usr/bin/env python3
"""Exact frozen GH shell-symbol checks; no curved IBVP theorem is inferred."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--inject-error", action="store_true")
    args = parser.parse_args()
    checks: list[str] = []
    controls: list[str] = []
    control_witnesses: dict[str, str] = {}
    lam, s, k, r = sp.symbols("lambda s k r", nonzero=True)
    eta = [-1, 1, 1, 1]
    q = [s, sp.I * k, sp.Integer(0), sp.Integer(0)]
    norm_subs = {
        s: lam * (1 - r**2) / (1 + r**2),
        k: 2 * lam * r / (1 + r**2),
    }
    a, w = sp.symbols("a w")
    b = sp.symbols("b0:4")
    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    cs = sp.symbols(" ".join(f"c{i}{j}" for i, j in pairs))
    c = sp.zeros(4)
    for (i, j), value in zip(pairs, cs):
        c[i, j] = c[j, i] = value
    trc = sum(eta[i] * c[i, i] for i in range(4))
    qb = sum(eta[i] * q[i] * b[i] for i in range(4))
    trace_factor = sp.Rational(1, 3) if args.inject_error else sp.Rational(1, 2)
    cn = lam * a + qb - trace_factor * lam * (a + trc)
    ca = [
        lam * b[i]
        + sum(eta[j] * q[j] * c[i, j] for j in range(4))
        - trace_factor * q[i] * (a + trc)
        for i in range(4)
    ]
    physical_metric = [
        lam * c[i, j] - q[i] * b[j] - q[j] * b[i]
        for i, j in pairs
    ]
    physical = physical_metric + [lam * w]
    variables = [a, *b, *cs, w]
    israel_solution = {
        c[i, j]: (q[i] * b[j] + q[j] * b[i]) / lam
        for i, j in pairs
    }

    def canonical(expr: sp.Expr) -> sp.Expr:
        return sp.factor(sp.cancel(expr.subs(norm_subs, simultaneous=True)))

    def zero(name: str, *expressions: sp.Expr) -> None:
        for expression in expressions:
            residual = canonical(expression)
            if residual != 0:
                raise RuntimeError(f"{name}: nonzero residual {residual}")
        checks.append(name)

    def equal_integer(name: str, actual: int, expected: int) -> None:
        if actual != expected:
            raise RuntimeError(f"{name}: found {actual}, expected {expected}")
        checks.append(name)

    def rejected(name: str, wrong_residual: sp.Expr) -> None:
        # A single exact nonzero witness refutes each deliberately false identity.
        # Avoid multivariate factorization of a known nonzero control polynomial.
        witness = wrong_residual.subs(norm_subs, simultaneous=True).subs({r: sp.Rational(1, 2), lam: 2})
        witness = sp.simplify(witness.subs({value: 1 for value in witness.free_symbols}))
        if witness == 0:
            raise RuntimeError(f"wrong-formula control unexpectedly passed: {name}")
        controls.append(name)
        control_witnesses[name] = str(witness)

    def matrix(rows: list[sp.Expr], columns: list[sp.Symbol]) -> sp.Matrix:
        return sp.Matrix(rows).jacobian(columns)

    receipt: dict[str, object] = {
        "scope": "frozen principal symbol only; no curved Einstein-scalar IBVP theorem",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version": sp.__version__,
    }
    try:
        zero("shell wave-root relation", -s**2 - k**2 + lam**2)
        zero(
            "normal harmonic trace-reversal identity after Israel",
            cn.subs(israel_solution) - lam * a / 2,
        )
        zero(
            "tangential harmonic trace-reversal identities after Israel",
            *[ca[i].subs(israel_solution) + q[i] * a / 2 for i in range(4)],
        )
        full_physical = matrix(physical, variables)
        full_naive = matrix([cn, *ca, *physical], variables)
        physical_pivot = full_physical[:, [variables.index(value) for value in [*cs, w]]]
        zero("full physical nonzero maximal minor", physical_pivot.det() - lam**11)
        equal_integer("full physical rank", physical_pivot.rank(), 11)
        # A determinant-one unknown change eliminates the Israel shift terms.
        # This avoids costly generic symbolic row reduction before q^2=-lambda^2.
        full_change = matrix([a, *b, *[value + israel_solution[value] for value in cs], w], variables)
        zero("full Israel elimination change has determinant one", full_change.det() - 1)
        full_naive_eliminated = (full_naive * full_change).applyfunc(canonical)
        equal_integer(
            "full physical-plus-all-harmonic rank defect",
            full_naive_eliminated.rank(),
            12,
        )
        zero("full naive boundary determinant vanishes", full_naive_eliminated.det())
        kernel_vectors = []
        for index in range(4):
            bk = [sp.Integer(int(j == index)) for j in range(4)]
            ck = [(q[i] * bk[j] + q[j] * bk[i]) / lam for i, j in pairs]
            kernel_vectors.append(sp.Matrix([0, *bk, *ck, 0]))
        zero(
            "four tangential harmonic pure-gauge kernel modes",
            *[entry for vector in kernel_vectors for entry in full_naive * vector],
        )
        equal_integer(
            "four kernel modes are independent",
            sp.Matrix.hstack(*kernel_vectors).rank(),
            4,
        )
        completed_full = matrix([cn, *b, *physical], variables)
        zero(
            "completed full candidate determinant",
            completed_full.det(method="domain-ge") - lam**12 / 2,
        )
        equal_integer(
            "completed full candidate rank",
            completed_full.subs({s: sp.Rational(3, 5), k: sp.Rational(4, 5), lam: 1}).rank(),
            16,
        )

        bt, bl, ctt, ctl, cll, ctrans = sp.symbols("bt bl ctt ctl cll ctrans")
        scalar_vars = [a, bt, bl, ctt, ctl, cll, ctrans, w]
        scalar_subs = {
            b[0]: bt, b[1]: bl, b[2]: 0, b[3]: 0,
            c[0, 0]: ctt, c[0, 1]: ctl,
            c[0, 2]: 0, c[0, 3]: 0,
            c[1, 1]: cll, c[1, 2]: 0, c[1, 3]: 0,
            c[2, 2]: ctrans, c[2, 3]: 0, c[3, 3]: ctrans,
        }
        scalar_physical = [
            (lam * c[0, 0] - 2 * q[0] * b[0]).subs(scalar_subs),
            (lam * c[0, 1] - q[0] * b[1] - q[1] * b[0]).subs(scalar_subs),
            (lam * c[1, 1] - 2 * q[1] * b[1]).subs(scalar_subs),
            lam * ctrans,
            lam * w,
        ]
        scalar_constraints = [cn.subs(scalar_subs), ca[0].subs(scalar_subs), ca[1].subs(scalar_subs)]
        zero(
            "transverse harmonic rows vanish in scalar block",
            ca[2].subs(scalar_subs), ca[3].subs(scalar_subs),
        )
        scalar_physical_matrix = matrix(scalar_physical, scalar_vars)
        scalar_naive = matrix([*scalar_constraints, *scalar_physical], scalar_vars)
        scalar_pivot = scalar_physical_matrix[:, [scalar_vars.index(value) for value in [ctt, ctl, cll, ctrans, w]]]
        zero("scalar physical nonzero maximal minor", scalar_pivot.det() - lam**5)
        equal_integer("scalar physical rank", scalar_pivot.rank(), 5)
        scalar_change = matrix([
            a, bt, bl, ctt + 2 * s * bt / lam,
            ctl + (s * bl + sp.I * k * bt) / lam,
            cll + 2 * sp.I * k * bl / lam, ctrans, w,
        ], scalar_vars)
        zero("scalar Israel elimination change has determinant one", scalar_change.det() - 1)
        scalar_naive_eliminated = (scalar_naive * scalar_change).applyfunc(canonical)
        equal_integer(
            "scalar physical-plus-all-harmonic rank defect",
            scalar_naive_eliminated.rank(),
            6,
        )
        zero("scalar naive determinant vanishes", scalar_naive_eliminated.det())
        scalar_completed = matrix([scalar_constraints[0], bt, bl, *scalar_physical], scalar_vars)
        zero(
            "completed scalar candidate determinant",
            scalar_completed.det(method="domain-ge") - lam**6 / 2,
        )
        equal_integer(
            "completed scalar candidate rank",
            scalar_completed.subs({s: sp.Rational(3, 5), k: sp.Rational(4, 5), lam: 1}).rank(),
            8,
        )

        an = sp.Symbol("an")
        bnn = sp.symbols("bnn0:4")
        cns = sp.symbols(" ".join(f"cn{i}{j}" for i, j in pairs))
        cnormal = sp.zeros(4)
        for (i, j), value in zip(pairs, cns):
            cnormal[i, j] = cnormal[j, i] = value
        trcn = sum(eta[i] * cnormal[i, i] for i in range(4))
        qsq = sum(eta[i] * q[i] ** 2 for i in range(4))
        cn_boundary = an / 2 + qb - trcn / 2
        ca_normal = [
            bnn[i] + sum(eta[j] * q[j] * cnormal[i, j] for j in range(4))
            - q[i] * (an + trcn) / 2
            for i in range(4)
        ]
        codazzi_twice = [
            sum(eta[j] * q[j] * cnormal[i, j] for j in range(4))
            - q[i] * trcn - qsq * b[i] + q[i] * qb
            for i in range(4)
        ]
        zero(
            "subsidiary normal-tangential identity modulo reduced mixed wave equation",
            *[
                codazzi_twice[i] - ca_normal[i] - q[i] * cn_boundary
                + bnn[i] + qsq * b[i]
                for i in range(4)
            ],
        )
        neumann_subs = {an: 0, **{value: 0 for value in cns}, **{value: 0 for value in b}, **{value: 0 for value in bnn}}
        zero(
            "candidate Dirichlet normal / Neumann tangential subsidiary principal boundary",
            cn_boundary.subs(neumann_subs),
            *[value.subs(neumann_subs) for value in ca_normal],
        )
        d, tension_prime, field_gradient = sp.symbols("d sigma_prime phi_gradient")
        zero(
            "variable-tension Codazzi scalar-normal cancellation",
            (1 - d) * tension_prime / (2 * (d - 1)) * field_gradient
            - (-tension_prime / 2) * field_gradient,
        )
        lapse, lapse_normal, umbilic = sp.symbols("N N_y Lambda", positive=True)
        absolute_gamma = lapse_normal / lapse**3 - d * umbilic / lapse
        zero(
            "absolute harmonic d-dependent normal metric coefficient",
            absolute_gamma.subs(lapse_normal, d * umbilic * lapse**2),
        )
        zero(
            "five-dimensional absolute harmonic coefficient is eight",
            (2 * lapse * d * umbilic * lapse**2).subs(d, 4) - 8 * umbilic * lapse**3,
        )
        zero(
            "Gaussian normal background has nonzero absolute harmonic normal vector",
            absolute_gamma.subs({lapse: 1, lapse_normal: 0, d: 4}) + 4 * umbilic,
        )
        tension, bg_umbilic, hh0trace = sp.symbols("sigma a_background tr_h_h0")
        relative_gamma = lapse_normal / lapse**3 - 2 * tension / (3 * lapse) + bg_umbilic * hh0trace
        zero(
            "exact relative-wave-map normal lapse boundary formula",
            relative_gamma.subs(lapse_normal, 2 * tension * lapse**2 / 3 - bg_umbilic * lapse**3 * hh0trace),
        )
        eps, ctrace = sp.symbols("epsilon ctrace")
        relative_linear = sp.diff(relative_gamma.subs({
            lapse: 1 + eps * a / 2,
            lapse_normal: eps * an / 2,
            tension: 6 * bg_umbilic + eps * tension_prime * w,
            hh0trace: 4 - eps * ctrace,
        }), eps).subs(eps, 0)
        physical_trace_normal = 2 * bg_umbilic * ctrace + 4 * bg_umbilic * a + 4 * tension_prime * w / 3
        zero(
            "relative normal gauge linearization agrees with covariant formula after Israel",
            relative_linear - ((an - physical_trace_normal) / 2 + 4 * bg_umbilic * a),
        )
        tangential_value, tangential_normal = sp.symbols("tangential_value tangential_normal")
        normal_covariant = tangential_normal - bg_umbilic * tangential_value
        tangent_normal_covariant = -bg_umbilic * tangential_value
        zero(
            "shell-preserving covector coordinate gauge Robin coefficient two",
            normal_covariant + tangent_normal_covariant
            - (tangential_normal - 2 * bg_umbilic * tangential_value),
        )
        zero(
            "tangential subsidiary covector Robin coefficient two",
            normal_covariant - bg_umbilic * tangential_value
            - (tangential_normal - 2 * bg_umbilic * tangential_value),
        )
        epsilon = sp.Symbol("epsilon", positive=True)
        glancing_root = sp.sqrt((epsilon + sp.I) ** 2 + 1)
        zero("glancing normal wave root tends to zero", sp.limit(glancing_root, epsilon, 0, dir="+"))

        rejected(
            "wrong harmonic trace coefficient",
            (lam * a + qb - sp.Rational(1, 3) * lam * (a + trc)).subs(israel_solution)
            - lam * a / 2,
        )
        rejected("wrong Israel shift sign", (lam * c[0, 1] + q[0] * b[1] + q[1] * b[0]).subs(israel_solution))
        rejected(
            "wrong subsidiary tangential normal-condition sign",
            codazzi_twice[0] - ca_normal[0] + q[0] * cn_boundary
            + bnn[0] + qsq * b[0],
        )
        rejected("wrong scalar determinant exponent", scalar_completed.det(method="domain-ge") - lam**8 / 2)
        receipt.update({
            "status": "passed",
            "positive_check_count": len(checks),
            "wrong_formula_controls_rejected": len(controls),
            "checks": checks,
            "wrong_formula_controls": controls,
            "wrong_formula_control_nonzero_witnesses": control_witnesses,
            "full": {"incoming_wave_fields": 16, "physical_rank": 11, "physical_and_all_harmonic_rank": 12, "tangential_gauge_kernel_dimension": 4, "candidate_completed_rank": 16, "candidate_completed_determinant": "lambda**12/2"},
            "scalar_l_ge_2": {"incoming_wave_fields": 8, "physical_rank": 5, "physical_and_all_harmonic_rank": 6, "tangential_gauge_kernel_dimension": 2, "candidate_completed_rank": 8, "candidate_completed_determinant": "lambda**6/2"},
            "uniformity_limit": "Neumann lambda symbol vanishes at glancing in the unweighted stable amplitude norm; the flat homogeneous energy remains conserved",
        })
    except Exception as exc:
        receipt.update({"status": "failed", "reason": str(exc), "completed_checks": checks, "deliberate_error_injected": args.inject_error})
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        print(str(exc), file=sys.stderr)
        return 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": receipt["status"], "positive_check_count": len(checks), "wrong_formula_controls_rejected": len(controls)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
