#!/usr/bin/env python3
"""Portable exact algebra checks, not verification of the analytic/physical proof.

Print a deterministic JSON receipt to stdout. Pass --output PATH to write the
same bytes to an explicitly chosen receipt file. No source-side receipt is
written by default. Failed identities or unrejected wrong formulas return 1;
CLI, dependency, or receipt-writing errors return 2. Checks remain active under
python -O.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys

try:
    import sympy as s
except ImportError:
    sys.stderr.write("ERROR: SymPy is required to run the exact algebra checks.\n")
    raise SystemExit(2)


def build_receipt():
    """Check zero residuals and reject wrong formulas by exact witnesses."""
    checks = {}
    residuals = {}
    controls = {}

    def identity(name, residual):
        reduced = s.simplify(residual)
        checks[name] = bool(reduced == 0)
        residuals[name] = s.sstr(reduced)
        return reduced

    def wrong_formula(name, residual, witness, formula, witness_description):
        reduced = s.simplify(residual)
        evaluated = s.simplify(witness)
        # A concrete exact nonzero value certifies that the residual is not
        # identically zero; mere failure of simplify() to return zero does not.
        rejected = bool(
            reduced != 0
            and evaluated.is_number is True
            and evaluated.is_zero is False
        )
        controls[name] = {
            "rejected": rejected,
            "symbolic_residual": s.sstr(reduced),
            "witness_residual": s.sstr(evaluated),
            "wrong_formula": formula,
            "witness": witness_description,
        }

    y = s.symbols("y", positive=True)
    rho = s.Function("rho")(y)
    p = s.Function("p")(y)
    F = s.Function("F")
    x = s.Function("x")(y)
    psi = F(x)
    psi_x = s.diff(psi, x)
    alpha = 1 / (rho**2 * p**2)
    w = 1 / (rho**4 * p**2)
    q = 2 / (3 * rho**2)
    # p has fixed sign on the interval. Replacing p by |p| in h changes
    # only its constant sign, leaving its logarithmic derivatives unchanged.
    h = 1 / (rho ** s.Rational(3, 2) * p)
    u = psi / h
    replace_x = {
        s.diff(x, y, 2): s.diff(rho, y) / rho**2,
        s.diff(x, y): -1 / rho,
    }
    tau = (-s.diff(alpha * s.diff(u, y), y) + q * u) / w
    transformed = s.simplify((h * tau).subs(replace_x))
    logder = s.simplify(-rho * s.diff(h, y) / h)
    V = s.simplify(q / w - rho * s.diff(logder, y) + logder**2)
    expected = -s.diff(psi, x, 2) + V * psi

    # Preserve all ten identities from the original authoring checker.
    identity("unitary_bulk_measure", h**2 - w * rho)
    identity("liouville_no_first_derivative", transformed - expected)
    identity(
        "liouville_log_derivative",
        logder - (s.Rational(3, 2) * s.diff(rho, y)
                  + rho * s.diff(p, y) / p),
    )
    flux = s.simplify((alpha * s.diff(u, y)).subs(replace_x))
    identity("shell_flux_sign", flux - h * (-psi_x + logder * psi))

    A, L, Ap, Lp, P, R, Upp = s.symbols(
        "a ell a_prime ell_prime p rho U_second", real=True
    )
    V_coeff = R**2 * (
        s.Rational(3, 4) * A**2 + 2 * A * L + L**2
        - s.Rational(3, 2) * Ap - Lp + s.Rational(2, 3) * P**2
    )
    V_background = s.expand(
        V_coeff.subs({Lp: Upp - 4 * Ap - 4 * A * L - L**2})
        .subs(Ap, -1 / R**2 - P**2 / 3)
    )
    identity(
        "background_potential_simplification",
        V_background - (-s.Rational(5, 2) + R**2 * (
            s.Rational(3, 4) * A**2 + 6 * A * L + 2 * L**2
            - Upp - P**2 / 6
        )),
    )
    Ac = s.symbols("A_c", nonzero=True, real=True)
    cone_V = s.simplify(V.subs({rho: y, p: Ac * y}).doit())
    identity("cone_essential_threshold", s.limit(cone_V, y, 0) - s.Rational(25, 4))
    identity("shifted_threshold", s.Rational(25, 4) - 4 - s.Rational(9, 4))
    U1 = s.symbols("U_prime", nonzero=True, real=True)
    identity("jacobi_normalization_cancellation", (-U1 / 60)**2 / (U1 / 5)**2 - s.Rational(1, 144))

    g = s.Function("g")(y)
    identity(
        "power_hardy_integration_identity",
        s.diff(y**3 * g**2, y) - 3 * y**2 * g**2
        - 2 * y**3 * g * s.diff(g, y),
    )
    P0 = s.Function("P0")(y)
    f = s.Function("f")(y)
    a = s.Function("alpha")(y)
    ww = s.Function("w")(y)
    qq = 4 * ww + s.diff(a * s.diff(P0, y), y) / P0
    lhs = a * s.diff(P0 * f, y)**2 + (qq - 4 * ww) * (P0 * f)**2
    primitive = a * P0 * s.diff(P0, y) * f**2
    rhs = a * P0**2 * s.diff(f, y)**2 + s.diff(primitive, y)
    identity("picone_integrand_identity", lhs - rhs)

    # Separately test the full spectral-parameter substitution, not only the
    # threshold's numerical subtraction.
    Lambda, mu_sq = s.symbols("Lambda mu_squared", real=True)
    spectral_residual = expected - Lambda * psi
    shifted_residual = -s.diff(psi, x, 2) + (V - 4) * psi - mu_sq * psi
    identity(
        "lambda_to_mu_squared_operator_identity",
        spectral_residual.subs(Lambda, mu_sq + 4) - shifted_residual,
    )

    # Integrating the checked total derivative gives J(T)-J(lower). Retain
    # the lower flux symbol explicitly; a cone-core vector sets it to zero.
    # Subtraction of 4 times the endpoint mass is a separate operation.
    aT, P0T, beta = s.symbols("alpha_T P0_T beta", positive=True)
    P0pT, fT, E_f, J_lower = s.symbols(
        "P0_prime_T f_T transformed_bulk_energy lower_flux", real=True
    )
    endpoint_u = P0T * fT
    endpoint_flux = aT * P0T * P0pT * fT**2
    integrated_shifted_form = E_f + endpoint_flux - J_lower - 4 * endpoint_u**2 / beta
    declared_boundary = aT * P0pT / P0T - 4 / beta
    declared_shifted_form = E_f - J_lower + declared_boundary * endpoint_u**2
    identity("picone_boundary_mass_identity", integrated_shifted_form - declared_shifted_form)

    # Independently integrate a generic Picone fixture. Its q=0 is deliberate:
    # this is not asserted to be a background satisfying the physical q law.
    fixture_alpha, fixture_w = y**-4, y**-6
    fixture_P0, fixture_f, fixture_beta = y**4, 1 + y, s.Integer(2)
    fixture_u = fixture_P0 * fixture_f
    fixture_q = s.simplify(
        4 * fixture_w
        + s.diff(fixture_alpha * s.diff(fixture_P0, y), y) / fixture_P0
    )
    fixture_energy = s.integrate(
        fixture_alpha * s.diff(fixture_u, y)**2 + fixture_q * fixture_u**2,
        (y, 0, 1),
    )
    fixture_norm = (
        s.integrate(fixture_w * fixture_u**2, (y, 0, 1))
        + fixture_u.subs(y, 1)**2 / fixture_beta
    )
    fixture_k = (
        fixture_alpha * fixture_P0 * s.diff(fixture_P0, y)
        - 4 * fixture_P0**2 / fixture_beta
    ).subs(y, 1)
    fixture_D = (
        s.integrate(fixture_alpha * fixture_P0**2 * s.diff(fixture_f, y)**2, (y, 0, 1))
        + fixture_k * fixture_f.subs(y, 1)**2
    )
    fixture_lower = s.limit(
        fixture_alpha * fixture_P0 * s.diff(fixture_P0, y) * fixture_f**2, y, 0
    )
    identity("picone_integrated_boundary_mass_fixture", fixture_energy - 4 * fixture_norm - fixture_D)

    bad_flux = s.simplify(flux - h * (psi_x + logder * psi))
    wrong_formula(
        "liouville_flux_derivative_sign",
        bad_flux,
        bad_flux.subs({psi_x: 1, rho: 1, p: 1}, simultaneous=True),
        "alpha*u_y = h*(+psi_x+s*psi)",
        "rho=1, p=1, psi_x=1 at a regular point",
    )
    wrong_V = q / w + rho * s.diff(logder, y) + logder**2
    bad_potential = s.simplify(V - wrong_V)
    wrong_formula(
        "liouville_potential_derivative_sign",
        bad_potential,
        bad_potential.subs({rho: y, p: y * (1 + y**2)}).doit().subs(y, 1),
        "V = q/w + rho*s_y + s**2 (instead of -rho*s_y)",
        "rho(y)=y, p(y)=y*(1+y**2), y=1; s is not constant",
    )
    bad_shift = s.simplify(shifted_residual - spectral_residual.subs(Lambda, mu_sq - 4))
    wrong_formula(
        "lambda_to_mu_squared_shift_sign",
        bad_shift,
        bad_shift.subs(psi, 1),
        "Lambda = mu_squared - 4",
        "psi=1; all differential and potential terms cancel in the residual",
    )
    wrong_boundary = aT * P0pT / P0T + 4 / beta
    bad_boundary = s.simplify(
        integrated_shifted_form - (E_f - J_lower + wrong_boundary * endpoint_u**2)
    )
    wrong_formula(
        "picone_boundary_mass_sign",
        bad_boundary,
        bad_boundary.subs({P0T: 1, fT: 2, beta: 2}),
        "r0 = alpha_T*P0_prime_T/P0_T + 4/beta",
        "P0_T=1, f_T=2, beta=2, as in the integrated Picone fixture",
    )

    passed = sum(checks.values())
    rejected = sum(control["rejected"] for control in controls.values())
    return {
        "schema_version": 1,
        "scope": (
            "Exact algebra and rejection of specified wrong formulas only; "
            "not formal proof verification, analytic domain verification, "
            "or validation of the full physical perturbation problem."
        ),
        "checks": checks,
        "identity_residuals": residuals,
        "passed": passed,
        "failed": len(checks) - passed,
        "wrong_formula_controls": controls,
        "wrong_formulas_rejected": rejected,
        "wrong_formulas_not_rejected": len(controls) - rejected,
        "all_passed": passed == len(checks) and rejected == len(controls),
        "cone_potential_exact_model": s.sstr(cone_V),
        "picone_integrated_fixture": {
            "scope": "Generic Picone algebra fixture; not a physical background coefficient model.",
            "data": "0<y<=1, alpha=y**-4, w=y**-6, P0=y**4, f=1+y, beta=2",
            "q": s.sstr(fixture_q),
            "energy": s.sstr(fixture_energy),
            "norm_squared": s.sstr(fixture_norm),
            "endpoint_penalty": s.sstr(fixture_k),
            "shifted_energy": s.sstr(fixture_energy - 4 * fixture_norm),
            "transformed_energy": s.sstr(fixture_D),
            "lower_boundary_flux": s.sstr(fixture_lower),
        },
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version": s.__version__,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, metavar="PATH",
        help="also write the deterministic JSON receipt to PATH (default: stdout only)",
    )
    args = parser.parse_args(argv)
    try:
        result = build_receipt()
        encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(encoded, encoding="utf-8")
        sys.stdout.write(encoded)
    except Exception as exc:
        sys.stderr.write(f"ERROR: exact algebra checker could not complete: {exc}\n")
        return 2
    if not result["all_passed"]:
        failed = [name for name, passed in result["checks"].items() if not passed]
        failed.extend(
            name for name, control in result["wrong_formula_controls"].items()
            if not control["rejected"]
        )
        sys.stderr.write("ERROR: failed explicit checks: " + ", ".join(failed) + "\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
