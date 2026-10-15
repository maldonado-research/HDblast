#!/usr/bin/env python3
"""Finite exact identity checks; these do not prove the PDE estimates."""
import json
import sys
from pathlib import Path

import sympy as sp


def demand(condition, label):
    if not condition:
        raise RuntimeError(label)


def fixture(seed, defect=None):
    n = 4
    q = sp.Rational
    gamma = sp.eye(n) + sp.Matrix([q(i + 1, seed + 9) for i in range(n)]) * sp.Matrix(
        [[q(i + 1, seed + 9) for i in range(n)]]
    )
    inverse = gamma.inv()
    kappa = sp.Matrix(n, n, lambda i, j: q((i + 1) * (j + 1) + seed, 43))
    alpha = q(seed + 4, seed + 3)
    beta = sp.Matrix([q((-1) ** i * (seed + i + 1), 29) for i in range(n)])
    alpha_space = [q(seed + j + 1, 31) for j in range(n)]
    beta_space = [sp.Matrix([q((j + 1) * (i + 2) - seed, 47) for i in range(n)]) for j in range(n)]
    gamma_space = [sp.Matrix(n, n, lambda i, j: q((i + 1) * (j + 1) + seed + d, 53)) for d in range(n)]
    christoffel_spatial = [
        sum(
            inverse[i, ell] * inverse[j, k]
            * (gamma_space[j][k, ell] + gamma_space[k][j, ell] - gamma_space[ell][j, k])
            for ell in range(n) for j in range(n) for k in range(n)
        ) / 2
        for i in range(n)
    ]
    metric = sp.zeros(n + 1)
    metric[0, 0] = -alpha ** 2 + (beta.T * gamma * beta)[0]
    for i in range(n):
        metric[0, i + 1] = metric[i + 1, 0] = (gamma * beta)[i]
        for j in range(n):
            metric[i + 1, j + 1] = gamma[i, j]
    metric_inverse = metric.inv()
    reference = [
        sp.Matrix(n + 1, n + 1, lambda b, c: q((a + 2) * (b + c + 3) + seed, 71))
        for a in range(n + 1)
    ]
    target = [sum(metric_inverse[b, c] * reference[a][b, c] for b in range(n + 1) for c in range(n + 1)) for a in range(n + 1)]
    trace_k = sp.trace(inverse * kappa)
    sign_k = -1 if defect == "wrong_extrinsic_sign" else 1
    sign_reference = 1 if defect == "wrong_reference_sign" else -1
    alpha_time = sum(beta[j] * alpha_space[j] for j in range(n)) + sign_k * alpha ** 2 * trace_k + sign_reference * alpha ** 3 * target[0]
    beta_time = []
    for i in range(n):
        lapse_factor = 2 if defect == "wrong_shift_lapse_factor" else 1
        beta_time.append(
            sum(beta[j] * beta_space[j][i] for j in range(n))
            + alpha ** 2 * (christoffel_spatial[i] - target[i + 1] - beta[i] * target[0])
            - lapse_factor * alpha * sum(inverse[i, j] * alpha_space[j] for j in range(n))
        )
    beta_time = sp.Matrix(beta_time)
    gamma_time = sp.Matrix(n, n, lambda i, j:
        2 * alpha * kappa[i, j]
        + sum(beta[k] * gamma_space[k][i, j] + gamma[k, j] * beta_space[i][k] + gamma[i, k] * beta_space[j][k] for k in range(n)))
    derivatives = []
    for d in range(n + 1):
        ga = gamma_time if d == 0 else gamma_space[d - 1]
        al = alpha_time if d == 0 else alpha_space[d - 1]
        be = beta_time if d == 0 else beta_space[d - 1]
        dg = sp.zeros(n + 1)
        dg[0, 0] = -2 * alpha * al + (beta.T * ga * beta)[0] + 2 * (beta.T * gamma * be)[0]
        for i in range(n):
            dg[0, i + 1] = dg[i + 1, 0] = (ga * beta + gamma * be)[i]
            for j in range(n):
                dg[i + 1, j + 1] = ga[i, j]
        derivatives.append(dg)
    actual = [
        sum(
            metric_inverse[a, d] * metric_inverse[b, c]
            * (derivatives[b][c, d] + derivatives[c][b, d] - derivatives[d][b, c])
            for d in range(n + 1) for b in range(n + 1) for c in range(n + 1)
        ) / 2
        for a in range(n + 1)
    ]
    residual = [sp.factor(actual[a] - target[a]) for a in range(n + 1)]
    return residual


def main():
    positive = []
    rejected = []
    for seed in (1, 2, 4):
        residual = fixture(seed)
        demand(all(value == 0 for value in residual), "exact ADM harmonic completion")
        positive.append(f"nonzero_shift_exact_harmonic_fixture_{seed}")
    for defect in ("wrong_extrinsic_sign", "wrong_reference_sign", "wrong_shift_lapse_factor"):
        residual = fixture(2, defect)
        demand(any(value != 0 for value in residual), f"failed to reject {defect}")
        rejected.append(defect)

    alpha, sigma, N, rho, cp, dp, Ca, Cay = sp.symbols("alpha sigma N rho sigma_prime Dphi Ca Cay", nonzero=True)
    # Four-dimensional shell: Codazzi derivative minus scalar stress.
    codazzi = (sp.Rational(1, 6) - sp.Rational(4, 6)) * cp * dp
    scalar_stress = -cp * dp / 2
    demand(sp.expand(codazzi - scalar_stress) == 0, "field-dependent-tension Codazzi cancellation")
    positive.append("variable_tension_codazzi_scalar_cancellation")
    wrong_codazzi = -scalar_stress
    demand(wrong_codazzi != 0, "failed to reject omitted tension derivative")
    rejected.append("omit_field_dependent_tension_derivative")
    K = sigma / 6
    covariant_robin = Cay / N - K * Ca - K * Ca
    coordinate_robin = (Cay - N * sigma * Ca / 3) / N
    demand(sp.expand(covariant_robin - coordinate_robin) == 0, "coordinate subsidiary Robin factor two")
    positive.append("coordinate_subsidiary_robin_factor_two")
    demand(sp.expand(covariant_robin - (Cay - N * sigma * Ca / 6) / N) != 0, "failed to reject halved Robin")
    rejected.append("halved_coordinate_robin")
    # Initial C=0: all tangential derivatives vanish, so trace E=-Xu.
    Xu = sp.symbols("Xu")
    traceE = -Xu
    Huu = Xu - sp.Rational(1, 2) * (-1) * traceE
    demand(sp.expand(Huu - Xu / 2) == 0, "initial trace-reversed Hamiltonian projection")
    positive.append("initial_trace_reversed_Huu")
    demand(sp.expand(Huu - Xu) != 0, "failed to reject untraced Einstein Hamiltonian")
    rejected.append("replace_Huu_with_Euu")
    lapse = sp.symbols("lapse", positive=True)
    beta = sp.symbols("b0:4")
    jac = sp.zeros(5)
    jac[0, 0] = -lapse ** -3
    for i in range(4):
        jac[i + 1, 0] = beta[i] * lapse ** -3
        jac[i + 1, i + 1] = -lapse ** -2
    demand(sp.factor(jac.det() + lapse ** -11) == 0, "harmonic timevelocity map invertible")
    positive.append("unique_harmonic_lapse_shift_velocities")
    output = {
        "status": "PASS",
        "optimized_python": not __debug__,
        "positive_checks": positive,
        "wrong_formula_controls_rejected": rejected,
        "positive_count": len(positive),
        "rejected_count": len(rejected),
        "scope": "Exact finite algebra only; no analytic PDE or existence theorem is certified by this script.",
    }
    destination = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("GAUGE_SUBSIDIARY_CHECK_RESULTS.json")
    destination.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: output[key] for key in ("status", "optimized_python", "positive_count", "rejected_count")}))


if __name__ == "__main__":
    main()
