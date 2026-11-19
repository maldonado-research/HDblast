"""First elastic change of the comoving one-particle moments in the dilute contact model."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import expit, gamma, roots_genlaguerre, roots_legendre, zeta


ARMS = ("quiet", "plus", "minus", "far")
BANDS = ("low", "high", "full")


def moment_vector(record):
    return np.array([record["one_sign_comoving_moments"][str(m)] for m in range(5)])


def incoming(moments):
    # The contact collision probability in each energy interval is proportional to k dN.
    return moments[2:5] / moments[1]


def event_change(moments_a, moments_b):
    """Symmetric average of the change for a single participant in one A-B collision."""
    a, b, c = incoming(moments_a)
    d, e, f = incoming(moments_b)
    second = -(b + e) / 6 + 5 * a * d / 12
    third = -(c + f) / 4 + 3 * (b * d + a * e) / 8
    return np.array([0.0, 0.0, second, third])


def positive_rule(values, weights, n=3):
    """Gauss rule from positive atomic input using orthonormal polynomials."""
    values = np.ravel(values)
    weights = np.ravel(weights)
    measure = weights / np.sum(weights)
    previous = np.zeros_like(values)
    current = np.ones_like(values)
    mat = np.zeros((n, n))
    beta = 0.
    for j in range(n):
        product = values * current - beta * previous
        alpha = np.sum(measure * current * product)
        mat[j, j] = alpha
        if j + 1 < n:
            residual = product - alpha * current
            beta = np.sqrt(np.sum(measure * residual**2))
            mat[j, j + 1] = mat[j + 1, j] = beta
            previous, current = current, residual / beta
    nodes, vec = np.linalg.eigh(mat)
    rule = vec[0]**2
    for m in range(2*n):
        np.testing.assert_allclose(np.sum(rule * nodes**m),
                                   np.sum(measure * values**m), rtol=2e-13)
    return nodes, rule


def projected_rule(q, omega, pair):
    """Direct Dirac two-body projection; no use of archived single-particle moments."""
    z, w = roots_legendre(4)
    k = (omega[:, :, None] + q[:, None, None] * z) / 2
    weights = pair[:, :, None] * w / 2
    if np.min(weights) < 0 or np.min(k) <= 0:
        raise ValueError("Nonpositive representative momenta or negative pair weight")
    return positive_rule(k, k * weights)


def numeric_event(rule_a, rule_b):
    """Independent quadrature in both incident and COM outgoing angles in one collision."""
    k_a, prob_a = rule_a
    k_b, prob_b = rule_b
    z, wz = roots_legendre(6)
    u, wu = roots_legendre(6)
    # Sample incoming directions uniformly and then weight by Møller flux times sigma.
    wu = wu * (1 - u)**2
    wu /= np.sum(wu)
    weights = (prob_a[:, None, None, None] * prob_b[None, :, None, None] *
               wu[None, None, :, None] * (wz[None, None, None, :] / 2))
    x = k_a[:, None, None, None]
    y = k_b[None, :, None, None]
    momentum = np.sqrt((x - y)**2 + 2*x*y*(1 + u[None, None, :, None]))
    out_a = (x + y + momentum*z[None, None, None, :]) / 2
    out_b = (x + y - momentum*z[None, None, None, :]) / 2
    return np.array([np.sum(weights * (out_a**m + out_b**m - x**m - y**m) / 2)
                     for m in range(4)])


def load_daughter_module(path):
    spec = importlib.util.spec_from_file_location("shape_daughter_projection", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def thermal_control():
    # Maxwell--Boltzmann, T=1, one-particle density proportional to k^2 e^-k.
    moments = np.array([1., 3., 12., 60., 360.])
    nodes, weights = roots_genlaguerre(8, 3.)
    numeric = numeric_event((nodes, weights/np.sum(weights)), (nodes, weights/np.sum(weights)))
    analytic = event_change(moments, moments)
    np.testing.assert_allclose(analytic, 0., atol=1e-13)
    np.testing.assert_allclose(numeric, 0., atol=3e-12)
    mono = np.ones(5)
    mono_numeric = numeric_event((np.ones(1), np.ones(1)), (np.ones(1), np.ones(1)))
    np.testing.assert_allclose(mono_numeric, [0., 0., 1/12, 1/4], atol=3e-14)
    np.testing.assert_allclose(mono_numeric, event_change(mono, mono), atol=3e-14)
    return dict(thermal_deltas=numeric.tolist(), monoenergetic_deltas=mono_numeric.tolist())


def chemical_comparison(raw):
    """Fixed-energy comparison restricted to one massless Dirac flavor, per sign with g=2."""
    g = 2
    mean_k = raw[1] / raw[0]
    # The common scale factor drops from all ratios if the emitted gas redshifts freely.
    theta_zero_unit = (240 * raw[1] / (7*g*np.pi**2))**0.25
    n_zero_unit = 3*g*zeta(3)*theta_zero_unit**3 / (4*np.pi**2)
    ratio_unit = n_zero_unit / raw[0]
    example_y = 0.1
    n_target = example_y**2 * raw[0]
    energy_target = example_y**2 * raw[1]

    def integral(power, eta):
        value, error = quad(lambda x: x**(power-1) * expit(eta-x), 0., np.inf,
                            epsabs=3e-14, epsrel=2e-12)
        if error > max(1e-11, abs(value)*3e-10):
            raise ArithmeticError("Fermi integral quadrature")
        return value / gamma(power)

    def density_at_mean(eta):
        f3 = integral(3, eta)
        theta = mean_k*f3 / (3*integral(4, eta))
        return g*theta**3*f3/np.pi**2

    eta = brentq(lambda value: density_at_mean(value)-n_target, -40., 0., xtol=1e-13)
    f3, f4, f5, f6 = (integral(power,eta) for power in (3, 4, 5, 6))
    theta = mean_k*f3/(3*f4)
    n_check = g*theta**3*f3 / np.pi**2
    rho_check = 3*g*theta**4*f4 / np.pi**2
    np.testing.assert_allclose([n_check,rho_check], [n_target,energy_target],rtol=3e-11)
    np.testing.assert_allclose([integral(3,0.),integral(4,0.)],
                               [3*zeta(3)/4,7*np.pi**4/720],rtol=3e-13)
    return dict(spin_states_per_sign=g,
                equilibrium_zero_mu_number_over_source_times_abs_y_to_half=float(ratio_unit),
                example_y=example_y,
                zero_mu_over_emitted_number_at_example_y=float(ratio_unit/np.sqrt(example_y)),
                elastic_number_preserving_FD_example={
                    "fugacity": float(np.exp(eta)), "mu_over_temperature": float(eta),
                    "comoving_temperature_H0": float(theta),
                    "R2_R3": [float(4*f5*f3/(3*f4**2)),float(20*f6*f3**2/(9*f4**3))],
                })


def study(repo, daughter_records, module):
    audits = thermal_control()
    results = {}
    for arm in ARMS:
        records = daughter_records["arms"][arm]["bands"]
        moments = {name: moment_vector(records[name]) for name in BANDS}
        rules = {}
        band_results = {}
        for band in BANDS:
            raw = moments[band]
            mean = raw[1] / raw[0]
            a, b, _ = incoming(raw)
            q, omega, pair = module.pair_nodes(repo, arm, band, frequency_order=96, time_order=5)
            rules[band] = projected_rule(q, omega, pair)
            direct = numeric_event(rules[band], rules[band]) / mean**np.arange(4)
            analytic = event_change(raw, raw) / mean**np.arange(4)
            np.testing.assert_allclose(direct, analytic, rtol=3e-7, atol=3e-9)
            ratios = raw[2:4] / (raw[0] * mean**np.arange(2, 4))
            band_results[band] = {
                "mean_k": float(mean), "initial_R2_R3": ratios.tolist(),
                "per_mean_collision_dR2_dR3": analytic[2:].tolist(),
                "angular_control_dR0_through_dR3": direct.tolist(),
                "size_biased_relative_variance": float(b/a**2 - 1),
            }
        # The labels describe source bands of the incoming particles. After scattering the
        # participants can exchange momentum. Their energies need not stay in either band.
        mean = moments["full"][1] / moments["full"][0]
        mixture = []
        sums = np.zeros(4)
        for band_a, band_b, multiplicity in [("low", "low", 1), ("low", "high", 2), ("high", "high", 1)]:
            first, second = moments[band_a], moments[band_b]
            frac = multiplicity * first[1] * second[1] / moments["full"][1]**2
            raw = event_change(first, second) / mean**np.arange(4)
            direct = numeric_event(rules[band_a], rules[band_b]) / mean**np.arange(4)
            np.testing.assert_allclose(direct, raw, rtol=3e-7, atol=3e-9)
            sums += frac*raw
            mixture.append(dict(incident=[band_a,band_b], fraction=float(frac),
                                weighted_dR2_dR3=(frac*raw[2:]).tolist()))
        np.testing.assert_allclose(sums[2:],band_results["full"]["per_mean_collision_dR2_dR3"],
                                   rtol=2e-7,atol=1e-9)
        np.testing.assert_allclose(sum(row["fraction"] for row in mixture),1.,rtol=2e-8)
        results[arm] = dict(bands=band_results, full_event_decomposition=mixture,
                            full_band_chemical_comparison=chemical_comparison(moments["full"]))
    return dict(normalization="Comoving one-sign R_m=<k^m>/<k>^m. Slopes are initial derivatives per mean contact depth, with identical particle and antiparticle marginals.", equilibrium_R2_R3=[4/3,20/9], analytic_controls=audits, arms=results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    project = Path(__file__).resolve().parent.parent
    parser.add_argument("--daughters", type=Path, default=project/"daughters"/"DISTRIBUTIONS.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    module = load_daughter_module(project/"daughters"/"daughters.py")
    result = study(args.repo,json.loads(args.daughters.read_text()),module)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
