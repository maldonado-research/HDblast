"""First-collision optical depth for the stipulated opposite-charge contact law."""
import argparse
import importlib.util
import io
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.special import roots_legendre

NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
QUANTUM = Path("research/HDBLAST_QUANTUM_MEMORY_20261119/quantum/SPECTRUM.json")
TIMES = (2.10, 2.16)
MOMENTA = (1., 5., 10., 25., 45., 60.)

def column(repo, arm, wall_end):
    """Read the measured proper clock and scale factor. Return dimensionless columns."""
    seed = json.loads((repo / NATIVE / 'FORK.json').read_text())['fine']
    gate = json.loads((repo / QUANTUM).read_text())['arms'][arm]['gate_proper'][1]
    with ZipFile(repo / NATIVE / 'WALLS.zip') as z:
        raw = z.read('fine/retain/' + arm + '/wall.csv')
    wall = np.genfromtxt(io.BytesIO(raw), delimiter=',', names=True)
    s = wall['H0tau']
    assert s[0] < gate < wall_end <= s[-1]
    a = CubicHermiteSpline(s, wall['A_b'] - seed['reference_A'], wall['H_over_H0'],
                           extrapolate=False)
    knots = np.r_[gate, s[(s > gate) & (s < wall_end)], wall_end]
    nodes, weights = roots_legendre(8)
    half = np.diff(knots)[:, None] / 2
    stamps = knots[:-1, None] + half * (1 + nodes)
    def integrate(power):
        return float(np.sum(half * weights * np.exp(-power * a(stamps))))
    J = integrate(5)
    J_trap = float(trapezoid(np.exp(-5 * a(knots)), knots))
    J_simpson = float(simpson(np.exp(-5 * a(knots)), x=knots))
    assert abs(J / J_simpson - 1) < 1e-8
    assert np.min(a(stamps, 1)) > 0
    return dict(gate_close=gate, end=wall_end, a_gate=float(np.exp(a(gate))),
                a_end=float(np.exp(a(wall_end))), H_gate=float(a(gate, 1)),
                H_end=float(a(wall_end, 1)), delta_eta=integrate(1),
                J=J, J_native_trapezoid=J_trap, J_native_simpson=J_simpson,
                fractional_J_difference=J_trap / J - 1,
                log_expansion=float(a(wall_end) - a(gate)))

def load_daughters():
    path = Path(__file__).resolve().parents[1] / 'daughters/daughters.py'
    spec = importlib.util.spec_from_file_location('hdblast_kinetic_daughters', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.pair_nodes

class Spectrum:
    """One species from each pair, with isotropic rest-frame decays."""
    def __init__(self, nodes, upper_frequency):
        q, omega, dN = [np.asarray(x) for x in nodes]
        self.width = np.broadcast_to(q[:, None], omega.shape)
        self.max_k_support = float((upper_frequency + np.max(q)) / 2)
        self.kmin = (omega - self.width) / 2
        self.mean_pair_k = omega / 2
        self.weights = dN
        self.number = float(dN.sum())
        self.energy = float(np.sum(dN * self.mean_pair_k))
        self.mean_k = self.energy / self.number
        assert np.min(dN) >= 0 and np.min(self.kmin) >= 0

    def collided(self, u):
        """Number-weighted 1-exp(-u k), integrated over each daughter energy box."""
        if u == 0.:
            return 0.
        x = u * self.width
        # 1 - (1 - exp(-x))/x. Evaluate the small-x limit without cancellation.
        b = np.empty_like(x)
        small = x < 1e-3
        z = x[small]
        b[small] = z * (0.5 + z * (-1/6 + z * (1/24 - z/120)))
        b[~small] = 1 + np.expm1(-x[~small]) / x[~small]
        v = u * self.kmin
        probability = -np.expm1(-v) + np.exp(-v) * b
        return float(np.sum(self.weights * probability) / self.number)

    def quadrature_collided(self, u):
        x, w = roots_legendre(24)
        k = self.kmin[:, :, None] + self.width[:, :, None] * (1 + x) / 2
        return float(np.sum((self.weights[:, :, None] * w / 2) *
                            (-np.expm1(-u * k))) / self.number)

def contact(tracers, targets, wall):
    # Frequencies are in H0, κhat=κ H0^4. λ=y² κhat. dN is per y² per comoving H0³.
    c = (8 / 3) * targets.energy
    coefficient = c * wall['J']
    tau_bar = coefficient * tracers.mean_k
    # A strictly isotropic one-channel s-wave contact law obeys σ≤16π/s over the spectrum.
    # Take the continuous frequency supremum at the largest native parent momentum.
    # Later elastic collisions can populate higher energies. A unitary completion
    # with σ≤κs keeps the reported contact first-collision fraction an upper bound.
    s_bound = 4 * targets.max_k_support ** 2 / wall['a_gate'] ** 2
    kappa_edge = 16 * np.pi / s_bound ** 2
    def F(lam):
        return tracers.collided(lam * coefficient)
    def quantile(f):
        high = 2 / tau_bar
        while F(high) < f:
            high *= 2
        return float(brentq(lambda t: F(t) - f, 0., high, xtol=1e-14))
    half = quantile(.5)
    samples = {}
    for key, lam in [('illustrative_y01_contact_ceiling', .01 * kappa_edge),
                     ('unit_y2_contact_ceiling', kappa_edge), ('lambda_0.01', .01),
                     ('half_collision_root', half)]:
        u = lam * coefficient
        prob = F(lam)
        numerical = tracers.quadrature_collided(u)
        assert abs(prob - numerical) < 1e-12 * max(prob, 1e-6)
        assert prob <= -np.expm1(-lam * tau_bar) + 1e-13
        samples[key] = dict(lambda_y2_kappahat=lam, mean_initial_tag_depth=lam * tau_bar,
                            fraction=prob, angular_quadrature_fraction=numerical)
    gate_rate = c * tracers.mean_k / (wall['a_gate']**5 * wall['H_gate'])
    end_rate = c * tracers.mean_k / (wall['a_end']**5 * wall['H_end'])
    return dict(one_species_target_conformal_energy=targets.energy,
                one_species_tracer_number=tracers.number, tracer_mean_comoving_k=tracers.mean_k,
                mean_optical_depth_per_lambda=tau_bar, tau_per_k_per_lambda=coefficient,
                gate_mean_Gamma_over_H_per_lambda=gate_rate,
                end_mean_Gamma_over_H_per_lambda=end_rate,
                lambda_for_gate_mean_Gamma_equal_H=1 / gate_rate,
                lambda_for_mean_one_collision=1 / tau_bar, lambda_for_half_colliding=half,
                lambda_for_90pct_colliding=quantile(.9),
                s_support_upper_at_gate=s_bound, one_channel_contact_kappahat_ceiling=kappa_edge,
                y_squared_necessary_for_half_at_ceiling=half / kappa_edge,
                tagged=[dict(k=k, tau_per_lambda=coefficient*k,
                             end_Gamma_over_H_per_lambda=c*k/(wall['a_end']**5*wall['H_end']))
                        for k in MOMENTA], samples=samples)

def experiment(repo, order, time_order):
    quadrature = load_daughters()
    mu, angular_weights = roots_legendre(30)
    relative_factor = float(np.sum(angular_weights * (1 - mu) ** 2) / 2)
    assert abs(relative_factor - 4 / 3) < 5e-15
    output = dict(n_frequency=order, n_time_local=time_order, flux_angular_integral=relative_factor,
                  tracer_and_target='one fermion species incident on its opposite species',
                  lambda_definition='y^2 times kappa H0^4', arms={})
    for arm in ('quiet', 'plus', 'minus', 'far'):
        spectra = {b: Spectrum(quadrature(repo, arm, band=b, frequency_order=order,
                                           time_order=time_order), upper_frequency=(60. if b == 'low' else 120.))
                   for b in ('full', 'high', 'low')}
        high, low, full = [spectra[b] for b in ('high', 'low', 'full')]
        split_number = high.number + low.number - full.number
        split_energy = high.energy + low.energy - full.energy
        assert abs(split_number) < 1e-5 * full.number and abs(split_energy) < 1e-5 * full.energy
        data = dict(split_number_residual=split_number, split_energy_residual=split_energy, readouts={})
        for end in TIMES:
            wall = column(repo, arm, end)
            cases = {name: contact(spectra[tr], spectra[tg], wall)
                     for name, tr, tg in [('full_on_full', 'full', 'full'),
                                          ('high_on_full', 'high', 'full'),
                                          ('high_on_high', 'high', 'high'),
                                          ('low_on_full', 'low', 'full'),
                                          ('low_on_low', 'low', 'low')]}
            u = cases['full_on_full']['tau_per_k_per_lambda'] * .01
            split_fraction = ((low.number * low.collided(u) + high.number * high.collided(u))
                              / (low.number + high.number))
            residual = split_fraction - full.collided(u)
            assert abs(residual) < 1e-6 * full.collided(u)
            data['readouts'][f'{end:.2f}'] = dict(wall=wall, cases=cases,
                                                  split_fraction_residual_lambda_001=residual)
        output['arms'][arm] = data
    return output

if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--repo', type=Path, required=True)
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--frequency-order', type=int, default=96)
    cli.add_argument('--time-order', type=int, default=5)
    args = cli.parse_args()
    results = experiment(args.repo, args.frequency_order, args.time_order)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
