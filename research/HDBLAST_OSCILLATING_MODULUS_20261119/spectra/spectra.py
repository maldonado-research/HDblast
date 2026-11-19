"""Collisionless two-body spectrum from a decaying pressureless parent in flat FRW."""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import zeta

PI = math.pi
BOSE_COEFFICIENT = zeta(3.0) / PI**2 * (30 / PI**2)**0.75
OCCUPATION_BOUND = (8 * PI**4 / 30)**0.75 / (6 * zeta(3.0))
QUANTILES = (0.1, 0.5, 0.9)


def history(alpha, mass, initial_hubble_mass, method):
    # Reduced Planck mass is unity. One real daughter, unit branching fraction.
    width = alpha**2 * mass**3 / (32 * PI)
    gamma = width / (mass * initial_hubble_mass)
    log_b_initial = (2 / 3) * math.log(gamma)

    def scale(n):
        return math.exp(log_b_initial + n)

    def ode(n, state):
        time, radiation = state
        b = scale(n)
        matter = math.exp(-time)
        dt_dn = b**2 / math.sqrt(b * matter + radiation)
        return [dt_dn, b * matter * dt_dn]

    def event(ratio):
        def reached(n, state):
            return state[1] - ratio * scale(n) * math.exp(-state[0])
        reached.direction = 1
        reached.terminal = ratio == 9
        return reached

    solution = solve_ivp(ode, (0, -log_b_initial + 4), [0., 0.],
                         method=method, dense_output=True, max_step=0.13,
                         events=[event(1), event(9)], rtol=2e-11, atol=2e-15)
    if not solution.success or any(len(events) != 1 for events in solution.t_events):
        raise RuntimeError(solution.message)

    def state_at(n):
        time, radiation = solution.sol(n)
        b = scale(n)
        matter = math.exp(-time)
        hubble = width * math.sqrt(b * matter + radiation) / b**2
        rho_m = 3 * width**2 * matter / b**3
        occupation = 32 * PI**2 * width * rho_m / (mass**4 * hubble)
        return time, radiation, hubble, occupation

    peak_occupation = 3 * PI * alpha**2 * initial_hubble_mass
    epochs, profile = {}, []
    for label, roots in zip(('equality', 'ninety_percent'), solution.t_events):
        n_star = roots[0]
        b_star = scale(n_star)
        u_star, p_star, h_star, _ = state_at(n_star)
        number = 6 * width**2 / mass * (-math.expm1(-u_star)) / b_star**3
        energy = 3 * width**2 * p_star / b_star**4
        matched = number / (BOSE_COEFFICIENT * energy**0.75)

        def momentum_fraction(n):
            return math.exp(n - n_star)  # p/(m/2) at observation

        def cdfs(n):
            time, radiation, _, _ = state_at(n)
            return math.expm1(-time) / math.expm1(-u_star), radiation / p_star

        def normalized_moment(n, power, target):
            physical_p = mass * momentum_fraction(n) / 2
            return physical_p**power * state_at(n)[3] / (2 * PI**2 * target)

        moment_errors = [abs(quad(normalized_moment, 0, n_star, args=(power, target),
                                  epsabs=2e-10, epsrel=2e-10)[0] - 1)
                         for power, target in [(3, number), (4, energy)]]
        quantiles = {}
        for index, kind in enumerate(('number', 'energy')):
            quantiles[kind] = [momentum_fraction(brentq(lambda n: cdfs(n)[index] - q, 0, n_star))
                               for q in QUANTILES]
        epochs[label] = dict(efolds_since_match=n_star, gamma_elapsed_time=u_star,
                             H_over_gamma=h_star / width, bose_count_fraction=matched,
                             mean_p_over_half_mass=2 * energy / (number * mass),
                             momentum_quantiles_over_half_mass=quantiles,
                             moment_relative_errors=moment_errors)
        for n in np.linspace(0, n_star, 129):
            cdf_n, cdf_rho = cdfs(n)
            profile.append((label, momentum_fraction(n), state_at(n)[3], cdf_n, cdf_rho))
        if max(moment_errors) > 2e-7 or matched > OCCUPATION_BOUND * peak_occupation**0.25:
            raise ValueError('Number/energy closure or occupation bound failed')
    return dict(mass_over_reduced_Mpl=mass, width_over_initial_hubble=gamma,
                H_initial_over_mass=initial_hubble_mass, peak_occupation=peak_occupation,
                all_time_bose_count_upper=OCCUPATION_BOUND * peak_occupation**0.25,
                epochs=epochs), profile


def run(ledger_path, output, background_path):
    ledger = json.loads(ledger_path.read_text())['B']
    result, rows, independent_error = {}, [], 0.
    resolved = json.loads(background_path.read_text())['branches']
    late_anchors = {}
    for minimum in (0.5, -0.5):
        source = ledger[f'phi_m={minimum}']
        branch = {}
        result[str(minimum)] = branch
        for detuning in ('0.001', '1e-10', '1e-20'):
            scale = source['scales']['t=' + detuning]['H_top_over_Mpl']
            mass = source['m_over_H_top'] * scale
            args = source['alpha'], mass, 1 / 50
            primary, profile = history(*args, 'DOP853')
            check, _ = history(*args, 'RK45')
            for epoch, measures in primary['epochs'].items():
                for key in ('bose_count_fraction', 'mean_p_over_half_mass', 'gamma_elapsed_time', 'H_over_gamma'):
                    err = abs(check['epochs'][epoch][key] / measures[key] - 1)
                    independent_error = max(independent_error, err)
            branch[detuning] = primary
            rows.extend((minimum, detuning, *row) for row in profile)
        undamped = resolved[f'phi_m={minimum}']['flat']['off']['flat_anchor']
        scale = source['scales']['t=0.001']['H_top_over_Mpl']
        anchor_h = undamped['h'] / source['m_over_H_top']
        late_anchors[str(minimum)], _ = history(source['alpha'], source['m_over_H_top'] * scale,
                                                anchor_h, 'DOP853')
    if independent_error > 3e-7:
        raise ValueError('Two ODE integrations disagree')
    output.mkdir(parents=True, exist_ok=True)
    payload = dict(quantile_probabilities=QUANTILES, bose_coefficient=BOSE_COEFFICIENT,
                   occupation_bound_prefactor=OCCUPATION_BOUND,
                   maximum_two_solver_relative_difference=independent_error, branches=result,
                   delayed_start_undamped_flat_at_detuning_0p001=late_anchors)
    (output / 'SPECTRA.json').write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')
    with (output / 'SPECTRA.csv').open('w', newline='') as stream:
        writer = csv.writer(stream, lineterminator='\n')
        writer.writerow(('phi_min', 'detuning', 'epoch', 'p_over_half_mass', 'occupation_per_state',
                         'cumulative_number_fraction', 'cumulative_energy_fraction'))
        writer.writerows(rows)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--background', type=Path,
                        default=Path(__file__).resolve().parents[1] / 'background' / 'BACKGROUND.json')
    args = parser.parse_args()
    run(args.ledger, args.output, args.background)
