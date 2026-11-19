"""Test a scalar-dependent Yukawa coupling calibrated only on the quiet trajectory."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre

HERE = Path(__file__).resolve().parent
BUDGET = importlib.util.spec_from_file_location('budget', HERE.parent / 'solver' / 'budget.py')
budget = importlib.util.module_from_spec(BUDGET)
BUDGET.loader.exec_module(budget)
WINDOW = np.array([1.35, 2.05])
CUT = 2.10


def gate(clock):
    phase = np.clip((clock - WINDOW[0]) / np.diff(WINDOW)[0], 0., 1.)
    return np.where((phase > 0) & (phase < 1), np.sin(np.pi * phase)**4, 0.)


def min_velocity(mass):
    # A piecewise cubic Hermite trajectory has quadratic velocity on each span.
    extrema = mass.derivative(2).roots(extrapolate=False)
    times = np.r_[mass.x, extrema[np.isfinite(extrema)]]
    value = np.min(mass(times, 1) / 100)
    if value <= 0:
        raise ValueError('The scalar does not give a one-to-one clock on the archived history')
    return float(value)


def quiet_time_at_mass(quiet, values):
    """Invert the quiet spline using vectorized bisection, clipping beyond the gate."""
    left = np.full_like(values, WINDOW[0])
    right = np.full_like(values, WINDOW[1])
    for _ in range(49):
        mid = (left + right) / 2
        advance = quiet(mid) < values
        left = np.where(advance, mid, left)
        right = np.where(advance, right, mid)
    return np.where(values <= quiet(WINDOW[0]), WINDOW[0],
                    np.where(values >= quiet(WINDOW[1]), WINDOW[1], (left + right) / 2))


def free_parent(q, state, background):
    start, loga, mass = background
    a0 = np.exp(loga(start))
    inertia = a0**3 * np.hypot(q / a0, mass(start))
    initial = 1 / np.sqrt(2 * inertia)
    X, C, P = state
    n = (inertia * X + P / inertia - 1) / 2
    c = (inertia * X - P / inertia) / 2 + 1j * C
    size = len(q)

    def rhs(t, y):
        a = np.exp(loga(t))
        return np.r_[y[size:2*size] / a**3, -a**3*((q/a)**2+mass(t)**2)*y[:size], 1/a]

    sol = solve_ivp(rhs, (start, CUT), np.r_[initial, -1j*inertia*initial, 0.],
                    method='DOP853', dense_output=True, rtol=1e-11, atol=5e-14, max_step=.004)
    if not sol.success:
        raise RuntimeError(sol.message)
    u, p = sol.y[:size], sol.y[size:2*size]
    error = float(np.max(abs(u*p.conj() - u.conj()*p - 1j)))
    if error > 3e-9:
        raise ValueError('Free canonical commutator failed')
    return sol.sol, n, c, error


def pair_energy(q, weights, n, c, modes, eta, scale, gate_weight, cutoff, order):
    nodes, wg = roots_legendre(order)
    total = 0.
    for k, wq, pop, anomalous, wave in zip(q, weights, n, c, modes, strict=True):
        omega = k + (cutoff - k)*(nodes + 1)/2
        dweight = wg*(cutoff-k)/2 * (omega**2-k**2)/(8*np.pi**2)
        phase = np.exp(1j*omega[:, None]*eta)
        A = np.einsum('ij,j->i', phase, wave*gate_weight)
        B = np.einsum('ij,j->i', phase.conj(), wave*gate_weight)
        probability = pop*abs(A)**2 + (pop+1)*abs(B)**2 + 2*(anomalous*A*B).real
        if np.min(probability) < -1e-13:
            raise ValueError('Negative inclusive spectral weight')
        total += wq * np.dot(dweight*omega/scale, probability)
    return float(total)


def run_arm(native, grid, arm, quiet):
    q, weights, states, background = budget.read_native(native, grid, arm)
    start, loga, mass = background
    lower_bound = min_velocity(mass)
    field_times = np.array([brentq(lambda s: mass(s) - quiet(point), start, CUT, xtol=1e-14)
                            for point in [WINDOW[0], WINDOW.mean(), WINDOW[1]]])
    free, n, c, wronskian = free_parent(q, states['retain'], background)
    diagnostic_grid = np.linspace(start, CUT, 6001)
    field_gate = gate(quiet_time_at_mass(quiet, mass(diagnostic_grid)))
    row = {'minimum_velocity_on_full_archived_wall':lower_bound,
           'last_archived_time':float(mass.x[-1]), 'field_gate_on_peak_off_times':field_times.tolist(),
           'sampled_max_abs_gate_change':float(np.max(abs(field_gate - gate(diagnostic_grid)))),
           'wronskian_error':wronskian, 'energy':{}}
    for cutoff in (120, 240, 480):
        result = row['energy'][str(cutoff)] = {}
        for clock, interval in [('proper_time', WINDOW), ('scalar', field_times[[0,2]])]:
            values = []
            settings = [(1024, 144), (2048, 192)] if cutoff == 480 else [(512, 96), (1024, 128)]
            for steps, order in settings:
                s = np.linspace(*interval, steps+1)
                times = s if clock == 'proper_time' else quiet_time_at_mass(quiet, mass(s))
                y = free(s)
                ds = (s[-1] - s[0])/steps
                ws = np.ones(len(s))
                ws[1:-1:2] = 4
                ws[2:-1:2] = 2
                values.append(pair_energy(q, weights, n, c, y[:len(q)], y[-1].real,
                                          np.exp(loga(CUT)), gate(times)*ws*ds/3, cutoff, order))
            result[clock] = values[-1]
            result[clock+'_refinement_change'] = values[-1] - values[0]
    return row


def calculate(native):
    report = {'proper_window':WINDOW.tolist(), 'readout_time':CUT,
              'energy_units':'coefficient of y^2 H0^4 per comoving volume', 'grids':{}}
    for grid in ('coarse', 'fine'):
        quiet = budget.read_native(native, grid, 'quiet')[-1][-1]
        arms = {name:run_arm(native, grid, name, quiet) for name in ('quiet','plus','minus','far')}
        report['grids'][grid] = {'quiet_field_on_peak_off':(quiet([WINDOW[0],WINDOW.mean(),WINDOW[1]])/100+.5).tolist(), 'arms':arms}
        for arm in arms.values():
            for cutoff, series in arm['energy'].items():
                base = arms['quiet']['energy'][cutoff]['proper_time']
                for clock in ('scalar','proper_time'):
                    series[clock+'_pulse_percent_of_quiet'] = 100*(series[clock]-arms['quiet']['energy'][cutoff][clock])/base
                series['change_of_pulse_contrast_pp'] = series['scalar_pulse_percent_of_quiet'] - series['proper_time_pulse_percent_of_quiet']
        if arms['quiet']['sampled_max_abs_gate_change'] > 2e-13:
            raise ValueError('Scalar-clock quiet calibration failed')
        for series in arms['quiet']['energy'].values():
            if abs(series['scalar'] - series['proper_time']) > 2e-11:
                raise ValueError('Quiet energies differ between clocks')
    return report


def verify(actual, expected, field=""):
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            raise ValueError('Unexpected fields in output')
        for key in expected:
            verify(actual[key], expected[key], key)
    elif isinstance(expected, str):
        if actual != expected:
            raise ValueError('Metadata differs')
    else:
        # Subtracting two energies and converting to percent changes the
        # absolute roundoff scale. The pulse energy is about 22 in these units.
        absolute = 3e-6 if field.endswith(("_percent_of_quiet", "_pp")) else 2e-8
        np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=absolute)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('native', type=Path)
    parser.add_argument('--check', action='store_true')
    options = parser.parse_args()
    result = calculate(options.native)
    path = HERE/'CLOCK.json'
    if options.check:
        verify(result, json.loads(path.read_bytes()))
        budgets = json.loads((HERE.parent / 'pulses' / 'PULSES.json').read_bytes())['runs']
        for grid, data in result['grids'].items():
            flux = budgets[grid + '_120']
            quiet_energy = flux['quiet']['bare']['daughter']
            for name, arm in data['arms'].items():
                pulse = 0. if name == 'quiet' else flux['pulse_minus_quiet'][name]['bare']['daughter']
                np.testing.assert_allclose(arm['energy']['120']['proper_time'], quiet_energy + pulse,
                                           rtol=1e-8, atol=1e-8)
        print('Quiet calibration and eight two-clock wall histories reproduced.')
    else:
        path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
