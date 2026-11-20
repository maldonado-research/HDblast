"""Bound closed-wall magnetic-to-electric conversion and read the archived wall history."""
import argparse
import json
import math
from bisect import bisect_right
from pathlib import Path

HERE = Path(__file__).resolve().parent
TARGET = 0.9
GROWTH = 2.0


def simpson(function, lo, hi, n):
    h = (hi-lo)/n
    terms = ((4 if i % 2 else 2) * function(lo+i*h) for i in range(1, n))
    return (function(lo) + function(hi) + math.fsum(terms)) * h/3


def available_phase(u0, u1, epsilon, n=256):
    """Twice conformal time, with a regular integrand up to turnaround."""
    if not 0 <= u0 <= u1 <= 1 or epsilon < 0:
        raise ValueError("The normalized radii must increase along the expanding branch.")
    def integrand(theta):
        x = math.sin(theta)**2
        return 1. if epsilon == 0. else x/math.sqrt(x*x + epsilon*(1+x))
    return 2*simpson(integrand, math.asin(u0), math.asin(u1), n)


def invert(fn, lo, hi):
    if not fn(lo) <= 0 <= fn(hi):
        raise ValueError("Root outside bracket.")
    for _ in range(65):
        mid = (lo+hi)/2
        if fn(mid) > 0:
            hi = mid
        else:
            lo = mid
    return (lo+hi)/2


def seed_bound(phase, target):
    return math.sin(max(0., math.asin(math.sqrt(target))-phase))**2


def loss_bound(phase, seed):
    return math.sin(min(math.pi/2, math.asin(math.sqrt(seed))+phase))**2


def radiation_match(row):
    a = row['a']
    energy = a**4 * row['e']
    constant = 2*energy + row['mass']
    root_function = lambda radius: radius**2 - constant - (energy / radius**2)**2
    hi = a
    while root_function(hi) < 0:
        hi *= 2
    maximum = invert(root_function, a, hi)
    u0, stop = a/maximum, 1/GROWTH
    epsilon = energy**2/maximum**6
    phase = available_phase(u0, stop, epsilon)
    seed = row['X']**2/(row['X']**2+row['B']**2)
    # A second integral uses the original Friedmann function in radius units.
    direct = lambda r: 2/math.sqrt(constant-r*r+(energy/r**2)**2)
    direct_phase = simpson(direct, a, stop*maximum, 512)
    if abs(phase-direct_phase) > 2e-10 or abs(root_function(maximum)) > 1e-10:
        raise ValueError("The independent phase quadratures or turning point disagree.")
    return dict(t=row['t'], a=a, e=row['e'], mass=row['mass'], comoving_total_radiation=energy,
                turnaround=maximum, u0=u0, epsilon=epsilon, phase_before_half_turn_radius=phase,
                radius_quadrature_difference=direct_phase-phase, initial_electric_share=seed,
                minimum_seed_electric_share=seed_bound(phase, TARGET),
                max_net_sheet_loss_before_doubling=loss_bound(phase, seed))


def archived_run(data):
    run = data['runs']['primary']
    rows = [dict(zip(data['trace_columns'], values, strict=True)) for values in run['trace']]
    t_exit = run['exit_complete_tau']
    current_on = run['threshold_tau']
    t_max = run['maximum_tau']
    anchor = next(row for row in rows if row['t'] >= t_exit)
    times = [row['t'] for row in rows]

    def interpolate(t, field):
        j = bisect_right(times, t)
        lower, upper = rows[j-1], rows[j]
        w = (t-lower['t'])/(upper['t']-lower['t'])
        return (1-w)*lower[field] + w*upper[field]

    def trace_phase(begin, end, stride=1):
        grid = [begin] + [t for t in times[::stride] if begin<t<end] + [end]
        # Integrate 1/a exactly for piecewise linear a(t).
        total = []
        for t0, t1 in zip(grid, grid[1:]):
            a, b = interpolate(t0, 'a'), interpolate(t1, 'a')
            z = (b-a)/a
            factor = math.log1p(z)/z if z else 1.
            total.append(2*(t1-t0)*factor/a)
        return math.fsum(total)

    def electric_angle(t):
        return math.atan2(-interpolate(t, 'X'), interpolate(t, 'B'))

    angle = electric_angle(anchor['t'])
    phase = trace_phase(anchor['t'], current_on)
    measured = electric_angle(current_on) - angle
    check_times = [t for t in times if anchor['t'] < t < current_on] + [current_on]
    worst_error = max(abs(electric_angle(t)-angle-trace_phase(anchor['t'], t)) for t in check_times)
    electric90 = invert(lambda t: electric_angle(t)-math.asin(math.sqrt(TARGET)),
                        anchor['t'], current_on)
    geometric = {event: dict(t=t, radius=interpolate(t, 'a'),
                              possible_growth=run['maximum_radius']/interpolate(t, 'a'))
                 for event, t in [('electric90_before_current', electric90), ('threshold', current_on)]}
    geometric['heat90'] = dict(t=run['first_ninety_intrinsic_time'],
                               radius=run['first_ninety_intrinsic_radius'],
                               possible_growth=run['maximum_radius']/run['first_ninety_intrinsic_radius'])
    geometric['exit'] = dict(t=t_exit, radius=run['exit_complete_radius'],
                             possible_growth=run['maximum_radius']/run['exit_complete_radius'])
    if worst_error > .002 or abs(math.sin(angle+phase)**2-run['threshold_electric_share']) > .001:
        raise ValueError("Neutral Maxwell prediction no longer matches the archived pre-charge history.")
    return dict(trajectory=dict(events=geometric, turnaround=run['maximum_radius'], turnaround_t=t_max,
                                first_constant_K_sample=anchor['t'], phase_at_first_sample=angle,
                                threshold_electric_share=run['threshold_electric_share'],
                                neutral_angle_prediction_at_threshold=angle+phase,
                                measured_angle_at_threshold=measured+angle, max_prethreshold_angle_error=worst_error,
                                geometric_phase_from_complete_exit_to_turn=trace_phase(t_exit, t_max),
                                geometric_phase_from_threshold_to_turn=trace_phase(current_on, t_max),
                                coarser_prethreshold_quadrature_change=trace_phase(anchor['t'], current_on, 2)-phase),
                radiation_match=radiation_match(anchor))


def coupled_match(path):
    document = json.loads(path.read_text())
    run = document['runs']['primary']
    trace = (dict(zip(document['trace_columns'],row,strict=True)) for row in run['trace'])
    row = next(row for row in trace if row['t'] >= run['exit_complete_tau'])
    match = radiation_match(row)
    return dict(radiation_match=match, actual_radius_fraction=row['a']/run['maximum_radius'],
                actual_turnaround=run['maximum_radius'], complete_exit_time=run['exit_complete_tau'])


def calculate(data):
    beta = math.asin(math.sqrt(TARGET))
    phase = available_phase(0., 1/GROWTH, 0.)
    # Closed, low-energy GR: a/A=sin(eta), with eta zero at the formal origin.
    if abs(phase-math.pi/3) > 1e-14 or abs(loss_bound(phase, 0)-.75) > 1e-14:
        raise ValueError("Low-energy analytic limit failed.")
    return dict(target_net_EM_fraction=TARGET, subsequent_radius_factor=GROWTH,
                limiting_low_energy_radiation=dict(
                    maximum_phase_before_half_turn_radius=phase, max_loss_from_pure_magnetic=loss_bound(phase, 0.),
                    max_growth_after_pure_magnetic_90=1/math.sin(beta/2),
                    max_growth_after_quarter_cycle=math.sqrt(2),
                    min_seed_electric_for_90_then_double=seed_bound(phase, TARGET)),
                **archived_run(data))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--decay', type=Path, default=HERE.parent/'prior/DECAY.json')
    parser.add_argument('--axion', type=Path, default=HERE.parent/'axion/GROWTH.json')
    parser.add_argument('--even', type=Path, default=HERE.parent/'exit/EXIT.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = calculate(json.loads(args.decay.read_text()))
    result['new_preparations'] = {'parity_even':coupled_match(args.even),'parity_odd':coupled_match(args.axion)}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
