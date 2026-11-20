"""Closed-wall continuation with a fixed bulk mass and a conserved wall fluid."""
import argparse
import json
from math import exp, log, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def bisect(fn, low, high):
    left, right = fn(low), fn(high)
    if left == 0:
        return low
    if right == 0:
        return high
    if (left > 0) == (right > 0):
        raise ValueError("Root is outside bracket")
    for _ in range(70):
        middle = (low + high) / 2
        if (fn(middle) > 0) == (left > 0):
            low = middle
        else:
            high = middle
    return (low + high) / 2


def matched_tail(state, exponent, vacuum_share):
    """Return H²(y) and its minimum on [0,1], y=a*/a."""
    a, e, mass = state['a'], state['e'], state['mass']
    curvature, dark_radiation = a**-2, mass / a**4
    vacuum = vacuum_share * e
    material = e - vacuum
    n = exponent

    def h2(y):
        density = vacuum + material * y**n
        return density * (2 + density) + dark_radiation * y**4 - curvature * y**2

    def gradient(y):
        return (2*n*(1 + vacuum + material*y**n)*material*y**(n-2)
                + 4*dark_radiation*y**2 - 2*curvature)

    stationary = 1. if gradient(1.) <= 0 else bisect(gradient, 0., 1.)
    return h2, stationary


def budget(state, n):
    def minimum(share):
        h2, location = matched_tail(state, n, share)
        return h2(location)

    critical = bisect(minimum, 0., 1.)
    h2, critical_y = matched_tail(state, n, critical)
    report = dict(critical_vacuum_share=critical, loitering_a_over_anchor=1/critical_y,
                  critical_h2=h2(critical_y))
    for tag, share in [('without_vacuum', 0.), ('ten_percent_vacuum', .1)]:
        h2, location = matched_tail(state, n, share)
        if h2(location) >= 0 or h2(1) <= 0:
            raise ValueError("Reference budget does not give a future turnaround")
        y = bisect(h2, location, 1.)
        report[tag] = dict(vacuum_share=share, first_maximum_a_over_anchor=1/y,
                           first_maximum_efolds=-log(y), minimum_h2=h2(location))
    # Compute again by a dense direct maximum over the vacuum required at each y.
    # e_req=-1+sqrt(1+C); use C/(1+sqrt(1+C)) to avoid cancellation.
    maximum = 0.
    a, mass, e = state['a'], state['mass'], state['e']
    for i in range(1, 10001):
        y = i / 10001
        c = (y/a)**2 - mass*(y/a)**4
        e_req = c/(1+(1+c)**.5)
        maximum = max(maximum, (e_req/e-y**n)/(1-y**n))
    report['grid_supremum_share_check'] = maximum
    if abs(maximum-critical) > 2e-7:
        raise ValueError("Stationary bound and direct grid disagree")
    return report


def radius_target(state, n, expansion_factor):
    """Vary the anchor radius at fixed initial e and bulk mass on the large-radius branch."""
    e_final = state['e'] / expansion_factor**n
    source = e_final * (2 + e_final)
    mass = state['mass']
    discriminant = 1 - 4*mass*source
    if discriminant <= 0 or state['a']**2 <= 2*mass:
        raise ValueError("Target is outside the large-radius regime considered here")
    inverse_radius_sq = 2*source / (1+sqrt(discriminant))
    minimum = 1/(sqrt(inverse_radius_sq)*expansion_factor)
    size = minimum * expansion_factor
    residual = source + mass / size**4 - 1/size**2
    if abs(residual) > 5e-14 or minimum < state['a']:
        raise ValueError("Target radius is inconsistent with the stopping condition")
    return dict(required_anchor_radius=minimum,
                factor_over_measured_anchor=minimum/state['a'])


def evaluate(path, name='primary'):
    doc = json.loads(path.read_bytes())
    model = doc['runs'][name]
    trace = (dict(zip(doc['trace_columns'], row, strict=True)) for row in model['trace'])
    state = next(d for d in trace if d['t'] >= model['exit_complete_tau'] and d['H'] > 0
                 and d['intrinsic_random'] >= .9*d['e'])
    h2, _ = matched_tail(state, 4, 0.)
    if state['mass'] < 0 or abs(h2(1.) - state['H']**2) > 2e-7:
        raise ValueError("Bulk mass or matched Friedmann constraint failed")
    result = {'anchor':{k:state[k] for k in ('t','a','e','H','mass')},
              'intrinsic_share':state['intrinsic_random']/state['e'],
              'friedmann_h2_residual':h2(1.) - state['H']**2,
              'coupled_turnaround_a_over_anchor':model['maximum_radius']/state['a'],
              'dust':budget(state, 3), 'radiation':budget(state, 4)}
    for label, n in [('dust', 3), ('radiation', 4)]:
        result[label]['radius_targets_without_vacuum'] = {
            key:radius_target(state, n, factor)
            for key, factor in [('double', 2.), ('one_efold', exp(1.)), ('three_efolds', exp(3.))]
        }
    if not (.1 < result['dust']['critical_vacuum_share'] < result['radiation']['critical_vacuum_share']):
        raise ValueError("Dust and radiation thresholds have unexpected ordering")
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--decay', type=Path, default=HERE.parent/'coupled/DECAY.json')
    p.add_argument('--previous', type=Path, default=HERE.parent/'prior/CHARGED.json')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = {'scalar_threshold':evaluate(args.decay), 'previous_clock':evaluate(args.previous)}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n')
