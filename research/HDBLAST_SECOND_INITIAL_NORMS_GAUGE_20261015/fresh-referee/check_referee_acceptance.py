#!/usr/bin/env python3
"""Separate exact collar and coarse initial-frame checks for the referee.

No project module is imported. This is a bounded supplement to the analytic
review, not a replacement for the full interval or constraint derivations.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

PINS = {
    'background/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json': 'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
    'background/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'background/ACTUAL_ANNULUS_STATIC_JETS.json': '3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22',
    'producer/SECOND_INITIAL_NORMS_NORMAL.json': 'e18b348a86a2702a6b520c0cc56e9205918cb725d0677eb970e9db68a129598d',
}


def need(condition, message):
    if not condition:
        raise ValueError(message)


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def times(a, b):
    values = [x*y for x in a for y in b]
    return (min(values), max(values))


def scale(c, a):
    return times((F(c), F(c)), a)


def square(a):
    lo = F(0) if a[0] <= 0 <= a[1] else min(a[0]**2, a[1]**2)
    return (lo, max(a[0]**2, a[1]**2))


def inverse(a):
    need(a[0] > 0 or a[1] < 0, 'zero denominator')
    return (1/a[1], 1/a[0])


def run(root):
    reviewed = json.loads(Path(__file__).with_name('REFEREE_ACCEPTANCE_PINS.json').read_text())
    need(reviewed['decision'] == 'ACCEPT_CONDITIONAL_INITIAL_CHECKPOINT', 'accepted referee pin declaration')
    for name, entry in reviewed['reviewed_files'].items():
        need(not name.startswith('fresh-referee/') and name not in ('MANIFEST.json', 'REPLAY_ROSTER.json', 'reproduce.py'),
             'nonrecursive referee pin set')
        raw = (root/name).read_bytes()
        need(len(raw) == entry['bytes'] and sha256(raw).hexdigest() == entry['sha256'],
             'reviewed candidate changed: '+name)
    inputs = {}
    for name, expected in PINS.items():
        raw = (root/name).read_bytes()
        need(sha256(raw).hexdigest() == expected, 'referee input changed: '+name)
        inputs[name] = json.loads(raw)
    old, shell, bg, prod = [inputs[name] for name in PINS]
    boxes = {k: tuple(map(F, old['background_intervals'][k]['exact']))
             for k in ('rho', 'eta', 'p', 'a')}
    r, e, p, a = [boxes[k] for k in ('rho', 'eta', 'p', 'a')]
    # U'(eta), evaluated independently by interval Horner arithmetic.
    u1_coefficients = [F(0), F(28, 9), F(50, 9), F(-2, 3), F(-20, 9), F(-4, 9)]
    u1 = (F(0), F(0))
    for c in reversed(u1_coefficients):
        u1 = add(times(u1, e), (c, c))
    rhs = {
        'rho': scale(-1, times(square(r), a)),
        'eta': scale(-1, times(r, p)),
        'p': scale(-1, times(r, add(u1, scale(-4, times(a, p))))),
        'a': add(inverse(r), scale(F(1, 3), times(r, square(p)))),
    }
    margins = {}
    margin_lower_bounds = {
        'rho': (F(9, 10), F(107, 100)),
        'eta': (F(27, 10**7), F(69, 10**8)),
        'p': (F(28, 10**7), F(41, 10**7)),
        'a': (F(377, 10**6), F(622, 10**6)),
    }
    length = F(8, 10**6)
    for name, shell_name in zip(('rho', 'eta', 'p', 'a'), ('rho_b', 'eta_b', 'p_b', 'a_b')):
        expected = tuple(map(F, old['background_intervals'][name+'_z']['exact']))
        need(rhs[name] == expected, 'independent outer ODE image: '+name)
        start = tuple(map(F, shell['boxes'][shell_name]['exact']))
        image = add(start, times((F(0), length), rhs[name]))
        row = (image[0]-boxes[name][0], boxes[name][1]-image[1])
        need(min(row) > 0, 'strict whole-collar first-exit inclusion: '+name)
        need(all(x > lower for x, lower in zip(row, margin_lower_bounds[name])),
             'simple exact first-exit margin lower bounds: '+name)
        margins[name] = list(map(str, row))

    def sup(name, j=0):
        return max(abs(F(x)) for x in bg['z_jets'][name][j]['exact'])

    rmin = min(map(F, bg['z_jets']['rho'][0]['exact']))
    need(rmin >= 129, 'coarse positive frame radius')
    qz = F(16*10**6)
    U = sup('C')
    Uy = (sup('C', 1)+U*qz)/rmin
    DU = 2*U
    Pi1 = sup('B')*(qz+sup('k')+1)/rmin
    need(U < F(1, 10**4) and Uy < 8 and DU < F(2, 10**4) and Pi1 < 1500,
         'independent coarse first-data majorants')
    need(sup('a') < F(14, 125), 'coarse connection coefficient')
    caps = {
        'gamma2_components': F('1.06e24'),
        'gamma2_full_spatial_covariant_components': F('4.25e29'),
        'gamma2_radial_covariant_components': F('4.25e29'),
        'gamma2_angular_tangent_covariant_components': F('1.34e22'),
        'K2_components': F('7.96e7'),
        'phi2': F('5.13e33'),
        'phi2_full_spatial_derivative': F('3.06e39'),
        'Pi2': F('2e11'),
    }
    for key, value in caps.items():
        need(F(prod['bounds'][key]['exact_upper']) < value, 'physical cap input: '+key)
    u, uy, du, pi1, amin, rlo = F('1e-4'), F(8), F('2e-4'), F(1500), F(14, 125), F(129)
    g, gy, gang, q, pp = [caps[k] for k in ('gamma2_components', 'gamma2_radial_covariant_components',
                                          'gamma2_angular_tangent_covariant_components', 'K2_components', 'Pi2')]
    frame = {
        'normal_00': 8*q+44*u*u/rlo,
        'normal_0y': 26*u*uy+F(3, 2)*gy+3*amin*g,
        'normal_0a': 22*u*du/rlo+gang/2,
        'normal_yy': 2*q+8*u*u/rlo,
        'normal_ab': 2*q+4*u*u/rlo,
        'additional_spatial_00_radial': 4*u*uy,
        'additional_spatial_00_angular': 4*u*du/rlo,
    }
    need(max(frame.values()) < F('6.4e29'), 'all added initial metric derivative classes')
    need(caps['gamma2_full_spatial_covariant_components'] < F('6.4e29'), 'inherited metric spatial derivatives')
    scalar_normal = pp+2*u*pi1
    need(scalar_normal < F('3.06e39'), 'added initial scalar time derivative')
    return {
        'status': 'PASS_INDEPENDENT_REFEREE_COLLAR_AND_COARSE_INITIAL_FRAME_CHECKS',
        'referee_decision': reviewed['decision'],
        'reviewed_files': reviewed['reviewed_files'],
        'input_sha256': PINS,
        'collar_length': str(length),
        'strict_first_exit_margins': margins,
        'strict_first_exit_margin_lower_bounds': {k: list(map(str, v)) for k, v in margin_lower_bounds.items()},
        'coarse_frame_derivative_upper_exact': {k: str(v) for k, v in frame.items()},
        'coarse_scalar_normal_upper_exact': str(scalar_normal),
        'scope': 'Exact input and coarse-bound checks supplement the analytic referee review; no evolution or exact nonlinear family.',
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', '--candidate-root', dest='candidate_root', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse to overwrite referee result')
    args.output.write_text(json.dumps(run(args.candidate_root), indent=2, sort_keys=True)+'\n')


if __name__ == '__main__':
    main()
