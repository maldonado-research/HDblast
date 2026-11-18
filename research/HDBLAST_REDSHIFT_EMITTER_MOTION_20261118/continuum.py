"""Bound an interval of source redshifts and its local Hubble-law limit."""
import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

PARENT = Path(__file__).resolve().parent.parent / 'HDBLAST_JOINT_CURVATURE_20261118'
PINS = {
    'inputs/FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'inputs/SHELL_BOX.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'METRIC.json': 'd7ce7f5c602bd5f7df159ccdd7e2e099c03c4b8740bfb4397c178ab3f2addf42',
}
ELL = F(1, 10**6)
A, B = F(49, 32), F(99, 64)
K = 1 + F(1, 10**11)
ANGULAR = F(18, 25)


def read(name):
    data = (PARENT / name).read_bytes()
    if hashlib.sha256(data).hexdigest() != PINS[name]:
        raise ValueError(name)
    return data.decode('ascii')


def fraction_pair(pair):
    return [str(value) for value in pair]


def evaluate():
    metric = json.loads(read('METRIC.json'))
    box = json.loads(read('inputs/SHELL_BOX.json'))['boxes']
    b_lo, b_hi = map(F, box['B_b']['exact'])
    r_lo, r_hi = map(F, box['rho_b']['exact'])
    scale = (ELL * b_lo / r_hi, ELL * b_hi / r_lo)
    error = F(metric['absolute_mean_error'])
    m_lo, m_hi = map(F, metric['uniform_metric_intervals'][1])
    if not (list(map(F, metric['windows'][1])) == [A, F(25,16)] and m_lo > 0
            and 1 / (1 - (ELL * B)**2) < K):
        raise ValueError('clock or uniform metric range')

    # Integrate (A-t)^2 * u(t)/2 before the reception time using the positive kernel.
    history = [F(0), F(0)]
    edge = F(1)
    for row in csv.DictReader(read('inputs/FUTURE_CELLS.csv').splitlines()):
        left = max(F(1), F(int(row['left_grid']), 4096))
        right = min(A, F(int(row['right_grid']), 4096))
        if right <= left:
            continue
        if left != edge:
            raise ValueError('incomplete waveform')
        weight = ((A - left)**3 - (A - right)**3) / 6
        history[0] += weight * (F(int(row['H_lo']), 2**112) - F(3, 10**6))
        history[1] += weight * (F(int(row['H_hi']), 2**112) + F(3, 10**6))
        edge = right
    history[0] -= (A - 1) * error
    history[1] += (A - 1) * error
    if edge != A or history[0] <= 0:
        raise ValueError('missing or nonpositive receiver history')

    # J(x)=integral_0^x M. The published lower bound on M is positive after A.
    j_max = history[1] + (B - A) * m_hi
    p_low = 1 / (K**5 * B)
    p_high = (1 + K**2 * ELL**2 * A * B) / A
    free_low = 2 * history[0] - K * (B - A) * m_hi
    free_high = 2 * K * j_max
    if free_low <= 0:
        raise ValueError('continuous free-fall sign test')
    free = (ANGULAR * scale[0] * p_low * free_low, ANGULAR * scale[1] * p_high * free_high)
    supported_floor = -ANGULAR * scale[1] * p_high * K * (B - A) * m_hi

    # In the local limit the free-fall axis response in units of ell*B_b/rho_b
    # is 2*J(A)*ell/(cosh(q0)*sinh(q0)).
    slope = (2 * ANGULAR * scale[0] * history[0] / (K**2 * A),
             2 * ANGULAR * scale[1] * history[1] / A)
    return {
        'source_digests': PINS,
        'source_x_range': fraction_pair((A, B)),
        'receiver_history_J': fraction_pair(history),
        'free_fall_cap_derivative_entire_range': fraction_pair(free),
        'supported_cap_derivative_floor': str(supported_floor),
        'free_fall_cap_reciprocal_Hubble_limit': fraction_pair(slope),
        'display_only': {'receiver_history': list(map(float, history)),
                         'free_fall_range': list(map(float, free)),
                         'supported_floor': float(supported_floor),
                         'free_fall_local_limit': list(map(float, slope))},
    }


if __name__ == '__main__':
    print(json.dumps(evaluate(), indent=2))
