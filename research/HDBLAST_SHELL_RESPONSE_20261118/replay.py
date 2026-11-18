"""Recompute certified cut statistics with exact rational arithmetic."""

import csv
import hashlib
import io
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUTS = ROOT / 'inputs'
PINS = {
    'FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'FUTURE_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PAST_TRANSFER.json': '0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_input(name):
    data = (INPUTS / name).read_bytes()
    require(hashlib.sha256(data).hexdigest() == PINS[name], f'{name}: source identity')
    return data


def square_interval(low, high):
    """Exact image of a closed real interval under squaring."""
    require(low <= high, 'inverted enclosure')
    near = 0 if low <= 0 <= high else min(low**2, high**2)
    return near, max(low**2, high**2)


def bounds_to_json(bounds):
    return [str(value) for value in bounds]


def calculate():
    future = json.loads(read_input('FUTURE_ENCLOSURE.json'))
    past = json.loads(read_input('PAST_TRANSFER.json'))
    cell_bytes = read_input('FUTURE_CELLS.csv')
    grid = future['cell_data']
    ell = Fraction(future['state']['ell'])
    value_scale = grid['value_denominator'] * ell  # CSV stores ell*S.
    clock_scale = grid['time_grid_denominator']   # x = q/ell.
    reflection = Fraction(past['difference']['safe_strict_global_upper'])
    require(grid['sha256'] == PINS['FUTURE_CELLS.csv'] and grid['bytes'] == len(cell_bytes),
            'report to CSV binding')
    require(ell == Fraction(past['ell']) == Fraction(1, 10**6), 'packet width')
    require(clock_scale == 4096 and grid['value_denominator'] == 2**112, 'grid scales')
    require(reflection == Fraction(1, 40), 'past-minus-future bound in S units')

    horizon = Fraction(25, 16)
    window_left, window_right = Fraction(5, 4), Fraction(11, 8)
    # Isotropic moments of a uniform unit vector in four dimensions.
    y_squared = 2 * (Fraction(3, 4 * 6) - Fraction(1, 4 * 6))
    require(y_squared == Fraction(1, 6), 'round quadrupole normalization')
    integral = [Fraction(0), Fraction(0)]
    total_width = Fraction(0)
    window_width = Fraction(0)
    window_bounds = []
    used_cells = 0
    last_edge = row_count = 0
    for row in csv.DictReader(io.StringIO(cell_bytes.decode('ascii'))):
        left, right = int(row['left_grid']), int(row['right_grid'])
        require(last_edge == left < right, 'noncontiguous clock cells')
        row_count += 1
        last_edge = right
        x_left, x_right = Fraction(left, clock_scale), Fraction(right, clock_scale)
        overlap = min(horizon, x_right) - x_left
        if overlap <= 0:
            continue
        low = int(row['ell_S_lo']) / value_scale - reflection
        high = int(row['ell_S_hi']) / value_scale + reflection
        sq_low, sq_high = square_interval(low, high)
        integral[0] += overlap * sq_low
        integral[1] += overlap * sq_high
        total_width += overlap
        used_cells += 1
        window_overlap = min(window_right, x_right) - max(window_left, x_left)
        if window_overlap > 0:
            window_width += window_overlap
            window_bounds.append((low, high))
    require(row_count == grid['rows'] == 4353 and last_edge == 3 * clock_scale,
            'complete source domain')
    require(total_width == horizon and used_cells == 2305, 'observer clock coverage')
    require(window_width == window_right - window_left and len(window_bounds) == 512,
            'signed interval coverage')
    hull = min(pair[0] for pair in window_bounds), max(pair[1] for pair in window_bounds)
    power = [value * y_squared / horizon for value in integral]
    require(-46653 < hull[0] <= hull[1] < -34748, 'signed whole-window trace')
    require(35373971 < power[0] <= power[1] < 43294110, 'mean square coefficient')

    witnesses = {}
    for label, position, direction in [('central_left', Fraction(74999, 50000), -1),
                                       ('central_right', Fraction(149999, 100000), 1)]:
        witness = future['sign_witnesses'][label]
        require(Fraction(witness['x']) == position < horizon, 'signed clock witness')
        a, b = [Fraction(x) / ell for x in witness['ell_S']]
        bounds = [a - reflection, b + reflection]
        retained = past['selected_point_transfers'][label]
        require(bounds == [Fraction(x) for x in retained['past_S']], 'point transfer')
        require(bounds[1] < 0 if direction < 0 else bounds[0] > 0, 'opposite signs')
        witnesses[str(position)] = bounds_to_json(bounds)
    return {
        'packet_width': str(ell), 'horizon_in_packet_widths': str(horizon),
        'source_cells': row_count, 'integrated_cells': used_cells,
        'window_in_packet_widths': bounds_to_json((window_left, window_right)),
        'window_cells': len(window_bounds), 'past_coefficient_window_hull': bounds_to_json(hull),
        'angular_second_moment_of_Y': str(y_squared),
        'time_averaged_variance_epsilon_squared_coefficient': bounds_to_json(power),
        'outward_integer_bounds': [math.floor(power[0]), math.ceil(power[1])],
        'past_coefficient_at_central_witnesses': witnesses,
    }


if __name__ == '__main__':
    print(json.dumps(calculate(), indent=2, sort_keys=True))
