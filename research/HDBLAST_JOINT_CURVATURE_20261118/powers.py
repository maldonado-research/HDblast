"""Finite-exposure bounds for angular and temporal gradient powers on actual shell cuts."""
import argparse
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path

PINS = {
    'FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'FUTURE_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PAST_TRANSFER.json': '0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564',
}


def calculate(source):
    data = {}
    for name, checksum in PINS.items():
        raw = (source / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != checksum:
            raise ValueError(name)
        data[name] = raw
    meta = json.loads(data['FUTURE_ENCLOSURE.json'])
    transfer = json.loads(data['PAST_TRANSFER.json'])
    grid = meta['cell_data']['time_grid_denominator']
    den = meta['cell_data']['value_denominator']
    ell = F(meta['state']['ell'])
    if (grid, den, ell) != (4096, 2**112, F(1, 1000000)):
        raise ValueError('registered grid')
    gain = tuple(map(F, meta['physical_gain_interval']))
    floor = 1 - F(1,10**11)
    radius = F(3,10**6)
    if not gain[0] < gain[1] < 0:
        raise ValueError('gain sign')
    if not (F(transfer['difference']['safe_strict_global_upper']) == F(1,40)
            and floor < 1 - F(3,4)*(ell*F(25,16))**2
            and ell*F(1,40)/(-gain[1]*floor) < radius):
        raise ValueError('past transfer and dilution')
    first, last = 5120, 5632  # [5/4,11/8].
    energy_square_lo = F(0)
    energy_square_hi = F(0)
    knots = range(first, last+1, 32)
    bounds = {i: [] for i in knots}
    count = 0
    for row in csv.DictReader(io.StringIO(data['FUTURE_CELLS.csv'].decode('ascii'))):
        a,b = map(int, (row['left_grid'],row['right_grid']))
        if b < first or a > last:
            continue
        lo,hi = F(int(row['H_lo']),den)-radius, F(int(row['H_hi']),den)+radius
        if not 0 < lo <= hi:
            raise ValueError('positive shell window')
        lo *= floor  # A=sech^(3/2)u.
        for edge in (a,b):
            if edge in bounds:
                bounds[edge].append((lo,hi))
        overlap = min(b,last) - max(a,first)
        if overlap > 0:
            if a != first+count or b != a+1:
                raise ValueError('continuous cells')
            count += 1
            energy_square_lo += lo**2 * F(overlap,last-first)
            energy_square_hi += hi**2 * F(overlap,last-first)
    if count != 512 or any(len(pair)!=2 for pair in bounds.values()):
        raise ValueError('closed endpoint enclosures')
    nodes = [(max(t[0] for t in bounds[i]), min(t[1] for t in bounds[i])) for i in knots]
    if any(a>b for a,b in nodes):
        raise ValueError('overlap of endpoint enclosures')
    slope_power = F(0)
    for left, right in zip(nodes,nodes[1:]):
        gap = max(F(0), right[0]-left[1], left[0]-right[1])
        slope_power += gap**2 / F(32,grid) / F(last-first,grid)

    # sech^2 >= sech^(3/2)^2 >= floor^2, using 0<sech<=1.
    mean_space_lo = F(16,9) / (131*ell)**2 * gain[1]**2 * floor**2 * energy_square_lo
    mean_space_hi = F(16,9) / (129*ell)**2 * gain[0]**2 * energy_square_hi
    mean_time_lo = gain[1]**2*slope_power/(3*131**2*ell**4)
    ratio_lo = F(3,16)*slope_power/(ell**2*energy_square_hi)
    mean_residual_hi = F(3,2)*mean_space_hi-mean_time_lo
    if not (ratio_lo>10**11 and mean_residual_hi<0):
        raise ValueError('gradient separation')
    results = {
        'window_x': ['5/4','11/8'], 'partition_intervals': len(nodes)-1,
        'input_sha256': PINS,
        'mean_A_squared': [str(energy_square_lo),str(energy_square_hi)],
        'mean_dA_dx_squared_lower': str(slope_power),
        'mean_spatial_projector_second_variation': [str(mean_space_lo),str(mean_space_hi)],
        'mean_temporal_projector_second_variation_lower': str(mean_time_lo),
        'ratio_temporal_to_spatial_second_variation_lower': str(ratio_lo),
        'mean_unweighted_residual_second_variation_upper': str(mean_residual_hi),
    }
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', type=Path, default=Path(__file__).resolve().parent/'inputs')
    print(json.dumps(calculate(parser.parse_args().source),indent=2))
