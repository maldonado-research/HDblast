"""Replay the nine-step bandwidth/preparation test from authenticated cell bounds."""
import csv
import hashlib
import io
import json
import sys
from fractions import Fraction as F
from pathlib import Path

PINS = {
    'FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'FUTURE_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PAST_TRANSFER.json': '0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564',
}


def run(directory):
    blobs = {}
    for name, pin in PINS.items():
        raw = (directory / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != pin:
            raise ValueError(f'input differs: {name}')
        blobs[name] = raw
    meta = json.loads(blobs['FUTURE_ENCLOSURE.json'])
    past = json.loads(blobs['PAST_TRANSFER.json'])
    ell = F(1, 10**6)
    if F(meta['state']['ell']) != ell or F(past['ell']) != ell:
        raise ValueError('packet scale')
    gain_lo, gain_hi = map(F, meta['physical_gain_interval'])
    delta_s = F(past['difference']['safe_strict_global_upper'])
    n, step = 9, F(1, 8)
    end = n * step
    dilution_floor = 1 - F(1, 10**11)
    transfer = F(3, 10**6)
    if not (gain_lo < gain_hi < 0 and dilution_floor < 1 - F(3, 4)*(ell*end)**2
            and ell * delta_s / (-gain_hi * dilution_floor) < transfer):
        raise ValueError('past-transfer enclosure')
    grid, denom = meta['cell_data']['time_grid_denominator'], meta['cell_data']['value_denominator']
    if grid != 4096 or denom != 2**112:
        raise ValueError('cell grid')
    rows = list(csv.DictReader(io.StringIO(blobs['FUTURE_CELLS.csv'].decode('ascii'))))
    adjacent = [r for r in rows if int(r['left_grid']) <= end * grid <= int(r['right_grid'])]
    if len(adjacent) != 2:
        raise ValueError('expected two closed enclosures at the registered point')
    lower = max(F(int(r['H_lo']), denom) for r in adjacent) - transfer
    upper = min(F(int(r['H_hi']), denom) for r in adjacent) + transfer
    eta = F(35, 8) * ell**4 * end**4
    curvature = F(1, 10**22)
    if not (F('0.365949') < lower <= upper and eta / (1 - eta) < curvature):
        raise ValueError('late enclosure or curvature tolerance')
    report = {'x': str(end), 'late_u_interval': [str(lower), str(upper)],
              'frozen_mode_relative_error': str(curvature), 'bounds': []}
    for error in [F(0), F(1,2000)]:
        for cutoff in [2,4]:
            loss = 2**n * error / dilution_floor
            lower_budget = (lower - loss) / ((cutoff * step)**n + 2**n * curvature)
            report['bounds'].append({'absolute_error_A': str(error), 'frequency_cap': cutoff,
                                     'budget_lower': str(lower_budget), 'decimal': float(lower_budget)})
            if error == F(1,2000) and lower_budget <= {2:F(28820), 4:F('56.29')}[cutoff]:
                raise ValueError('reported threshold')
    width, stride = F(1,64), F(9,64)
    n = 8
    left, right = n * stride, n * stride + width
    target_mean = coverage = F(0)
    for row in rows:
        overlap = min(right, F(int(row['right_grid']), grid)) - max(left, F(int(row['left_grid']), grid))
        if overlap > 0:
            low_u = F(int(row['H_lo']), denom) - transfer
            target_mean += overlap * min(low_u, dilution_floor * low_u) / width
            coverage += overlap
    eta = F(35,8) * ell**4 * right**4
    physical_curvature = F(2,10**11)
    if not (coverage == width and target_mean > F('0.51393') and
            dilution_floor < 1 - F(3,4)*(ell*right)**2 and
            1 - dilution_floor + eta / (1 - eta) < physical_curvature):
        raise ValueError('finite exposure geometry and mean')
    exposure = {'windows': [[str(k*stride), str(k*stride+width)] for k in range(n+1)],
                'last_mean_A_lower': str(target_mean), 'bounds': []}
    for error in [F(0), F(1,2000)]:
        for cutoff in [2,4]:
            lower_budget = (target_mean - 2**n * error) / ((cutoff*stride)**n + 2**n * physical_curvature)
            if error == F(1,2000) and lower_budget <= {2:F(9856), 4:F('38.50')}[cutoff]:
                raise ValueError('finite exposure threshold')
            exposure['bounds'].append({'absolute_error_mean_A': str(error), 'frequency_cap': cutoff,
                                       'budget_lower': str(lower_budget), 'decimal': float(lower_budget)})
    report['exposures'] = exposure
    return report


if __name__ == '__main__':
    if len(sys.argv) > 2:
        raise SystemExit('usage: bandwidth.py [capsule-input-directory]')
    directory = Path(sys.argv[1]) if len(sys.argv) == 2 else Path(__file__).resolve().parent / 'inputs'
    print(json.dumps(run(directory), indent=2))
