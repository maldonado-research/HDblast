"""Certified proper-clock comparisons using time-cell intervals and exact rationals."""
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PINS = {
    'FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'FUTURE_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PAST_TRANSFER.json': '0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564',
}
ELL = F(1, 10**6)
Q = F(25, 16)
COS_FLOOR = 1 - F(1, 10**11)  # Lower bound for sech^(3/2)(ell*x).
TRANSFER = F(3, 10**6)         # Absolute bound in normalized H units.
W1 = (F(5, 4), F(11, 8))
W2 = (F(49, 32), Q)


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def read(name):
    data = (ROOT / 'inputs' / name).read_bytes()
    require(hashlib.sha256(data).hexdigest() == PINS[name], name)
    return data


def product(a, b):
    corners = [x*y for x in a for y in b]
    return min(corners), max(corners)


def quotient(a, b):
    require(b[0]*b[1] > 0, 'denominator may vanish')
    return product(a, (1/b[1], 1/b[0]))


def load_cells():
    old = json.loads(read('FUTURE_ENCLOSURE.json'))
    past = json.loads(read('PAST_TRANSFER.json'))
    grid = old['cell_data']
    gain_left, gain_right = map(F, old['physical_gain_interval'])
    bound = F(past['difference']['safe_strict_global_upper'])
    require(F(old['state']['ell']) == F(past['ell']) == ELL, 'packet width')
    require(gain_left < gain_right < 0 and bound == F(1, 40), 'signed gain and transfer')
    require(1 - 3*(ELL*Q)**2/4 > COS_FLOOR, 'time-dependent readout floor')
    require(ELL*bound/(-gain_right*COS_FLOOR) < TRANSFER, 'transfer in H units')
    require(grid['sha256'] == PINS['FUTURE_CELLS.csv'], 'report/cell file identity')
    require(grid['time_grid_denominator'] == 4096 and grid['value_denominator'] == 2**112,
            'time and value grid')
    rows = []
    edge = 0
    for row in csv.DictReader(io.StringIO(read('FUTURE_CELLS.csv').decode('ascii'))):
        left, right = int(row['left_grid']), int(row['right_grid'])
        require(left == edge < right, 'contiguous positive-length cells')
        edge = right
        low, high = (F(int(row[k]), 2**112) for k in ('H_lo', 'H_hi'))
        require(low <= high, 'actual master enclosure')
        rows.append((F(left,4096), F(right,4096), low-TRANSFER, high+TRANSFER))
    require(len(rows) == grid['rows'] == 4353 and edge == 3*4096, 'complete source table')
    return rows


def overlaps(rows, start, end):
    for left, right, low, high in rows:
        lft, rgt = max(left,start), min(right,end)
        if lft < rgt:
            yield lft, rgt, low, high


def window_mean(rows, limits):
    start, end = limits
    lo = hi = width = F(0)
    count = 0
    for left, right, low, high in overlaps(rows, start, end):
        weight = right-left
        a,b = product((low,high), (COS_FLOOR,F(1)))
        lo += weight*a
        hi += weight*b
        width += weight
        count += 1
    require(width == end-start, 'mean window width')
    return (lo/width, hi/width), count


def integrated_spline(t, h):
    # Antiderivative of the unnormalized three-box convolution starting at zero.
    return sum(c * max(t-i*h,0)**3 for i,c in enumerate((1,-3,3,-1))) / 6


def mass_test(rows, start, h, error):
    den_lo = den_hi = num_lo = num_hi = weight_sum = F(0)
    count = 0
    for left,right,lo,hi in overlaps(rows, start,start+3*h):
        lo -= error
        hi += error
        require(lo > 0, 'sign on the complete spline support')
        sector = int((left-start)//h)
        require(right <= start+(sector+1)*h, 'cell crosses spline segment')
        weight = integrated_spline(right-start,h) - integrated_spline(left-start,h)
        require(weight >= 0, 'positive weight')
        den_lo += weight*lo
        den_hi += weight*hi
        weight_sum += weight
        signed_weight = (-1, 2, -1)[sector] * (right - left)
        a, b = sorted((signed_weight * lo, signed_weight * hi))
        num_lo += a
        num_hi += b
        count += 1
    require(weight_sum == h**3 and count == 384, 'spline mass or coverage')
    a,b = quotient((num_lo,num_hi),(den_lo,den_hi))
    # k=35/4 sech^2-9/4 belongs to [-9/4,13/2]; positive weight & measure.
    return a-ELL**2*F(13,2), b+ELL**2*F(9,4)


def run():
    rows = load_cells()
    a, ca = window_mean(rows,W1)
    b, cb = window_mean(rows,W2)
    require((ca,cb) == (512,128) and 0 < a[0] < a[1] < 5 and b[1] < 0,
            'exposures and signs in normalized units')
    ratio = quotient(b,a)
    require(F(-36754, 100000) < ratio[0] <= ratio[1] < F(-36459, 100000),
            'rounded ratio in the paper')
    out = {'scaled_window_means': [a,b], 'ratio': ratio, 'windows_cells': [ca,cb]}
    for key,err in [('waveform',F(0)), ('enlarged',F(1,2000))]:
        first = mass_test(rows,F(5,4),F(1,32),err)
        second = mass_test(rows,F(11,8),F(1,32),err)
        require(first[0]>second[1], 'independent constant-mass intervals overlap')
        gap = first[0] - second[1]
        required_gap = F(1374, 100) if err == 0 else F(1237, 100)
        require(gap > required_gap, 'mass spread times cancellation bound')
        out[key] = {'mass_squared_intervals':[first,second],
                    'minimum_mass_spread_times_cancellation':gap,
                    'minimum_relative_weak_residual':(first[0]-second[1])/2}
    return out


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, default=str))
