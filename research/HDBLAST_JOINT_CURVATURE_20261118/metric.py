"""Reproduce finite exposures of the intrinsic trace-free Ricci perturbation."""
import csv
import hashlib
import io
import json
from fractions import Fraction as F
from pathlib import Path
import sys

PINS = {
    'FUTURE_CELLS.csv': '2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
    'FUTURE_ENCLOSURE.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    'PAST_TRANSFER.json': '0955ee307c1df84f6cfd7868a0e9ef5e043a176b8285b83ca52a0c48cd9a5564',
}
ELL = F(1, 10**6)
END = F(25, 16)
FRONT = F(1)
WINDOWS = [(F(5,4), F(11,8)),(F(49,32), F(25,16))]
FLOOR = 1 - F(1,10**11)
TRANSFER = F(3,10**6)


def verify(flag, message):
    if not flag:
        raise ValueError(message)


def kernel_weight(a, b, start, stop):
    def tail(c):
        return (max(c-start,F(0))**3 - max(c-stop,F(0))**3)/6
    return (tail(b)-tail(a))/(b-a)


def calculate(directory):
    def read(name):
        raw = (directory/name).read_bytes()
        verify(hashlib.sha256(raw).hexdigest() == PINS[name], name)
        return raw
    report = json.loads(read('FUTURE_ENCLOSURE.json'))
    past = json.loads(read('PAST_TRANSFER.json'))
    grid = report['cell_data']['time_grid_denominator']
    denom = report['cell_data']['value_denominator']
    verify(grid == 4096 and denom == 2**112 and report['cell_data']['sha256'] == PINS['FUTURE_CELLS.csv'], 'grid')
    verify(F(report['state']['ell']) == F(past['ell']) == ELL, 'clock')
    gain = F(report['physical_gain_interval'][1])
    delta = F(past['difference']['safe_strict_global_upper'])
    verify(-gain>0 and delta == F(1,40) and 1-F(3,4)*(ELL*END)**2>FLOOR, 'physical enclosure')
    verify(ELL*delta/(-gain*FLOOR) < TRANSFER, 'transfer radius')
    cells = []
    edge = F(0)
    for row in csv.DictReader(io.StringIO(read('FUTURE_CELLS.csv').decode('ascii'))):
        l, r = F(int(row['left_grid']),grid), F(int(row['right_grid']),grid)
        verify(edge == l < r, 'coverage')
        edge = r
        s, t = max(l, FRONT), min(r, END)
        if s < t:
            low, high = F(int(row['H_lo']),denom)-TRANSFER, F(int(row['H_hi']),denom)+TRANSFER
            verify(low <= high,'cell order')
            cells.append((s,t,low,high))
    verify(edge == 3 and sum(t-s for s,t,*_ in cells)==END-FRONT,'covered domain')
    source_bound = max(max(abs(low),abs(high)) for _,_,low,high in cells)
    length = END - FRONT
    # 0<5/2 - 35/4 tanh(ell*x)^2 <=5/2 on this interval.
    verify(F(5,2)-F(35,4)*(ELL*END)**2>0,'time coefficient')
    theta = F(5,2)*ELL**2 * length**2/2
    z_bound = source_bound*length**2/(2*(1-theta))
    flat_error = theta*z_bound
    error = flat_error + (1-FLOOR)*z_bound
    means=[]
    pointwise=[]
    null_windows=[]
    # sech^2(ell*x) = 1 - tanh^2(ell*x) >= 1 - (ell*x)^2.
    spatial_frame_floor = 1 - (ELL*END)**2
    def point_kernel(c, s, t):
        return (max(c-s,F(0))**2-max(c-t,F(0))**2)/2
    for a,b in WINDOWS:
        low,high=F(0),F(0)
        point_low,point_high=F(0),F(0)
        for s,t,u,v in cells:
            w=kernel_weight(a,b,s,t)
            verify(w>=0,'positive weight')
            low+=u*w
            high+=v*w
            left, right = point_kernel(a,s,t), point_kernel(b,s,t)
            point_low+=min(u*left,u*right)
            point_high+=max(v*left,v*right)
        means.append([low-error,high+error])
        pointwise.append([point_low-error,point_high+error])
        verify(point_low-error>0, "pointwise positive metric")
        null_windows.append([(low-error)*spatial_frame_floor,high+error])
    verify(means[0][0]>0 and means[1][0]>0,'signed metric response')
    ratio=[means[1][0]/means[0][1],means[1][1]/means[0][0]]
    verify(all(isinstance(v,F) for group in means+pointwise+null_windows+[ratio] for v in group),
           "exact Fraction endpoints")
    noise = F(1,2000)
    full_theta = F(5,2)*ELL**2 * END**2/2
    transfer_noise = noise*END**2 / (2*FLOOR*(1-full_theta))
    verify(0<transfer_noise<F(611,10**6), 'absolute sensitivity to source error')
    return {'inputs':PINS,'normalized_geometry':'M = mu*rho_b^2*U/ell = sech(ell*x)^(3/2)*z',
            'input_error_in_A':noise,'metric_output_error_if_transfer_fixed':transfer_noise,
            'cell_count':len(cells),'source_bound':source_bound,'theta':theta,'z_bound':z_bound,
            'absolute_mean_error':error,'windows':WINDOWS,'metric_means':means,
            'uniform_metric_intervals':pointwise,'uniform_null_abs_intervals':[[a*spatial_frame_floor,b] for a,b in pointwise],
            'null_abs_exposures':null_windows,'ratio_second_over_first':ratio}

if __name__=='__main__':
    source=Path(sys.argv[1]) if len(sys.argv)==2 else Path(__file__).resolve().parent/'inputs'
    if len(sys.argv)>2: raise SystemExit('usage: metric.py [input directory]')
    print(json.dumps(calculate(source),indent=2,default=str))
