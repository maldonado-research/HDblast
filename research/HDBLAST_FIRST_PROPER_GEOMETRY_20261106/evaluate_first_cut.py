"""Exact arithmetic for the evaluated first proper-cut bound.

This consumes a retained actual first-master certificate; it does not
rerun that certificate or prove the retained continuous/PDE/gauge premises.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def need(ok, why):
    if not ok:
        raise ValueError(why)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def evaluate(root):
    pins = {
        'HIGHER_JETS_NORMAL.json': '92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
        'OLD_FUTURE_REPORT.json': 'd2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
    }
    inputs = {}
    for name, pin in pins.items():
        path = root/'inputs'/name
        need(path.is_file() and not path.is_symlink(), 'Regular pinned input required')
        raw = path.read_bytes()
        need(sha(raw) == pin, 'Input identity differs: '+name)
        inputs[name] = json.loads(raw)
    old = inputs['OLD_FUTURE_REPORT.json']
    higher = inputs['HIGHER_JETS_NORMAL.json']
    ell = F(1, 10**6)
    need(old['status'] == higher['status'] == 'PASS', 'Retained report status differs')
    need(F(old['state']['ell']) == ell and old['state']['initial_velocity'] == 'exp(-z) Q(z/ell)', 'Registered state differs')
    need(higher['input_pins']['inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json'] == pins['OLD_FUTURE_REPORT.json'], 'Actual first-trace input join differs')
    need(F(higher['source_bound_adapter']['ell']) == ell, 'Higher certificate width differs')
    need(higher['source_bound_adapter']['domain'] == '0<=x<=3, r>=0, x+r<=8', 'Actual certificate domain differs')
    y = F(higher['exact']['boundary']['y'])
    passed = []

    def check(name, ok):
        need(ok, name)
        passed.append(name)

    check('actual_future_trace_simple_upper_bound', 0 < y < 189*ell)
    check('same_state_trace_transfer_coefficient', F(2,125) < 1)
    ypast = 190*ell
    check('past_trace_triangle_bound', 189*ell+ell == ypast)
    rho_min = F(129)
    u_raw = ypast/rho_min**2
    ustar = F(23,2*10**9)
    check('normalized_actual_shell_metric_coefficient', u_raw < ustar)
    target = 25*ell/16
    active = target-ell
    check('target_inside_retained_trace_domain', 0 < ell < target < 3*ell)
    check('first_tangent_active_duration', active == 9*ell/16)
    imax = active*ustar
    jmax = active**2*ustar/2
    nmax = 2*jmax
    wmax = 2*imax
    check('proper_clock_first_variation', imax == F(207,32*10**15))
    check('angular_velocity_first_variation', wmax == F(207,16*10**15))
    check('angular_position_first_variation', nmax == F(1863,512*10**21) and nmax < F(364,10**23))
    g1 = ustar*(2+2*target*active+4*active**2)
    area1 = ustar*(3+3*target*active+4*active**2)
    check('full_first_cut_metric_upper_bound', g1 < F(24,10**9))
    check('full_first_area_density_upper_bound', area1 < F(35,10**9))
    return {
        'schema':'HDBLAST_FIRST_PROPER_CUT_BOUND_V1',
        'status':'PASS', 'named_checks':passed, 'named_check_count':len(passed),
        'input_sha256':pins,
        'exact':{k:str(v) for k,v in {
            'ell':ell, 'future_actual_y_bound':y,
            'past_trace_simple_bound':ypast, 'rho_b_lower_bound':rho_min,
            'u_raw_upper':u_raw, 'u_star':ustar,
            'target_q':target, 'first_tangent_active_q':active,
            'I_upper':imax, 'J_upper':jmax,
            'longitudinal_t1_upper':imax, 'longitudinal_w1_upper':wmax,
            'longitudinal_n1_upper':nmax,
            'actual_G1_upper':g1, 'actual_j1_upper':area1,
            'display_G1_strict_upper':F(24,10**9),
            'display_j1_strict_upper':F(35,10**9),
        }.items()},
        'scope':'Retained actual first-trace input and exact rational transfer/constants only. '
                'Continuous first geodesic, physical pullback, state/gauge and PDE premises '
                'are separate. Bounds are expansion-point data, not uniform nonlinear '
                'P1/P2, harmonic U1, amplitude radius or interacting nonlinear zero.',
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    report=evaluate(root)
    need(not args.output.exists(), 'Fresh output required')
    with args.output.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'checks':report['named_check_count']}))


if __name__=='__main__':
    main()
