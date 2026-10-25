#!/usr/bin/env python3
"""Join exact shell candidate intervals to their complete error ledger.

The anchored volume result is a conditional kinematical second jet, with
smooth mean histories and the parent projection/subsidiary assumptions.
This script neither integrates the bulk equations nor constructs a physical
finite-amplitude family. All decisions use exact rational arithmetic.
"""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

PINS = {
    'numeric/COMPONENT_LEDGER_TEN_INTEGRATED_NORMAL.json':
        '12dfb1247dba4ca1a9ed59f49e3c67209ff8c839d0d554180b5382d4b126664d',
    'numeric/SECOND_FEEDBACK_NORMAL.json':
        '23883937b7e5a55f1155e8b5c795644253ed5cf0abd51e0c40a09ece50220974',
    'balanced/b2/EXACT_B2_NORMAL.json':
        '3db587c30dca20f01e4964ac87f8f5ba457d55d15f7720d09eb840ed3e0e4935',
}
FIRST_PATH = ('research/HDBLAST_COMMON_GAUGE_FORCING_20261016/'
              'shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json')
FIRST_PIN = '559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def integer_outward(interval):
    lo, hi = interval
    return [lo.__floor__(), hi.__ceil__()]


def first_quadratic_bound(caps):
    expected = {'X': F(2987, 2000000000), 'Z': F(1119, 100000000000),
                'V_gauge': F(11117, 2500000000000000),
                'd': F(33351, 10000000000000000000000),
                'Z_t': F(28811, 1000000), 'd_t': F(11117, 5000000000000000)}
    for name, value in expected.items():
        need(F(caps[name]['strict_readable_cap']['exact']) == value, 'FIRST_CAP:' + name)
        need(F(caps[name]['exact_upper']) < value, 'FIRST_STRICT_CAP:' + name)
    ell, x, delta = F(1, 10**6), F(9, 8), F(1, 8)
    X, Z, V, d, Zx, dx = [expected[k] for k in ('X', 'Z', 'V_gauge', 'd', 'Z_t', 'd_t')]
    Zx, dx = ell*Zx, ell*dx
    cbar, a, t = delta*X, ell**2*x, ell*x
    shift, spatial = V+ell*cbar, Z+a*cbar
    groups = {
        'first_clock_history': a*(delta*X**2/2 + 4*delta*shift**2 + cbar*X),
        'first_evaluation': cbar*(3*Zx+6*a*Z+8*dx+16*a*d)/3,
        'background_evaluation': ell**2*cbar**2*(1+t**2)/2,
        'shift': 8*ell*cbar*V/3,
        'clock_gradient': 4*ell**2*cbar**2/3,
        'spatial_determinant': spatial**2/2+F(8, 3)*spatial*d+F(16, 3)*d**2,
    }
    total = sum(groups.values(), F(0))
    need(total == F(6686222274308036835249874398211279,
                   1228800000000000000000000000000000000000000000000),
         'INDEPENDENT_FIRST_QUADRATIC_SUM')
    need(total < F(55, 10**16), 'QUADRATIC_READABLE_BOUND')
    return groups, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--capsule-root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    root = args.capsule_root.resolve()
    need(root not in args.output.resolve().parents, 'OUTPUT_INSIDE_CAPSULE')
    data = {}
    for path, pin in PINS.items():
        member = root/path
        need(member.is_file() and not member.is_symlink() and sha(member) == pin,
             'JOIN_INPUT_HASH:' + path)
        data[path] = json.loads(member.read_text())
        need(data[path]['complete_closures_authenticated_before_and_after'] is True,
             'INPUT_AUTHENTICATION:' + path)
    component, error, candidate = [data[k] for k in PINS]
    need(component['authentication'] == error['authentication'] == candidate['authentication'],
         'SAME_SEVEN_FROZEN_CLOSURES')
    need(error['no_old_u_du_p1_defect'] is True, 'CANDIDATE_SCOPE_EXACT_U')
    need(candidate['arithmetic']['fixed_subintervals'] == 256, 'FIXED_BUDGET')
    need(candidate['mean_order'] == component['mean_order'] == ['K0', 'L0', 'N0', 'J0', 'w0'],
         'MEAN_COMPONENT_ORDER')
    samples = []
    for x in (F(9, 8), F(5, 4)):
        cv = next(r for r in candidate['samples'] if F(r['x']) == x)
        er = next(r for r in error['results'] if F(r['time']) == x)
        need(list(map(F, cv['evaluation']['three_term'][1])) == [0, 0], 'DIRICHLET_SHELL')
        means = {}
        for i, name in enumerate(candidate['mean_order']):
            lo, hi = map(F, cv['evaluation']['three_term'][i])
            radius = F(er['complete_shell_value_error_upper'][i])
            need(lo <= hi and radius >= 0, 'ORDERED_INTERVAL_AND_NONNEGATIVE_ERROR')
            full = [lo-radius, hi+radius]
            means[name] = {'candidate_interval_including_quadrature': [lo, hi],
                           'actual_response_error_upper': radius,
                           'complete_actual_response_interval': full,
                           'integer_outward_interval': integer_outward(full),
                           'strictly_negative': full[1] < 0,
                           'strictly_positive': full[0] > 0}
        old = next(r for r in component['shell_sample_intervals'] if F(r['x']) == x)
        old_n = list(map(F, old['traces']['N0']['conditional_response_value_interval']))
        new_n = means['N0']['complete_actual_response_interval']
        n_intersection = [max(old_n[0], new_n[0]), min(old_n[1], new_n[1])]
        need(n_intersection[0] <= n_intersection[1] < 0, 'COMPATIBLE_TWO_N_ENCLOSURES')
        need(means['J0']['strictly_negative'], 'COMPLETE_J_NEGATIVE')
        generic_j_radius = F(er['generic_value_derivative_error'][3])
        j_candidate = means['J0']['candidate_interval_including_quadrature']
        generic_j = [j_candidate[0]-generic_j_radius, j_candidate[1]+generic_j_radius]
        need(generic_j_radius >= 0 and generic_j[1] < 0, 'COARSER_GENERIC_J_NEGATIVE')
        samples.append({'x': x, 'r': F(0), 'means': means,
                        'N_interval_intersecting_two_valid_bounds': n_intersection,
                        'coarser_generic_J_interval': generic_j,
                        'coarser_generic_J_integer_outward_interval': integer_outward(generic_j)})

    first_path = args.repo_root.resolve()/FIRST_PATH
    need(first_path.is_file() and not first_path.is_symlink() and sha(first_path) == FIRST_PIN,
         'PRIMARY_FIRST_TRACE_HASH')
    first = json.loads(first_path.read_text())
    groups, qfirst = first_quadratic_bound(first['bounds']['actual_field_amplitude_caps'])
    ten = next(r for r in component['value_derivative_ten_state_comparison'] if F(r['time']) == F(9, 8))
    need(ten['exact_positive_inverse_and_identity_verified'] is True, 'POSITIVE_N_CAP_PREMISE')
    ncap = F(ten['full_response_value_derivative_upper'][2])
    need(0 <= ncap < F(1213)*10**19, 'DIRECT_N_COMPONENT_VALUE_CAP')
    ell, x, delta = F(1, 10**6), F(9, 8), F(1, 8)
    nclock = F(3, 2)*ell**2*x*delta*ncap
    need(nclock < 256*10**7, 'CLOCK_READABLE_CAP')
    jlo, jhi = samples[0]['means']['J0']['complete_actual_response_interval']
    volume = [F(3, 2)*jlo-nclock-qfirst, F(3, 2)*jhi+nclock+qfirst]
    need(volume[1] < 0, 'CONDITIONAL_ANCHORED_VOLUME_JET_NEGATIVE')
    result = {
        'schema': 'HDBLAST_SHELL_PRECISION_JOIN_V1',
        'status': 'PASS_CONDITIONAL_J_SAMPLE_AND_ANCHORED_VOLUME_JET_SIGNS',
        'input_sha256': {**PINS, FIRST_PATH: FIRST_PIN},
        'join_source_sha256': sha(Path(__file__)),
        'arithmetic': 'Exact Fraction decisions; candidate intervals already include quadrature, charged once.',
        'shell_samples': samples,
        'anchored_volume_second_jet_at_x_9_8': {
            'physical_proper_time': 'tau=rho_b*ell*9/8; background-normalized clock held fixed in epsilon',
            'formula': 'v2=V2/V0=(3/2)[J0+ell*tanh(ell*x)*integral_0^x N0 da]+Q_first',
            'direct_N_value_upper_on_history': ncap,
            'entire_N_clock_term_upper': nclock,
            'first_quadratic_group_uppers': groups,
            'first_quadratic_total_upper': qfirst,
            'v2_interval': volume,
            'v2_integer_outward_interval': integer_outward(volume),
            'volume_radius_second_jet_R2_over_R0_interval': [v/3 for v in volume],
            'strictly_negative': True,
        },
        'premises': [
            'All inherited regular-branch, first-jet, mean projection and subsidiary assumptions remain conditional.',
            'Exact smooth mean histories and the known smooth first-order induced tensor on a finite open shell tube.',
            'The physical initial anchor and future orientation are fixed and carried under gauge changes.',
            'The exact shell history is zero for x<=1; partial_x=ell*partial_t.',
            'Missing second-order angular tensors may be completed kinematically; volume depends only on the stated second-order means.',
        ],
        'limits': [
            'Coordinate J and N sample signs are representation dependent.',
            'The volume sign concerns a parameter second derivative of an anchored kinematical observable.',
            'No finite-amplitude error/remainder or epsilon range is supplied.',
            'No full angular second-order Einstein solution or exact nonlinear physical family is constructed.',
            'No measured expansion rate, physical force or evidence for a higher-dimensional Big Bang origin follows.',
        ],
    }
    with args.output.open('x') as handle:
        json.dump(encode(result), handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'output': str(args.output), 'sha256': sha(args.output),
                      'J_intervals': {str(r['x']): r['means']['J0']['integer_outward_interval'] for r in samples},
                      'v2_interval': integer_outward(volume)}, sort_keys=True))


if __name__ == '__main__':
    main()
