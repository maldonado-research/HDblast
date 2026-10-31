#!/usr/bin/env python3
"""Actual first-geodesic variation bound and complete proper-clock coefficient.

The numerical decisions are rational outward interval calculations. Analytic
family, geodesic and curvature-map premises are retained and documented.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PINS_SHA = '3bf9ae71f20d406c2620af085327c4564281231877f80e218ec0a2677674cca2'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def regular(path):
    return path.is_file() and not path.is_symlink() and all(not x.is_symlink() for x in path.absolute().parents)


def read(path):
    return json.loads(Path(path).read_bytes())


def authenticate(root, pin):
    need(root.is_dir() and not root.is_symlink(), 'CAPSULE_ROOT')
    manifest = root / 'MANIFEST.json'
    need(regular(manifest) and sha(manifest) == pin, 'CAPSULE_MANIFEST:' + root.name)
    seen = {'MANIFEST.json'}
    for row in read(manifest)['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'SAFE_UNIQUE_MEMBER')
        path = root / rel
        need(regular(path) and path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],
             'CAPSULE_MEMBER:' + root.name + '/' + str(rel))
        seen.add(str(rel))
    actual = set()
    for path in root.rglob('*'):
        need(not path.is_symlink(), 'CAPSULE_SYMLINK')
        if path.is_file():
            actual.add(str(path.relative_to(root)))
        else:
            need(path.is_dir(), 'CAPSULE_SPECIAL_FILE')
    need(actual == seen, 'CAPSULE_CLOSED_INVENTORY')
    return {'manifest_sha256': pin, 'members_including_manifest': len(seen)}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def encode(value):
    if hasattr(value, 'lo') and hasattr(value, 'hi'):
        return [str(value.lo), str(value.hi)]
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


def decimal(n, places=12):
    scale = 10**places
    sign = '-' if n < 0 else ''
    n = abs(n)
    return sign + str(n // scale) + '.' + str(n % scale).zfill(places)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--background-report', required=True, type=Path)
    parser.add_argument('--tt-report', required=True, type=Path)
    parser.add_argument('--formula-api', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve(strict=True)
    out = args.output.absolute()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT')
    need(all(not p.is_symlink() for p in out.parents), 'OUTPUT_ANCESTOR_SYMLINK')
    need(repo != out.resolve() and repo not in out.resolve().parents
         and HERE != out.resolve() and HERE not in out.resolve().parents, 'OUTPUT_OUTSIDE_SOURCE')
    pins_path = HERE / 'SOURCE_PINS.json'
    need(regular(pins_path) and sha(pins_path) == PINS_SHA, 'SOURCE_PINS')
    pins = read(pins_path)
    own_pin = sha(__file__)
    before = {label: authenticate(repo / 'research' / label, pin)
              for label, pin in sorted(pins['old_capsules'].items())}
    supplied = {'background_report': args.background_report, 'tt_report': args.tt_report,
                'formula_api': args.formula_api, 'freeze': HERE / 'FREEZE.json'}
    for key, path in supplied.items():
        need(regular(path) and sha(path) == pins['new_inputs'][key], 'NEW_INPUT:' + key)
    bg = read(args.background_report)
    tt = read(args.tt_report)
    freeze = read(HERE / 'FREEZE.json')
    need(bg['status'].startswith('PASS') or bg['status'].startswith('CERTIFIED'), 'BACKGROUND_STATUS')
    need(tt['status'] == 'CERTIFIED_CONDITIONAL_ACTUAL_COMBINED_TT_JET_ONLY', 'TT_STATUS')
    need(freeze['time_cells'] == 64 and freeze['scaled_x'] == ['0', '5/4']
         and freeze['scaled_r'] == ['1249/10000', '1251/10000']
         and freeze['no_adaptive_refinement'] is True, 'FROZEN_CAMPAIGN')
    need(tt['point'] == ['5/4', '1/8'], 'TT_EVENT')

    # Scientific imports occur only after all selected old/new inputs authenticate.
    oldcap = repo / 'research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031'
    native = load('proper_clock_actual_input_context', oldcap / 'main/producer.py')
    ctx = native.Context(repo)
    formula = load('proper_clock_curvature_map', args.formula_api)
    I = ctx.I
    def interval(row):
        return I(Q(row[0]), Q(row[1]))

    ell = Q(1, 10**6)
    tstar = ell * Q(5, 4)
    T = Q(bg['shoot']['tau_dagger'][1])
    restricted_tau = Q(1, 5500)
    family_time = Q(3, 2) * ell
    rho_lower = Q(129)
    tube_chi = max(abs(Q(x)) for x in bg['compact_tube']['rapidity'])
    need(tube_chi**2 < 1, 'TUBE_RAPIDITY_FINITE')
    tube_cosh_upper = 1 / (1 - tube_chi**2 / 2)
    restricted_t = restricted_tau * tube_cosh_upper / rho_lower
    need(T < restricted_tau and restricted_t < family_time, 'EXACT_FAMILY_CLOCK_DOMAIN')
    need(Q(bg['actual_background_local']['rho_a'][0]) > 0
         and Q(bg['actual_background_local']['rho_a'][1]) < 20, 'CENTER_MONOTONE_RADIAL_SPEED')
    rlo, rhi = Q(1249, 10000), Q(1251, 10000)
    need(20 * tstar**2 / (2 * ell) < Q(1, 10000), 'NARROW_CENTER_PATH')
    need(rlo < Q(bg['shoot']['scaled_z_i'][0])
         <= Q(bg['shoot']['scaled_z_i'][1]) <= Q(1, 8) < rhi, 'LAUNCH_INSIDE_CENTER_PATH')
    need(Q(bg['compact_tube']['scaled_initial_z'][1]) < Q(1, 2), 'EXACT_ZERO_INITIAL_COLLAR')
    c = I(1, Q(bg['boost']['cosh_chi_dagger'][1]))
    s = I(Q(bg['boost']['sinh_chi_dagger'][0]), 0)
    r = I(rlo, rhi)
    N = 64
    raw_keys = ['A', 'B', 'D', 'At', 'Ay', 'Bt', 'By', 'Dt', 'Dy']
    rows = []
    source_integral = Q(0)
    M = K = as_cap = Q(0)
    for j in range(N):
        xl, xh = Q(5, 4) * j / N, Q(5, 4) * (j + 1) / N
        m, _, b = ctx.actual(I(xl, xh), r)
        rho, aa, ay = b['rho'], b['a'], b['ay']
        SP = aa*m['B'] - (c*c*m['Dt'] - s*s*m['At'])/(2*rho) - c*s*m['Dy'] - s*s*m['By']
        SQ = c*c*(m['Bt']/rho - aa*(m['A']+m['D']) - m['Dy']/2) \
            + c*s*(m['At']/rho - 2*aa*m['B']) + s*s*m['Ay']/2
        source_cap = max(SP.mag(), SQ.mag())
        weight = ell * (xh-xl) * rho.hi
        source_integral += weight * source_cap
        mcap = max((aa*s).mag() + 2*(aa*c).mag(), 2*(aa*c).mag())
        kcap = max(2*(ay*c*s).mag(), ((2*aa*aa+ay)*c*c).mag())
        M, K, as_cap = max(M, mcap), max(K, kcap), max(as_cap, (aa*s).mag())
        rows.append({'x': [str(xl), str(xh)], 'r': r.row(),
                     'metric_m1': {key: m[key] for key in raw_keys},
                     'background': {'rho': rho, 'a': aa, 'ay': ay},
                     'SP': SP, 'SQ': SQ, 'source_cap': source_cap,
                     'proper_length_upper': weight, 'M_cell': mcap, 'K_cell': kcap})
    den = 1 - M*T - K*T*T/2
    denx = 1 - as_cap*T
    need(den > 0 and denx > 0, 'VARIATION_DENOMINATORS_POSITIVE')
    W = source_integral / den
    Y = T * W
    X = T * W / denx
    m, _, b = ctx.actual(I(Q(5, 4)), I(Q(1, 8)))
    invariants = formula.first_weyl_invariants(m, b)
    Ft, Fy = formula.time_quadratic_F(m, b), formula.radial_quadratic_F(m, b)
    cd = interval(bg['boost']['cosh_chi_dagger'])
    sd = interval(bg['boost']['sinh_chi_dagger'])
    R = b['rho'] * b['c']
    quadratic = 4 * (cd*Ft + sd*Fy) / R
    tilt = -4 * (I(-X, X)*(cd*invariants['At'] + sd*invariants['K'])
                 + I(-Y, Y)*(sd*invariants['Ay'] + cd*invariants['K'])) / R
    linear = interval(tt['enclosure']['metric_linear_second_TT_coefficient'])
    total = linear + quadratic + tilt
    need(total.hi < 0, 'FULL_PROPER_CLOCK_SIGN_INCONCLUSIVE')
    scale = 10**12
    lower_integer = (total.lo * scale).__floor__()
    upper_integer = -(-total.hi * scale).__floor__()
    outward = [decimal(lower_integer), decimal(upper_integer)]
    need(Q(outward[0]) <= total.lo <= total.hi <= Q(outward[1]) < 0, 'OUTWARD_SIGN_DISPLAY')

    ctx.finish()
    after = {label: authenticate(repo / 'research' / label, pin)
             for label, pin in sorted(pins['old_capsules'].items())}
    need(before == after and sha(__file__) == own_pin and sha(pins_path) == PINS_SHA, 'FROZEN_SOURCE_POSTCHECK')
    need(all(sha(path) == pins['new_inputs'][key] for key, path in supplied.items()), 'NEW_INPUTS_UNCHANGED')
    report = encode({
        'schema': 'HDBLAST_FULL_PROPER_CLOCK_COEFFICIENT_V1',
        'status': 'PASS_CONDITIONAL_FULL_PROPER_CLOCK_SCREEN_COEFFICIENT',
        'source_authentication': before, 'new_input_sha256': pins['new_inputs'],
        'source_pins_sha256': PINS_SHA, 'producer_source_sha256': own_pin,
        'point': {'scaled_x': '5/4', 'scaled_r': '1/8'}, 'ell': ell,
        'cells': rows, 'proper_T_upper': T, 'source_integral': source_integral,
        'M': M, 'K': K, 'as': as_cap, 'den': den, 'den_X': denx,
        'velocity_cap_W': W, 'position_X_cap': X, 'position_Y_cap': Y,
        'initial_variation_exact_zero_from_open_collar': True,
        'common_cosh_cap_to_target': c, 'common_sinh_cap_to_target': s,
        'restricted_tau_upper': restricted_tau, 'family_time_upper': family_time,
        'tube_rho_lower': rho_lower, 'tube_chi_abs_upper': tube_chi,
        'tube_cosh_upper': tube_cosh_upper, 'restricted_background_t_upper': restricted_t,
        'family_time_margin': family_time - restricted_t,
        'proper_clock_headroom': restricted_tau - T,
        'apex_metric_m2': m, 'apex_background': b,
        'point_weyl_invariants': invariants, 'Ft': Ft, 'Fy': Fy,
        'background_cosh_at_target': cd, 'background_sinh_at_target': sd,
        'boosted_quadratic': quadratic, 'full_geodesic_tilt_interval': tilt,
        'metric_linear_TT': linear, 'full_coefficient': total,
        'opposite_orientation_coefficient': -total,
        'outward_decimal': outward, 'qualitative_sign_upper': Q(outward[1]),
        'ineffective_quadratic_sign_constant': -Q(outward[1])/4,
        'actual_oracle_calls': N+1, 'frozen_campaign': freeze,
        'all_inputs_authenticated_before_and_after': True,
        'scope': {'conditional_exact_family_and_actual_first_field_premises': True,
                  'proper_time_gaussian_label_observer': True,
                  'full_radial_quadratic_and_angular_tilt_retained': True,
                  'effective_parameter_radius': False,
                  'detector_or_cosmological_origin_claim': False}})
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'outward_decimal': outward,
                      'actual_oracle_calls': N+1, 'sha256': sha(args.output)}))


if __name__ == '__main__':
    main()
