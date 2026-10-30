#!/usr/bin/env python3
"""Compare independent reconstruction enclosures and transfer the direct C1 bound."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read(path, digest):
    raw = path.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'INPUT_HASH:' + str(path))
    return json.loads(raw)


def interval(row):
    need(len(row) == 2, 'interval shape')
    lo, hi = map(F, row)
    need(lo <= hi, 'interval order')
    return lo, hi


def sqrtup(value):
    need(value >= 0, 'nonnegative norm square')
    grid = 1 << 64
    n = math.isqrt(value.numerator * grid * grid // value.denominator)
    if n*n*value.denominator < value.numerator*grid*grid:
        n += 1
    result = F(n, grid)
    need(result*result >= value, 'outward square root')
    return result


def widen(row, radius):
    need(radius >= 0, 'nonnegative error radius')
    return row[0] - radius, row[1] + radius


def inside(inner, outer):
    return outer[0] <= inner[0] <= inner[1] <= outer[1]


def encode(row):
    return list(map(str, row))


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
for name in ('resummed', 'coarse', 'direct', 'transfer'):
    ap.add_argument('--' + name, type=Path)
args = ap.parse_args()
root = args.repo_root.resolve()
ref = Path(__file__).resolve().parents[1] / 'reference'
pins = {
    'resummed': 'bc7a0ef6545ddf09e54f16a71f4ac3dc0cf8f1a0a85d37d78a44bf1845ed1f69',
    'coarse': 'a3c5d2a8c20b6b62870c2b0d969ac19a73c8168b0b75d571c8c77eea4956b103',
    'direct': '4506b4104f6d80e7b46b527b968ab7046765f6235cf605a1f4643d799f884b54',
    'transfer': '7d816712ba63e808d951bfe1fd2167f24664ac9c13e99aed13bb83303b0c7a52'}
data = {name: read(getattr(args, name) or ref / (name.upper() + '.json'), pin)
        for name, pin in pins.items()}
d, c, e, t = (data[name] for name in ('resummed', 'coarse', 'direct', 'transfer'))
need(all(z['status'] == 'PASS' for z in data.values()), 'successful reports')
need(d['point'] == c['point'] == ['9/8', '0'], 'same active shell point')
base = root / 'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
p = read(base / 'compute/response-normal/RESPONSE_PROFILES.json',
         '00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
r = read(base / 'compute/response-normal/RESPONSE_REPORT.json',
         'ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
labels, weights = p['labels'], list(map(F, p['gram']))
need(len(labels) == len(weights) == 23 and min(weights) > 0 and
     labels == d['labels'] == t['labels'] == r['labels'] and
     weights == list(map(F, r['gram'])), 'registered positive modal space')
E = F(1735*10**24)
need(F(e['direct_C1_X_error_cap']) <
     F(e['direct_C1_X_error_using_retained_contraction_cap']) < E,
     'conservative readable direct error')
need(F(e['inherited_prior_C1_X_error']) == F(t['inherited_prior_exact_C1_X_error'])
     and E < F(t['inherited_prior_exact_C1_X_error'])/12,
     'upper certificate over twelvefold tighter')
need(F(e['boundary_correction_whole_cone_X_cap']) ==
     F(t['global_retained_history']['correction_C1_X_cap']),
     'same independently checked global boundary correction')
M = [[interval(z) for z in row] for row in r['actual_coefficient_boxes']['M']]
need(len(M) == 23 and all(len(row) == 23 for row in M), 'actual M shape')
rows = []
checks = 0
need(len(d['full23_point_table']) == len(c['full23_point_table']) == 23,
     'complete reconstructed point tables')
for i, (dr, cr) in enumerate(zip(d['full23_point_table'], c['full23_point_table'])):
    need(dr['label'] == cr['label'] == labels[i] and dr['parity'] == cr['parity'],
         'modal row identity')
    factor = sqrtup(1/weights[i])
    candidate = {key: interval(dr['reconstructed_boundary_candidate'][key])
                 for key in ('value', 'dx', 'dr')}
    actual = {}
    for key in candidate:
        need(inside(candidate[key], interval(cr['reconstructed_boundary_candidate'][key])),
             'candidate outside independent coarse enclosure: ' + labels[i] + ' ' + key)
        checks += 1
        actual[key] = widen(candidate[key], factor*E*(F(1, 8) if key == 'value' else 1))
    if dr['parity'] < 0:
        need(candidate['value'] == candidate['dx'] == (0, 0), 'Dirichlet candidate trace')
        actual['value'] = actual['dx'] = (F(0), F(0))
        radial_error = None
        sharpened_radial = actual['dr']
    else:
        # U and W share the same actual M and G, so (U-W)_r=M(U-W).
        # Cauchy-Schwarz in H and the shell value error E/8 bound each row.
        rownorm2 = sum(max(abs(lo), abs(hi))**2 / weights[j]
                       for j, (lo, hi) in enumerate(M[i]))
        radial_error = sqrtup(rownorm2)*E/8
        boundary = widen(candidate['dr'], radial_error)
        sharpened_radial = (max(actual['dr'][0], boundary[0]),
                             min(actual['dr'][1], boundary[1]))
        need(sharpened_radial[0] <= sharpened_radial[1], 'compatible shared Robin error')
        need(inside(sharpened_radial,
                    interval(cr['conditional_actual_formal_shared_Robin_radial'])),
             'refined actual radial box outside independent coarse transfer: ' + labels[i])
        checks += 1
        need(dr['exact_affine_Robin_residual_by_analytic_identity'] == ['0', '0'],
             'reported analytic boundary identity')
    for key in actual:
        need(inside(actual[key], interval(cr['conditional_actual_formal_response'][key])),
             'new actual box outside independent coarse transfer: ' + labels[i] + ' ' + key)
        checks += 1
    if i in (21, 22):
        need(radial_error == 0, 'autonomous zero-M tensor radial error')
    rows.append({'label': labels[i], 'parity': dr['parity'],
                 'reconstructed_candidate': {key: encode(z) for key, z in candidate.items()},
                 'conditional_actual_formal_response': {key: encode(z) for key, z in actual.items()},
                 'shared_Robin_actual_radial': encode(sharpened_radial),
                 'shared_Robin_radial_error_abs_upper':
                     str(radial_error) if radial_error is not None else None})
primary = rows[12]
need(primary['label'] == 'w:0', 'primary registered scalar mean')
for key in ('value', 'dx'):
    need(interval(primary['conditional_actual_formal_response'][key])[0] > 0,
         'positive conditional scalar mean ' + key)
need(interval(primary['shared_Robin_actual_radial'])[0] > 0,
     'positive conditional actual inward derivative')
need(interval(d['full23_point_table'][12]['trace_correction_from_old_candidate'])[0] > 0,
     'positive actual-law reconstruction trace correction')
for key, lo, hi in [('value', F(82901*10**26), F(82907*10**26)),
                    ('dx', F(197460*10**27), F(197465*10**27))]:
    row = interval(primary['conditional_actual_formal_response'][key])
    need(lo < row[0] <= row[1] < hi, 'readable actual mean ' + key + ' bounds')
row = interval(primary['shared_Robin_actual_radial'])
need(F(215401*10**22) < row[0] <= row[1] < F(215414*10**22),
     'readable shared Robin radial bounds')
print(json.dumps({'status': 'PASS', 'point': ['9/8', '0'],
    'cross_method_interval_checks': checks, 'modal_rows': 23,
    'global_C1_X_error_readable_upper': str(E),
    'prior_upper_certificate_to_new_readable_upper_ratio':
        str(F(t['inherited_prior_exact_C1_X_error'])/E),
    'all_three_registered_scalar_mean_first_jet_signs_positive': True,
    'full23_point_table': rows,
    'boundary_identity_semantics':
        'exact integral-defined W shares fixed actual M and G with U; finite series and source bounds enclose W; zeros are analytic identities, not residual interval subtraction',
    'scope': 'conditional formal response error bound and one coordinate modal sample; no numerical C2 error, full bulk self-consistency, pointwise Einstein residual, nonlinear family or physical observable'}, sort_keys=True))
