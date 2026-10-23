#!/usr/bin/env python3
"""Independent geometry, interval joins and selected actual-oracle checks.

This is a bounded internal numerical review, not a reproof of every inherited
analytic premise. It authenticates all retained input bytes before importing
their code. Polygon areas and moments use vertical slices independently of the
producer's clipping/shoelace routines. Source maps use the reference Fraction
interval engine rather than the producer's integer-tick evaluator.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import stat
import sys

OLD_SHA = '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'
PREVIOUS_SHA = '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'
A, B, D = F(3, 2), F(1, 4), F(5, 4)
FH, GH = F(17387513156573), F(3851400)
SIGMA, ROBIN_M = F(203, 1000000), F(1, 2000)
E_TRACE = F(14866325944167915, 3992782)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique(pairs):
    out = {}
    for key, value in pairs:
        need(key not in out, 'DUPLICATE_JSON_KEY')
        out[key] = value
    return out


def decode(raw):
    return json.loads(raw, object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('NONFINITE_JSON')))


def authenticate(root, expected, count):
    need(root.is_dir() and not root.is_symlink(), 'SOURCE_ROOT')
    manifest_raw = (root / 'MANIFEST.json').read_bytes()
    need(sha(manifest_raw) == expected, 'SOURCE_MANIFEST_PIN')
    manifest = decode(manifest_raw)
    expected_rows = {}
    for row in manifest['files']:
        name = row['path']; pure = PurePosixPath(name)
        need(name == pure.as_posix() and name not in ('', '.', 'MANIFEST.json')
             and not pure.is_absolute() and '..' not in pure.parts and '\\' not in name
             and name not in expected_rows, 'SOURCE_MANIFEST_MEMBER')
        expected_rows[name] = row
    found, modes, total = {}, {}, 0
    for path in root.rglob('*'):
        need(not path.is_symlink() and (path.is_file() or path.is_dir()), 'SOURCE_NONREGULAR')
        name = path.relative_to(root).as_posix()
        modes[name] = stat.S_IMODE(path.stat().st_mode)
        if path.is_file():
            raw = path.read_bytes(); total += len(raw)
            found[name] = {'bytes': len(raw), 'sha256': sha(raw)}
    need(set(found) == set(expected_rows) | {'MANIFEST.json'} and len(found) == count,
         'SOURCE_EXACT_CLOSURE')
    need((root / 'MANIFEST.json').read_bytes() == manifest_raw, 'SOURCE_MANIFEST_CHANGED')
    for name, row in expected_rows.items():
        need(found[name] == {'bytes': row['bytes'], 'sha256': row['sha256']}, 'SOURCE_MEMBER_PIN')
    return {'manifest_sha256': expected, 'files_including_manifest': count,
            'total_input_bytes': total, 'complete_source_authenticated': True}, modes


def interval(row):
    need(isinstance(row, list) and len(row) == 2, 'INTERVAL_SHAPE')
    lo, hi = map(F, row)
    need(lo <= hi, 'INTERVAL_ORDER')
    return lo, hi


def plus(a, b):
    return a[0] + b[0], a[1] + b[1]


def scaled(a, c):
    return (c * a[0], c * a[1]) if c >= 0 else (c * a[1], c * a[0])


def equality(saved, expected, message):
    need(interval(saved) == expected, message)


def polynomial_primitive(coefficients, x):
    return sum((c * x ** (k + 1) / (k + 1) for k, c in enumerate(coefficients)), F(0))


def integrate_poly(coefficients, lo, hi):
    return polynomial_primitive(coefficients, hi) - polynomial_primitive(coefficients, lo)


def polygon_moments(polygon):
    """Exact area, integral s, integral u from independent vertical slices."""
    need(len(polygon) >= 3 and len(set(polygon)) == len(polygon), 'POLYGON_VERTICES')
    turns = []
    for a, b, c in zip(polygon, polygon[1:] + polygon[:1], polygon[2:] + polygon[:2]):
        turns.append((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]))
    need(all(value >= 0 for value in turns) and any(value > 0 for value in turns),
         'CONVEX_COUNTERCLOCKWISE_POLYGON')
    xs = sorted({p[0] for p in polygon})
    totals = [F(0), F(0), F(0)]
    for lo, hi in zip(xs, xs[1:]):
        mid = (lo + hi) / 2
        lines = set()
        for left, right in zip(polygon, polygon[1:] + polygon[:1]):
            if left[0] == right[0]:
                continue
            if min(left[0], right[0]) < mid < max(left[0], right[0]):
                slope = (right[1] - left[1]) / (right[0] - left[0])
                intercept = left[1] - slope * left[0]
                lines.add((intercept, slope))
        need(len(lines) == 2, 'CONVEX_VERTICAL_SLICE')
        lower, upper = sorted(lines, key=lambda line: line[0] + line[1] * mid)
        c, k = upper[0] - lower[0], upper[1] - lower[1]
        need(c + k * mid > 0, 'POSITIVE_SLICE')
        totals[0] += integrate_poly((c, k), lo, hi)
        totals[1] += integrate_poly((0, c, k), lo, hi)
        totals[2] += integrate_poly(((upper[0] ** 2 - lower[0] ** 2) / 2,
                                     upper[0] * upper[1] - lower[0] * lower[1],
                                     (upper[1] ** 2 - lower[1] ** 2) / 2), lo, hi)
    need(totals[0] > 0, 'POSITIVE_POLYGON_AREA')
    return tuple(totals)


def whole_domain_moments():
    # Independent integration on vertical slices, subtracting lower triangles.
    def below(c, k, lo, hi):
        c, k, lo, hi = map(F, (c, k, lo, hi))
        return (integrate_poly((c, k), lo, hi),
                integrate_poly((0, c, k), lo, hi),
                integrate_poly((c * c / 2, c * k, k * k / 2), lo, hi))
    upper = below(A + B, -1, 0, A)
    lower = below(-D, 1, D, A)
    plus_domain = tuple(x - y for x, y in zip(upper, lower))
    minus_domain = (D * D / 2, D ** 3 / 6, D ** 3 / 6)
    zero = (F(1, 2), F(1, 6), F(1, 6))
    return {'Dminus': tuple(x - y for x, y in zip(minus_domain, zero)),
            'strip': tuple(x - y for x, y in zip(plus_domain, minus_domain))}


def check_geometry(report):
    n = report['partition']['n']; h = A / n
    need(type(n) is int and n in (4, 8, 16) and F(report['partition']['h']) == h, 'MESH_SPECIFICATION')
    sums = {'Dminus': [F(0), F(0), F(0)], 'strip': [F(0), F(0), F(0)]}
    keys = set(); eligible = 0
    for index, row in enumerate(report['actual_bulk_source_cell_enclosures']):
        need(row['index'] == index, 'BULK_ROW_INDICES')
        poly = tuple(tuple(map(F, p)) for p in row['polygon'])
        need(all(len(p) == 2 for p in poly), 'POLYGON_POINT_DIMENSION')
        region = row['region']; need(region in sums, 'FOLDED_REGION')
        split = row['characteristic_split']
        need(split in ('eta_nonpositive', 'eta_nonnegative'), 'CHARACTERISTIC_SPLIT')
        i, j = row['mesh_index']; need(type(i) is int and type(j) is int and i >= 0 and j >= 0,
                                       'MESH_INDEX')
        key = i, j, region, split
        need(key not in keys, 'DUPLICATE_SPLIT_CELL'); keys.add(key)
        for s, u in poly:
            need(i * h <= s <= min((i + 1) * h, A)
                 and j * h <= u <= min((j + 1) * h, A + B), 'POLYGON_MESH_BOX')
            need(s >= 0 and u >= 0 and F(1) <= s + u <= A + B and s - u <= D,
                 'POLYGON_ACTIVE_CONE')
            need(s + u <= D if region == 'Dminus' else s + u >= D, 'POLYGON_FOLDED_REGION')
            need(s - u <= 0 if split == 'eta_nonpositive' else s - u >= 0,
                 'POLYGON_CHARACTERISTIC_HALF')
        area, ms, mu = polygon_moments(poly)
        need(F(row['exact_area']) == area, 'INDEPENDENT_VERTICAL_SLICE_AREA')
        for k, value in enumerate((area, ms, mu)): sums[region][k] += value
        bbox = min(s for s, u in poly), max(s for s, u in poly), min(u for s, u in poly), max(u for s, u in poly)
        need(tuple(map(F, row['bbox_s'] + row['bbox_u'])) == bbox, 'EXACT_VERTEX_BBOX')
        slo, shi, ulo, uhi = bbox
        need(0 <= slo <= shi <= A and 0 <= ulo <= uhi and shi + uhi < 8, 'WHOLE_ORACLE_DOMAIN')
        is_exterior = ulo >= shi
        need(row['bbox_exterior_error_eligible'] is is_exterior
             and row['actual_error_regime'] == ('exterior' if is_exterior else 'global'),
             'WHOLE_BBOX_ERROR_REGIME')
        eligible += is_exterior
        need(F(row['wave_kernel_weight']) == (1 if region == 'Dminus' else F(1, 2)),
             'REFLECTED_ROBIN_WEIGHT')
    totals = {k: tuple(v) for k, v in sums.items()}
    need(totals == whole_domain_moments(), 'GAP_FREE_DISJOINT_TOTALS')
    # Each positive-area cell is contained in exactly one mesh rectangle,
    # folded region and eta half-plane. Unique keys therefore imply disjoint
    # interiors. Subset membership + disjointness + exact total area proves
    # complete coverage up to measure-zero boundaries; moment totals add a
    # separate independent geometric check.
    part = report['partition']
    need(part['active_bulk_cells'] == len(keys) and part['bbox_exterior_error_eligible_cells'] == eligible,
         'PARTITION_COUNTS')
    need({k: F(v) for k, v in part['active_area'].items()} == {k: v[0] for k, v in totals.items()},
         'ACTIVE_AREA_METADATA')
    need({k: F(v) for k, v in part['zero_wedge_area'].items()} == {'Dminus': F(1, 2), 'strip': F(0)},
         'EXACT_ZERO_WEDGE_AREA')
    need(F(part['active_Robin_wave_kernel_weighted_area']) == F(5, 8), 'ACTIVE_ROBIN_WEIGHTED_AREA')
    return {'method': 'Independent vertical-slice polynomial integration; mesh/region/eta uniqueness proves disjoint interiors',
            'active_area_first_s_first_u': {k: list(map(str, v)) for k, v in totals.items()},
            'positive_bulk_cells': len(keys), 'whole_bbox_exterior_cells': eligible}


def load_reference(source):
    sys.path.insert(0, str(source / 'numeric'))
    sys.path.insert(0, str(source / 'projection'))
    import actual_first_jets as jets
    import angular_generic as angular
    need(tuple(angular.inner(s, angular.SCALARS[0]) for s in angular.SCALARS)
         == (F(1), F(0), F(0)), 'YY_CONSTANT_IS_PROBABILITY_MEAN')
    return jets


def parse_inputs(encoded, I):
    need(set(encoded) == {'m', 's', 'background'}, 'SOURCE_INPUT_GROUPS')
    return {group: {k: I(*interval(v)) if isinstance(v, list) else F(v)
                    for k, v in values.items()} for group, values in encoded.items()}


def encode_inputs(rows):
    return {group: {k: v.row() if hasattr(v, 'row') else str(v) for k, v in values.items()}
            for group, values in rows.items()}


def evaluate_reference(terms, inputs, I):
    values = {k: v for group in inputs.values() for k, v in group.items()}
    need(sum(len(group) for group in inputs.values()) == len(values), 'SOURCE_SYMBOL_OVERLAP')
    ans = I(0)
    for term in terms:
        monomial = I(1)
        for name, power in sorted(term['powers'].items()):
            monomial *= values[name] ** power
        ans += F(term['coefficient']) * monomial
    return ans


def reference_shell(inputs, data, I):
    m, s, bg = inputs['m'], inputs['s'], inputs['background']
    H, ss, aa, pp = m['A'], s['s'], bg['a'], bg['p']
    sigma1 = -2 * pp
    sigma2 = 4 * (1 + I(*map(F, data['shell']['boxes']['eta_b']['exact'])))
    affine = aa * H ** 2 / 4 - sigma1 * H * ss / 6 - sigma2 * ss ** 2 / 6
    linear = aa * H + sigma1 * ss / 3
    q00 = F(1, 6) * (-affine - linear * m['D'])
    qab = F(1, 6) * (affine - linear * m['L']) + F(-4, 9) * (-linear * m['M'])
    return q00, qab


def check_numeric(report, source, jets, oracle_spots):
    maps = decode((source / 'projection/GEOMETRIC_COEFFICIENT_MAPS.json').read_bytes())
    fterms = maps['maps']['metric_forcing']['yy'][0]
    qterms = maps['maps']['QC']['y'][0]
    I = jets.I
    data, errors, _ = jets.inputs()
    actual = report['actual_uncertainty']
    need(actual['retained_actual_master_C3_error_budget'] == errors['master_error_budget_by_scaled_total_order']
         and actual['retained_actual_gauge_C3_error_caps'] == errors['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'],
         'RETAINED_ACTUAL_C3_UNCERTAINTY')
    need(actual['exterior_improvement_only_when_bbox_u_lo_ge_s_hi'] is True
         and actual['inherited_actual_uncertainty_in_source_intervals_once'] is True
         and actual['within_box_variation_and_outward_arithmetic_in_source_intervals'] is True
         and actual['separate_point_quadrature_error_required'] is False
         and actual['all23_feedback_couplings_in_analytic_allowance'] is True, 'SOURCE_UNCERTAINTY_CONVENTIONS')
    need(actual['exact_source_map_terms_evaluated'] == {'metric_forcing_yy0': len(fterms), 'QC_y0': len(qterms)},
         'EXACT_FULL_SOURCE_MAP_SELECTION')
    totals = {'Dminus': (F(0), F(0)), 'strip': (F(0), F(0))}
    bulk = report['actual_bulk_source_cell_enclosures']
    for row in bulk:
        inputs = parse_inputs(row['actual_source_inputs'], I)
        value = evaluate_reference(fterms, inputs, I)
        equality(row['actual_Fyy_probability_mean_source_interval'], (value.lo, value.hi), 'REFERENCE_FRACTION_BULK_SOURCE')
        contribution = scaled((value.lo, value.hi), F(row['exact_area']) * F(row['wave_kernel_weight']))
        equality(row['signed_bulk_contribution_interval'], contribution, 'BULK_CELL_SIGNED_INTEGRAL')
        totals[row['region']] = plus(totals[row['region']], contribution)
    shells = report['actual_boundary_source_interval_enclosures']
    starts, stops = [], []
    gint = F(0), F(0)
    for index, row in enumerate(shells):
        need(row['index'] == index and row['r'] == ['0', '0'] and row['actual_error_regime'] == 'global',
             'BOUNDARY_ROW_DOMAIN')
        lo, hi = map(F, row['s']); need(1 <= lo < hi <= D and F(row['length']) == hi - lo,
                                       'BOUNDARY_INTERVAL_LENGTH')
        starts.append(lo); stops.append(hi)
        inputs = parse_inputs(row['actual_source_inputs'], I)
        q00, qab = reference_shell(inputs, data, I)
        qc = evaluate_reference(qterms, inputs, I)
        equality(row['one_copy_outward_QK00_0'], (q00.lo, q00.hi), 'ONE_COPY_OUTWARD_QK00')
        equality(row['one_copy_outward_QKab_0'], (qab.lo, qab.hi), 'ONE_COPY_OUTWARD_QKAB')
        equality(row['varying_metric_lowered_QC_y_0'], (qc.lo, qc.hi), 'FULL_VARYING_METRIC_QCY')
        bracket = -q00 + 3 * qab + qc
        equality(row['lorentzian_trace_plus_QC_interval'], (bracket.lo, bracket.hi), 'LORENTZIAN_SPHERE_TRACE_AND_QC')
        bg = inputs['background']
        need(F(row['ell']) == F(1, 1000000) == bg['ell'], 'ELL_NORMALIZATION')
        equality(row['rho'], (bg['rho'].lo, bg['rho'].hi), 'BOUNDARY_RHO_INPUT')
        gyy = 2 * I.cv(bg['ell']) * bg['rho'] * bracket
        equality(row['actual_Gyy_mean_interval'], (gyy.lo, gyy.hi), 'BOUNDARY_GYY_NORMALIZATION')
        gc = scaled((gyy.lo, gyy.hi), hi - lo)
        equality(row['G_integral_contribution_interval'], gc, 'BOUNDARY_G_INTEGRAL')
        equality(row['signed_negative_LG_contribution_interval'], scaled(gc, -1), 'NEGATIVE_LIFT_SIGN')
        gint = plus(gint, gc)
    need(bool(shells) and starts[0] == 1 and stops[-1] == D
         and all(x == y for x, y in zip(stops, starts[1:])), 'GAP_FREE_BOUNDARY_HISTORY')
    need(report['partition']['boundary_cells'] == len(shells)
         and F(report['partition']['boundary_positive_history_length']) == D - 1, 'BOUNDARY_PARTITION_METADATA')
    kf = plus(totals['Dminus'], totals['strip']); lg = scaled(gint, -1); ss = plus(kf, lg)
    ledger = report['signed_integral_ledger']
    for name, expected in {'bulk_Dminus_weight_1': totals['Dminus'], 'bulk_strip_weight_1_over_2': totals['strip'],
                           'KF_mean': kf, 'G_history_integral_mean': gint, 'LG_negative_history_mean': lg,
                           'S_equals_KF_plus_LG_mean': ss}.items():
        equality(ledger[name], expected, 'SIGNED_LEDGER_' + name)
    fb = report['full_feedback_ledger']; qd = A * SIGMA + 3 * ROBIN_M
    bd = (A * FH + GH) / (1 - qd); eb = F(5, 8) * SIGMA * bd
    er = ROBIN_M * (D - 1) ** 2 * bd / 2
    values = {'A_operator_bound': F(1, 10000), 'B_operator_bound': F(1, 10000),
              'D0_operator_bound': F(1, 1000000), 'sigma': SIGMA, 'M_operator_bound': ROBIN_M,
              'F_full23_H_bound': FH, 'G_full23_H_bound': GH, 'q_D': qd, 'W_X_domain_bound': bd,
              'bulk_value_allowance': eb, 'Robin_trace_value_allowance': er,
              'total_point_mean_allowance': E_TRACE}
    need(eb + er == E_TRACE and all(F(fb[k]) == v for k, v in values.items()), 'COMPLETE_FULL_STATE_FEEDBACK_LEDGER')
    ww = ss[0] - E_TRACE, ss[1] + E_TRACE
    equality(report['actual_Wyy_probability_mean_interval'], ww, 'ACTUAL_AUXILIARY_MEAN_JOIN')
    need(report['strict_positive'] is (ww[0] > 0) and report['strict_negative'] is (ww[1] < 0)
         and report['contains_zero'] is (ww[0] <= 0 <= ww[1]), 'HONEST_SIGN_DECISION')
    spots = []
    if oracle_spots:
        chosen = [next(row for row in reversed(bulk) if row['region'] == 'Dminus' and row['actual_error_regime'] == 'global'),
                  next(row for row in bulk if row['region'] == 'strip' and row['actual_error_regime'] == 'exterior')]
        for row in chosen:
            m, s, bg = jets.actual_jet(I(*map(F, row['bbox_s'])), I(*map(F, row['bbox_u'])), mode='actual')
            fresh = {'m': m, 's': s, 'background': bg}
            need(encode_inputs(fresh) == row['actual_source_inputs'], 'INDEPENDENT_ACTUAL_BULK_ORACLE_SPOT')
            spots.append({'kind': 'bulk', 'index': row['index'], 'region': row['region'],
                          'error_regime': row['actual_error_regime'], 'all_ordinary_jet_and_background_intervals_equal': True})
        row = shells[-1]
        m, s, bg = jets.actual_jet(I(*map(F, row['s'])), I(0), mode='actual')
        need(encode_inputs({'m': m, 's': s, 'background': bg}) == row['actual_source_inputs'],
             'INDEPENDENT_ACTUAL_BOUNDARY_ORACLE_SPOT')
        spots.append({'kind': 'boundary', 'index': row['index'], 'error_regime': 'global',
                      'all_ordinary_jet_and_background_intervals_equal': True})
    return {'all_bulk_reference_fraction_sources_checked': len(bulk),
            'all_boundary_reference_fraction_sources_checked': len(shells),
            'fresh_actual_oracle_spots': spots,
            'signed_integral_ledger': {k: list(map(str, interval(v))) for k, v in ledger.items()},
            'full_feedback_allowance': str(E_TRACE), 'actual_Wyy_probability_mean_interval': list(map(str, ww)),
            'contains_zero': ww[0] <= 0 <= ww[1], 'S_strict_positive': ss[0] > 0,
            'LG_strict_positive': lg[0] > 0}


def review(report, source, previous, producer, oracle_spots):
    before_old, modes_old = authenticate(source, OLD_SHA, 169)
    before_previous, modes_previous = authenticate(previous, PREVIOUS_SHA, 49)
    need(report['status'] == 'CONDITIONAL_ACTUAL_BOUNDARY_AUXILIARY_MEAN_ENCLOSURE', 'PRODUCER_STATUS')
    need(report['producer_sha256'] == sha(producer.read_bytes()), 'PRODUCER_PROGRAM_PIN')
    need(report['input_authentication_before_numeric_import']['geometric169'] == before_old
         and report['input_authentication_before_numeric_import']['forced49'] == before_previous
         and report['input_authentication_before_numeric_import'] == report['input_authentication_after_evaluation']
         and report['input_bytes_unchanged'] is True, 'SOURCE_AUTHENTICATION_REPORT')
    domain = report['domain']
    need((F(domain['a']), F(domain['b']), F(domain['d'])) == (A, B, D)
         and domain['backward_closed'] is True and domain['shell_history'] == ['0', '5/4']
         and domain['shell_history_exact_zero_interval'] == ['0', '1'], 'TARGET_CAUSAL_DOMAIN')
    geometric = check_geometry(report)
    jets = load_reference(source)  # all numerical imports follow both authentications.
    numeric = check_numeric(report, source, jets, oracle_spots)
    after_old, after_modes_old = authenticate(source, OLD_SHA, 169)
    after_previous, after_modes_previous = authenticate(previous, PREVIOUS_SHA, 49)
    need(before_old == after_old and modes_old == after_modes_old
         and before_previous == after_previous and modes_previous == after_modes_previous,
         'ALL_RETAINED_BYTES_AND_MODES_UNCHANGED')
    return {'status': 'PASS_INDEPENDENT_BOUNDARY_NUMERICAL_REVIEW', 'partition_n': report['partition']['n'],
            'producer_sha256': report['producer_sha256'], 'geometry': geometric, 'numeric': numeric,
            'retained_source_authentication': {'geometric169': before_old, 'forced49': before_previous},
            'retained_bytes_and_modes_unchanged': True,
            'scope': 'Every reported geometric and source/ledger arithmetic join plus specified fresh actual-oracle spots; no reproof of all inherited analytic premises or entire fresh oracle campaign.'}


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode, 'USE_PYTHON_I_S_B')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--previous-source-root', type=Path, required=True)
    ap.add_argument('--producer-script', type=Path, required=True)
    ap.add_argument('--report', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--oracle-spots', action='store_true')
    args = ap.parse_args()
    need(args.output.is_absolute() and not args.output.exists(), 'NEW_ABSOLUTE_REVIEW_OUTPUT')
    own_pin = sha(Path(__file__).read_bytes())
    raw = args.report.read_bytes(); report = decode(raw)
    result = review(report, args.source_root.resolve(), args.previous_source_root.resolve(), args.producer_script.resolve(), args.oracle_spots)
    need(sha(Path(__file__).read_bytes()) == own_pin and args.report.read_bytes() == raw,
         'REVIEW_PROGRAM_OR_REPORT_CHANGED_DURING_REVIEW')
    result.update(producer_report_sha256=sha(raw), reviewer_source_sha256=own_pin)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'partition_n': result['partition_n'],
                      'oracle_spots': len(result['numeric']['fresh_actual_oracle_spots']),
                      'review_sha256': sha(args.output.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit('BOUNDARY_NUMERICAL_REVIEW_FAILED:' + str(error))
