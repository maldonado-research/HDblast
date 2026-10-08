"""Exact outer review of fixed-probe U/W rectangles and separate error budgets.

This checks serialized intervals and registered bookkeeping; the frozen proof
and independent review establish scientific coverage. It does not infer a
continuum error or an original stored binary80-state error from sampled values.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json
import re
import sys

sys.dont_write_bytecode = True
from registration_guard import MOMENTA, OUTPUT_SCOPE, SOURCES, load, rational, require

GATE = Q(1, 10**20)
KEYS = {(source, k) for source in SOURCES for k in MOMENTA}
CONFIGURATIONS = {'primary': 'PRIMARY_ARB512_SOURCE112_PHASE128', 'independent': 'INDEPENDENT_DYADIC512_SOURCE112_ODE160'}

def geometric_parents():
    left = Q(1, 128)
    result = []
    while left < Q(1, 2):
        right = min(3 * left / 2, Q(1, 2))
        width = right - left
        count_ratio = width * 64
        children = (count_ratio.numerator + count_ratio.denominator - 1) // count_ratio.denominator
        cx = (left + right) / 2
        radius = left / 2
        rho = 1 - cx + radius
        separation = 5 - cx - radius
        D = 1 - rho**2
        mb = 2 * (4 / separation**2 + 4 * rho / (separation * D**2) + 2 / D**2 + 8 * rho**2 / D**3 + 4 * rho**2 / D**4)
        mz = rho * mb + 2 * (2 / separation + 4 * rho / D**2)
        result.append({'center': (left + right) / 2 - 5, 'half': width / 2, 'radius': left / 2,
                       'interval': (left - 5, right - 5), 'children': children,
                       'majorants': {'positive_B': mb, 'signed_uB': mz}})
        left = right
    require(len(result) == 11, 'Fixed geometry incomplete')
    return result

def source_models(data, budget, route, fabricated):
    rows = data.get('source_model_rows') if route == 'primary' else budget.get('source_model_rows')
    require(type(rows) is list and len(rows) == 22, 'Exactly 22 source-parent model rows required')
    geometries = geometric_parents()
    seen = set()
    model_radii = {source: {'U': Q(0), 'W': Q(0)} for source in SOURCES}
    for row in rows:
        require(type(row) is dict and row.get('source') in SOURCES and type(row.get('parent')) is int and 0 <= row['parent'] < 11, 'Invalid source-parent model key')
        key = (row['source'], row['parent'])
        require(key not in seen, 'Duplicate source-parent model evidence')
        seen.add(key)
        geom = geometries[row['parent']]
        M = geom['majorants'][row['source']]
        tail = M / 2**112
        require(rational(row.get('center')) == geom['center'] and rational(row.get('half_width')) == geom['half'], 'Source-parent geometry differs')
        require(type(row.get('children')) is int and row['children'] == geom['children'], 'Source-parent child partition differs')
        if route == 'primary':
            require(rational(row.get('disk_radius')) == geom['radius'], 'Primary source disk differs')
            require(type(row.get('coefficient_count')) is int and row['coefficient_count'] == 113, 'Primary source coefficient count differs')
            coefficients = row.get('coefficients')
            require(type(coefficients) is list and len(coefficients) == 113, 'Primary source coefficients missing')
            for coefficient in coefficients:
                require(type(coefficient) is dict and set(coefficient) == {'lo', 'hi'}, 'Exact real source coefficient interval required')
                require(rational(coefficient['lo']) <= rational(coefficient['hi']), 'Inverted source coefficient interval')
            for field in ('source_majorant', 'source_uniform_tail'):
                value = rational(row.get(field))
                require(value >= 0 and (fabricated or value > 0), 'Invalid primary analytic source bound')
            require(rational(row['source_majorant']) == M, 'Primary exact source majorant differs')
            require(rational(row['source_uniform_tail']) == tail, 'Primary exact source Cauchy tail differs')
            error = tail
        else:
            interval = row.get('interval')
            require(type(interval) is list and len(interval) == 2 and tuple(rational(v) for v in interval) == geom['interval'], 'Independent source-parent coverage differs')
            require(rational(row.get('analytic_radius')) == geom['radius'] and rational(row.get('disk_bound')) > 0, 'Independent source disk differs')
            require(rational(row['disk_bound']) == M, 'Independent exact source majorant differs')
            require(type(row.get('chosen_coefficients_sha256')) is str and re.fullmatch('[0-9a-f]{64}', row['chosen_coefficients_sha256']), 'Independent chosen coefficient identity absent')
            for field in ('analytic_tail', 'coefficient_error', 'uniform_source_error'):
                require(rational(row.get(field)) >= 0, 'Negative independent source error')
            if not fabricated:
                require(rational(row['analytic_tail']) > 0 and rational(row['uniform_source_error']) > 0, 'Actual independent source-model term omitted')
                require(rational(row['analytic_tail']) == tail, 'Independent exact source Cauchy tail differs')
            else:
                require(rational(row['analytic_tail']) == 0, 'Exact finite-polynomial fabricated source must have zero analytic tail')
            error = rational(row['uniform_source_error'])
            coefficient = rational(row['coefficient_error'])
            analytic = rational(row['analytic_tail'])
            ulp = Q(1, 2**512)
            for value in (error, coefficient):
                require((value * 2**512).denominator == 1, 'Independent source error must be outward 512-bit dyadic')
            require(error >= analytic and max(Q(0), analytic + coefficient - ulp) <= error <= analytic + coefficient + ulp, 'Independent source/coefficient outward-rounding relation differs')
        width = 2 * geom['half']
        model_radii[row['source']]['W'] += width * error
        model_radii[row['source']]['U'] += width * (Q(-9, 2) - geom['center']) * error
    require(seen == {(source, parent) for source in SOURCES for parent in range(11)}, 'Source-parent universe incomplete')
    return model_radii

def independent_mode_trace(budget):
    geometry = geometric_parents()
    expected = {(source, k, parent, child) for source in SOURCES for k in MOMENTA
                for parent, g in enumerate(geometry) for child in range(g['children'])}
    rows = budget.get('mode_error_rows')
    require(type(rows) is list and len(rows) == len(expected) == 666, 'Complete 666-row independent defect trace required')
    seen = set()
    by_key = {}
    fields = {'source_W', 'source_U', 'defect_W', 'defect_U', 'point_shift_W', 'point_shift_U', 'radius_round_W', 'radius_round_U'}
    for row in rows:
        require(type(row) is dict and row.get('source') in SOURCES and row.get('momentum') in MOMENTA, 'Invalid independent trace source/probe')
        require(type(row.get('parent')) is int and type(row.get('child')) is int, 'Exact trace parent/child required')
        key = (row['source'], row['momentum'], row['parent'], row['child'])
        require(key in expected and key not in seen, 'Duplicate or nonregistered independent trace child')
        seen.add(key)
        by_key[key] = row
        g = geometry[row['parent']]
        require(rational(row.get('width')) == 2 * g['half'] / g['children'], 'Independent trace child width differs')
        require(fields <= set(row), 'Independent trace error component omitted')
        for field in fields:
            require(rational(row[field]) >= 0, 'Negative/noncanonical independent trace error')
    require(seen == expected, 'Independent defect trace incomplete')
    ulp = Q(1, 2**512)
    parent_error = {(row['source'], row['parent']): rational(row['uniform_source_error']) for row in budget['source_model_rows']}
    whole = indexed_rows(budget['whole_rows'], 'independent trace aggregate')
    category_fields = {'source': ('source', 'source_model_disk_radius'), 'defect': ('defect', 'ode_defect_disk_radius'),
                       'point': ('point_shift', 'point_round_disk_radius'), 'radius': ('radius_round', 'radius_round_disk_radius')}
    for source in SOURCES:
        for k in MOMENTA:
            totals = {category: {'U': Q(0), 'W': Q(0)} for category in category_fields}
            slack = {category: {'U': Q(0), 'W': Q(0)} for category in category_fields}
            for parent, g in enumerate(geometry):
                for child in range(g['children']):
                    row = by_key[(source, k, parent, child)]
                    width = 2 * g['half'] / g['children']
                    E = parent_error[(source, parent)]
                    for category, (prefix, unused_field) in category_fields.items():
                        increment = {channel: rational(row[prefix + '_' + channel]) for channel in ('U', 'W')}
                        increment_slack = {'U': Q(0), 'W': Q(0)}
                        if category == 'source':
                            raw = {'U': E * width**2 / 2, 'W': E * width}
                            for channel in ('U', 'W'):
                                value = raw[channel] * 2**512
                                ceiling = Q((value.numerator + value.denominator - 1) // value.denominator, 2**512)
                                require(increment[channel] == ceiling, 'Independent source trace does not equal exact directed model increment')
                            # The source raw increments are known exactly; their
                            # displayed upward rounding adds no uncertainty here.
                            increment = raw
                        elif category in ('defect', 'point'):
                            for channel in ('U', 'W'):
                                value = increment[channel]
                                require((value * 2**512).denominator == 1, 'Independent defect/point trace must be outward dyadic512')
                                increment_slack[channel] = ulp if value > 0 else Q(0)
                        else:
                            require(all(0 <= value < ulp for value in increment.values()), 'Independent exact radius-round increment exceeds one ulp')
                        totals[category]['U'] += width * totals[category]['W'] + increment['U']
                        slack[category]['U'] += width * slack[category]['W'] + increment_slack['U']
                        totals[category]['W'] += increment['W']
                        slack[category]['W'] += increment_slack['W']
            for category, (unused_prefix, field) in category_fields.items():
                recorded = pair(whole[(source, k)].get(field), 'independent trace whole category')
                for channel in ('U', 'W'):
                    upper = totals[category][channel]
                    allowance = slack[category][channel]
                    require(recorded[channel] <= upper, 'Independent whole category exceeds trace telescope')
                    if allowance:
                        require(upper - recorded[channel] < allowance, 'Independent whole category differs beyond strict one-ulp trace slack')
                    else:
                        require(recorded[channel] == upper, 'Independent exact/zero trace category differs')

def indexed_rows(value, label):
    require(type(value) is list and len(value) == 18, label + ': exactly 18 whole rows required')
    result = {}
    for row in value:
        require(type(row) is dict and row.get('source') in SOURCES, label + ': unknown source')
        k = row.get('momentum')
        rational(k)
        require(k in MOMENTA, label + ': nonregistered momentum')
        key = (row['source'], k)
        require(key not in result, label + ': duplicate source/probe row')
        result[key] = row
    require(set(result) == KEYS, label + ': incomplete source/probe universe')
    return result

def pair(value, label):
    require(type(value) is dict and set(value) == {'U', 'W'}, label + ': U/W pair required')
    result = {channel: rational(value[channel]) for channel in ('U', 'W')}
    require(all(v >= 0 for v in result.values()), label + ': negative error budget')
    return result

def rectangle(value, label):
    require(type(value) is dict and set(value) == {'real', 'imag'}, label + ': real/imag rectangle required')
    result = {}
    for component in ('real', 'imag'):
        entry = value[component]
        require(type(entry) is dict and set(entry) == {'lo', 'hi'}, label + ': exact endpoints required')
        lo, hi = rational(entry['lo']), rational(entry['hi'])
        require(lo <= hi, label + ': inverted interval')
        result[component] = (lo, hi)
    radius = sum(((hi - lo) / 2 for lo, hi in result.values()), Q(0))
    return result, radius

def cap_bounds(source, k, route):
    delta = Q(1, 128)
    D = delta * (2 - delta)
    ell = 1 / (5 - delta)
    H = Q(3, 8)**63
    if route == 'primary':
        if source == SOURCES[0]:
            g = delta * H * (4 * ell**2 + 4 * ell / D**2 + 4 / D**4)
        else:
            g = delta * H * (4 * ell**2 + 2 * ell + (4 * ell + 4) / D**2 + 4 / D**4)
        return {'U': g / 2, 'W': g}
    p = 2 * (1 - delta) / D**2
    rho = 0 if source == SOURCES[0] else 1
    return {'U': H * ((p + rho + 2 * ell) / 2 + 1 + delta * (3 * ell**2 + 2 * k + 2 * ell)),
            'W': H * (p + rho + 2 * ell + 2 * k + delta * (4 * k**2 + 4 * k * ell + 6 * ell**2))}

def validate_output(data, budget, route, fabricated=False, enforce_gate=None):
    require(route in CONFIGURATIONS, 'Known route required')
    require(type(data) is dict and type(budget) is dict, 'Data and budget objects required')
    require(type(data.get('schema_version')) is int and data['schema_version'] == 1, 'Exact data schema_version 1 required')
    require(data.get('configuration') == CONFIGURATIONS[route], 'Route configuration differs')
    scope = 'FABRICATED_ONLY' if fabricated else OUTPUT_SCOPE
    require(data.get('scope') == scope, 'Output target scope differs')
    require(type(data.get('archive_arrays_decoded')) is int and data['archive_arrays_decoded'] == 0, 'Retained-array decode forbidden')
    require(data.get('original_binary80_state_error') == 'NOT_ENCLOSED' and data.get('full_twelve_case_certificate') == 'UNRESOLVED', 'Scientific status promotion forbidden')
    require(type(data.get('source_callbacks')) is int and data['source_callbacks'] == (0 if fabricated else 22), 'Physical source callback count differs')
    require(type(budget.get('physical_source_evaluations')) is int and budget['physical_source_evaluations'] == data['source_callbacks'], 'Mandatory budget callback counter differs')
    decode_field = 'archive_arrays_decoded' if route == 'primary' else 'retained_arrays_decoded'
    require(type(budget.get(decode_field)) is int and budget[decode_field] == 0, 'Mandatory budget decode counter absent or nonzero')
    model_radii = source_models(data, budget, route, fabricated)
    if route == 'independent':
        independent_mode_trace(budget)
    if 'scope' in budget:
        require(budget['scope'] == scope, 'Budget target scope differs')
    for field in ('original_binary80_state_error', 'full_twelve_case_certificate'):
        if field in budget:
            require(budget[field] == data[field], 'Budget scientific status differs')
    if 'retained_arrays_decoded' in budget:
        require(type(budget['retained_arrays_decoded']) is int and budget['retained_arrays_decoded'] == 0, 'Budget retained-array decoding forbidden')
    if 'physical_source_evaluations' in budget:
        require(type(budget['physical_source_evaluations']) is int and budget['physical_source_evaluations'] == data['source_callbacks'], 'Budget callback count differs')
    if enforce_gate is None:
        enforce_gate = not fabricated
    require(type(enforce_gate) is bool, 'Boolean gate policy required')
    rows = indexed_rows(data.get('whole_rows'), 'data')
    budgets = indexed_rows(budget.get('whole_rows'), 'budget')
    reviewed = {}
    max_radius = Q(0)
    for key in sorted(KEYS):
        row, b = rows[key], budgets[key]
        k = rational(key[1])
        factor = Q(1 if k == 0 else 2)
        cap = pair(b.get('cap_disk_radius'), 'cap')
        source = pair(b.get('source_model_disk_radius'), 'source model')
        require(source == model_radii[key[0]], 'Whole source-model budget does not equal parent proof transport')
        excess = pair(b.get('export_excess_L1_radius'), 'export excess')
        complete = pair(b.get('complete_output_L1_radius'), 'complete radius')
        recorded = pair(row.get('total_absolute_radii'), 'recorded exported radius')
        require(pair(row.get('cap_error'), 'data cap') == cap, 'Data/budget cap differs')
        if not fabricated:
            require(cap == cap_bounds(key[0], k, route), 'Cap differs from exact registered route formula')
            require(all(value > 0 for value in cap.values()) and all(value > 0 for value in source.values()), 'Actual cap/source-model term omitted')
        if 'inflation_factor' in b:
            if type(b['inflation_factor']) is dict:
                require(all(v == factor for v in pair(b['inflation_factor'], 'inflation factor').values()), 'Disk/L1 inflation factor differs')
            else:
                require(rational(b['inflation_factor']) == factor, 'Disk/L1 inflation factor differs')
        if 'output_real_dimensions' in b:
            require(type(b['output_real_dimensions']) is int and b['output_real_dimensions'] == int(factor), 'Exact-real projection dimension differs')
        if route == 'primary':
            phase = pair(b.get('phase_model_disk_radius'), 'phase model')
            baseline = pair(b.get('baseline_L1_radius'), 'arithmetic baseline')
            require(pair(b.get('coefficient_and_arithmetic_baseline_L1_radius'), 'arithmetic alias') == baseline, 'Arithmetic baseline alias differs')
            inflation = pair(b.get('inflation_disk_radius'), 'inflation disk')
            require(all(inflation[channel] == cap[channel] + source[channel] + phase[channel] for channel in ('U', 'W')), 'Primary omitted inflation budget')
            if k == 0:
                require(all(value == 0 for value in phase.values()), 'k0 phase model must use its exact-zero proof')
            elif not fabricated:
                require(all(value > 0 for value in phase.values()), 'Actual nonzero-k finite phase remainder omitted')
            expected = {channel: baseline[channel] + factor * inflation[channel] + excess[channel] for channel in ('U', 'W')}
        else:
            defect = pair(b.get('ode_defect_disk_radius'), 'ODE defect')
            point = pair(b.get('point_round_disk_radius'), 'point rounding')
            radius_round = pair(b.get('radius_round_disk_radius'), 'radius rounding')
            expected = {channel: factor * (cap[channel] + source[channel] + defect[channel] + point[channel] + radius_round[channel]) + excess[channel] for channel in ('U', 'W')}
        require(complete == expected, 'Route budget/export identity differs')
        reviewed[key] = {}
        for channel in ('U', 'W'):
            box, actual = rectangle(row.get(channel), channel)
            require(actual == recorded[channel] == complete[channel], 'Recorded radius differs from actual exported endpoints')
            if k == 0:
                require(box['imag'] == (Q(0), Q(0)), 'k0 imaginary projection is not exact zero')
            if enforce_gate:
                require(actual <= GATE, 'Actual exported L1 radius exceeds new 1e-20 gate')
            max_radius = max(max_radius, actual)
            reviewed[key][channel] = box
    return {'status': 'PASS_FABRICATED_OUTPUT_STRUCTURE' if fabricated else 'PASS_REGISTERED_BD_PREHISTORY_OUTPUT',
            'route': route, 'scope': scope, 'whole_rows': 18, 'complex_rectangles': 36,
            'max_actual_exported_L1_radius': f'{max_radius.numerator}/{max_radius.denominator}',
            'width_gate_applied': enforce_gate, 'width_gate': '1/100000000000000000000',
            'original_binary80_state_error': 'NOT_ENCLOSED', 'full_twelve_case_certificate': 'UNRESOLVED'}, reviewed

def compare_outputs(primary_data, primary_budget, independent_data, independent_budget, fabricated=False, enforce_gate=None):
    p, pb = validate_output(primary_data, primary_budget, 'primary', fabricated, enforce_gate)
    i, ib = validate_output(independent_data, independent_budget, 'independent', fabricated, enforce_gate)
    overlap = []
    for key in sorted(KEYS):
        for channel in ('U', 'W'):
            components = {}
            for component in ('real', 'imag'):
                a, b = pb[key][channel][component], ib[key][channel][component]
                lo, hi = max(a[0], b[0]), min(a[1], b[1])
                require(lo <= hi, 'Independent boxes do not intersect: ' + '/'.join((*key, channel, component)))
                components[component] = {'lo': f'{lo.numerator}/{lo.denominator}', 'hi': f'{hi.numerator}/{hi.denominator}'}
            overlap.append({'source': key[0], 'momentum': key[1], 'channel': channel, 'intersection': components})
    return {'status': 'PASS_FABRICATED_SHARED_FIXTURE_INTERSECTIONS' if fabricated else 'PASS_36_INDEPENDENT_BD_PREHISTORY_RECTANGLE_INTERSECTIONS',
            'primary_review': p, 'independent_review': i, 'overlap_rectangles': 36, 'component_intersections': 72,
            'comparison': 'Separate route validation and exact interval intersections; no averaging', 'rows': overlap,
            'original_binary80_state_error': 'NOT_ENCLOSED', 'full_twelve_case_certificate': 'UNRESOLVED'}

def read_limited(path):
    p = Path(path)
    require(p.is_file() and not p.is_symlink() and p.stat().st_size <= 20 * 1024 * 1024, 'Regular bounded JSON output required')
    return load(p.read_bytes())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--primary-data', required=True)
    parser.add_argument('--primary-budget', required=True)
    parser.add_argument('--independent-data', required=True)
    parser.add_argument('--independent-budget', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--fabricated-shared-fixture', action='store_true')
    args = parser.parse_args()
    result = compare_outputs(read_limited(args.primary_data), read_limited(args.primary_budget), read_limited(args.independent_data), read_limited(args.independent_budget), args.fabricated_shared_fixture, enforce_gate=True)
    with Path(args.output).open('x') as handle:
        handle.write(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(result['status'])

if __name__ == '__main__':
    main()
