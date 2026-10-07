"""Fabricated-only outer guard mutation controls; never calls either real route."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import copy
import hashlib
import json
import subprocess
import sys

sys.dont_write_bytecode = True
from registration_guard import load, validate_contract, require, zero_counts
from validate_outputs import CONFIGURATIONS, MOMENTA, SOURCES, compare_outputs, geometric_parents, validate_output

def qs(value):
    return f'{value.numerator}/{value.denominator}'

def fixture(route):
    data = {'schema_version': 1, 'configuration': CONFIGURATIONS[route], 'scope': 'FABRICATED_ONLY',
            'archive_arrays_decoded': 0, 'source_callbacks': 0,
            'original_binary80_state_error': 'NOT_ENCLOSED', 'full_twelve_case_certificate': 'UNRESOLVED', 'whole_rows': []}
    budget = {'physical_source_evaluations': 0, 'whole_rows': []}
    budget['archive_arrays_decoded' if route == 'primary' else 'retained_arrays_decoded'] = 0
    models = []
    source_radii = {source: {'U': Q(0), 'W': Q(0)} for source in SOURCES}
    for source in SOURCES:
        for parent, g in enumerate(geometric_parents()):
            model = {'source': source, 'parent': parent, 'center': qs(g['center']), 'half_width': qs(g['half']), 'children': g['children']}
            if route == 'primary':
                M = g['majorants'][source]
                error = M / 2**112
                model.update(disk_radius=qs(g['radius']), source_majorant=qs(M), source_uniform_tail=qs(error), coefficient_count=113, coefficients=[{'lo': '0/1', 'hi': '0/1'} for unused in range(113)])
            else:
                error = Q(0)
                model.update(interval=[qs(v) for v in g['interval']], analytic_radius=qs(g['radius']), disk_bound=qs(g['majorants'][source]), chosen_coefficients_sha256=hashlib.sha256(b'fabricated zero coefficient fixture').hexdigest(), analytic_tail='0/1', coefficient_error='0/1', uniform_source_error='0/1')
            source_radii[source]['W'] += 2 * g['half'] * error
            source_radii[source]['U'] += 2 * g['half'] * (Q(-9, 2) - g['center']) * error
            models.append(model)
    (data if route == 'primary' else budget)['source_model_rows'] = models
    if route == 'independent':
        budget['mode_error_rows'] = []
        for source in SOURCES:
            for k in MOMENTA:
                for parent, g in enumerate(geometric_parents()):
                    for child in range(g['children']):
                        budget['mode_error_rows'].append({'source': source, 'momentum': k, 'parent': parent, 'child': child,
                            'width': qs(2 * g['half'] / g['children']), **{field: '0/1' for field in ['source_W', 'source_U', 'defect_W', 'defect_U', 'point_shift_W', 'point_shift_U', 'radius_round_W', 'radius_round_U']}})
        # A known dyadic synthetic defect is placed only in the last child,
        # producing exactly the advertised whole U/W category without a later
        # W-to-U prefix contribution. No physical equation is evaluated.
        final_parent = len(geometric_parents()) - 1
        final_child = geometric_parents()[-1]['children'] - 1
        for row in budget['mode_error_rows']:
            if row['parent'] == final_parent and row['child'] == final_child:
                row['defect_U'] = row['defect_W'] = qs(Q(1, 2**150))
    zero = {'U': '0/1', 'W': '0/1'}
    for source in SOURCES:
        for k in MOMENTA:
            factor = 1 if Q(k) == 0 else 2
            half = Q(1, 2**150)
            radii = {channel: qs(factor * (half + source_radii[source][channel])) for channel in ('U', 'W')}
            data_row = {'source': source, 'momentum': k, 'total_absolute_radii': radii, 'cap_error': copy.deepcopy(zero)}
            for channel in ('U', 'W'):
                complete_half = Q(radii[channel]) / factor
                real = {'lo': qs(Q(k) - complete_half), 'hi': qs(Q(k) + complete_half)}
                imag = {'lo': '0/1', 'hi': '0/1'} if factor == 1 else {'lo': qs(-complete_half), 'hi': qs(complete_half)}
                data_row[channel] = {'real': real, 'imag': imag}
            data['whole_rows'].append(data_row)
            source_pair = {channel: qs(source_radii[source][channel]) for channel in ('U', 'W')}
            row = {'source': source, 'momentum': k, 'cap_disk_radius': copy.deepcopy(zero), 'source_model_disk_radius': source_pair, 'complete_output_L1_radius': radii, 'export_excess_L1_radius': copy.deepcopy(zero), 'inflation_factor': f'{factor}/1'}
            if route == 'primary':
                baseline = {channel: qs(factor * half) for channel in ('U', 'W')}
                row.update(phase_model_disk_radius=copy.deepcopy(zero), baseline_L1_radius=baseline, coefficient_and_arithmetic_baseline_L1_radius=copy.deepcopy(baseline), inflation_disk_radius=copy.deepcopy(source_pair))
            else:
                row.update(ode_defect_disk_radius={'U': qs(half), 'W': qs(half)}, point_round_disk_radius=copy.deepcopy(zero), radius_round_disk_radius=copy.deepcopy(zero), output_real_dimensions=factor)
            budget['whole_rows'].append(row)
    return data, budget

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--contract', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    contract = load(Path(args.contract).read_bytes())
    validate_contract(contract)
    p, b = fixture('primary')
    i, ib = fixture('independent')
    passing = []
    for route, data, budget in [('primary', p, b), ('independent', i, ib)]:
        validate_output(data, budget, route, True, True)
        passing.append('complete exact synthetic ' + route)
    overlap = compare_outputs(p, b, i, ib, True, True)
    require(overlap['overlap_rectangles'] == 36, 'All fabricated box comparisons required')
    passing.append('36 exact shared-fixture intersections')
    rejected = []

    def reject(label, callback):
        try:
            callback()
        except (ValueError, KeyError, TypeError, ZeroDivisionError):
            rejected.append(label)
        else:
            raise RuntimeError('Mutation unexpectedly accepted: ' + label)

    def mutate(label, callback):
        dp, bp = copy.deepcopy(p), copy.deepcopy(b)
        callback(dp, bp)
        reject(label, lambda: validate_output(dp, bp, 'primary', True, True))

    mutate('missing whole row', lambda d, e: d['whole_rows'].pop())
    mutate('duplicate whole row', lambda d, e: d['whole_rows'].__setitem__(1, copy.deepcopy(d['whole_rows'][0])))
    mutate('extra whole row', lambda d, e: d['whole_rows'].append(copy.deepcopy(d['whole_rows'][0])))
    mutate('wrong source', lambda d, e: d['whole_rows'][0].__setitem__('source', 'unknown'))
    mutate('wrong probe', lambda d, e: d['whole_rows'][0].__setitem__('momentum', '2/1'))
    mutate('noncanonical probe', lambda d, e: d['whole_rows'][0].__setitem__('momentum', '0/2'))
    mutate('float endpoint', lambda d, e: d['whole_rows'][0]['U']['real'].__setitem__('lo', 0.0))
    mutate('nonfinite endpoint', lambda d, e: d['whole_rows'][0]['U']['real'].__setitem__('lo', 'NaN'))
    mutate('noncanonical endpoint', lambda d, e: d['whole_rows'][0]['U']['real'].__setitem__('lo', '0/2'))
    mutate('inverted endpoint', lambda d, e: d['whole_rows'][0]['U']['real'].__setitem__('lo', '1/1'))
    mutate('recorded exported radius mismatch', lambda d, e: d['whole_rows'][0]['total_absolute_radii'].__setitem__('U', '0/1'))
    mutate('missing model budget', lambda d, e: e['whole_rows'][0].pop('source_model_disk_radius'))
    mutate('negative cap budget', lambda d, e: e['whole_rows'][0]['cap_disk_radius'].__setitem__('U', '-1/1'))
    mutate('wrong inflation factor', lambda d, e: e['whole_rows'][0].__setitem__('inflation_factor', '2/1'))
    mutate('arithmetic alias drift', lambda d, e: e['whole_rows'][0]['coefficient_and_arithmetic_baseline_L1_radius'].__setitem__('U', '0/1'))
    mutate('omitted inflation term', lambda d, e: e['whole_rows'][1]['source_model_disk_radius'].__setitem__('U', '1/1000000000000000000000000000000'))
    mutate('binary80 error promoted', lambda d, e: d.__setitem__('original_binary80_state_error', 'ENCLOSED'))
    mutate('full certificate promoted', lambda d, e: d.__setitem__('full_twelve_case_certificate', 'PASS'))
    mutate('retained array decode', lambda d, e: d.__setitem__('archive_arrays_decoded', 1))
    mutate('false fabricated callback', lambda d, e: d.__setitem__('source_callbacks', 22))
    mutate('missing budget callback counter', lambda d, e: e.pop('physical_source_evaluations'))
    mutate('missing budget decode counter', lambda d, e: e.pop('archive_arrays_decoded'))
    mutate('missing source-parent models', lambda d, e: d.pop('source_model_rows'))
    mutate('duplicate source-parent model', lambda d, e: d['source_model_rows'].__setitem__(1, copy.deepcopy(d['source_model_rows'][0])))
    mutate('altered source-parent center', lambda d, e: d['source_model_rows'][0].__setitem__('center', '-5/1'))
    mutate('source coefficient absent', lambda d, e: d['source_model_rows'][0]['coefficients'].pop())
    mutate('wrong source child partition', lambda d, e: d['source_model_rows'][0].__setitem__('children', 99))
    mutate('wrong exact source majorant', lambda d, e: d['source_model_rows'][0].__setitem__('source_majorant', '1/1'))
    mutate('wrong exact source tail', lambda d, e: d['source_model_rows'][0].__setitem__('source_uniform_tail', '0/1'))
    mutate('budget probe missing', lambda d, e: e['whole_rows'].pop())

    wide_data, wide_budget = copy.deepcopy(p), copy.deepcopy(b)
    wide_data['whole_rows'][0]['U']['real'] = {'lo': '-1/1', 'hi': '1/1'}
    wide_data['whole_rows'][0]['total_absolute_radii']['U'] = '1/1'
    for field in ['complete_output_L1_radius', 'baseline_L1_radius', 'coefficient_and_arithmetic_baseline_L1_radius']:
        wide_budget['whole_rows'][0][field]['U'] = '1/1'
    base_wide = Q(1) - Q(wide_budget['whole_rows'][0]['inflation_disk_radius']['U'])
    for field in ['baseline_L1_radius', 'coefficient_and_arithmetic_baseline_L1_radius']:
        wide_budget['whole_rows'][0][field]['U'] = qs(base_wide)
    reject('actual endpoint radius exceeds gate despite consistent records', lambda: validate_output(wide_data, wide_budget, 'primary', True, True))
    shifted = copy.deepcopy(i)
    for component in ['lo', 'hi']:
        value = Q(shifted['whole_rows'][1]['U']['real'][component]) + 1
        shifted['whole_rows'][1]['U']['real'][component] = qs(value)
    reject('nonintersecting independent rectangle', lambda: compare_outputs(p, b, shifted, ib, True, True))
    for label, action in [('missing independent trace', lambda e: e.pop('mode_error_rows')),
                          ('duplicate independent trace child', lambda e: e['mode_error_rows'].__setitem__(1, copy.deepcopy(e['mode_error_rows'][0]))),
                          ('missing independent trace component', lambda e: e['mode_error_rows'][0].pop('defect_U')),
                          ('wrong independent dyadic source relation', lambda e: e['source_model_rows'][0].__setitem__('uniform_source_error', '1/1'))]:
        changed = copy.deepcopy(ib)
        action(changed)
        reject(label, lambda changed=changed: validate_output(i, changed, 'independent', True, True))
    for label, action in [('zeroed independent positive defect trace', lambda e: [row.__setitem__(field, '0/1') for row in e['mode_error_rows'] for field in ['defect_U', 'defect_W']]),
                          ('non-dyadic independent defect trace', lambda e: e['mode_error_rows'][0].__setitem__('defect_W', '1/3')),
                          ('overwide exact radius-round trace', lambda e: e['mode_error_rows'][0].__setitem__('radius_round_W', '1/1')),
                          ('wrong exact directed source trace', lambda e: e['mode_error_rows'][0].__setitem__('source_W', '1/1'))]:
        changed = copy.deepcopy(ib)
        action(changed)
        reject(label, lambda changed=changed: validate_output(i, changed, 'independent', True, True))
    for key, wrong in [('target', 'ROUNDED_ARRAY_TARGET'), ('represented_epsilon', '1/10000'), ('represented_Pi', '3/1'), ('source_degree', 111), ('cap_delta', '1/64'), ('input_scope', 'ARRAYS_ALLOWED'), ('zero_anchor_eta', '-5/1')]:
        changed = copy.deepcopy(contract)
        changed[key] = wrong
        reject('contract changed ' + key, lambda changed=changed: validate_contract(changed))
    changed = copy.deepcopy(contract)
    changed['cap_formulae']['primary'] = 'arbitrary bound'
    reject('arbitrary cap formula prose', lambda: validate_contract(changed))
    for key in contract['new_target_evaluations_before_freeze']:
        changed = copy.deepcopy(contract)
        changed['new_target_evaluations_before_freeze'].pop(key)
        reject('missing preexecution counter ' + key, lambda changed=changed: validate_contract(changed))
    reject('duplicate JSON keys', lambda: load('{"x":1,"x":2}'))
    reject('nonfinite JSON number', lambda: load('{"x":NaN}'))
    result = {'status': 'PASS_FABRICATED_OUTER_GUARD_CONTROLS', 'optimization': sys.flags.optimize,
              'passing_checks': passing, 'passing_checks_count': len(passing),
              'rejected_mutations': rejected, 'rejected_mutations_count': len(rejected),
              'registered_source_callbacks': 0, 'retained_array_decodes': 0, 'physical_trajectories': 0,
              'stored_state_comparisons': 0, 'likelihood_evaluations': 0,
              'scope': 'FABRICATED_ONLY; structural guards do not establish physical output readiness'}
    with Path(args.output).open('x') as handle:
        handle.write(json.dumps(result, sort_keys=True, indent=2) + '\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ['passing_checks', 'rejected_mutations']}))

if __name__ == '__main__':
    main()
