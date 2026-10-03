#!/usr/bin/env python3
"""Exact-decimal cross-route validation; negative science is a valid outcome."""
from __future__ import annotations
import argparse
from decimal import Decimal
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

TOL = Q('1e-12')
SCIENCE = Q('2e-7')
ATTRIBUTION = Q('2e-6')
SOURCES = ('positive_B', 'signed_uB')
SETTINGS = ('coarse', 'fine')
KS = (64, 128, 256)
COMMON = ('times R P F stored_R stored_P stored_F I_ab E_flow triangle_bound S_ab '
          'S_direct_reset DeltaR S_global_minus_reset D_S D_cont E_Q decomposition_error '
          'flow_operator_closure R_profile_max_error P_profile_max_error F_profile_max_error').split()
FINE = ('S_ab_double_global', 'S_native_minus_double')
FLOW_FIELDS = {'k', 'delta_u_real', 'delta_u_imag', 'delta_w_real', 'delta_w_imag',
               'weighted_signed_flow', 'weighted_triangle', 'weighted_direct_ledger'}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def number(value):
    require(isinstance(value, str), 'Scientific values must be decimal strings')
    decimal = Decimal(value)
    require(decimal.is_finite(), 'Nonfinite scientific decimal')
    return Q(decimal)

def gap(a, b):
    if isinstance(a, list):
        require(isinstance(b, list) and len(a) == len(b), 'Profile shapes differ')
        return max((gap(x, y) for x, y in zip(a, b)), default=Q(0))
    return abs(number(a) - number(b))

def close(a, b, name):
    require(abs(a-b) <= TOL, 'Identity or arithmetic gate failed: ' + name)

def scientific_outcome(rows):
    failures, witnesses = [], []
    for source, setting, row in rows:
        identity = {'source': source, 'setting': setting, 'K': row['K']}
        delta = number(row['stored_R'][-1])-number(row['stored_R'][0])
        actual = {key+'_profile_max_error':max(abs(number(a)-number(b)) for a,b in
                   zip(row[key],row['stored_'+key])) for key in ('R','P')}
        actual.update(D_cont=delta-number(row['I_ab']), E_flow=number(row['E_flow']))
        for key, value in actual.items():
            if abs(value) > SCIENCE:
                failures.append({**identity, 'gate': key, 'exact_rational_value':str(value)})
        ds = abs(delta-number(row['S_ab']))
        if ds > ATTRIBUTION and abs(actual['D_cont']) <= ds/10 \
                and abs(number(row['I_ab'])-number(row['S_ab'])) >= 9*ds/10:
            witnesses.append(identity)
    if failures:
        return 'CONSISTENCY_FAILURE', failures, []
    return ('LEDGER_ERROR_DEMONSTRATED' if witnesses else 'NO_GATE_SCALE_ATTRIBUTION'), [], witnesses

def check_row(row, setting, route):
    required = set(COMMON) | ({*FINE} if setting == 'fine' else set()) | {'K'}
    require(required <= row.keys(), 'Required common scientific fields missing')
    n, divisor = (129, 128) if setting == 'coarse' else (257, 256)
    for key in ('times', 'R', 'P', 'F', 'stored_R', 'stored_P', 'stored_F'):
        require(isinstance(row[key], list) and len(row[key]) == n, 'Incomplete profile: ' + key)
        for value in row[key]: number(value)
    require([number(x) for x in row['times']] == [Q(-5, 2)+Q(j, divisor) for j in range(n)],
            'Fixed inclusive interval/grid changed')
    v = {key: number(value) for key, value in row.items() if key != 'K' and not isinstance(value, list)}
    for component in ('R', 'P', 'F'):
        maximum = max(abs(number(a)-number(b)) for a,b in zip(row[component], row['stored_'+component]))
        close(v[component+'_profile_max_error'], maximum, component+' profile maximum')
    close(v['DeltaR'], number(row['stored_R'][-1])-number(row['stored_R'][0]), 'stored endpoint DeltaR')
    for field, expression in {
        'D_S': v['DeltaR']-v['S_ab'], 'D_cont': v['DeltaR']-v['I_ab'],
        'E_Q': v['I_ab']-v['S_ab'],
        'decomposition_error': v['D_S']-v['D_cont']-v['E_Q'],
        'flow_operator_closure': v['D_cont']-v['E_flow'],
        'S_global_minus_reset': v['S_ab']-v['S_direct_reset'],
    }.items(): close(v[field], expression, field)
    require(abs(v['decomposition_error']) <= TOL, 'Signed ledger decomposition failed')
    require(abs(v['D_S']-v['D_cont']-v['E_Q']) <= TOL, 'Actual signed ledger decomposition failed')
    require(v['triangle_bound'] >= 0 and abs(v['E_flow']) <= v['triangle_bound']+TOL,
            'Measured flow projection exceeds its triangle bound')
    if setting == 'fine':
        close(v['S_native_minus_double'], v['S_ab']-v['S_ab_double_global'], 'same-trajectory Simpson difference')
    if route == 'primary':
        close(v['I_ab'], v['independent_phase_primitive'], 'primary canonical direct phase')
        close(v['recurrence_primitive_gap'], v['recurrence_phase_primitive']-v['I_ab'], 'primary recurrence difference')
        require(abs(v['recurrence_primitive_gap']) <= TOL, 'Primary recurrence/direct discrepancy')
        require(abs(v['recurrence_phase_primitive']-v['I_ab']) <= TOL, 'Actual primary recurrence/direct discrepancy')
        require(v['quantization_only_RP_bound'] >= 0, 'Negative conditional quantization bound')
    else:
        close(v['I_recurrence_minus_direct'], v['I_recurrence']-v['I_ab'], 'independent recurrence difference')
        require(abs(v['I_recurrence_minus_direct']) <= TOL, 'Independent recurrence/direct discrepancy')
        require(abs(v['I_recurrence']-v['I_ab']) <= TOL, 'Actual independent recurrence/direct discrepancy')

def check_result(result, route):
    require(result['schema_version'] == 1 and result['fixed_cases'] == 12, 'Wrong case count/schema')
    require(result['status'] == 'COMPLETED_SAVED_DATA_DIAGNOSTIC' and result['fatal_failures'] == [],
            'Incomplete or internally failed route')
    require(result['old_metric_status'] == 'FAIL', 'Earlier metric status changed')
    expected_route = 'primary_integer_three_moment' if route == 'primary' else 'independent_mpf_bare_bilinear_fourier'
    require(result['route'] == expected_route, 'Wrong independent route identity')
    resources = result['resources']
    require(0 <= resources['seconds'] <= 900 and 0 < resources['peak_rss_kib'] <= 262144,
            'Reported full-route budget exceeded')
    require(resources['limit_seconds'] == 900 and resources['limit_peak_rss_kib'] == 262144,
            'Resource limits changed')
    records = result['records']
    require([(r['source'],r['setting']) for r in records] == [(s,t) for s in SOURCES for t in SETTINGS],
            'Exactly four ordered source/resolution records required')
    maximum, selected, lookup = Q(0), [], {}
    for record in records:
        source, setting = record['source'], record['setting']
        require(set(record['levels']) == {'80','100'}, 'Precision levels changed')
        for level in ('80','100'):
            rows = record['levels'][level]
            require(tuple(row['K'] for row in rows) == KS and all(type(row['K']) is int for row in rows),
                    'Fixed ordered cutoffs changed')
            for row in rows:
                check_row(row, setting, route)
                lookup[(source, setting, level, row['K'])] = row
                if level == '100': selected.append((source, setting, row))
        record_gap = Q(0)
        for lo, hi in zip(record['levels']['80'], record['levels']['100']):
            require(lo.keys() == hi.keys(), 'Precision field universe changed')
            for key in lo:
                if key != 'K': record_gap = max(record_gap, gap(lo[key], hi[key]))
        require(record_gap <= TOL, 'Full reported precision gap exceeds threshold')
        require(record_gap <= Q(record['precision_gap_exact_rational']) <= TOL,
                'Reported precision gap understates field universe')
        maximum = max(maximum, record_gap)
        if route == 'primary':
            require(set(record['serialization_max_error']) == {'80','100'}, 'Missing serialization evidence')
            require(all(0 <= number(x) <= TOL for x in record['serialization_max_error'].values()),
                    'Serialization gate failed')
    require(maximum <= Q(result['precision_gap_max_exact_rational']) <= TOL, 'Overall precision gap differs')
    classification, failures, witnesses = scientific_outcome(selected)
    require(result['classification'] == classification and result['attribution_cases'] == witnesses,
            'Classification not supported by fixed numerical criteria')
    reported = {(x['source'],x['setting'],x['K'],x['gate']) for x in result['failures']}
    actual = {(x['source'],x['setting'],x['K'],x['gate']) for x in failures}
    require(reported == actual and len(result['failures']) == len(actual), 'Scientific failure list differs')
    for failure in result['failures']:
        row = lookup[(failure['source'],failure['setting'],'100',failure['K'])]
        require(number(failure['value']) == number(row[failure['gate']])
                and number(failure['threshold']) == SCIENCE, 'Scientific failure evidence differs from row')
    return lookup, {'classification':classification, 'failures':failures, 'attribution_cases':witnesses,
                    'profile_precision_gap_exact_rational':str(maximum)}

def check_execution(receipts):
    physical = [x for x in receipts if x.get('physical_route') is True]
    require([x['name'] for x in physical] == ['primary','independent'], 'Two complete physical executions required')
    for r in physical:
        require(r['status'] == 'PASS_EXECUTION' and type(r['exit_code']) is int and r['exit_code'] == 0
                and r['timed_out'] is False and r['memory_exceeded'] is False and 'exception' not in r,
                'Execution exit/resource failure')
        require(r['timeout_seconds'] == 900 and r['memory_limit_kib'] == 262144,
                'Outer execution limits changed')
        require(r['process_creation_limit'] == 0, 'Single-process execution containment missing')
        require(0 <= r['elapsed_seconds'] <= 900 and 0 < r['peak_rss_kib'] <= 262144,
                'Actual complete-route resource budget exceeded')
    return [{key:r[key] for key in ('name','elapsed_seconds','peak_rss_kib','exit_code')} for r in physical]

def cross_validate(primary, independent, receipts):
    p, ps = check_result(primary, 'primary')
    i, ins = check_result(independent, 'independent')
    require(ps['classification'] == ins['classification'] and ps['attribution_cases'] == ins['attribution_cases'],
            'Routes disagree on scientific classification')
    maximum = Q(0)
    for identity, row in p.items():
        other = i[identity]
        for key in COMMON + (list(FINE) if identity[1] == 'fine' else []):
            maximum = max(maximum, gap(row[key],other[key]))
    require(maximum <= TOL, 'Cross-route common scientific fields disagree')
    return {'schema_version':1, 'status':'VERIFIED_COMPLETE_DIAGNOSTIC',
            'classification':ps['classification'], 'old_metric_status':'FAIL', 'fixed_cases':12,
            'cross_route_max_gap_exact_rational':str(maximum),
            'routes':{'primary':ps,'independent':ins}, 'executions':check_execution(receipts),
            'interpretation':'Saved-data post-support diagnostic only; no metric calibration or higher-dimensional evidence.'}

def check_flow_files(result, directory, read_json, file_sha):
    expected = {f'flow_{s}_{t}_{d}.jsonl' for s in SOURCES for t in SETTINGS for d in (80,100)}
    entries = result['flow_files']
    require(len(entries) == 8 and {x['path'] for x in entries} == expected, 'Streamed flow membership differs')
    for item in entries:
        path = directory/item['path']
        require(path.is_file() and not path.is_symlink() and path.stat().st_size == item['bytes']
                and file_sha(path) == item['sha256'], 'Streamed flow bytes changed')
    maximum = Q(0)
    for record in result['records']:
        source, setting = record['source'],record['setting']
        n = 8192 if setting == 'coarse' else 16384
        sums = {d:{k:Q(0) for k in ('weighted_signed_flow','weighted_triangle','weighted_direct_ledger')}
                for d in ('80','100')}
        count = 0
        with (directory/f'flow_{source}_{setting}_80.jsonl').open() as a, \
             (directory/f'flow_{source}_{setting}_100.jsonl').open() as b:
            from itertools import zip_longest
            for count, pair in enumerate(zip_longest(a,b),1):
                require(None not in pair, 'Stream precision lengths differ')
                records = [json.loads(line) for line in pair]
                for d, row in zip(('80','100'),records):
                    require(set(row) == FLOW_FIELDS|{'node_index'} and type(row['node_index']) is int
                            and row['node_index'] == count-1, 'Incomplete/reordered endpoint defects')
                    values = {k:number(row[k]) for k in FLOW_FIELDS}
                    require(values['weighted_triangle'] >= abs(values['weighted_signed_flow'])-TOL,
                            'Per-node triangle bound failed')
                    for key in sums[d]: sums[d][key] += values[key]
                    if count in (n//4,n//2,n):
                        row_index = (n//4,n//2,n).index(count)
                        target = record['levels'][d][row_index]
                        for key,outkey in (('weighted_signed_flow','E_flow'),('weighted_triangle','triangle_bound'),
                                           ('weighted_direct_ledger','I_ab')):
                            close(sums[d][key],number(target[outkey]),'streamed '+outkey+' prefix')
                maximum = max(maximum, *(gap(records[0][key],records[1][key]) for key in FLOW_FIELDS))
        require(count == n, 'Not all endpoint node defects retained')
    require(maximum <= TOL and maximum <= Q(result['precision_gap_max_exact_rational']),
            'Streamed endpoint precision gap failed')
    return str(maximum)

def main():
    parser = argparse.ArgumentParser()
    for name in ('checkpoint-root','primary','independent','execution-receipts','output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--freeze-commit',required=True)
    args = parser.parse_args()
    root = args.checkpoint_root.absolute()
    # Authenticate the helper before importing it. The frozen driver authenticates this entry.
    raw = (root/'FULL_REGISTRATION.json').read_bytes()
    require(hashlib.sha256(raw).hexdigest() == args.registration_sha256, 'Registration hash differs')
    registration = json.loads(raw)
    helper_path = root/'code/ledger_integrity.py'
    require(not helper_path.is_symlink() and hashlib.sha256(helper_path.read_bytes()).hexdigest()
            == registration['files']['code/ledger_integrity.py'], 'Unregistered integrity helper')
    spec = importlib.util.spec_from_file_location('ledger_integrity',helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    evidence = helper.verify_frozen(root,args.registration_sha256,args.freeze_commit)
    require(helper.file_sha(Path(__file__)) == evidence['frozen_files']['code/validate_ledger.py'],
            'Validator source differs')
    output = helper.fresh_external(root,args.output_dir)
    output.mkdir(parents=True)
    try:
        primary,independent = helper.read_json(args.primary),helper.read_json(args.independent)
        for result in (primary,independent):
            require(result['input_manifest_sha256'] == evidence['frozen_files']['inputs/INPUT_MANIFEST.json'],
                    'Wrong immutable input lineage')
            require(result['input_lineage_public_commit'] == '71d00cc423e9049c8166b7ee7afbefbac28d8a18',
                    'Wrong source science commit')
        require(primary['registration_sha256'] == args.registration_sha256 and primary['freeze_commit'] == args.freeze_commit,
                'Primary provenance differs')
        require(independent['provenance']['registration_sha256'] == args.registration_sha256
                and independent['provenance']['public_freeze_commit'] == args.freeze_commit, 'Independent provenance differs')
        report = cross_validate(primary,independent,helper.read_json(args.execution_receipts))
        report['streamed_flow_precision_gap_exact_rational'] = check_flow_files(independent,args.independent.parent,
                                                                              helper.read_json,helper.file_sha)
        report.update(registration_sha256=args.registration_sha256,public_freeze_commit=args.freeze_commit,
                      result_sha256={'primary':helper.file_sha(args.primary),'independent':helper.file_sha(args.independent)},
                      python_optimization=sys.flags.optimize)
        require(helper.verify_frozen(root,args.registration_sha256,args.freeze_commit) == evidence,
                'Frozen files changed during validation')
        helper.write_json(output/'VALIDATION.json',report)
        print(json.dumps({'status':report['status'],'classification':report['classification']}))
        return 0
    except Exception as error:
        helper.write_json(output/'FAILED_VALIDATION.json',{'error_type':type(error).__name__,'error':str(error),
                                                         'old_metric_status':'FAIL'})
        raise

if __name__ == '__main__':
    sys.exit(main())
