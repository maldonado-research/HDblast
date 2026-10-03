#!/usr/bin/env python3
"""Synthetic schema/mutation tests. No physical files, arrays, or computations."""
from __future__ import annotations
import argparse
import copy
from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

DEFAULT_SOURCE = Path('/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/code/validate_ledger.py')

def dec(x):
    return format(Decimal(str(x)), 'f')

def synthetic_fixture():
    """Fabricate all full profile/schema fields from constants and a rational grid."""
    results = []
    for route, identity in (('primary', 'primary_integer_three_moment'),
                            ('independent', 'independent_mpf_bare_bilinear_fourier')):
        records = []
        for source in ('positive_B', 'signed_uB'):
            for setting in ('coarse', 'fine'):
                divisor = 128 if setting == 'coarse' else 256
                levels = {}
                for level in ('80', '100'):
                    rows = []
                    for k in (64, 128, 256):
                        row = {field: '0' for field in (
                            'I_ab E_flow triangle_bound S_ab S_direct_reset DeltaR S_global_minus_reset '
                            'D_S D_cont E_Q decomposition_error flow_operator_closure '
                            'R_profile_max_error P_profile_max_error F_profile_max_error').split()}
                        row['K'] = k
                        row['times'] = [format(Decimal(-5)/2+Decimal(j)/divisor, 'f')
                                        for j in range(divisor+1)]
                        for field in ('R','P','F','stored_R','stored_P','stored_F'):
                            row[field] = ['0']*(divisor+1)
                        if setting == 'fine':
                            row.update(S_ab_double_global='0', S_native_minus_double='0')
                        if route == 'primary':
                            row.update(independent_phase_primitive='0', recurrence_phase_primitive='0',
                                       recurrence_primitive_gap='0', quantization_only_RP_bound='0')
                        else:
                            row.update(I_recurrence='0', I_recurrence_minus_direct='0')
                        rows.append(row)
                    levels[level] = rows
                record = {'source': source, 'setting': setting, 'levels': levels,
                          'precision_gap_exact_rational':'0'}
                if route == 'primary':
                    record['serialization_max_error'] = {'80':'0','100':'0'}
                records.append(record)
        results.append({'schema_version':1, 'fixed_cases':12,
                        'status':'COMPLETED_SAVED_DATA_DIAGNOSTIC', 'fatal_failures':[],
                        'old_metric_status':'FAIL','route':identity,
                        'resources':{'seconds':1.0,'peak_rss_kib':1024,
                                     'limit_seconds':900,'limit_peak_rss_kib':262144},
                        'records':records,'precision_gap_max_exact_rational':'0',
                        'classification':'NO_GATE_SCALE_ATTRIBUTION',
                        'failures':[],'attribution_cases':[]})
    receipts = [{'name':r, 'physical_route':True, 'status':'PASS_EXECUTION',
                 'exit_code':0,'timed_out':False,'memory_exceeded':False,
                 'timeout_seconds':900,'memory_limit_kib':262144,
                 'process_creation_limit':0,
                 'elapsed_seconds':2.0,'peak_rss_kib':2048}
                for r in ('primary','independent')]
    return results+[receipts]

def rows(result):
    for record in result['records']:
        for level in ('80','100'):
            yield from record['levels'][level]

def first(result, level='80', record=0):
    return result['records'][record]['levels'][level][0]

def outcome(result):
    """Test oracle expressed independently using finite Decimal comparisons."""
    failures, witnesses = [], []
    for record in result['records']:
        for row in record['levels']['100']:
            ident = {'source':record['source'],'setting':record['setting'],'K':row['K']}
            for gate in ('R_profile_max_error','P_profile_max_error','D_cont','E_flow'):
                if abs(Decimal(row[gate])) > Decimal('0.0000002'):
                    failures.append(dict(ident,gate=gate,value=row[gate],threshold='0.0000002'))
            ds = abs(Decimal(row['D_S']))
            if (ds > Decimal('0.000002') and abs(Decimal(row['D_cont'])) <= ds/10
                    and abs(Decimal(row['E_Q'])) >= ds*Decimal('0.9')):
                witnesses.append(ident)
    result['classification'] = ('CONSISTENCY_FAILURE' if failures else
                                'LEDGER_ERROR_DEMONSTRATED' if witnesses else
                                'NO_GATE_SCALE_ATTRIBUTION')
    result['failures'] = failures
    result['attribution_cases'] = [] if failures else witnesses

def apply_both(data, f):
    for result in data[:2]:
        f(result)

def ledger_row(row, amount):
    """Self-consistent synthetic reset/global ledger with direct integral zero."""
    row.update(S_ab=dec(-Decimal(amount)), S_direct_reset=dec(-Decimal(amount)),
               D_S=dec(amount), E_Q=dec(amount))
    if 'S_ab_double_global' in row:
        row['S_ab_double_global'] = dec(-Decimal(amount))

def ledger_fixture(data, amount):
    for result in data[:2]:
        for row in rows(result):
            ledger_row(row, amount)
        outcome(result)

def failure_fixture(data, amount='0.0000003'):
    for result in data[:2]:
        for row in rows(result):
            row['R'][1] = amount
            row['R_profile_max_error'] = amount
        outcome(result)

def precision_fixture(data, amount, reported=None):
    result=data[0]
    first(result,'100')['quantization_only_RP_bound']=amount
    value=str(Fraction(Decimal(amount))) if reported is None else reported
    result['records'][0]['precision_gap_exact_rational']=value
    result['precision_gap_max_exact_rational']=value

def cross_profile_fixture(data, amount):
    for row in rows(data[1]):
        row['R'][1]=amount
        row['stored_R'][1]=amount

def shared_biased_recurrence(data):
    for ri,result in enumerate(data[:2]):
        for row in rows(result):
            row.update(I_ab='0.00000001', D_cont='-0.00000001', E_Q='0.00000001',
                       E_flow='-0.00000001', triangle_bound='0.00000001')
            if ri == 0:
                row['recurrence_phase_primitive']='0.00000001'
            else:
                row['I_recurrence']='0.00000001'

def make_tests():
    tests=[]
    def add(name, accepted, mutate=lambda d:None):
        tests.append((name,accepted,mutate))
    add('complete_no_attribution_fixture',True)
    add('complete_ledger_attribution_fixture',True,lambda d:ledger_fixture(d,'0.000003'))
    add('complete_consistency_failure_fixture',True,failure_fixture)
    add('arithmetic_tolerance_exact_boundary',True,lambda d:apply_both(d,lambda r:[row.__setitem__('decomposition_error','0.000000000001') for row in rows(r)]))
    add('arithmetic_tolerance_overrun',False,lambda d:first(d[0]).__setitem__('decomposition_error','0.000000000001000001'))
    add('science_gate_exact_boundary',True,lambda d:failure_fixture(d,'0.0000002'))
    add('science_gate_overrun_honest_classification',True,lambda d:failure_fixture(d,'0.000000200000000001'))
    def false_science(d):
        failure_fixture(d,'0.000000200000000001')
        for r in d[:2]:r.update(classification='NO_GATE_SCALE_ATTRIBUTION',failures=[])
    add('science_gate_overrun_false_classification',False,false_science)
    add('attribution_exact_boundary',True,lambda d:ledger_fixture(d,'0.000002'))
    add('attribution_over_boundary_honest',True,lambda d:ledger_fixture(d,'0.000002000000000001'))
    def false_attribution_boundary(d):
        ledger_fixture(d,'0.000002')
        for r in d[:2]:
            r['classification']='LEDGER_ERROR_DEMONSTRATED'
            r['attribution_cases']=[{'source':'positive_B','setting':'coarse','K':64}]
    add('attribution_boundary_false_witness',False,false_attribution_boundary)
    add('full_precision_exact_boundary',True,lambda d:precision_fixture(d,'0.000000000001'))
    add('full_precision_overrun',False,lambda d:precision_fixture(d,'0.000000000001000001'))
    add('record_precision_underreported',False,lambda d:precision_fixture(d,'0.0000000000001','0'))
    def overall_gap_lie(d):
        precision_fixture(d,'0.0000000000001')
        d[0]['precision_gap_max_exact_rational']='0'
    add('overall_precision_underreported',False,overall_gap_lie)
    add('cross_profile_exact_boundary',True,lambda d:cross_profile_fixture(d,'0.000000000001'))
    add('cross_profile_overrun_same_reported_maxima',False,lambda d:cross_profile_fixture(d,'0.000000000001000001'))
    add('shared_biased_recurrence_direct_primitive_guard',False,shared_biased_recurrence)
    add('independent_recurrence_direct_discrepancy',False,lambda d:first(d[1]).update(I_recurrence='0.00000001',I_recurrence_minus_direct='0.00000001'))
    add('independent_recurrence_gap_falsely_zero',False,lambda d:first(d[1]).__setitem__('I_recurrence','0.00000001'))
    add('primary_recurrence_gap_falsely_zero',False,lambda d:first(d[0]).__setitem__('recurrence_phase_primitive','0.00000001'))
    add('record_missing',False,lambda d:d[0]['records'].pop())
    add('record_reordered',False,lambda d:d[0]['records'].reverse())
    add('record_duplicate',False,lambda d:d[0]['records'].__setitem__(1,copy.deepcopy(d[0]['records'][0])))
    add('cutoff_missing',False,lambda d:d[0]['records'][0]['levels']['80'].pop())
    add('cutoff_reordered',False,lambda d:d[0]['records'][0]['levels']['80'].reverse())
    add('cutoff_duplicate',False,lambda d:first(d[0]).__setitem__('K',128))
    add('cutoff_boolean',False,lambda d:first(d[0]).__setitem__('K',True))
    add('precision_level_missing',False,lambda d:d[0]['records'][0]['levels'].pop('80'))
    add('precision_level_extra',False,lambda d:d[0]['records'][0]['levels'].__setitem__('120',[]))
    add('fine_field_missing',False,lambda d:first(d[0],record=1).pop('S_native_minus_double'))
    add('full_profile_member_missing',False,lambda d:first(d[0])['R'].pop())
    add('inclusive_grid_shifted',False,lambda d:first(d[0])['times'].__setitem__(1,'-2.49218749'))
    add('full_profile_maximum_underreported',False,lambda d:first(d[0])['R'].__setitem__(1,'0.00000001'))
    def profile_gate_laundering(d):
        for result in d[:2]:
            for row in rows(result):
                row['R'][1]='0.0000002000005'
                row['R_profile_max_error']='0.0000002'
    add('true_profile_science_overrun_laundered_by_reported_maximum',False,profile_gate_laundering)
    def direct_defect_gate_laundering(d):
        for ri,result in enumerate(d[:2]):
            for row in rows(result):
                row.update(I_ab='-0.0000002000005', D_cont='0.0000002', E_flow='0.0000002',
                           triangle_bound='0.0000002', E_Q='-0.0000002000005',
                           decomposition_error='0.0000000000005')
                if ri == 0:
                    row.update(independent_phase_primitive='-0.0000002000005',
                               recurrence_phase_primitive='-0.0000002000005')
                else:row['I_recurrence']='-0.0000002000005'
    add('true_direct_defect_science_overrun_laundered_by_reported_value',False,direct_defect_gate_laundering)
    def recurrence_laundering(d,route):
        result=d[route]
        for row in rows(result):
            if route == 0:
                row.update(recurrence_phase_primitive='0.000000000002',recurrence_primitive_gap='0.000000000001')
            else:
                row.update(I_recurrence='0.000000000002',I_recurrence_minus_direct='0.000000000001')
    add('true_primary_recurrence_overrun_laundered_by_reported_gap',False,lambda d:recurrence_laundering(d,0))
    add('true_independent_recurrence_overrun_laundered_by_reported_gap',False,lambda d:recurrence_laundering(d,1))
    def decomposition_laundering(d):
        for result in d[:2]:
            for row in rows(result):
                row.update(S_ab='-0.000000000001',S_direct_reset='-0.000000000001',
                           D_S='0.000000000002',decomposition_error='0.000000000001')
                if 'S_ab_double_global' in row:row['S_ab_double_global']='-0.000000000001'
    add('true_signed_decomposition_overrun_laundered_by_reported_residual',False,decomposition_laundering)
    add('stored_endpoint_delta_underreported',False,lambda d:first(d[0])['stored_R'].__setitem__(-1,'0.00000001'))
    add('flow_closure_lie',False,lambda d:first(d[0]).__setitem__('flow_operator_closure','0.00000001'))
    add('negative_triangle_bound',False,lambda d:first(d[0]).__setitem__('triangle_bound','-0.00000001'))
    add('negative_quantization_bound',False,lambda d:first(d[0]).__setitem__('quantization_only_RP_bound','-0.00000001'))
    add('serialization_gap_overrun',False,lambda d:d[0]['records'][0]['serialization_max_error'].__setitem__('80','0.000000000001000001'))
    add('serialization_level_missing',False,lambda d:d[0]['records'][0]['serialization_max_error'].pop('80'))
    add('false_classification',False,lambda d:d[0].__setitem__('classification','LEDGER_ERROR_DEMONSTRATED'))
    add('false_witness',False,lambda d:d[0]['attribution_cases'].append({'source':'positive_B','setting':'coarse','K':64}))
    def duplicate_witness(d):
        ledger_fixture(d,'0.000003')
        d[0]['attribution_cases'].append(copy.deepcopy(d[0]['attribution_cases'][0]))
    add('duplicate_witness',False,duplicate_witness)
    def failure_list_missing(d):
        failure_fixture(d)
        d[0]['failures'].pop()
    add('scientific_failure_missing',False,failure_list_missing)
    def failure_list_duplicate(d):
        failure_fixture(d)
        d[0]['failures'].append(copy.deepcopy(d[0]['failures'][0]))
    add('scientific_failure_duplicate',False,failure_list_duplicate)
    def failure_value_lie(d):
        failure_fixture(d)
        d[0]['failures'][0]['value']='123'
    add('scientific_failure_value_lie',False,failure_value_lie)
    def failure_threshold_lie(d):
        failure_fixture(d)
        d[0]['failures'][0]['threshold']='0.0000003'
    add('scientific_failure_threshold_lie',False,failure_threshold_lie)
    def failure_nonfinite_value(d):
        failure_fixture(d)
        d[0]['failures'][0]['value']='NaN'
    add('scientific_failure_nonfinite_value',False,failure_nonfinite_value)
    add('route_identity_primary_wrong',False,lambda d:d[0].__setitem__('route','independent_mpf_bare_bilinear_fourier'))
    add('route_identity_independent_wrong',False,lambda d:d[1].__setitem__('route','primary_integer_three_moment'))
    add('fixed_case_count_wrong',False,lambda d:d[0].__setitem__('fixed_cases',11))
    add('old_metric_status_upgraded',False,lambda d:d[0].__setitem__('old_metric_status','PASS'))
    add('route_status_incomplete',False,lambda d:d[0].__setitem__('status','RUNNING'))
    add('route_fatal_failure_retained',False,lambda d:d[0]['fatal_failures'].append({'error':'synthetic'}))
    def resource_boundary(d):
        for r in d[:2]:r['resources'].update(seconds=900,peak_rss_kib=262144)
        for r in d[2]:r.update(elapsed_seconds=900,peak_rss_kib=262144)
    add('route_and_outer_exact_resource_boundary',True,resource_boundary)
    add('route_time_budget_overrun',False,lambda d:d[0]['resources'].__setitem__('seconds',900.000001))
    add('route_memory_budget_overrun',False,lambda d:d[0]['resources'].__setitem__('peak_rss_kib',262145))
    add('route_memory_zero',False,lambda d:d[0]['resources'].__setitem__('peak_rss_kib',0))
    add('route_resource_limit_changed',False,lambda d:d[0]['resources'].__setitem__('limit_seconds',901))
    add('outer_time_budget_overrun',False,lambda d:d[2][0].__setitem__('elapsed_seconds',900.000001))
    add('outer_memory_budget_overrun',False,lambda d:d[2][0].__setitem__('peak_rss_kib',262145))
    add('outer_resource_limit_changed',False,lambda d:d[2][0].__setitem__('memory_limit_kib',262145))
    add('outer_process_creation_limit_missing',False,lambda d:d[2][0].pop('process_creation_limit'))
    add('outer_process_creation_limit_changed',False,lambda d:d[2][0].__setitem__('process_creation_limit',1))
    add('outer_missing_independent',False,lambda d:d[2].pop())
    add('outer_duplicate_primary',False,lambda d:d[2].__setitem__(1,copy.deepcopy(d[2][0])))
    add('outer_routes_reordered',False,lambda d:d[2].reverse())
    add('outer_failed_status',False,lambda d:d[2][0].__setitem__('status','FAIL_EXECUTION'))
    add('outer_nonzero_exit',False,lambda d:d[2][0].__setitem__('exit_code',1))
    add('outer_boolean_exit',False,lambda d:d[2][0].__setitem__('exit_code',False))
    add('outer_timeout',False,lambda d:d[2][0].__setitem__('timed_out',True))
    add('outer_memory_exceeded',False,lambda d:d[2][0].__setitem__('memory_exceeded',True))
    add('outer_exception_present',False,lambda d:d[2][0].__setitem__('exception','synthetic'))
    # Timers/RSS measurements have different inner and outer scopes; budgets are enforced
    # separately. Their cross-inequalities are informational, not specified gates.
    add('separate_resource_scope_elapsed_allowed_within_budget',True,lambda d:d[0]['resources'].__setitem__('seconds',3.0))
    add('separate_resource_scope_rss_allowed_within_budget',True,lambda d:d[0]['resources'].__setitem__('peak_rss_kib',4096))
    for value,label in [(0.0,'float'),(True,'boolean'),('NaN','nan'),('Infinity','infinity'),
                        ('-Infinity','negative_infinity'),('1/2','fraction'),('garbage','text')]:
        add('scientific_scalar_'+label,False,lambda d,v=value:first(d[0]).__setitem__('I_ab',v))
    add('scientific_profile_nondecimal',False,lambda d:first(d[0])['P'].__setitem__(10,0.0))
    add('scientific_profile_nonfinite',False,lambda d:first(d[0])['P'].__setitem__(10,'NaN'))
    return tests

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=DEFAULT_SOURCE)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    before=hashlib.sha256(args.source.read_bytes()).hexdigest()
    # Prevent import from creating __pycache__ inside the sole-editor checkout.
    sys.dont_write_bytecode=True
    spec=importlib.util.spec_from_file_location('ledger_validator_synthetic',args.source)
    validator=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    checks=[]
    for name,expected,mutate in make_tests():
        fixture=synthetic_fixture()
        mutate(fixture)
        serialized=json.dumps(fixture,sort_keys=True,separators=(',',':')).encode()
        try:
            output=validator.cross_validate(*fixture)
            accepted=True
            detail={'classification':output['classification'],
                    'cross_route_max_gap_exact_rational':output['cross_route_max_gap_exact_rational']}
        except Exception as error:
            accepted=False
            detail={'error_type':type(error).__name__,'error':str(error)}
        checks.append({'name':name,'expected_accept':expected,'actual_accept':accepted,
                       'status':'PASS' if accepted == expected else 'FAIL',
                       'synthetic_fixture_sha256':hashlib.sha256(serialized).hexdigest(),**detail})
    after=hashlib.sha256(args.source.read_bytes()).hexdigest()
    failed=[x for x in checks if x['status'] != 'PASS']
    result={'schema_version':1,'scope':'SYNTHETIC_ONLY_cross_validate_SCHEMA_AND_MUTATIONS',
            'physical_files_loaded':False,'physical_route_executed':False,
            'source':str(args.source),'source_sha256_before':before,'source_sha256_after':after,
            'source_unchanged':before == after,
            'test_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'python_executable':sys.executable,'python_version':sys.version,
            'python_optimization':sys.flags.optimize,
            'status':'PASS' if not failed and before == after else 'BLOCKED_MUTATION_GAPS',
            'counts':{'total':len(checks),'pass':len(checks)-len(failed),'fail':len(failed)},
            'failed_checks':failed,'checks':checks}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'counts':result['counts'],'source_sha256':before,
                      'failures':[x['name'] for x in failed]},sort_keys=True))
    return 0 if result['status'] == 'PASS' else 1

if __name__ == '__main__':
    raise SystemExit(main())
