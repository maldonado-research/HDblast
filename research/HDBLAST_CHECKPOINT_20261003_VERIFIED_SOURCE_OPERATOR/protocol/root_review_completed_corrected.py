"""Strict stdlib review of archived guarded source/operator executions.

No source or numerical method is imported. Mathematical validity remains
conditional on the frozen analytic proofs and numerical library contracts.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import factorial, isfinite
from pathlib import Path
import re
import sys

from registration_guard import authenticate, require, load

SOURCES=('positive_B','signed_uB')
MOMENTA=('0/1','1/1099511627776','1/4096','1/4','1/1','16/1','64/1','128/1','256/1')
CHANNELS=('M0','Mexp','Mu')
SOURCE_TAIL=Q(1,15*2**90)
WORK_TAIL=SOURCE_TAIL/2
GATE=Q(1,10**26)
OUTER_KEYS={'status','route','command','uid','fabricated_only','exit_code','wait_status',
    'timed_out','wall_seconds','peak_rss_kib_wait4','cpu_seconds_wait4','outer_limits',
    'registration_sha256','remote_receipt_sha256','child_log_sha256',
    'source_attempt_journal_sha256','source_attempt_journal_bytes','output_sha256','output_bytes'}
ENVELOPE_KEYS={'schema_version','route','payload','method_evidence','semantic_validation','execution'}
EXECUTION_KEYS={'uid','python','optimization','wall_seconds_internal','peak_rss_kib_internal',
    'registration_sha256','remote_receipt_sha256','freeze_commit','source_chronology'}
MODE_BOUNDS={'source_ME_error','source_Mu_error','local_ME_defect_error','local_Mu_defect_error',
    'prefix_ME_defect_increment','prefix_Mu_defect_increment','inherited_w_radius',
    'inherited_u_radius','chosen_w_point_shift','chosen_u_point_shift',
    'w_radius_upward_rounding','u_radius_upward_rounding'}


def exact_int(value, minimum=0):
    return type(value) is int and value>=minimum


def measured(value, positive=True):
    return type(value) in (int,float) and isfinite(value) and (value>0 if positive else value>=0)


def rational(text):
    require(type(text) is str and len(text)<=10000 and
        re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*',text),
        'Canonical bounded rational evidence required')
    # Match the frozen child's explicit decimal integer conversion contract.
    # This does not relax outward rounding or allow unbounded output length.
    if sys.get_int_max_str_digits()<10000:
        sys.set_int_max_str_digits(10000)
    q=Q(text)
    require(text==str(q.numerator)+'/'+str(q.denominator),'Reduced rational evidence required')
    return q


def nonnegative(text):
    q=rational(text);require(q>=0,'Negative error/radius evidence');return q


def hexadecimal(value):
    return type(value) is str and re.fullmatch('[0-9a-f]{64}',value) is not None


def read_file(directory,name):
    path=directory/name
    require(path.is_file() and not path.is_symlink() and
        all(not p.is_symlink() for p in path.parents),'Missing/symlink execution artifact')
    require(path.stat().st_size<=64*1024*1024,'Oversized execution artifact')
    return path.read_bytes()


def canonical(q):
    return str(q.numerator)+'/'+str(q.denominator)


def expected_events():
    return [{'source':source,'center':str(Q(-9,2)+Q(2*j+1,128))}
        for source in SOURCES for j in range(64)]


def validate_chronology(envelope,outer,journal,route,registration_sha,receipt_sha,commit):
    execution=envelope['execution']
    require(type(execution) is dict and set(execution)==EXECUTION_KEYS,
        'Exact child execution schema required')
    require(exact_int(execution['uid'],1) and execution['uid']==outer['uid'] and
        type(execution['python']) is str and execution['python'].startswith('3.12.14 ') and
        type(execution['optimization']) is int and execution['optimization'] in (0,1),
        'Child uid/Python/optimization evidence differs')
    require(measured(execution['wall_seconds_internal']) and
        execution['wall_seconds_internal']<=outer['wall_seconds'] and
        exact_int(execution['peak_rss_kib_internal'],1) and
        execution['peak_rss_kib_internal']<=outer['peak_rss_kib_wait4'],
        'Child resource evidence differs from authoritative outer scope')
    require(execution['registration_sha256']==registration_sha and
        execution['remote_receipt_sha256']==receipt_sha and execution['freeze_commit']==commit,
        'Child registration/readback/freeze pins differ')
    chronology=execution['source_chronology']
    require(type(chronology) is dict and set(chronology)=={
        'route','source_bundle_constructions','archive_arrays_decoded','events'} and
        chronology['route']==route and type(chronology['source_bundle_constructions']) is int and
        chronology['source_bundle_constructions']==128 and
        type(chronology['archive_arrays_decoded']) is int and chronology['archive_arrays_decoded']==0,
        'Actual source chronology schema/count differs')
    events=expected_events()
    require(chronology['events']==events,'Exact ordered128 source event universe differs')
    require(journal.endswith(b'\n'),'Source attempt journal must finish a complete line')
    lines=journal.splitlines()
    require(len(lines)==128 and all(lines),'Exact128 source attempt records required')
    for index,(raw,event) in enumerate(zip(lines,events),1):
        row=load(raw)
        require(type(row) is dict and set(row)=={'index','route','source','center'} and
            type(row['index']) is int and row['index']==index and row['route']==route and
            row['source']==event['source'] and row['center']==event['center'],
            'Persisted source attempt journal differs from source chronology')


def analytic_model_expected():
    return {'source_degree':24,'phase_degree':96,'source_complex_disk_radius':'1/8',
        'panel_half_width':'1/128','source_majorant':'64/1','Lg_majorant':'32/1',
        'source_uniform_tail':canonical(SOURCE_TAIL),'Lg_uniform_tail':canonical(WORK_TAIL),
        'local_phase_cap':'4/1','integrand_phase_cap':'8/1','exponential_majorant':'4096/1',
        'phase_pointwise_E_tail':canonical(Q(4096*8**97,factorial(97))),
        'phase_pointwise_Q_tail':canonical(Q(2*4096*8**97,factorial(98))),
        'phase_integrated_E_tail':canonical(Q(2*4096*8**97,factorial(97))),
        'phase_integrated_Q_tail':canonical(Q(4*4096*8**97,factorial(98))),
        'zero_momentum_identity':'All source moments are real at k=0; exact imaginary zero is projected only for the real-source moment target, never for incoming complex states.'}


def primary_evidence(evidence,payload):
    required={'schema_version','configuration','panel_rows','whole_rows',
        'source_work_panel_rows','source_work_whole_rows','budget_semantics','whole_gate',
        'analytic_model','maximum_whole_complete_L1_radius','status'}
    require(type(evidence) is dict and set(evidence)==required,
        'Complete primary model/error evidence required')
    require(type(evidence['schema_version']) is int and evidence['schema_version']==1 and
        evidence['configuration']==payload['configuration']=='PRIMARY_ARB256_SOURCE24_PHASE96',
        'Primary evidence configuration/schema differs')
    expected=analytic_model_expected()
    require(evidence['analytic_model']==expected and
        type(evidence['analytic_model']['source_degree']) is int and
        type(evidence['analytic_model']['phase_degree']) is int,
        'Primary analytic source/phase model differs')
    require(type(evidence['budget_semantics']) is str and evidence['budget_semantics'] and
        evidence['whole_gate']=='1e-26 complete complex L1 radius from outward exact rational endpoints',
        'Primary error scope/gate differs')
    max_radius=Q(0)
    for field in ('panel_rows','whole_rows','source_work_panel_rows','source_work_whole_rows'):
        rows=evidence[field];p_rows=payload[field]
        require(type(rows) is list and len(rows)==len(p_rows),'Missing primary budget rows')
        moment_rows=field in ('panel_rows','whole_rows')
        whole='whole' in field
        for row,p_row in zip(rows,p_rows):
            metadata={'source','interval'}|({'panel'} if not whole else set())|({'momentum'} if moment_rows else set())
            budget_fields={'source_model_disk_radius','phase_model_disk_radius',
                'coefficient_and_arithmetic_baseline_L1_radius','complete_output_L1_radius'} if moment_rows else {
                'source_model_radius','phase_model_radius','coefficient_and_arithmetic_baseline_radius','complete_output_radius'}
            require(type(row) is dict and set(row)==metadata|budget_fields,
                'Primary budget row schema differs')
            require(all(row[k]==p_row[k] for k in metadata) and
                (whole or type(row['panel']) is int),'Primary budget row target/order differs')
            if moment_rows:
                for name in budget_fields:
                    require(type(row[name]) is dict and set(row[name])==set(CHANNELS),
                        'Primary complete moment budget inventory differs')
                    for value in row[name].values():nonnegative(value)
                require(row['complete_output_L1_radius']==p_row['total_absolute_radii'],
                    'Primary error evidence/output complete radii differ')
                width=Q(1) if whole else Q(1,64)
                required_source={'M0':width*SOURCE_TAIL,'Mexp':width*SOURCE_TAIL,
                    'Mu':width**2*SOURCE_TAIL/2}
                require(all(rational(row['source_model_disk_radius'][k])==v
                    for k,v in required_source.items()),'Primary source model propagation differs')
                require(row['phase_model_disk_radius']['M0']=='0/1' and
                    (row['momentum']!='0/1' or all(v=='0/1'
                        for v in row['phase_model_disk_radius'].values())),
                    'Primary zero-phase error convention differs')
                if whole:max_radius=max(max_radius,*(rational(v)
                    for v in row['complete_output_L1_radius'].values()))
            else:
                for name in budget_fields:nonnegative(row[name])
                require(row['complete_output_radius']==p_row['total_absolute_radius'] and
                    rational(row['source_model_radius'])==(WORK_TAIL if whole else WORK_TAIL/64) and
                    row['phase_model_radius']=='0/1','Primary work error/output target differs')
                if whole:max_radius=max(max_radius,rational(row['complete_output_radius']))
    require(nonnegative(evidence['maximum_whole_complete_L1_radius'])==max_radius,
        'Primary maximum radius differs from exact outputs')
    expected_status='PASS_CERTIFIED_OPERATOR_WIDTH_GATE' if max_radius<=GATE else 'UNRESOLVED_CERTIFICATE'
    require(evidence['status']==expected_status,'Primary width status differs from complete output radii')


def independent_evidence(evidence):
    required={'method','source_degree','mode_degree','scalar_exp_degree','dyadic_bits',
        'source_model_rows','mode_error_rows','physical_source_evaluations','fabricated_source_provider',
        'retained_arrays_decoded','primary_helper_imports','global_state_origin',
        'source_forcing_disk_bound','Lg_disk_bound','evidence_error_values',
        'exported_radius_definition','exported_endpoint_rounding_per_component_upper',
        'resource','construction_benchmark'}
    require(type(evidence) is dict and set(evidence)==required,
        'Complete independent source/defect evidence required')
    for name,value in [('source_degree',24),('mode_degree',96),('scalar_exp_degree',200),
            ('dyadic_bits',512),('physical_source_evaluations',128),
            ('retained_arrays_decoded',0),('primary_helper_imports',0)]:
        require(type(evidence[name]) is int and evidence[name]==value,
            'Independent degree/precision/source counter differs: '+name)
    require(evidence['method']=='independent exact polynomial ODE defect and rational source model' and
        evidence['fabricated_source_provider'] is False and evidence['construction_benchmark'] is None and
        evidence['global_state_origin']=='exact zero before the first source panel' and
        evidence['source_forcing_disk_bound']=='64/1' and evidence['Lg_disk_bound']=='32/1' and
        evidence['evidence_error_values']=='each error component is independently rounded upward on512bit grid' and
        evidence['exported_radius_definition']=='L1 sum of actual real/imag rectangle halfwidths' and
        evidence['exported_endpoint_rounding_per_component_upper']=='1/'+str(2**512),
        'Independent method/source bound/zero-anchor scope differs')
    internal=evidence['resource']
    require(type(internal) is dict and set(internal)=={'wall_seconds','peak_rss_kib'} and
        measured(internal['wall_seconds']) and exact_int(internal['peak_rss_kib'],1),
        'Independent internal resource metadata malformed')
    rows=evidence['source_model_rows']
    require(type(rows) is list and len(rows)==128,'Complete128 source model rows required')
    forcing_errors={}
    for index,(source,center) in enumerate((s,Q(-9,2)+Q(2*j+1,128))
            for s in SOURCES for j in range(64)):
        row=rows[index];panel=index%64
        require(type(row) is dict and set(row)=={'source','panel','center','forcing','Lg'} and
            row['source']==source and type(row['panel']) is int and row['panel']==panel and
            row['center']==canonical(center),'Independent source model row/order differs')
        for name,tail in [('forcing',SOURCE_TAIL),('Lg',WORK_TAIL)]:
            part=row[name]
            require(type(part) is dict and set(part)=={'analytic_tail','coefficient_error_uniform',
                'uniform_error','left_chosen_model_sha256'} and hexadecimal(part['left_chosen_model_sha256']),
                'Independent chosen source model metadata differs')
            require(nonnegative(part['analytic_tail'])==tail and
                nonnegative(part['uniform_error'])==tail+nonnegative(part['coefficient_error_uniform']),
                'Independent complete source error decomposition differs')
        forcing_errors[(source,panel)]=rational(row['forcing']['uniform_error'])
    rows=evidence['mode_error_rows']
    require(type(rows) is list and len(rows)==1152,'Complete1152 defect rows required')
    expected=((source,k,panel) for source in SOURCES for k in MOMENTA for panel in range(64))
    for row,(source,k,panel) in zip(rows,expected):
        require(type(row) is dict and set(row)==MODE_BOUNDS|{'source','panel','momentum'} and
            row['source']==source and type(row['panel']) is int and row['panel']==panel and
            row['momentum']==k,'Independent defect row/order differs')
        for name in MODE_BOUNDS:nonnegative(row[name])
        source_error=forcing_errors[(source,panel)]
        for name,target in [('source_ME_error',source_error/64),
                ('source_Mu_error',source_error/(2*64**2))]:
            reported=rational(row[name])
            require(target<=reported<target+Q(1,2**512),
                'Independent source-moment error rounding/propagation differs')
        if panel==0:require(row['inherited_w_radius']==row['inherited_u_radius']=='0/1',
            'Independent initial inherited radius is not exact zero')


def stable_frame(envelope):
    evidence=dict(envelope['method_evidence'])
    evidence.pop('resource',None)
    evidence.pop('construction_benchmark',None)
    return {'payload':envelope['payload'],'method_evidence':evidence,
        'semantic_validation':envelope['semantic_validation']}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def review(root,receipt,receipt_sha,registration_sha,run_dirs):
    verified=authenticate(root,receipt,receipt_sha,registration_sha)
    for artifact in (verified['receipt'],verified['registration']):
        for marker in ('fabricated_fixture','not_a_real_remote_freeze'):
            require(marker not in artifact or artifact[marker] is False,
                'Fabricated authentication metadata cannot support an actual source review')
    root=Path(root).absolute()
    sys.path.insert(0,str(root/'protocol'))
    from operator_probe_contract import validate_pair,validate_payload
    outputs={};receipts={}
    limits=verified['contract']['resources']
    require(exact_int(limits['each_route_wall_seconds'],1) and
        exact_int(limits['each_route_peak_rss_kib'],1),'Exact positive frozen resource limits required')
    for route in ('primary','independent'):
        directory=Path(run_dirs[route]).absolute()
        archived_mode=directory==root/'results'/route
        allowed_location=(root not in directory.parents and directory!=root) or archived_mode
        require(directory.is_dir() and directory.resolve()==directory and allowed_location,
            'Real external or exact archived route-results directory required')
        outer=load(read_file(directory,'EXECUTION.json'))
        raw=read_file(directory,'OUTPUT.json');envelope=load(raw)
        log=read_file(directory,'child.log');journal=read_file(directory,'SOURCE_ATTEMPTS.jsonl')
        require(type(outer) is dict and set(outer)==OUTER_KEYS and
            type(envelope) is dict and set(envelope)==ENVELOPE_KEYS,
            'Exact execution/envelope artifact schemas required')
        require(type(envelope['schema_version']) is int and envelope['schema_version']==1 and
            outer['status']=='PASS_BOUNDED_ROUTE_EXECUTION' and outer['fabricated_only'] is False and
            type(outer['exit_code']) is int and outer['exit_code']==0 and
            type(outer['wait_status']) is int and outer['wait_status']==0 and outer['timed_out'] is False and
            exact_int(outer['uid'],1) and outer['route']==envelope['route']==route,
            'Actual bounded execution state differs')
        require(outer['registration_sha256']==registration_sha and
            outer['remote_receipt_sha256']==receipt_sha and hexadecimal(outer['output_sha256']) and
            outer['output_sha256']==hashlib.sha256(raw).hexdigest() and
            exact_int(outer['output_bytes'],1) and outer['output_bytes']==len(raw) and
            hexadecimal(outer['child_log_sha256']) and outer['child_log_sha256']==hashlib.sha256(log).hexdigest() and
            hexadecimal(outer['source_attempt_journal_sha256']) and
            outer['source_attempt_journal_sha256']==hashlib.sha256(journal).hexdigest() and
            exact_int(outer['source_attempt_journal_bytes'],1) and outer['source_attempt_journal_bytes']==len(journal),
            'Execution output/log/source-journal hashes or sizes differ')
        require(measured(outer['wall_seconds']) and outer['wall_seconds']<=limits['each_route_wall_seconds'] and
            exact_int(outer['peak_rss_kib_wait4'],1) and outer['peak_rss_kib_wait4']<=limits['each_route_peak_rss_kib'] and
            measured(outer['cpu_seconds_wait4'],positive=False) and
            outer['outer_limits']=={'wall_seconds':limits['each_route_wall_seconds'],
                'peak_rss_kib':limits['each_route_peak_rss_kib']} and
            type(outer['outer_limits']['wall_seconds']) is int and
            type(outer['outer_limits']['peak_rss_kib']) is int,
            'Authoritative outer resource evidence differs')
        command=outer['command']
        require(type(command) is list and all(type(x) is str for x in command) and len(command) in (15,16),
            'Recorded bounded command malformed')
        expected=[str(root/'execution/run_registered.py'),'--root',str(root),'--receipt',
            str(Path(receipt).absolute()),'--receipt-sha256',receipt_sha,'--registration-sha256',
            registration_sha,'--route',route,'--output',str(directory/'OUTPUT.json')]
        optimized=command[1:3]==['-B','-O']
        require(command[1]=='-B' and len(command)==(16 if optimized else 15),
            'Recorded Python flag contract differs')
        arguments=command[3:] if optimized else command[2:]
        require(envelope['execution']['optimization']==(1 if optimized else 0),
            'Recorded command/child optimization differs')
        if archived_mode:
            original_root=Path(arguments[2]);original_output=Path(arguments[12])
            require(original_root.is_absolute() and Path(arguments[4]).is_absolute() and
                original_output.is_absolute() and original_output.name=='OUTPUT.json' and
                arguments[0]==str(original_root/'execution/run_registered.py') and
                [arguments[i] for i in (1,3,5,7,9,11)]==['--root','--receipt',
                    '--receipt-sha256','--registration-sha256','--route','--output'] and
                arguments[6]==receipt_sha and arguments[8]==registration_sha and arguments[10]==route,
                'Archived original command relationships/pins differ')
        else:
            require(arguments==expected,'Recorded fresh-external child command differs')
        validate_chronology(envelope,outer,journal,route,registration_sha,receipt_sha,
            verified['receipt']['freeze_commit'])
        numerical_route=validate_payload(envelope['payload'])
        require(envelope['semantic_validation']==numerical_route,
            'Recorded child semantic validation differs from recomputation')
        if route=='primary':primary_evidence(envelope['method_evidence'],envelope['payload'])
        else:independent_evidence(envelope['method_evidence'])
        outputs[route]=envelope
        receipts[route]={k:outer[k] for k in ('wall_seconds','peak_rss_kib_wait4','output_bytes',
            'output_sha256','child_log_sha256','source_attempt_journal_sha256','source_attempt_journal_bytes')}
    numerical=validate_pair(outputs['primary']['payload'],outputs['independent']['payload'])
    status=('PASS_UNIFORM_MODEL_AND_REGISTERED_PROBE_CERTIFICATE'
        if numerical['status']=='PASS_BOTH_REGISTERED_PROBE_OUTPUTS' else
        'CERTIFICATE_CONSISTENCY_FAILURE' if numerical['status']=='CERTIFICATE_CONSISTENCY_FAILURE' else
        'UNRESOLVED_CERTIFICATE')
    require(authenticate(root,receipt,receipt_sha,registration_sha)==verified,
        'Frozen source/readback/configuration changed during completed review')
    result={'schema_version':1,'status':status,'numerical_validation':numerical,
        'freeze_commit':verified['receipt']['freeze_commit'],'registration_sha256':registration_sha,
        'remote_receipt_sha256':receipt_sha,
        'row_counts_per_route':{'moment_panel':1152,'moment_whole':18,'Lg_panel':128,'Lg_whole':2},
        'resources':receipts,'stable_scientific_frame_sha256':{r:digest(stable_frame(v)) for r,v in outputs.items()},
        'uniform_scope':'Analytic source/panel/phase model bounds on64panels and real k∈[0,256]',
        'computed_width_scope':'All1300registered rows per route; nine exact rational momenta only',
        'proof_premises':['archived analytic Cauchy and Duhamel proofs','reviewed source/geometry identities',
            'Arb and exact-rational arithmetic contracts','custodian public freeze/readback chronology'],
        'limitations':{'twelve_case_action_pressure_contact_certificate':'UNRESOLVED',
            'metric_endpoint_29_and_refinement_30':'FAIL_PRESERVED',
            'retained_quantum_state_and_momentum_integration_error':'UNRESOLVED',
            'nonlinear_gravity_and_hot_radiation_era':'NOT_ESTABLISHED',
            'higher_dimensional_Big_Bang_cause':'NOT_ESTABLISHED',
            'external_mathematical_novelty':'NOT_ASSESSED','external_peer_review':'NONE'}}
    return result,outputs


def main():
    parser=argparse.ArgumentParser()
    for name in ('root','receipt','receipt-sha256','registration-sha256','primary','independent','output'):
        parser.add_argument('--'+name,required=True)
    args=parser.parse_args();output=Path(args.output)
    require(not output.exists() and not output.is_symlink(),'Fresh review output required')
    result,_=review(args.root,args.receipt,args.receipt_sha256,args.registration_sha256,
        {'primary':args.primary,'independent':args.independent})
    with output.open('x') as handle:
        json.dump(result,handle,sort_keys=True,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps(result,allow_nan=False),flush=True)


if __name__=='__main__':main()
