"""Strict completed-review tests; fake authentication and fake envelopes only."""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

BASE=Path('/workspace/hdblast-research-work/verified-integration-20261003')
P=BASE/'protocol'
FORBIDDEN={'flint','numpy','scipy','sympy','mpmath','matplotlib',
    'generic_operator','source_models','verified_moments','route'}


class NoNumericalImports:
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in FORBIDDEN:raise RuntimeError('Numerical/source import forbidden')
        return None


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def sha(raw):return hashlib.sha256(raw).hexdigest()
def canonical(q):q=Q(q);return str(q.numerator)+'/'+str(q.denominator)
def dump(path,value):path.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')
def ceil512(q):return Q(-((-Q(q)*2**512)//1),2**512)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args();out=args.output_directory.absolute()
    if out.exists():raise RuntimeError('Fresh strict review test directory required')
    out.mkdir(parents=True);sys.dont_write_bytecode=True
    if FORBIDDEN.intersection(sys.modules):raise RuntimeError('Numerical imports before fixture')
    sys.meta_path.insert(0,NoNumericalImports())
    sys.path.insert(0,str(P));sys.path.insert(0,str(BASE/'execution-design'))
    source=P/'root_review_completed_corrected.py';reviewer=module(source,'corrected_fake_only_reviewer')
    contractmod=module(P/'operator_probe_contract.py','strict_fixture_contract')
    fixturemod=module(P/'test_operator_contract_fabricated.py','strict_fixture_provider')
    contract=json.loads((P/'EXPERIMENT.json').read_text())
    root=out/'checkpoint';(root/'protocol').mkdir(parents=True)
    (root/'protocol/operator_probe_contract.py').write_bytes((P/'operator_probe_contract.py').read_bytes())
    regsha='a'*64;receiptsha='b'*64;commit='c'*40;receipt=out/'FAKE_READBACK_NOT_READ.json'
    verified={'contract':contract,'registration':{},'receipt':{'freeze_commit':commit}}
    reviewer.authenticate=lambda *a,**kw:copy.deepcopy(verified)
    frame=fixturemod.fixture();width=Q(1,2**90)
    for field in ('panel_rows','whole_rows'):
        for row in frame[field]:
            for name,box in row['moments'].items():
                box['real']=contractmod.encode(-width,width)
                if name!='M0' and row['momentum']!='0/1':box['imag']=contractmod.encode(-width,width)
                row['total_absolute_radii'][name]=canonical(contractmod.absolute_radius(box))
    for field in ('source_work_panel_rows','source_work_whole_rows'):
        for row in frame[field]:
            row['moment']['real']=contractmod.encode(-width,width)
            row['total_absolute_radius']=canonical(contractmod.absolute_radius(row['moment']))

    def primary_budget(payload):
        evidence={'schema_version':1,'configuration':payload['configuration'],
            'analytic_model':reviewer.analytic_model_expected(),
            'budget_semantics':'Explicit fabricated schema fixture; no computation of the actual source.',
            'whole_gate':'1e-26 complete complex L1 radius from outward exact rational endpoints'}
        maximum=Q(0)
        for field in ('panel_rows','whole_rows','source_work_panel_rows','source_work_whole_rows'):
            whole='whole' in field;moment=field in ('panel_rows','whole_rows');rows=[]
            for p_row in payload[field]:
                metadata={k:v for k,v in p_row.items() if k not in ('moments','total_absolute_radii','moment','total_absolute_radius')}
                if moment:
                    w=Q(1) if whole else Q(1,64)
                    r=reviewer.SOURCE_TAIL
                    row={**metadata,'source_model_disk_radius':{'M0':canonical(w*r),'Mexp':canonical(w*r),'Mu':canonical(w*w*r/2)},
                        'phase_model_disk_radius':{k:'0/1' for k in reviewer.CHANNELS},
                        'coefficient_and_arithmetic_baseline_L1_radius':{k:'0/1' for k in reviewer.CHANNELS},
                        'complete_output_L1_radius':copy.deepcopy(p_row['total_absolute_radii'])}
                    if whole:maximum=max(maximum,*(Q(v) for v in p_row['total_absolute_radii'].values()))
                else:
                    row={**metadata,'source_model_radius':canonical(reviewer.WORK_TAIL if whole else reviewer.WORK_TAIL/64),
                        'phase_model_radius':'0/1','coefficient_and_arithmetic_baseline_radius':'0/1',
                        'complete_output_radius':p_row['total_absolute_radius']}
                    if whole:maximum=max(maximum,Q(p_row['total_absolute_radius']))
                rows.append(row)
            evidence[field]=rows
        evidence['maximum_whole_complete_L1_radius']=canonical(maximum)
        evidence['status']='PASS_CERTIFIED_OPERATOR_WIDTH_GATE' if maximum<=reviewer.GATE else 'UNRESOLVED_CERTIFICATE'
        return evidence

    def independent_evidence():
        rows=[]
        for source in reviewer.SOURCES:
            for panel in range(64):
                center=Q(-9,2)+Q(2*panel+1,128)
                row={'source':source,'panel':panel,'center':canonical(center)}
                for name,tail in [('forcing',reviewer.SOURCE_TAIL),('Lg',reviewer.WORK_TAIL)]:
                    row[name]={'analytic_tail':canonical(tail),'coefficient_error_uniform':'0/1',
                        'uniform_error':canonical(tail),'left_chosen_model_sha256':'f'*64}
                rows.append(row)
        errors=[]
        for source in reviewer.SOURCES:
            for k in reviewer.MOMENTA:
                for panel in range(64):
                    row={name:'0/1' for name in reviewer.MODE_BOUNDS}
                    row.update(source=source,momentum=k,panel=panel,
                        source_ME_error=canonical(ceil512(reviewer.SOURCE_TAIL/64)),
                        source_Mu_error=canonical(ceil512(reviewer.SOURCE_TAIL/(2*64**2))))
                    errors.append(row)
        return {'method':'independent exact polynomial ODE defect and rational source model',
            'source_degree':24,'mode_degree':96,'scalar_exp_degree':200,'dyadic_bits':512,
            'source_model_rows':rows,'mode_error_rows':errors,'physical_source_evaluations':128,
            'fabricated_source_provider':False,'retained_arrays_decoded':0,'primary_helper_imports':0,
            'global_state_origin':'exact zero before the first source panel','source_forcing_disk_bound':'64/1',
            'Lg_disk_bound':'32/1','evidence_error_values':'each error component is independently rounded upward on512bit grid',
            'exported_radius_definition':'L1 sum of actual real/imag rectangle halfwidths',
            'exported_endpoint_rounding_per_component_upper':'1/'+str(2**512),
            'resource':{'wall_seconds':0.5,'peak_rss_kib':10000},'construction_benchmark':None}

    base={}
    for route in ('primary','independent'):
        payload=copy.deepcopy(frame)
        if route=='independent':payload['configuration']=contractmod.CONFIGURATIONS[1]
        base[route]={'schema_version':1,'route':route,'payload':payload,
            'method_evidence':primary_budget(payload) if route=='primary' else independent_evidence(),
            'semantic_validation':contractmod.validate_payload(payload),
            'execution':{'uid':1000,'python':'3.12.14 EXPLICIT FAKE ENVELOPE','optimization':0,
                'wall_seconds_internal':0.75,'peak_rss_kib_internal':11000,
                'registration_sha256':regsha,'remote_receipt_sha256':receiptsha,'freeze_commit':commit,
                'source_chronology':{'route':route,'source_bundle_constructions':128,
                    'archive_arrays_decoded':0,'events':reviewer.expected_events()}}}
    records=[]
    def attempt(name,change=None,outer_change=None,journal_change=None,log_change=None,
            verified_change=None,expected='REJECT',archived=False):
        original_verified=copy.deepcopy(verified)
        if verified_change:verified_change(verified)
        case=out/name;case.mkdir();dirs={}
        for route in reviewer.SOURCES[:0]+('primary','independent'):
            directory=(root/'results'/route) if archived else case/route
            directory.mkdir(parents=True)
            dirs[route]=directory;envelope=copy.deepcopy(base[route])
            if change:change(envelope,route)
            raw=json.dumps(envelope,sort_keys=True,allow_nan=False).encode()+b'\n'
            (directory/'OUTPUT.json').write_bytes(raw)
            log=directory/'child.log';log.write_text('EXPLICIT FABRICATED CHILD LOG; NO EXECUTION.\n')
            journal=directory/'SOURCE_ATTEMPTS.jsonl'
            lines=[json.dumps({'index':i,'route':route,**event},sort_keys=True) for i,event in enumerate(reviewer.expected_events(),1)]
            journal.write_text('\n'.join(lines)+'\n')
            recorded_root=Path('/previous-host/frozen-checkpoint') if archived else root
            recorded_receipt=Path('/previous-host/remote-readback.json') if archived else receipt
            recorded_output=Path('/previous-host/raw-run')/route/'OUTPUT.json' if archived else directory/'OUTPUT.json'
            command=['/usr/local/bin/python','-B',str(recorded_root/'execution/run_registered.py'),
                '--root',str(recorded_root),'--receipt',str(recorded_receipt),'--receipt-sha256',receiptsha,
                '--registration-sha256',regsha,'--route',route,'--output',str(recorded_output)]
            outer={'status':'PASS_BOUNDED_ROUTE_EXECUTION','route':route,'command':command,'uid':1000,
                'fabricated_only':False,'exit_code':0,'wait_status':0,'timed_out':False,
                'wall_seconds':1.0,'peak_rss_kib_wait4':12000,'cpu_seconds_wait4':0.5,
                'outer_limits':{'wall_seconds':120,'peak_rss_kib':131072},
                'registration_sha256':regsha,'remote_receipt_sha256':receiptsha,
                'child_log_sha256':sha(log.read_bytes()),'source_attempt_journal_sha256':sha(journal.read_bytes()),
                'source_attempt_journal_bytes':journal.stat().st_size,'output_sha256':sha(raw),'output_bytes':len(raw)}
            if outer_change:outer_change(outer,route)
            dump(directory/'EXECUTION.json',outer)
            if journal_change:journal_change(journal,route)
            if log_change:log_change(log,route)
        try:result,_=reviewer.review(root,receipt,receiptsha,regsha,dirs)
        except (ValueError,KeyError,TypeError,FileNotFoundError) as error:
            observed='REJECT';diagnostic=str(error)
        else:
            observed=result['status'];diagnostic=None;dump(case/'REVIEW_RESULT.json',result)
        verified.clear();verified.update(original_verified)
        record={'case':name,'observed':observed,'expected':expected,'diagnostic':diagnostic,'passed':observed==expected}
        records.append(record)
        if observed!=expected:raise RuntimeError('Strict reviewer fixture mismatch: '+json.dumps(record))

    passed='PASS_UNIFORM_MODEL_AND_REGISTERED_PROBE_CERTIFICATE'
    attempt('valid_complete_fresh_envelopes',expected=passed)
    for name,field,value in [('negative_wall','wall_seconds',-1),('negative_rss','peak_rss_kib_wait4',-1),
        ('boolean_uid','uid',True),('boolean_exit','exit_code',False),('boolean_wait','wait_status',False),
        ('cpu_boolean','cpu_seconds_wait4',False),('outer_infinite_wall','wall_seconds',float('inf'))]:
        # Infinity is passed in memory to the JSON writer below; separately
        # test measured() without emitting a prohibited nonfinite JSON artifact.
        if name=='outer_infinite_wall':
            if reviewer.measured(value):raise RuntimeError('Infinity resource survived')
            records.append({'case':name,'observed':'REJECT','expected':'REJECT','passed':True});continue
        attempt(name,outer_change=lambda o,r,f=field,v=value:o.update({f:v}))
    attempt('arbitrary128_events',change=lambda e,r:e['execution']['source_chronology'].update(events=[{}]*128))
    attempt('boolean_archive_count',change=lambda e,r:e['execution']['source_chronology'].update(archive_arrays_decoded=False))
    attempt('inner_registration_pin',change=lambda e,r:e['execution'].update(registration_sha256='d'*64))
    attempt('inner_receipt_pin',change=lambda e,r:e['execution'].update(remote_receipt_sha256='d'*64))
    attempt('inner_freeze_pin',change=lambda e,r:e['execution'].update(freeze_commit='d'*40))
    attempt('stale_semantic_verdict',change=lambda e,r:e.update(semantic_validation={'status':'FAIL_EXECUTION'}))
    attempt('log_hash_mismatch',log_change=lambda p,r:p.write_text('MUTATED FAKE LOG\n'))
    attempt('journal_hash_mismatch',journal_change=lambda p,r:p.write_text('MUTATED FAKE JOURNAL\n'))
    attempt('missing_primary_model',change=lambda e,r:e['method_evidence'].pop('analytic_model') if r=='primary' else None)
    attempt('missing_primary_budget_row',change=lambda e,r:e['method_evidence']['whole_rows'].pop() if r=='primary' else None)
    attempt('primary_negative_budget',change=lambda e,r:e['method_evidence']['whole_rows'][0]['phase_model_disk_radius'].update(Mu='-1/1') if r=='primary' else None)
    attempt('primary_radius_mismatch',change=lambda e,r:e['method_evidence']['whole_rows'][0]['complete_output_L1_radius'].update(Mu='0/1') if r=='primary' else None)
    attempt('primary_stale_width_status',change=lambda e,r:e['method_evidence'].update(status='UNRESOLVED_CERTIFICATE') if r=='primary' else None)
    attempt('missing_independent_source_row',change=lambda e,r:e['method_evidence']['source_model_rows'].pop() if r=='independent' else None)
    attempt('missing_independent_defect_row',change=lambda e,r:e['method_evidence']['mode_error_rows'].pop() if r=='independent' else None)
    attempt('independent_negative_bound',change=lambda e,r:e['method_evidence']['mode_error_rows'][0].update(local_ME_defect_error='-1/1') if r=='independent' else None)
    attempt('independent_noncanonical_bound',change=lambda e,r:e['method_evidence']['mode_error_rows'][0].update(local_ME_defect_error='2/4') if r=='independent' else None)
    attempt('independent_understated_source_bound',change=lambda e,r:e['method_evidence']['mode_error_rows'][0].update(source_ME_error='0/1') if r=='independent' else None)
    attempt('independent_wrong_precision',change=lambda e,r:e['method_evidence'].update(dyadic_bits=256) if r=='independent' else None)
    attempt('independent_nonzero_initial_anchor',change=lambda e,r:e['method_evidence']['mode_error_rows'][0].update(inherited_w_radius='1/1') if r=='independent' else None)
    attempt('independent_fake_provider',change=lambda e,r:e['method_evidence'].update(fabricated_source_provider=True) if r=='independent' else None)
    attempt('receipt_fake_marker',verified_change=lambda v:v['receipt'].update(fabricated_fixture=True))
    attempt('registration_fake_marker',verified_change=lambda v:v['registration'].update(not_a_real_remote_freeze=True))
    attempt('outer_fake_only',outer_change=lambda o,r:o.update(fabricated_only=True))
    attempt('wrong_recorded_route',outer_change=lambda o,r:o['command'].__setitem__(-3,'other'))
    attempt('command_optimization_mismatch',outer_change=lambda o,r:o['command'].insert(2,'-O'))
    def wide(e,route):
        row=e['payload']['whole_rows'][1];width=Q(1,2**80)
        row['moments']['Mexp']['real']=contractmod.encode(-width,width)
        row['total_absolute_radii']['Mexp']=canonical(contractmod.absolute_radius(row['moments']['Mexp']))
        e['semantic_validation']=contractmod.validate_payload(e['payload'])
        if route=='primary':e['method_evidence']=primary_budget(e['payload'])
    attempt('genuine_valid_wide_scientific_outcome',change=wide,expected='UNRESOLVED_CERTIFICATE')
    def disjoint(e,route):
        if route!='independent':return
        row=e['payload']['whole_rows'][1];row['moments']['Mexp']['real']=contractmod.encode(Q(1),Q(1))
        row['total_absolute_radii']['Mexp']=canonical(contractmod.absolute_radius(row['moments']['Mexp']))
        e['semantic_validation']=contractmod.validate_payload(e['payload'])
    attempt('genuine_valid_discordant_scientific_outcome',change=disjoint,expected='CERTIFICATE_CONSISTENCY_FAILURE')
    attempt('archived_original_absolute_paths_preserved',archived=True,expected=passed)
    if FORBIDDEN.intersection(sys.modules):raise RuntimeError('Numerical import occurred')
    report={'status':'PASS_STRICT_COMPLETED_REVIEW_FABRICATED_CONTROLS','checks_passed':len(records),
        'records':records,'reviewer_sha256':sha(source.read_bytes()),
        'authentication_replaced_with_explicit_fabricated_stub':True,
        'real_source_calls':0,'archived_array_decodes':0,'numerical_imports':[],
        'remote_operations':0,'warning':'All envelopes and authorization are fabricated; these are not real-source or remote-freeze results.'}
    dump(out/'TEST_RESULT.json',report)
    print(json.dumps({'checks':len(records),'reviewer_sha256':report['reviewer_sha256'],'status':report['status']},sort_keys=True))


if __name__=='__main__':main()
