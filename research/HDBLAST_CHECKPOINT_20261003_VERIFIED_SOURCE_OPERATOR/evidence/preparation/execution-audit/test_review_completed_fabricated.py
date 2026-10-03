"""Adversarial stdlib-only audit of actual staged final result reviewer.

Authenticate is explicitly replaced by a fabricated stub. No real source,
numerical core, array reader, remote receipt or network operation is invoked.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

BASE=Path('/workspace/hdblast-research-work/verified-integration-20261003')
PROTOCOL=BASE/'protocol'
FORBIDDEN={'flint','numpy','scipy','sympy','mpmath','matplotlib','generic_operator','source_models','verified_moments','route'}


class NoNumericalImports:
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in FORBIDDEN:
            raise RuntimeError('Numerical/physical import forbidden in audit: '+fullname)
        return None


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def sha(raw):return hashlib.sha256(raw).hexdigest()
def dump(path,value):path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-directory',type=Path,required=True)
    args=parser.parse_args()
    out=args.output_directory.absolute()
    if out.exists():raise RuntimeError('Fresh audit output required')
    out.mkdir(parents=True)
    if FORBIDDEN.intersection(sys.modules):raise RuntimeError('Unexpected numerical module before audit')
    sys.meta_path.insert(0,NoNumericalImports())
    sys.dont_write_bytecode=True
    sys.path.insert(0,str(BASE/'execution-design'))
    sys.path.insert(0,str(PROTOCOL))
    fixture_module=module(PROTOCOL/'test_operator_contract_fabricated.py','fabricated_contract_provider')
    contract_module=module(PROTOCOL/'operator_probe_contract.py','fabricated_contract_reviewer')
    source=BASE/'execution-design/review_completed.py'
    rawsource=source.read_bytes()
    (out/'REVIEWER_SOURCE_SNAPSHOT.py').write_bytes(rawsource)
    reviewer=module(source,'actual_reviewer_under_fabricated_authentication')
    regsha='a'*64;receiptsha='b'*64;commit='c'*40
    contract=json.loads((PROTOCOL/'EXPERIMENT.json').read_text())
    verified={'contract':contract,'receipt':{'freeze_commit':commit}}
    reviewer.authenticate=lambda *a,**kw:copy.deepcopy(verified)
    root=out/'checkpoint';(root/'protocol').mkdir(parents=True)
    (root/'protocol/operator_probe_contract.py').write_bytes((PROTOCOL/'operator_probe_contract.py').read_bytes())
    payload=fixture_module.fixture()
    # A zero-source fixture with conservative exact rectangles. This is not
    # a source-model calculation or evidence of the true source's value.
    width=Q(1,2**90)
    for field in ('panel_rows','whole_rows'):
        for row in payload[field]:
            for name,box in row['moments'].items():
                box['real']=contract_module.encode(-width,width)
                if name!='M0' and row['momentum']!='0/1':box['imag']=contract_module.encode(-width,width)
                row['total_absolute_radii'][name]=contract_module.rat(contract_module.absolute_radius(box))
    for field in ('source_work_panel_rows','source_work_whole_rows'):
        for row in payload[field]:
            row['moment']['real']=contract_module.encode(-width,width)
            row['total_absolute_radius']=contract_module.rat(contract_module.absolute_radius(row['moment']))
    events=[{'source':s,'center':str(Q(-9,2)+Q(2*j+1,128))}
            for s in ('positive_B','signed_uB') for j in range(64)]
    base_envelopes={}
    for route in ('primary','independent'):
        value=copy.deepcopy(payload)
        if route=='independent':value['configuration']=contract_module.CONFIGURATIONS[1]
        evidence={'status':'PASS_CERTIFIED_OPERATOR_WIDTH_GATE'} if route=='primary' else {
            'physical_source_evaluations':128,'fabricated_source_provider':False,
            'retained_arrays_decoded':0,'primary_helper_imports':0}
        base_envelopes[route]={'schema_version':1,'route':route,'payload':value,
            'method_evidence':evidence,'semantic_validation':contract_module.validate_payload(value),
            'execution':{'uid':1000,'python':'3.12.14 FABRICATED ENVELOPE',
                'optimization':0,'wall_seconds_internal':1.0,'peak_rss_kib_internal':12000,
                'registration_sha256':regsha,'remote_receipt_sha256':receiptsha,
                'freeze_commit':commit,'source_chronology':{'route':route,
                    'source_bundle_constructions':128,'archive_arrays_decoded':0,'events':copy.deepcopy(events)}}}
    cases=[]
    def attempt(name,envelope_change=None,outer_change=None,log_change=None,expected_reject=True):
        case=out/name;case.mkdir()
        dirs={}
        for route in ('primary','independent'):
            directory=case/route;directory.mkdir();dirs[route]=directory
            envelope=copy.deepcopy(base_envelopes[route])
            if envelope_change:envelope_change(envelope,route)
            raw=json.dumps(envelope,sort_keys=True,allow_nan=False).encode()+b'\n'
            (directory/'OUTPUT.json').write_bytes(raw)
            logfile=directory/'child.log';logfile.write_text('FABRICATED CHILD LOG. NO COMMAND WAS EXECUTED.\n')
            outer={'status':'PASS_BOUNDED_ROUTE_EXECUTION','route':route,'uid':1000,
                'fabricated_only':False,'exit_code':0,'wait_status':0,'timed_out':False,
                'wall_seconds':1.0,'peak_rss_kib_wait4':12000,'cpu_seconds_wait4':0.5,
                'outer_limits':{'wall_seconds':120,'peak_rss_kib':131072},
                'registration_sha256':regsha,'remote_receipt_sha256':receiptsha,
                'child_log_sha256':sha(logfile.read_bytes()),'output_sha256':sha(raw),'output_bytes':len(raw)}
            if outer_change:outer_change(outer,route)
            dump(directory/'EXECUTION.json',outer)
            if log_change:log_change(logfile,route)
        try:
            result,_=reviewer.review(root,'FABRICATED_RECEIPT_NOT_READ',receiptsha,regsha,dirs)
        except (ValueError,KeyError,TypeError,FileNotFoundError) as error:
            record={'case':name,'accepted':False,'diagnostic':str(error),'expected_reject':expected_reject}
        else:
            dump(case/'REVIEW_RESULT.json',result)
            record={'case':name,'accepted':True,'review_status':result['status'],'expected_reject':expected_reject}
        record['expectation_met']=record['accepted'] is not expected_reject
        cases.append(record)

    attempt('baseline_minimal_evidence',expected_reject=False)
    attempt('negative_wall',outer_change=lambda o,r:o.update(wall_seconds=-1))
    attempt('negative_rss',outer_change=lambda o,r:o.update(peak_rss_kib_wait4=-1))
    attempt('boolean_uid',outer_change=lambda o,r:o.update(uid=True))
    attempt('boolean_archive_count',envelope_change=lambda e,r:e['execution']['source_chronology'].update(archive_arrays_decoded=False))
    attempt('arbitrary_128_events',envelope_change=lambda e,r:e['execution']['source_chronology'].update(events=[{}]*128))
    attempt('wrong_inner_registration_pin',envelope_change=lambda e,r:e['execution'].update(registration_sha256='d'*64))
    attempt('wrong_inner_receipt_pin',envelope_change=lambda e,r:e['execution'].update(remote_receipt_sha256='d'*64))
    attempt('stale_semantic_validation',envelope_change=lambda e,r:e.update(semantic_validation={'status':'FAIL_EXECUTION'}))
    attempt('changed_child_log',log_change=lambda p,r:p.write_text('MUTATED FABRICATED CHILD LOG\n'))
    attempt('missing_primary_model_error_rows',envelope_change=lambda e,r:e.update(method_evidence={'status':'PASS_CERTIFIED_OPERATOR_WIDTH_GATE'}) if r=='primary' else None)
    attempt('missing_independent_model_error_rows',envelope_change=lambda e,r:e.update(method_evidence={'physical_source_evaluations':128,'fabricated_source_provider':False,'retained_arrays_decoded':0,'primary_helper_imports':0}) if r=='independent' else None)
    attempt('unresolved_method_evidence_ignored',envelope_change=lambda e,r:e['method_evidence'].update(status='UNRESOLVED_CERTIFICATE') if r=='primary' else None)
    attempt('malformed_payload_missing_row',envelope_change=lambda e,r:e['payload']['whole_rows'].pop())
    attempt('known_fabricated_receipt_rejected',outer_change=lambda o,r:o.update(fabricated_only=True))
    forbidden_loaded=sorted(FORBIDDEN.intersection(sys.modules))
    if forbidden_loaded:raise RuntimeError('Numerical module imported during fabricated audit')
    report={'status':'FABRICATED_ADVERSARIAL_REVIEW_AUDIT','reviewer_sha256':sha(rawsource),
        'authentication_replaced_with_explicit_fabricated_stub':True,'real_source_calls':0,
        'archived_array_decodes':0,'numerical_imports':forbidden_loaded,
        'cases':cases,'mutation_acceptance_gaps':[r['case'] for r in cases
            if r['expected_reject'] and r['accepted']],
        'warning':'Fake envelopes and fake authentication are never scientific execution or remote-freeze evidence.'}
    dump(out/'AUDIT_RESULT.json',report)
    print(json.dumps({'reviewer_sha256':report['reviewer_sha256'],'cases':len(cases),
        'acceptance_gaps':report['mutation_acceptance_gaps'],'numerical_imports':forbidden_loaded},sort_keys=True))


if __name__=='__main__':main()
