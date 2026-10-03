"""Narrow schema/guard controls on fabricated zero-source interval output.

This is not an integrator workload benchmark or a physical source evaluation.
"""
from __future__ import annotations
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

from operator_probe_contract import (
    CHANNELS,CONFIGURATIONS,MOMENTA,SCOPE,SOURCES,absolute_radius,encode,fraction,interval,
    load_json_no_duplicates,panel_interval,rat,validate_chronology,validate_payload,validate_pair)


def fabricated_row(source,momentum,panel=None,width=Q(1,2**120)):
    box = encode(-width,width)
    zero = encode(Q(0),Q(0))
    real_only = {"real":box,"imag":zero}
    complex_box = {"real":box,"imag":zero if momentum == 0 else box}
    moments = {"M0":copy.deepcopy(real_only),
               "Mexp":copy.deepcopy(complex_box),"Mu":copy.deepcopy(complex_box)}
    row = {"source":source,"momentum":rat(momentum),
           "interval":panel_interval(panel) if panel is not None else ["-9/2","-7/2"],
           "moments":moments,"total_absolute_radii":{k:rat(absolute_radius(v))
                                                      for k,v in moments.items()}}
    if panel is not None:
        row["panel"] = panel
    return row


def fixture():
    def work(source,panel=None):
        w=Q(1,2**120)
        moment={"real":encode(-w,w),"imag":encode(Q(0),Q(0))}
        row={"source":source,"interval":panel_interval(panel) if panel is not None else ["-9/2","-7/2"],
             "moment":moment,"total_absolute_radius":rat(w)}
        if panel is not None:row['panel']=panel
        return row
    return {"schema_version":1,"scope":SCOPE,"configuration":CONFIGURATIONS[0],
            "archive_arrays_decoded":0,
            "panel_rows":[fabricated_row(s,k,p) for s in SOURCES for p in range(64) for k in MOMENTA],
            "whole_rows":[fabricated_row(s,k) for s in SOURCES for k in MOMENTA],
            "source_work_panel_rows":[work(s,p) for s in SOURCES for p in range(64)],
            "source_work_whole_rows":[work(s) for s in SOURCES]}


def run():
    checks=[]
    def check(name,condition):
        if not condition:raise RuntimeError("Fabricated control failed: "+name)
        checks.append(name)
    def reject(name,callback):
        try:callback()
        except (ValueError,KeyError,TypeError):check(name,True)
        else:raise RuntimeError("Expected rejection: "+name)
    root=Path(__file__).parent
    registration=load_json_no_duplicates((root/'EXPERIMENT.json').read_text())
    validate_chronology(registration['chronology'])
    check('canonical_marker_matches_narrow_registration',True)
    check('all64panels_all9momenta_bothsources_fixed',
          registration['panel_count']==64 and len(registration['probe_momenta'])==9 and
          registration['sources']==list(SOURCES))
    check('full_twelve_case_integral_remains_uncomputed',
          registration['larger_study_status']['original_twelve_case_continuity_integral']=='UNCOMPUTED')
    check('resource_budget_fixed_before_real_source_evaluation',
          registration['resources']['each_route_wall_seconds']==120 and
          registration['resources']['each_route_peak_rss_kib']==131072)
    check('unverified_registration_cannot_authorize_real_source_evaluation',
          registration['freeze_authorization']=='NOT_AUTHORIZED_BEFORE_REMOTE_BYTE_VERIFICATION')

    counts={'fabricated_decode_stub':0,'fabricated_source_stub':0}
    def fake_entry(value):
        validate_chronology(value['chronology'])
        if value['archive_array_decoding_allowed'] is not False:
            raise ValueError('Archived arrays forbidden')
        counts['fabricated_decode_stub']+=1
        counts['fabricated_source_stub']+=1
        return 'FABRICATED_STUB_ONLY'
    renamed=copy.deepcopy(registration)
    marker=renamed['chronology'].pop('new_target_evaluations_before_public_verification')
    renamed['chronology']['physical_diagnostic_evaluations_before_freeze']=marker
    reject('renamed_execution_marker_rejected',lambda:fake_entry(renamed))
    check('renamed_marker_fails_before_fabricated_stubs',counts=={
        'fabricated_decode_stub':0,'fabricated_source_stub':0})
    for wrong,name in [(False,'boolean_zero'),('0','string_zero'),(1,'prior_evaluation'),(-1,'negative_count')]:
        bad=copy.deepcopy(registration['chronology'])
        bad['new_target_evaluations_before_public_verification']['primary']=wrong
        reject(name+'_rejected',lambda b=bad:validate_chronology(b))
    check('same_finalized_marker_passes_entry_contract',fake_entry(registration)=='FABRICATED_STUB_ONLY')

    payload=fixture()
    validated=validate_payload(payload)
    check('exact1170moment_and130work_rows_pass_validation',
          validated['status']=='PASS_REGISTERED_PROBE_OUTPUT_VALIDATION')
    check('contract_validation_does_not_claim_uniform_proof',
          validated['uniform_analytic_proof_validated_by_this_module'] is False)
    check('contract_validation_does_not_authorize_freeze',validated['production_freeze_authorized'] is False)
    for field in ('panel_rows','whole_rows','source_work_panel_rows','source_work_whole_rows'):
        bad=copy.deepcopy(payload);bad[field].pop()
        reject('missing_'+field+'_rejected',lambda b=bad:validate_payload(b))
        bad=copy.deepcopy(payload);bad[field][1]=copy.deepcopy(bad[field][0])
        reject('duplicate_'+field+'_rejected',lambda b=bad:validate_payload(b))
    bad=copy.deepcopy(payload);bad['panel_rows'][0]['panel']=False
    reject('boolean_panel_index_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload);bad['archive_arrays_decoded']=False
    reject('boolean_decode_count_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload);bad['archive_arrays_decoded']=1
    reject('any_archived_array_decode_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload);bad['whole_rows'][0]['interval']=['-6/1','-3/2']
    reject('changed_target_interval_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload);bad['panel_rows'][2]['momentum']='1/2048'
    reject('unregistered_momentum_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload);bad['whole_rows'][2]['total_absolute_radii']['Mexp']='0/1'
    reject('understated_complex_radius_rejected',lambda:validate_payload(bad))
    bad=copy.deepcopy(payload)
    wider=encode(-Q(1,2**80),Q(1,2**80))
    bad['whole_rows'][8]['moments']['Mexp']={'real':wider,'imag':wider}
    bad['whole_rows'][8]['total_absolute_radii']['Mexp']=rat(Q(1,2**79))
    check('valid_wide_output_is_unresolved',
          validate_payload(bad)['status']=='UNRESOLVED_PROBE_RADIUS')
    component_width=Q(1,2**87)
    box={'real':encode(-component_width,component_width),'imag':encode(-component_width,component_width)}
    check('complex_L1_radius_combines_both_components',absolute_radius(box)==2*component_width)
    # Each component is below1e-26, yet their total is above it.
    component_width=Q(1,2**87)
    check('component_gate_alone_can_understate_absolute_radius',
          component_width < Q(1,10**26) < 2*component_width)
    bad=copy.deepcopy(payload)
    shifted=encode(Q(1),Q(1))
    bad['whole_rows'][0]['moments']['Mexp']['real']=shifted
    bad['whole_rows'][0]['total_absolute_radii']['Mexp']='0/1'
    check('zero_momentum_identity_failure_retained',
          validate_payload(bad)['status']=='CERTIFICATE_CONSISTENCY_FAILURE')
    bad=copy.deepcopy(payload)
    bad['whole_rows'][0]['moments']['M0']['real']=shifted
    bad['whole_rows'][0]['total_absolute_radii']['M0']='0/1'
    check('incompatible_global_additivity_retained',
          validate_payload(bad)['status']=='CERTIFICATE_CONSISTENCY_FAILURE')
    check('general_exact_rational_endpoints_accepted',interval({'lo':'1/3','hi':'1/2'})==(Q(1,3),Q(1,2)))
    reject('midpoint_string_not_interval',lambda:interval('0 +/- 1e-26'))
    reject('inverted_interval_rejected',lambda:interval({'lo':'1/1','hi':'0/1'}))
    reject('noncanonical_rational_rejected',lambda:fraction('2/4'))
    reject('duplicate_JSON_rejected',lambda:load_json_no_duplicates('{"x":0,"x":1}'))
    reject('nonfinite_JSON_rejected',lambda:load_json_no_duplicates('{"x":Infinity}'))

    # Pure proof arithmetic, not real source evaluation. The mathematical
    # complex-domain premise M_g<=64 still requires independent proof review.
    from math import factorial
    source_bound=Q(1,15*2**90)
    phase_pointwise=Q(4096*8**97,factorial(97))
    kernel_pointwise=Q(2*4096*8**97,factorial(98))
    check('source_geometric_tail_exact',Q(64)*Q(1,16)**25/(1-Q(1,16))==source_bound)
    check('pointwise_phase_tail_is_below_1e60',phase_pointwise<Q(1,10**60))
    check('centered_kernel_integral_requires_factor2',
          2*phase_pointwise>phase_pointwise and 2*kernel_pointwise>kernel_pointwise)
    check('uniform_source_tail_below_narrow_radius_gate',source_bound<Q(1,10**26))
    serialized=json.dumps(payload,sort_keys=True,allow_nan=False)
    check('full_shape_exact_roundtrip',load_json_no_duplicates(serialized)==payload)
    independent=copy.deepcopy(payload);independent['configuration']=CONFIGURATIONS[1]
    both=validate_pair(payload,independent)
    check('both_full_route_output_shapes_validated',both['status']=='PASS_BOTH_REGISTERED_PROBE_OUTPUTS')
    check('both_route_validation_still_requires_outer_proof_authentication',
          both['outer_bootstrap_proof_and_resource_checks_required'] is True)
    bad_pair=copy.deepcopy(independent)
    bad_pair['whole_rows'][1]['moments']['Mexp']['real']=encode(Q(1),Q(1))
    bad_pair['whole_rows'][1]['total_absolute_radii']['Mexp']=rat(
        absolute_radius(bad_pair['whole_rows'][1]['moments']['Mexp']))
    check('disjoint_both_route_target_enclosures_fail',
          validate_pair(payload,bad_pair)['status']=='CERTIFICATE_CONSISTENCY_FAILURE')
    p_hull=copy.deepcopy(payload);i_hull=copy.deepcopy(independent)
    radius=Q(5,2**89)
    p_hull['whole_rows'][1]['moments']['Mexp']['real']=encode(-radius,radius)
    i_hull['whole_rows'][1]['moments']['Mexp']['real']=encode(Q(0),2*radius)
    for value in (p_hull,i_hull):
        value['whole_rows'][1]['total_absolute_radii']['Mexp']=rat(
            absolute_radius(value['whole_rows'][1]['moments']['Mexp']))
    check('overlap_intersection_does_not_replace_wider_cross_route_hull',
          validate_payload(p_hull)['status']=='PASS_REGISTERED_PROBE_OUTPUT_VALIDATION' and
          validate_payload(i_hull)['status']=='PASS_REGISTERED_PROBE_OUTPUT_VALIDATION' and
          validate_pair(p_hull,i_hull)['status']=='UNRESOLVED_PROBE_RADIUS')
    return {'status':'PASS_FABRICATED_NARROW_PROTOCOL_CONTRACT','checks':checks,
            'checks_passed':len(checks),'panel_rows':1152,'whole_rows':18,
            'source_work_panel_rows':128,'source_work_whole_rows':2,
            'fabricated_payload_sha256':hashlib.sha256(serialized.encode()).hexdigest(),
            'real_source_callbacks':0,'archive_arrays_decoded':0,
            'actual_integrator_whole_route_benchmark':False,
            'uniform_complex_domain_proof_established_by_this_fixture':False,
            'production_freeze_authorized':False,
            'scope':'Schema/guard/radius controls on fabricated zero-source output only.'}


if __name__=='__main__':
    result=run()
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(sys.argv)==2:Path(sys.argv[1]).write_text(raw)
    else:print(raw,end='')
