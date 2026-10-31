#!/usr/bin/env python3
"""Separate exact Cauchy-packet conditional response; exact rational ledger."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    p=Path(path)
    need(p.is_file() and not p.is_symlink() and all(not a.is_symlink() for a in p.absolute().parents),'regular symlink-free source: '+str(p))
    return hashlib.sha256(p.read_bytes()).hexdigest()

def encode(o):
    if isinstance(o,Q): return str(o)
    if isinstance(o,dict): return {k:encode(v) for k,v in o.items()}
    if isinstance(o,(list,tuple)): return [encode(v) for v in o]
    return o

def read(p): return json.loads(Path(p).read_text())

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    for role in ('same-response','absolute','velocity','independent-velocity'):
        ap.add_argument('--'+role+'-dir',type=Path,required=True)
    ap.add_argument('--freeze',type=Path,default=Path(__file__).with_name('FREEZE.json'))
    ap.add_argument('--freeze-sha256',required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--fault',choices=['reuse_old_velocity_absolute','reuse_old_velocity_errors','omit_f00','omit_r00','omit_full_Robin_source','omit_original_Fx','omit_derivative_commutator','lose_r1_initial_cancellation','omit_feedback_inverse'])
    ap.add_argument('--verify-report',type=Path)
    args=ap.parse_args()
    need(not args.output.exists(),'fresh output required')
    need(sha(args.freeze)==args.freeze_sha256,'frozen Cauchy response inputs/recipe')
    freeze=read(args.freeze)
    need(sha(__file__)==freeze['adapter_sha256'],'frozen Cauchy recipe source')
    roots={role:getattr(args,role+'_dir') for role in ('same_response','absolute','velocity','independent_velocity')}
    before={str(Path(__file__).resolve()):sha(__file__),str(args.freeze.resolve()):args.freeze_sha256}
    portable={'response_cauchy.py':sha(__file__),'FREEZE.json':args.freeze_sha256}
    def pin(p,expected,label):
        need(sha(p)==expected,'selected source/report/table pin: '+label)
        before[str(p.resolve())]=expected; portable[label]=expected
    for role,pins in freeze['input_pins'].items():
        for rel,expected in pins.items(): pin(roots[role]/rel,expected,role+'/'+rel)
    original_freeze=read(roots['same_response']/'checks/FREEZE.json')
    for rel,expected in original_freeze['source_pins'].items(): pin(args.repo_root/rel,expected,rel)
    native=read(roots['velocity']/'INITIAL_NORMAL.json')
    independent=read(roots['independent_velocity']/'INITIAL_NORMAL.json')
    regularity=read(roots['independent_velocity']/'REGULARITY_FINAL_NORMAL.json')
    native_audit=read(roots['independent_velocity']/'AUDIT_PORTABLE_FINAL_NORMAL.json')
    # Independently read back every selected original source authenticated by
    # either velocity proof before importing the local rational-only helper.
    for report,key in [(native,'source_authentication_before_after'),(independent,'source_sha256_before_and_after'),(regularity,'source_sha256_before_and_after')]:
        for rel,expected in report[key].items():
            need(rel.startswith('research/') and '..' not in Path(rel).parts,'portable original source path')
            pin(args.repo_root/rel,expected,rel)
    helper=roots['same_response']/'checks/response_adapter.py'
    spec=importlib.util.spec_from_file_location('frozen_rational_response_recipe',helper)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    reports={k:read(args.repo_root/rel) for k,rel in m.SOURCES.items()}
    published=reports['baseline']; old_campaign=reports['initial']['campaigns'][0]
    old_data=m.caps(old_campaign['candidate_initial_absolute_H_caps'],m.CAPS)
    old_error=m.caps(old_campaign['initial_error_H_caps'],m.ERRORS)
    k=m.constants(reports)
    baseline=m.reconstruct(k,old_data,old_error)
    nodes=m.exact_equal(baseline,published['ledger'],'all PR109 baseline ledger nodes')
    m.exact_equal(k,published['constants'],'all unchanged constants')
    need(native['candidate']==independent['candidate']==regularity['candidate']=='W1_exact_Cauchy_packet','distinct new candidate identity')
    need(native['target_F_Fx_G_M_H_domain_unchanged'] is True and native['old_candidate_functions_unchanged'] is True,'unchanged target and older candidates')
    passed=[c for c in native['campaigns'] if c['status']=='INITIAL_AND_LOCAL_GATES_PASS']
    need(len(passed)==1 and native['status']=='INITIAL_AND_LOCAL_GATES_PASS','native complete initial/local PASS')
    c=passed[0]
    need(c['initial_gate_pass'] is True and c['unresolved_cells']==[],'native initial complete coverage')
    need(c['labels']==published['labels'] and c['Gram']==published['Gram'],'native full23 modal norm')
    local=c['local']; low=Q(local['actual_H_norm_interval'][0]); upper=Q(local['new_error_H_upper']); oldlow=Q(local['old_error_H_lower'])
    need(local['point']==['0','65/64'] and local['observable']=='initial_velocity_r_full23_H' and local['old_point_table_reused'] is False,'fresh velocity local role')
    need(low>0 and 10*upper<=low and 2*upper<oldlow and local['pass'] is True,'native fresh local true-error inequalities')
    need(independent['status']=='PASS_INDEPENDENT_EXACT_VELOCITY_INITIAL_GATES' and all(independent['coverage']['initial_gates_pass']) and independent['local']['gates_pass'] is True,'independent complete initial/local PASS')
    il=independent['local']; need(Q(il['actual_v_r_H_interval'][0])>0 and 10*Q(il['new_error_v_r_H_interval'][1])<=Q(il['actual_v_r_H_interval'][0]) and 2*Q(il['new_error_v_r_H_interval'][1])<Q(il['old_error_v_r_H_interval'][0]),'independent fresh local inequalities')
    need(independent['coverage']['exact_zero_region']==[0,1],'flat support collar')
    need(native_audit['status']=='PASS_INDEPENDENT_VELOCITY_REPORT_AND_PROVENANCE_AUDIT','independent native-table audit')
    need(regularity['status']=='PASS_INDEPENDENT_ANALYTIC_REGULARITY_AUDIT' and regularity['required_q_order_for_analytic_velocity_r']==4,'analytic velocity derivative inventory')
    labels=published['labels']
    traces=c['endpoint_componentwise_affine_matching']
    need([row['label'] for row in traces]==labels and all(row['exterior_minus_interior_exact_jet']==[0,0] for row in traces),'all23 exact completed velocity value/r outer traces')
    independent_traces=regularity['complete_component_velocity_affine_extension_jet_differences_through1']
    need(set(independent_traces)==set(labels) and all(jets==['0','0'] for jets in independent_traces.values()),'independent C1 completed velocity outer traces')
    need(set(regularity['velocity_parity'])==set(labels) and all(abs(v)==1 for v in regularity['velocity_parity'].values()),'H-isometric all23 parity')
    for pair in regularity['nominal_denominator_enclosures_on_0_to_7over4'].values():
        lo,hi=map(Q,pair); need(lo<=hi and not lo<=0<=hi,'regular nominal denominators')
    native_data=m.caps(c['candidate_initial_absolute_H_caps'],m.CAPS)
    error=m.caps(c['initial_error_H_caps'],m.ERRORS)
    refined=read(roots['absolute']/'REFINEMENT.json')
    data=m.caps(refined['candidate_initial_absolute_H_caps'],m.CAPS)
    need(m.caps(refined['initial_error_H_caps'],m.ERRORS)==old_error,'absolute refinement target errors remain frozen')
    need(all(native_data[key]==old_data[key] and data[key]<=old_data[key] for key in ('u','u_r','u_rr')),'same completed displacement and no-worse absolute allowances')
    for key in ('e0','e0_r','e0_rr'):
        need(error[key]==old_error[key],'exactly frozen displacement errors')
    for key in ('v','v_r'): data[key]=native_data[key]
    need(all(error[key]<=old_error[key] for key in ('v0','v0_r')),'fresh complete velocity error gates')
    result=m.reconstruct(k,data,error)
    gates={s:2*result['second_coordinate_error_H_caps'][s]<=Q(published['ledger']['second_coordinate_error_H_caps'][s]) for s in ('xx','xr','rr')}
    gates['C1_no_increase']=result['complete_C1_error_X_cap']<=Q(published['ledger']['complete_C1_error_X_cap'])
    need(all(gates.values()),'separate new candidate full response gates')
    faults=('reuse_old_velocity_absolute','reuse_old_velocity_errors','omit_f00','omit_r00','omit_full_Robin_source','omit_original_Fx','omit_derivative_commutator','lose_r1_initial_cancellation','omit_feedback_inverse')
    witnesses={}
    for fault in faults:
        d,e=dict(data),dict(error)
        if fault=='reuse_old_velocity_absolute':
            for key in ('v','v_r'):d[key]=old_data[key]
            wrong=m.reconstruct(k,d,e)
        elif fault=='reuse_old_velocity_errors':
            for key in ('v0','v0_r'):e[key]=old_error[key]
            wrong=m.reconstruct(k,d,e)
        else: wrong=m.reconstruct(k,d,e,fault)
        delta=wrong['second_coordinate_error_H_caps']['rr']-result['second_coordinate_error_H_caps']['rr']
        need(delta!=0,'effective new candidate corruption: '+fault)
        witnesses[fault]={'wrong_minus_correct_rr':delta,'rejected':True}
        if args.fault==fault: m.exact_equal(wrong,result,'CORRUPTION_REJECTED:'+fault)
    if args.verify_report:
        check=read(args.verify_report)
        m.exact_equal(check['ledger'],result,'verified Cauchy report exact nodes')
        m.exact_equal(check['constants'],k,'verified Cauchy report unchanged constants')
        need(m.caps(check['candidate_initial_absolute_H_caps'],m.CAPS)==data,'verified Cauchy exact absolute caps')
        need(check['candidate']=='W1_exact_Cauchy_packet' and check['same_function_as_W1_exact_packet'] is False,'verified distinct candidate role')
        need(m.caps(check['initial_error_H_caps'],m.ERRORS)==error,'verified fresh Cauchy errors')
    comparisons={node:{'before':value,'after':m.flatten(result)[node],
        'old_over_new':value/m.flatten(result)[node] if m.flatten(result)[node] else None} for node,value in m.flatten(baseline).items()}
    for path,expected in before.items(): need(sha(path)==expected,'selected input unchanged on post-readback: '+path)
    out={'schema':'HDBLAST_SEPARATE_EXACT_CAUCHY_RESPONSE_V1','status':'PASS_SEPARATE_EXACT_CAUCHY_CONDITIONAL_RESPONSE',
        'candidate':'W1_exact_Cauchy_packet','same_function_as_W1_exact_packet':False,
        'older_candidate_functions_and_certificates_preserved':True,
        'fresh_analytic_velocity_absolute_and_paired_error_caps_used':True,
        'old_P1_velocity_caps_not_used_for_new_candidate':True,
        'actual_target_original_F_Fx_G_M_operators_domain_H_unchanged':True,
        'baseline_exact_node_count':len(nodes),'constants':k,'domain':published['domain'],'labels':labels,'Gram':published['Gram'],
        'candidate_initial_absolute_H_caps':data,'initial_error_H_caps':error,
        'field_roles':{'u/u_r/u_rr':'same displacement function, refined absolute caps','e0/e0_r/e0_rr':'exact frozen published PR109 displacement errors','v/v_r':'fresh native analytic nominal velocity absolute caps','v0/v0_r':'fresh native paired actual-minus-analytic-nominal velocity errors'},
        'initial_cancellation':'r1(0)=0 from common first Cauchy jet; f00/r00 initial acceleration charges retained',
        'regularity':'Completed displacement C2; new completed velocity C1 globally with exact component affine outer extension and H-isometric parity; flat through1. Conditional physical-side C2 response under retained coefficient/source/corner premises.',
        'full_response_gates_against_PR109':gates,'ledger':result,'node_comparisons':comparisons,
        'formula_and_stale_velocity_corruption_witnesses':witnesses,
        'native_and_independent_initial_local_and_regularity_gates_verified':True,
        'no_new_evolved_point_table_natural_constraint_residual_or_observable':True,
        'no_global_true_error_factor_claim':True,
        'selected_input_sha256_before':portable,'selected_input_sha256_after':portable,'selected_inputs_unchanged':True,
        'scope':'Distinct new exact-Cauchy candidate. Global comparisons are conditional coordinate upper certificates; the separately fresh local initial velocity witness supplies the local true-error statement. No inherited old point table, nonlinear shell-contact result, invariant observable, or cosmological claim.'}
    if args.verify_report:
        need(read(args.verify_report)==encode(out),'verified Cauchy complete canonical report identity')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(out),indent=2,sort_keys=True)+'\n')

if __name__=='__main__': main()
