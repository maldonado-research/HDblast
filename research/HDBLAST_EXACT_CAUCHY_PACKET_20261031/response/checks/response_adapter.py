#!/usr/bin/env python3
"""Portable exact rational response certificate adapter; no scientific imports.

Reconstruct all published baseline nodes before accepting a refinement of the
same candidate's absolute displacement allowances. Target-error inputs never
come from a refined report. Explicit exceptions survive python -O.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

CAPS = ('u', 'u_r', 'u_rr', 'v', 'v_r')
ERRORS = ('e0', 'e0_r', 'e0_rr', 'v0', 'v0_r')
BASE = 'research/HDBLAST_EXACT_PACKET_INITIAL_20261031/'
SOURCES = {
    'initial': BASE+'reference/INITIAL.json',
    'baseline': BASE+'reference/RESPONSE.json',
    'transfer': 'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TRANSFER.json',
    'time': 'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json',
    'response': 'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json',
    'absolute': 'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json',
    'fx': 'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json',
}

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def digest(path):
    p = Path(path)
    need(p.is_file() and not p.is_symlink(), 'regular input required: '+str(p))
    need(all(not parent.is_symlink() for parent in p.absolute().parents),
         'symlink-free input ancestry: '+str(p))
    return hashlib.sha256(p.read_bytes()).hexdigest()

def encoded(obj):
    if isinstance(obj, Q):
        return str(obj)
    if isinstance(obj, dict):
        return {k: encoded(v) for k,v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [encoded(v) for v in obj]
    return obj

def canonical(obj):
    return json.dumps(encoded(obj), indent=2, sort_keys=True)+'\n'

def caps(obj, names):
    need(set(obj) == set(names), 'exact cap schema')
    out = {k: Q(obj[k]) for k in names}
    need(all(v >= 0 for v in out.values()), 'nonnegative cap values')
    return out

def flatten(obj, prefix=''):
    out = {}
    for k,v in obj.items():
        key = prefix+'/'+k if prefix else k
        if isinstance(v, dict):
            out.update(flatten(v, key))
        else:
            out[key] = Q(v)
    return out

def exact_equal(got, expected, label):
    a,b = flatten(got), flatten(expected)
    need(set(a) == set(b), label+': complete node schema')
    for node in a:
        need(a[node] == b[node], label+': exact node '+node)
    return sorted(a)

def constants(reports):
    r = reports['response']; t = reports['transfer']; dt = reports['time']
    a = {k:Q(v) for k,v in r['actual_coefficient_H_norm_caps'].items()}
    f = {k:Q(v) for k,v in t['frozen_H_operator_caps'].items()}
    d = {k:Q(v) for k,v in r['actual_coefficient_error_H_norm_caps'].items()}
    dx = {k:Q(v) for k,v in dt['global']['H_operator_caps'].items()}
    need(all(a[s] == f[s]+d[s] for s in a), 'actual=frozen+delta')
    need(dx['B_x'] == dx['M_x'] == 0, 'B_x=M_x=0')
    k = dict(actual=a, frozen=f, delta=d, derivative=dx, T=Q(3,2),
        q=Q(r['retained_coupled_X_operator_cap']),
        F=Q(reports['absolute']['source_caps']['F_H']),
        Fx=Q(reports['fx']['genuine_scaled_F_x_H_cap']),
        G=Q(r['actual_affine_boundary_G_H_cap']),
        Gx=Q(reports['absolute']['source_caps']['G_x_H']),
        KF=Q(r['affine_bulk_response_KF_X_cap']))
    for name,op in [('c',a), ('c0',f), ('dc',d)]:
        k[name] = op['A']+op['B']+Q(3,4)*op['D']
    k['comm'] = dx['A_x']+dx['B_x']+Q(3,4)*dx['D_x']
    k['bdry'] = a['M']/4
    k['Pinv'] = Q(1)/(1-k['bdry'])
    k['betaR'] = (k['q']-k['bdry'])/(1-k['bdry'])
    need(k['T']*k['c']+k['bdry'] <= k['q'] < 1, 'causal contraction')
    need((k['F'],k['Fx'],k['G'],k['Gx']) == (17387513156573,905405545954284,487122,58807538), 'original source values')
    need(k['Pinv'] == Q(t['P_boundary_inverse_X_cap']), 'full Robin inverse')
    need(k['comm'] == Q(dt['coefficient_commutator_X_to_H_upper']), 'full derivative commutator')
    need(reports['fx']['scope']['actual_original_full23_F_x_bound'] is True, 'original full23 Fx')
    return k

def reconstruct(k, absolute, error, fault=None):
    """Nonnegative grouped equations, independent of the published oracle."""
    u,ur,urr,v,vr = (absolute[s] for s in CAPS)
    eu,eur,eurr,ev,evr = (error[s] for s in ERRORS)
    A,B,D,M = (k['actual'][s] for s in ('A','B','D','M'))
    A0,B0,D0 = (k['frozen'][s] for s in ('A','B','D'))
    dA,dB,dD = (k['delta'][s] for s in ('A','B','D'))
    T,P,c,c0,dc,comm = (k[s] for s in ('T','Pinv','c','c0','dc','comm'))
    if fault == 'absolute_error_curvature_interchange':
        urr = eurr
    Bs = max(u/3+T*v/3, ur+v)
    I = max(eu/3+T*ev/3, eur+ev)
    correction0 = P*(T*c0*Bs+k['bdry']*Bs+k['G'])
    mismatch, feedback = dc*Bs, c*correction0
    residual = T*mismatch+T*feedback+k['KF']
    correction = P*residual
    E = (I+k['betaR']*residual)/(1-k['q'])
    free2 = urr+vr
    f00 = A0*v+B0*ur+D0*u
    f0x = A0*free2+B0*free2+D0*Bs
    Robin = M*Bs+k['Gx']
    if fault == 'omit_full_Robin_source':
        Robin = Q(0)
    acceleration0 = Q(0) if fault == 'omit_f00' else f00
    correction0x = P*(acceleration0+T*f0x+Robin)
    r00 = dA*v+dB*ur+dD*u+k['F']
    Fx = Q(0) if fault == 'omit_original_Fx' else k['Fx']
    r0x = ((A+B)*correction0x+(D+comm)*correction0+
           (dA+dB)*free2+(dD+comm)*Bs+Fx)
    acceleration1 = Q(0) if fault == 'omit_r00' else r00
    correctionx = P*(acceleration1+T*r0x)
    defect0 = c*correction
    defectx = (A+B)*correctionx+(D+comm)*correction
    err_acc = eurr+A*ev+B*eur+D*eu
    initx = max(ev/3+T*err_acc/3, evr+err_acc)
    denominator = Q(1) if fault == 'omit_feedback_inverse' else 1-k['q']
    if fault == 'omit_derivative_commutator':
        defectx = (A+B)*correctionx+D*correction
        comm = Q(0)
    # r1(0)=0 exactly: the correction has zero first Cauchy jet.
    invented_r10 = r00 if fault == 'lose_r1_initial_cancellation' else Q(0)
    Z = (initx+invented_r10+T*defectx+T*comm*E)/denominator
    W0 = Bs+correction0
    W1 = W0+correction
    f1init = A*v+B*ur+D*u+k['F']
    f1x = (A+B)*(free2+correction0x)+(D+k['comm'])*W0+k['Fx']
    acc = urr+f1init
    Jabs = max(v/3+T*acc/3, vr+acc)
    W1x = P*(Jabs+T*f1x+k['Gx'])
    return {
        'free_wave_X_cap':Bs, 'initial_free_wave_X_error_cap':I,
        'W0_minus_Bs_X_cap':correction0, 'W0_X_cap':W0,
        'coefficient_mismatch_residual_H_cap':mismatch,
        'feedback_residual_H_cap':feedback,
        'bulk_residual_response_X_cap':residual,
        'bulk_correction_X_cap':correction, 'W1_X_cap':W1,
        'complete_C1_error_X_cap':E,
        'free_wave_second_jet_H_cap':free2,
        'frozen_initial_forcing_H_cap':f00,
        'frozen_forcing_x_H_cap':f0x,
        'W0_minus_Bs_time_derivative_X_cap':correction0x,
        'r0_initial_H_cap':r00, 'r0_time_derivative_H_cap':r0x,
        'correction_time_derivative_X_cap':correctionx,
        'reduced_defect_H_cap':defect0,
        'reduced_defect_time_derivative_H_cap':defectx,
        'initial_error_acceleration_H_cap':err_acc,
        'differentiated_initial_free_wave_X_cap':initx,
        'second_coordinate_error_H_caps':{'xx':Z,'xr':Z,'rr':Z+c*E+defect0},
        'absolute_W1_H_caps':{'value':Q(3,4)*W1,'x':W1,'r':W1,
            'xx':W1x,'xr':W1x,'rr':W1x+c*W0+k['F']},
        'absolute_W1_prescribed_forcing_ledger':{
            'f1_initial_H_cap':f1init,'f1_time_derivative_H_cap':f1x,
            'initial_acceleration_H_cap':acc,
            'differentiated_free_wave_X_cap':Jabs,
            'W1_time_derivative_X_cap':W1x},
    }

def sensitivity(k, data, error, baseline):
    # Zero is a hypothetical ledger input; it is not a certified candidate cap.
    zero = dict(data, u_rr=Q(0)); one = dict(data, u_rr=Q(1))
    floor = reconstruct(k,zero,error)
    step = reconstruct(k,one,error)
    out = {}
    for name in ('xx','xr','rr'):
        before = baseline['second_coordinate_error_H_caps'][name]
        low = floor['second_coordinate_error_H_caps'][name]
        slope = step['second_coordinate_error_H_caps'][name]-low
        need(slope>0 and before == low+slope*data['u_rr'], 'exact affine rr sensitivity')
        out[name] = {'current':before,'zero_u_rr_hypothetical_floor':low,
            'slope_per_unit_u_rr':slope,'fraction_due_to_u_rr':(before-low)/before,
            'required_u_rr_for_half_current':(before/2-low)/slope,
            'half_current_possible_by_u_rr_alone':low<before/2}
    return {'status':'EXACT_LEDGER_SENSITIVITY_ONLY','not_a_candidate_or_true_error_lower_bound':True,
        'fixed_absolute_u_u_r_velocity_and_all_initial_errors':True,'second_errors':out}

FAULTS = ('absolute_error_curvature_interchange','omit_full_Robin_source',
          'omit_f00','omit_r00','omit_original_Fx','omit_feedback_inverse',
          'omit_derivative_commutator','lose_r1_initial_cancellation')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--freeze',type=Path,default=Path(__file__).with_name('FREEZE.json'))
    ap.add_argument('--freeze-sha256',required=True)
    ap.add_argument('--caps',type=Path)
    ap.add_argument('--caps-sha256')
    ap.add_argument('--proof-input',action='append',default=[],help='path=sha256')
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--verify-report',type=Path)
    ap.add_argument('--fault',choices=FAULTS)
    args=ap.parse_args()
    need(not args.output.exists(),'fresh output required')
    need(digest(args.freeze)==args.freeze_sha256,'frozen recipe hash')
    freeze=json.loads(args.freeze.read_text())
    need(digest(__file__)==freeze['adapter_sha256'],'frozen ledger recipe source')
    selected={str(Path(__file__).resolve()):digest(__file__),str(args.freeze.resolve()):args.freeze_sha256}
    portable={'response_adapter.py':digest(__file__),'FREEZE.json':args.freeze_sha256}
    for relative,expected in freeze['source_pins'].items():
        p=args.repo_root/relative
        need(digest(p)==expected,'selected source pin: '+relative)
        selected[str(p.resolve())]=expected; portable[relative]=expected
    reports={k:json.loads((args.repo_root/path).read_text()) for k,path in SOURCES.items()}
    initial=reports['initial']; published=reports['baseline']
    need(initial['candidate']==published['candidate']=='W1_exact_packet','same candidate identity')
    campaigns=[c for c in initial['campaigns'] if c['status']=='INITIAL_AND_LOCAL_GATES_PASS']
    need(len(campaigns)==1,'one frozen initial campaign')
    campaign=campaigns[0]
    data=caps(campaign['candidate_initial_absolute_H_caps'],CAPS)
    error=caps(campaign['initial_error_H_caps'],ERRORS)
    need(encoded(data)==published['candidate_initial_absolute_H_caps'],'baseline absolute inputs')
    need(encoded(error)==published['initial_error_H_caps'],'baseline errors fixed from INITIAL')
    k=constants(reports)
    exact_equal(k,published['constants'],'all frozen constants')
    baseline=reconstruct(k,data,error)
    nodes=exact_equal(baseline,published['ledger'],'every baseline C1/C2 node')
    refined=dict(data)
    cap_metadata=None
    if args.caps:
        need(args.caps_sha256 is not None and digest(args.caps)==args.caps_sha256,'refined caps report pin')
        selected[str(args.caps.resolve())]=args.caps_sha256; portable['refined_caps_report']=args.caps_sha256
        report=json.loads(args.caps.read_text())
        need(report['candidate']=='W1_exact_packet','refined same candidate identity')
        refined=caps(report['candidate_initial_absolute_H_caps'],CAPS)
        if 'initial_error_H_caps' in report:
            need(caps(report['initial_error_H_caps'],ERRORS)==error,'refined report cannot change target-error caps')
        need(all(refined[s]<=data[s] for s in ('u','u_r','u_rr')),'no worse absolute displacement charges')
        need(all(refined[s]==data[s] for s in ('v','v_r')),'unchanged absolute velocity caps')
        need(any(refined[s]<data[s] for s in ('u','u_r','u_rr')),'strict refinement required')
        cap_metadata={'sha256':args.caps_sha256,'proof_authority':'externally authenticated same-candidate absolute-cap report; response adapter does not reprove interval leaves'}
    else:
        need(args.caps_sha256 is None,'no caps hash without report')
    for item in args.proof_input:
        name,sep,expected=item.rpartition('=')
        need(sep and len(expected)==64,'proof-input path=sha256')
        p=Path(name); need(digest(p)==expected,'proof source/table/report pin: '+str(p))
        selected[str(p.resolve())]=expected; portable['proof-input:'+p.name]=expected
    result=reconstruct(k,refined,error)
    if args.verify_report:
        checked=json.loads(args.verify_report.read_text())
        exact_equal(checked['ledger'],result,'verified report every exact node')
        need(caps(checked['frozen_initial_error_H_caps'],ERRORS)==error,'verified report unchanged errors')
        need(caps(checked['candidate_initial_absolute_H_caps'],CAPS)==refined,'verified report exact absolute caps')
        need(checked['same_function_true_error_unchanged'] is True,'verified report same true error')
    if args.fault:
        bad=reconstruct(k,refined,error,args.fault)
        exact_equal(bad,result,'CORRUPTION_REJECTED:'+args.fault)
    witnesses={}
    for fault in FAULTS:
        bad=reconstruct(k,refined,error,fault)
        delta=bad['second_coordinate_error_H_caps']['rr']-result['second_coordinate_error_H_caps']['rr']
        need(delta!=0,'meaningful formula corruption: '+fault)
        witnesses[fault]={'wrong_minus_correct_rr':delta,'rejected':True}
    comparisons={}
    for name,new in flatten(result).items():
        old=flatten(baseline)[name]
        need(new<=old,'no worse node: '+name)
        comparisons[name]={'before':old,'after':new,'decrease':old-new,
            'old_over_new':old/new if new else None,'strictly_smaller':new<old}
    for path,expected in selected.items():
        need(digest(path)==expected,'post-read selected source unchanged: '+path)
    output={'schema':'HDBLAST_SAME_FUNCTION_RESPONSE_REFINEMENT_V1',
        'status':'PASS_SAME_FUNCTION_CONDITIONAL_CERTIFICATE_REFINEMENT' if args.caps else 'PASS_ALL_EXACT_BASELINE_NODES',
        'candidate':'W1_exact_packet','same_function_true_error_unchanged':True,
        'same_candidate_absolute_certificate_refinement_only':True,
        'actual_target_original_F_Fx_G_M_operators_domain_H_unchanged':True,
        'domain':published['domain'],'labels':published['labels'],'Gram':published['Gram'],
        'finite_definition':published['finite_definition'],
        'initial_cancellation':'r1(0)=0; common first Cauchy jet retained; f00/r00 initial acceleration charges retained',
        'boundary_error':published['boundary_error'],
        'physical_side_estimates_conditional_on_retained_source_coefficient_and_corner_premises':True,
        'no_new_waveform_point_table_residual_observable_or_nonlinear_claim':True,
        'frozen_initial_error_H_caps':error,
        'candidate_initial_absolute_H_caps':refined,'baseline_absolute_H_caps':data,
        'constants':k,'baseline_exact_node_count':len(nodes),'baseline_exact_nodes':nodes,
        'ledger':result,'node_comparisons':comparisons,
        'sensitivity':sensitivity(k,data,error,baseline),
        'formula_corruption_witnesses':witnesses,'refined_cap_metadata':cap_metadata,
        'selected_input_sha256_before':portable,'selected_input_sha256_after':portable,
        'selected_inputs_unchanged':True,
        'scope':'Same candidate and same true error. Only conditional coordinate upper certificates improve. No inherited old point/quadrature or natural/constraint residual table.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(canonical(output))

if __name__=='__main__':
    main()
