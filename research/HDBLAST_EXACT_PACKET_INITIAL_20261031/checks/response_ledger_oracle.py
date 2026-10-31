#!/usr/bin/env python3
"""Exact original-source finite-response majorants; no scientific imports.

The interval producer is an authenticated input. This is an independent
arithmetic reconstruction of the response theorem, not of its interval jets.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

PINS = {
 'transfer': ('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TRANSFER.json', '4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3'),
 'time': ('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json', '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d'),
 'response': ('HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json', 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84'),
 'absolute': ('HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json', '44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f'),
 'original_c2': ('HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/C2.json', '64967ab9afda32150153c934522bfb3d9784a6a96515fec6ef34343057d71a58'),
 'fx': ('HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json', '63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60'),
 'combined': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/COMBINED.json', '3e4fb0174df7e88c566346a1b68aedbe4e731fbd062de11eafa28a2cae0a4989'),
 'old_c2': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json', '85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923'),
 'normalization': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/NORMALIZATION.json', 'e63711829cad2e3dea6ac7e133f754d0e25d04fc40c460c681c7ed1ef3a67a5a'),
}

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def encode(x):
    if isinstance(x, Q): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)): return [encode(v) for v in x]
    return x

def read_inputs(root):
    data, paths, portable = {}, {}, {}
    for name, (relative, digest) in PINS.items():
        path = Path(root)/'research'/relative
        need(path.is_file() and not path.is_symlink(), 'regular pinned input: '+name)
        need(sha(path) == digest, 'input hash: '+name)
        data[name] = json.loads(path.read_text())
        paths[str(path.resolve())] = digest
        portable['research/'+relative] = digest
    return data, paths, portable

def constants(inputs):
    r, t, dt, p, fx = (inputs[k] for k in ('response','transfer','time','absolute','fx'))
    op = {k: Q(v) for k,v in r['actual_coefficient_H_norm_caps'].items()}
    fr = {k: Q(v) for k,v in t['frozen_H_operator_caps'].items()}
    de = {k: Q(v) for k,v in r['actual_coefficient_error_H_norm_caps'].items()}
    der = {k: Q(v) for k,v in dt['global']['H_operator_caps'].items()}
    need(all(op[k] == fr[k]+de[k] for k in op), 'actual=frozen+discrepancy')
    need(der['B_x'] == der['M_x'] == 0, 'actual B_x=M_x=0')
    need(fx['scope']['actual_original_full23_F_x_bound'] is True and
         fx['source_coordinate'] == 'x=t/ell' and fx['status'].startswith('PASS'),
         'original full23 source derivative scope')
    out = {'actual':op,'frozen':fr,'delta':de,'derivative':der,'T':Q(3,2),
           'q':Q(r['retained_coupled_X_operator_cap']),
           'F':Q(p['source_caps']['F_H']), 'Fx':Q(fx['genuine_scaled_F_x_H_cap']),
           'G':Q(r['actual_affine_boundary_G_H_cap']), 'Gx':Q(p['source_caps']['G_x_H']),
           'KF':Q(r['affine_bulk_response_KF_X_cap'])}
    out['c'] = op['A']+op['B']+Q(3,4)*op['D']
    out['c0'] = fr['A']+fr['B']+Q(3,4)*fr['D']
    out['dc'] = de['A']+de['B']+Q(3,4)*de['D']
    out['comm'] = der['A_x']+der['B_x']+Q(3,4)*der['D_x']
    out['bdry'] = op['M']/4
    out['Pinv'] = 1/(1-out['bdry'])
    out['betaR'] = (out['q']-out['bdry'])/(1-out['bdry'])
    need(out['q'] == Q(3409,8000000) and 0 < out['bdry'] < 1 and
         out['T']*out['c']+out['bdry'] <= out['q'] < 1,
         'same verified causal contractions')
    need(out['F'] == 17387513156573 and out['Fx'] == 905405545954284 and
         out['G'] == 487122 and out['Gx'] == 58807538,
         'genuine unfiltered source caps')
    need(out['KF'] == Q(52162547172519,2)+out['G'], 'source-response units')
    need(out['Pinv'] == Q(t['P_boundary_inverse_X_cap']) and
         out['comm'] == Q(dt['coefficient_commutator_X_to_H_upper']),
         'same actual Robin inverse and full commutator')
    return out

def allcaps(values, keys, name):
    need(set(values) == set(keys), 'complete cap schema: '+name)
    result = {k:Q(values[k]) for k in keys}
    need(min(result.values()) >= 0, 'nonnegative caps: '+name)
    return result

def differentiated(k, data, error, Bs, D0X, CX, E, fault=None):
    a,b,d,m = (k['actual'][s] for s in ('A','B','D','M'))
    a0,b0,d0 = (k['frozen'][s] for s in ('A','B','D'))
    da,db,dd = (k['delta'][s] for s in ('A','B','D'))
    U,U1,U2,V,V1 = (data[s] for s in ('u','u_r','u_rr','v','v_r'))
    eu,eur,eurr,ev,evr = (error[s] for s in ('e0','e0_r','e0_rr','v0','v0_r'))
    T,P,comm = k['T'],k['Pinv'],k['comm']
    if fault == 'absolute_error_jet_interchange': U2 = eurr
    Qfree = U2+V1
    f00 = a0*V+b0*U1+d0*U
    f0x = (a0+b0)*Qfree+d0*Bs
    D0xX = P*((0 if fault=='omit_frozen_initial_acceleration' else f00)+T*f0x+
               (0 if fault=='omit_Robin_free_wave_boundary_source' else m*Bs)+k['Gx'])
    r00 = da*V+db*U1+dd*U+k['F']
    Fx = k['Fx']
    if fault == 'omit_original_Fx': Fx = Q(0)
    if fault == 'filtered_Fx_substitution': Fx = Q(66768050521240320)
    if fault == 'lose_W0_minus_Bs_cancellation':
        # A valid but unnecessarily broad absolute derivative route differs.
        r0x = (a+b)*(Qfree+D0xX)+(d+comm)*(Bs+D0X)+(a0+b0)*Qfree+d0*Bs+Fx
    else:
        r0x = (a+b)*D0xX+(d+comm)*D0X+(da+db)*Qfree+(dd+comm)*Bs+Fx
    CxX = P*((0 if fault=='omit_bulk_initial_acceleration' else r00)+T*r0x)
    F1 = (a+b)*CxX+(d+comm)*CX
    if fault == 'omit_complete_commutator':
        F1 = (a+b)*CxX+d*CX
        comm = Q(0)
    w0 = eurr+a*ev+b*eur+d*eu
    if fault == 'filtered_initial_error_forcing': w0 += k['F']
    J = max((ev+T*w0)/3,(0 if fault=='omit_velocity_gradient_error' else evr)+w0)
    denominator = 1 if fault == 'omit_full_feedback_inverse' else 1-k['q']
    Z = (J+T*(F1+comm*E))/denominator
    F0 = k['c']*CX
    return {'free_wave_second_jet_H_cap':Qfree,
            'frozen_initial_forcing_H_cap':f00,'frozen_forcing_x_H_cap':f0x,
            'W0_minus_Bs_time_derivative_X_cap':D0xX,
            'r0_initial_H_cap':r00,'r0_time_derivative_H_cap':r0x,
            'correction_time_derivative_X_cap':CxX,
            'reduced_defect_H_cap':F0,'reduced_defect_time_derivative_H_cap':F1,
            'initial_error_acceleration_H_cap':w0,'differentiated_initial_free_wave_X_cap':J,
            'second_coordinate_error_H_caps':{'xx':Z,'xr':Z,'rr':Z+k['c']*E+F0}}

def replay_retained(inputs, k):
    t,p,old = (inputs[s] for s in ('transfer','absolute','old_c2'))
    dr = t['direct_reproof']
    data = allcaps(p['initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'old absolute')
    error = allcaps(old['initial_error_H_caps'],('e0','e0_r','e0_rr','v0','v0_r'),'old error')
    Bs = Q(dr['new_smooth_Bs_X_cap'])
    D0X = Q(dr['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
    CX = Q(dr['new_actual_bulk_correction_X_cap'])
    E = Q(t['best_reproved_global_C1_error_X_cap'])
    out = differentiated(k,data,error,Bs,D0X,CX,E)
    need(out['second_coordinate_error_H_caps'] ==
         {s:Q(v) for s,v in old['original_target_second_coordinate_error_H_caps'].items()},
         'independent replay of every exact retained C2 upper')
    need(out['reduced_defect_time_derivative_H_cap'] ==
         Q(old['original_reduced_defect_time_derivative_H_cap']),
         'retained exact differentiated-defect replay')
    need(out['reduced_defect_H_cap'] == Q(old['complete_original_reduced_defect_H_cap']),
         'retained exact reduced-defect replay')
    need(out['initial_error_acceleration_H_cap'] == Q(old['initial_differentiated_acceleration_H_cap']) and
         out['differentiated_initial_free_wave_X_cap'] == Q(old['differentiated_initial_free_wave_X_cap']),
         'retained exact differentiated initial problem')
    return out

def local_gate(campaign, old):
    labels, gram = old['labels'], list(map(Q,old['Gram']))
    local = campaign['local']
    need(local['point'] == ['0','65/64'] and len(local['rows']) == 23,
         'same complete local point')
    A2lo=A2hi=E2=Q(0)
    for label,w,row in zip(labels,gram,local['rows']):
        need(row['label'] == label and Q(row['Gram']) == w, 'local all23 Gram role')
        lo,hi=map(Q,row['actual_u_rr']); el,eh=map(Q,row['paired_actual_minus_candidate_u_rr'])
        need(lo<=hi and el<=eh, 'ordered local component intervals')
        amin=Q(0) if lo<=0<=hi else min(abs(lo),abs(hi))
        A2lo+=w*amin*amin; A2hi+=w*max(abs(lo),abs(hi))**2
        E2+=w*max(abs(el),abs(eh))**2
    low,up=map(Q,local['actual_H_norm_interval']); err=Q(local['new_error_H_upper'])
    oldlower=Q(local['old_error_H_strict_lower'])
    need(0<low and low*low<=A2lo and up*up>=A2hi and err*err>=E2,
         'independent local Gram-square aggregation')
    need(oldlower == Q(839330547758611430698479,1208925819614629174706176),
         'same retained old true-error lower witness')
    need(10*err <= low and 2*err < oldlower, 'local true-error acceptance inequalities')
    return {'target_H_interval':[low,up],'new_error_H_upper':err,
            'target_H_squared_lower':A2lo,'target_H_squared_upper':A2hi,
            'new_error_H_squared_upper':E2,'old_error_H_strict_lower':oldlower,
            'new_relative_error_H_upper':err/low,'old_relative_error_H_strict_lower':oldlower/up,
            'ten_percent_gate':True,'old_lower_half_gate':True}

def fresh_majorants(k,data,error):
    T=k['T']
    Bs=max((data['u']+T*data['v'])/3,data['u_r']+data['v'])
    I=max((error['e0']+T*error['v0'])/3,error['e0_r']+error['v0'])
    D0X=k['Pinv']*((T*k['c0']+k['bdry'])*Bs+k['G'])
    Rcoeff=k['dc']*Bs; Rfeedback=k['c']*D0X
    R=T*(Rcoeff+Rfeedback)+k['KF']
    CX=k['Pinv']*R
    E=(I+k['betaR']*R)/(1-k['q'])
    result=differentiated(k,data,error,Bs,D0X,CX,E)
    result.update({'free_wave_X_cap':Bs,'initial_free_wave_X_error_cap':I,
       'W0_minus_Bs_X_cap':D0X,'W0_X_cap':Bs+D0X,
       'coefficient_mismatch_residual_H_cap':Rcoeff,'feedback_residual_H_cap':Rfeedback,
       'bulk_residual_response_X_cap':R,'bulk_correction_X_cap':CX,
       'W1_X_cap':Bs+D0X+CX,'complete_C1_error_X_cap':E})
    # Direct prescribed-forcing absolute candidate bounds. Distinct from errors.
    a,b,d=(k['actual'][s] for s in ('A','B','D'))
    f1init=a*data['v']+b*data['u_r']+d*data['u']+k['F']
    f1x=(a+b)*(result['free_wave_second_jet_H_cap']+result['W0_minus_Bs_time_derivative_X_cap'])+(d+k['comm'])*result['W0_X_cap']+k['Fx']
    acc=data['u_rr']+f1init
    Jabs=max((data['v']+T*acc)/3,data['v_r']+acc)
    W1x=k['Pinv']*(Jabs+T*f1x+k['Gx'])
    result['absolute_W1_H_caps']={'value':Q(3,4)*result['W1_X_cap'],
      'x':result['W1_X_cap'],'r':result['W1_X_cap'], 'xx':W1x,'xr':W1x,
      'rr':W1x+k['c']*result['W0_X_cap']+k['F']}
    result['absolute_W1_prescribed_forcing_ledger']={'f1_initial_H_cap':f1init,
      'f1_time_derivative_H_cap':f1x,'initial_acceleration_H_cap':acc,
      'differentiated_free_wave_X_cap':Jabs,'W1_time_derivative_X_cap':W1x}
    return result

FAULTS=('absolute_error_jet_interchange','omit_frozen_initial_acceleration',
 'omit_Robin_free_wave_boundary_source','lose_W0_minus_Bs_cancellation',
 'omit_bulk_initial_acceleration','omit_original_Fx','filtered_Fx_substitution',
 'omit_complete_commutator','omit_velocity_gradient_error',
 'filtered_initial_error_forcing','omit_full_feedback_inverse')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--candidate-report',type=Path,required=True)
    ap.add_argument('--candidate-sha256',required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--fault',choices=FAULTS)
    args=ap.parse_args()
    need(not args.output.exists(), 'fresh absent output required')
    me=Path(__file__).resolve(); before={str(me):sha(me)}
    inputs,pins,portable=read_inputs(args.repo_root);before.update(pins)
    portable['response_ledger_oracle.py']=before[str(me)]
    need(sha(args.candidate_report)==args.candidate_sha256,'candidate report hash')
    before[str(args.candidate_report.resolve())]=args.candidate_sha256
    portable['candidate_initial_report']=args.candidate_sha256
    report=json.loads(args.candidate_report.read_text())
    need(report['candidate']=='W1_exact_packet' and report['target_and_original_F_Fx_G_M_unchanged'] is True,
         'frozen new candidate and unchanged actual target/source')
    successful=[c for c in report['campaigns'] if c['status']=='INITIAL_AND_LOCAL_GATES_PASS']
    need(len(successful)==1,'one successful initial campaign before response')
    campaign=successful[0]
    old=inputs['old_c2'];k=constants(inputs);replay=replay_retained(inputs,k)
    need(campaign['labels']==old['labels'] and campaign['Gram']==old['Gram'] and
         len(old['labels'])==23 and min(map(Q,old['Gram']))>0,'complete same modal space')
    data=allcaps(campaign['candidate_initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'new absolute')
    error=allcaps(campaign['initial_error_H_caps'],('e0','e0_r','e0_rr','v0','v0_r'),'new error')
    olderr={s:Q(v) for s,v in old['initial_error_H_caps'].items()}
    initial_gate={s:error[s]<=olderr[s] for s in ('e0','e0_r')}
    initial_gate['e0_rr_half']=2*error['e0_rr']<=Q(inputs['combined']['refined_uniform_initial_u_rr_H_upper'])
    initial_gate.update({s:error[s]==olderr[s] for s in ('v0','v0_r')})
    need(all(initial_gate.values()),'exact complete initial acceptance inequalities')
    need(campaign['unresolved_cells']==[] and campaign['initial_gate_pass'] is True,
         'initial campaign no unresolved cells')
    local=local_gate(campaign,old)
    result=fresh_majorants(k,data,error)
    # Fault executions must fail before writing a success artifact. Each has a
    # numerical witness against the complete original-source expression.
    witnesses={}
    for fault in FAULTS:
        wrong=differentiated(k,data,error,result['free_wave_X_cap'],
          result['W0_minus_Bs_X_cap'],result['bulk_correction_X_cap'],
          result['complete_C1_error_X_cap'],fault)
        delta=wrong['second_coordinate_error_H_caps']['rr']-result['second_coordinate_error_H_caps']['rr']
        need(delta!=0,'effective nonzero formula corruption: '+fault)
        witnesses[fault]={'wrong_minus_correct_rr_majorant':delta,'rejected_by_exact_identity':True}
        if args.fault==fault:
            need(delta==0,'CORRUPTION_REJECTED:'+fault+':nonzero_exact_residual='+str(delta))
    oldcaps={s:Q(v) for s,v in old['original_target_second_coordinate_error_H_caps'].items()}
    full_gate={s:2*result['second_coordinate_error_H_caps'][s]<=oldcaps[s] for s in oldcaps}
    full_gate['C1_no_increase']=result['complete_C1_error_X_cap']<=Q(old['complete_original_C1_error_X_cap'])
    scale=Q(3249*10**30)
    need(scale==Q('3.249e33'),'same global actual-target lower scale')
    output={'status':'PASS_CONDITIONAL_EXACT_PACKET_RESPONSE_CERTIFICATE' if all(full_gate.values()) else 'FULL_RESPONSE_CERTIFICATE_GATE_UNRESOLVED',
      'capsule':'HDBLAST_EXACT_PACKET_INITIAL_20261031','candidate':'W1_exact_packet',
      'target':old['target'],'domain':old['domain'],'labels':old['labels'],'Gram':old['Gram'],
      'actual_target_original_F_Fx_G_M_unchanged':True,
      'finite_definition':{'W0':'P[Bs_new+K L0 Bs_new+LG]','W1':'P[Bs_new+K(La W0_new+F)+LG]'},
      'regularity':'completed displacement C2; unchanged smoothed P1 velocity C4; physical-side C2 response with one-sided shell traces under retained source/coefficient and compatible corner premises',
      'initial_cancellation':'r1(0)=0 from C=W1-W0 zero initial value, C_x and C_r; f0(0) and r0(0) acceleration charges retained',
      'boundary_error':{'g':'0','g_x':'0','M_x':'0'},
      'new_point_or_quadrature_values_produced':False,
      'input_sha256':portable,'candidate_initial_absolute_H_caps':data,'initial_error_H_caps':error,
      'constants':k,'initial_certificate_gates':initial_gate,'local_gate':local,
      'full_response_certificate_gates':full_gate,'ledger':result,
      'old_exact_second_error_H_caps':oldcaps,
      'old_exact_C1_error_X_cap':Q(old['complete_original_C1_error_X_cap']),
      'global_normalization':{'actual_target_strict_lower_scale':scale,
        'initial_rr_upper_ratio':error['e0_rr']/scale,
        'evolved_rr_upper_ratio':result['second_coordinate_error_H_caps']['rr']/scale,
        'evolved_max_second_seminorm_upper_ratio':max(result['second_coordinate_error_H_caps'].values())/scale},
      'retained_exact_C2_replayed':replay,'formula_corruption_witnesses':witnesses,
      'scope':'Conditional absolute coordinate response error bounds from certified new initial inputs. Local true-error improvement separately proved by initial local gate; global factor-two compares upper certificates only. No new point values, natural/constraint residual, invariant observable or nonlinear claim.'}
    need(all(sha(path)==digest for path,digest in before.items()),'every read input and executed helper unchanged on readback')
    output['authentication_before_after_equal']=True
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(output),indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    main()
