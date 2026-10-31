#!/usr/bin/env python3
"""Independent exact rational equation ledger; imports only Python stdlib.

Interval expression evaluation has its own independent audit. This checker
authenticates that audit's producer input, checks every positive-Gram leaf
aggregation/coverage and local witness, then builds response majorants as a
sum of named nonnegative terms. No retained scientific program is imported.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


PIN = {
    'transfer': ('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TRANSFER.json', '4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3'),
    'derivatives': ('HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json', '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d'),
    'operators': ('HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json', 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84'),
    'shared_absolute': ('HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/CANDIDATE.json', '44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f'),
    'original_Fx': ('HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json', '63ed8e3a7c0544db93de0d3bcbdba4652eb4d9a6bcd8f434327e785cd9986e60'),
    'old_combined': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/COMBINED.json', '3e4fb0174df7e88c566346a1b68aedbe4e731fbd062de11eafa28a2cae0a4989'),
    'old_C2': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json', '85bf4133d931d93401b339dc508553ae488bf91f43f75bd7e736432cc856a923'),
    'normalization': ('HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/NORMALIZATION.json', 'e63711829cad2e3dea6ac7e133f754d0e25d04fc40c460c681c7ed1ef3a67a5a'),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    return value


def rational_dict(data):
    return {key: F(value) for key, value in data.items()}


def scalar_caps(data, keys, role):
    require(set(data) == set(keys), role + ': exact key set')
    result = rational_dict(data)
    require(min(result.values()) >= 0, role + ': nonnegative caps')
    return result


def inputs(root, record):
    result = {}
    for role, (name, expected) in PIN.items():
        path = root / 'research' / name
        require(path.is_file() and not path.is_symlink(), 'regular retained input: ' + role)
        require(digest(path) == expected, 'retained hash mismatch: ' + role)
        record[path] = expected
        result[role] = json.loads(path.read_text())
    return result


def constants(old):
    a = rational_dict(old['operators']['actual_coefficient_H_norm_caps'])
    z = rational_dict(old['transfer']['frozen_H_operator_caps'])
    delta = rational_dict(old['operators']['actual_coefficient_error_H_norm_caps'])
    dt = rational_dict(old['derivatives']['global']['H_operator_caps'])
    require(all(a[s] == z[s] + delta[s] for s in a), 'operator cap split')
    require(dt['B_x'] == dt['M_x'] == 0, 'constant B and actual Robin M')
    proof = old['original_Fx']
    require(proof['status'].startswith('PASS') and proof['source_coordinate'] == 'x=t/ell'
            and proof['scope']['actual_original_full23_F_x_bound'] is True,
            'original unfiltered scaled source derivative role')
    T = F(3, 2)
    c = sum((a['A'], a['B'], F(3, 4) * a['D']))
    cz = sum((z['A'], z['B'], F(3, 4) * z['D']))
    cd = sum((delta['A'], delta['B'], F(3, 4) * delta['D']))
    comm = sum((dt['A_x'], dt['B_x'], F(3, 4) * dt['D_x']))
    h = a['M'] / 4
    q = F(old['operators']['retained_coupled_X_operator_cap'])
    require(q == F(3409, 8000000) and T*c+h <= q < 1 and 0 < h < 1,
            'retained coupled contraction and support factors')
    inv = 1 / (1-h)
    require(inv == F(old['transfer']['P_boundary_inverse_X_cap']), 'same actual Robin inverse')
    require(comm == F(old['derivatives']['coefficient_commutator_X_to_H_upper']),
            'complete coefficient commutator')
    source = old['shared_absolute']['source_caps']
    F0, Fx = F(source['F_H']), F(proof['genuine_scaled_F_x_H_cap'])
    G, Gx = F(old['operators']['actual_affine_boundary_G_H_cap']), F(source['G_x_H'])
    KF = F(old['operators']['affine_bulk_response_KF_X_cap'])
    require((F0, Fx, G, Gx) == (17387513156573, 905405545954284, 487122, 58807538),
            'fixed original F Fx G Gx')
    require(KF == F(52162547172519, 2) + G, 'source response cap has X units')
    return dict(actual=a, frozen=z, difference=delta, derivatives=dt, T=T,
                c=c, c0=cz, dc=cd, comm=comm, boundary=h, q=q, P=inv,
                beta=(q-h)*inv, F=F0, Fx=Fx, G=G, Gx=Gx, KF=KF)


def response_equations(k, u, e, inherited=None):
    """Named positive terms expose mandatory charges and independent grouping."""
    A,B,D,M = (k['actual'][s] for s in ('A','B','D','M'))
    A0,B0,D0 = (k['frozen'][s] for s in ('A','B','D'))
    DA,DB,DD = (k['difference'][s] for s in ('A','B','D'))
    T,P,Qinv = k['T'],k['P'],1/(1-k['q'])
    terms, v = {}, {}
    def put(name, **parts):
        require(min(parts.values()) >= 0, 'nonnegative equation terms: '+name)
        terms[name] = parts
        v[name] = sum(parts.values())
        return v[name]
    if inherited is None:
        v['free_wave_X_cap'] = max((u['u']+T*u['v'])/3,u['u_r']+u['v'])
        v['initial_free_wave_X_error_cap'] = max((e['e0']+T*e['v0'])/3,e['e0_r']+e['v0'])
        Bs = v['free_wave_X_cap']
        dfield = put('W0_minus_Bs_X_cap', frozen_bulk=P*T*k['c0']*Bs,
                     Robin_free_wave=P*k['boundary']*Bs, affine_shell=P*k['G'])
        put('coefficient_mismatch_residual_H_cap', coefficient_difference=k['dc']*Bs)
        put('feedback_residual_H_cap', feedback=k['c']*dfield)
        residual = put('bulk_residual_response_X_cap', coefficient=T*v['coefficient_mismatch_residual_H_cap'],
                       feedback=T*v['feedback_residual_H_cap'], affine_bulk=k['KF'])
        correction = put('bulk_correction_X_cap', Robin_inverted_residual=P*residual)
        E = put('complete_C1_error_X_cap', initial=Qinv*v['initial_free_wave_X_error_cap'],
                one_update_residual=Qinv*k['beta']*residual)
        v['W0_X_cap'], v['W1_X_cap'] = Bs+dfield, Bs+dfield+correction
    else:
        Bs,dfield,correction,E = inherited
    second = put('free_wave_second_jet_H_cap', displacement_curvature=u['u_rr'], velocity_gradient=u['v_r'])
    fzero = put('frozen_initial_forcing_H_cap', velocity=A0*u['v'], displacement_gradient=B0*u['u_r'], displacement=D0*u['u'])
    fx = put('frozen_forcing_x_H_cap', derivative=(A0+B0)*second, lower_order=D0*Bs)
    dx = put('W0_minus_Bs_time_derivative_X_cap', initial_acceleration=P*fzero,
             differentiated_bulk=P*T*fx, Robin_free_wave_derivative=P*M*Bs, affine_shell_derivative=P*k['Gx'])
    rzero = put('r0_initial_H_cap', velocity_difference=DA*u['v'], gradient_difference=DB*u['u_r'],
                value_difference=DD*u['u'], actual_source=k['F'])
    rx = put('r0_time_derivative_H_cap', correction_second=(A+B)*dx,
             correction_first=(D+k['comm'])*dfield, free_second_difference=(DA+DB)*second,
             free_first_difference=(DD+k['comm'])*Bs, original_source_derivative=k['Fx'])
    cx = put('correction_time_derivative_X_cap', initial_acceleration=P*rzero, differentiated_residual=P*T*rx)
    defect = put('reduced_defect_H_cap', first_order_correction=k['c']*correction)
    defectx = put('reduced_defect_time_derivative_H_cap', correction_second=(A+B)*cx,
                 correction_first=(D+k['comm'])*correction)
    w = put('initial_error_acceleration_H_cap', curvature=e['e0_rr'], velocity=A*e['v0'],
            gradient=B*e['e0_r'], value=D*e['e0'])
    j = max((e['v0']+T*w)/3,e['v0_r']+w)
    v['differentiated_initial_free_wave_X_cap'] = j
    z = put('second_error_common', initial=Qinv*j, differentiated_defect=Qinv*T*defectx,
            coefficient_commutator=Qinv*T*k['comm']*E)
    v['second_coordinate_error_H_caps'] = dict(xx=z,xr=z,rr=z+k['c']*E+defect)
    if inherited is None:
        forcing_initial = put('absolute_f1_initial', velocity=A*u['v'], gradient=B*u['u_r'],
                              value=D*u['u'], actual_source=k['F'])
        forcing_derivative = put('absolute_f1_time_derivative', principal=(A+B)*(second+dx),
                                 lower_and_commutator=(D+k['comm'])*v['W0_X_cap'], source_derivative=k['Fx'])
        acceleration = u['u_rr']+forcing_initial
        absj = max((u['v']+T*acceleration)/3,u['v_r']+acceleration)
        w1x = P*(absj+T*forcing_derivative+k['Gx'])
        v['absolute_W1_H_caps'] = dict(value=F(3,4)*v['W1_X_cap'],x=v['W1_X_cap'],r=v['W1_X_cap'],
                                     xx=w1x,xr=w1x,rr=w1x+k['c']*v['W0_X_cap']+k['F'])
        v['absolute_W1_prescribed_forcing_ledger'] = dict(f1_initial_H_cap=forcing_initial,
            f1_time_derivative_H_cap=forcing_derivative,initial_acceleration_H_cap=acceleration,
            differentiated_free_wave_X_cap=absj,W1_time_derivative_X_cap=w1x)
    return v, terms


def numeric_tree(value):
    if isinstance(value, dict):
        return {k:numeric_tree(v) for k,v in value.items()}
    return F(value)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--freeze',type=Path,default=Path(__file__).resolve().parent.parent/'CONSERVATIVE_RESPONSE_FREEZE.json');ap.add_argument('--input',action='append',default=[],help='portable role=path override; frozen SHA remains required');ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    require(not args.output.exists(),'FRESH_OUTPUT');freeze=json.loads(args.freeze.read_bytes());pre={Path(__file__):digest(__file__),args.freeze:digest(args.freeze)};reports={}
    overrides={}
    for item in args.input:
        require('=' in item,'INPUT_ROLE_PATH');role,path=item.split('=',1);require(role in freeze['inputs'] and role not in overrides,'INPUT_ROLE');overrides[role]=Path(path)
    for role,row in freeze['inputs'].items():
        p=overrides.get(role,Path(row['path']));require(p.is_file() and not p.is_symlink() and digest(p)==row['sha256'],'FROZEN_INPUT:'+role);pre[p]=row['sha256'];reports[role]=json.loads(p.read_bytes())
    retained=inputs(args.source,pre);k=constants(retained)
    base=reports['baseline'];old=reports['baseline_initial']['campaigns'][-1]
    absolute=scalar_caps(old['candidate_initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'old absolute')
    error=scalar_caps(old['initial_error_H_caps'],('e0','e0_r','e0_rr','v0','v0_r'),'old error')
    original,_=response_equations(k,absolute,error)
    nodes=[]
    for key,value in original.items():
        if key=='second_error_common' or key.startswith('absolute_f1_'):continue
        require(numeric_tree(base['ledger'][key])==value,'BASELINE_EXACT_NODE:'+key);nodes.append(key)
    own=reports['own_initial'];reg=reports['own_regularity']
    require(own['status']=='PASS_INDEPENDENT_EXACT_VELOCITY_INITIAL_GATES' and reg['status']=='PASS_INDEPENDENT_ANALYTIC_REGULARITY_AUDIT' and reports['own_audit']['status']=='PASS_INDEPENDENT_VELOCITY_REPORT_AND_PROVENANCE_AUDIT','INDEPENDENT_INITIAL_REGULARITY_AUDIT')
    require(reg['candidate']=='W1_exact_Cauchy_packet' and all(j==['0','0'] for j in reg['complete_component_velocity_affine_extension_jet_differences_through1'].values()),'COMPLETED_C1_VELOCITY')
    for key in own['source_sha256_before_and_after']:
        p=args.source/key;h=own['source_sha256_before_and_after'][key];require(digest(p)==h,'COMPLETE_SOURCE_CLOSURE:'+key);pre[p]=h
    data=scalar_caps(reports['absolute']['candidate_initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'refined absolute')
    data['v'],data['v_r']=map(F,own['coverage']['nominal_absolute_H_caps'])
    error['v0'],error['v0_r']=map(F,own['coverage']['error_H_caps'])
    values,terms=response_equations(k,data,error)
    gates={key:2*value<=F(base['ledger']['second_coordinate_error_H_caps'][key]) for key,value in values['second_coordinate_error_H_caps'].items()}
    gates['C1_no_increase']=values['complete_C1_error_X_cap']<=F(base['ledger']['complete_C1_error_X_cap'])
    require(all(gates.values()),'CONSERVATIVE_COMPLETE_RESPONSE_GATES')
    native=reports['native_response'];native_values=numeric_tree(native['ledger'])
    for key in ('u','u_r','u_rr'):require(data[key]==F(native['candidate_initial_absolute_H_caps'][key]),'SAME_DISPLACEMENT_ABSOLUTE:'+key)
    for key in ('e0','e0_r','e0_rr'):require(error[key]==F(native['initial_error_H_caps'][key]),'FROZEN_DISPLACEMENT_ERROR:'+key)
    for key in ('v','v_r'):require(data[key]>=F(native['candidate_initial_absolute_H_caps'][key]),'VALID_LARGER_INDEPENDENT_VELOCITY_ABSOLUTE:'+key)
    for key in ('v0','v0_r'):require(error[key]>=F(native['initial_error_H_caps'][key]),'VALID_LARGER_INDEPENDENT_VELOCITY_ERROR:'+key)
    require(values['second_coordinate_error_H_caps']['rr']>=native_values['second_coordinate_error_H_caps']['rr'],'CONSERVATIVE_MONOTONICITY')
    require(values['frozen_initial_forcing_H_cap']>0 and values['r0_initial_H_cap']>0,'BOTH_INITIAL_ACCELERATION_CHARGES_RETAINED')
    require(all(digest(p)==h for p,h in pre.items()),'SOURCE_AND_SELF_CLOSURE_AFTER')
    out={'status':'PASS_INDEPENDENT_CONSERVATIVE_EXACT_CAUCHY_RESPONSE','candidate':'W1_exact_Cauchy_packet','no_scientific_imports':True,'source_formula':'Independent PR109 named nonnegative response equations','baseline_exact_nodes':nodes,'constants':k,'candidate_initial_absolute_H_caps':data,'initial_error_H_caps':error,'values':values,'named_nonnegative_terms':terms,'gates_against_PR109':gates,'old_over_new_second_upper_ratios':{s:F(base['ledger']['second_coordinate_error_H_caps'][s])/v for s,v in values['second_coordinate_error_H_caps'].items()},'C1_old_over_new_upper_ratio':F(base['ledger']['complete_C1_error_X_cap'])/values['complete_C1_error_X_cap'],'uses_own_independent_larger_fresh_velocity_caps':True,'r1_initial_exact_cancellation':True,'all_f00_r00_full_Robin_original_Fx_commutator_feedback_inverse_terms_retained':True,'authenticated_files_before_after':{str(p):h for p,h in pre.items()},'scope':'Conditional coordinate upper certificate only, no global true-error factor, evolved point table, natural residual, physical observable or nonlinear conclusion.'}
    args.output.write_text(json.dumps(encode(out),indent=2,sort_keys=True)+'\n');print(json.dumps({'status':out['status'],'second_errors':{s:float(v) for s,v in values['second_coordinate_error_H_caps'].items()},'C1':float(values['complete_C1_error_X_cap']),'old_over_new':{s:float(v) for s,v in out['old_over_new_second_upper_ratios'].items()}}))
if __name__=='__main__':main()


