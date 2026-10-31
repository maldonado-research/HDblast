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


def check_initial(initial, retained):
    require(initial['schema'] == 'HDBLAST_EXACT_PACKET_INITIAL_EXPERIMENT_V1'
            and initial['candidate'] == 'W1_exact_packet'
            and initial['status'] == 'INITIAL_AND_LOCAL_GATES_PASS'
            and initial['target_and_original_F_Fx_G_M_unchanged'] is True,
            'frozen exact packet initial input role')
    successful = [c for c in initial['campaigns'] if c['status'] == 'INITIAL_AND_LOCAL_GATES_PASS']
    require(len(successful) == 1, 'exactly one accepted initial campaign')
    c = successful[0]
    old = retained['old_C2']
    labels, weights = old['labels'], list(map(F, old['Gram']))
    require(c['labels'] == labels and c['Gram'] == old['Gram']
            and len(labels) == 23 and min(weights) > 0, 'complete positive all23 Gram')
    require(c['unresolved_cells'] == [] and c['initial_gate_pass'] is True,
            'no unresolved initial cells')
    require(c['bits'] in (128, 192) and c['leaf_count'] == len(c['leaves']) <= 4096,
            'declared finite campaign budget')
    absolute = scalar_caps(c['candidate_initial_absolute_H_caps'], ('u','u_r','u_rr','v','v_r'), 'absolute')
    error = scalar_caps(c['initial_error_H_caps'], ('e0','e0_r','e0_rr','v0','v0_r'), 'error')
    require(c['paired_channel_names'] == ['background_coefficient_difference','exact_shared_q_E_zero','actual_ODE_remainder'],
            'separate coefficient shared analytic and ODE channels')
    absolute_max = [F(0)] * 3
    paired_max = [[F(0)] * 3 for _ in range(3)]
    cursor = F(0)
    for leaf in c['leaves']:
        lo, hi = map(F, leaf['r'])
        require(lo == cursor and lo < hi and leaf['depth'] <= 16,
                'complete ordered leaf cover including tail')
        cursor = hi
        require(leaf['threshold_satisfied'] is True, 'every initial leaf accepted')
        for order in range(3):
            components = list(map(F, leaf['absolute_component_caps'][order]))
            cap = F(leaf['absolute_H_caps'][order])
            require(len(components) == 23 and min(components) >= 0 and
                    cap >= 0 and cap*cap >= sum(w*v*v for w,v in zip(weights,components)),
                    'leaf absolute positive Gram aggregation')
            absolute_max[order] = max(absolute_max[order], cap)
            for channel in range(3):
                vector = list(map(F, leaf['paired_channel_component_caps'][channel][order]))
                pcap = F(leaf['paired_channel_H_caps'][channel][order])
                require(len(vector) == 23 and min(vector) >= 0 and pcap >= 0 and
                        pcap*pcap >= sum(w*v*v for w,v in zip(weights,vector)),
                        'leaf paired positive Gram aggregation')
                if channel == 1:
                    require(pcap == 0 and all(v == 0 for v in vector),
                            'functional q and E identity retained exactly')
                paired_max[channel][order] = max(paired_max[channel][order], pcap)
    require(cursor == F(7,4), 'complete physical initial footprint')
    require([absolute[s] for s in ('u','u_r','u_rr')] == absolute_max,
            'fresh absolute displacement caps aggregated from full leaves')
    require(c['paired_channel_H_caps'] == encode(paired_max), 'global paired channel maxima')
    require([error[s] for s in ('e0','e0_r','e0_rr')] ==
            [sum(paired_max[ch][order] for ch in range(3)) for order in range(3)],
            'all displacement errors retain every background channel')
    olderr = rational_dict(old['initial_error_H_caps'])
    velocity = c['velocity']
    require(error['v0'] == olderr['v0'] and error['v0_r'] == olderr['v0_r'] and
            velocity['error_ledger_retained_exactly'] is True and
            scalar_caps(velocity['error_H_caps'], ('v0','v0_r'), 'velocity errors') ==
            {s:error[s] for s in ('v0','v0_r')}, 'unchanged complete velocity error charges')
    require(scalar_caps(velocity['absolute_H_caps'], ('v','v_r'), 'absolute velocity') ==
            {s:absolute[s] for s in ('v','v_r')} and F(velocity['support_zero_through']) > 1,
            'unchanged absolute velocity and support collar')
    require(error['e0'] <= olderr['e0'] and error['e0_r'] <= olderr['e0_r'] and
            2*error['e0_rr'] <= F(retained['old_combined']['refined_uniform_initial_u_rr_H_upper']),
            'all exact retained initial gates')
    local = c['local']
    require(local['point'] == ['0','65/64'] and len(local['rows']) == 23,
            'same local complete target error witness')
    sqlo = sqhi = sqerror = F(0)
    for label, w, row in zip(labels, weights, local['rows']):
        require(row['label'] == label and F(row['Gram']) == w, 'local complete Gram slot')
        lo, hi = map(F, row['actual_u_rr'])
        elo, ehi = map(F, row['paired_actual_minus_candidate_u_rr'])
        require(lo <= hi and elo <= ehi, 'ordered actual and paired local intervals')
        lower = F(0) if lo <= 0 <= hi else min(abs(lo),abs(hi))
        sqlo += w*lower**2
        sqhi += w*max(abs(lo),abs(hi))**2
        sqerror += w*max(abs(elo),abs(ehi))**2
    dl, du = map(F, local['actual_H_norm_interval'])
    numerator = F(local['new_error_H_upper'])
    oldlower = F(retained['normalization']['local_same_H_calibration']['actual_error_H_norm_strict_lower'])
    require(0 < dl <= du and dl*dl <= sqlo and du*du >= sqhi and
            numerator >= 0 and numerator*numerator >= sqerror,
            'independently aggregate local actual lower and new error upper')
    require(F(local['old_error_H_strict_lower']) == oldlower and 10*numerator <= dl
            and 2*numerator < oldlower, 'local true error improvement exact inequalities')
    require(F(local['new_relative_error_H_upper']) == numerator/dl,
            'local denominator from same point complete H norm')
    return absolute, error, dict(target_lower=dl, target_upper=du,
        error_upper=numerator, old_error_strict_lower=oldlower,
        new_relative_upper=numerator/dl, actual_squared_lower=sqlo,
        actual_squared_upper=sqhi, error_squared_upper=sqerror)


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



FREEZE_SHA='6cb41096a36598a130dda1d19bec14c53c33ac723aec2a1e6638389d4308d751'
PRODUCER_SHA='449d1d0293cfe2f1b8620ea8827a6ffa4343731c60c62b26f5acb94dcbf0639c'
BASELINE_INITIAL_SHA='692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3'
BASELINE_RESPONSE_SHA='991ee9e04d17a879e4e65f3f68ad8f9e76929b4c10d4698d2af06dd6588f3273'
EQUATIONS_SHA='48b112a26788f5e715ff439ad10b2166ddebae9c1e65a25c3626d1d3b1c6ad4c'
PACK='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/'


def sha_regular(path):
    p=Path(path)
    require(p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.parents),'regular symlink-free input: '+p.name)
    return hashlib.sha256(p.read_bytes()).hexdigest()


def local_squares(local, labels, gram, own=False):
    require(local['point']==['0','65/64'] and len(local['rows'])==23,'fresh all23 velocity local point')
    targetlo=targethi=errorhi=oldlo=F(0)
    for name,g,row in zip(labels,gram,local['rows']):
        require(row['label']==name and F(row['Gram'])==g,'fresh local all23 positive Gram')
        def span(key):
            lo,hi=map(F,row[key]);require(lo<=hi,'fresh local interval order')
            return (0 if lo<=0<=hi else min(abs(lo),abs(hi))),max(abs(lo),abs(hi))
        al,ah=span('actual_v_r');el,eh=span('paired_error_v_r')
        ol,oh=span('actual_minus_old_smoothed_v_r' if own else 'actual_minus_old_v_r')
        targetlo+=g*al**2;targethi+=g*ah**2;errorhi+=g*eh**2;oldlo+=g*ol**2
    if own:
        low,up=map(F,local['actual_v_r_H_interval']);new=F(local['new_error_v_r_H_interval'][1]);old=F(local['old_error_v_r_H_interval'][0])
    else:
        low,up=map(F,local['actual_H_norm_interval']);new=F(local['new_error_H_upper']);old=F(local['old_error_H_lower'])
    require(0<low<=up and low**2<=targetlo and up**2>=targethi and new>=0 and new**2>=errorhi and old>0 and old**2<=oldlo,'fresh local independent square aggregation')
    require(10*new<=low and 2*new<old,'fresh local true-error upper versus lower gates')
    return dict(actual_target_H_interval=[low,up],new_error_H_upper=new,old_error_H_lower=old,
      new_over_target_upper=new/low,actual_squared_lower=targetlo,actual_squared_upper=targethi,
      new_error_squared_upper=errorhi,old_error_squared_lower=oldlo)


def public_values(values):
    return {name:value for name,value in values.items() if name!='second_error_common' and not name.startswith('absolute_f1_')}


def flat_public(value,prefix=''):
    out={}
    for name,v in value.items():
        key=prefix+'/'+name if prefix else name
        if isinstance(v,dict):out.update(flat_public(v,key))
        else:out[key]=F(v)
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-root',required=True,type=Path)
    for role in ('same-response','absolute','velocity','independent-velocity'):
        ap.add_argument('--'+role+'-root',required=True,type=Path)
    ap.add_argument('--producer-source',required=True,type=Path)
    ap.add_argument('--producer-freeze',required=True,type=Path)
    ap.add_argument('--producer-report',required=True,type=Path)
    ap.add_argument('--producer-report-sha256',required=True)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    require(not args.output.exists(),'fresh independent Cauchy output')
    before={Path(__file__).resolve():sha_regular(__file__)}
    portable={'checker.py':before[Path(__file__).resolve()]}
    def pin(path,expected,label):
        rel=Path(label);require(not rel.is_absolute() and '..' not in rel.parts,'portable selected source role')
        require(sha_regular(path)==expected,'selected source hash: '+label)
        before[Path(path)]=expected;portable[label]=expected
    pin(args.producer_source,PRODUCER_SHA,'producer/response_cauchy.py')
    pin(args.producer_freeze,FREEZE_SHA,'producer/FREEZE.json')
    pin(args.producer_report,args.producer_report_sha256,'producer/RESPONSE.json')
    freeze=json.loads(args.producer_freeze.read_text())
    require(freeze['candidate']=='W1_exact_Cauchy_packet' and freeze['adapter_sha256']==PRODUCER_SHA,'separate Cauchy frozen role')
    expected_producer_sources={'response_cauchy.py':PRODUCER_SHA,'FREEZE.json':FREEZE_SHA}
    roots={role:getattr(args,role+'_root') for role in ('same_response','absolute','velocity','independent_velocity')}
    for role,pins in freeze['input_pins'].items():
        require(role in roots,'known proof root role')
        for rel,expected in pins.items():
            r=Path(rel);require(not r.is_absolute() and '..' not in r.parts,'proof member relative')
            pin(roots[role]/r,expected,role+'/'+rel)
            expected_producer_sources[role+'/'+rel]=expected
    oldfreeze=json.loads((roots['same_response']/'checks/FREEZE.json').read_text())
    for rel,expected in oldfreeze['source_pins'].items():
        pin(args.source_root/rel,expected,rel);expected_producer_sources[rel]=expected
    for relative,expected in [(PACK+'reference/INITIAL.json',BASELINE_INITIAL_SHA),(PACK+'reference/RESPONSE.json',BASELINE_RESPONSE_SHA),(PACK+'independent/check_response_independent.py',EQUATIONS_SHA)]:
        pin(args.source_root/relative,expected,relative)
    for relative,expected in PIN.values():pin(args.source_root/'research'/relative,expected,'research/'+relative)
    native=json.loads((roots['velocity']/'INITIAL_NORMAL.json').read_text())
    independent=json.loads((roots['independent_velocity']/'INITIAL_NORMAL.json').read_text())
    regularity=json.loads((roots['independent_velocity']/'REGULARITY_FINAL_NORMAL.json').read_text())
    for report,key in [(native,'source_authentication_before_after'),(independent,'source_sha256_before_and_after'),(regularity,'source_sha256_before_and_after')]:
        for rel,expected in report[key].items():
            require(rel.startswith('research/'),'original scientific source role')
            pin(args.source_root/rel,expected,rel);expected_producer_sources[rel]=expected
    # All selected scientific sources and both proof programs are authenticated.
    # No scientific module imports occur: the distinct named-term equations
    # are retained textually from the pinned independent source above.
    retained=inputs(args.source_root,before);k=constants(retained)
    oldinitial=json.loads((args.source_root/PACK/'reference/INITIAL.json').read_text())
    oldresponse=json.loads((args.source_root/PACK/'reference/RESPONSE.json').read_text())
    olddata,olderror,oldlocal=check_initial(oldinitial,retained)
    baseline_values,baseline_terms=response_equations(k,olddata,olderror)
    require(numeric_tree(oldresponse['ledger'])==public_values(baseline_values),'every PR109 baseline ledger node exactly reproduced')
    expected_constants={name:k[name] for name in ('actual','frozen','T','q','F','Fx','G','Gx','KF','c','c0','dc','comm')}
    expected_constants.update(delta=k['difference'],derivative=k['derivatives'],bdry=k['boundary'],Pinv=k['P'],betaR=k['beta'])
    require(numeric_tree(oldresponse['constants'])==expected_constants,'every baseline constant exactly reproduced')
    passed=[c for c in native['campaigns'] if c['status']=='INITIAL_AND_LOCAL_GATES_PASS']
    require(native['candidate']==independent['candidate']==regularity['candidate']=='W1_exact_Cauchy_packet' and len(passed)==1,'distinct native and independent velocity identity')
    c=passed[0]
    require(native['status']=='INITIAL_AND_LOCAL_GATES_PASS' and c['initial_gate_pass'] is True and c['unresolved_cells']==[],'native complete initial/local gates')
    require(independent['status']=='PASS_INDEPENDENT_EXACT_VELOCITY_INITIAL_GATES' and independent['global_coverage_performed'] is True and all(independent['coverage']['initial_gates_pass']),'independent whole initial gates')
    require(native['target_F_Fx_G_M_H_domain_unchanged'] is True and native['old_candidate_functions_unchanged'] is True,'same actual target and older candidate retained')
    require(c['labels']==oldresponse['labels'] and c['Gram']==oldresponse['Gram'],'complete unchanged all23 modal norm')
    labels=c['labels'];gram=list(map(F,c['Gram']))
    native_local=local_squares(c['local'],labels,gram)
    independent_local=local_squares(independent['local'],labels,gram,True)
    require(c['local']['old_point_table_reused'] is False and c['local']['fresh_direct_actual_reconstruction'] is True,'fresh local witness reconstruction')
    native_data=scalar_caps(c['candidate_initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'fresh native velocity absolute')
    error=scalar_caps(c['initial_error_H_caps'],('e0','e0_r','e0_rr','v0','v0_r'),'fresh native paired error')
    require(all(error[name]==olderror[name] for name in ('e0','e0_r','e0_rr')),'exactly frozen PR109 displacement errors')
    require(all(error[name]<=olderror[name] for name in ('v0','v0_r')),'fresh velocity error gates')
    refinement=json.loads((roots['absolute']/'REFINEMENT.json').read_text())
    data=scalar_caps(refinement['candidate_initial_absolute_H_caps'],('u','u_r','u_rr','v','v_r'),'refined displacement')
    for name in ('v','v_r'):data[name]=native_data[name]
    require(all(data[name]<=olddata[name] for name in ('u','u_r','u_rr')),'same displacement no worse refined allowances')
    if 'COMBINED_INITIAL_NORMAL.json' in freeze['input_pins']['velocity']:
        combined=json.loads((roots['velocity']/'COMBINED_INITIAL_NORMAL.json').read_text())
        require(numeric_tree(combined['candidate_initial_absolute_H_caps'])==data and numeric_tree(combined['initial_error_H_caps'])==error,'fresh combined initial provenance')
    audit_name=next(name for name in freeze['input_pins']['independent_velocity'] if 'AUDIT' in name and name.endswith('.json'))
    audit=json.loads((roots['independent_velocity']/audit_name).read_text())
    require(audit['status']=='PASS_INDEPENDENT_VELOCITY_REPORT_AND_PROVENANCE_AUDIT','native initial table independently audited')
    require(regularity['status']=='PASS_INDEPENDENT_ANALYTIC_REGULARITY_AUDIT' and regularity['required_q_order_for_analytic_velocity_r']==4,'complete analytic velocity regularity inventory')
    require([row['label'] for row in c['endpoint_componentwise_affine_matching']]==labels and all(row['exterior_minus_interior_exact_jet']==[0,0] for row in c['endpoint_componentwise_affine_matching']),'all23 native affine value/r matching')
    require(set(regularity['complete_component_velocity_affine_extension_jet_differences_through1'])==set(labels) and all(v==['0','0'] for v in regularity['complete_component_velocity_affine_extension_jet_differences_through1'].values()),'all23 independent affine value/r matching')
    require(independent['coverage']['exact_zero_region']==[0,1] and set(regularity['velocity_parity'])==set(labels) and all(abs(v)==1 for v in regularity['velocity_parity'].values()),'support collar and all23 H-isometric parity')
    values,terms=response_equations(k,data,error)
    producer=json.loads(args.producer_report.read_text())
    require(producer['schema']=='HDBLAST_SEPARATE_EXACT_CAUCHY_RESPONSE_V1' and producer['status']=='PASS_SEPARATE_EXACT_CAUCHY_CONDITIONAL_RESPONSE','advertised response schema and status')
    require(producer['candidate']=='W1_exact_Cauchy_packet' and producer['same_function_as_W1_exact_packet'] is False and producer['actual_target_original_F_Fx_G_M_operators_domain_H_unchanged'] is True,'distinct Cauchy response candidate role')
    require(numeric_tree(producer['constants'])==expected_constants,'every Cauchy operator/source/commutator constant')
    require(producer['selected_input_sha256_before']==producer['selected_input_sha256_after']==expected_producer_sources,'every advertised selected source pin')
    for name in ('older_candidate_functions_and_certificates_preserved','fresh_analytic_velocity_absolute_and_paired_error_caps_used','old_P1_velocity_caps_not_used_for_new_candidate','native_and_independent_initial_local_and_regularity_gates_verified','no_new_evolved_point_table_natural_constraint_residual_or_observable','no_global_true_error_factor_claim','selected_inputs_unchanged'):
        require(producer[name] is True,'advertised Cauchy boolean flag: '+name)
    require(producer['labels']==labels and producer['Gram']==c['Gram'] and producer['domain']==oldresponse['domain'] and producer['baseline_exact_node_count']==35,'advertised response domain/full23 norm/baseline count')
    expected_roles={'u/u_r/u_rr':'same displacement function, refined absolute caps','e0/e0_r/e0_rr':'exact frozen published PR109 displacement errors','v/v_r':'fresh native analytic nominal velocity absolute caps','v0/v0_r':'fresh native paired actual-minus-analytic-nominal velocity errors'}
    require(producer['field_roles']==expected_roles,'every advertised initial field role')
    require(producer['initial_cancellation']=='r1(0)=0 from common first Cauchy jet; f00/r00 initial acceleration charges retained','advertised common first-jet cancellation')
    require(producer['regularity']=='Completed displacement C2; new completed velocity C1 globally with exact component affine outer extension and H-isometric parity; flat through1. Conditional physical-side C2 response under retained coefficient/source/corner premises.','advertised conditional completed-component regularity')
    require(producer['scope']=='Distinct new exact-Cauchy candidate. Global comparisons are conditional coordinate upper certificates; the separately fresh local initial velocity witness supplies the local true-error statement. No inherited old point table, nonlinear shell-contact result, invariant observable, or cosmological claim.','advertised scientific scope')
    require(numeric_tree(producer['candidate_initial_absolute_H_caps'])==data,'fresh analytic absolute velocity inputs required')
    require(numeric_tree(producer['initial_error_H_caps'])==error,'fresh paired analytic velocity errors required')
    require(numeric_tree(producer['ledger'])==public_values(values),'every Cauchy ledger node independently reproduced')
    baseline_flat=flat_public(public_values(baseline_values));new_flat=flat_public(public_values(values))
    comparisons={name:dict(before=old,after=new_flat[name],old_over_new=old/new_flat[name] if new_flat[name] else None) for name,old in baseline_flat.items()}
    require(producer['node_comparisons']==encode(comparisons),'every advertised baseline-to-new node comparison')
    def gates(v):
        out={name:2*value<=F(oldresponse['ledger']['second_coordinate_error_H_caps'][name]) for name,value in v['second_coordinate_error_H_caps'].items()}
        out['C1_no_increase']=v['complete_C1_error_X_cap']<=F(oldresponse['ledger']['complete_C1_error_X_cap']);return out
    actual_gates=gates(values);require(all(actual_gates.values()) and producer['full_response_gates_against_PR109']==actual_gates,'all second errors half PR109 and C1 no worse')
    # A separate robustness proof uses the own independent velocity enclosures,
    # retaining the same refined displacement inputs and frozen displacement errors.
    own_data,own_error=dict(data),dict(error)
    for name,val in zip(('v','v_r'),independent['coverage']['nominal_absolute_H_caps']):own_data[name]=F(val)
    for name,val in zip(('v0','v0_r'),independent['coverage']['error_H_caps']):own_error[name]=F(val)
    own_values,own_terms=response_equations(k,own_data,own_error)
    own_gates=gates(own_values);require(all(own_gates.values()),'independent velocity caps robustly pass full response gates')
    # Mandatory stale-velocity controls have exact nonzero numeric witnesses.
    stale={}
    for role in ('absolute','error'):
        d,e=dict(data),dict(error)
        if role=='absolute':
            for name in ('v','v_r'):d[name]=olddata[name]
        else:
            for name in ('v0','v0_r'):e[name]=olderror[name]
        wrong,_=response_equations(k,d,e);delta=wrong['second_coordinate_error_H_caps']['rr']-values['second_coordinate_error_H_caps']['rr']
        require(delta!=0,'effective stale velocity '+role);stale[role]={'wrong_minus_correct_rr':delta,'rejected':True}
    A,B=k['actual']['A'],k['actual']['B'];P,T,q,comm=k['P'],k['T'],k['q'],k['comm']
    expected_witnesses={
      'reuse_old_velocity_absolute':dict(wrong_minus_correct_rr=stale['absolute']['wrong_minus_correct_rr'],rejected=True),
      'reuse_old_velocity_errors':dict(wrong_minus_correct_rr=stale['error']['wrong_minus_correct_rr'],rejected=True),
      'omit_f00':dict(wrong_minus_correct_rr=-T*T*P*P*(A+B)**2*values['frozen_initial_forcing_H_cap']/(1-q),rejected=True),
      'omit_r00':dict(wrong_minus_correct_rr=-T*P*(A+B)*values['r0_initial_H_cap']/(1-q),rejected=True),
      'omit_full_Robin_source':dict(wrong_minus_correct_rr=-T*T*P*P*(A+B)**2*(k['actual']['M']*values['free_wave_X_cap']+k['Gx'])/(1-q),rejected=True),
      'omit_original_Fx':dict(wrong_minus_correct_rr=-T*T*P*(A+B)*k['Fx']/(1-q),rejected=True),
      'omit_derivative_commutator':dict(wrong_minus_correct_rr=-T*comm*(values['bulk_correction_X_cap']+values['complete_C1_error_X_cap'])/(1-q),rejected=True),
      'lose_r1_initial_cancellation':dict(wrong_minus_correct_rr=values['r0_initial_H_cap']/(1-q),rejected=True),
      'omit_feedback_inverse':dict(wrong_minus_correct_rr=-q*values['second_coordinate_error_H_caps']['xx'],rejected=True)}
    require(all(row['wrong_minus_correct_rr']!=0 for row in expected_witnesses.values()),'effective independent formula corruption witnesses')
    require(producer['formula_and_stale_velocity_corruption_witnesses']==encode(expected_witnesses),'every advertised formula/stale-velocity corruption witness')
    require(all(sha_regular(path)==expected for path,expected in before.items()),'all selected sources reports checker postauthentication')
    output=dict(status='PASS_INDEPENDENT_EXACT_CAUCHY_RESPONSE',candidate='W1_exact_Cauchy_packet',assertion_free=True,no_scientific_imports=True,
      selected_source_sha256=portable,selected_source_count=len(before),all_selected_sources_authenticated_before_after=True,
      prior_independent_equation_source_sha256=EQUATIONS_SHA,every_PR109_baseline_ledger_node_reproduced=True,
      every_Cauchy_ledger_node_reproduced=True,public_ledger_node_count=35,constants=expected_constants,
      initial_absolute_H_caps=data,initial_error_H_caps=error,independent_values=values,named_nonnegative_terms=terms,
      complete_response_gates=actual_gates,native_local_true_error_gate=native_local,independent_local_true_error_gate=independent_local,
      own_independent_velocity_robustness=dict(initial_absolute_H_caps=own_data,initial_error_H_caps=own_error,
        independent_values=own_values,named_nonnegative_terms=own_terms,complete_response_gates=own_gates),
      stale_velocity_exact_corruption_witnesses=stale,formula_corruption_witnesses=expected_witnesses,
      common_first_Cauchy_jet_implies_r1_initial_zero=True,boundary_error=dict(g='0',g_x='0',M_x='0'),
      scope='Distinct exact Cauchy velocity candidate with same displacement and original actual target/source/operators/domain/H. Global factor comparisons concern coordinate upper certificates only. Fresh local velocity upper-versus-lower comparisons retain their own positive target denominators. Conditional exact initial/precontact hypotheses retained; nonlinear shell contact separate.')
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(encode(output),sort_keys=True,indent=2)+'\n')


if __name__=='__main__':main()
