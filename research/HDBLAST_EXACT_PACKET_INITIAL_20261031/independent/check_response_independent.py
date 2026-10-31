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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--initial', required=True, type=Path)
    parser.add_argument('--initial-sha256', required=True)
    parser.add_argument('--producer-response', required=True, type=Path)
    parser.add_argument('--producer-response-sha256', required=True)
    parser.add_argument('--initial-source', required=True, type=Path)
    parser.add_argument('--initial-freeze', required=True, type=Path)
    parser.add_argument('--response-source', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    require(not args.output.exists(), 'fresh absent review output required')
    pre = {Path(__file__).resolve():digest(__file__)}
    canonical_pins = {'checker_self':digest(__file__)}
    old = inputs(args.repo_root, pre)
    canonical_pins.update({'research/'+relative:expected for relative,expected in PIN.values()})
    for path, expected in [(args.initial,args.initial_sha256),(args.producer_response,args.producer_response_sha256)]:
        require(digest(path) == expected, 'external report hash mismatch: '+path.name)
        pre[path] = expected
    canonical_pins.update(initial_report=args.initial_sha256,producer_response=args.producer_response_sha256)
    initial = json.loads(args.initial.read_text())
    producer = json.loads(args.producer_response.read_text())
    require(digest(args.initial_source) == initial['producer_sha256'], 'executed initial producer source hash')
    require(digest(args.initial_freeze) == initial['freeze_sha256'], 'initial frozen definition map hash')
    response_digests = [value for name,value in producer['input_sha256'].items() if name.endswith('response_ledger_oracle.py')]
    require(len(response_digests) == 1 and digest(args.response_source) == response_digests[0],
            'executed response producer source hash')
    pre[args.initial_source],pre[args.initial_freeze],pre[args.response_source] = \
        initial['producer_sha256'],initial['freeze_sha256'],response_digests[0]
    canonical_pins.update(initial_source=initial['producer_sha256'],initial_freeze=initial['freeze_sha256'],
                          response_source=response_digests[0])
    for relative, expected in initial['authentication_before_after'].items():
        path = args.repo_root / relative
        require(digest(path) == expected, 'initial selected source closure hash: '+relative)
        pre[path] = expected
        canonical_pins[relative] = expected
    k = constants(old)
    absolute,error,local = check_initial(initial,old)
    require(producer['candidate'] == 'W1_exact_packet' and
            producer['actual_target_original_F_Fx_G_M_unchanged'] is True and
            producer['boundary_error'] == {'g':'0','g_x':'0','M_x':'0'},
            'same original candidate source and affine error boundary law')
    require(producer['labels'] == old['old_C2']['labels'] and producer['Gram'] == old['old_C2']['Gram'],
            'same response modal norm')
    expected_constants = {name:k[name] for name in ('actual','frozen','T','q','F','Fx','G','Gx','KF','c','c0','dc','comm')}
    expected_constants.update(delta=k['difference'],derivative=k['derivatives'],bdry=k['boundary'],Pinv=k['P'],betaR=k['beta'])
    require(numeric_tree(producer['constants']) == expected_constants,
            'every reported operator source contraction and commutator constant')
    require(producer['new_point_or_quadrature_values_produced'] is False and
            not any(name in producer for name in ('full23_point_table','point_table','quadrature_point_values')),
            'no copied old point or quadrature table in ledger-only report')
    require(numeric_tree(producer['candidate_initial_absolute_H_caps']) == absolute and
            numeric_tree(producer['initial_error_H_caps']) == error, 'response takes fresh distinct absolute and error jets')
    values,terms = response_equations(k,absolute,error)
    for name,value in values.items():
        if name == 'second_error_common' or name.startswith('absolute_f1_'):
            continue
        require(numeric_tree(producer['ledger'][name]) == value, 'independent exact equation mismatch: '+name)
    dr = old['transfer']['direct_reproof']
    baseline,_ = response_equations(k,rational_dict(old['shared_absolute']['initial_absolute_H_caps']),
        rational_dict(old['old_C2']['initial_error_H_caps']),
        inherited=(F(dr['new_smooth_Bs_X_cap']),F(dr['residual_majorants']['W_minus_B0_whole_cone_X_cap']),
                   F(dr['new_actual_bulk_correction_X_cap']),F(old['old_C2']['complete_original_C1_error_X_cap'])))
    require(baseline['second_coordinate_error_H_caps'] ==
            rational_dict(old['old_C2']['original_target_second_coordinate_error_H_caps']),
            'independent old specialization reproduces every retained exact C2 cap')
    oldC1 = F(old['old_C2']['complete_original_C1_error_X_cap'])
    oldC2 = rational_dict(old['old_C2']['original_target_second_coordinate_error_H_caps'])
    gates = dict(C1_no_increase=values['complete_C1_error_X_cap'] <= oldC1)
    gates.update({name:2*value <= oldC2[name] for name,value in values['second_coordinate_error_H_caps'].items()})
    require(all(gates.values()) and producer['full_response_certificate_gates'] == gates,
            'exact predeclared complete response gates')
    scale = F(old['normalization']['target_scale_strict_lower'])
    require(scale == F(3249*10**30), 'retained global endpoint lower scale')
    normalization = dict(actual_target_strict_lower_scale=scale,
        initial_rr_upper_ratio=error['e0_rr']/scale,
        evolved_rr_upper_ratio=values['second_coordinate_error_H_caps']['rr']/scale,
        evolved_max_second_seminorm_upper_ratio=max(values['second_coordinate_error_H_caps'].values())/scale)
    require(numeric_tree(producer['global_normalization']) == normalization,
            'global ratios secondary and based on retained global endpoint scale')
    require(values['frozen_initial_forcing_H_cap'] > 0 and values['r0_initial_H_cap'] > 0,
            'nonzero correction acceleration charges retained')
    require(all(digest(path) == expected for path,expected in pre.items()),
            'selected sources maps reports and checker unchanged after evaluation')
    report = dict(status='PASS_INDEPENDENT_EXACT_PACKET_RESPONSE_LEDGER',
        no_scientific_imports=True,assertion_free=True,
        authenticated_file_count=len(pre),authentication_before_after_equal=True,
        source_input_sha256=canonical_pins,
        independent_values=values,named_nonnegative_terms=terms,local_true_error_gate=local,
        initial_absolute_caps=absolute,paired_initial_error_caps=error,
        exact_complete_response_gates=gates,old_exact_C1=oldC1,old_exact_C2=oldC2,
        global_normalization=normalization,retained_C2_specialization_reproduced=True,
        candidate='W1_exact_packet',
        boundary_error=dict(g='0',g_x='0',M_x='0'),
        common_first_Cauchy_jet_implies_r1_initial_zero=True,
        point_tables_inherited=False,
        scope='Conditional physical-side absolute coordinate coefficient estimates. Local true-error improvement uses its own positive all23 target denominator. Global factor-two compares complete upper certificates only; no new physical-constraint, shell-contact nonlinear existence or finite-amplitude claim.')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(report),sort_keys=True,indent=2)+'\n')


if __name__ == '__main__':
    main()
