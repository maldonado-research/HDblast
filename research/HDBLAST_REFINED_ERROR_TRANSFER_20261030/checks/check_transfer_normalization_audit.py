#!/usr/bin/env python3
"""Independent combined-ledger, affine-transfer and coordinate-ratio audit.

No producer or transfer checker code is imported. New raw derivative bounds
are premises supplied by separately audited reports; this checker rebuilds
their complete combined initial error and checks the independent C2 delta.
"""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


PINS = {
    'higher_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'endpoint_audit': 'e21d2cb172893495b557bbcdcee198d33d21a7ba287c5215d7319eff3d9dc9cc',
    'PR103_initial': '915acf8abf3c3bb2c843eac509315ebdbd568d7c8b10ff2682c900d672e7d4d7',
    'PR103_C2': 'fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0',
    'response': 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
    'PR103_joins': 'e1d0d6e3d8650cbf402ccb4f5a87fd20c20db651e3bd323c68f0bdea0219dd35',
}


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    for name in ('paired', 'nominal', 'initial', 'transfer'):
        ap.add_argument('--'+name+'-report', type=Path, required=True)
        ap.add_argument('--'+name+'-sha256', required=True)
    ap.add_argument('--normalization-report', type=Path)
    ap.add_argument('--normalization-sha256')
    a = ap.parse_args()
    checks, controls, files, hashes = [], [], {}, {}

    def need(ok, why):
        if not ok:
            raise RuntimeError(why)
        checks.append(why)

    def load(path, role, digest):
        raw = Path(path).read_bytes()
        need(hashlib.sha256(raw).hexdigest() == digest, 'authenticate '+role)
        files[role] = (Path(path), raw)
        hashes[role] = digest
        return json.loads(raw)

    def reject(name, fn):
        n = len(checks)
        caught = False
        try:
            fn()
        except (RuntimeError, KeyError, TypeError, ValueError, IndexError):
            caught = True
        del checks[n:]
        need(caught, 'reject '+name)
        controls.append(name)

    root = a.repo_root
    hi = root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference'
    oldroot = root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference'
    h = load(hi/'INITIAL.json', 'higher_initial', PINS['higher_initial'])
    endpoint = load(hi/'ENDPOINT.json', 'endpoint_audit', PINS['endpoint_audit'])
    old_i = load(oldroot/'INITIAL.json', 'PR103_initial', PINS['PR103_initial'])
    old = load(oldroot/'C2.json', 'PR103_C2', PINS['PR103_C2'])
    joins = load(oldroot/'JOINS.json', 'PR103_joins', PINS['PR103_joins'])
    response = load(root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json',
                    'response', PINS['response'])
    reports = {name: load(getattr(a, name+'_report'), name, getattr(a, name+'_sha256'))
               for name in ('paired', 'nominal', 'initial', 'transfer')}
    p, n, i, t = [reports[k] for k in ('paired', 'nominal', 'initial', 'transfer')]
    labels, gram = old['labels'], old['Gram']
    weights = list(map(Q, gram))
    need(len(labels) == 23 and len(gram) == 23 and min(weights) > 0,
         'complete positive full23 comparison norm')
    need(labels == h['labels'] == old_i['labels'] and gram == h['gram'] == old_i['Gram'],
         'same historical initial and evolved norm')
    for role, d in reports.items():
        need(d['status'].startswith('PASS') and d['labels'] == labels and d['Gram'] == gram,
             role+' same successful complete norm')
        need(d['candidate_changed'] is False, role+' unchanged candidate')
    need(t['actual_sources_changed'] is False and i['initial_data_unchanged'] is True and
         i['target_unchanged'] is True and i['F_G_and_coefficients_unchanged'] is True,
         'same actual target source and initial functions')
    need(i['refined_u_rr_H_cap_is_upper_bound'] is True and
         i['bound_quantity'] == 'sup_H_norm_actual_minus_smoothed_initial_u_rr',
         'combined input is a uniform upper initial error bound')
    need(n['nominal_M3_upper'] is True and n['actual_order_six_bound_claimed'] is False,
         'nominal derivative bounds require no actual sixth-order premise')
    need(i['input_sha256']['paired'] == a.paired_sha256 and
         i['input_sha256']['nominal'] == a.nominal_sha256 and
         t['input_sha256']['initial_refinement'] == a.initial_sha256,
         'combined and transfer reports retain exact new input pins')
    du = Q(1, 2**30)
    need(Q(i['fixed_displacement_kernel_width']) == Q(p['fixed_displacement_kernel_width']) ==
         Q(n['fixed_displacement_kernel_width']) == du and Q(i['ell']) == Q(h['ell']) == Q(1, 10**6),
         'same smoothing width and scaled coordinate convention')
    m3 = list(map(Q, n['refined_nominal_M3u_component_caps']))
    need(len(m3) == 23 and min(m3) >= 0 and
         m3 == list(map(Q, i['refined_nominal_M3u_component_caps'])), 'same refined nominal derivative vector')
    jumps = h['curvature_jumps']
    need([Q(j['r']) for j in jumps] == [Q(65,64), Q(66,64), Q(7,4)] and
         all(Q(v) == 0 for v in jumps[0]['jump_abs_component_caps']), 'same complete joins and zero first nominal jump')
    inner = [max(Q(j['jump_abs_component_caps'][k]) for j in jumps[:2]) for k in range(23)]
    need(inner == list(map(Q, i['retained_inner_jump_component_caps'])), 'same retained inner jump allowances')
    smooth = [Q(35,128)*du*d+j/2 for d,j in zip(m3,inner)]
    need(smooth == list(map(Q, i['refined_curvature_smoothing_component_caps'])) ==
         list(map(Q,n['refined_curvature_smoothing_component_caps'])), 'independent component convolution ledger')
    smooth_h = Q(i['refined_curvature_smoothing_H_cap'])
    need(smooth_h == Q(n['refined_curvature_smoothing_H_cap']) and
         smooth_h**2 >= sum((w*s*s for w,s in zip(weights,smooth)), Q(0)),
         'positive full-Gram smoothing norm enclosed outward')
    channels = [[Q(v) for v in ch] for ch in p['paired_channel_H_caps']]
    need(len(channels) == 3 and all(len(ch) == 2 and min(ch) >= 0 for ch in channels) and
         channels == [[Q(v) for v in ch] for ch in i['paired_channel_H_caps']], 'all three paired channel suprema retained')
    paired_h = sum((ch[0] for ch in channels), Q(0))
    need(paired_h == Q(p['paired_u_rr_H_upper']) == Q(i['paired_u_rr_H_upper']),
         'paired channel sum independently reconstructed')
    u2 = paired_h+smooth_h
    need(u2 == Q(i['refined_uniform_initial_u_rr_H_upper']) and
         0 < u2 < Q(old_i['refined_uniform_initial_u_rr_H_upper']), 'complete combined initial error bound improves')
    need(Q(i['unchanged_full_initial_v_r_error_H_upper']) == Q(old_i['unchanged_full_initial_v_r_error_H_upper']),
         'complete original velocity-gradient error retained')
    need(t['candidate'] == old['candidate'] and t['target'] == old['target'] and
         t['domain'] == old['domain'] and t['norm'] == old['norm'], 'same target domain norm and numerical candidate')
    olddata = {k:Q(v) for k,v in old['initial_error_H_caps'].items()}
    data = {k:Q(v) for k,v in t['initial_error_H_caps'].items()}
    need(data['e0_rr'] == u2 and all(data[k] == olddata[k] for k in olddata if k != 'e0_rr'),
         'only initial curvature error cap replaced')
    q = Q(response['retained_coupled_X_operator_cap'])
    need(q == Q(3409,8000000), 'same retained feedback contraction')
    acceleration = Q(old['initial_differentiated_acceleration_H_cap'])-olddata['e0_rr']+u2
    J = max((data['v0']+Q(3,2)*acceleration)/3, data['v0_r']+acceleration)
    delta = (J-Q(old['differentiated_initial_free_wave_X_cap']))/(1-q)
    need(acceleration == Q(t['initial_differentiated_acceleration_H_cap']) and
         J == Q(t['differentiated_initial_free_wave_X_cap']), 'independent differentiated initial max ledger')
    evolved = {}
    for key in ('xx','xr','rr'):
        value = Q(old['original_target_second_coordinate_error_H_caps'][key])+delta
        need(value == Q(t['original_target_second_coordinate_error_H_caps'][key]) and value > 0,
             'independent exact affine C2 delta '+key)
        need(value < Q(t['original_target_second_coordinate_error_H_readable_strict_uppers'][key]),
             'outward readable new C2 bound '+key)
        evolved[key] = value
    retained = ('actual_original_bulk_Fx_H_cap','actual_Fx_contribution_to_second_error_upper',
                'complete_original_C1_error_X_cap','complete_original_reduced_defect_H_cap',
                'original_reduced_defect_time_derivative_H_cap','boundary_error')
    need(all(t[k] == old[k] for k in retained), 'all source correction C1 and boundary ledgers retained exactly')
    need(labels[12] == 'w:0' and weights[12] == 1, 'actual scalar mean has unit Gram weight')
    scale = Q('3.249e33')
    actual = list(map(Q,h['endpoint']['rows'][12]['actual_u_rr']))
    independent = list(map(Q,endpoint['rows'][12]['actual_u_rr']))

    def endpoint_scale(d, label):
        need(len(d) == 2 and d[0] <= d[1] < -scale, label)

    endpoint_scale(actual,'paired actual endpoint establishes strict scale')
    endpoint_scale(independent,'independent actual ODE endpoint establishes strict scale')
    need(max(actual[0],independent[0]) <= min(actual[1],independent[1]), 'independent actual endpoint intervals overlap')
    need(h['endpoint']['rows'][12]['label'] == endpoint['rows'][12]['label'] == 'w:0' and
         Q(h['endpoint']['r']) == Q(7,4), 'denominator uses actual target endpoint row')
    x, r = Q(0), Q(7,4)
    need(0 <= x <= Q(3,2) and max(0,x-Q(5,4)) <= r <= Q(7,4)-x,
         'actual initial endpoint lies in declared physical Dplus')
    ratios = {'initial_rr_H_sup':u2/scale, 'evolved_rr_H_sup':evolved['rr']/scale,
              'evolved_max_second_derivative_H_sup_seminorm':max(evolved.values())/scale}
    baseline = {'initial_rr_H_sup':Q(old_i['refined_uniform_initial_u_rr_H_upper'])/scale,
        'evolved_rr_H_sup':Q(old['original_target_second_coordinate_error_H_caps']['rr'])/scale,
        'evolved_max_second_derivative_H_sup_seminorm':max(map(Q,old['original_target_second_coordinate_error_H_caps'].values()))/scale}
    need(all(0 < v < baseline[k] for k,v in ratios.items()), 'each exact coordinate global ratio cap improves')
    need(baseline['initial_rr_H_sup'] < Q('3.648e-6') and
         baseline['evolved_max_second_derivative_H_sup_seminorm'] < Q('3.654e-6'),
         'old frozen baseline normalization arithmetic independently verified')
    join = joins['joins'][0]
    row = join['rows'][12]
    need(Q(join['r']) == Q(65,64) and row['label'] == 'w:0' and
         joins['initial_candidate_unchanged'] is True and joins['actual_target_unchanged'] is True,
         'local calibration uses the same fixed candidate actual target and exact first join')
    local_target = list(map(Q,row['actual_u_rr']))
    local_error = list(map(Q,row['actual_minus_smoothed_u_rr']))
    need(0 < local_target[0] <= local_target[1] and
         0 < local_error[0] <= local_error[1], 'local scalar component denominators and error exclude zero')
    local_lower = local_error[0]/local_target[1]
    need(local_lower > Q(999,1000), 'local component relative error exceeds 999 per thousand')
    need(joins['labels'] == labels and joins['positive_Gram'] == gram and len(join['rows']) == 23,
         'local full-H calibration has exactly the same complete norm')
    actual_caps, error_mins = [], []
    for k,row_i in enumerate(join['rows']):
        values = list(map(Q,row_i['actual_u_rr']))
        need(row_i['label'] == labels[k] and values[0] <= values[1], 'ordered local actual interval '+labels[k])
        actual_caps.append(max(map(abs,values)))
        ev = list(map(Q,row_i['actual_minus_smoothed_u_rr']))
        need(ev[0] <= ev[1], 'ordered local error interval '+labels[k])
        error_mins.append(Q(0) if ev[0] <= 0 <= ev[1] else min(map(abs,ev)))
    actual_h_square_upper = sum((w*b*b for w,b in zip(weights,actual_caps)), Q(0))
    error_h_square_lower = sum((w*b*b for w,b in zip(weights,error_mins)), Q(0))
    need(error_h_square_lower == Q(join['fixed_width_error_H_norm_squared_lower']),
         'complete local H-error squared lower independently reconstructed')
    error_h_strict_lower = Q(join['fixed_width_error_H_norm_strict_lower'])
    need(error_h_strict_lower > 0 and
         error_h_strict_lower**2 < Q(join['fixed_width_error_H_norm_squared_lower']),
         'authenticated strict local full-H error lower bound')
    need(error_h_strict_lower**2 > Q(999,1000)**2*actual_h_square_upper,
         'same full-H pointwise relative error exceeds 999 per thousand')
    # Fixtures exercise the actual role/scale acceptance predicates.
    reject('old affine error interval used as actual denominator',
           lambda:endpoint_scale(list(map(Q,h['endpoint']['rows'][12]['actual_minus_smoothed_u_rr_at_fixed_width'])), 'wrong denominator role'))
    for name, wrong in (('omit full feedback inverse',J-Q(old['differentiated_initial_free_wave_X_cap'])),
                        ('ignore changed initial majorant',Q(0))):
        reject(name,lambda wrong=wrong:need(wrong == delta,'require exact feedback delta'))
    reject('omit paired initial comparison',lambda:need(u2 == smooth_h,'full initial comparison required'))
    reject('omit smoothing charge',lambda:need(u2 == paired_h,'full convolution charge required'))
    if a.normalization_report is not None:
        need(a.normalization_sha256 is not None,'normalization report has explicit digest')
        norm = load(a.normalization_report,'normalization',a.normalization_sha256)
        allowed = {'global_coordinate_rr_sup_ratio':True,
            'global_coordinate_maximum_second_derivative_seminorm_ratio':True,
            'uniform_pointwise_relative_error_upper_bound':False,
            'one_initial_point_same_H_relative_error_lower_bound':True,
            'one_initial_component_relative_error_lower_bound':True,
            'separate_relative_xx_or_xr_sup_ratios':False,
            'full_C2_norm_relative_error':False,'physical_covariant_or_observational_relative_accuracy':False,
            'finite_amplitude_nonlinear_or_cosmological_realization':False}
        required = {'same_actual_target_initial_trace_Urr_equals_U0rr':True,
            'classical_physical_side_C2_with_one_sided_corner_trace':True,
            'registered_coordinate_and_positive_full23_H_norm':True}

        def validate_normalization(d):
            need(d['status'].startswith('PASS'),'successful normalization producer report')
            need(d['normalization_quantity'] == 'global_rr_sup_and_max_second_derivative_seminorm',
                 'normalization quantity is the explicit global coordinate seminorm')
            need(d['claimed_conclusions'] == allowed,'only justified coordinate conclusions are claimed')
            need(d['required_premises'] == required,'required classical target trace and norm premises explicit')
            need(Q(d['target_scale_strict_lower']) == scale,'same actual target denominator scale')
            need(d['candidate_changed'] is False and d['actual_sources_changed'] is False and
                 d['target'] == old['target'] and d['domain'] == old['domain'],
                 'normalization retains the same conditional target and domain')
            need({k:Q(v) for k,v in d['global_coordinate_ratios'].items()} == ratios,
                 'all producer global coordinate ratios match independent arithmetic')
            need({k:Q(v) for k,v in d['baseline_PR103_global_coordinate_ratios'].items()} == baseline,
                 'all producer baseline coordinate ratios match independent arithmetic')
            need(d['input_sha256']['combined'] == a.initial_sha256 and
                 d['input_sha256']['c2'] == a.transfer_sha256,'normalization retains exact initial and transfer input pins')
            for k,v in ratios.items():
                need(v < Q(d['global_coordinate_ratio_readable_strict_uppers'][k]),
                     'safe readable coordinate ratio '+k)

        validate_normalization(norm)
        calibration = norm['local_component_calibration']
        need(list(map(Q,calibration['actual_target_interval'])) == local_target and
             list(map(Q,calibration['actual_minus_candidate_interval'])) == local_error,
             'producer local component calibration retains authenticated actual intervals')
        need(Q(calibration['component_relative_absolute_error_lower_bound']) == local_lower and
             Q(calibration['component_relative_absolute_error_readable_strict_lower']) == Q(999,1000),
             'producer local component quotient lower bound independently reproduced')
        local_h = norm['local_same_H_calibration']
        need(Q(local_h['actual_H_squared_upper_from_all23_intervals']) == actual_h_square_upper and
             Q(local_h['error_H_squared_lower_from_all23_intervals']) == error_h_square_lower,
             'producer local H actual/error squared bounds independently reproduced')
        h_actual_upper = Q(local_h['actual_target_H_norm_upper'])
        need(h_actual_upper > 0 and h_actual_upper**2 >= actual_h_square_upper and
             Q(local_h['actual_error_H_norm_strict_lower']) == error_h_strict_lower,
             'producer local H actual upper and strict error lower valid')
        need(Q(local_h['pointwise_H_relative_error_strict_lower']) == error_h_strict_lower/h_actual_upper > Q(999,1000) and
             Q(local_h['pointwise_H_relative_error_readable_strict_lower']) == Q(999,1000),
             'producer same-H pointwise relative error lower independently reproduced')
        mutations = (
            ('NOT_PASS normalization status',lambda d:d.update(status='NOT_PASS')),
            ('uniform pointwise interpretation',lambda d:d['claimed_conclusions'].update(uniform_pointwise_relative_error_upper_bound=True)),
            ('separate xx or xr denominator',lambda d:d['claimed_conclusions'].update(separate_relative_xx_or_xr_sup_ratios=True)),
            ('full C2 norm relabeling',lambda d:d['claimed_conclusions'].update(full_C2_norm_relative_error=True)),
            ('physical accuracy relabeling',lambda d:d['claimed_conclusions'].update(physical_covariant_or_observational_relative_accuracy=True)),
            ('missing classical initial trace premise',lambda d:d['required_premises'].update(same_actual_target_initial_trace_Urr_equals_U0rr=False)),
            ('wrong denominator scale',lambda d:d.update(target_scale_strict_lower='1.6245e33')),
        )
        for label,mutate in mutations:
            bad = copy.deepcopy(norm)
            mutate(bad)
            reject(label,lambda bad=bad:validate_normalization(bad))
    for role,(path,raw) in files.items():
        need(path.read_bytes() == raw,'readback '+role)
    result = {'status':'PASS_INDEPENDENT_COMBINED_TRANSFER_COORDINATE_NORMALIZATION_AUDIT',
        'input_sha256':hashes,'combined_initial_u_rr_H_cap':u2,'independent_initial_acceleration_H_cap':acceleration,
        'independent_differentiated_initial_free_wave_X_cap':J,'each_C2_upper_cap_delta_from_PR103':delta,
        'original_target_second_coordinate_error_H_caps':evolved,
        'actual_endpoint_interval':actual,'independent_actual_endpoint_interval':independent,
        'strict_actual_target_global_coordinate_H_sup_lower':scale,
        'global_coordinate_ratios':ratios,'baseline_PR103_global_coordinate_ratios':baseline,
        'local_component_calibration':{'initial_point':['0','65/64'],'component':'w:0',
            'actual_target_u0_rr_interval':local_target,'actual_minus_candidate_u0_rr_interval':local_error,
            'component_relative_error_lower_bound':local_lower,
            'readable_component_relative_error_strict_lower':'999/1000'},
        'local_same_H_calibration':{'initial_point':['0','65/64'],
            'actual_target_H_norm_squared_upper':actual_h_square_upper,
            'actual_error_H_norm_strict_lower':error_h_strict_lower,
            'pointwise_H_relative_error_squared_strict_lower':error_h_strict_lower**2/actual_h_square_upper,
            'pointwise_H_relative_error_readable_strict_lower':'999/1000'},
        'target_trace_premise':'Retained physical-side C2 trace and zero-Cauchy integral lifts imply U_rr(0,7/4)=actual U0_rr(7/4)',
        'coordinate_semantics':{'global_coordinate_ratio':True,'uniform_pointwise_relative_error_upper_bound':False,
            'one_initial_point_same_H_relative_error_lower_bound':True,
            'separate_xx_xr_relative_error':False,'full_C2_norm':False,'physical_accuracy':False},
        'named_checks':checks,'wrong_formula_controls_rejected':controls,
        'scope':'Conditional global ratios in the fixed scaled-coordinate H rr supremum and maximum-second-derivative seminorm only. A large initial endpoint can dominate the denominator. No pointwise, separate xx/xr, full C2 norm, covariant or physical accuracy; no nonlinear or cosmological realization.'}
    print(json.dumps(encode(result),indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
