#!/usr/bin/env python3
"""Independent exact arithmetic replay of the anchored-volume join.

No producer/join module is imported. The first-order error bound is derived
directly in dimensionless ordinary time t=s/rho, rather than copied in x.
All sign decisions and comparisons use Fraction. Floats are displays only.
"""
from fractions import Fraction as Q
import argparse
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
JOIN_REPORT_PIN = 'f37b2ce9ebb45ae96e86ed715180f752266b2013b42d59a37265c71196ecf77f'
JOIN_SOURCE_PIN = '589ee7a7efab0acb33030562c554f315f7217a580d53a899932ea182cc3fc2ef'
PINS = {
    'numeric/COMPONENT_LEDGER_TEN_INTEGRATED_NORMAL.json': '12dfb1247dba4ca1a9ed59f49e3c67209ff8c839d0d554180b5382d4b126664d',
    'numeric/SECOND_FEEDBACK_NORMAL.json': '23883937b7e5a55f1155e8b5c795644253ed5cf0abd51e0c40a09ece50220974',
    'balanced/b2/EXACT_B2_NORMAL.json': '3db587c30dca20f01e4964ac87f8f5ba457d55d15f7720d09eb840ed3e0e4935',
    'repo/research/HDBLAST_COMMON_GAUGE_FORCING_20261016/shell-producer/COMMON_HARMONIC_SHELL_NORMAL.json': '559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a',
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def need(p, msg):
    if not p:
        raise AssertionError(msg)


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, Path):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def stats(iv):
    a, b = iv
    return {'exact': iv, 'integer_outward': [a.__floor__(), b.__ceil__()],
            'display_only': [float(a), float(b)], 'strictly_negative': b < 0,
            'strictly_positive': a > 0}


def inflate(iv, r):
    a, b = map(Q, iv)
    need(a <= b and r >= 0, 'ordered interval / nonnegative radius')
    return [a-r, b+r]


def first_bound(caps, ell, x, delta, derivative_multiplier=Q(1)):
    # Set rho=1 in the dimensionless ordinary-time formula. Every rho
    # cancels. C=integral X dt; |H*rho|<=t; rho/A<=1.
    get = lambda k: Q(caps[k]['strict_readable_cap']['exact'])
    X, Z, V, d = map(get, ['X', 'Z', 'V_gauge', 'd'])
    Zt, dt = [get(k)*derivative_multiplier for k in ['Z_t', 'd_t']]
    t, active_length = ell*x, ell*delta
    C = active_length*X
    E = Z+t*C
    shift_in_clock = V+C
    groups = {
        'first_clock_history': t*(active_length*X*X/Q(2)
                                +4*active_length*shift_in_clock**2+C*X),
        'first_evaluation': C*(3*Zt+6*t*Z+8*dt+16*t*d)/3,
        'background_evaluation': C*C*(1+t*t)/2,
        'shift': Q(8, 3)*C*V,
        'clock_gradient': Q(4, 3)*C*C,
        'spatial_determinant': E*E/2+Q(8, 3)*E*d+Q(16, 3)*d*d,
    }
    return groups, sum(groups.values(), Q(0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--capsule-root', type=Path, default=BASE)
    parser.add_argument('--report', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    root, repo = args.capsule_root.resolve(), args.repo_root.resolve()
    report_path = args.report.resolve() if args.report is not None else root/'results/JOINED_NORMAL.json'
    output = args.output.resolve()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT_REQUIRED')
    need(root not in output.parents, 'OUTPUT_INSIDE_CAPSULE')
    need(repo not in output.parents, 'OUTPUT_INSIDE_REPO')

    def member(logical):
        return repo/logical[5:] if logical.startswith('repo/') else root/logical

    hashes = {}
    data = {}
    for rel, pin in PINS.items():
        p = member(rel)
        hashes[rel] = sha(p)
        need(hashes[rel] == pin, 'input pin: '+rel)
        data[rel] = json.loads(p.read_text())
    component, error, candidate, first = data.values()
    report = json.loads(report_path.read_text())
    logical_sources = [
        'repo/research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/observable/SHELL_PROPER_TIME_OBSERVABLE.md',
        'analysis/join_results.py', 'analysis/FINITE_TUBE_OBSERVABLE.md',
        'numeric/component_ledger.py', 'numeric/check_second_feedback.py',
        'balanced/b2/exact_b2.py', 'balanced/b2/check_b2.py',
    ]
    source_paths = {logical:member(logical) for logical in logical_sources}
    source_paths['results/JOINED_NORMAL.json'] = report_path
    source_paths['review-volume/check_volume_join.py'] = Path(__file__).resolve()
    for logical, p in source_paths.items():
        hashes[logical] = sha(p)
    need(sha(report_path) == JOIN_REPORT_PIN, 'final join report pin')
    need(sha(root/'analysis/join_results.py') == JOIN_SOURCE_PIN, 'final join source pin')
    need(report['join_source_sha256'] == sha(root/'analysis/join_results.py'), 'join source identity')
    for rel, pin in PINS.items():
        report_rel = rel[5:] if rel.startswith('repo/') else rel
        need(report['input_sha256'][report_rel] == pin, 'report input pin: '+rel)
    need(component['authentication'] == error['authentication'] == candidate['authentication'], 'same seven frozen closures')
    need(all(d['complete_closures_authenticated_before_and_after'] is True for d in [component,error,candidate]), 'source authentication')
    order = ['K0','L0','N0','J0','w0']
    need(candidate['mean_order'] == component['mean_order'] == order, 'mean order')
    need(error['no_old_u_du_p1_defect'] is True, 'exact-u error scope')
    need(candidate['arithmetic']['fixed_subintervals'] == 256, 'frozen budget')
    samples = []
    for x in [Q(9,8),Q(5,4)]:
        cv = next(r for r in candidate['samples'] if Q(r['x']) == x)
        er = next(r for r in error['results'] if Q(r['time']) == x)
        actual = next(r for r in report['shell_samples'] if Q(r['x']) == x)
        need(er['state_order'] == order+order, 'error state order')
        need(list(map(Q,cv['evaluation']['three_term'][1])) == [0,0], 'Dirichlet shell candidate')
        means = {}
        for i,name in enumerate(order):
            cand = list(map(Q,cv['evaluation']['three_term'][i]))
            radius = Q(er['complete_shell_value_error_upper'][i])
            value = inflate(cand,radius)
            need(value == list(map(Q,actual['means'][name]['complete_actual_response_interval'])), 'actual interval '+name)
            need(radius == Q(actual['means'][name]['actual_response_error_upper']), 'radius '+name)
            need(cand == list(map(Q,actual['means'][name]['candidate_interval_including_quadrature'])), 'candidate once '+name)
            quad = Q(cv['quadrature_error_component_uppers']['three_term'][i])
            means[name] = {'candidate_including_quadrature':stats(cand),
                           'actual_response_error_upper':radius,
                           'actual_interval':stats(value),
                           'quadrature_upper_for_wrong_double_charge_control':quad,
                           'wrong_double_charged_interval':stats(inflate(value,quad))}
        old = next(r for r in component['shell_sample_intervals'] if Q(r['x']) == x)
        oldn = list(map(Q,old['traces']['N0']['conditional_response_value_interval']))
        newn = means['N0']['actual_interval']['exact']
        intersection = [max(oldn[0],newn[0]), min(oldn[1],newn[1])]
        need(intersection[0] <= intersection[1] < 0, 'compatible negative N intervals')
        need(intersection == list(map(Q,actual['N_interval_intersecting_two_valid_bounds'])), 'N intersection equality')
        need(means['J0']['actual_interval']['strictly_negative'], 'J negative')
        generic_radius = Q(er['generic_value_derivative_error'][3])
        generic_j = inflate(means['J0']['candidate_including_quadrature']['exact'],generic_radius)
        need(generic_j == list(map(Q,actual['coarser_generic_J_interval'])), 'generic J interval equality')
        need(generic_j[1] < 0, 'generic J negative')
        need(stats(generic_j)['integer_outward'] == actual['coarser_generic_J_integer_outward_interval'], 'generic J outward equality')
        samples.append({'x':x,'means':means,'N_intersection':stats(intersection),
                        'coarser_generic_J_error_upper':generic_radius,
                        'coarser_generic_J_interval':stats(generic_j)})

    ell, x, delta = Q(1,10**6), Q(9,8), Q(1,8)
    caps = first['bounds']['actual_field_amplitude_caps']
    for k in ['X','Z','V_gauge','d','Z_t','d_t']:
        need(Q(caps[k]['exact_upper']) < Q(caps[k]['strict_readable_cap']['exact']), 'strict first cap '+k)
    groups, qfirst = first_bound(caps,ell,x,delta)
    joined = report['anchored_volume_second_jet_at_x_9_8']
    need(groups == {k:Q(v) for k,v in joined['first_quadratic_group_uppers'].items()}, 'all six Q groups')
    need(qfirst == Q(joined['first_quadratic_total_upper']), 'Q total')
    need(qfirst < Q(55,10**16), 'readable Q bound')
    ten = next(r for r in component['value_derivative_ten_state_comparison'] if Q(r['time']) == x)
    need(ten['value_then_derivative_order'] == order+order, 'value component selected')
    need(ten['exact_positive_inverse_and_identity_verified'] is True, 'positive history comparison')
    ncap = Q(ten['full_response_value_derivative_upper'][2])
    need(ncap == Q(joined['direct_N_value_upper_on_history']) and 0 <= ncap < Q(1213)*10**19, 'direct N cap')
    # Ordinary time t makes both powers of ell explicit: H*rho<=ell*x,
    # and active t-history has length ell*(x-1).
    nclock = Q(3,2)*(ell*x)*(ell*delta)*ncap
    need(nclock == Q(joined['entire_N_clock_term_upper']) and nclock < 256*10**7, 'full N clock cap')
    j = samples[0]['means']['J0']['actual_interval']['exact']
    value = [Q(3,2)*j[0]-nclock-qfirst,Q(3,2)*j[1]+nclock+qfirst]
    need(value == list(map(Q,joined['v2_interval'])), 'volume exact interval')
    need(value[1] < 0, 'volume negative')
    need([v/3 for v in value] == list(map(Q,joined['volume_radius_second_jet_R2_over_R0_interval'])), 'radius derivative convention')
    generic_j = samples[0]['coarser_generic_J_interval']['exact']
    generic_volume = [Q(3,2)*generic_j[0]-nclock-qfirst,Q(3,2)*generic_j[1]+nclock+qfirst]
    need(generic_volume[1] < 0, 'coarser generic volume also negative')
    wrong_derivative_groups, wrong_q = first_bound(caps,ell,x,delta,1/ell)
    qj = samples[0]['means']['J0']['quadrature_upper_for_wrong_double_charge_control']
    controls = {
        'omit_first_time_chain_rule_ell': {'description':'Using ordinary-time derivative caps as x derivative caps multiplies evaluation derivative terms by 1/ell.',
                                         'wrong_Q_upper':wrong_q,'correct_Q_upper':qfirst,'exact_bound_changes':wrong_q != qfirst},
        'double_charge_candidate_quadrature': {'description':'The B2 evaluation already contains quadrature radii; adding the metadata cap again widens J and V spuriously.',
                                              'wrong_v2_interval':stats([value[0]-Q(3,2)*qj,value[1]+Q(3,2)*qj]),
                                              'spurious_volume_radius':Q(3,2)*qj},
        'omit_trace_factor_three': {'description':'Replacing tr(n_ij)=3J by J gives J/2 and a different interval despite retaining the same sign.',
                                   'wrong_v2_interval':stats([j[0]/2-nclock-qfirst,j[1]/2+nclock+qfirst])},
        'omit_both_clock_ell_factors': {'description':'The wrong bound (3/2)*x*delta*Ncap loses the two ordinary-time scale factors and cannot certify this sign.',
                                       'wrong_clock_bound':nclock/ell**2,
                                       'wrong_v2_interval':stats([Q(3,2)*j[0]-nclock/ell**2-qfirst,Q(3,2)*j[1]+nclock/ell**2+qfirst])},
    }
    need(controls['omit_first_time_chain_rule_ell']['exact_bound_changes'], 'chain rule control differs')
    need(not controls['omit_both_clock_ell_factors']['wrong_v2_interval']['strictly_negative'], 'clock scale control loses sign')
    result = {
        'schema':'HDBLAST_INDEPENDENT_ANCHORED_VOLUME_JOIN_REVIEW_V1',
        'status':'PASS_INDEPENDENT_EXACT_INTERVAL_REPLAY_AND_FROZEN_FORMULA_REVIEW',
        'input_and_source_sha256':hashes,
        'method':'No join/producer modules imported. Fraction arithmetic; Q first groups derived in ordinary dimensionless t=s/rho with active length ell/8.',
        'samples':samples,
        'anchored_volume_at_x_9_8':{'ell':ell,'x':x,'active_x_length':delta,'direct_N_history_value_cap':ncap,
                                   'entire_N_clock_upper':nclock,'Q_first_groups':groups,'Q_first_upper':qfirst,
                                   'v2_interval':stats(value),'R2_over_R0_interval':stats([v/3 for v in value]),
                                   'coarser_generic_v2_interval':stats(generic_volume),
                                   'negative_margin_to_zero':-value[1]},
        'wrong_formula_controls':controls,
        'manual_source_findings':[
            'Frozen equations (6),(8): Bclock=1/2 integral N ds + first-order quadratic terms; tr(n_spatial)=3J gives +(3/2)J. Thus both J and lapse-history signs and factors match.',
            's=rho ell x gives H integral N ds = ell tanh(ell x) integral N dx; rho cancels. Exact history vanishes for x<=1.',
            'The N bound is a direct value cap on the shell-history causal triangle, selected from index 2 of value-then-derivative order. It is not a derivative cap multiplied by a history length.',
            'B2 source builds moment intervals center+quadrature error, then evaluation; the separate second-feedback ledger contains response defects and tail, and explicitly requires candidate quadrature separately. Joining evaluation plus response radius charges quadrature once.',
            'All six first-quadratic groups bound the complete frozen (8) after separating Bclock N history; Hprime+2H²=(1+tanh² t)/rho² and 1/cosh t<=1 suffice.',
            'A smooth family on a common open finite buffer gives a compact normal tube for some unspecified small epsilon by ODE dependence and the parameter-extended inverse function theorem; no measured caustic or epsilon bound is asserted.',
            'The smooth kinematic completion n=N0 ds²+J0 A² kappa has the required mean nss and spatial trace. Nonmean second angular tensors and nsi do not enter this integrated second variation.',
        ],
        'limits':[
            'Assumes the previously independently reviewed pinned component, second-feedback, candidate and first-trace inputs and their inherited conditional response premises.',
            'A kinematical anchored proper-time parameter second derivative only. No finite-amplitude interval, nonlinear Einstein family, force, expansion-rate sign or cosmology conclusion.',
            'Negative J and N samples alone do not establish the invariant-volume sign; the complete clock and first-order terms are bounded here.',
        ],
        'blockers':[],
    }
    with args.output.open('x') as handle:
        handle.write(json.dumps(encode(result),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),'status':result['status'],
                      'J_9_8':stats(j)['integer_outward'],'v2':stats(value)['integer_outward'],
                      'Nclock_display':float(nclock),'Qfirst_display':float(qfirst)},sort_keys=True))


if __name__ == '__main__':
    main()
