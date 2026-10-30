#!/usr/bin/env python3
"""Independent expanded-ledger audit of the conditional W1q natural bound.

Reuses the authenticated older geometric identities, not producer code or
the older source-averaged candidate's numerical jets.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def load(path, digest):
    raw = Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == digest, 'authenticated input '+str(path.name))
    return json.loads(raw)


def enc(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: enc(v) for k, v in x.items()}
    if isinstance(x, list):
        return [enc(v) for v in x]
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--c2-report', type=Path, required=True)
    ap.add_argument('--c2-sha256', required=True)
    ap.add_argument('--natural-report', type=Path, required=True)
    ap.add_argument('--natural-sha256', required=True)
    a = ap.parse_args()
    c = load(a.c2_report, a.c2_sha256)
    n = load(a.natural_report, a.natural_sha256)
    t = load(a.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json',
             '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d')
    old = load(a.repo_root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/C2.json',
               'fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0')
    retained = load(a.repo_root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/NATURAL.json',
                    '5cc576f03fc1d2245249d9d9b4b31a955e5d6fe817733930e80687a796904d6d')
    checks = []

    def ck(ok, reason):
        need(ok, reason)
        checks.append(reason)

    ck(c['status'] == 'PASS_COMBINED_INITIAL_ORIGINAL_TARGET_C2_ERROR' and
       n['status'] == 'PASS_CONDITIONAL_W1Q_NATURAL_ERROR_TRANSFER', 'successful roles')
    ck(n['candidate'] == c['candidate'] == old['candidate'] and
       n['target'] == c['target'] == old['target'] and
       n['domain'] == c['domain'] == old['domain'], 'same original W1q target and domain')
    ck(n['candidate_changed'] is False and n['actual_sources_changed'] is False and
       n['actual_target_unchanged'] is True, 'fixed mathematical functions')
    ck(n['candidate_source_is_filtered'] is False and n['filtered_source_charge_added'] is False,
       'original unsmoothed source and original reduced defect')
    ck(n['input_sha256']['new_C2'] == a.c2_sha256, 'same new error input')
    p = n['additional_conditional_premises']
    ck(p['first_harmonic_neighbourhood_assumed'] is True and
       p['formal_subsidiary_completion_assumed'] is True and
       p['physical_initial_shell_and_subsidiary_premises_independently_reproved_here'] is False,
       'additional identities explicitly remain assumptions')
    ck(c['labels'] == old['labels'] == t['labels'] and c['Gram'] == old['Gram'] == t['Gram'] and
       len(c['labels']) == 23 and min(map(Q, c['Gram'])) > 0, 'complete positive input H norm')
    ck(set(n['complete_shared_source_metric_identities']) == set(retained['complete_metric_map_identities']) and
       len(n['complete_shared_source_metric_identities']) == 20, 'same retained full metric identity closure')
    E1 = Q(c['complete_original_C1_error_X_cap'])
    F0 = Q(c['complete_original_reduced_defect_H_cap'])
    ck(E1 == Q(old['complete_original_C1_error_X_cap']) and
       F0 == Q(old['complete_original_reduced_defect_H_cap']), 'unchanged first-error and defect caps')
    e = {k: Q(v) for k, v in c['original_target_second_coordinate_error_H_caps'].items()}
    ck(set(e) == {'xx', 'xr', 'rr'} and min(e.values()) >= 0, 'all three second errors')

    # Independent degree-by-degree weighted Frobenius arithmetic for the
    # retained full23 ->8 covector value operator. Covectors have no metric
    # mixed-component factor two.
    alpha2 = 3*(16+1+3)+3*Q(25, 2)+2*Q(25, 2)
    eta2 = 3*Q(9, 2)+3*(9+3)+2*8
    zeta2 = sum((lam+Q(7, 12)*lam+Q(2, 3)*(lam-3) for lam in (8, 24)), Q(0))
    ck((alpha2, eta2, zeta2) == (Q(245, 2), Q(131, 2), Q(68)), 'independent primitive covector norm squares')
    ck(12**2 > alpha2 and 9**2 > eta2 and 9**2 > zeta2, 'retained rational primitive caps')
    ell = Q(1, 10**6)
    low, high = map(Q, t['global']['rho_interval'])
    aa = max(abs(Q(v)) for v in t['global']['a_interval'])
    pp = max(abs(Q(v)) for v in t['global']['p_interval'])
    ck(low > 0 and high >= low, 'positive actual frame scale')
    alpha, eta = ell*high*aa, Q(3, 2)*ell**2
    alpha_r = ell**2*(1+high**2*(aa**2+pp**2/3))
    z0 = 12*alpha+9*eta+9*ell
    zx, zr = 9*ell**2+Q(27, 2)*ell**3, 12*alpha_r
    connection = Q(15, 4)*alpha+5*ell+Q(7, 4)*eta
    ck(n['Zmat_caps'] == {'value': str(z0), 'x_derivative': str(zx), 'r_derivative': str(zr)},
       'independent value-operator derivatives')
    ck(Q(n['gradient_DC_connection_scaled_cap']) == connection, 'independent combined frame connection charge')
    q = Q(15, 8)
    ck(q*q > Q(7, 2) and Q(n['joint_derivative_constant_upper']) == q,
       'joint derivative map rational constant')
    V0 = 3*E1
    ck({k: Q(v) for k, v in n['error_caps'].items()} ==
       dict(X=E1, value=V0, x=E1, r=E1, **e), 'X value and complete error jets retained')
    denominator = ell**2*low**2
    # Expand the sum directly rather than relying on the producer's P/Px/Pr
    # intermediates. The radial rho product derivative and all connection
    # charges are retained in the authenticated geometric inequality.
    principal = q*(max(e['xx'], e['xr'])+max(e['xr'], e['rr']))
    lower = 2*z0*E1+(zx+zr)*V0+connection*(q*E1+z0*V0)
    gradient = (principal+lower)/denominator
    metric = (F0/2+principal+lower)/denominator
    scalar = F0/denominator
    ck(Q(n['physical_denominator_lower']) == denominator, 'ordinary physical second-derivative conversion')
    ck(Q(n['physical_frame_gradient_DC_error_L2_cap']) == gradient, 'expanded harmonic-error gradient ledger')
    ck(Q(n['retained_original_reduced_defect_H_cap']) == F0, 'original defect without averaging charge')
    expected = {'trace_reversed_Einstein': metric, 'scalar': scalar}
    ck({k: Q(v) for k, v in n['natural_second_parameter_residual_L2_caps'].items()} == expected,
       'independent natural coefficient bounds')
    for k, value in expected.items():
        ck(0 < value < Q(n['natural_second_parameter_residual_readable_strict_uppers'][k]),
           'strict outward readable bound '+k)
    controls = []
    for name, wrong, correct in (
        ('missing frame scale', F0/2+principal+lower, metric),
        ('metric scalar factor confused', F0/denominator+gradient, metric),
        ('coefficient versus second parameter derivative', metric/2, metric),
        ('omit first-derivative and connection terms', principal/denominator, gradient),
        ('omit radial rho product charge', gradient-alpha*(q*E1+z0*V0)/denominator, gradient),
    ):
        ck(wrong != correct, 'reject '+name)
        controls.append(name)
    result = {'status': 'PASS_INDEPENDENT_CONDITIONAL_NATURAL_ERROR_AUDIT',
        'input_sha256': {'new_C2': a.c2_sha256, 'natural': a.natural_sha256},
        'expanded_principal_gradient_numerator': principal, 'expanded_lower_order_gradient_numerator': lower,
        'natural_second_parameter_residual_L2_caps': expected,
        'additional_target_subsidiary_identity_is_assumed': True,
        'named_checks': checks, 'wrong_formula_controls_rejected': controls,
        'scope': 'Independent arithmetic audit of a conditional second-parameter natural equation residual in a fixed positive physical-frame comparison norm. No independent subsidiary proof or physical relative accuracy.'}
    print(json.dumps(enc(result), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
