#!/usr/bin/env python3
"""Conditional natural-equation coefficient residual for unchanged W1_quadratic.

The additional first-harmonic and formal-target subsidiary premises must be
named explicitly. They are not established by numerical enclosures. This
checker imports no scientific producer code and never consumes filtered jets.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
from pathlib import Path

PINS = {
    'natural_bridge': '5cc576f03fc1d2245249d9d9b4b31a955e5d6fe817733930e80687a796904d6d',
    'time': '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
    'formulas': 'f3f037f890655ab05f936491d378a494b56e10f0c8a3d32dda0f1b8b2ea7556f',
    'source_map': '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615',
    'prior_C2': 'fa13efe588fa19ab6a09c517b67bc277cd3dcaeee9d428143bf479ef98128bb0',
}


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def load(path, digest):
    data = Path(path).read_bytes()
    need(hashlib.sha256(data).hexdigest() == digest, 'INPUT_HASH:'+str(path))
    return json.loads(data)


def poly(terms):
    result = {}
    for term in terms:
        key = tuple(sorted((name, int(power)) for name, power in term['powers'].items() if power))
        result[key] = result.get(key, Q(0))+Q(term['coefficient'])
    return {key: value for key, value in result.items() if value}


def add(*pairs):
    result = {}
    for factor, terms in pairs:
        for key, value in terms.items():
            result[key] = result.get(key, Q(0))+factor*value
    return {key: value for key, value in result.items() if value}


def shifted(terms, factor, changes):
    result = {}
    for key, value in terms.items():
        powers = dict(key)
        for name, amount in changes.items():
            powers[name] = powers.get(name, 0)+amount
        result[tuple(sorted((name, power) for name, power in powers.items() if power))] = factor*value
    return result


def sqrt_integer_upper(value):
    need(value >= 0, 'nonnegative square')
    k = isqrt(value.numerator//value.denominator)
    if Q(k*k) < value:
        k += 1
    need(Q(k*k) >= value, 'outward integer square root')
    return Q(k)


def upper4(value):
    need(value > 0, 'positive residual cap')
    scale = Q(1)
    while value >= 10000*scale:
        scale *= 10
    while value < 1000*scale:
        scale /= 10
    return (value//scale+1)*scale


def encode(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--c2-report', type=Path, required=True)
    ap.add_argument('--c2-sha256', required=True)
    ap.add_argument('--assume-first-harmonic-neighbourhood', action='store_true')
    ap.add_argument('--assume-formal-subsidiary-completion', action='store_true')
    ap.add_argument('--inject-error', choices=('symmetrization', 'metric-normalization',
                    'radial-product-sign', 'filtered-source-charge', 'value-bound'))
    a = ap.parse_args()
    checked, controls = [], []

    def ck(ok, reason):
        need(ok, reason)
        checked.append(reason)

    ck(a.assume_first_harmonic_neighbourhood,
       'additional retained first-harmonic neighbourhood premise must be explicit')
    ck(a.assume_formal_subsidiary_completion,
       'additional retained formal subsidiary completion premise must be explicit')
    root = a.repo_root
    paths = {
        'natural_bridge': root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/NATURAL.json',
        'time': root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference/TIME.json',
        'formulas': root/'research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/reference/FORMULAS.json',
        'source_map': root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json',
        'prior_C2': root/'research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/reference/C2.json',
    }
    inherited = {role: load(path, PINS[role]) for role, path in paths.items()}
    bridge, time, formulas, maps, previous = [inherited[key] for key in
        ('natural_bridge', 'time', 'formulas', 'source_map', 'prior_C2')]
    c2 = load(a.c2_report, a.c2_sha256)
    ck(c2['status'] == 'PASS_COMBINED_INITIAL_ORIGINAL_TARGET_C2_ERROR',
       'new report is a complete original-target absolute C2 error certificate')
    ck(c2['candidate'] == previous['candidate'] and c2['candidate'].startswith('W1_quadratic'),
       'unchanged original-source W1_quadratic candidate')
    ck(c2['candidate_changed'] is False and c2['actual_sources_changed'] is False and
       c2['actual_target_unchanged'] is True, 'same candidate target and actual source')
    ck(c2['target'] == previous['target'] and c2['domain'] == previous['domain'] and
       c2['norm'] == previous['norm'], 'same original target domain and coordinate norm')
    labels, gram = time['labels'], time['Gram']
    ck(c2['labels'] == labels == previous['labels'] == formulas['input_labels'] and
       c2['Gram'] == gram == previous['Gram'] == formulas['input_Gram'], 'same complete full23 H norm')
    ck(len(labels) == 23 and min(map(Q, gram)) > 0, 'positive complete angular Gram')
    ck(c2['boundary_error'] == previous['boundary_error'] == {'g': '0', 'g_x': '0', 'M_x': '0'},
       'same homogeneous physical-side error boundary')
    E1 = Q(c2['complete_original_C1_error_X_cap'])
    F0 = Q(c2['complete_original_reduced_defect_H_cap'])
    ck(E1 == Q(previous['complete_original_C1_error_X_cap']) and
       F0 == Q(previous['complete_original_reduced_defect_H_cap']) and E1 > 0 and F0 > 0,
       'retained C1 and original reduced-defect caps')
    jets = {key: Q(value) for key, value in c2['original_target_second_coordinate_error_H_caps'].items()}
    ck(set(jets) == {'xx', 'xr', 'rr'} and min(jets.values()) >= 0,
       'all three nonnegative absolute second-error caps')
    ck(all(jets[key] <= Q(previous['original_target_second_coordinate_error_H_caps'][key]) for key in jets),
       'new coordinate error caps do not exceed retained predecessor caps')

    # The quadratic natural/reduced difference is the symmetric derivative of
    # QC only under the explicitly retained C1=0 neighbourhood premise.
    identities = []
    sym_factor = Q(2) if a.inject_error == 'symmetrization' else Q(1)
    for field, rows in maps['maps']['QE'].items():
        for index, row in enumerate(rows):
            qe = poly(row)
            reduced = poly(maps['maps']['QEred'][field][index])
            gradient = poly(maps['maps']['gauge_gradient'][field][index])
            ck(add((1, qe), (-1, reduced)) == add((sym_factor, gradient)),
               'exact shared-source one-half symmetric gradient '+field+'/'+str(index))
            identities.append(field+'/'+str(index))
    ck(len(identities) == 20, 'all20 metric source identities')
    metric_factor = Q(-1) if a.inject_error == 'metric-normalization' else Q(-2)
    for field, rows in maps['maps']['QEred'].items():
        for index, row in enumerate(rows):
            ck(poly(maps['maps']['metric_forcing'][field][index]) ==
               shifted(poly(row), metric_factor, {'ell': 2, 'rho': 2}),
               'exact original metric forcing normalization '+field+'/'+str(index))
    for index, row in enumerate(maps['maps']['Qphi']):
        ck(poly(maps['maps']['scalar_forcing'][index]) ==
           shifted(poly(row), Q(1), {'ell': 2, 'rho': 2}),
           'exact original natural scalar forcing normalization '+str(index))

    # Rebuild the physical derivative charges from authenticated background
    # and exact harmonic-operator constants, rather than old absolute Z jets.
    ell = Q(1, 10**6)
    qder = Q(15, 8)
    ck(Q(formulas['sharp_scaled_derivative_cap_squared']) == Q(7, 2) and
       qder*qder > Q(7, 2) and qder == Q(bridge['joint_derivative_constant_upper']),
       'retained joint harmonic derivative constant')
    squares = {'alpha': Q(245, 2), 'eta': Q(131, 2), 'zeta': Q(68)}
    primitive = {key: sqrt_integer_upper(value) for key, value in squares.items()}
    ck(squares == {key: Q(value) for key, value in bridge['value_matrix_primitive_Frobenius_squares'].items()} and
       primitive == {key: Q(value) for key, value in bridge['value_matrix_primitive_operator_caps'].items()},
       'retained weighted primitive value-operator bounds')
    background = time['global']
    rho_lo, rho_hi = map(Q, background['rho_interval'])
    aa = max(map(abs, map(Q, background['a_interval'])))
    pp = max(map(abs, map(Q, background['p_interval'])))
    ck(rho_lo > 0 and rho_hi >= rho_lo, 'positive actual physical scale')
    alpha = ell*rho_hi*aa
    eta = Q(3, 2)*ell**2
    zeta = ell
    alpha_r = ell**2*(1+rho_hi**2*(aa**2+pp**2/3))
    eta_x, zeta_x = ell**2, Q(3, 2)*ell**3
    z0 = primitive['alpha']*alpha+primitive['eta']*eta+primitive['zeta']*zeta
    zx = primitive['eta']*eta_x+primitive['zeta']*zeta_x
    zr = primitive['alpha']*alpha_r
    connection = 2*alpha+5*ell+Q(7, 4)*(alpha+eta)
    ck({'value': z0, 'x_derivative': zx, 'r_derivative': zr} ==
       {key: Q(value) for key, value in bridge['Zmat_caps'].items()}, 'reconstruct exact value-matrix derivative constants')
    ck(connection == Q(bridge['gradient_DC_connection_scaled_cap']) and
       rho_lo == Q(bridge['background_caps']['rho_lower']), 'reconstruct exact physical connection and scale')

    # rho_r=-alpha rho gives ey DC=-(P_r+alpha P)/(ell^2 rho^2).
    # A polynomial product-rule fixture retains this sign independently.
    radial_numerator = Q(3, 7)+(Q(-2, 5)*Q(11, 13))
    if a.inject_error == 'radial-product-sign':
        radial_numerator = Q(3, 7)-(Q(-2, 5)*Q(11, 13))
    ck(radial_numerator == Q(3, 7)+(Q(-2, 5)*Q(11, 13)), 'retain radial physical product-rule sign')
    V0 = (Q(0) if a.inject_error == 'value-bound' else Q(3))*E1
    ck(V0 == 3*E1, 'value bound uses retained X definition without additional zero-wedge improvement')
    P0 = qder*E1+z0*V0
    Px = qder*max(jets['xx'], jets['xr'])+z0*E1+zx*V0
    Pr = qder*max(jets['xr'], jets['rr'])+z0*E1+zr*V0
    denominator = ell**2*rho_lo**2
    gradient = (Px+Pr+connection*P0)/denominator
    original_defect = F0 if a.inject_error != 'filtered-source-charge' else F0+Q(2)
    ck(original_defect == F0, 'original W1q reduced defect without filtered-source charge')
    einstein = original_defect/(2*denominator)+gradient
    scalar = original_defect/denominator
    ck(einstein > gradient > 0 and scalar > 0, 'positive conditional natural residual bounds')
    controls.extend(['source_symmetrization_factor', 'metric_forcing_factor',
                     'radial_product_rule_sign', 'filtered_source_charge', 'X_value_bound'])
    for role, path in paths.items():
        load(path, PINS[role])
    load(a.c2_report, a.c2_sha256)
    result = {
        'status': 'PASS_CONDITIONAL_W1Q_NATURAL_ERROR_TRANSFER',
        'candidate': c2['candidate'], 'target': c2['target'], 'domain': c2['domain'],
        'candidate_changed': False, 'actual_sources_changed': False, 'actual_target_unchanged': True,
        'additional_conditional_premises': {
            'first_harmonic_neighbourhood_assumed': True,
            'formal_subsidiary_completion_assumed': True,
            'first_harmonic_identity_on_neighbourhood': 'C1=0 as an actual field identity, not interval point selections',
            'formal_target_subsidiary_completion': 'DC[U]+QC=0 as an actual C1 covector identity on the domain',
            'physical_initial_shell_and_subsidiary_premises_independently_reproved_here': False,
        },
        'input_sha256': dict(PINS, new_C2=a.c2_sha256),
        'natural_equations': 'E_ab=Ric_ab-dphi_a*dphi_b-(2/3)Vg_ab; S=Box_g phi-Vprime(phi)',
        'parameter_convention': 'g=g0+lambda*q+(lambda^2/2)*W; output is second parameter derivative, lambda^2 coefficient is half',
        'error_orientation': 'e=W1_quadratic-U',
        'exact_bridge': 'E2[W]=R_metric[W]/(2ell^2rho^2)+sym(nabla DC[e]); S2[W]=-R_scalar[W]/(ell^2rho^2)',
        'complete_shared_source_metric_identities': identities,
        'candidate_source_is_filtered': False, 'filtered_source_charge_added': False,
        'retained_original_reduced_defect_H_cap': F0,
        'error_caps': {'X': E1, 'value': V0, 'x': E1, 'r': E1, **jets},
        'error_value_bound_origin': 'value<=3X from definition; no new zero-wedge or .75X improvement used',
        'background_caps': bridge['background_caps'],
        'joint_derivative_constant_upper': qder,
        'Zmat_caps': {'value': z0, 'x_derivative': zx, 'r_derivative': zr},
        'gradient_DC_connection_scaled_cap': connection,
        'scaled_harmonic_error_caps': {'P': P0, 'P_x': Px, 'P_r': Pr},
        'physical_denominator_lower': denominator,
        'physical_frame_gradient_DC_error_L2_cap': gradient,
        'natural_second_parameter_residual_L2_caps': {'trace_reversed_Einstein': einstein, 'scalar': scalar},
        'natural_second_parameter_residual_readable_strict_uppers': {
            'trace_reversed_Einstein': upper4(einstein), 'scalar': upper4(scalar)},
        'norm': 'uniform in (x,r), positive background-orthonormal-frame Frobenius S3 probability L2; no Lorentz invariant scalar norm',
        'named_checks': checked, 'wrong_formula_controls_rejected': controls,
        'scope': 'Conditional natural-equation second-coefficient residual bound for unchanged W1q, additional formal physical/subsidiary completion premise explicit. No independent complete constraint proof, small residual, relative physical accuracy, finite-amplitude nonlinear remainder, observable or cosmological evidence.',
    }
    print(json.dumps(encode(result), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
