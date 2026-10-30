#!/usr/bin/env python3
"""Conditional original-target C2 gate, with actual unfiltered F_x input.

No scientific imports.  This is an error bound for U-W1_quadratic, not for
the causally averaged auxiliary candidate introduced in PR101.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


PINS = {
    'higher_initial': 'bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207',
    'initial_repair': '4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002',
    'quadratic_transfer': '4d0b83994270f6e2d3aa5433887a0d4b978acbac372bf28197050345de2881c3',
    'actual_time_coefficients': '6bf1b2fbb2b6f6a9df2f2a24bcac142bfe881437fa5638bcff20b20c7b3b509d',
    'refined_response': 'ef097736c337003f0e2d578a830d14824772065cc36f950932819f524c6fbe84',
    'prior_candidate_ledger': '44d05eadf32770e152ba15e8482f2d290a8baafdcf59f20ac281841960b76c1f',
}
checks = []
controls = []


def need(ok, label):
    if not ok:
        raise RuntimeError(label)
    checks.append(label)


def load(path, role, digest=None):
    raw = Path(path).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    need(actual == (digest or PINS[role]), 'authenticate ' + role)
    return json.loads(raw)


def reject(name, wrong, correct):
    need(wrong != correct, 'reject ' + name)
    controls.append(name)


def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def readable_strict_upper(value):
    need(value > 0, 'positive coordinate error upper bound')
    scale = Q(1)
    while value >= 10000*scale:
        scale *= 10
    while value < 1000*scale:
        scale /= 10
    return (value//scale+1)*scale


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--higher-capsule', type=Path)
ap.add_argument('--causal-capsule', type=Path)
ap.add_argument('--refined-response', type=Path)
ap.add_argument('--fx-report', type=Path)
ap.add_argument('--fx-sha256')
args = ap.parse_args()
root = args.repo_root
hi = args.higher_capsule or root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
ca = args.causal_capsule or root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030'
initial = load(hi/'reference/INITIAL.json', 'higher_initial')
repair = load(hi/'reference/REPAIR.json', 'initial_repair')
t = load(hi/'reference/TRANSFER.json', 'quadratic_transfer')
dt = load(hi/'reference/TIME.json', 'actual_time_coefficients')
r = load(args.refined_response or root/'research/HDBLAST_INITIAL_CUTOFF_ENCLOSURES_20261030/reference/RESPONSE.json', 'refined_response')
prior = load(ca/'reference/CANDIDATE.json', 'prior_candidate_ledger')
need(t['status'] == r['status'] == 'PASS', 'successful retained response inputs')
need(t['labels'] == initial['labels'] == repair['labels'] == dt['labels'] and
     t['Gram'] == initial['gram'] == repair['gram'] == dt['Gram'], 'same all23 positive norm')
need(len(t['labels']) == 23 and min(map(Q, t['Gram'])) > 0, 'complete positive full23 H')
need(t['direct_reproof']['extension_component_curvature_caps_unchanged'] and
     repair['velocity_and_its_initial_error_bound_unchanged'], 'same PR100 candidate initial data')

T = Q(3, 2)
op = {k: Q(v) for k, v in r['actual_coefficient_H_norm_caps'].items()}
frozen = {k: Q(v) for k, v in t['frozen_H_operator_caps'].items()}
delta = {k: Q(v) for k, v in r['actual_coefficient_error_H_norm_caps'].items()}
der = {k: Q(v) for k, v in dt['global']['H_operator_caps'].items()}
need(all(op[k] == frozen[k] + delta[k] for k in op), 'actual and frozen coefficient discrepancy caps')
need(der['B_x'] == der['M_x'] == 0, 'constant actual B and Robin M in time')
c = op['A'] + op['B'] + Q(3, 4)*op['D']
comm = der['A_x'] + der['B_x'] + Q(3, 4)*der['D_x']
b = op['M']/4
q = Q(r['retained_coupled_X_operator_cap'])
P = 1/(1-b)
need(T*c == Q(t['k_actual_bulk_X_operator_cap']) and b == Q(t['b_actual_boundary_X_operator_cap']) and
     P == Q(t['P_boundary_inverse_X_cap']), 'same actual physical-side Robin inverse')
need(0 < b < 1 and T*c+b <= q < 1 and q == Q(3409, 8000000), 'retained full causal contraction')
need(comm == Q(dt['coefficient_commutator_X_to_H_upper']), 'complete actual coefficient time derivative ledger')

data = {k: Q(v) for k, v in prior['initial_absolute_H_caps'].items()}
Bs = Q(t['direct_reproof']['new_smooth_Bs_X_cap'])
E0 = Q(t['direct_reproof']['residual_majorants']['W_minus_B0_whole_cone_X_cap'])
C1 = Q(t['direct_reproof']['new_actual_bulk_correction_X_cap'])
E = Q(t['best_reproved_global_C1_error_X_cap'])
F0 = Q(t['direct_reproof']['new_reduced_forcing_defect_H_cap'])
FH = Q(prior['source_caps']['F_H'])
Gx = Q(prior['source_caps']['G_x_H'])
need(FH == 17387513156573 and Gx == 58807538, 'genuine actual source and shell derivative caps')
need(Bs == Q(prior['free_wave_X_cap']) and Bs+E0 == Q(prior['W0_X_cap']), 'same W0_quadratic and free-wave ledger')
need(F0 == c*C1, 'same original unfiltered reduced-defect magnitude')

# E0field=W0-Bs has zero Cauchy data. Its x derivative has zero
# displacement and initial velocity f0(0), plus boundary source M Bs_x+G_x.
Qfree = data['u_rr'] + data['v_r']
f00 = frozen['A']*data['v'] + frozen['B']*data['u_r'] + frozen['D']*data['u']
f0x = (frozen['A']+frozen['B'])*Qfree + frozen['D']*Bs
need(f00 == Q(prior['W0_initial_forcing_H_cap']) and
     f0x == Q(prior['frozen_forcing_x_H_cap']), 'same frozen forcing and its true initial acceleration')
E0x = P*(f00 + T*f0x + op['M']*Bs + Gx)

# r0=La(W0-Bs)+(La-L0)Bs+F. Differentiate this cancellation-
# preserving form, not the sum of two large absolute W0/Bs second jets.
r0x_const = ((op['A']+op['B'])*E0x + (op['D']+comm)*E0 +
             (delta['A']+delta['B'])*Qfree + (delta['D']+comm)*Bs)
r00 = delta['A']*data['v'] + delta['B']*data['u_r'] + delta['D']*data['u'] + FH
# C=W1-W0=P K r0 has homogeneous boundary and zero Cauchy data.
# C_x has initial velocity r0(0), not zero.  Its free-wave X cap is r00.
Cx_const = P*(r00+T*r0x_const)
Cx_slope = P*T
F1_const = (op['A']+op['B'])*Cx_const + (op['D']+comm)*C1
F1_slope = (op['A']+op['B'])*Cx_slope
need(F1_slope > 0, 'original actual F_x enters the differentiated defect')

# The same PR100 global C1 initial error bounds were reconstructed there.
eu, eur, ev = map(Q, r['initial_smoothed_data_H_error_caps'])
need(max((eu+T*ev)/3, eur+ev) == Q(t['direct_reproof']['new_initial_free_wave_X_error_cap']),
     'same individual PR100 smoothed initial error caps')
eurr = Q(repair['new_uniform_initial_u_rr_H_upper'])
evr = Q(initial['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
need(eurr < Q('6.219e33') and evr < Q('9.672e24'), 'continuous higher initial errors for original candidate')
w0 = eurr + op['A']*ev + op['B']*eur + op['D']*eu
J = max((ev+T*w0)/3, evr+w0)
Z_const = (J+T*(F1_const+comm*E))/(1-q)
Z_slope = T*F1_slope/(1-q)
rr_const = Z_const+c*E+F0
need(Z_slope == T*T*(op['A']+op['B'])/((1-b)*(1-q)), 'exact affine original-source derivative sensitivity')
need(J == evr+w0, 'derivative branch controls the actual initial free-wave majorant')

# General nonzero-initial-forcing fixture: w=x^2/2 has principal forcing1,
# zero Cauchy value/velocity, but (w_x)_x(0)=1.  Both differentiated
# finite corrections must retain their forcing initial velocity.
reject('zero initial frozen correction acceleration', Q(0), Q(1))
reject('zero initial actual bulk correction acceleration', Q(0), Q(1))
reject('omit actual Robin inverse', f00+T*f0x+op['M']*Bs+Gx, E0x)
reject('omit complete actual coefficient time commutator',
       (op['A']+op['B'])*Cx_const+op['D']*C1, F1_const)
reject('omit original source derivative sensitivity', Q(0), Z_slope)
reject('omit initial higher derivative error', w0-eurr, w0)
reject('omit full feedback inverse', J+T*(F1_const+comm*E), Z_const)

report = {
    'status': 'PASS_AFFINE_ORIGINAL_TARGET_C2_GATE_FX_INPUT_OPEN',
    'candidate': 'W1_quadratic (W1_dagger), the unchanged PR100 original-source finite candidate',
    'target': 'same original actual U with full unsmoothed F, actual coefficients and actual shell G/M',
    'domain': 'D+:0<=x<=3/2,max(0,x-5/4)<=r<=7/4-x',
    'norm': 'positive full23 H; scaled coordinate derivatives, one-sided at physical shell',
    'input_sha256': PINS,
    'initial_error_H_caps': {'e0': eu, 'e0_r': eur, 'v0': ev, 'e0_rr': eurr, 'v0_r': evr},
    'initial_original_defect': 'r1(0,r)=0, because W1_quadratic-W0_quadratic has zero Cauchy first jets',
    'initial_differentiated_acceleration_H_cap': w0,
    'differentiated_initial_free_wave_X_cap': J,
    'boundary_error': {'g': '0', 'g_x': '0', 'M_x': '0'},
    'complete_original_C1_error_X_cap': E,
    'complete_original_reduced_defect_H_cap': F0,
    'W0_minus_Bs_X_cap': E0,
    'W0_minus_Bs_time_derivative_X_cap': E0x,
    'r0_initial_H_cap': r00,
    'r0_time_derivative_H_affine_in_Fx': {'constant': r0x_const, 'Fx_coefficient': Q(1)},
    'bulk_correction_time_derivative_X_affine_in_Fx': {'constant': Cx_const, 'Fx_coefficient': Cx_slope},
    'reduced_defect_time_derivative_H_affine_in_Fx': {'constant': F1_const, 'Fx_coefficient': F1_slope},
    'error_xx_xr_H_affine_in_Fx': {'constant': Z_const, 'Fx_coefficient': Z_slope},
    'error_rr_H_affine_in_Fx': {'constant': rr_const, 'Fx_coefficient': Z_slope},
    'Fx_threshold_for_all_second_errors_below_6_222e33': (Q('6.222e33')-rr_const)/Z_slope,
    'new_Fx_input_supplied': False,
    'scope': 'A sufficient conditional original-target scaled-coordinate C2 error gate; until actual full F_x is supplied the affine output is not a numerical C2 error certificate. No filtered source/candidate substitution, physical relative accuracy, nonlinear realization, full constraint or cosmological conclusion.',
    'named_checks': checks,
    'wrong_formula_controls_rejected': controls,
}
if args.fx_report is not None:
    need(args.fx_sha256 is not None, 'explicit actual original source derivative report digest')
    fx = load(args.fx_report, 'actual_original_bulk_Fx', args.fx_sha256)
    need(fx['status'].startswith('PASS') and fx['scope']['actual_original_full23_F_x_bound'] is True and
         fx['source_coordinate'] == 'x=t/ell', 'complete original scaled-time source derivative scope')
    Fxc = Q(fx['genuine_scaled_F_x_H_cap'])
    need(Fxc >= 0, 'nonnegative actual raw F_x H cap')
    caps = {'xx': Z_const+Z_slope*Fxc, 'xr': Z_const+Z_slope*Fxc,
            'rr': rr_const+Z_slope*Fxc}
    report.update({'status': 'PASS_CONDITIONAL_ORIGINAL_TARGET_C2_ERROR', 'new_Fx_input_supplied': True,
        'actual_original_bulk_Fx_report_sha256': args.fx_sha256, 'actual_original_bulk_Fx_H_cap': Fxc,
        'original_target_second_coordinate_error_H_caps': caps,
        'original_reduced_defect_time_derivative_H_cap': F1_const+F1_slope*Fxc,
        'scope': 'Conditional continuous original-target C2 coordinate error bounds for unchanged PR100 W1_quadratic with full original sources. Absolute bounds only; no relative physical accuracy, across-doubled-shell C2, full constraint propagation, nonlinear realization or cosmological evidence.'})
    report['original_target_second_coordinate_error_H_readable_strict_uppers'] = {
        k: readable_strict_upper(v) for k, v in caps.items()}
    initial_fraction = J/(J+T*(F1_const+F1_slope*Fxc+comm*E))
    need(initial_fraction > Q('0.999999999'), 'initial term dominates this specific error majorant')
    report['initial_free_wave_fraction_of_differentiated_majorant_numerator'] = initial_fraction
    report['actual_Fx_contribution_to_second_error_upper'] = Z_slope*Fxc
    report['named_checks'] = checks
print(json.dumps(encode(report), indent=2, sort_keys=True))
