#!/usr/bin/env python3
"""Independent exact bridge/parity/half-line audit; no physical theorem is inferred.

Positive checks require zero rational residual. Negative controls require a
nonzero residual and have explicit names. No Python assert determines outcomes.
"""
from __future__ import annotations

import hashlib
import json
import platform
import argparse
from fractions import Fraction
from pathlib import Path

import sympy as s

HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--basis', type=Path, default=PACKAGE_ROOT / 'basis',
                    help='Dated proof inputs; default is the sibling basis directory.')
parser.add_argument('--output', type=Path, default=HERE / 'SCALAR_EXACT_CHECKS.json',
                    help='JSON receipt destination; parent directories are created.')
arguments = parser.parse_args()
BASIS_EXPECTED = {'static/theory/JACOBI_ENDPOINT_MONOTONICITY.md': 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8', 'static/bounds/RESCALED_BOUND_CERTIFICATE.json': '22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d', 'bridge/PROSPECTIVE_ZERO_MODE_STRATEGY.md': '6d816800055083af4d0637eb725d032e01e2f4e03dece2e6a76519970e7b5cef'}
source_pins = {}
for logical_name, expected_sha in BASIS_EXPECTED.items():
    path = arguments.basis / logical_name
    if not path.is_file():
        raise FileNotFoundError('Missing dated basis input: ' + logical_name)
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    if observed != expected_sha:
        raise ValueError('Dated basis SHA256 mismatch: ' + logical_name)
    source_pins['basis/' + logical_name] = observed

rho, a, p, h, q, v, b, U1, U2, U3, S2, lam = s.symbols(
    'rho a p h q v b U1 U2 U3 S2 Lambda', nonzero=True)
H = rho**-2
flow = {rho: a*rho, a: -H-p*p/3, p: U1-4*a*p,
        h: q, q: U2*h-4*a*q-4*p*b, v: b,
        b: 2*H*v-2*p*q/3, U1: U2*p, U2: U3*p}

def derivative(expr, extra=None):
    f = dict(flow)
    if extra:
        f.update(extra)
    return s.expand(sum(s.diff(expr, x)*dx for x, dx in f.items()))

# Eliminate b, unlike the archived transform checker which eliminates v.
constraint_b = (p*q-U1*h-12*H*v)/(12*a)
checks = {}
controls = {}

def residual(expr, constraint=False):
    if constraint:
        expr = expr.subs(b, constraint_b)
    return s.factor(s.cancel(expr))

def check(name, expr, constraint=False):
    r = residual(expr, constraint)
    checks[name] = {'passed': r == 0, 'residual': str(r)}

def wrong(name, expr, constraint=False):
    r = residual(expr, constraint)
    controls[name] = {'detected': r != 0, 'residual': str(r)}

j = b+p*h/3
D = h*flow[p]-p*q
K = rho**2*j
X0 = h+p*K
Psi0 = v+a*K
g = flow[p]/p
B = g+S2/2
J = q/h+S2/2
M0 = derivative(X0)+S2*X0/2+2*p*Psi0
check('Jacobi_j_damping', derivative(j)+4*a*j+2*H*v, True)
check('Jacobi_Wronskian_evolution', derivative(D)+4*a*D-4*p*(H*h+p*j))
check('K_damping', derivative(K)+2*a*K+2*v, True)
check('zero_mode_momentum', derivative(Psi0)+2*a*Psi0+p*X0/3, True)
check('zero_mode_Gauss', p*derivative(X0)-flow[p]*X0-(12*H-2*p*p)*Psi0, True)
check('zero_mode_Psi_Wronskian', Psi0+rho**2*D/12, True)
check('B_minus_J_Wronskian', B-J-D/(p*h))
check('zero_mode_full_shell_cancellation', M0-h*J-B*p*K, True)
check('zero_mode_XZ_mismatch', M0-B*X0-12*H*Psi0/p, True)
check('static_junction_determinant', j*p*B+H*h*J-H*M0, True)
check('constant_odd_continuity_determinant', a*h-p*v-(a*X0-p*Psi0))

wrong('omit_Jacobi_K_correction_breaks_momentum', derivative(v)+2*a*v+p*h/3, True)
wrong('flip_K_sign_breaks_momentum', derivative(v-a*K)+2*a*(v-a*K)+p*(h-p*K)/3, True)
wrong('omit_gravitational_shell_term_breaks_cancellation',
      derivative(X0)+S2*X0/2-h*J-B*p*K, True)
wrong('reverse_scalar_shell_tension_sign_breaks_cancellation',
      derivative(X0)-S2*X0/2+2*p*Psi0-h*J-B*p*K, True)

# Reconstruct the second-order expression from generic first-order equations.
X, Psi = s.symbols('X Psi')
gi = {X: g*X+(3*lam*H/p-2*p)*Psi, Psi: -2*a*Psi-p*X/3}
P = rho**2*Psi
pSL = 1/(rho**2*p*p)
qSL = 2/(3*rho**2)
wSL = 1/(rho**4*p*p)
Pprime = derivative(P, gi)
flux = pSL*Pprime
beta = rho**4*p*p*B
M = B*X+3*lam*H*Psi/p
check('SL_flux_from_momentum', flux+X/(3*p))
check('SL_expression_from_Gauss', -derivative(flux, gi)+qSL*P-lam*wSL*P)
check('symmetric_boundary_orientation', flux-lam*P/beta+M/(3*p*B))
P0 = rho**2*Psi0
r0 = pSL*derivative(P0)/P0-4/beta
check('reference_boundary_positive_margin', r0+M0/(3*rho**2*p*B*Psi0), True)
wrong('reverse_SL_boundary_sign', flux+lam*P/beta+M/(3*p*B))
wrong('use_rho2_instead_of_rho4_boundary_beta',
      flux-lam*P/(rho**2*p*p*B)+M/(3*p*B))
wrong('replace_qSL_by_negative_value', -derivative(flux, gi)-qSL*P-lam*wSL*P)

# Generic full two-copy shell constraints. The summed trace-free equation and
# induced-field continuity split into even and odd parts.
ze, zo, Xe, Xo, Pe, Po = s.symbols('zeta_even zeta_odd X_even X_odd Psi_even Psi_odd')
zp, zm = ze+zo, ze-zo
Xp, Xm, Psip, Psim = Xe+Xo, Xe-Xo, Pe+Po, Pe-Po
check('two_copy_TF_sum_selects_even_bending', (zp+zm)/2-ze)
check('two_copy_scalar_continuity_selects_odd_field',
      ((Xp+p*zp)-(Xm+p*zm))/2-(Xo+p*zo))
check('two_copy_metric_continuity_selects_odd_field',
      ((Psip+a*zp)-(Psim+a*zm))/2-(Po+a*zo))
odd_boundary = {X: p*Psi/a}
co = 1/(3*a*rho**2)
check('odd_Robin_orientation', (flux+co*P).subs(odd_boundary))
check('odd_shell_volume_kinetic_suggestion',
      (6*a)*rho*rho*(P/(a*rho*rho))**2-18*co*P*P)
wrong('odd_Robin_wrong_sign', (flux-co*P).subs(odd_boundary))
wrong('odd_continuity_wrong_scalar_sign', (flux+co*P).subs(X, -p*Psi/a))
dplus, dminus = s.symbols('Delta_plus Delta_minus')
ell1_matrix = s.Matrix([[p, -p], [p*B, p*B]])
check('ell1_two_copy_scalar_constraint_determinant', ell1_matrix.det()-2*p*p*B)

# Local Picone identity, independent abstract coefficient names and test field.
ps, dps, qs, ws, ref, dref, test, dtest, l0 = s.symbols(
    'pSL dpSL qSL wSL reference dreference test dtest Lambda0', nonzero=True)
ddref = ((qs-l0*ws)*ref-dps*dref)/ps
picone_total = derivative(ps*dref/ref*test*test,
    {ps:dps, ref:dref, dref:ddref, test:dtest})
picone_square = ps*ref*ref*(dtest/ref-test*dref/ref**2)**2
check('local_Picone_divergence_orientation',
      ps*dtest*dtest+qs*test*test-l0*ws*test*test-picone_square-picone_total)
ref1, ref2, f1, f2, l1, l2 = s.symbols('u1 u2 flux1 flux2 Lambda1 Lambda2')
wronskian_derivative = derivative(ref1*f2-f1*ref2,
 {ref1:f1/ps, ref2:f2/ps, f1:(qs-l1*ws)*ref1, f2:(qs-l2*ws)*ref2})
check('SL_Wronskian_reality_bulk_orientation',
      wronskian_derivative-(l1-l2)*ws*ref1*ref2)
wrong('wrong_Picone_endpoint_orientation',
      ps*dtest*dtest+qs*test*test-l0*ws*test*test-picone_square+picone_total)

# Half-line x = integral_y^T dt/rho(t), x increases toward the cone.
# A signed Liouville factor is allowed because p has constant nonzero sign.
F, Ft, Ftt = s.symbols('F F_x F_xx')
scale = rho**s.Rational(-3,2)/p
half_flow = {F:-Ft/rho, Ft:-Ftt/rho}
P_liouville = F/scale
fprime = derivative(P_liouville, half_flow)
L = -derivative(pSL*fprime, half_flow)+qSL*P_liouville-lam*wSL*P_liouville
gprime = derivative(g)
V_before = rho**2*(s.Rational(3,4)*a*a+2*a*g+g*g-s.Rational(3,2)*flow[a]-gprime+2*p*p/3)
V = rho**2*(s.Rational(3,4)*a*a+6*a*g+2*g*g-U2-p*p/6)-s.Rational(5,2)
check('half_line_potential_background_simplification', V_before-V)
check('half_line_transformed_differential_expression', scale*L/wSL+Ftt-(V-lam)*F)
kappa = rho*(s.Rational(3,2)*a+g)
check('half_line_even_boundary_orientation',
      (pSL*fprime-lam*P_liouville/beta)*rho**3*p*p*scale-
      (-Ft+kappa*F-lam*F/(rho*B)))
check('half_line_odd_boundary_orientation',
      (pSL*fprime+co*P_liouville)*rho**3*p*p*scale-
      (-Ft+(kappa+rho*p*p/(3*a))*F))
check('half_line_L2_unitarity', wSL*rho/scale**2-1)
check('half_line_even_endpoint_mass', 1/(beta*scale**2)-1/(rho*B))
wrong('half_line_reverse_x_orientation',
      (pSL*fprime-lam*P_liouville/beta)*rho**3*p*p*scale-
      (Ft+kappa*F-lam*F/(rho*B)))
wrong('half_line_omit_pprime_from_Liouville_factor',
      V-rho**2*(s.Rational(3,4)*a*a-s.Rational(3,2)*flow[a]+2*p*p/3))

# Kinetic comparison is a bulk algebraic identity, NOT an action reduction.
check('bulk_kinetic_equals_18_SL_energy_density',
      2*rho*rho*(X*X+6*Psi*Psi)-18*(pSL*Pprime**2+qSL*P**2))

y, Ac, Uc = s.symbols('y A_center U_center', positive=True)
check('cone_half_line_potential_threshold',
      s.limit(V.subs({rho:y,a:1/y,p:Ac*y,U1:5*Ac,U2:Uc}),y,0)-s.Rational(25,4))
gamma = s.symbols('gamma')
cone_indicial = -(s.diff(y**-4*s.diff(y**gamma,y),y))+lam*y**(gamma-6)
check('cone_SL_indicial_polynomial',
      s.expand(cone_indicial/y**(gamma-6))+gamma*gamma-5*gamma-lam)
check('cone_reference_weight_cancellation',
      (1/(Ac*Ac*y**6))*(-Ac*y**4/12)**2-y*y/144)
check('cone_reference_energy_cancellation',
      (1/(Ac*Ac*y**4))*(-Ac*y**3/3)**2-y*y/9)
eta, c, delta = s.symbols('eta c delta')
phi = eta+1
check('exact_shift_W', 1-phi+phi**3/3-(s.Rational(1,3)+eta*eta+eta**3/3))
check('exact_shift_detuning', 1+c*phi-(1+c+c*eta))

# Independent covariant graph-boundary calculation. For n^2 dy^2+r^2 gamma,
# keeping two powers of tangential graph derivatives gives
# sqrt(h_graph)K_graph-sqrt(h_coordinate)K_y
# = (n*r*r_y-n_y*r^2)(D f)^2 - n^-1 D_mu(n^2*r^2 D^mu f).
# Expand at y=T+epsilon*zeta*Q, with n=exp(epsilon*N*Q),
# r=rho*exp(epsilon*psi*Q). A homogeneous de Sitter chart realizes
# (DQ)^2=-Q_t^2 and Box Q=-(Q_tt+3Q_t). This implementation starts from
# the covariant density, independently of the operator level-set code.
Nb, psib, zeta, Qamp, Qdot, Qddot = s.symbols(
    'N_boundary psi_boundary zeta Q Q_t Q_tt')
def time_derivative_coefficient(density):
    expanded = s.expand(density)
    coefficient = s.diff(expanded, Qddot)
    if coefficient.has(Qddot, Qdot):
        raise ValueError('Unexpected time-derivative structure in graph density')
    integrated = s.expand(expanded-coefficient*Qddot-
        (s.diff(coefficient,Qamp)*Qdot+3*coefficient)*Qdot)
    return s.simplify(s.diff(integrated,Qdot,2)/2)

graph_minus_radial_density = rho**2*(
    -a*zeta**2*Qdot**2
    +zeta*(Nb+2*psib)*Qamp*(Qddot+3*Qdot)
    +2*zeta*(Nb+psib)*Qdot**2
    +2*a*zeta**2*Qamp*(Qddot+3*Qdot))
graph_before_lapse_flux = time_derivative_coefficient(graph_minus_radial_density)
lapse_moving_domain_flux = -rho**2*Nb*zeta*Qdot**2
graph_complete = time_derivative_coefficient(
    graph_minus_radial_density+lapse_moving_domain_flux)
moving_intrinsic_density = 3*rho**2*zeta*psib*Qamp*(Qddot+3*Qdot)
moving_intrinsic = time_derivative_coefficient(moving_intrinsic_density)
sigmab = s.symbols('sigma_boundary')
shell_kinetic = sigmab*rho**2*zeta**2/2
total_boundary_time_kinetic = 2*(graph_complete+moving_intrinsic)+shell_kinetic
odd_boundary_time_kinetic = s.factor(total_boundary_time_kinetic.subs(
    {psib:-a*zeta,sigmab:6*a}))
check('graph_GHY_radial_difference_lapse_term',
      graph_before_lapse_flux-rho**2*(-3*a*zeta**2+Nb*zeta))
check('graph_EH_lapse_divergence_cancellation',
      graph_complete+3*a*rho**2*zeta**2)
check('moving_intrinsic_linear_curvature_time_kinetic',
      moving_intrinsic+3*rho**2*zeta*psib)
check('full_odd_graph_shell_time_kinetic',
      odd_boundary_time_kinetic-3*a*rho**2*zeta**2)
check('full_odd_boundary_K_matches_Robin_energy_trace',
      (2*odd_boundary_time_kinetic).subs(zeta,-P/(a*rho**2))-18*co*P**2)
Ebulk, Wbulk = s.symbols('odd_bulk_energy odd_bulk_norm')
check('full_odd_total_K_matches_spectral_bulk_norm',
      (18*Ebulk+(2*odd_boundary_time_kinetic).subs(
          zeta,-P/(a*rho**2))).subs(Ebulk,lam*Wbulk-co*P**2)-18*lam*Wbulk)
wrong('omit_EH_lapse_flux_breaks_graph_kinetic',
      graph_before_lapse_flux-graph_complete)
wrong('omit_moving_domain_curvature_breaks_odd_kinetic',
      (2*graph_complete+shell_kinetic).subs(
          {psib:-a*zeta,sigmab:6*a})-odd_boundary_time_kinetic)
wrong('reverse_shell_action_sign_breaks_odd_kinetic',
      (2*(graph_complete+moving_intrinsic)-shell_kinetic).subs(
          {psib:-a*zeta,sigmab:6*a})-odd_boundary_time_kinetic)
wrong('reverse_odd_metric_matching_sign_breaks_kinetic',
      total_boundary_time_kinetic.subs(
          {psib:a*zeta,sigmab:6*a})-odd_boundary_time_kinetic)

# A nontrivial exact Hardy witness: g(y)=y^(-7/5)-1 on (0,1].
# Both weighted integrals are finite, g(1)=0, and I/J=5/12. A corrupted
# Hardy constant of 2/5 therefore fails on this actual function; this control
# checks an inequality witness rather than merely subtracting two constants.
power = s.Rational(7,5)
hardy_numerator = 1/(3-2*power)-2/(3-power)+s.Rational(1,3)
hardy_denominator = power**2/(3-2*power)
hardy_ratio = s.factor(hardy_numerator/hardy_denominator)
check('Hardy_exact_integral_witness_ratio',hardy_ratio-s.Rational(5,12))
checks['Hardy_correct_constant_accepts_actual_witness'] = {
    'passed':bool(hardy_numerator<=s.Rational(4,9)*hardy_denominator),
    'residual':str(s.Rational(4,9)*hardy_denominator-hardy_numerator)}
controls['Hardy_corrupted_constant_rejected_by_actual_witness'] = {
    'detected':bool(hardy_numerator>s.Rational(2,5)*hardy_denominator),
    'residual':str(hardy_numerator-s.Rational(2,5)*hardy_denominator)}

# Map exact registered model and stored binary64 values separately.
cexact = 2/Fraction('1.0357712571566784')-Fraction(4,3)
cbinary = Fraction.from_float(2.0/1.0357712571566784-4.0/3.0)
f0stored = Fraction.from_float(1.0+float(cbinary))
models = {
 'registered_exact_decimal': {'k':'1/9', 'c':str(cexact),'f0':str(1+cexact)},
 'binary64_c_only_as_exact_model': {'k':'1/9','c':str(cbinary),'f0':str(1+cbinary)},
 'stored_binary64_coefficients_as_exact_model': {
    'k':str(Fraction.from_float(1/9)),'c':str(cbinary),'f0':str(f0stored)},
 'exact_c_difference_binary_minus_registered':str(cbinary-cexact),
 'exact_f0_stored_minus_exact_addition':str(f0stored-(1+cbinary)),
 'warning':'All three mathematical coefficient models differ from executing binary64 arithmetic at every ODE step; no numerical trajectory is certified.'}

passed = all(item['passed'] for item in checks.values())
detected = all(item['detected'] for item in controls.values())
report = {
 'basis_date_utc':'2026-10-11',
 'continuation_date_utc':'2026-10-13',
 'scope':'Independent exact algebra and conditional parity/domain audit; no complete physical scalar stability or mass-gap certificate.',
 'python':platform.python_version(),'sympy':s.__version__,
 'positive_check_count':len(checks),'positive_checks_passed':passed,
 'deliberate_wrong_control_count':len(controls),'all_wrong_controls_detected':detected,
 'checks':checks,'wrong_controls':controls,'coefficient_models':models,
 'source_sha256':source_pins,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'graph_boundary_time_kinetic':{
    'GHY_without_EH_lapse_flux':str(graph_before_lapse_flux),
    'EH_plus_GHY_graph':str(graph_complete),
    'moving_intrinsic':str(moving_intrinsic),
    'odd_total_boundary':str(odd_boundary_time_kinetic)},
 'Hardy_actual_counterexample_to_corrupted_constant':{
    'function':'y^(-7/5)-1, 0<y<=1',
    'weighted_L2_integral':str(hardy_numerator),
    'weighted_derivative_integral':str(hardy_denominator),
    'ratio':str(hardy_ratio),'correct_constant':'4/9','corrupted_constant':'2/5'},
 'claim_limits':[
  'Full two-copy odd junction equations require independently audited orientation and matching; algebra alone does not derive them.',
  'Finite norm at the cone is a specified state class, not a universal boundary condition for all possible extensions.',
  'Local odd action completion is checked here and in the operator note; global gauge, cone symplectic flux and physical state-space completeness remain conditional.',
  'Closure, self-adjointness, exceptional-harmonic completeness, and global gauge admissibility require mathematical/physical argument.',
  'The reviewed domain proof derives fixed-background scalar gaps; no detuning-uniform dimensionless gap or uniform positive physical mass floor is proved.']}
destination = arguments.output
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
print(json.dumps({k:report[k] for k in ('positive_check_count','positive_checks_passed','deliberate_wrong_control_count','all_wrong_controls_detected')}))
if not passed or not detected:
    print(json.dumps({'failed_checks':[k for k,v in checks.items() if not v['passed']],
                      'undetected_controls':[k for k,v in controls.items() if not v['detected']]}))
raise SystemExit(0 if passed and detected else 1)
