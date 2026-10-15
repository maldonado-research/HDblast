"""Exact pullback, gauge and normalization checks; not an analytic PDE proof."""
import argparse
import json
from pathlib import Path
import sys
import sympy as sp

FAULTS = ['drop_embedding', 'wrong_gauge_sign', 'half_hessian',
          'wrong_angular_norm', 'missing_proper_time', 'second_jet_leak', 'double_shell']
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--fault', choices=['none', *FAULTS], default='none')
parser.add_argument('--output', type=Path)
args = parser.parse_args()
checks = []


def equal(name, actual, expected):
    residual = sp.simplify(sp.expand(actual - expected))
    checks.append({'check': name, 'status': 'PASS' if residual == 0 else 'FAIL',
                   'residual': str(residual)})


e, dy, dt = sp.symbols('epsilon dy dt')
phib, p, py = sp.symbols('phi_b p p_y', real=True)
chi, chiy, chit, psi2 = sp.symbols('chi chi_y chi_t psi_2', real=True)
zeta, zeta2, clock1 = sp.symbols('zeta zeta_2 clock_1', real=True)
mu0, mu1, mu2 = sp.symbols('mu_0 mu_1 mu_2', real=True)
bulk = phib + p*dy + py*dy**2/2 + e*(chi + chiy*dy + chit*dt) + e**2*psi2/2
pulled = sp.expand(bulk.subs({dy: e*zeta + e**2*zeta2/2, dt: e*clock1}))
first = sp.diff(pulled, e).subs(e, 0)
candidate = chi if args.fault == 'drop_embedding' else chi + p*zeta
equal('physical scalar first variation from an explicit moving pullback', candidate, first)

gauge = sp.symbols('gauge_normal', real=True)
gauge_zeta = zeta + gauge if args.fault == 'wrong_gauge_sign' else zeta - gauge
transformed_first = chi + p*gauge + p*gauge_zeta
equal('bulk normal gauge shift cancels the shell embedding shift', transformed_first, first)
equal('first clock shift disappears because the background shell scalar is constant',
      sp.diff(first, clock1), 0)

measure = mu0 + e*mu1 + e**2*mu2/2
integrand = sp.expand((pulled - phib)**2 * measure)
hessian = sp.diff(integrand, e, 2).subs(e, 0)
candidate_hessian = 2*mu0*first**2
if args.fault == 'half_hessian':
    candidate_hessian /= 2
if args.fault == 'second_jet_leak':
    candidate_hessian += mu0*psi2
equal('second derivative of squared physical response includes its factor two',
      candidate_hessian, hessian)
equal('quadratic coefficient is independent of the second scalar/embedding/measure jets',
      sum(sp.diff(hessian, variable)**2 for variable in [psi2, zeta2, mu1, mu2]), 0)

# Hopf coordinates: n1=sqrt(u)cos(theta), n2=sqrt(u)sin(theta),
# n3=sqrt(1-u)cos(v), n4=sqrt(1-u)sin(v), dOmega=du dtheta dv/2.
u, theta, v = sp.symbols('u theta v', real=True)
angular = sp.integrate(u**2*sp.cos(2*theta)**2/2,
                       (u, 0, 1), (theta, 0, 2*sp.pi), (v, 0, 2*sp.pi))
candidate_angular = 2*sp.pi**2/3 if args.fault == 'wrong_angular_norm' else sp.pi**2/3
equal('actual unit S3 quadrupole norm from Hopf-coordinate integration', candidate_angular, angular)

rho, t, S = sp.symbols('rho_b t S', positive=True)
section_measure = rho**3*sp.cosh(t)**3
proper_time_jacobian = rho
window_density = angular*section_measure*S**2*proper_time_jacobian
candidate_window = sp.pi**2* rho**4 * sp.cosh(t)**3 * S**2 / 3
if args.fault == 'missing_proper_time':
    candidate_window /= rho
if args.fault == 'double_shell':
    candidate_window *= 2
equal('one physical shell and proper-time integral give the rho_b power four',
      candidate_window, window_density)

success = all(row['status'] == 'PASS' for row in checks)
result = {'schema': 'hdblast-physical-shell-observable-checks-v1',
          'status': 'PASS' if success else 'FAIL', 'fault': args.fault,
          'python_optimized': not __debug__, 'sympy_version': sp.__version__,
          'finite_positive_check_groups': len(checks), 'checks': checks,
          'scope': 'Finite symbolic pullback/gauge/normalization identities. These do not certify analytic existence, activation or nonlinear remainder bounds.'}
if args.output:
    args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
sys.exit(0 if success else 1)
