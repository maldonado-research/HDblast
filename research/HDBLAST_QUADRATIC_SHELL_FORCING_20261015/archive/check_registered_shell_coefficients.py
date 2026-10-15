"""Independent exact audit of the registered shell tension and Taylor sources.

Only copied, authenticated public inputs are read. No historical program is
imported or executed. Symbolic equalities check derivative coefficients, not
coefficients of epsilon squared without the Taylor factorial.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import sys

import sympy as s


ROOT = Path(__file__).resolve().parent.parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-error', choices=['third-sign', 'material-factor',
                                                  'omit-normal', 'metric-factor'])
    args = parser.parse_args()
    groups = []

    def check(name, condition):
        require(condition, name)
        groups.append(name)

    def equal(a, b):
        return s.expand(a-b) == 0

    pins_raw = (ROOT/'inputs/INPUT_PINS.json').read_bytes()
    pins = json.loads(pins_raw)
    require(pins['status'] == 'PASS', 'Authenticated public input pins required')
    for row in pins['files']:
        relative = PurePosixPath(row['path'])
        require(not relative.is_absolute() and '..' not in relative.parts
                and relative.as_posix() == row['path'], 'Canonical copied input path')
        path = ROOT/row['path']
        require(not any(p.is_symlink() for p in [path, *path.parents])
                and path.is_file(), 'Real copied input required')
        raw = path.read_bytes()
        require(len(raw) == row['bytes'] and sha256(raw).hexdigest() == row['sha256'],
                'Copied public input differs: '+row['path'])
    check('all_nine_selected_public_inputs_match_pins', len(pins['files']) == 9)
    box = json.loads((ROOT/'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_bytes())
    delta = F(box['parameters_exact']['delta'])
    c = F(box['parameters_exact']['c'])
    check('same_registered_exact_delta_and_model_c',
          delta == F(1,1000) and c == F(96714321425692,161839258930731))

    reference = (ROOT/'inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md').read_text()
    boundary = (ROOT/'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md').read_text()
    check('actual_registered_cubic_tension_and_one_copy_convention_retained',
          r'+\mathfrak d(1+c_*+c_*z)' in reference and
          r"\sigma'=4z+2z^2+\mathfrak dc_*" in reference and
          r"p=-\sigma'_0/2" in boundary and
          r'\tfrac12\varepsilon^2h_2' in boundary and
          r"S_p=-4pu^2+2\sigma''uv-2v^2" in boundary)

    eta, d, model_c = s.symbols('eta d model_c')
    u, v, a, p, sp, s2, s3, uy, vy, eps, w = s.symbols(
        'u v a p sp s2 s3 uy vy eps w')
    sigma = 2*(s.Rational(1,3)+eta**2+eta**3/3)+d*(1+model_c+model_c*eta)
    sigma1 = 4*eta+2*eta**2+d*model_c
    sigma2 = 4+4*eta
    sigma3 = s.Integer(4)
    expected3 = -sigma3 if args.inject_error == 'third-sign' else sigma3
    check('registered_polynomial_three_derivatives',
          equal(s.diff(sigma,eta),sigma1) and
          equal(s.diff(sigma,eta,2),sigma2) and
          equal(s.diff(sigma,eta,3),expected3))

    lo, hi = map(F,box['boxes']['eta_b']['exact'])
    check('same_near_positive_phi_branch_and_strict_sigma_prime_monotonicity',
          -1 < lo < hi < 0 and 4+4*lo > 0)
    sigma_prime = lambda e: 4*e+2*e*e+delta*c
    sharp1 = [sigma_prime(lo),sigma_prime(hi)]
    sharp2 = [4+4*lo,4+4*hi]
    retained1 = list(map(F,box['boxes']['sigma_prime_b']['exact']))
    retained2 = list(map(F,box['boxes']['sigma_second_b']['exact']))
    pb = list(map(F,box['boxes']['p_b']['exact']))
    check('correlated_monotone_sigma_prime_enclosure_refines_retained_box',
          retained1[0] < sharp1[0] < sharp1[1] < retained1[1])
    check('sigma_second_exact_affine_enclosure_matches_retained_box', sharp2 == retained2)
    check('retained_shell_identity_sigma_prime_equals_minus_two_p',
          retained1 == [-2*pb[1],-2*pb[0]] and sharp1[0] > 0 and pb[1] < 0)

    # One component of the tangential metric suffices because the first gauge
    # is conformal there. Its background derivative is 2 a h0, and u_y=uy.
    lapse_inverse = (1-4*eps*u)**s.Rational(-1,2)
    curvature = (a*(1+2*eps*u)+eps*uy)*lapse_inverse
    normal_scalar = (p+eps*vy)*lapse_inverse
    normal_K2 = s.diff(curvature,eps,2).subs(eps,0)
    normal_p2 = s.diff(normal_scalar,eps,2).subs(eps,0)
    check('unit_normal_second_metric_derivative',
          equal(normal_K2,4*u*uy+20*a*u**2))
    check('unit_normal_second_scalar_derivative',
          equal(normal_p2,4*u*vy+12*p*u**2))

    first_rows = {uy:-2*a*u-p*v/3,vy:-2*p*u-s2*v/2,sp:-2*p}
    metric_source = (2*sp*u*v/3+s2*v**2/6-normal_K2).subs(first_rows)
    metric_expected = s2*v**2/6-12*a*u**2
    if args.inject_error == 'metric-factor':
        metric_expected = s2*v**2/3-12*a*u**2
    check('complete_Eulerian_metric_source_after_first_shell_rows',
          equal(metric_source,metric_expected))
    scalar_source = (-s3*v**2/2-normal_p2).subs(first_rows).subs(s3,4)
    scalar_expected = -4*p*u**2+2*s2*u*v-2*v**2
    if args.inject_error == 'omit-normal':
        scalar_expected = -2*v**2
    check('complete_Eulerian_scalar_source_after_first_shell_rows',
          equal(scalar_source,scalar_expected))

    phi_path = eta+eps*v+eps**2*w/2
    scalar_tension_path = sigma1.subs(eta,phi_path)
    second_tension = s.diff(scalar_tension_path,eps,2).subs(eps,0)
    check('factorial_second_scalar_tension_chain_rule',
          equal(second_tension,sigma2*w+sigma3*v**2))
    material = -sigma3*v**2/2
    material_expected = -v**2 if args.inject_error == 'material-factor' else -2*v**2
    check('one_copy_material_scalar_source_has_correct_factor_and_sign',
          equal(material,material_expected) and equal(2*material,-4*v**2))

    wrong = [
        ('quadratic_tension_wrong_third_derivative',sigma3,s.Integer(0)),
        ('wrong_sign_third_derivative',sigma3,s.Integer(-4)),
        ('doubled_third_derivative',sigma3,s.Integer(8)),
        ('wrong_material_scalar_sign',material,2*v**2),
        ('lost_Taylor_factor_in_material_source',material,-v**2),
        ('scalar_source_omits_normal_terms',scalar_source,-2*v**2),
        ('wrong_Eulerian_metric_v_squared_factor',metric_source,s2*v**2/3-12*a*u**2),
        ('wrong_Eulerian_scalar_u_squared_sign',scalar_source,4*p*u**2+2*s2*u*v-2*v**2),
    ]
    for name, correct, incorrect in wrong:
        require(not equal(correct,incorrect), 'Wrong formula unexpectedly passed: '+name)
    check('eight_embedded_factor_sign_and_wrong_model_controls_rejected', len(wrong) == 8)

    def interval(values):
        return {'exact':[str(x) for x in values],
                'display_only':[float(x) for x in values]}

    result = {
        'schema':'hdblast-registered-shell-coefficient-audit-v1','status':'PASS',
        'scope':'Independent exact rational and symbolic coefficient checks on authenticated public inputs. Retains conditional actual background theorem; no second-order evolution, nonlinear remainder, physical calibration or external peer review is established.',
        'input_pins_sha256':sha256(pins_raw).hexdigest(),
        'parameters_exact':{'delta':str(delta),'c':str(c)},
        'tension_in_eta':'2/3 + 2 eta^2 + 2 eta^3/3 + delta (1+c+c eta)',
        'sigma_prime':'4 eta + 2 eta^2 + delta c',
        'sigma_second':'4+4 eta', 'sigma_third':'4',
        'eta_b':interval([lo,hi]),
        'sigma_prime_b_correlated_monotone':interval(sharp1),
        'sigma_prime_b_retained':interval(retained1),
        'sigma_second_b':interval(sharp2),
        'p_b_correlated_monotone':interval([-sharp1[1]/2,-sharp1[0]/2]),
        'material_single_copy_source':'-2 psi_1^2',
        'material_two_copy_source':'-4 psi_1^2',
        'Eulerian_metric_source':'sigma_second v^2/6 - 12 a u^2',
        'Eulerian_scalar_source':'-4 p u^2 + 2 sigma_second u v - 2 v^2',
        'gauge_scope':'The physical material first trace is psi_1=S(t)Y in the retained straight-shell observable construction. It must not be identified silently with a wave-gauge coordinate scalar trace. The -2 psi_1^2 term is the material tension contribution, not the complete coordinate scalar source or a solved second-order response.',
        'positive_check_groups':len(groups),'groups':groups,
        'embedded_wrong_formula_controls_rejected':len(wrong),
        'embedded_control_names':[x[0] for x in wrong],
        'python_optimized':not __debug__,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','positive_check_groups':len(groups),
                      'embedded_controls_rejected':len(wrong)}))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print('REJECTED: '+str(exc),file=sys.stderr)
        raise
