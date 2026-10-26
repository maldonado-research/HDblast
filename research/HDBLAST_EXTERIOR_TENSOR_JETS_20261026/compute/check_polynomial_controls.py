#!/usr/bin/env python3
"""Exact polynomial controls for cone/line weights and all six two-jet slots."""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from fractions import Fraction as F
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--evaluator', type=Path, required=True)
    parser.add_argument('--evaluator-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT_REQUIRED')
    need(not args.evaluator.is_symlink() and sha(args.evaluator) == args.evaluator_sha256, 'EVALUATOR_EXTERNAL_PIN')
    own_pin = sha(Path(__file__))
    spec = importlib.util.spec_from_file_location('externally_pinned_TT_twojet', args.evaluator)
    numeric = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(numeric)
    repo = args.repo_root.resolve()
    numeric.guard_output(repo, args.output)
    before = numeric.authenticate(repo)
    frozen = numeric.load('frozen_TT_for_polynomial_controls', repo/'research'/numeric.LABELS['F']/'tensor/tt_response.py')
    base = frozen.Context(repo)
    I = base.I

    class PolynomialContext:
        def __init__(self, coefficients):
            self.coefficients, self.base = coefficients, base

        def forcing_with_radial(self, x, r, audit=False):
            a, b, c, d, e = self.coefficients
            return a+b*x+c*r+d*x*r+e*r*r, c+d*x+2*e*r

    op = {'q': F(0), 'q_r': F(0)}
    apex = {
        'coefficient_intervals': {key: I(0) for key in ('A', 'B', 'D')},
        'background_value_intervals': {'rho': I(2), 'a': I(3), 'ell': base.ell, 'L': I(1)/(2*base.ell)},
    }
    checks, retained = [], {}
    for coefficients in ((1, 0, 0, 0, 0), (1, 2, 3, 4, 5), (0, 0, 0, 0, -1)):
        ctx = PolynomialContext(coefficients)
        a, b, c, d, e = map(F, coefficients)
        T, R = numeric.T, numeric.R
        C, H = a+c*R+e*R*R, b+d*R
        exact = {
            'u': C*T*T/2+H*T**3/6+e*T**4/12,
            'u_x': C*T+H*T*T/2+e*T**3/3,
            'u_r': (c+2*e*R)*T*T/2+d*T**3/6,
            'u_rr': e*T*T,
            'u_xr': (c+2*e*R)*T+d*T*T/2,
            'u_xx': C+H*T+e*T*T,
        }
        for n, segments in ((4, 8), (8, 16)):
            with contextlib.redirect_stdout(io.StringIO()):
                row = numeric.campaign_run(ctx, n, segments, op, apex, ctx.forcing_with_radial(I(T), I(R)))
            need(all(row['total_apex_twojet'][key].lo <= value <= row['total_apex_twojet'][key].hi for key, value in exact.items()), 'EXACT_POLYNOMIAL_TWOJET_ENCLOSURE')
            # Also verify the ordinary-coordinate conversion, including a L ur.
            ell, L, av = base.ell, 1/(2*base.ell), F(3)
            exact_natural = {'u': exact['u'], 'u_t': exact['u_x']/ell, 'u_y': -L*exact['u_r'], 'u_tt': exact['u_xx']/ell**2, 'u_ty': -L*exact['u_xr']/ell, 'u_yy': L**2*exact['u_rr']+av*L*exact['u_r']}
            need(all(row['ordinary_t_y_coefficient_twojet'][key].lo <= value <= row['ordinary_t_y_coefficient_twojet'][key].hi for key, value in exact_natural.items()), 'EXACT_POLYNOMIAL_NATURAL_COORDINATE_CONVERSION')
            checks.append({'polynomial_coefficients': coefficients, 'cone_n': n, 'segments_per_characteristic': segments, 'all_six_exact_derivatives': exact, 'computed_intervals': row['total_apex_twojet'], 'ordinary_t_y_exact_derivatives': exact_natural, 'all_exact_derivatives_enclosed': True})
            retained[coefficients, n] = row
    constant = retained[(1, 0, 0, 0, 0), 4]
    polynomial = retained[(1, 2, 3, 4, 5), 4]
    wrong_ur = constant['total_apex_twojet']['u_x']
    wrong_uxx = constant['total_apex_twojet']['u_rr']
    lines = polynomial['characteristic_quadrature']['half_integrals']
    wrong_urr = lines['Fr_plus']+lines['Fr_minus']
    expected_urr = F(5)*numeric.T**2
    controls = {
        'wrong_characteristic_difference_rejected': not (wrong_ur.lo <= 0 <= wrong_ur.hi),
        'omitted_apex_source_in_PDE_rejected': not (wrong_uxx.lo <= 1 <= wrong_uxx.hi),
        'wrong_second_radial_characteristic_sign_rejected': not (wrong_urr.lo <= expected_urr <= wrong_urr.hi),
        'boundary_triangles_as_full_squares_rejected': numeric.T**2/2 != (numeric.T**2+4*(numeric.T/4)**2)/2,
    }
    need(all(controls.values()), 'POLYNOMIAL_CORRUPTION_CONTROL')
    base.finish()
    need(numeric.authenticate(repo) == before, 'FROZEN_CLOSURES_CHANGED')
    need(sha(args.evaluator) == args.evaluator_sha256 and sha(Path(__file__)) == own_pin, 'CONTROL_SOURCES_CHANGED')
    result = {'schema': 'HDBLAST_TT_TWOJET_ANALYTIC_POLYNOMIAL_CONTROLS_V1', 'status': 'PASS_EXACT_POLYNOMIAL_GEOMETRY_AND_TWOJET_CONTROLS', 'checker_sha256': own_pin, 'evaluator_sha256': args.evaluator_sha256, 'authentication': before, 'forcing_family': 'F(s,z)=a+b*s+c*z+d*s*z+e*z²; Fr=c+d*s+2e*z; exact zero-data scalar wave solution with A=B=D=0', 'exact_solution': 'u(x,r)=(a+c*r+e*r²)*x²/2+(b+d*r)*x³/6+e*x⁴/12', 'controls': checks, 'deliberate_algebra_corruptions_rejected': controls, 'actual_scientific_campaign_replayed': False}
    with args.output.open('x') as handle:
        json.dump(numeric.encode(result), handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(result['status'])


if __name__ == '__main__':
    main()
