#!/usr/bin/env python3
"""Exact S3 integration and quadrupole symmetry conditions, without numerics."""
import argparse
import json
from pathlib import Path
import sympy as s


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-failure', action='store_true',
                        help='Deliberately reverse the first rotation identity; must fail without a result file')
    args = parser.parse_args()
    require(not args.output.exists(), 'Output must be fresh')
    require(s.__version__ == '1.14.0', 'Registered SymPy 1.14.0 required')
    x = s.Matrix(s.symbols('x1:5', real=True))
    coordinates = tuple(x)

    def average(expression):
        """Exact monomial average over the unit S3; its volume is factored out."""
        total = s.S.Zero
        for powers, coefficient in s.Poly(s.expand(expression), *coordinates).terms():
            if any(p % 2 for p in powers):
                continue
            degree = sum(powers)
            numerator = s.prod(s.factorial2(p - 1) if p else 1 for p in powers)
            denominator = s.prod(4 + 2*j for j in range(degree // 2))
            total += coefficient * s.Rational(numerator, denominator)
        return s.expand(total)

    def matrix(prefix):
        a,b,c,d,e,f,g,h,j = s.symbols(prefix+'0:9', real=True)
        return s.Matrix([[a,d,e,f],[d,b,g,h],[e,g,c,j],[f,h,j,-a-b-c]])

    A, B = matrix('a'), matrix('b')
    HA, HB = (x.T*A*x)[0], (x.T*B*x)[0]
    checks, controls = [], []

    def zero(name, expression):
        residual = s.expand(expression)
        require(residual == 0, 'Exact identity failed: '+name)
        checks.append({'name': name, 'residual': '0'})

    def wrong(name, residual):
        reduced = s.expand(residual)
        require(reduced != 0, 'Wrong-formula control was not distinguished: '+name)
        controls.append({'name': name, 'incorrect_formula_rejected': True})

    zero('sphere_fourth_moments_normalize_radius_four', average((x.dot(x))**2)-1)
    zero('harmonic_A_trace_free', sum(s.diff(HA,z,2) for z in x))
    zero('harmonic_B_trace_free', sum(s.diff(HB,z,2) for z in x))
    zero('harmonic_A_zero_mean', average(HA))
    zero('harmonic_B_zero_mean', average(HB))
    zero('arbitrary_quadrupole_bilinear_sphere_integral', average(HA*HB)-s.trace(A*B)/12)
    rotation_rows = []
    for i in range(4):
        for j in range(i+1,4):
            M = s.zeros(4)
            M[i,j], M[j,i] = -1, 1
            JHA = x[i]*s.diff(HA,x[j])-x[j]*s.diff(HA,x[i])
            expected = (x.T*(A*M-M*A)*x)[0]
            if args.inject_failure and (i,j) == (0,1):
                expected = -expected
            zero(f'rotation_{i+1}{j+1}_polynomial_covariance', JHA-expected)
            integrated = average(HB*JHA)
            commutator = (A*B-B*A)[i,j]
            zero(f'rotation_{i+1}{j+1}_integral_is_commutator_component', integrated+commutator/6)
            rotation_rows.append({'generator':f'J_{i+1}{j+1}',
                                  'integral_divided_by_sphere_volume':str(s.factor(integrated)),
                                  'commutator_component':str(s.expand(commutator))})

    u,v,w = s.symbols('u v w', real=True)
    for i in range(4):
        # The three independent radial coefficients are kept arbitrary.
        # Radial differentiation does not alter these angular parities.
        angular = u*HA*x[i]*HA + v*HA*s.diff(HA,x[i]) + w*x[i]*HB*HB
        zero(f'boost_{i+1}_all_radial_coefficients_integrate_to_zero', average(angular))

    t = s.symbols('t', real=True)
    zero('common_harmonic_matrices_commute', sum(e**2 for e in (A*(t*A)-(t*A)*A)))

    H = x[0]**2-x[1]**2
    J = lambda value: x[0]*s.diff(value,x[1])-x[1]*s.diff(value,x[0])
    JH = J(H)
    zero('explicit_quadrupole_generator', JH+4*x[0]*x[1])
    zero('explicit_quadrupole_second_generator', J(JH)+4*H)
    zero('explicit_quadrupole_norm', average(H*H)-s.Rational(1,6))
    zero('explicit_rotated_quadrupole_norm', average(JH*JH)-s.Rational(2,3))
    zero('explicit_orthogonal_displacement_velocity_harmonics', average(H*JH))
    zero('explicit_doubled_Hessian_coefficient', 4*average(JH*JH)-s.Rational(8,3))

    wrong('wrong_sphere_trace_normalization', average(HA*HB)-s.trace(A*B)/24)
    wrong('reversed_rotation_commutator_sign', JH-4*x[0]*x[1])
    wrong('unrestricted_quadrupole_zero_claim', 4*average(JH*JH))
    wrong('one_interior_used_as_doubled_Hessian', 2*average(JH*JH)-s.Rational(8,3))

    result = {
        'status':'PASS_EXACT_QUADRUPOLE_SYMMETRY_CHECKS',
        'sympy_version':s.__version__,
        'sphere_dimension':3,
        'sphere_volume':'2*pi**2',
        'integration_convention':'Displayed exact moments are divided by Vol(S3)',
        'positive_checks':len(checks),
        'wrong_formula_controls':len(controls),
        'checks':checks, 'controls':controls,
        'rotation_components':rotation_rows,
        'witness':{'H':'n1**2-n2**2','J12_H':'-4*n1*n2',
                   'norm_H_divided_by_volume':'1/6',
                   'norm_JH_divided_by_volume':'2/3',
                   'doubled_Hessian_divided_by_omega_Krad_volume':'8/3'},
        'scope':'Exact harmonic/matrix algebra only; the nonlinear constraint identity, PDE existence and classical-family hypotheses require the accompanying proofs.',
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(result,stream,indent=2,sort_keys=True)
        stream.write('\n')
    print(f"PASS {len(checks)} exact identities; {len(controls)} wrong-formula controls")


if __name__ == '__main__':
    main()
