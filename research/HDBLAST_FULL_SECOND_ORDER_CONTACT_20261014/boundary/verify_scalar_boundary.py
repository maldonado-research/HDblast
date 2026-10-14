#!/usr/bin/env python3
"""Independent exact checks of the Eulerian scalar shell boundary rows.

This checks local geometry/Taylor/gauge identities, not a curved IBVP theorem.
No prior scientific code is imported. All guards survive Python optimization.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inject-error', action='store_true')
    args = parser.parse_args()
    checks: dict[str, str] = {}
    controls: dict[str, str] = {}

    def zero(name: str, value: s.Expr) -> None:
        result = s.simplify(s.expand(value))
        if result != 0:
            raise RuntimeError(f'{name}: nonzero exact residual {result}')
        checks[name] = 'PASS'

    def wrong(name: str, value: s.Expr) -> None:
        if s.simplify(s.expand(value)) == 0:
            raise RuntimeError(f'{name}: corrupt formula unexpectedly accepted')
        controls[name] = 'REJECTED'

    eps = s.Symbol('eps')
    a, p, u, uy, v, vy, s2, r, z = s.symbols('a p u uy v vy sigma2 rho zeta')
    sp = -2*p
    inv_lapse = (1-4*eps*u)**s.Rational(-1, 2)
    exact_K = inv_lapse * (a + eps*(uy+2*a*u))
    exact_p = inv_lapse * (p + eps*vy)
    K1 = s.diff(exact_K, eps).subs(eps, 0)
    K2 = s.diff(exact_K, eps, 2).subs(eps, 0)
    p1 = s.diff(exact_p, eps).subs(eps, 0)
    p2 = s.diff(exact_p, eps, 2).subs(eps, 0)
    zero('exact_first_extrinsic_curvature', K1-(uy+4*a*u))
    K2_claim = 4*u*uy + (18 if args.inject_error else 20)*a*u**2
    zero('quadratic_extrinsic_curvature_factorial', K2-K2_claim)
    zero('exact_first_normal_scalar', p1-(vy+2*p*u))
    zero('quadratic_normal_scalar_factorial', p2-(4*u*vy+12*p*u**2))
    zero('first_metric_shell_row', K1-2*a*u-sp*v/6-(uy+2*a*u+p*v/3))
    zero('first_scalar_shell_row', p1+s2*v/2-(vy+2*p*u+s2*v/2))
    first_rows = {uy: -2*a*u-p*v/3, vy: -2*p*u-s2*v/2}
    raw_K = 2*sp*u*v/3+s2*v**2/6-K2
    raw_p = -2*v**2-p2
    S_K = s2*v**2/6-12*a*u**2
    S_p = -4*p*u**2+2*s2*u*v-2*v**2
    zero('metric_forcing_on_first_shell_solution', raw_K.subs(first_rows)-S_K)
    zero('scalar_forcing_on_first_shell_solution', raw_p.subs(first_rows)-S_p)
    wrong('missing_K_factorial', (K2/2)-K2)
    wrong('wrong_scalar_cubic_sign', 2*v**2-(-2*v**2))
    wrong('discard_first_lapse_quadratic', raw_p-(-2*v**2))
    wrong('double_shell_tension_normalization', (-4*v**2)-(-2*v**2))
    phi, detuning, cstar = s.symbols('phi detuning cstar')
    qphi = phi-1
    sigma = 2*(s.Rational(1,3)+qphi**2+qphi**3/3)+detuning*(1+cstar+cstar*qphi)
    zero('registered_tension_sigma_third', s.diff(sigma,phi,3)-4)
    zero('registered_tension_sigma_second', s.diff(sigma,phi,2)-4*phi)

    # Independent chain-rule check in the family-adapted normal coordinate.
    # The coordinate acceleration derivative is NOT the generator derivative.
    n2, d2 = s.symbols('N2 coordinate_acceleration_y')
    jacobian = 1+2*eps*u+eps**2*d2/2
    transformed_yy = (1-4*eps*u+eps**2*n2/2)*jacobian**2
    d2_solution = s.solve(s.diff(transformed_yy,eps,2).subs(eps,0),d2)[0]
    zero('material_coordinate_acceleration_jet', d2_solution-(-n2/2+12*u**2))
    zero('material_generator_jet', (d2_solution-4*u**2)-(-n2/2+8*u**2))
    chain_metric_quadratic = 4*u*(2*uy+4*a*u)+2*a*(d2_solution+n2/2)
    zero('material_tangential_metric_normal_jet', chain_metric_quadratic-(8*u*uy+40*a*u**2))
    chain_scalar_quadratic = 4*u*vy+p*(d2_solution+n2/2)
    zero('material_scalar_normal_jet', chain_scalar_quadratic-(4*u*vy+12*p*u**2))
    wrong('confuse_coordinate_acceleration_with_generator', d2_solution-(-n2/2+8*u**2))
    ap = -1/r**2-p**2/3
    zero('graph_Israel_background_coefficient', ap-sp*p/6+1/r**2)

    # Exact degree-four harmonic decomposition in R^4 and restriction to S^3.
    xx = s.symbols('x1:5')
    rr = sum(x*x for x in xx)
    Y = xx[0]**2-xx[1]**2
    J = xx[0]**2+xx[1]**2
    Q4 = Y**2-rr*J/2+rr**2/12
    zero('F4_homogeneous_harmonic', sum(s.diff(Q4,x,2) for x in xx))
    zero('F4_Euler_degree_four', sum(x*s.diff(Q4,x) for x in xx)-4*Q4)
    zero('Y_squared_Q_coefficient_one', Y**2-(Q4+rr*(J-rr/2)/2+rr**2/6))
    zero('F4_sphere_eigenvalue', 4*(4+2)-24)

    # Complete scalar component rows and second-gauge covariance.  Christoffel
    # terms of gamma=-dt^2+cosh(t)^2 kappa are retained explicitly.
    yy, tt = s.symbols('y t', real=True)
    rho = s.Function('rho')(yy)
    pp = s.Function('p')(yy)
    aa = s.diff(rho,yy)/rho
    cc = s.cosh(tt)
    hh = s.tanh(tt)
    lam = s.Symbol('lambda', positive=True)
    sig2 = s.Symbol('sigma2')
    sig1 = -2*pp
    pb = s.diff(pp,yy)+sig2*pp/2
    field = {name: s.Function(name)(yy,tt) for name in ('N','B','C','T','J','K','G','w','Z')}

    def rows(f: dict[str,s.Expr]) -> dict[str,s.Expr]:
        N,B,C,T,J,K,G,w,Z = [f[n] for n in ('N','B','C','T','J','K','G','w','Z')]
        return {
            'tt': s.diff(T,yy)/2-s.diff(B,tt)/rho**2+aa*N/2+sig1*w/6-(s.diff(Z,tt,2)-Z)/rho**2,
            'ti': s.diff(J,yy)/2-(s.diff(C,tt)+B-2*hh*C)/(2*rho**2)-(s.diff(Z,tt)-hh*Z)/rho**2,
            'angular_TF': s.diff(G,yy)/2-(C+Z)/(rho**2*cc**2),
            'angular_trace': s.diff(K,yy)/2+lam*(C+Z)/(3*rho**2*cc**2)+hh*(B+s.diff(Z,tt))/rho**2-aa*N/2-sig1*w/6-Z/rho**2,
            'scalar': s.diff(w,yy)-pp*N/2+sig2*w/2+pb*Z,
        }

    alpha = s.Function('alpha')(yy,tt)
    beta = s.Function('beta')(yy,tt)
    ell = s.Function('ell')(yy,tt)
    gauge = {
        'N': 2*s.diff(alpha,yy),
        'B': s.diff(alpha,tt)-rho**2*s.diff(beta,yy),
        'C': alpha+rho**2*cc**2*s.diff(ell,yy),
        'T': -2*aa*alpha-2*s.diff(beta,tt),
        'J': cc**2*s.diff(ell,tt)-beta,
        'K': 2*aa*alpha+2*hh*beta-2*lam*ell/3,
        'G': 2*ell,
        'w': pp*alpha,
        'Z': -alpha,
    }
    background = {s.diff(rho,yy,2): (s.diff(rho,yy)**2-1)/rho-rho*pp**2/3}
    for name, expression in rows(gauge).items():
        zero('second_gauge_covariance_'+name, s.simplify(expression).subs(background))

    rw = dict(field)
    rw.update({n:s.Integer(0) for n in ('C','J','G','Z','B')})
    rwrows = rows(rw)
    zero('RW_tt_row', rwrows['tt']-(s.diff(field['T'],yy)/2+aa*field['N']/2+sig1*field['w']/6))
    zero('RW_trace_row', rwrows['angular_trace']-(s.diff(field['K'],yy)/2-aa*field['N']/2-sig1*field['w']/6))
    zero('RW_scalar_row', rwrows['scalar']-(s.diff(field['w'],yy)-pp*field['N']/2+sig2*field['w']/2))
    zero('RW_tt_plus_trace_neumann_sum', rwrows['tt']+rwrows['angular_trace']-s.diff(field['T']+field['K'],yy)/2)
    wrong('delete_graph_Hessian', -s.diff(alpha,tt,2)/rho**2)
    wrong('omit_shift_covariant_friction', 2*hh*field['C']/rho**2)
    wrong('unrescaled_metric_neumann', aa*field['K'])

    # Join to the independently derived bulk angular-TF algebraic relation.
    # This verifies the join algebra; it does not re-derive that bulk equation.
    Csource = 16*u**2-4*v**2/3
    Cy = s.diff(Csource,u)*uy+s.diff(Csource,v)*vy
    EN = -16*a*u**2-16*p*u*v/3+2*s2*v**2/3
    zero('angular_TF_boundary_join_source', (Cy-4*S_K).subs(first_rows)-EN)
    M = s.Matrix([[2*a,-4*p/3],[-p/2,s2/2]])
    D = s.diag(3,8)
    zero('Robin_pair_positive_weight_symmetrization', (D*M-(D*M).T).norm()**2)
    zero('Robin_pair_weighted_determinant', (D*M).det()-(24*a*s2-16*p**2))
    wrong('Robin_pair_wrong_mixing_factor', (s.diag(1,1)*M-M.T).norm()**2)

    # Independent Taylor differentiation of the exact boundary-orthogonal
    # relative-GH laws, for arbitrary transformed first tangential fields.
    n1,n2,psi1,psi2,tr1,tr2,trA2 = s.symbols('n1 n2 psi1 psi2 tr1 tr2 trA2')
    lapse = 1+eps*n1+eps**2*n2/2
    scalar = psi1*eps+psi2*eps**2/2
    tension = 6*a+sp*scalar+s2*scalar**2/2+4*scalar**3/6
    tension_prime = sp+s2*scalar+4*scalar**2/2
    m0,m1,m2 = s.symbols('h0_component h1_component h2_component')
    hm = m0+eps*m1+eps**2*m2/2
    twice = lambda expression: s.diff(expression,eps,2).subs(eps,0)
    gh_metric_rhs = (2*sp*psi1/3+4*a*n1)*m1+(s2*psi1**2/3+2*sp*psi1*n1/3)*m0
    zero('relative_GH_metric_quadratic_boundary_source', twice(tension*lapse*hm/3)-2*a*m2-(sp*psi2/3+2*a*n2)*m0-gh_metric_rhs)
    gh_scalar_rhs = -2*psi1**2-s2*psi1*n1
    zero('relative_GH_scalar_quadratic_boundary_source', twice(-tension_prime*lapse/2)+s2*psi2/2+sp*n2/2-gh_scalar_rhs)
    trace_h0 = 4-eps*tr1+eps**2*(2*trA2-tr2)/2
    exact_normal_rhs = 2*tension*lapse**2/3-a*lapse**3*trace_h0
    gh_normal_rhs = 2*s2*psi1**2/3+8*sp*psi1*n1/3-16*a*n1**2-2*a*trA2+6*a*n1*tr1
    zero('relative_GH_lapse_quadratic_boundary_source', twice(exact_normal_rhs)+4*a*n2-a*tr2-2*sp*psi2/3-gh_normal_rhs)
    zero('relative_GH_background_normal_condition', exact_normal_rhs.subs(eps,0))
    zero('metric_lapse_factorial_conversion', twice(lapse**2)-(2*n2+2*n1**2))
    zero('nonconstant_tension_Codazzi_cancellation', -3*sp/6-p)
    wrong('GH_omit_first_tension_lapse_cross', gh_scalar_rhs-(-2*psi1**2))
    # Character check of the complete finite angular representation set.
    sx,sy = s.symbols('SU2_x SU2_y', nonzero=True)
    character = lambda j,q: sum(q**(2*m) for m in range(-j,j+1))
    first_character = character(1,sx)*character(1,sy)
    symmetric_square = (first_character**2+first_character.subs({sx:sx**2,sy:sy**2},simultaneous=True))/2
    irreps = [(0,0),(1,1),(2,2),(2,0),(0,2)]
    zero('SO4_complete_quadratic_character_decomposition', symmetric_square-sum(character(j,sx)*character(k,sy) for j,k in irreps))
    zero('SO4_quadratic_representation_dimension', sum((2*j+1)*(2*k+1) for j,k in irreps)-9*10//2)

    receipt = {
        'status':'PASS',
        'positive_checks':len(checks),
        'wrong_formula_controls':len(controls),
        'checks':checks,
        'controls':controls,
        'normalization':'g=g0+eps*h1+eps^2*h2/2; phi likewise; one shell tension',
        'scope':'Exact local Taylor, harmonic and gauge/boundary identities; not a curved coupled IBVP existence theorem.',
        'sympy_version':s.__version__,
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','positive_checks':len(checks),'wrong_formula_controls':len(controls)}))


if __name__ == '__main__':
    main()
