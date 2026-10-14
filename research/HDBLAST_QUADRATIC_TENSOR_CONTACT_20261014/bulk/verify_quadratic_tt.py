#!/usr/bin/env python3
"""Independent exact checks for the bulk quadratic tensor source.

Only SymPy is imported; no inherited scientific implementation is used.
The direct Ricci test uses the complete five-dimensional diagonal metric
g(eps)=diag(1-4 eps Psi, rho^2(1+2 eps Psi) gamma), with a flat dS chart.
The eps^2 coefficient, not the second derivative, is called R2.
"""
import argparse
import json
from pathlib import Path
import sympy as s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    ap.add_argument('--inject-failure', action='store_true')
    args = ap.parse_args()
    checks = {}
    controls = {}

    def zero(name, value):
        value = s.simplify(s.expand(value))
        if value != 0:
            raise RuntimeError(f'{name}: {value}')
        checks[name] = 'passed'

    def nonzero(name, value):
        value = s.simplify(s.expand(value))
        if value == 0:
            raise RuntimeError(f'{name}: mutation was not detected')
        controls[name] = str(value)

    y, t, x, z, w = coords = s.symbols('y t x z w', real=True)
    rho = s.Function('rho')(y)
    psi = s.Function('Psi')(*coords)
    gam = [-s.Integer(1), s.exp(2*t), s.exp(2*t), s.exp(2*t)]
    g = [s.Integer(1)] + [rho**2*v for v in gam]
    q = [-4*psi] + [2*psi]*4
    h = [g[i]*q[i] for i in range(5)]
    G = [{}, {}, {}]
    for A in range(5):
        for B in range(5):
            for C in range(5):
                dg = ((s.diff(g[A],coords[B]) if A == C else 0)
                      +(s.diff(g[A],coords[C]) if A == B else 0)
                      -(s.diff(g[B],coords[A]) if B == C else 0))
                dh = ((s.diff(h[A],coords[B]) if A == C else 0)
                      +(s.diff(h[A],coords[C]) if A == B else 0)
                      -(s.diff(h[B],coords[A]) if B == C else 0))
                G[0][A,B,C] = s.simplify(dg/(2*g[A]))
                G[1][A,B,C] = s.simplify(dh/(2*g[A])-q[A]*G[0][A,B,C])
                G[2][A,B,C] = s.simplify(q[A]**2*G[0][A,B,C]-q[A]*dh/(2*g[A]))

    def ric2(A,B):
        ans = sum(s.diff(G[2][C,A,B],coords[C])
                  -s.diff(G[2][C,A,C],coords[B]) for C in range(5))
        for C in range(5):
            for D in range(5):
                for k in range(3):
                    ans += (G[k][C,C,D]*G[2-k][D,A,B]
                            -G[k][C,B,D]*G[2-k][D,A,C])
        return s.expand(ans)

    # Every tangential component modulo gamma trace, not just one witness.
    diffdiag = []
    for A in range(1,5):
        for B in range(A,5):
            Hess = (s.diff(psi**2,coords[A],coords[B])
                    -sum(G[0][C,A,B]*s.diff(psi**2,coords[C]) for C in range(1,5)))
            target = 6*Hess-6*s.diff(psi,coords[A])*s.diff(psi,coords[B])
            difference = s.simplify(ric2(A,B)-target)
            if A == B:
                diffdiag.append(difference/gam[A-1])
            else:
                zero(f'direct_Ricci2_offdiag_{A}{B}',difference)
    for j in range(1,4):
        zero(f'direct_Ricci2_trace_difference_{j}',diffdiag[j]-diffdiag[0])

    # Independent radial-lapse/conformal calculation in d tangential dimensions.
    e, P, Pi, Pj, Pij, d = s.symbols('e P Pi Pj Pij d')
    omega = s.log(1+2*e*P)/2
    nu = s.log(1-2*(d-2)*e*P)/2
    D2 = lambda f: s.diff(f,P)*Pij+s.diff(f,P,2)*Pi*Pj
    anis = (-(d-2)*D2(omega)-D2(nu)
            +(d-2)*s.diff(omega,P)**2*Pi*Pj
            -s.diff(nu,P)**2*Pi*Pj
            +2*s.diff(omega,P)*s.diff(nu,P)*Pi*Pj)
    a2 = s.expand(s.diff(anis,e,2).subs(e,0)/2)
    zero('general_dimension_anisotropy',a2-(d-1)*(d-2)*(2*P*Pij+Pi*Pj))
    zero('dimension_four_anisotropy',a2.subs(d,4)-(12*P*Pij+6*Pi*Pj))

    # Angular TT projection kills Hess(Psi^2), and has gradient coefficient -4/3.
    u, v = s.symbols('u v',real=True)
    R2tt = (-6*u**2)*(-s.Rational(4,3))
    E2tt = R2tt-v**2*(-s.Rational(4,3))
    zero('full_quadratic_residual_T_coefficient',E2tt-s.Rational(4,3)*(6*u**2+v**2))
    zero('second_derivative_source_T_coefficient',-2*E2tt+s.Rational(8,3)*(6*u**2+v**2))
    nonzero('delete_geometric_source',R2tt)
    nonzero('reverse_matter_sign',-2*v**2*(-s.Rational(4,3)))
    nonzero('miss_second_derivative_factor',E2tt)

    # Exact change of variables for the principal tensor operator.
    # T=r sinh(tau), R=r cosh(tau), Omega=Omega(r).
    r, tau = s.symbols('r tau',positive=True,real=True)
    Om = s.Function('Omega')(r)
    Z = s.Function('Z')(r,tau)
    DT = lambda f: -s.sinh(tau)*s.diff(f,r)+s.cosh(tau)/r*s.diff(f,tau)
    DR = lambda f: s.cosh(tau)*s.diff(f,r)-s.sinh(tau)/r*s.diff(f,tau)
    E = lambda f: r*s.diff(f,r)
    ambient = (-DT(DT(Z))+DR(DR(Z))+3/(r*s.cosh(tau))*DR(Z)
               +3*s.diff(Om,r)/(r*Om)*E(Z)-8/(r**2*s.cosh(tau)**2)*Z)
    radial = (s.diff(Z,r,2)+(4/r+3*s.diff(Om,r)/Om)*s.diff(Z,r)
              -(s.diff(Z,tau,2)+3*s.tanh(tau)*s.diff(Z,tau)
                +8*Z/s.cosh(tau)**2)/r**2)
    zero('ambient_warped_tensor_operator',s.trigsimp(ambient-radial))

    # Gaussian radial conversion d_y=Omega^-1 d_r, rho=r Omega.
    dy = lambda f: s.diff(f,r)/Om
    a = dy(r*Om)/(r*Om)
    zero('radial_warp_factor',Om**2*(dy(dy(Z))+4*a*dy(Z))
         -s.diff(Z,r,2)-(4/r+3*s.diff(Om,r)/Om)*s.diff(Z,r))
    nonzero('wrong_tensor_angular_eigenvalue',2*Z/(r**2*s.cosh(tau)**2))
    nonzero('wrong_conformal_drift',3*s.diff(Om,r)/Om*s.diff(Z,r))

    if args.inject_failure:
        zero('injected_matter_only_claim',R2tt)
    result = {'status':'passed','convention':'g=g0+eps*h1+eps^2*k2/2',
              'checks':checks,'negative_controls':controls,
              'limits':['No shell junction condition or IBVP estimate is checked.',
                        'S3 TT decomposition is independently documented in the derivation; the coefficient is supplied to this source check.',
                        'The tensor operator is derived geometrically in the note; this script checks its coordinate conversion.']}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
