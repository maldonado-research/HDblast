#!/usr/bin/env python3
"""Exact reconstruction/residual checks; no substitute for the PDE theorem."""
import argparse
import json
from pathlib import Path

import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--inject-failure', action='store_true')
    args = parser.parse_args()
    checks, controls = {}, {}

    def zero(name, expr):
        value = s.simplify(s.expand(expr))
        if value != 0:
            raise RuntimeError(f'{name}: {value}')
        checks[name] = 'passed'

    def reject(name, expr):
        value = s.simplify(expr)
        if value == 0:
            raise RuntimeError(f'{name}: corruption was not detected')
        controls[name] = str(value)

    # Full five-dimensional Christoffel variation, no harmonic ansatz.
    y, t, x1, x2, x3 = coords = s.symbols('y t x1 x2 x3', real=True)
    rho = s.Function('rho')(y)
    psi = s.Function('Psi')(*coords)
    chi = s.Function('chi')(*coords)
    phi = s.Function('phi')(y)
    a = s.diff(rho, y)/rho
    gam = [-1, s.exp(2*t), s.exp(2*t), s.exp(2*t)]
    g = [s.Integer(1)] + [rho**2*v for v in gam]
    h = [-4*psi] + [2*psi*v for v in g[1:]]
    Gamma, dGamma = {}, {}
    for A in range(5):
        for B in range(5):
            for C in range(5):
                dg = ((s.diff(g[A], coords[B]) if A == C else 0)
                      + (s.diff(g[A], coords[C]) if A == B else 0)
                      - (s.diff(g[B], coords[A]) if B == C else 0))
                dh = ((s.diff(h[A], coords[B]) if A == C else 0)
                      + (s.diff(h[A], coords[C]) if A == B else 0)
                      - (s.diff(h[B], coords[A]) if B == C else 0))
                Gamma[A,B,C] = s.simplify(dg/(2*g[A]))
                dGamma[A,B,C] = s.simplify(dh/(2*g[A])
                                           -h[A]*Gamma[A,B,C]/g[A])

    def box_gamma(f):
        return -s.diff(f,t,2)-3*s.diff(f,t)+s.exp(-2*t)*sum(
            s.diff(f,x,2) for x in (x1,x2,x3))

    radial_ric = -4*s.diff(psi,y,2)-16*a*s.diff(psi,y)+2*box_gamma(psi)/rho**2
    slice_ric = -box_gamma(psi)-rho**2*(s.diff(psi,y,2)+10*a*s.diff(psi,y)
                  +6*(s.diff(a,y)+4*a**2)*psi)
    for A in range(5):
        for B in range(A,5):
            val = sum(s.diff(dGamma[C,A,B],coords[C])
                      -s.diff(dGamma[C,A,C],coords[B]) for C in range(5))
            for C in range(5):
                for D in range(5):
                    val += (dGamma[C,C,D]*Gamma[D,A,B]
                            +Gamma[C,C,D]*dGamma[D,A,B]
                            -dGamma[C,B,D]*Gamma[D,A,C]
                            -Gamma[C,B,D]*dGamma[D,A,C])
            if A == B == 0:
                target = radial_ric
            elif A == 0:
                target = -3*s.diff(s.diff(psi,y)+2*a*psi,coords[B])
            elif A == B:
                target = gam[A-1]*slice_ric
            else:
                target = 0
            zero(f'direct_Ricci_{A}{B}',val-target)

    p = s.diff(phi,y)
    scalar = 0
    for A in range(5):
        scalar += s.diff(chi,coords[A],2)/g[A]
        scalar -= sum(Gamma[D,A,A]*s.diff(chi,coords[D]) for D in range(5))/g[A]
        background_hessian = s.diff(phi,coords[A],2)-sum(
            Gamma[D,A,A]*s.diff(phi,coords[D]) for D in range(5))
        scalar -= h[A]*background_hessian/g[A]**2
        scalar -= sum(dGamma[D,A,A]*s.diff(phi,coords[D]) for D in range(5))/g[A]
    scalar_target = (s.diff(chi,y,2)+4*a*s.diff(chi,y)+box_gamma(chi)/rho**2
                     +6*p*s.diff(psi,y)+4*(s.diff(p,y)+4*a*p)*psi)
    zero('direct_scalar_wave_variation',scalar-scalar_target)

    # Background identities and arbitrary (not separated) tangential boxes.
    A, H, p, py, pyy, U, U1, U2 = s.symbols('a H p py pyy U U1 U2')
    P, Py, Pyy, X, Xy, Xyy, Bpsi, BX = s.symbols('Psi Psiy Psiyy X Xy Xyy BoxPsi BoxX')
    ay = -H-p**2/3
    mom = -2*A*P-p*X/3
    mom2 = -2*ay*P-2*A*mom-py*X/3-p*Xy/3
    bg = {U: 6*H+p**2/2-6*A**2, U1: py+4*A*p}
    C = p*Xy-py*X-3*H*(Bpsi+4*P)+2*p**2*P
    Eyy = -4*mom2-16*A*mom+2*H*Bpsi-2*p*Xy-2*U1*X/3+8*U*P/3
    Etrace = (-H*Bpsi-mom2-10*A*mom-6*(ay+4*A**2)*P
              -4*U*P/3-2*U1*X/3)
    zero('Einstein_radial_factor',(Eyy+2*C/3).subs(bg))
    zero('Einstein_tangential_factor',(Etrace-C/3).subs(bg))
    # Differentiate C, commuting the tangential box with y derivatives.
    Cprime = (p*Xyy-pyy*X+6*A*H*(Bpsi+4*P)
              -3*H*(-2*A*(Bpsi+4*P)-p*(BX+4*X)/3)
              +4*p*py*P+2*p**2*mom)
    S = Xyy+4*A*Xy+H*BX-U2*X+6*p*mom+4*U1*P
    zero('scalar_residual_factor',
         (p*S-Cprime-4*A*C).subs({pyy:U2*p-4*ay*p-4*A*py, U1:py+4*A*p}))
    reject('omit_Einstein_scalar_cross_term',2*p*Xy)
    reject('wrong_tangential_residual_sign',2*C/3)

    # Scalar master M[P] and its regular ambient equation.
    r = s.symbols('r', positive=True)
    O = s.Function('Omega')(r)
    K = s.Function('Phi_s')(r)
    F = s.Function('F')(r)
    Q = s.Function('wave_residual')(r)
    boxF = s.symbols('Box_gamma_F')
    dy = lambda value: s.diff(value,r)/O
    rho = r*O
    p = 2*r*K/O
    aa = dy(rho)/rho
    L = (s.diff(O,r)/O+2*s.diff(K,r)/K)/(2*r)
    masterP = r**4*F
    M = (dy(dy(masterP))-(2*aa+2*dy(p)/p)*dy(masterP)
         +(r**4*boxF+4*masterP)/rho**2-2*p**2*masterP/3)
    waveF = (s.diff(F,r,2)+4*s.diff(F,r)/r+boxF/r**2
             -2*r*L*s.diff(F,r)-(8*L+8*r**2*K**2/3)*F)
    zero('exact_ambient_master',M-r**4*waveF/O**2)
    reconstruction = -3*dy(masterP)/(rho**2*p)
    zero('exact_ambient_scalar_reconstruction',reconstruction+3*(r*s.diff(F,r)+4*F)/(2*O**2*K))
    psi_ambient = masterP/rho**2
    zero('exact_ambient_metric_eta_coefficient',2*O**2*psi_ambient-2*r**2*F)
    zero('exact_ambient_metric_XX_coefficient',-6*O**2*psi_ambient/r**2+6*F)
    Cdirect = (p*dy(reconstruction)-dy(p)*reconstruction
               -3*(r**4*boxF+4*masterP)/rho**4+2*p**2*psi_ambient)
    zero('radial_constraint_master_relation',Cdirect+3*M/rho**2)
    scalar_factor = -3*(dy(r**4*Q/O**2)+2*aa*r**4*Q/O**2)/(rho**2*p)
    zero('ambient_scalar_residual',scalar_factor+3*(r*s.diff(Q,r)+6*Q)/(2*O**4*K))
    # In E_AB the radial/tangential decomposition is (M/rho^2)(3dy^2-g).
    metric_radial_coeff = (r**4*Q/O**2)/rho**2*3*O**2/r**2
    metric_eta_coeff = -(r**4*Q/O**2)/rho**2*O**2
    zero('ambient_Einstein_XX_coefficient',metric_radial_coeff-3*Q/O**2)
    zero('ambient_Einstein_eta_coefficient',metric_eta_coeff+r**2*Q/O**2)
    # Corrupt P=r^4 F to P=r^3 F. In the algebraic test O=K=F=1,
    # the supposed regular wave equation misses a genuine -2r residual.
    badP = r**3
    badM = s.diff(badP,r,2)-4*s.diff(badP,r)/r+4*badP/r**2-8*r**2*badP/3
    reject('wrong_master_radial_power',(badM-badP*(-8*r**2/3)).subs(r,1))
    reject('omit_master_drift',2*r*L*s.diff(F,r))

    # Relational scalar Z=delta K_background-(K_s/Phi_s) delta phi.
    z = s.symbols('s')
    O = s.Function('Omega')(z)
    K = s.Function('Phi_s')(z)
    f, ef, eef = s.symbols('F E_F E2_F')
    L = s.diff(O,z)/O+2*s.diff(K,z)/K
    coeff = -3/(2*O**2*K)
    chi = coeff*(ef+4*f)
    echi = 2*z*s.diff(coeff,z)*(ef+4*f)+coeff*(eef+4*ef)
    bgK = 4*z*K**2/O**2
    Z = 16*z**2*K**2*f/O**4+4*K*echi/O**2-s.diff(bgK,z)*chi/K
    Ztarget = -6/O**4*(eef+(3-2*z*L)*ef-(4+8*z*L)*f-8*z**2*K**2*f/3)
    zero('relational_scalar_invariant',Z-Ztarget)

    # Cauchy jets proving a physical pair with positive small-cone Wronskian.
    time, v, d, e, L0, K0, lam = s.symbols('T v d e L0 Phi_s0 lambda')
    spatial = s.symbols('X1:5')
    R2 = sum(q*q for q in spatial)
    cubic = (4*e-5*L0*d)/3
    F2 = d*time+e*time*R2+cubic*time**3
    box = -s.diff(F2,time,2)+sum(s.diff(F2,q,2) for q in spatial)
    Euler = time*s.diff(F2,time)+sum(q*s.diff(F2,q) for q in spatial)
    linear_residual = s.expand(box-2*L0*Euler-8*L0*F2)
    # Its homogeneous degree-one part is the only condition on the cubic jet.
    scale = s.symbols('scale')
    scaled = linear_residual.subs({time:scale*time, **{q:scale*q for q in spatial}}, simultaneous=True)
    zero('Cauchy_cubic_master_jet',s.expand(scaled).coeff(scale,1))
    F2cone = d*v/2+(7*e-5*L0*d)*v**3/24
    chi2 = -3*(v*s.diff(F2cone,v)+4*F2cone)/(2*K0)
    # Simultaneous substitution then d replacement avoids an ordering convention.
    normalized = chi2.subs(e,(5*L0*d+lam)/7).subs(d,-2*K0/15)
    zero('normalized_cone_scalar_jet',normalized-(v/2-7*lam*v**3/(16*K0)))
    Z2 = 4*K0*(v*s.diff(normalized,v)-normalized)
    zero('physical_second_solution_witness',Z2+7*lam*v**3/2)
    c2 = s.symbols('chi1_quadratic')
    chi1 = 1+c2*v**2
    W = s.expand(chi1*s.diff(normalized,v)-normalized*s.diff(chi1,v))
    zero('nonzero_flux_leading_wronskian',W.subs(v,0)-s.Rational(1,2))
    reject('both_solutions_are_pure_gauge',Z2.subs({lam:1,v:1}))
    reject('zero_local_characteristic_flux',W.subs(v,0))

    if args.inject_failure:
        zero('injected_failure',s.Integer(1))
    result = {'scope':'Exact local differential identities and jets; PDE existence uses the stated normally hyperbolic Cauchy theorem, and no shell solution is checked.',
              'positive_count':len(checks),'positive_checks':checks,
              'wrong_formula_count':len(controls),'wrong_formula_controls':controls}
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text,end='')


if __name__ == '__main__':
    main()
