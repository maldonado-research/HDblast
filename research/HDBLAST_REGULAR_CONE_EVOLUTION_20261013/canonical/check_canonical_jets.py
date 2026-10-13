#!/usr/bin/env python3
"""Independent exact background, regular-master and physical-jet checks."""
import argparse
import json
from pathlib import Path
import sympy as S


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--inject-failure', action='store_true')
    args = parser.parse_args()
    checks, controls = {}, {}

    def zero(name, expr):
        residual = S.factor(S.cancel(S.together(expr)))
        if residual != 0:
            raise RuntimeError(f'{name}: {residual}')
        checks[name] = 'passed'

    def nonzero(name, expr):
        residual = S.factor(expr)
        if residual == 0:
            raise RuntimeError(f'{name}: wrong formula escaped')
        controls[name] = str(residual)

    # In five ambient dimensions, Box H(s)=10H_s+4sH_ss and E H=2sH_s.
    s = S.symbols('s', real=True)
    O, f, psi = [S.Function(x)(s) for x in ('Omega', 'Phi_s', 'psi')]
    ell, j = S.diff(O,s)/O, S.diff(f,s)/f
    L = ell+2*j
    weight = 1/(O*f**2)
    zero('divergence_form_weight', S.diff(weight,s)/weight+L)
    # On either fixed sign branch of f, this signed multiplier differs from
    # the positive square-root conjugation only by a constant sign.
    multiplier = S.sqrt(O)*f
    sigma_s = S.diff(multiplier,s)/multiplier
    zero('canonical_multiplier', 2*sigma_s-L)
    box_radial = lambda x: 10*S.diff(x,s)+4*s*S.diff(x,s,2)
    op = lambda x: box_radial(x)-4*s*L*S.diff(x,s)-(8*L+S.Rational(8,3)*s*f**2)*x
    V = 3*L-2*s*S.diff(L,s)+s*L**2+S.Rational(8,3)*s*f**2
    zero('canonical_operator_all_radial_jets', op(multiplier*psi)/multiplier-box_radial(psi)+V*psi)
    # Nonradial cross-derivative coefficient cancels as well.
    zero('canonical_operator_nonradial_drift', 4*sigma_s-2*L)
    nonzero('wrong_canonical_mass_sign', 2*V)

    # Solve the lowest regular background coefficients independently.
    U, U1, U2, k, rho3, phi4 = S.symbols('U U1 U2 k rho3 phi4')
    background = {k:U1/5, rho3:-U/36,
                  phi4:(U1/5)*(U2/S.Integer(56)+U/S.Integer(126))}
    zero('radial_constraint_cubic_coefficient', (6*rho3+U/6).subs(background))
    zero('scalar_leading_coefficient', (5*k-U1).subs(background))
    zero('scalar_next_coefficient', (28*phi4-k*(U2/2-8*rho3)).subs(background))
    ell0 = S.Rational(3,2)*rho3
    j0 = (k*rho3+2*phi4)/(k/2)
    L0 = (ell0+2*j0).subs(background)
    zero('regular_master_cone_coefficient', L0-(U2/7-5*U/56))
    zero('canonical_potential_cone_coefficient', 3*L0-(3*U2/7-15*U/56))
    # Cross-check the independently found minimally coupled curved master.
    Vg0 = -3*U2/7-5*U/14
    zero('curved_master_potential_consistency', 15*ell0.subs(background)-Vg0-3*L0)
    nonzero('omit_conformal_coordinate_correction', ell0-rho3)

    # Exact jets at T=R=t0 for any real S^3 harmonic H, lambda=l(l+2).
    t, a, b, C, Lc, lam, H, f0 = S.symbols('t a b C L0 lambda H f0')
    F1tt = -(lam/t**2+8*Lc)*a*H
    F2tt = -2*Lc*t*b*H
    x1 = -4*C*a*H
    x2 = -C*t*b*H
    x1v = -C*t*F1tt/2
    x2v = -C*(5*b*H+t*F2tt)/2
    W = S.expand(x1*x2v-x2*x1v)
    zero('physical_null_Wronskian_all_harmonics', W-C**2*a*b*(10+lam/2)*H**2)
    zero('Wronskian_background_independence', S.diff(W,Lc))
    zero('l2_Wronskian', W.subs(lam,8)-14*C**2*a*b*H**2)
    Z1 = 4*f0*(2*t*x1v-x1)
    Z2 = 4*f0*(2*t*x2v-x2)
    zero('first_relational_invariant', Z1-4*f0*C*a*(lam+4+8*t**2*Lc)*H)
    zero('second_relational_invariant', Z2-4*f0*C*t*b*(-4+2*t**2*Lc)*H)
    # Choose |t^2 L0|<1. l=2 then first bracket >4, second <-2.
    zero('first_invariant_lower_endpoint', (8+4+8*S.Integer(-1))-4)
    zero('second_invariant_upper_endpoint', (-4+2*S.Integer(1))+2)
    nonzero('wrong_null_derivative_euler_factor', C**2*a*b*H**2)
    nonzero('wrong_angular_flux_coefficient', W-C**2*a*b*10*H**2)
    # Integrate a unit-normalized l=2 harmonic over S^3: per-copy density
    # v^3/8, doubled v^3/4 at v0=2t. This is bilinear, not energy.
    doubled = (2*t)**3/S.Integer(4)*14*C**2*a*b
    zero('doubled_integrated_flux_density_at_initial_cone', doubled-28*t**3*C**2*a*b)
    nonzero('single_copy_mislabeled_doubled', doubled-doubled/2)
    if args.inject_failure:
        zero('injected_failure', S.Integer(1))
    result = {'positive_count':len(checks), 'positive_checks':checks,
              'wrong_formula_count':len(controls), 'wrong_formula_controls':controls,
              'scope':'Exact algebra; no numerical evolution, no global shell theorem, no observational claim.'}
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text,end='')


if __name__ == '__main__':
    main()
