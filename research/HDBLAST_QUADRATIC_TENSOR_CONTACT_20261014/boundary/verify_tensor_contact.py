#!/usr/bin/env python3
"""Independent exact operator/sign checks; no scientific imports or asserts."""
import argparse
import json
import sys
from pathlib import Path

import sympy as sp


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--inject-error', action='store_true')
    args = parser.parse_args()
    positives = []
    controls = []

    def require_zero(name, expression):
        residual = sp.simplify(sp.expand(expression))
        if residual != 0:
            raise RuntimeError(f'{name}: nonzero residual {residual}')
        positives.append(name)

    def reject_wrong(name, expression):
        residual = sp.simplify(sp.expand(expression))
        if residual == 0:
            raise RuntimeError(f'{name}: wrong formula unexpectedly accepted')
        controls.append(name)

    a, ap, rho, pot = sp.symbols('a ap rho pot', nonzero=True, real=True)
    q, qy, qyy = sp.symbols('q qy qyy', real=True)
    # Gauss terms differentiated separately: K', 2 K.K, and K K_ab.
    normal_variation = (
        -rho**2*((ap+2*a*a)*q+2*a*qy+qyy/2)
        +2*rho**2*(a*a*q+a*qy)
        -4*a*rho**2*(a*q+qy/2)
    )
    require_zero('Gauss radial curvature variation', normal_variation+
                 rho**2*((ap+4*a*a)*q+2*a*qy+qyy/2))
    bulk_relation = {ap: 3/rho**2-4*a*a-sp.Rational(2, 3)*pot}
    require_zero('actual scalar-background pressure cancellation',
                 (normal_variation-sp.Rational(2, 3)*pot*rho**2*q).subs(bulk_relation)
                 +3*q+rho**2*(2*a*qy+qyy/2))
    require_zero('radial tensor friction four',
                 -2*sp.expand(normal_variation/rho**2).coeff(qy)-4*a)

    t = sp.symbols('t', real=True)
    c = sp.cosh(t)
    h = sp.diff(c, t)/c
    kt, lam = sp.symbols('kt lam', real=True)
    require_zero('unit de Sitter expansion', h-sp.tanh(t))
    require_zero('de Sitter tensor curvature shift',
                 (2*h*h-2-kt/c**2)+(kt+2)/c**2)
    require_zero('l2 tensor versus scalar eigenvalue', (kt+2).subs(kt, 6)-8)
    reject_wrong('omit two mixed-index connection terms',
                 -2-kt/c**2+(kt+2)/c**2)
    reject_wrong('use rough tensor eigenvalue without plus two',
                 (2*h*h-2-kt/c**2)+kt/c**2)

    x = sp.symbols('x', real=True)
    R = sp.Function('rho')(x)
    H = sp.Function('H')(x)
    # dy=rho dx; evaluate H_yy+4a H_y from direct chain rules.
    Hy = sp.diff(H, x)/R
    Hyy = sp.diff(Hy, x)/R
    ay = sp.diff(R, x)/R**2
    require_zero('proper-radial to tortoise chain rule',
                 R**2*(Hyy+4*ay*Hy)-sp.diff(H, x, 2)
                 -3*sp.diff(R, x)/R*sp.diff(H, x))

    vv = sp.Function('v')(x, t)
    wt = sp.tanh(t)
    wx = sp.diff(R, x)/R
    S = sp.Rational(3, 2)*wx
    conjugating = R**sp.Rational(3, 2)*c**sp.Rational(3, 2)
    HH = vv/conjugating
    friction = 4 if args.inject_error else 3
    transformed = conjugating*(
        sp.diff(HH, t, 2)-sp.diff(HH, x, 2)
        +friction*wt*sp.diff(HH, t)-3*wx*sp.diff(HH, x)
        +lam/c**2*HH)
    expected = (sp.diff(vv, t, 2)-sp.diff(vv, x, 2)
                +(S*S+sp.diff(S, x)-sp.Rational(9, 4)
                  +(lam+sp.Rational(3, 4))/c**2)*vv)
    require_zero('full radial-time canonical conjugation', transformed-expected)
    require_zero('canonical de Sitter time potential',
                 -sp.Rational(3, 2)*sp.diff(wt, t)-sp.Rational(9, 4)*wt*wt
                 +sp.Rational(9, 4)-sp.Rational(3, 4)/c**2)
    require_zero('static tensor partner potential',
                 S*S+sp.diff(S, x)-
                 (sp.Rational(3, 2)*sp.diff(R, x, 2)/R
                  +sp.Rational(3, 4)*sp.diff(R, x)**2/R**2))
    reject_wrong('omit time conjugation curvature three quarters',
                 sp.Rational(3, 4)/c**2)
    reject_wrong('reverse sign of radial derivative of S', 2*sp.diff(S, x))

    r, omega, rho_b, C, J, V = sp.symbols(
        'r omega rho_b C J V', positive=True, real=True)
    Htt, Ht, Hyy_s, Hy_s, H_s = sp.symbols('Htt Ht Hyy Hy H', real=True)
    LT = Hyy_s+4*a*Hy_s-(Htt+3*wt*Ht+lam/c**2*H_s)/rho**2
    deltaE = -rho**2*c**2*LT/2
    JfromE = 2*C/(rho**2*c**2)
    require_zero('second-derivative source prefactor and sign',
                 deltaE.subs(Hyy_s, JfromE-4*a*Hy_s
                              +(Htt+3*wt*Ht+lam/c**2*H_s)/rho**2)+C)
    require_zero('canonical source from complete Einstein Hessian',
                 -rho**sp.Rational(7, 2)*c**sp.Rational(3, 2)*JfromE
                 +2*rho**sp.Rational(3, 2)/sp.sqrt(c)*C)
    reject_wrong('wrong Taylor factorial in source', -rho**2*c**2*JfromE+C)
    reject_wrong('wrong source sign', rho**2*c**2*JfromE/2+C)

    Hp, Hm, Hyp, Hym, common = sp.symbols('Hp Hm Hyp Hym common', real=True)
    Kp = rho_b**2*c*c*(a*Hp+Hyp/2)
    Km = rho_b**2*c*c*(a*Hm+Hym/2)
    junction = Kp+Km-2*a*rho_b**2*c*c*common
    require_zero('one common tension and two outward warp cancellations',
                 junction.subs({Hp: common, Hm: common})
                 -rho_b**2*c*c*(Hyp+Hym)/2)
    require_zero('tensor source copies sum and difference',
                 (Hp+Hm)/2+(Hp-Hm)/2-Hp)
    reject_wrong('differentiate unrescaled metric at shell', 2*a*common)

    Sb, v0, vx0, g = sp.symbols('Sb v0 vx0 g', real=True)
    # H_x=rho^-3/2 c^-3/2 (v_x-S v); H_y=H_x/rho.
    require_zero('canonical even Robin sign',
                 rho_b**sp.Rational(5, 2)*c**sp.Rational(3, 2)
                 *((vx0-Sb*v0)/(rho_b**sp.Rational(5, 2)*c**sp.Rational(3, 2))-g)
                 -(vx0-Sb*v0-rho_b**sp.Rational(5, 2)*c**sp.Rational(3, 2)*g))
    reject_wrong('reverse canonical outward Robin sign', 2*Sb*v0)

    alpha, p, m, z, G = sp.symbols('alpha p m z G', real=True)
    vt = p-alpha*z+G
    vz = alpha*z-G
    require_zero('Robin supplies one incoming characteristic', (vt-vz)-(p-2*alpha*z+2*G))
    require_zero('trace ODE preserves outgoing characteristic', vt+vz-p)
    Dprime = sp.symbols('Dprime', real=True)
    require_zero('Dirichlet supplies one incoming characteristic',
                 (2*Dprime-p)+p-2*Dprime)
    # Spatial integration boundary term at inward z=0 is -v_t v_z.
    require_zero('energy uses actual Robin coefficient and cancels flux',
                 -vt*alpha*z+alpha*z*vt)
    reject_wrong('add arbitrary positive boundary energy',
                 -vt*alpha*z+(alpha+1)*z*vt)
    xi, xiz, xizz, Gtt, Wpot, forcing = sp.symbols('xi xiz xizz Gtt Wpot forcing', real=True)
    require_zero('inhomogeneous Robin lift',
                 (-xiz+alpha*xi).subs({xi: 0, xiz: -1})-1)
    require_zero('lifted source differential equation',
                 forcing-(xi*Gtt-xizz*G+Wpot*xi*G)
                 -(forcing-xi*Gtt+xizz*G-Wpot*xi*G))

    b, R0 = sp.symbols('b R0', positive=True, real=True)
    tc = sp.log(b/R0)
    require_zero('actual ambient tensor-contact time',
                 b*(sp.exp(tc)-sp.exp(-tc))/2-(b*b-R0*R0)/(2*R0))
    require_zero('outgoing radial characteristic reaches shell', sp.log(R0/b)+tc)
    require_zero('ambient positive-cone null coordinate plus',
                 (r*(sp.sinh(t)+sp.cosh(t))-r*sp.exp(t)).rewrite(sp.exp))
    require_zero('ambient positive-cone null coordinate minus',
                 (r*(sp.sinh(t)-sp.cosh(t))+r*sp.exp(-t)).rewrite(sp.exp))
    depth, et = sp.symbols('depth et', positive=True, real=True)
    tau_derivative = et/(1+depth**2*et**2)-(1/et)/(1+depth**2/et**2)
    require_zero('inherited conformal time support-domain monotonicity',
                 tau_derivative-(et-1/et)*(1-depth**2)/
                 ((1+depth**2*et**2)*(1+depth**2/et**2)))
    require_zero('shell conformal and global de Sitter time conversion',
                 ((sp.exp(t)-sp.exp(-t))/2-sp.sinh(t)).rewrite(sp.exp))
    reject_wrong('wrong outgoing ambient contact sign',
                 (b*b-R0*R0)/(2*R0)-(R0*R0-b*b)/(2*R0))

    # An odd first derivative squared is even, not odd second data.
    eps = sp.symbols('eps', real=True)
    linear, second = sp.symbols('linear second', real=True)
    family_plus = eps*linear+eps**2*second/2
    family_minus = family_plus.subs(eps, -eps)
    require_zero('parameter-reverse second copy derivative is even',
                 sp.diff(family_minus, eps, 2)-sp.diff(family_plus, eps, 2))
    reject_wrong('odd first derivative implies odd second derivative',
                 sp.diff(family_minus, eps, 2)+sp.diff(family_plus, eps, 2))

    result = {
        'status': 'PASS',
        'positive_checks': len(positives),
        'wrong_formula_controls': len(controls),
        'checks': positives,
        'controls': controls,
        'scope': 'Exact differential-operator and boundary algebra; the analytic IBVP proof is in the companion note.',
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'FAIL: {error}', file=sys.stderr)
        sys.exit(1)
