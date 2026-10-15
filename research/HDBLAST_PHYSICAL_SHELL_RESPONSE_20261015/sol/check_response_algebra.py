#!/usr/bin/env python3
"""Exact local identities for the coupled shell-response note.

This program does not verify any PDE existence theorem or its hypotheses.
No prior scientific checker is imported. Assertions are deliberately avoided.
"""
import argparse
import json
import sys

import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inject-error", choices=["potential", "shell-sign", "ambient-velocity"])
    args = parser.parse_args()
    names = []

    def check(name, expression):
        value = s.simplify(s.expand(expression))
        if value != 0:
            raise RuntimeError(f"{name}: nonzero residual {value}")
        names.append(name)

    t, z = s.symbols("t z", real=True)
    lam, kap, cb, db = s.symbols("lambda kappa Cb Db", real=True)
    C = s.Function("C")(z)
    D = s.Function("D")(z)
    alpha = s.Function("alpha")(z)
    w = s.Function("w")(t, z)
    f = s.cosh(t) ** (-s.Rational(3, 2))
    P = alpha * f * w
    original = (s.diff(P, t, 2)-s.diff(P, z, 2)
                +3*s.tanh(t)*s.diff(P, t)-C*s.diff(P, z)
                +(lam/s.cosh(t)**2-4+D)*P)
    a_rules = {s.diff(alpha, z, 2): (C*C/4-s.diff(C, z)/2)*alpha,
               s.diff(alpha, z): -C*alpha/2}
    transformed = s.expand(original).subs(a_rules) / (alpha*f)
    time_shift = s.Rational(3, 4)
    if args.inject_error == "potential":
        time_shift = s.Rational(1, 4)
    potential = ((lam+time_shift)/s.cosh(t)**2-s.Rational(25, 4)
                 +s.diff(C, z)/2+C*C/4+D)
    check("full_drift_conjugation", transformed-s.diff(w,t,2)+s.diff(w,z,2)-potential*w)
    check("l2_potential_coefficient", (lam+s.Rational(3,4)).subs(lam,8)-s.Rational(35,4))

    beta = s.Function("beta")(t)
    wn = s.Function("wn")(t)
    J = s.Function("J")(t)
    boundary_original = (s.diff(P,t,2)+3*s.tanh(t)*s.diff(P,t)
                         +(lam/s.cosh(t)**2-4)*P-kap*s.diff(P,z))
    boundary = s.expand(boundary_original).subs(a_rules)/(alpha*f)
    trace_rules = {w: beta, s.diff(w,t): s.diff(beta,t),
                   s.diff(w,t,2): s.diff(beta,t,2), s.diff(w,z): wn,
                   C: cb}
    H = (lam+s.Rational(3,4))/s.cosh(t)**2-s.Rational(25,4)+kap*cb/2
    target = s.diff(beta,t,2)-kap*wn+H*beta
    if args.inject_error == "shell-sign":
        target = s.diff(beta,t,2)+kap*wn+H*beta
    check("shell_conjugation_sign_and_mass", boundary.subs(trace_rules)-target)
    check("shell_characteristic_ODE", target.subs(wn,J-s.diff(beta,t))
          -(s.diff(beta,t,2)+kap*s.diff(beta,t)+H*beta-kap*J))
    check("scalar_transfer_characteristic", (wn-cb*beta/2).subs(wn,J-s.diff(beta,t))
          -(J-s.diff(beta,t)-cb*beta/2))

    rho, a, p, py = s.symbols("rho a p py", nonzero=True)
    log_alpha_z = s.Rational(3,2)*(-rho*a)-rho*py/p
    check("closed_radial_conjugation", log_alpha_z+rho*(3*a+2*py/p)/2)
    Pz, Pzz, Pt, Ptt = s.symbols("Pz Pzz Pt Ptt")
    box_P = -Ptt-3*s.tanh(t)*Pt-lam/s.cosh(t)**2*s.Symbol("P")
    normal_P = -Pz/rho
    check("inward_shell_normal_sign", -rho*normal_P-Pz)
    check("physical_scalar_normalization", -3*normal_P/(rho*rho*p)-3*Pz/(rho**3*p))

    r, nu, G = s.symbols("r nu G", positive=True)
    ambient_FT = nu*G*r*r
    tangent_factor = r
    if args.inject_error == "ambient-velocity":
        tangent_factor = s.Integer(1)
    global_Pt = r**4*tangent_factor*ambient_FT
    check("ambient_to_global_initial_velocity", global_Pt-nu*r*(r**6*G))

    d, ell = s.symbols("d ell", positive=True)
    bump_log = 4-ell**2/((z-d)*(d+ell-z))
    derivative_numerator = ell**2*(2*d+ell-2*z)
    check("bump_outer_half_derivative", s.diff(bump_log,z)
          -derivative_numerator/((z-d)**2*(d+ell-z)**2))
    check("bump_peak_normalization", bump_log.subs(z,d+ell/2))
    eps, M, c0, tt = s.symbols("eps M c0 tt", positive=True)
    aa = M*tt*eps/2
    Q = s.symbols("Q", positive=True)
    Sup_bound = Q*(1+c0*eps/2)/(1-aa)
    lower = Q/2-aa*Sup_bound
    check("Duhamel_lower_bound_algebra", lower-Q*(1-aa*(3+c0*eps))/(2*(1-aa)))
    ss = s.symbols("ss", nonnegative=True)
    total_area = s.integrate(eps,(ss,0,tt-eps/2))+s.integrate(2*(tt-ss),(ss,tt-eps/2,tt))
    check("near_front_source_area", total_area-(tt*eps-eps**2/4))

    # Moments of the uniform probability measure on S^3 are exact.
    n14, n12n22 = s.Rational(1,8), s.Rational(1,24)
    Y2mean = 2*n14-2*n12n22
    check("quadrupole_L2_mean", Y2mean-s.Rational(1,6))
    check("quadrupole_angular_integral", 2*s.pi**2*Y2mean-s.pi**2/3)
    amplitude = 3*Pz/(rho**3*p)
    Cstar = s.pi**2/3*rho**3*s.cosh(t)**3*amplitude**2
    check("physical_measure_response_coefficient", Cstar-3*s.pi**2*s.cosh(t)**3*Pz**2/(rho**3*p**2))
    check("Q_second_variation_factor", 2*Cstar-6*s.pi**2*s.cosh(t)**3*Pz**2/(rho**3*p**2))

    print(json.dumps({"status":"pass", "positive_identity_count":len(names), "identities":names}, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(json.dumps({"status":"fail", "error":str(exc)},sort_keys=True),file=sys.stderr)
        sys.exit(1)
