#!/usr/bin/env python3
"""Exact algebra checks for RESIDUAL_TO_PHYSICAL_RESPONSE.md.

These finite identities are not formal PDE or interval-background validation.
Use --inject-error {boundary-source,moving-flux,physical-prefactor} to test
that selected incorrect formulas are rejected, including under python -O.
"""

import argparse
import json
import sys

import sympy as s


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inject-error", choices=(
        "boundary-source", "moving-flux", "physical-prefactor"))
    args = parser.parse_args()
    checked = []

    def equal(name, actual, expected):
        difference = s.simplify(s.expand(actual - expected))
        if difference != 0:
            raise RuntimeError(f"{name}: nonzero discrepancy {difference}")
        checked.append(name)

    u, v, r, acc, mu, beta, nu, lam, F, G, J, W = s.symbols(
        "u v r acc mu beta nu lam F G J W", real=True)
    source_sign = -1 if args.inject_error == "boundary-source" else 1
    equal("boundary flux and inertia cancellation",
          -v * (mu * acc + beta * u - G) + mu * v * acc + nu * u * v,
          (nu - beta) * u * v + source_sign * G * v)

    moving_sign = 1 if args.inject_error == "moving-flux" else -1
    equal("shrinking right boundary flux",
          v * r - (v**2 + r**2 + lam*u**2)/2,
          moving_sign * ((v-r)**2 + lam*u**2)/2)

    equal("incoming shell ODE sign",
          -(J-v) + mu*acc + beta*u - G,
          mu*acc + v + beta*u - J - G)

    utz, utt, uzz = s.symbols("utz utt uzz")
    equal("incoming characteristic derivative",
          (utt+utz)-(utz+uzz), utt-uzz)
    equal("incoming characteristic source",
          (utt-uzz).subs(utt, uzz-W*u+F), F-W*u)

    sq_lam = s.symbols("sqrt_lambda", positive=True)
    equal("bulk energy cross normalization",
          v**2 + sq_lam**2*u**2 - 2*sq_lam*u*v,
          (v-sq_lam*u)**2)
    sq_mu, sq_nu = s.symbols("sqrt_mu sqrt_nu", positive=True)
    equal("shell energy cross normalization",
          sq_mu**2*v**2 + sq_nu**2*u**2 - 2*sq_mu*sq_nu*u*v,
          (sq_mu*v-sq_nu*u)**2)
    c = s.symbols("c", real=True)
    equal("shell observation Cauchy Schwarz remainder",
          (1/sq_mu**2+c**2/sq_nu**2)*(sq_mu**2*v**2+sq_nu**2*u**2)
          -(v+c*u)**2,
          (sq_nu*u/sq_mu-c*sq_mu*v/sq_nu)**2)

    wt, wz, uu, wh, mh, bh, rr, rb = s.symbols(
        "trial_tt trial_zz trial_value What muhat betahat R rb")
    exact_bulk_residual = wt-wz+W*uu
    equal("coefficient-aware bulk error source",
          -exact_bulk_residual,
          -(wt-wz+wh*uu)-(W-wh)*uu)
    exact_boundary_residual = -r+mu*wt+beta*uu
    equal("coefficient-aware boundary error source",
          -exact_boundary_residual,
          -(-r+mh*wt+bh*uu)-(mu-mh)*wt-(beta-bh)*uu)

    rho, p, H, D = s.symbols("rho p cosh_t D", positive=True)
    S = 3*D/(rho**3*p*H**s.Rational(3,2))
    equal("physical weighted norm cancellation", H**3*S**2,
          9*D**2/(rho**6*p**2))
    wrong_factor = 9 if args.inject_error == "physical-prefactor" else 3
    equal("physical coefficient normalization",
          s.pi**2*rho**4*H**3*S**2/3,
          wrong_factor*s.pi**2*D**2/(rho**2*p**2))

    t = s.symbols("t", real=True)
    q = s.Function("q")(t)
    equal("Dirichlet comparator cross term",
          2*s.exp(-t)*q*s.diff(q,t),
          s.diff(s.exp(-t)*q**2,t)+s.exp(-t)*q**2)

    x = s.symbols("x", positive=True)
    profile = s.exp(4-1/(x*(1-x)))
    equal("dimensionless bump derivative", s.diff(profile,x),
          profile*(1-2*x)/(x**2*(1-x)**2))
    y = s.symbols("y", positive=True)
    equal("bump derivative envelope monotonicity",
          s.diff(y**2*s.exp(4-y),y), y*(2-y)*s.exp(4-y))

    # For ell=x^2<=1 the explicit polynomial has only powers >=1.
    M, b, root_m = s.symbols("M b root_m", positive=True)
    T = 3*x**2
    f_bound = 4*s.sqrt(8)*M*x
    d_bound = f_bound+17/(root_m*x**2)
    expression = T*f_bound + M*d_bound*3*T**2/2 + b*d_bound*3*T**s.Rational(3,2)/s.sqrt(3)
    polynomial = (12*s.sqrt(8)*M*x**3+54*s.sqrt(8)*M**2*x**5
                  +s.Rational(459,2)*M*x**2/root_m
                  +36*s.sqrt(8)*M*b*x**4+153*b*x/root_m)
    equal("short pulse explicit bound polynomial", expression, polynomial)
    K = (12*s.sqrt(8)*M+54*s.sqrt(8)*M**2+s.Rational(459,2)*M/root_m
         +36*s.sqrt(8)*M*b+153*b/root_m)
    equal("short pulse K coefficient sum", polynomial.subs(x,1), K)

    print(json.dumps({"status":"pass", "checks":len(checked),
                      "names":checked}, sort_keys=True, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"CHECK_FAILED: {exc}", file=sys.stderr)
        sys.exit(1)
