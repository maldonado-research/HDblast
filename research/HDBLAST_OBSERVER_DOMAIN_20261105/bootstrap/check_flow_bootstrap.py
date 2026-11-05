#!/usr/bin/env python3
"""Exact finite algebra checks; no evaluated physical observer inputs.

Both stdout and --output receive identical canonical JSON. Existing files
and symlink parents are rejected. No assert is used, including under -O.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + message)


def add(a, b):
    return [a[i] + b[i] for i in range(3)]


def mul(a, b):
    return [sum((a[j] * b[i-j] for j in range(i+1)), Q(0))
            for i in range(3)]


def power(a, n):
    out = [Q(1), Q(0), Q(0)]
    for _ in range(n):
        out = mul(out, a)
    return out


def scale(a, c):
    return [Q(c) * value for value in a]


def monomial_derivative(i, j, de, dx, e, x):
    if i < de or j < dx:
        return Q(0)
    return (Q(math.factorial(i), math.factorial(i-de))
            * Q(math.factorial(j), math.factorial(j-dx))
            * e**(i-de) * x**(j-dx))


def exp_coefficient(k, rate):
    return rate**k / math.factorial(k)


def phi_coefficient(k, rate):
    return Q(0) if k == 0 else rate**(k-1) / math.factorial(k)


def run_checks():
    checks = []

    def check(condition, name):
        require(condition, name)
        checks.append(name)

    # Fundamental theorem on both complete convex segments. Exact integral
    # of d/dtheta F(e,x0+theta*(x-x0)), followed by parameter segment.
    samples = [(Q(-2,3), Q(5,7), Q(-1,4)),
               (Q(2,5), Q(-3,2), Q(1,3)),
               (Q(0), Q(1), Q(1))]
    for i in range(6):
        for j in range(6-i):
            for k, (e,x,x0) in enumerate(samples):
                delta = x-x0
                state_integral = Q(0)
                if j:
                    state_integral = e**i * j * delta * sum(
                        (Q(math.comb(j-1,m), m+1) * x0**(j-1-m)
                         * delta**m for m in range(j)), Q(0))
                parameter_integral = e**i*x0**j if i else Q(0)
                total = e**i*x**j - (Q(0) if i else x0**j)
                check(total == state_integral + parameter_integral,
                      f"segment_FTC_{i}_{j}_{k}")

    # Exact formal Taylor coefficients of R'=LR+P1 and phi'=L phi+1.
    for l in [Q(0), Q(1,3), Q(7,2)]:
        for k in range(10):
            eco = exp_coefficient(k,l)
            pco = phi_coefficient(k,l)
            check(eco >= 0 and pco >= 0,
                  f"nonnegative_series_{l}_{k}")
            check((k+1)*exp_coefficient(k+1,l) == l*eco,
                  f"exp_ode_{l}_{k}")
            check((k+1)*phi_coefficient(k+1,l)
                  == l*pco + (1 if k == 0 else 0),
                  f"phi_ode_{l}_{k}")
            v,p = Q(5,11),Q(7,13)
            rco = eco*v+pco*p
            next_rco = exp_coefficient(k+1,l)*v+phi_coefficient(k+1,l)*p
            check((k+1)*next_rco == l*rco+(p if k == 0 else 0),
                  f"inhomogeneous_majorant_{l}_{k}")
            if l == 0:
                check(pco == (Q(1) if k == 1 else Q(0)),
                      f"zero_L_phi_{k}")

    # Chain-rule identities through second order; mixed term must be 2.
    e,x,v,w = Q(-2,3),Q(3,5),Q(7,11),Q(-5,13)
    for i in range(6):
        for j in range(6-i):
            coefficients = mul(power([e,Q(1),Q(0)],i),
                               power([x,v,w/2],j))
            d = lambda de,dx: monomial_derivative(i,j,de,dx,e,x)
            first = d(1,0)+d(0,1)*v
            second = d(2,0)+2*d(1,1)*v+d(0,2)*v*v+d(0,1)*w
            check(coefficients[1] == first, f"first_parameter_chain_{i}_{j}")
            check(2*coefficients[2] == second,
                  f"second_parameter_chain_{i}_{j}")
    # Transition Hessian term is the parameter-free case.
    for j in range(2,7):
        coefficients = power([x,v,w/2],j)
        check(2*coefficients[2] == j*x**(j-1)*w
              +j*(j-1)*x**(j-2)*v*v, f"transition_Hessian_{j}")

    # Clock rescaling is exact coefficient by coefficient, including L=0.
    for rho in [Q(1,7),Q(1),Q(17,3)]:
        for l in [Q(0),Q(2,5)]:
            p = Q(3,11)
            for k in range(10):
                check(exp_coefficient(k,l)*rho**k
                      ==exp_coefficient(k,rho*l),
                      f"rescaled_exp_{rho}_{l}_{k}")
                check(phi_coefficient(k,l)*rho**k*p
                      ==phi_coefficient(k,rho*l)*rho*p,
                      f"rescaled_forcing_{rho}_{l}_{k}")

    # Units are exponents of (scaled state, parameter, clock).
    state=(1,0,0); par=(0,1,0); clock=(0,0,1)
    unit_add=lambda *items: tuple(sum(item[k] for item in items) for k in range(3))
    l=(0,0,-1); p1=(1,-1,-1); phi=(0,0,1)
    v0=(1,-1,0); w0=(1,-2,0)
    h=(-1,0,-1); c=(0,-1,-1); p2=(1,-2,-1)
    check(unit_add(l,clock)==(0,0,0), "dimensionless_LT")
    check(unit_add(phi,p1)==v0, "first_forcing_units")
    check(unit_add(par,v0)==state, "radius_displacement_units")
    for label,forcing in [("H_V_squared",unit_add(h,v0,v0)),
                          ("C_V",unit_add(c,v0)),("P2",p2)]:
        check(unit_add(phi,forcing)==w0, "second_forcing_units_"+label)

    # L=0 affine ODE has exact radius majorant; equality can touch boundary.
    a,v0,p,t = Q(1,5),Q(2,3),Q(4,7),Q(3,2)
    bound=a*(v0+t*p)
    check(bound < bound+Q(1,17), "strict_guard_positive_margin_example")
    check(a*v0 < bound, "initial_interior_from_positive_duration_example")
    check(bound == a*v0+a*p*t, "affine_zero_L_exact_solution")
    check(not bound < bound, "equality_guard_not_strict")
    check(Q(0)*v0+Q(0)*p == 0, "zero_time_zero_displacement")
    check(exp_coefficient(0,Q(3))*v0+phi_coefficient(0,Q(3))*p==v0,
          "zero_time_majorant")
    r=bound+Q(1,17); intermediate=(r+bound)/2
    check(bound < intermediate < r, "compact_subradius_margin")

    # Exact retained horizon values. These are geometric bookkeeping, not
    # an evaluated coefficient or a physical nonlinear amplitude radius.
    ell=Q(1,10**6); tq=Q(25,16)*ell; tc=Q(13,8)*ell
    coarse_max=Q(6145,4096)*ell; coarse_min=Q(6143,4096)*ell
    guard=ell/32
    check(tc-tq==ell/16, "background_horizon_margin_ell_over_16")
    check(tc-tq-guard==ell/32, "perturbed_lower_margin_half_guard")
    check(tq-coarse_max==Q(255,4096)*ell, "coarse_root_horizon_separation")
    check(tc-coarse_max==Q(511,4096)*ell, "old_coarse_slab_separation")
    check(Q(0)<coarse_min<coarse_max<tq<tc, "complete_connecting_interval")
    check(tq+guard<tc, "time_projection_guard_fits_slab")
    check(tq/Q(50)==ell/32, "velocity_guard_one_over_50_gives_position_guard")
    check(1-Q(1,50)==Q(49,50), "velocity_guard_strict_clock_margin")
    check(tq==Q(1,640000), "Tq_exact_coordinate_number")
    check(tc==Q(13,8000000), "Tcoord_exact_coordinate_number")

    digests = {}
    for name in ["FLOW_BOOTSTRAP.md", "SOURCE_AUDIT.md", "check_flow_bootstrap.py"]:
        digests[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return {
        "schema":"hdblast.flow_bootstrap.exact_checks.v1",
        "status":"pass", "check_count":len(checks), "checks":checks,
        "claim_scope":"finite exact algebra only; analytic ODE argument is in FLOW_BOOTSTRAP.md",
        "physical_inputs_evaluated":False,
        "effective_amplitude_claimed":False,
        "source_sha256":digests,
        "exact_clock_values":{
            "ell":str(ell),"T_q":str(tq),"T_coord":str(tc),
            "background_margin":str(tc-tq),
            "guard_consumed_at_most":str(guard),
            "strict_remaining_margin_greater_than":str(tc-tq-guard),
            "coarse_root_to_horizon":str(tq-coarse_max)},
    }


def output_path(text):
    path=Path(text)
    if not path.is_absolute():
        raise argparse.ArgumentTypeError("--output must be absolute")
    if not path.parent.is_dir():
        raise argparse.ArgumentTypeError("output parent must already be a directory")
    if path.parent.resolve()!=path.parent:
        raise argparse.ArgumentTypeError("output parent must be canonical and contain no symlink")
    if os.path.lexists(path):
        raise argparse.ArgumentTypeError("output must be fresh; existing paths and symlinks are rejected")
    return path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=output_path)
    args=parser.parse_args()
    payload=(json.dumps(run_checks(),sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode("ascii")
    flags=os.O_WRONLY|os.O_CREAT|os.O_EXCL|getattr(os,"O_NOFOLLOW",0)
    descriptor=os.open(args.output,flags,0o600)
    with os.fdopen(descriptor,"wb") as handle:
        handle.write(payload)
    sys.stdout.buffer.write(payload)


if __name__=="__main__":
    main()
