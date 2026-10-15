#!/usr/bin/env python3
"""Exact algebra/arithmetic audit of UNIFORM_POINTWISE_WAVEFORM.md.

This is not a formal PDE or existence proof. No checks depend on `assert`.
Run with normal Python and `python -O`; the JSON receipts must agree.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

import sympy as s


class AuditFailure(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise AuditFailure(message)


def zero(expression):
    require(s.simplify(expression) == 0, "symbolic identity rejected")


def rational(value):
    require(isinstance(value, (int, F)), "only integers and Fractions accepted")
    return F(value)


def majorant(ell, T, m, M, b, g):
    ell, T, m, M, b, g = map(rational, (ell, T, m, M, b, g))
    require(ell > 0 and T > 0 and m > 0, "nonpositive length/time/inertia")
    require(M >= 0 and b >= 0 and g >= 0, "negative coefficient bound")
    U = 1 + T
    A = (1 + ell) / m
    k = M * T**2 / 2
    c = k / m
    d = b * T**2 / (2 * m)
    den = (1-k)*(1-d)-c*T
    require(k < 1, "bulk denominator is not positive")
    require(den > 0, "coupled denominator is not positive")
    V = (A*(1-k)+c*U) / den
    E = (T*V+k*U) / (1-k)
    eps = ell*(M*T*(U+E)+(1+g*T)*V)
    require(E == T*V+k*(U+E), "E does not solve majorant equality")
    require(V == A+c*(U+E)+d*V, "V does not solve majorant equality")
    return dict(ell=ell, T=T, m=m, M=M, b=b, g=g, U=U, A=A, k=k,
                c=c, d=d, denominator=den, V_star=V, E_star=E,
                epsilon_star=eps)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    checks = []
    controls = []

    def check(name, fn):
        fn()
        checks.append(name)

    def reject(name, fn):
        try:
            fn()
        except AuditFailure:
            controls.append(name)
        else:
            raise AuditFailure(f"wrong-formula control accepted: {name}")

    t, z = s.symbols("t z", real=True)
    e = s.Function("e")(t, z)
    et, ez = s.diff(e,t), s.diff(e,z)
    check("incoming_characteristic_sign", lambda: zero(
        s.diff(et+ez,t)-s.diff(et+ez,z)-s.diff(e,t,2)+s.diff(e,z,2)))
    check("outgoing_characteristic_sign", lambda: zero(
        s.diff(et-ez,t)+s.diff(et-ez,z)-s.diff(e,t,2)+s.diff(e,z,2)))
    reject("wrong_incoming_transport_direction", lambda: zero(
        s.diff(et+ez,t)+s.diff(et+ez,z)-s.diff(e,t,2)+s.diff(e,z,2)))

    mu, v, a, beta, qb, J, G, gamma, dhat, r = s.symbols(
        "mu v a beta qb J G gamma dhat r", nonzero=True)
    true_boundary = -(J-v)+mu*a+beta*qb-G
    check("boundary_ode_has_positive_damping_and_forcing", lambda: zero(
        true_boundary-(mu*a+v+beta*qb-J-G)))
    reject("negative_damping", lambda: zero(
        true_boundary-(mu*a-v+beta*qb-J-G)))
    reject("wrong_boundary_forcing_sign", lambda: zero(
        true_boundary-(mu*a+v+beta*qb-J+G)))
    check("dirichlet_comparator_residual_and_error_source", lambda: zero(
        (-r).subs(r,-dhat)-dhat))
    reject("discarded_comparator_boundary_residual", lambda: zero(dhat))
    check("observed_trace_includes_velocity_and_gamma", lambda: zero(
        ((J-v)-gamma*qb)-(J-v-gamma*qb)))
    reject("discarded_boundary_velocity", lambda: zero(
        ((J-v)-gamma*qb)-(J-gamma*qb)))
    reject("discarded_gamma_trace", lambda: zero(
        ((J-v)-gamma*qb)-(J-v)))
    reject("full_C_b_instead_of_half", lambda: zero(
        ((J-v)-gamma*qb)-(J-v-2*gamma*qb)))

    # I=integral_0^t exp(-(t-s)/mu)q(s)ds solves I'=q-I/mu.
    q=s.Function("q")(t)
    I=s.Function("I")(t)
    integ_derivative = q-I/mu
    damped_primitive = q-I/mu
    primitive_ode = s.diff(damped_primitive,t)+damped_primitive/mu-s.diff(q,t)
    check("integration_by_parts_negative_convolution", lambda: zero(
        primitive_ode.subs(s.diff(I,t),integ_derivative)))
    bad_primitive=q+I/mu
    reject("positive_convolution_in_integration_by_parts", lambda: zero(
        (s.diff(bad_primitive,t)+bad_primitive/mu-s.diff(q,t)).subs(
            s.diff(I,t),integ_derivative)))
    reject("drop_convolution_in_integration_by_parts", lambda: zero(q/mu))
    s_var=s.symbols("s", real=True)
    kernel=s.exp(-(t-s_var)/mu)
    check("kernel_derivative_with_respect_to_source_time", lambda: zero(
        s.diff(kernel,s_var)-kernel/mu))

    ell,x,Q,Qp=s.symbols("ell x Q Qp", positive=True)
    scaled=ell*(Qp/ell+s.exp(-ell*x)*Q)
    correct=Qp+ell*s.exp(-ell*x)*Q
    check("dimensionless_waveform_scaling", lambda: zero(scaled-correct))
    reject("missing_ell_on_velocity_data", lambda: zero(
        scaled-(Qp+s.exp(-ell*x)*Q)))
    reject("missing_ell_on_observation", lambda: zero(
        scaled/ell-correct))

    result=majorant(F(1,10**10),F(3,10**10),F(1,1000),10**6,1000,250)
    eps_expected=F(8000003000800720779892143999973,
                   79999999999992799998920000162000000000)
    target=F(1000001,10**13)
    check("exact_rational_epsilon", lambda: require(
        result["epsilon_star"]==eps_expected,"exact epsilon mismatch"))
    check("advertised_strict_uniform_bound", lambda: require(
        result["epsilon_star"]<target,"target bound not proved"))
    clock_correction=F(3,4)*result["T"]**2*(16+result["ell"])
    check("signed_physical_normalization_same_rounded_bound", lambda: require(
        result["epsilon_star"]+clock_correction<target,
        "clock correction exceeds target"))
    ell0,m0=result["ell"],result["m"]
    Vtail=(ell0/m0**2+ell0/m0+result["c"]*(result["U"]+result["E_star"])
           +result["d"]*result["V_star"])
    etail=ell0*(result["M"]*result["T"]*(result["U"]+result["E_star"])
                +Vtail+result["g"]*result["T"]*result["V_star"])
    check("post_pulse_tail_exact_fraction", lambda: require(
        etail==F(190040072038996407198648650081,
                 3999999999999639999946000008100000000000000),
        "tail fraction mismatch"))
    check("post_pulse_tail_strict_bound", lambda: require(
        etail<F(4752,10**17),"tail target not proved"))
    check("support_aware_L2_bound", lambda: require(
        result["epsilon_star"]**2+etail**2<F(10000005,10**14)**2,
        "support-aware L2 bound not proved"))
    reject("float_inputs", lambda: majorant(1e-10,F(3,10**10),F(1,1000),
                                           10**6,1000,250))
    reject("negative_inertia", lambda: majorant(F(1,10**10),F(3,10**10),
                                               F(-1,1000),10**6,1000,250))
    reject("nonpositive_majorant_denominator", lambda: majorant(
        F(1,10**10),1,F(1,1000),0,1000,250))
    reject("claim_too_small_epsilon", lambda: require(
        result["epsilon_star"]<F(1,10**7),"too-small target rejected"))

    receipt={
        "status":"PASS",
        "scope":"exact algebra and rational arithmetic conditional on analytic proof",
        "model":"actual-background positive-inertia first variation; frozen source retained",
        "checks":checks,
        "wrong_formula_controls_rejected":controls,
        "values":{k:str(v) for k,v in result.items()},
        "strict_rounded_epsilon":str(target),
        "strict_margin":str(target-result["epsilon_star"]),
        "physical_clock_correction_bound":str(clock_correction),
        "V_tail":str(Vtail),
        "epsilon_tail":str(etail),
        "strict_tail_target":str(F(4752,10**17)),
        "strict_support_aware_L2_target":str(F(10000005,10**14)),
        "p_b_sign":"negative on retained actual branch; raw S sign is opposite H",
        "unproved_by_this_program":[
            "analytic PDE characteristic argument and smooth existence",
            "actual background and exact-family hypotheses",
            "nonlinear remainder or finite amplitude radius",
            "physical scale calibration",
        ],
    }
    output=json.dumps(receipt,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output,end="")


if __name__=="__main__":
    main()
