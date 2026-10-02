#!/usr/bin/env python3
"""Pure symbolic implementation audit; never samples a physical source."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback

import sympy as S
from source_jet import POLYNOMIALS


def audit():
    checks = []
    def equal(name, left, right):
        if S.simplify(left-right) != 0:
            raise RuntimeError("Pure symbolic implementation mismatch: "+name)
        checks.append({"name": name, "passed": True})
    u = S.symbols("u", real=True)
    d = 1-u*u
    B = S.exp(-u*u/d)
    polynomials = [sum(S.Integer(c)*u**(len(coefficients)-index-1)
                       for index,c in enumerate(coefficients)) for coefficients in POLYNOMIALS]
    for n in range(6):
        expression = B*polynomials[n]/d**(2*n)
        equal("bump_jet_"+str(n), expression, S.diff(B,u,n))
        signed = u*expression+(n*B*polynomials[n-1]/d**(2*n-2) if n else 0)
        equal("signed_jet_"+str(n), signed, S.diff(u*B,u,n))
        if n < 5:
            equal("integer_polynomial_recurrence_"+str(n), polynomials[n+1],
                  S.expand(d*d*S.diff(polynomials[n],u)+4*n*u*d*polynomials[n]-2*u*polynomials[n]))
    K,M,L = S.symbols("K M L", positive=True)
    w = S.sqrt(K*K+M*M)
    v = K/w
    A = S.asinh(K/M)-v
    D = lambda expression: L*M*S.diff(expression,M)+L*L*S.diff(expression,L)
    equal("finite_local_A_first_derivative", D(A), -L*v**3)
    equal("finite_local_A_second_derivative", D(D(A)), L*L*v**3*(2-3*v*v))
    vv = S.symbols("v", positive=True)
    J5 = vv**3/6
    J7 = (vv**3/3-vv**5/5)/4
    J9 = (vv**3/3-2*vv**5/5+vv**7/7)/8
    f0,f1,f2,f3 = S.symbols("f0 f1 f2 f3", real=True)
    z0,z1,z2 = S.symbols("z0 z1 z2", real=True)
    Q0,Q0p = S.symbols("Q0 Q0p", real=True)
    rho_contact = (3*L*L*f0-L*f1)*vv**3/(96*S.pi**2)
    p_contact = (70*L*L*f0*J9-(30*L*L*f0+10*L*f1)*J7+(f2-L*f1-9*L*L*f0)*J5)/(48*S.pi**2)
    anomaly = ((f2-12*L*L*f0)*J5-(30*L*L*f0+10*L*f1)*J7+70*L*L*f0*J9)/(16*S.pi**2)
    equal("direct_pressure_contact_continuum",p_contact.subs(vv,1),(f2-3*L*f1-11*L*L*f0)/(288*S.pi**2))
    equal("direct_density_contact_continuum",rho_contact.subs(vv,1),(3*L*L*f0-L*f1)/(96*S.pi**2))
    equal("finite_trace_contact_cross_check",-rho_contact+3*p_contact,anomaly)
    equal("continuum_anomaly_contact",anomaly.subs(vv,1),(f2-2*L*f1-14*L*L*f0)/(96*S.pi**2))
    baseline_polynomial = vv**2/(2*(1+vv))-vv**3/48-vv**5/16
    equal("baseline_polynomial_derivative",S.diff(baseline_polynomial,vv),vv*(2+vv)/(2*(1+vv)**2)-vv**2/16-5*vv**4/16)
    equal("reference_Q0_continuum",baseline_polynomial.subs(vv,1)/(2*S.pi**2),1/(12*S.pi**2))
    w_ref = S.sqrt(K*K+2*L*L)
    v_ref = K/w_ref
    U = -L*L/w_ref-3*L**4/(2*w_ref**3)+5*L**6/(2*w_ref**5)
    raw_reference_integrand = K*K/(2*S.pi**2*L*L)*(1/(2*K)-1/(2*w_ref)+U/(2*w_ref**2))
    closed_reference = baseline_polynomial.subs(vv,v_ref)/(2*S.pi**2)
    equal("closed_Q0K_matches_complete_subtraction_integrand",S.diff(closed_reference,K),raw_reference_integrand)
    equal("closed_Q0K_lower_endpoint",S.limit(closed_reference,K,0),0)
    rho = (3*L*L*z0-L*z1)/2+L*L*Q0*f0/2+rho_contact
    p = (z2-3*L*z1-3*L*L*z0)/6-L*L*Q0*f0/6+p_contact
    density_differential = (L*L*S.diff(rho,L)+f1*S.diff(rho,f0)+f2*S.diff(rho,f1)
                            +z1*S.diff(rho,z0)+z2*S.diff(rho,z1)+Q0p*S.diff(rho,Q0)
                            -L*vv*(1-vv*vv)*S.diff(rho,vv)-4*L*rho)
    Dcontact = ((6*L**3*f0+2*L*L*f1-L*f2)*vv**3
                -(3*L*L*f0-L*f1)*3*L*vv**3*(1-vv*vv))/(96*S.pi**2)
    implemented_density_prime = (-3*L**3*z0+3*L*L*z1-L*z2/2
                                  +L*L*((Q0p/2-L*Q0)*f0+Q0*f1/2)+Dcontact-4*L*rho_contact)
    equal("direct_density_derivative",density_differential,implemented_density_prime)
    equal("finite_trace_variance_and_mass_terms",-rho+3*p,
          (z2-2*L*z1-6*L*L*z0)/2-L*L*Q0*f0+anomaly)
    q0ofv = baseline_polynomial/(2*S.pi**2)
    q0prime = S.diff(q0ofv,vv)*(-L*vv*(1-vv*vv))
    ward = (implemented_density_prime+3*L*(rho+p)-L*L*Q0*(f1-2*L*f0)/2).subs({Q0:q0ofv,Q0p:q0prime})
    equal("independently_defined_density_pressure_Ward_algebra",ward,0)
    return {"passed": True, "check_count": len(checks), "checks": checks,
            "scope": "Exact symbolic identities only. No source values, physical responses, quadratures, or modes evaluated."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("Refusing to overwrite symbolic development evidence")
    provenance = {"verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  "source_jet_sha256": hashlib.sha256(Path(__file__).with_name("source_jet.py").read_bytes()).hexdigest(),
                  "sympy": S.__version__, "python_optimization": sys.flags.optimize}
    try:
        result = audit()
        result["provenance"] = provenance
    except BaseException as exc:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({"passed": False, "exception": repr(exc), "traceback": traceback.format_exc(),
                                           "provenance": provenance}, indent=2)+"\n")
        raise
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print("Pure symbolic primary implementation audit passed:",result["check_count"],"identities")


if __name__ == "__main__":
    main()
