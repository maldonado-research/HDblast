#!/usr/bin/env python3
"""Portable exact stress-response proof; no numerical source or mode evaluation.

Run with --repo-root CHECKOUT --output NEW_JSON. All checks use explicit
exceptions and remain active under Python optimization. Inherited sources are
read and hashed, never imported or executed. Output files cannot be overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import sympy as S


PINS = {
    "research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM/theory/PROSPECTIVE_CAUSAL_RESPONSE_PROTOCOL.md": "2bcc72013110d14f4c076f8270ee6b83057a960405f47a27b7d8c4445f82083e",
    "research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/COMMON_ACTION_SMOOTH_FRW.md": "4896a6f83c4d63ac8300608bba923c7a7e45b0459d70aad621b400c001ee302f",
    "research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/code/derive_local_trace.py": "25a9703f720fe068e97b1fd6731808751d66d90960fd40b0607e79a48276e0f0",
    "research/HDBLAST_CHECKPOINT_20261002_SMOOTH_FRW/code/verify_local_identities.py": "d392740bc1926474c33fe37503c9c1c3406b63b0f7907218857c5d6af940daf0",
    "research/HDBLAST_CHECKPOINT_20261002_DESITTER/theory/DESITTER_COMMON_ACTION_DERIVATION.md": "fe84d5c108dc688b8d261f078786ac7a1483510381052c95ecc827013b89214b",
    "research/HDBLAST_CHECKPOINT_20261002_DESITTER/independent/exact_mode_bridge.py": "501a2dc5edbf3aed7bbf64d407c2acc8174b3479a555e344bb61508b85f333aa",
    "research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM/theory/STATIONARY_COMMON_ACTION_MODEL.md": "84364b8271fc73ced365202e820908bca2d23ab0760a6ea067e7c775cf2dddf5",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prove(report):
    exact, mutations = report["exact_identities"], report["detected_mutations"]

    def equal(name, lhs, rhs=0):
        defect = S.simplify(S.expand(lhs-rhs))
        require(defect == 0, f"{name}: nonzero exact defect {defect}")
        exact.append(name)

    def wrong(name, defect):
        residual = S.factor(S.simplify(defect))
        require(residual != 0, f"{name}: wrong formula was not detected")
        mutations.append({"name": name, "exact_residual": str(residual)})

    # Bare physical stress at the conformal-frequency reference.
    a, L, k = S.symbols("a L k", positive=True)
    A, Ap, source, sourcep = S.symbols("A Ap s sp", real=True)
    App = -4*k*k*A-source/k
    kinetic = (L*L-k*k)*A-L*Ap
    eb = (kinetic+(k*k+2*L*L)*A+source/(2*k))/(2*a**4)
    pb = (kinetic-(k*k/3+2*L*L)*A-source/(2*k))/(2*a**4)
    qb = A/a**2
    qbp = (Ap-2*L*A)/a**2
    qbpp = (App-4*L*Ap+2*L*L*A)/a**2
    q0b = 1/(2*k*a**2)
    d = source/a**2
    dp = (sourcep-2*L*source)/a**2
    equal("bare_minimal_density", eb, (3*L*L*A-L*Ap+source/(2*k))/(2*a**4))
    equal("bare_minimal_pressure", pb, (-(4*k*k/3+L*L)*A-L*Ap-source/(2*k))/(2*a**4))
    equal("bare_density_variance_reduction", eb, (L*L*qb-L*qbp)/(2*a**2)+q0b*d/2)
    equal("bare_pressure_variance_reduction", pb, (qbpp+L*qbp-3*L*L*qb)/(6*a**2)-q0b*d/6)
    equal("bare_minimal_trace", -eb+3*pb, -2*L*L*qb/a**2-q0b*d+(qbpp+2*L*qbp)/(2*a**2))
    Deta = lambda expr: (S.diff(expr,a)*L*a+S.diff(expr,L)*L*L
                        +S.diff(expr,A)*Ap+S.diff(expr,Ap)*App
                        +S.diff(expr,source)*sourcep)
    equal("bare_source_Ward_identity", Deta(eb)+3*L*(eb+pb), q0b*dp/2)
    wrong("omit_explicit_stress_mass_contacts", Deta(eb-q0b*d/2)+3*L*(eb+pb)-q0b*dp/2)

    # Direct W2/W4 stress subtraction at unit observation scale factor.
    # H and w are positive; source jets are arbitrary real symbols.
    w, r, H = S.symbols("w r H", positive=True)
    dx, dx1, dx2 = S.symbols("dx dx_dot dx_ddot", real=True)
    sp = dx1+2*H*dx
    spp = dx2+5*H*dx1+6*H*H*dx
    wp = r*H/w
    wpp = 3*r*H*H/w-r*r*H*H/w**3
    U = -H*H/w-wpp/(4*w*w)+3*wp*wp/(8*w**3)
    u = dx/(2*w)
    up = sp/(2*w)-dx*wp/(2*w*w)
    upp = spp/(2*w)-sp*wp/(w*w)-dx*wpp/(2*w*w)+dx*wp*wp/w**3
    v = (-U*u/w-upp/(4*w*w)+wpp*u/(4*w**3)
         +3*wp*up/(4*w**3)-3*wp*wp*u/(4*w**4))
    j4 = H*up/w**2-2*H*wp*u/w**3+wp*up/(2*w**3)-3*wp*wp*u/(4*w**4)
    b = -(w*w+2*r)/3
    c = 1-b/(w*w)
    esub_raw = (dx/w+2*U*u/w-dx*U/(w*w)-H*H*u/(w*w)+j4)/4
    esub = (dx/w-H*H*u/(w*w)+j4)/4
    psub = (c*u-dx/w+c*v+2*b*U*u/w**3+dx*U/w**2-H*H*u/w**2+j4)/4
    equal("full_fourth_order_density_variation_cancellation", esub_raw, esub)
    q = -dx/(4*w**3)
    qp = -sp/(4*w**3)+H*dx/(2*w**3)+3*dx*wp/(4*w**4)
    qpp = (-spp/(4*w**3)+H*sp/w**3-H*H*dx/(2*w**3)
           +3*sp*wp/(2*w**4)-3*H*dx*wp/w**4
           +3*dx*wpp/(4*w**4)-3*dx*wp*wp/w**5)
    q0 = 1/(2*w)-U/(2*w*w)
    eremainder = S.expand((H*H*q-H*qp)/2+q0*dx/2-esub)
    premainder_general = S.expand((qpp+H*qp-3*H*H*q)/6-q0*dx/6-psub)
    equal("direct_density_local_integrand", eremainder, r*H*(H*dx-dx1)/(16*w**5))
    premainder = S.expand(premainder_general.subs(r,2*H*H))
    pexpected = H*H*(70*H**6*dx-50*H**4*dx*w*w-10*H**3*dx1*w*w
                    -5*H*H*dx*w**4+4*H*dx1*w**4+dx2*w**4)/(24*w**9)
    equal("direct_pressure_local_integrand", premainder, pexpected)
    wrong("omit_reference_specialization_r_equals_2H_squared",
          premainder_general.coeff(w,-3))

    def integrate_rational_mode(expr, mass2):
        # Exact beta-function moment, excluding the common 1/(2 pi^2).
        result = 0
        for term in S.Add.make_args(S.expand(expr)):
            power = term.as_powers_dict().get(w,S.Integer(0))
            n = -power
            require(n.is_integer is True and n>3,
                    f"Nonintegrable rational subtraction mismatch: {term}")
            moment = (S.sqrt(S.pi)*S.gamma((n-3)/2)*mass2**((3-n)/2)
                      /(4*S.gamma(n/2)))
            result += term/w**power*moment
        return S.simplify(result)

    econ = integrate_rational_mode(eremainder,r)/(2*S.pi**2)
    pcon = integrate_rational_mode(premainder,2*H*H)/(2*S.pi**2)
    equal("integrated_density_contact", econ, (H*H*dx-H*dx1)/(96*S.pi**2))
    equal("integrated_pressure_contact", pcon, (dx2+2*H*dx1-11*H*H*dx)/(288*S.pi**2))

    # Independent cosmic-time jets check the final continuum identities.
    Q, Q1, Q2 = S.symbols("deltaQ deltaQ_dot deltaQ_ddot", real=True)
    Q0 = H*H/(12*S.pi**2)
    rho = H*H*Q/2-H*Q1/2+Q0*dx/2+econ
    pressure = Q2/6+H*Q1/3-H*H*Q/2-Q0*dx/6+pcon
    Dt = lambda expr: (S.diff(expr,Q)*Q1+S.diff(expr,Q1)*Q2
                      +S.diff(expr,dx)*dx1+S.diff(expr,dx1)*dx2)
    anomaly = -H*H*dx/(8*S.pi**2)+(dx2+3*H*dx1)/(96*S.pi**2)
    trace = -2*H*H*Q-Q0*dx+(Q2+3*H*Q1)/2+anomaly
    equal("matched_continuum_Ward_identity", Dt(rho)+3*H*(rho+pressure), Q0*dx1/2)
    equal("matched_continuum_trace_identity", -rho+3*pressure, trace)
    wrong("drop_density_contact", Dt(rho-econ)+3*H*(rho-econ+pressure)-Q0*dx1/2)
    wrong("drop_pressure_contact", Dt(rho)+3*H*(rho+pressure-pcon)-Q0*dx1/2)
    wrong("reverse_boxQ_sign", -rho+3*pressure-(-2*H*H*Q-Q0*dx-(Q2+3*H*Q1)/2+anomaly))
    wrong("omit_explicit_Q0_density_term", Dt(rho-Q0*dx/2)+3*H*(rho-Q0*dx/2+pressure)-Q0*dx1/2)

    # Finite-cutoff contacts use independently tabulated inherited moments.
    t = S.symbols("cutoff_ratio", nonnegative=True)
    moments = {}
    for n in (5,7,9):
        degree = (n-5)//2
        primitive = sum((-1)**j*S.binomial(degree,j)*t**(2*j+3)/S.Integer(2*j+3)
                        for j in range(degree+1))
        equal(f"finite_cutoff_moment_primitive_n{n}", S.diff(primitive,t), t*t*(1-t*t)**degree)
        moments[n] = r**S.Rational(3-n,2)*primitive
    anomalyK = (r*(-6*H*H*dx+5*H*dx1+dx2)*moments[5]/16
                -r*r*(50*H*H*dx+10*H*dx1)*moments[7]/32
                +70*r**3*H*H*dx*moments[9]/64)/(2*S.pi**2)
    equal("finite_cutoff_anomaly_continuum_limit", anomalyK.subs(t,1), anomaly)
    equal("finite_cutoff_density_contact", r*H*(H*dx-dx1)*moments[5]/(32*S.pi**2), econ*t**3)
    wrong("use_continuum_anomaly_at_finite_cutoff", anomalyK-anomaly)

    # Separate past-infinite stationary control, fixed r throughout.
    Qx = -(2*S.EulerGamma+S.log(2))/(16*S.pi**2)
    Qz = (-S.Rational(2,3)+4*S.EulerGamma+2*S.log(2))/(16*S.pi**2)
    rhox_closed = S.diff(rho.subs({Q:Qx*dx,Q1:0,dx1:0}),dx)
    rhox_action = Q0/2-H*H*Qz/4
    rhox_target = H*H*(S.Rational(5,3)-2*S.EulerGamma-S.log(2))/(32*S.pi**2)
    equal("stationary_rhox_common_action", rhox_closed, rhox_action)
    equal("stationary_rhox_closed_value", rhox_closed, rhox_target)
    px_closed = S.diff(pressure.subs({Q:Qx*dx,Q1:0,Q2:0,dx1:0,dx2:0}),dx)
    equal("stationary_pressure_deSitter_invariance", px_closed, -rhox_closed)
    wrong("moving_reference_static_derivative", H*H/(32*S.pi**2))
    report["derived_contacts"] = {"delta_rho":str(econ), "delta_p":str(pcon)}
    report["stationary_rho_x"] = str(rhox_closed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), f"Refusing to overwrite {args.output}")
    root = args.repo_root.resolve()
    report = {"status":"RUNNING", "scope":"Exact algebra only; no numerical source/response/mode evaluation",
              "python":platform.python_version(), "sympy":S.__version__,
              "executable":sys.executable, "optimization_level":sys.flags.optimize,
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "repo_root":str(root), "inherited_source_pins":[],
              "exact_identities":[], "detected_mutations":[]}
    error = None
    try:
        for relative, expected in PINS.items():
            path = root/relative
            actual = hashlib.sha256(path.read_bytes()).hexdigest()
            report["inherited_source_pins"].append({"path":relative,"resolved_path":str(path),
                                                     "expected_sha256":expected,"actual_sha256":actual})
            require(actual==expected, f"Inherited source hash mismatch: {relative}")
        prove(report)
        report["status"] = "PASS"
    except Exception as exc:
        error = exc
        report["status"] = "FAIL"
        report["failure"] = traceback.format_exc()
    report["exact_identity_count"] = len(report["exact_identities"])
    report["detected_mutation_count"] = len(report["detected_mutations"])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report,handle,indent=2)
        handle.write("\n")
    print(json.dumps(report,indent=2))
    if error is not None:
        raise RuntimeError("Exact proof failed; preserved output contains the original exception") from error


if __name__ == "__main__":
    main()
