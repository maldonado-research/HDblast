#!/usr/bin/env python3
"""Exact algebra supporting stress UV bounds; no numerical source evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback
import sympy as S


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prove(report):
    identities, mutations = report["identities"], report["mutations"]
    def eq(name, actual, expected=0):
        residual=S.factor(S.simplify(actual-expected))
        require(residual==0,f"{name}: {residual}")
        identities.append(name)
    def wrong(name, residual):
        residual=S.factor(S.simplify(residual))
        require(residual!=0,f"Mutation undetected: {name}")
        mutations.append({"name":name,"residual":str(residual)})
    v,L,M,K,k,H,r,a=S.symbols("v L M K k H r a",positive=True)
    f,f1,f2,q,q1,q2=S.symbols("f f1 f2 q q1 q2",real=True)
    C=1/(8*S.pi**2)
    dv=-L*v*(1-v*v)
    A1=-L*v**3
    A2=S.diff(A1,L)*L*L+S.diff(A1,v)*dv
    eq("cutoff_second_scale_derivative",A2,L*L*v**3*(2-3*v*v))
    eq("cutoff_first_scale_derivative",S.diff(S.asinh(K/M)-K/S.sqrt(K*K+M*M),M)*L*M,-L*K**3/(K*K+M*M)**S.Rational(3,2))
    eq("first_derivative_contact_difference",-L-A1,-L*(1-v**3))
    eq("second_derivative_contact_difference",-L*L-A2,-L*L*(1+2*v**3-3*v**5))
    eq("positive_second_contact_factor",1+2*v**3-3*v**5,(1-v)*(3*v**4+3*v**3+v*v+v+1))
    eq("positive_first_contact_factor",1-v**3,(1-v)*(1+v+v*v))
    eq("momentum_tail_factor",S.integrate(k**-3,(k,K,S.oo)),1/(2*K*K))
    w=S.symbols("w",positive=True)
    eq("subtraction_inequality_certificate",w**-3-1+S.Rational(3,2)*(w*w-1),(w-1)**2*(3*w**3+6*w*w+4*w+2)/(2*w**3))
    # w=sqrt(1+M^2/k^2)>=1; the factored expression is nonnegative.
    eq("stable_one_minus_cutoff_ratio",1-K/S.sqrt(K*K+M*M),M*M/(S.sqrt(K*K+M*M)*(S.sqrt(K*K+M*M)+K)))
    Q0=H*H/(12*S.pi**2)
    Q0K=H*H/(4*S.pi**2)*(v*v/(1+v)-v**3/24-v**5/8)
    DQ0=H*H*(1-v)**2*(v+2)*(3*v**3+3*v*v+10*v+4)/(96*S.pi**2*(v+1))
    eq("positive_reference_variance_tail",Q0-Q0K,DQ0)
    eq("reference_variance_continuum_limit",Q0K.subs(v,1),Q0)
    eq("reference_variance_time_derivative",S.diff(Q0K,v)*(-H*v*(1-v*v)),-H**3*v*v*(1-v)**2*(5*v**4+15*v**3+21*v*v+23*v+16)/(32*S.pi**2*(v+1)))
    moments={}
    tail_factors={
        5:(1-v)*(1+v+v*v)/3,
        7:(1-v)**2*(3*v**3+6*v*v+4*v+2)/15,
        9:(1-v)**3*(15*v**4+45*v**3+48*v*v+24*v+8)/105,
    }
    for n in (5,7,9):
        degree=(n-5)//2
        primitive=sum((-1)**j*S.binomial(degree,j)*v**(2*j+3)/S.Integer(2*j+3) for j in range(degree+1))
        eq(f"moment_{n}_primitive",S.diff(primitive,v),v*v*(1-v*v)**degree)
        eq(f"moment_{n}_positive_tail",primitive.subs(v,1)-primitive,tail_factors[n])
        moments[n]=primitive
    anomalyK=((f2-12*L*L*f)*moments[5]/16-(30*L*L*f+10*L*f1)*moments[7]/32+70*L*L*f*moments[9]/64)/(2*S.pi**2)
    anomalyInf=(f2-2*L*f1-14*L*L*f)/(96*S.pi**2)
    anomalyTail=((f2-12*L*L*f)*tail_factors[5]/16-(30*L*L*f+10*L*f1)*tail_factors[7]/32+70*L*L*f*tail_factors[9]/64)/(2*S.pi**2)
    eq("normalized_anomaly_continuum",anomalyK.subs(v,1),anomalyInf)
    eq("normalized_anomaly_tail",anomalyInf-anomalyK,anomalyTail)
    densityContact=L*(3*L*f-f1)/(96*S.pi**2)
    pressureContact=(f2-3*L*f1-11*L*L*f)/(288*S.pi**2)
    eq("normalized_pressure_contact",(anomalyInf+densityContact)/3,pressureContact)
    eq("combined_pressure_contact_tail",pressureContact-(anomalyK+densityContact*v**3)/3,(anomalyTail+densityContact*(1-v**3))/3)
    # Independent physical-variance conversion to normalized stress variables.
    physicalQ=q/a**2
    physicalQ1=(q1-2*L*q)/a**2
    physicalQ2=(q2-4*L*q1+2*L*L*q)/a**2
    eq("normalized_density_variance_part",a**2*(L*L*physicalQ-L*physicalQ1)/2,(3*L*L*q-L*q1)/2)
    eq("normalized_pressure_variance_part",a**2*(physicalQ2+L*physicalQ1-3*L*L*physicalQ)/6,(q2-3*L*q1-3*L*L*q)/6)
    # Exact finite-comoving-K Ward cancellation, with physical cosmic jets.
    dx,dx1,dx2,Q,Q1,Q2=S.symbols("dx dx_dot dx_ddot Q Q_dot Q_ddot",real=True)
    rhoLocal=(H*H*dx-H*dx1)/(96*S.pi**2)
    rhoK=H*H*Q/2-H*Q1/2+Q0K*dx/2+rhoLocal*v**3
    AK=((-6*H*H*dx+5*H*dx1+dx2)*moments[5]/16-(50*H*H*dx+10*H*dx1)*moments[7]/32+70*H*H*dx*moments[9]/64)/(2*S.pi**2)
    traceK=(Q2+3*H*Q1)/2-2*H*H*Q-Q0K*dx+AK
    pK=(traceK+rhoK)/3
    def Dt(expr):
        return (S.diff(expr,Q)*Q1+S.diff(expr,Q1)*Q2+S.diff(expr,dx)*dx1
                +S.diff(expr,dx1)*dx2-S.diff(expr,v)*H*v*(1-v*v))
    ward=Dt(rhoK)+3*H*(rhoK+pK)-Q0K*dx1/2
    eq("exact_finite_cutoff_Ward",ward)
    wrong("continuum_reference_in_finite_cutoff_Ward",Dt(rhoK)+3*H*(rhoK+pK)-Q0*dx1/2)
    Ainf=-H*H*dx/(8*S.pi**2)+(dx2+3*H*dx1)/(96*S.pi**2)
    wrong("continuum_anomaly_in_finite_cutoff_pressure",H*(Ainf-AK))
    wrong("omit_moving_scale_first_derivative_contact",C*L*f*(1-v**3))
    wrong("omit_moving_scale_second_derivative_contact",C*L*L*f*(1+2*v**3-3*v**5))
    wrong("omit_density_contact_cutoff_factor",densityContact*(1-v**3))
    wrong("omit_reference_tail_from_current",f*DQ0/4)
    report["proof_scope"]="Exact coefficient, factorization and finite-cutoff Ward identities; functional IBP bound and triangle inequalities are documented in STRESS_TAIL_DERIVATION.md; pulse root completeness is checked separately."


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),"Refusing to overwrite proof evidence")
    pins_path=Path(__file__).with_name("SOURCE_PINS.json")
    report={"status":"RUNNING","scope":"Pure analytic preparation; no source/response/mode numerical evaluation",
            "python":platform.python_version(),"sympy":S.__version__,"executable":sys.executable,
            "optimization":sys.flags.optimize,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "pins_sha256":hashlib.sha256(pins_path.read_bytes()).hexdigest(),"identities":[],"mutations":[],"input_pins":[]}
    failure=None
    try:
        for item in json.loads(pins_path.read_text())["files"]:
            path=args.repo_root/item["repository_relative_path"]
            actual=hashlib.sha256(path.read_bytes()).hexdigest()
            require(actual==item["sha256"],f"Inherited source mismatch: {path}")
            report["input_pins"].append(dict(item,actual_sha256=actual))
        prove(report)
        report["status"]="PASS"
    except Exception as exc:
        failure=exc
        report["status"]="FAIL"
        report["failure"]=traceback.format_exc()
    report["identity_count"]=len(report["identities"])
    report["mutation_count"]=len(report["mutations"])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open("x") as handle:
        json.dump(report,handle,indent=2)
        handle.write("\n")
    print(json.dumps(report,indent=2))
    if failure is not None:
        raise RuntimeError("Exact tail proof failed; output preserved") from failure


if __name__=="__main__":
    main()
