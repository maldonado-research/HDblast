"""Independent exact algebra only: no field evolution or numerical experiment."""
import json
from pathlib import Path
import hashlib

import sympy as s


def verify():
    H,d,e,f,h,z,zz=s.symbols("H d e f h z zz")
    V0,V1,V2,F0,F1,alpha=s.symbols("V0 V1 V2 F0 F1 alpha")
    D=lambda expr: (s.diff(expr,H)*d+s.diff(expr,d)*e+s.diff(expr,e)*f
                    +s.diff(expr,h)*z+s.diff(expr,z)*zz)
    R=6*(d+2*H**2)
    boxR=-6*f-24*d**2-42*H*e-72*H**2*d
    boxx=-zz-3*H*z
    riemann_minus_ricci=-12*H**2*(d+H**2)
    rhoR=36*d**2-216*H**2*d-72*H*e
    pR=108*d**2+216*H**2*d+144*H*e+24*f
    V=V0+V1*h+V2*h**2/2
    F=F0+F1*h
    rho=V-6*F*H**2-6*H*F1*z+alpha*rhoR
    pressure=-V+2*F*(2*d+3*H**2)+2*F1*zz+4*H*F1*z+alpha*pR
    Q=2*(V1+V2*h)-2*F1*R
    assert s.expand(D(rho)+3*H*(rho+pressure)-z*Q/2)==0
    assert s.expand(-rhoR+3*pR+12*boxR)==0
    # Covariant heat-kernel invariants independently expanded to cosmic-time jets.
    a2=(h-R/6)**2/2+riemann_minus_ricci/180+boxR/30-boxx/6
    claimed=(58*H**4-14*H**2*d-60*H**2*h-42*H*e+15*H*z
             -9*d**2-30*d*h-6*f+15*h**2+5*zz)/30
    assert s.expand(a2-claimed)==0
    assert s.expand(a2.subs({H:0,d:0,e:0,f:0})-h**2/2-zz/6)==0
    # Exact finite-cutoff primitives after p -> sqrt(r) t/sqrt(1-t^2).
    # This also avoids cancellation-prone subtraction of infinite moments.
    t=s.symbols("t",real=True)
    moment_checks=[]
    for n in range(5,16,2):
        m=(n-5)//2
        polynomial=sum((-1)**j*s.binomial(m,j)*t**(2*j+3)/s.Integer(2*j+3)
                       for j in range(m+1))
        assert s.expand(s.diff(polynomial,t)-t**2*(1-t**2)**m)==0
        beta_limit=s.sqrt(s.pi)*s.gamma(s.Rational(n-3,2))/(4*s.gamma(s.Rational(n,2)))
        assert s.simplify(polynomial.subs(t,1)-beta_limit)==0
        moment_checks.append(n)
    x,r=s.symbols("x r",positive=True)
    Vflat=(x**2*s.log(x/r)-s.Rational(3,2)*x**2+2*r*x-r**2/2)/(64*s.pi**2)
    Qflat=2*s.diff(Vflat,x)
    assert s.simplify(-4*Vflat+x*Qflat-(x-r)**2/(32*s.pi**2))==0
    return dict(status="PASS",scope="Exact symbolic identities; no numerical field evolution",
                action_exchange=True,regular_H_zero_pressure=True,
                curvature_squared_trace_shift="-12 alpha box R",
                heat_kernel_FRW_expansion=True,flat_trace_reference_shift=True,
                exact_finite_cutoff_moment_powers=moment_checks,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=="__main__":
    print(json.dumps(verify(),sort_keys=True))
