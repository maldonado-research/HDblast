"""Exact independent grading expansion and finite-cutoff trace integrand check.

No primary source imports, physical mode solves, or numerical quadrature.
"""
import hashlib
import json
from pathlib import Path

import sympy as s


def verify():
    eps=s.symbols("eps")
    w,r,h,H,d,e,f,z,zz=s.symbols("w r h H d e f z zz", nonzero=True)
    U,V,WP,UP=s.symbols("U V WP UP")
    W=w+eps**2*U+eps**4*V
    Wprime=eps*WP+eps**3*UP
    geometric=eps*H*Wprime/W**2+Wprime**2/(4*W**3)
    B=-(w**2+2*r)/3
    rho_formal=(W+(w**2+eps**2*(h+H**2))/W+geometric)/4
    p_formal=(W+(B+eps**2*(H**2-h))/W+geometric)/4
    rho=s.series(rho_formal,eps,0,5).removeO().subs(eps,1)
    pressure=s.series(p_formal,eps,0,5).removeO().subs(eps,1)
    D=lambda v: (s.diff(v,w)*(-H*(w*w-r)/w)+s.diff(v,H)*d
                 +s.diff(v,d)*e+s.diff(v,e)*f+s.diff(v,h)*z+s.diff(v,z)*zz)
    # W2=a*u; W2'=a^2*(H*u+D*u); W2''=a^3*(2H*up+D*up).
    wp=r*H/w
    wpp=r*(d+3*H**2)/w-r**2*H**2/w**3
    u=(h-d-2*H**2)/(2*w)-wpp/(4*w**2)+3*wp**2/(8*w**3)
    up=H*u+D(u)
    upp=2*H*up+D(up)
    v=(-u**2/(2*w)-upp/(4*w**2)+wpp*u/(4*w**3)
       +3*wp*up/(4*w**3)-3*wp**2*u/(4*w**4))
    replace={U:u,V:v,WP:wp,UP:up}
    rho=rho.subs(replace)
    pressure=pressure.subs(replace)
    q=1/(2*w)-u/(2*w**2)
    kinetic=(D(D(q))-3*H*D(q)-3*d*q)/2
    anomaly=s.expand(kinetic-(h+r)*q+rho-3*pressure)
    # Registered finite-K coefficient representation, independently checked here.
    coefficients={
        5:r*(-18*H**2*d-6*H**2*h-9*H*e+5*H*z-3*d**2-4*d*h-f+3*h**2+zz)/16,
        7:-r**2*(-40*H**4-24*H**2*d+50*H**2*h+5*H*e+10*H*z+10*d*h+f)/32,
        9:7*r**3*(55*H**4+48*H**2*d+10*H**2*h+4*H*e+3*d**2)/64,
        11:-231*H**2*r**4*(3*H**2+d)/64,
        13:1155*H**4*r**5/256,
    }
    assert s.expand(anomaly-sum(c/w**n for n,c in coefficients.items()))==0
    continuum=sum(c*s.sqrt(s.pi)*s.gamma(s.Rational(n-3,2))
                  *r**s.Rational(3-n,2)/(4*s.gamma(s.Rational(n,2)))
                  for n,c in coefficients.items())
    curvature=6*(d+2*H**2)
    boxR=-6*f-24*d**2-42*H*e-72*H**2*d
    a2=(h-curvature/6)**2/2-H**2*(d+H**2)/15+boxR/30+(zz+3*H*z)/6
    assert s.simplify(continuum-a2/8)==0
    # Full fourth-order stress is required: deliberate p4 omission is detected.
    p4=s.expand(pressure).coeff(f)
    assert p4 != 0
    return dict(status="PASS",scope="Exact formal algebra, no field evolution",
                complete_graded_WKB_expansion=True,
                finite_cutoff_trace_coefficients=True,
                covariant_continuum_anomaly=True,
                fourth_order_pressure_required=True,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=="__main__":
    print(json.dumps(verify(),sort_keys=True))
