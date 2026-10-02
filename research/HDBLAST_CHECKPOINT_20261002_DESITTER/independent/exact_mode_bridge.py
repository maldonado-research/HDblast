#!/usr/bin/env python3
"""Exact x=r=2H² plane-wave check of variance, energy and pressure.

Analytic diagnostic added after the proper-time registration. No evolution or
new physical parameter run. H=1 below sets units; results scale as H² or H⁴.
The WKB frequency is derived directly from its Riccati equation and each
stress component is integrated separately, without using the trace identity.
"""
import argparse
import json
from pathlib import Path
import sympy as s


def require(value, message):
    if not value:
        raise RuntimeError(message)


eta=s.symbols("eta",negative=True)
k,w=s.symbols("k w",positive=True)
v=s.symbols("v",nonnegative=True)
omega=s.sqrt(k*k+2/eta**2)
delta=-2/eta**2
u2=delta/(2*omega)-s.diff(omega,eta,2)/(4*omega**2)+3*s.diff(omega,eta)**2/(8*omega**3)
u4=(-u2**2/(2*omega)-s.diff(u2,eta,2)/(4*omega**2)
    +s.diff(omega,eta,2)*u2/(4*omega**3)
    +3*s.diff(omega,eta)*s.diff(u2,eta)/(4*omega**3)
    -3*s.diff(omega,eta)**2*u2/(4*omega**4))

def at_unit_scale(expr):
    return s.simplify(expr.subs(eta,-1).subs(k**2,w*w-2))

u,up,uu,wp=map(at_unit_scale,[u2,s.diff(u2,eta),u4,s.diff(omega,eta)])
m1=1+wp/(2*w)
m3=up/(2*w)-wp*u/(2*w*w)
shear2=m1*m1/w
shear4=2*m1*m3/w-m1*m1*u/w**2
c=w*w/3+s.Rational(4,3)
rho_sub=w/2+shear2/4+(u*u/w+shear4)/4
p_sub=(w-c/w)/4+(u+c*u/w**2+shear2)/4+(uu+c*uu/w**2-c*u*u/w**3+shear4)/4
q_sub=1/(2*w)-u/(2*w*w)
bare={"rho":k/2+s.Rational(3,4)/k,
      "p":k/6-s.Rational(1,4)/k,"Q":1/(2*k)}
sub={"rho":rho_sub,"p":p_sub,"Q":q_sub}
predicted={"rho":s.Rational(11,960),"p":-s.Rational(11,960),"Q":s.Rational(1,12)}
answers={}
for name in bare:
    # p=sqrt(2)*v/sqrt(1-v²), dp=sqrt(2)/(1-v²)^(3/2).
    # Include momentum measure p²dp/(2 pi²). Store answer times pi².
    integrand=s.factor((bare[name]-sub[name]).subs({k:s.sqrt(2)*v/s.sqrt(1-v*v),w:s.sqrt(2)/s.sqrt(1-v*v)})
                       *s.sqrt(2)*v*v/(1-v*v)**s.Rational(5,2))
    integrand=s.cancel(integrand)
    result=s.simplify(s.integrate(integrand,(v,0,1)))
    require(result==predicted[name],f"Direct {name} mismatch: {result}")
    answers[name]={"integrand_times_pi_squared":str(integrand),
                   "integral_times_pi_squared":str(result)}
document={"status":"PASS","derivation":"Separate exact-mode minus order-four stress/order-two variance integrals",
          "state":"minimal scalar x=r=2H², Euclidean/Bunch-Davies exact conformal-frequency mode",
          "results":answers,"pressure_inferred_from_trace":False,
          "negative_pressure_omission_detected":predicted["p"]!=0,
          "symbolic_only":True,"sympy_version":s.__version__}
parser=argparse.ArgumentParser()
parser.add_argument("--output",type=Path,default=Path(__file__).with_name("EXACT_MODE_BRIDGE.json"))
path=parser.parse_args().output
path.write_text(json.dumps(document,indent=2)+"\n")
print(json.dumps(document,indent=2))
