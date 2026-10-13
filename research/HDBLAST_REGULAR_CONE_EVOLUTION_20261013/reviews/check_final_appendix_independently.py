#!/usr/bin/env python3
"""Independent exact checks for Appendix A; no physical norm identification."""
import argparse
import json
import sympy as S

parser=argparse.ArgumentParser()
parser.add_argument('--inject-failure',action='store_true')
args=parser.parse_args()
checks={}
controls={}
def zero(name,value):
    value=S.factor(S.cancel(S.together(value)))
    if value != 0:
        raise RuntimeError(f'{name}: {value}')
    checks[name]='passed'
def reject(name,value):
    value=S.factor(S.cancel(S.together(value)))
    if value == 0:
        raise RuntimeError(f'{name}: wrong formula escaped')
    controls[name]=str(value)

r,rho,pbg,B,a,nF,F,boxh=S.symbols('r rho p_background B a nF F Box_h_F',nonzero=True)
alpha=1/(rho**2*pbg**2)
beta=rho**4*pbg**2*B
co=1/(3*a*rho**2)
Py=r**4*(nF+4*F/rho)
DP=r**4*(rho**2*boxh+4*F)
even=(alpha*Py-DP/beta)/(alpha*r**4)
odd=(alpha*Py+co*r**4*F)/(alpha*r**4)
zero('even_shell_conversion',even-(nF-boxh/B+(4/rho-4/(rho**2*B))*F))
zero('odd_shell_conversion',odd-(nF+(4/rho+pbg**2/(3*a))*F))
reject('wrong_even_tangential_sign',even-(nF+boxh/B+(4/rho-4/(rho**2*B))*F))

# Coordinate conversion of inherited physical time current on T=t=0:
# dy=Omega dr, dt-derivative=r dT-derivative, rho=Omega r.
O,f=S.symbols('Omega f',nonzero=True)
measure=O*(O*r)**2*r
zero('actual_time_current_measure',measure-O**3*r**3)
scalar_coefficient=measure*(-3/(2*O**2*f))**2
metric_coefficient=6*measure*(r**2/O**2)**2
zero('physical_scalar_weight',scalar_coefficient-9*r**3/(4*O*f**2))
zero('physical_metric_weight',metric_coefficient-6*r**7/O)
reject('wrong_metric_radial_power',metric_coefficient-6*r**5/O)

# Check the crucial derivative of E on an arbitrary function, without a
# field equation or a separated harmonic.
T=S.symbols('T')
X=S.symbols('X1:5')
fun=S.Function('F')(T,*X)
E=lambda h:T*S.diff(h,T)+sum(x*S.diff(h,x) for x in X)
zero('Euler_time_commutator',S.diff(E(fun),T)-E(S.diff(fun,T))-S.diff(fun,T))
zero('D_plus_five_velocity',S.diff(E(fun)+4*fun,T).subs(T,0)
     -sum(x*S.diff(fun,T,x).subs(T,0) for x in X)-5*S.diff(fun,T).subs(T,0))
reject('D_plus_four_is_wrong',S.diff(E(fun)+4*fun,T).subs(T,0)
     -sum(x*S.diff(fun,T,x).subs(T,0) for x in X)-4*S.diff(fun,T).subs(T,0))

# Weight derivative and exact integration-by-parts boundary identity.
z=S.symbols('s')
Oz=S.Function('Omega')(z)
fz=S.Function('f')(z)
Lz=S.diff(Oz,z)/Oz+2*S.diff(fz,z)/fz
mr=9*r**3/(4*Oz.subs(z,r**2)*fz.subs(z,r**2)**2)
zero('radial_log_weight',r*S.diff(mr,r)/mr-(3-2*r**2*Lz.subs(z,r**2)))
q=S.Function('q')(r)
Dq=r*S.diff(q,r)
raw=mr*(Dq+4*q)*(Dq+5*q)
reduced=mr*(Dq**2+(2+9*r**2*Lz.subs(z,r**2))*q**2)
boundary_derivative=S.diff(S.Rational(9,2)*mr*r*q**2,r)
zero('exact_IBP_density_identity',raw-reduced-boundary_derivative)
reject('wrong_IBP_potential',raw-(reduced+mr*q**2)-boundary_derivative)

if args.inject_failure:
    zero('injected_failure',1)
print(json.dumps({'positive_count':len(checks),'positive_checks':checks,
                  'wrong_formula_count':len(controls),'wrong_formula_controls':controls,
                  'scope':'Appendix A exact algebra; smooth endpoint decay and kernel proof checked analytically.'},indent=2,sort_keys=True))
