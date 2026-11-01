#!/usr/bin/env python3
"""Independent exact jet algebra. No source/producer imports; no numeric oracle."""
from fractions import Fraction as F
import json
from pathlib import Path
import sympy as s

e=s.symbols('e')
def tr2(x):
    return s.expand(x).series(e,0,3).removeO().expand()
def need(x, label):
    if not x:
        raise RuntimeError(label)

# A fully generic symmetric first/second Lorentzian metric perturbation.
H0=s.diag(-1,1,1,1)
H1=s.zeros(4); H2=s.zeros(4)
for i in range(4):
    for j in range(i,4):
        H1[i,j]=H1[j,i]=s.Symbol(f'h{i}{j}')
        H2[i,j]=H2[j,i]=s.Symbol(f'k{i}{j}')
G0=H0
G1=-G0*H1*G0
G2=2*G0*H1*G0*H1*G0-G0*H2*G0
need(H0*G0==s.eye(4),'inverse at order zero')
need(H0*G1+H1*G0==s.zeros(4),'inverse at order one')
need(s.simplify(H0*G2+2*H1*G1+H2*G0)==s.zeros(4),'inverse at order two')
a=s.Matrix(s.symbols('a0:4')); b=s.Matrix(s.symbols('b0:4'))
G=G0+e*G1+e**2*G2/2
dphi=e*a+e**2*b/2
# No need to expand high powers: coefficient extraction is exact polynomial algebra.
q=(dphi.T*G*dphi)[0]
need(s.expand(q).coeff(e,0)==0,'Q0')
need(s.expand(q).coeff(e,1)==0,'Q1')
need(s.expand(s.expand(q).coeff(e,2)-(a.T*G0*a)[0])==0,'Q2 completion independence')

u0=s.Matrix(s.symbols('u0:4')); u1=s.Matrix(s.symbols('w0:4'))
u=u0+e*u1
k=(u.dot(dphi))**2
need(s.expand(s.expand(k).coeff(e,2)-u0.dot(a)**2)==0,'K2 observer variation independence')
sp=(dphi.T*(G+u*u.T)*dphi)[0]
need(s.expand(s.expand(sp).coeff(e,2)-(a.T*G0*a)[0]-u0.dot(a)**2)==0,'S2 projector')

# Generic event pullback and normalized moving measure.
x,v,w,q20,q21,q22,n0,n1,n2,m0,m1,m2=s.symbols('x v w q20 q21 q22 n0 n1 n2 m0 m1 m2',nonzero=True)
x_event=x+e*v+e**2*w/2
q_event=e**2*(q20+q21*(x_event-x)+q22*(x_event-x)**2/2)/2
need(tr2(q_event)==e**2*q20/2,'event variations')
# N(e)=integral Q(e) dmu(e) starts at order two; D(e)=integral dmu(e).
mean=tr2((e**2*n2/2)/(m0+e*m1+e**2*m2/2))
need(s.simplify(mean-e**2*n2/(2*m0))==0,'normalized measure variations')

# Independent moment rule for uniform S3: E[prod n_i^(2k_i)].
def rising(x,k):
    out=F(1)
    for j in range(k): out*=x+j
    return out
def moment(ks):
    out=F(1)
    for ki in ks: out*=rising(F(1,2),ki)
    return out/rising(F(2),sum(ks))
y2=moment([2,0,0,0])+moment([0,2,0,0])-2*moment([1,1,0,0])
z=moment([1,0,0,0])+moment([0,1,0,0])
dy2=4*(z-y2)
need(y2==F(1,6) and dy2==F(4,3),'quadrupole moments')
A,B,R,C=s.symbols('A B R C',positive=True)
mean_k=2*A**2*s.Rational(y2.numerator,y2.denominator)/R**2
mean_sp=2*B**2*s.Rational(dy2.numerator,dy2.denominator)/(R**2*C**2)
need(s.simplify(mean_k-A**2/(3*R**2))==0,'mean K2 normalization')
need(s.simplify(mean_sp-8*B**2/(3*R**2*C**2))==0,'mean S2 normalization')

# Counterexample: background observer changes K2 even with fixed d eta.
# Minkowski frame, deta=dt, u=(5/3,4/3,0,0) is future unit.
boost=s.Matrix([s.Rational(5,3),s.Rational(4,3),0,0])
need((boost.T*H0*boost)[0]==-1,'boost unit')
need(2*boost[0]**2==s.Rational(50,9),'changed background K2')
# Counterexample: background measure changes quadrupole means.
# Weight 1+Y^2, normalized: E[Y^2]_w=(E[Y^2]+E[Y^4])/(1+E[Y^2]).
y4=moment([4,0,0,0])-4*moment([3,1,0,0])+6*moment([2,2,0,0])-4*moment([1,3,0,0])+moment([0,4,0,0])
weighted=(y2+y4)/(1+y2)
need(weighted!=y2,'changed background angular measure')
# Counterexample to finite-radius uniformity: same first tangent, arbitrary second completion.
# phi_e=e*t+e^2*c*x in a flat shell gives Q=-e^2+e^4*c^2 at (t,x) arbitrary.
# For any proposed e0>0, choose c>1/e0: Q(e0)>0 although Q2=-2.
report={
 'status':'PASS',
 'scope':'independent generic exact jet algebra, angular moments and necessary-hypothesis counterexamples; no source imports or actual-jet numerical evaluation',
 'checks':['generic symmetric first/second metric inverse','Q0=Q1=0 and Q2=2h0_inverse(deta,deta)','K2=2(u0 eta)^2','S2=Q2+K2','event pullback invariance','normalized moving measure invariance with fixed background measure','uniform S3 moments','second-derivative normalization'],
 'moments':{'Y2':str(y2),'Y4':str(y4),'DY2':str(dy2),'weighted_Y2':str(weighted)},
 'counterexamples':{'different_background_observer_K2':'50/9 versus 2 for deta=dt','different_background_measure_Y2':str(weighted)+' versus '+str(y2),'no_uniform_finite_parameter_radius':'phi_e=e*t+e^2*c*x yields Q=-e^2+e^4*c^2, Q2=-2 for every c'},
}
print(json.dumps(report,indent=2,sort_keys=True))
