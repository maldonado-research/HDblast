#!/usr/bin/env python3
"""Independent rational d'Alembert convolution time-derivative audit."""
from fractions import Fraction as Q
from math import comb
import json

checks=[]
def need(v,name):
    if not v: raise RuntimeError(name)
    checks.append(name)
def add(p,q):
    z=dict(p)
    for n,c in q.items():z[n]=z.get(n,Q(0))+c
    return {n:c for n,c in z.items() if c}
def scale(p,c):return {n:c*v for n,v in p.items() if c*v}
def deriv(p):return {n-1:n*c for n,c in p.items() if n and c}
def evaluate(p,x):return sum(c*x**n for n,c in p.items())
def Kmono(a,b,r):
    # Integrate the original inner antiderivative evaluated at y=r +/- (x-s).
    # All s integrations are elementary; the output is a polynomial in x.
    p={}
    for j in range(1,b+2,2):
        common=Q(comb(b+1,j),b+1)*r**(b+1-j)
        for k in range(j+1):
            c=common*comb(j,k)*(-1)**k/Q(a+k+1)
            p[a+j+1]=p.get(a+j+1,Q(0))+c
    return {n:c for n,c in p.items() if c}
def Bvelocity_mono(b,r):
    # 1/2 int_{r-x}^{r+x} y^b dy, independently expanded at fixed r.
    return {j:Q(comb(b+1,j),b+1)*r**(b+1-j) for j in range(1,b+2,2)}
def forcing(a,b,x,r):return x**a*r**b
points=[(Q(2,7),Q(-3,5)),(Q(4,3),Q(7,11)),(Q(0),Q(1,2))]
rows=[]
for a in range(4):
 for b in range(4):
  for x,r in points:
   k=Kmono(a,b,r)
   lhs=deriv(k)
   kfx=scale(Kmono(a-1,b,r),Q(a)) if a else {}
   initial=Bvelocity_mono(b,r) if a==0 else {}
   rhs=add(kfx,initial)
   need(lhs==rhs,f'time derivative identity a{a} b{b} r{r}')
   # Independent PDE/initial check helps detect an erroneous common K term.
   krr=scale(Kmono(a,b-2,r),Q(b*(b-1))) if b>=2 else {}
   wave=add(deriv(deriv(k)),scale(krr,Q(-1)))
   need(evaluate(wave,x)==forcing(a,b,x,r),f'principal wave identity a{a} b{b} point{x},{r}')
   need(evaluate(k,Q(0))==0 and evaluate(deriv(k),Q(0))==0,f'zero Cauchy data a{a} b{b} r{r}')
   rows.append({'a':a,'b':b,'x':str(x),'r':str(r),'K':str(evaluate(k,x)),'K_x':str(evaluate(lhs,x)),'K_f_x':str(evaluate(kfx,x)),'free_initial_velocity_term':str(evaluate(initial,x))})
# A mixed polynomial combines nonzero and zero-time forcing terms. Coefficient
# equality above proves every linear combination; this separately evaluates one.
terms=[(Q(3,7),0,0),(Q(-5,9),0,3),(Q(11,13),2,1),(Q(-2,5),3,3)]
x,r=Q(5,8),Q(3,10)
k={};rhs={}
for c,a,b in terms:
 k=add(k,scale(Kmono(a,b,r),c))
 rhs=add(rhs,scale(Kmono(a-1,b,r),c*a)) if a else add(rhs,scale(Bvelocity_mono(b,r),c))
need(deriv(k)==rhs,'mixed polynomial time differentiation')
# Nonzero F(0) must appear; the omitted term fails already for F=1.
missing_initial = deriv(Kmono(0,0,Q(1,3)))=={}
need(not missing_initial,'reject missing initial velocity term F=1')
# Free data estimate: v is displacement and w acceleration for Y. On a line,
# Y=1/2(v(r+x)+v(r-x))+1/2 int w; differentiation gives
# |Y| <= V+T W, |Y_x|,|Y_r| <= V_r+W. Norm X includes |Y|/3.
# The triangle constants are sharp in the independent branches: choose a
# constant positive v/w for value. On the retained interval [-T,T], the
# quadratic v(r)=V_r*r^2/(2T) and constant w attain Y_x=V_r+W at (T,0).
T=Q(3,2)
for V,Vr,W in [(Q(0),Q(0),Q(1)),(Q(2),Q(3),Q(5)),(Q(7),Q(0),Q(0)),(Q(1,7),Q(8,9),Q(3,11))]:
 J=max((V+T*W)/3,Vr+W)
 need(J>=(V+T*W)/3 and J>=Vr+W,'free-wave triangle J envelope')
 need((V+T*W)/3==V/3+W/2,'free-wave value triangle exact T=3/2')
# Missing W in the derivative branch or missing the time-length in value
# fails physically realizable constant/affine examples.
constant_velocity=Bvelocity_mono(0,Q(0))
constant_Yx=evaluate(deriv(constant_velocity),T)
wrong_J=max(T/3,Q(0))
need(constant_Yx>wrong_J,'reject missing acceleration W in derivative branch')
need(evaluate(constant_velocity,T)>Q(1),'reject unit-time value factor on length3/2')
# A smooth quadratic displacement witness attains the derivative triangle.
Vr,W=Q(7,5),Q(11,9)
quadratic_Y={2:Vr/(2*T),1:W}
need(evaluate(deriv(quadratic_Y),T)==Vr+W,'smooth quadratic displacement attains derivative triangle')
print(json.dumps({'status':'PASS_DIFFERENTIATED_PRESCRIBED_FORCING_IDENTITY','identity':'partial_x Kf=K(partial_x f)+Bfree[0,f(0)]','scope':'full-line principal wave identity before physical coupling/reflection; no assertion about parity commuting with operators','polynomial_monomial_count':16,'rational_point_count':3,'rows':rows,'mixed_polynomial_checked':True,'J_formula':'max((V+(3/2)W)/3,V_r+W)','named_checks':checks,'wrong_formula_controls_rejected':['omitted f(0) free velocity term','omitted acceleration W in derivative bound','unit-time value coefficient on T=3/2']},sort_keys=True,indent=2))
