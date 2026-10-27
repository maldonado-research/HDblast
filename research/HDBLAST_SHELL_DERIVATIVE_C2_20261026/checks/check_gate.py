#!/usr/bin/env python3
"""Rational synthetic check of the differentiated mixed-boundary response gate."""
from fractions import Fraction as F
from math import comb
import json

def req(v,m):
    if not v: raise RuntimeError(m)
def add(*ps):
    o={}
    for p in ps:
        for k,v in p.items():o[k]=o.get(k,F(0))+v
    return {k:v for k,v in o.items() if v}
def scale(p,c):return {k:F(c)*v for k,v in p.items() if c*v}
def mul(p,q):
    o={}
    for (i,j),a in p.items():
        for (k,l),b in q.items():o[i+k,j+l]=o.get((i+k,j+l),F(0))+a*b
    return {k:v for k,v in o.items() if v}
def der(p,axis):
    o={}
    for (i,j),v in p.items():
        n=(i,j)[axis]
        if n:o[(i-1,j) if axis==0 else (i,j-1)]=v*n
    return o
def ev(p,x,r):return sum((v*x**i*r**j for (i,j),v in p.items()),F(0))
def pw(p,n):
    out={(0,0):F(1)}
    for _ in range(n):out=mul(out,p)
    return out
x={(1,0):F(1)}; rr={(0,1):F(1)}; one={(0,0):F(1)}; z=add(x,rr,scale(one,-1))
e=pw(z,4);ex=der(e,0);er=der(e,1);exx=der(ex,0);err=der(er,1);v=ex
alpha=F(1,15000); beta=F(1,15000); delta=F(1,1500000); robin=F(1,2000)
A=scale(x,alpha);B=scale(x,beta);D=scale(x,delta)
f=add(exx,scale(err,-1),scale(mul(A,ex),-1),scale(mul(B,er),-1),scale(mul(D,e),-1))
fx=der(f,0)
newf=add(fx,mul(der(A,0),ex),mul(der(B,0),er),mul(der(D,0),e))
vd=add(der(der(v,0),0),scale(der(der(v,1),1),-1),scale(mul(A,der(v,0)),-1),scale(mul(B,der(v,1)),-1),scale(mul(D,v),-1))
req(newf==vd,'differentiated PDE')
wrong=add(fx,mul(der(B,0),er),mul(der(D,0),e))
req(wrong!=vd and ev(add(vd,scale(wrong,-1)),F(1),F(3,4))!=0,'A_x omission control')
g=add(er,scale(e,-robin));gx=der(g,0)
# Restrict each side to shell by evaluating r=0 for x>1.
for t in (F(9,8),F(5,4)):
    lhs=ev(der(v,1),t,F(0))
    rhs=robin*ev(v,t,F(0))+ev(gx,t,F(0))
    req(lhs==rhs,'differentiated Robin')
# Independent initial-acceleration control with a nonzero initial f.
e2=mul(add(one,x),e)
e2x=der(e2,0);e2r=der(e2,1);e2xx=der(e2x,0);e2rr=der(e2r,1)
f2=add(e2xx,scale(e2rr,-1),scale(mul(A,e2x),-1),scale(mul(B,e2r),-1),scale(mul(D,e2),-1))
init_rhs=add(e2rr,mul(A,e2x),mul(B,e2r),mul(D,e2),f2)
for r0 in (F(5,4),F(3,2),F(7,4)):
    req(ev(e2xx,F(0),r0)==ev(init_rhs,F(0),r0),'initial acceleration')
    req(ev(f2,F(0),r0)!=0 and ev(e2xx,F(0),r0)!=ev(e2rr,F(0),r0),'missing initial f control')
# C3 zero-wedge e=(x+r−1)_+^4; the interface derivatives through order 3 vanish.
for dx in range(4):
    p=e
    for _ in range(dx):p=der(p,0)
    for s in (F(1,4),F(1,2),F(3,4)):
        req(ev(p,s,1-s)==0,'zero wedge join')
# Reflected d'Alembert test with nontrivial compatible initial displacement.
def initial(s,odd=False):
    z=max(abs(s)-1,F(0))
    return z**4 * ((-1 if s<0 else 1) if odd else 1)
def initial_prime(s,odd=False):
    z=max(abs(s)-1,F(0))
    return 4*z**3 * (1 if odd else (-1 if s<0 else 1))
t=F(5,4)
free=lambda odd:(initial(t,odd)+initial(-t,odd))/2
free_r=lambda odd:(initial_prime(t,odd)+initial_prime(-t,odd))/2
req(free(True)==0 and free_r(False)==0,'odd Dirichlet/even Neumann reflection')
req(free(False)!=0 and free_r(True)!=0,'swapped parity negative controls')
# D+ has z<=3/4. These bounds are exact, with z=3/4 attained.
T=F(3,2); H=F(3,4); E=4*H**3; actual=12*H**2
req(H**4<=F(3,4)*E and E==F(27,16) and actual==F(27,4),'X(e) and X(v)')
a0=alpha*T;b0=beta*T;d0=delta*T
# e_x(0,r)=4(r−1)_+^3 and v_x(0,r)=12(r−1)_+^2.
v0=4*H**3;v0r=12*H**2;w0=12*H**2
J=max((v0+T*w0)/3,v0r+w0)
# f = −4(alpha+beta)x z³−delta x z⁴; bound |f| and |f_x| by triangle inequality.
F0=4*(alpha+beta)*T*H**3+delta*T*H**4
F1=4*(alpha+beta)*(H**3+3*T*H**2)+delta*(H**4+4*T*H**3)
G1=12*F(1,4)**2+4*robin*F(1,4)**3
q=F(3409,8000000)
Z=(J+T*(F1+(alpha+beta+F(3,4)*delta)*E)+G1)/(1-q)
req(0<q<1 and Z>=actual,'differentiated contraction gate')
req(Z+(a0+b0+F(3,4)*d0)*E+F0>=actual,'radial recovery gate')
print(json.dumps({'status':'PASS','example':'C3 zero-wedge scalar with nonzero time-dependent A,B,D0 and Robin feedback','X_error':str(E),'actual_X_time_derivative':str(actual),'free_wave_gate':str(J),'C2_gate_Z':str(Z),'prior_contraction_q':str(q),'controls':['omitted A_x term rejected','differentiated Robin sign','nonzero initial forcing in acceleration','mixed reflection parity','C3 wedge join','full radial recovery'],'scope':'exact synthetic algebra and bound; no actual HDBLAST numerical C2 estimate'},sort_keys=True))
