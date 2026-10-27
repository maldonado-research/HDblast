#!/usr/bin/env python3
"""Exact scalar-jet tests of gauge changes and an identified moving event.
A toy finite-dimensional model, not a spacetime field-equation proof.
"""
from fractions import Fraction as F
import json

def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def j(a=0,b=0,c=0): return (F(a),F(b),F(c)) # coefficients of 1, λ, λ²
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def neg(a): return tuple(-x for x in a)
def sub(a,b): return add(a,neg(b))
def mul(a,b): return tuple(sum(a[k]*b[i-k] for k in range(i+1)) for i in range(3))
def s(a,k): return tuple(F(k)*x for x in a)
L=j(0,1); HALF_L2=j(0,0,F(1,2))
def vals(v):return tuple(v[i][0] for i in range(2))
def xi(x,y): return (add(j(1),sub(s(x,2),y)),add(j(-2),add(x,y)))
def eta(x,y): return (sub(j(3),x),add(j(1),s(y,2)))
def jac_xi(v): return (2*v[0]-v[1], v[0]+v[1])
def liephi(x,y):
    xx=xi(x,y); ee=eta(x,y); dx=(sub(s(xx[0],2),xx[1]),add(xx[0],xx[1]))
    return (add(x,add(mul(L,xx[0]),mul(HALF_L2,add(ee[0],dx[0])))), add(y,add(mul(L,xx[1]),mul(HALF_L2,add(ee[1],dx[1])))))
def f0(x,y):return add(add(mul(x,x),mul(x,y)),add(s(mul(y,y),2),sub(s(x,3),y)))
def f1(x,y):return add(sub(add(mul(x,x),s(mul(x,y),2)),s(y,3)),j(4))
def f2(x,y):return add(sub(s(x,3),s(y,2)),mul(x,y))
def full(x,y):return add(f0(x,y),add(mul(L,f1(x,y)),mul(HALF_L2,f2(x,y))))
q0=(F(2),F(-1)); v=(F(3,2),F(-2,3)); w=(F(4,3),F(5,2))
xi0=vals(xi(j(q0[0]),j(q0[1]))); eta0=vals(eta(j(q0[0]),j(q0[1])))
q1=tuple(v[i]-xi0[i] for i in range(2)); dxq1=jac_xi(q1); dxxi=jac_xi(xi0)
q2=tuple(w[i]-2*dxq1[i]-eta0[i]-dxxi[i] for i in range(2))
p=(j(q0[0],v[0],w[0]/2),j(q0[1],v[1],w[1]/2))
q=(j(q0[0],q1[0],q2[0]/2),j(q0[1],q1[1],q2[1]/2))
check(liephi(*q)==p,'inverse event jets fail')
check(full(*liephi(*q))==full(*p),'anchored scalar must agree')
# The test checks first and second diffeomorphism generators directly.
x,y=map(j,q0); transformed=full(*liephi(x,y)); s0=f0(x,y)[0]; s1=f1(x,y)[0]; s2=f2(x,y)[0]
g0=(2*q0[0]+q0[1]+3,q0[0]+4*q0[1]-1); g1=(2*q0[0]+2*q0[1],2*q0[0]-3)
dot=lambda a,b:sum(a[i]*b[i] for i in range(2))
hessquad=2*xi0[0]**2+2*xi0[0]*xi0[1]+4*xi0[1]**2
lie1=dot(g0,xi0); lie2_s0=dot(g0,eta0)+dot(g0,dxxi)+hessquad
expected=(s0,s1+lie1,(s2+2*dot(g1,xi0)+lie2_s0)/2)
check(transformed==expected,'second order pullback / 1/2 convention')
check(dot(g0,eta0)!=0 and dot(g0,dxxi)!=0 and hessquad!=0 and dot(g1,xi0)!=0,'nonzero generator controls')
# With S0=0 and S1=x-x0, S1 vanishes at the reference event but not as a field.
# Changing the event without a matching pullback moves the second coefficient.
def h(x,y):return add(mul(L,sub(x,j(q0[0]))),mul(HALF_L2,add(j(1),s(y,2))))
original=h(*p); moving_unmatched=h(*q)
check(original[1]==0,'first-order value at base event vanishes')
check(original[2]!=moving_unmatched[2],'moving-event negative control')
check(h(*liephi(*q))==original,'matched moving event')
print(json.dumps({'status':'PASS','normalization':'lambda^2/2 for second derivatives', 'coordinate_gauge_jets':'S1+Lxi S0; S2+2 Lxi S1+(Leta+Lxi^2)S0','xi1':list(map(str,xi0)),'eta2':list(map(str,eta0)),'old_event_first_jet':list(map(str,v)),'new_event_first_jet':list(map(str,q1)),'nonzero_second_generator_term':str(dot(g0,eta0)),'nonzero_first_generator_squared_term':str(dot(g0,dxxi)+hessquad),'nonzero_first_field_cross_term':str(2*dot(g1,xi0)),'matched_anchor_second_coefficient':str(original[2]),'unmatched_anchor_second_coefficient':str(moving_unmatched[2]),'scope':'toy scalar fields and a two-dimensional diffeomorphism; no Einstein-scalar field equations'},sort_keys=True))
