#!/usr/bin/env python3
"""Exact five-label implicit-event and gauge-jet check at a background event.
This checks the local Jacobian algebra, not existence of the physical collar.
"""
from fractions import Fraction as F
import json

def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def mv(M,v): return [sum(M[i][j]*v[j] for j in range(5)) for i in range(5)]
def plus(a,b): return [a[i]+b[i] for i in range(5)]
def minus(a,b): return [a[i]-b[i] for i in range(5)]
def scale(a,s): return [s*x for x in a]
def jac(p):return [[(p if i==0 else F(1)) if i==j else F(0) for j in range(5)] for i in range(5)]
def implicit(p,z):
    check(p!=0,'degenerate scalar background cannot be inverted')
    return [-z[0]/p]+[-z[i] for i in range(1,5)]
# Labels are (phi,tau,alpha1,alpha2,alpha3); coordinates are (y,t,a1,a2,a3).
p=F(7,3); J=jac(p)
z=[F(5,7),F(-2,3),F(4,5),F(-8,9),F(11,12)]
xi=[F(-3,2),F(2,7),F(1,4),F(-5,8),F(7,6)]
x=implicit(p,z)
check(plus(mv(J,x),z)==[0]*5,'implicit first event equation')
z_new=plus(z,mv(J,xi))
x_new=implicit(p,z_new)
check(x_new==minus(x,xi),'transformed event must subtract generator')
check(plus(mv(J,x_new),z_new)==[0]*5,'transformed label equation')
check(x_new!=plus(x,xi),'wrong-sign generator must fail')
# p=0: the first row of the Jacobian vanishes, so scalar label can obstruct
# any solution; when z_phi=0, every radial displacement is in the kernel.
J0=jac(F(0)); bad=plus(mv(J0,x),z)
check(bad[0]==z[0] and bad[0]!=0,'degenerate inconsistent scalar label')
z_good=[F(0)]+z[1:]
x_a=[F(0)]+[-v for v in z_good[1:]]
x_b=[F(11,4)]+x_a[1:]
check(plus(mv(J0,x_a),z_good)==[0]*5,'degenerate compatible event A')
check(plus(mv(J0,x_b),z_good)==[0]*5 and x_a!=x_b,'degenerate non-unique event B')
print(json.dumps({'status':'PASS','coordinate_order':['y','t','a1','a2','a3'],'label_order':['phi','tau','alpha1','alpha2','alpha3'],'p':str(p),'Z1':list(map(str,z)),'xi':list(map(str,xi)),'X1':list(map(str,x)),'Z1_prime':list(map(str,z_new)),'X1_prime':list(map(str,x_new)),'p_zero_examples':['nonzero scalar perturbation is inconsistent at linear order','zero scalar perturbation leaves radial displacement undetermined'],'scope':'local first-jet Jacobian algebra; p=0 only shows failure of this inverse-function argument, not all possible anchoring schemes'},sort_keys=True))
