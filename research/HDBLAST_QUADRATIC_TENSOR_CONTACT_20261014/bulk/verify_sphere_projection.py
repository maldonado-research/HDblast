#!/usr/bin/env python3
"""Direct Hopf-coordinate verification of the complete scalar-gradient TT split."""
import argparse
import json
from pathlib import Path
import sympy as s

ap=argparse.ArgumentParser()
ap.add_argument('--output',type=Path)
args=ap.parse_args()
th,ph,ps=s.symbols('th ph ps',real=True)
coords=[th,ph,ps]
sn,cs=s.sin(th),s.cos(th)
g=[s.Integer(1),sn**2,cs**2]
T=[s.Rational(1,6),sn**2*(s.Rational(1,6)-sn**2/2),
   cs**2*(s.Rational(1,6)-cs**2/2)]
Y=sn**2*s.cos(2*ph)
F2=sn**2-s.Rational(1,2)
F4=Y**2-sn**2/2+s.Rational(1,12)
G={}
checks={}

def zero(name,v):
    v=s.simplify(s.trigsimp(s.expand_trig(v)))
    if v != 0:
        v=s.simplify(v.rewrite(s.exp))
    if v != 0:
        raise RuntimeError(f'{name}: {v}')
    checks[name]='passed'

for i in range(3):
    for j in range(3):
        for k in range(3):
            G[k,i,j]=s.simplify(((s.diff(g[k],coords[i]) if k==j else 0)
                                 +(s.diff(g[k],coords[j]) if k==i else 0)
                                 -(s.diff(g[i],coords[k]) if i==j else 0))/(2*g[k]))

def Hess(f,i,j):
    return s.diff(f,coords[i],coords[j])-sum(G[k,i,j]*s.diff(f,coords[k]) for k in range(3))

Q={}
for i in range(3):
    for j in range(i,3):
        Q[i,j]=(s.diff(Y,coords[i])*s.diff(Y,coords[j])-Hess(F2,i,j)/2
                -Hess(F4,i,j)/3+(4*T[i]/3 if i==j else 0))
        if i!=j:
            zero(f'TT_split_offdiagonal_{i}{j}',Q[i,j])
for i in [1,2]:
    zero(f'TT_split_trace_difference_{i}',Q[i,i]/g[i]-Q[0,0])

# x=sin(theta)^2 is uniform under the normalized Hopf sphere measure.
x=s.symbols('x',real=True)
norm=s.Rational(1,36)+(s.Rational(1,6)-x/2)**2+(s.Rational(1,6)-(1-x)/2)**2
stress_moment=s.Rational(2,3)*x-s.Rational(4,3)*x**2
zero('tensor_norm',s.integrate(norm,(x,0,1))-s.Rational(1,12))
zero('stress_contraction',s.integrate(stress_moment,(x,0,1))+s.Rational(1,9))
zero('projection_coefficient',s.integrate(stress_moment,(x,0,1))/s.integrate(norm,(x,0,1))+s.Rational(4,3))
out={'status':'passed','checks':checks,'TT_gradient_coefficient':'-4/3',
     'averaged_tensor_norm':'1/12','averaged_tensor_stress_contraction':'-1/9'}
if args.output:
    args.output.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
