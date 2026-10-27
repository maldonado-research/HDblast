#!/usr/bin/env python3
"""Exact compact-kernel smoothing of frozen full23 P1 profiles, with a new candidate."""
import argparse, hashlib, importlib.util, json
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
from pathlib import Path


def need(ok, message):
    if not ok:
        raise ValueError(message)

D=F(1,2048)
N=23

def moment(k):
    if k%2: return F(0)
    return F(35,32)*sum(F(2*c,2*j+k+1) for j,c in enumerate((1,-3,3,-1)))


def hinge(z,n):
    """Convolution of (z)_+^n/n! with normalized C² compact polynomial.

    n=0 is the step function, with its point value irrelevant under integral.
    All outputs are exact for rational z, D.
    """
    z=F(z)
    need(n>=0,'negative primitive order')
    if z<=-D: return F(0)
    if z>=D:
        return sum(F(comb(n,k),factorial(n))*z**(n-k)*(-D)**k*moment(k) for k in range(n+1))
    s=z/D
    ans=F(0)
    for j,c in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            p=2*j+k
            ans += c*comb(n,k)*s**(n-k)*(-1)**k*(s**(p+1)-(-1)**(p+1))/F(p+1)
    return F(35,32)*D**n*ans/F(factorial(n))


def interior_polynomial(n):
    """Coefficients in s=z/delta, with the common delta**n included."""
    coeff=[F(0)]*(n+8)
    for j,c in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            m=2*j+k
            scale=F(35*c*comb(n,k)*(-1)**k,32*factorial(n)*(m+1))
            coeff[n-k+m+1] += scale
            coeff[n-k] -= scale*(-1)**(m+1)
    return coeff

def deriv(p): return [F(i)*p[i] for i in range(1,len(p))]
def at(p,z): return sum((a*z**i for i,a in enumerate(p)), F(0))
def trim(p):
    while p and p[-1]==0:p.pop()
    return p

# Exact polynomial identities—not finite differences—certify convolution,
# the derivative/primitive recurrence and C4 join for the smoothed hinge.
for n in range(1,5):
    lhs=trim(deriv(interior_polynomial(n)))
    rhs=trim(interior_polynomial(n-1))
    need(lhs==rhs,'exact convolution derivative recurrence')
p=interior_polynomial(1)
for order in range(5):
    need(at(p,-1)==0,'left C4 join order '+str(order))
    need(at(p,1)==(F(1) if order in (0,1) else F(0)),'right C4 join order '+str(order))
    p=deriv(p)

class Smooth:
    def __init__(self,nodes, values,parity, radial_derivative=False):
        self.nodes=tuple(map(F,nodes)); self.parity=tuple(parity);self.radial_derivative=radial_derivative
        self.values=tuple(tuple(map(F,row)) for row in values)
        need(len(self.nodes)==len(self.values) and self.nodes[0]==0,'shape and zero origin')
        need(all(len(v)==N for v in self.values),'dimension')
        need(all(self.nodes[i]<self.nodes[i+1] for i in range(len(self.nodes)-1)),'node order')
        need(all(v==0 for v in self.values[0]+self.values[1]),'open zero collar')
        slopes=[tuple((self.values[j+1][i]-self.values[j][i])/(self.nodes[j+1]-self.nodes[j]) for i in range(N)) for j in range(len(self.nodes)-1)]
        need(all(v==0 for v in slopes[0]),'first slope')
        self.jumps=[(self.nodes[j],tuple(slopes[j][i]-slopes[j-1][i] for i in range(N))) for j in range(1,len(slopes))]
        # Continue the final P1 segment affinely past the endpoint; no endpoint hinge.
        self.jumps=[(a,c) for a,c in self.jumps if any(c)]
    def primitive(self,t,order=0):
        t=F(t); need(0<=order<=4 and abs(t)<=self.nodes[-1],'primitive range')
        negative=t<0
        t=abs(t)
        if self.radial_derivative:
            n=order if order>0 else 0
        else: n=order+1
        terms=[(hinge(t-a,n),v) for a,v in self.jumps if t>a-D]
        out=tuple(sum((h*v[i] for h,v in terms),F(0)) for i in range(N))
        if negative: out=tuple((-1)**order*self.parity[i]*out[i] for i in range(N))
        return out

parser=argparse.ArgumentParser();parser.add_argument('--repo-root',type=Path,required=True);a=parser.parse_args()
base=a.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
def load(path,sha):
    b=path.read_bytes();need(hashlib.sha256(b).hexdigest()==sha,'source pin mismatch: '+str(path));return b
p=json.loads(load(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
runner=base/'compute/angular_ridge.py'
load(runner,'d6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f')
spec=importlib.util.spec_from_file_location('pinned_angular_ridge',runner); ridge=importlib.util.module_from_spec(spec);spec.loader.exec_module(ridge)
parity=ridge.PARITY
U=Smooth(p['nodes'],p['profiles']['u'],parity)
V=Smooth(p['nodes'],p['profiles']['v'],parity)
UR=Smooth(p['nodes'],p['profiles']['u'],tuple(-s for s in parity),True)
need(moment(0)==1 and moment(1)==0,'kernel mass/evenness')
need(hinge(-D,1)==0 and hinge(D,1)==D and hinge(0,1)==F(35,256)*D,'hinge values')
need(F(65,64)-D>1,'zero wedge')
for n in range(1,5):
    need(hinge(-D,n)==0 and hinge(D,n)==sum(F(comb(n,k),factorial(n))*D**(n-k)*(-D)**k*moment(k) for k in range(n+1)),'hinge endpoint')
for r in (F(0),F(1),F(65,64)-D):
    need(U.primitive(r)==V.primitive(r)==UR.primitive(r)==ridge.ZERO,'zero collar')
for i in (21,22):
    need(all(F(row[i])==0 for row in p['profiles']['u']+p['profiles']['v']),'TT vanishing')
need(all(u==v for u,v in zip(UR.primitive(F(3,2),1),U.primitive(F(3,2),0))),'integrated derivative compatibility')
obj=ridge.CoupledRidges(p['matrices']['A'],p['matrices']['B'],p['matrices']['D'],U,UR,V)
def free(x,r):
    a,b=r+x,r-x
    return {'value':ridge.scale(F(1,2),ridge.va(ridge.va(U.primitive(a),U.primitive(b)),ridge.vs(V.primitive(a,1),V.primitive(b,1)))),
            'dx':ridge.scale(F(1,2),ridge.va(ridge.vs(UR.primitive(a),UR.primitive(b)),ridge.va(V.primitive(a),V.primitive(b)))),
            'dr':ridge.scale(F(1,2),ridge.va(ridge.va(UR.primitive(a),UR.primitive(b)),ridge.vs(V.primitive(a),V.primitive(b))))}

# Check the mixed boundary conditions without commuting any coupled matrix with
# parity: the frozen kernel reflects after the full matrix action.
shellx=F(6,5)
bzero=free(shellx,F(0)); bone=obj.feedback(shellx,F(0),p['matrices']['M'])
need(all(bzero['value'][i]==bone['value'][i]==0 for i in range(N) if parity[i]<0),'Dirichlet boundary')
need(all(bzero['dr'][i]==0 for i in range(N) if parity[i]>0),'free Neumann boundary')
expected=ridge.mv(ridge.matrix(p['matrices']['M']),tuple(bzero['value'][i] if parity[i]>0 else F(0) for i in range(N)))
need(all(bone['dr'][i]==expected[i] for i in range(N) if parity[i]>0),'one-feedback Robin identity')
# The approximation intentionally omits affine G; it does not solve that row.

# Quantitative spot at earlier nonsmooth point; exact feedback uses physical reflection.
x,r=F(1,1000),F(1499,1000)
f=obj.feedback(x,r,p['matrices']['M'])
b0=free(x,r)
out={k:ridge.va(b0[k],f[k]) for k in ('value','dx','dr')}
need(all(out[k][i]==0 for k in out for i in (21,22)),'TT candidate zero')
print(json.dumps({'status':'PASS','kernel':'35(1-(t/delta)^2)^3/(32 delta), |t|<=delta','delta':str(D),'regularity':'profiles C4; candidate C2 subject to one-sided mixed-boundary interpretation','zero_collar':str(F(65,64)-D),'boundary':'Dirichlet exact; first Robin feedback matches M times free trace, affine source omitted', 'profiles':'smoothed P1 u and v; du is derivative of smoothed u','old_kink_point':[str(x),str(r)],'K0_at_point':{k:str(out[k][0]) for k in out},'scope':'new nominal approximation; no actual-field C2 error or pointwise Einstein residual'},sort_keys=True))
