#!/usr/bin/env python3
"""Exact 5D Weyl test for an umbilic codimension-one leaf of constant curvature."""
from fractions import Fraction as F
from itertools import product, combinations
import json
N=5
eta=(-1,1,1,1,1)
tan=(0,2,3,4)
def check(x,msg):
    if not x:raise AssertionError(msg)
def g(i,j):return F(eta[i] if i==j else 0)
def newmatrix():return [[F(0) for _ in range(N)] for _ in range(N)]
def riemann(c,v,J):
    # On-leaf R_ijkl=c(h_ik h_jl-h_il h_jk).
    # Codazzi R_nijk=h_ik v_j-h_ij v_k; sign of v is conventional.
    # R_ninj=J_ij.  Remaining entries supplied by pair symmetries.
    R={x:F(0) for x in product(range(N),repeat=4)}
    def put(i,j,k,l,z):
        for a,b,d,e,s in [(i,j,k,l,1),(j,i,k,l,-1),(i,j,l,k,-1),(j,i,l,k,1),(k,l,i,j,1),(l,k,i,j,-1),(k,l,j,i,-1),(l,k,j,i,1)]:
            old=R[a,b,d,e]
            check(old==0 or old==s*z,'inconsistent pair symmetry')
            R[a,b,d,e]=s*z
    for i,j,k,l in product(tan,repeat=4):R[i,j,k,l]=c*(g(i,k)*g(j,l)-g(i,l)*g(j,k))
    for i,j,k in product(tan,repeat=3):put(1,i,j,k,g(i,k)*v[j]-g(i,j)*v[k])
    for i,j in product(tan,repeat=2):put(1,i,1,j,J[i][j])
    for i,j,k,l in product(range(N),repeat=4):check(R[i,j,k,l]+R[i,k,l,j]+R[i,l,j,k]==0,'first Bianchi')
    return R
def weyl(R):
    Ric=[[sum(eta[i]*R[i,j,i,l] for i in range(N)) for l in range(N)] for j in range(N)]
    scalar=sum(eta[i]*Ric[i][i] for i in range(N))
    C={}
    for i,j,k,l in product(range(N),repeat=4):
        sub=g(i,k)*Ric[j][l]-g(i,l)*Ric[j][k]-g(j,k)*Ric[i][l]+g(j,l)*Ric[i][k]
        C[i,j,k,l]=R[i,j,k,l]-F(1,3)*sub+F(1,12)*scalar*(g(i,k)*g(j,l)-g(i,l)*g(j,k))
    return C
zero_v=[F(0)]*N
cases=[(F(1),zero_v,newmatrix())]
for i in tan:
    v=[F(0)]*N;v[i]=F(1);cases.append((F(0),v,newmatrix()))
for i in tan:
    for j in tan:
        if j<i:continue
        J=newmatrix();J[i][j]=J[j][i]=F(1);cases.append((F(0),zero_v,J))
for c,v,J in cases:
    C=weyl(riemann(c,v,J)); E=[[C[1,i,1,j] for j in range(N)] for i in range(N)]
    check(sum(eta[i]*E[i][i] for i in tan)==0,'E trace')
    for i,j,k in product(tan,repeat=3):check(C[1,i,j,k]==0,'mixed Weyl')
    for i,j,k,l in product(tan,repeat=4):
        wed=g(i,k)*E[j][l]-g(i,l)*E[j][k]-g(j,k)*E[i][l]+g(j,l)*E[i][k]
        check(C[i,j,k,l]==-F(1,2)*wed,'tangential Schouten form')
print(json.dumps({'cases':len(cases),'mixed_weyl_component_checks':len(cases)*4**3,'tangential_weyl_component_checks':len(cases)*4**4,'E_trace_checks':len(cases),'interpretation':'formal Gauss-Codazzi basis, not actual field evaluation'},sort_keys=True,indent=2))
