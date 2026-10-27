#!/usr/bin/env python3
"""Exact algebraic-basis test of screen-H under arbitrary first Lorentz tilt.
No use is made of the frozen HDBLAST curvature implementations.
"""
from fractions import Fraction as F
from itertools import product, combinations
import json
N=5
eta=[-1,1,1,1,1]
screen=(2,3,4)
tangent=(0,2,3,4)
def g(i,j): return eta[i] if i==j else 0
def eps(i,j,k):
    if len({i,j,k})<3:return 0
    p=[screen.index(i),screen.index(j),screen.index(k)]
    return (-1)**sum(p[a]>p[b] for a in range(3) for b in range(a+1,3))
def wedge(A,B,i,j,k,l):return A[i][k]*B[j][l]-A[i][l]*B[j][k]-A[j][k]*B[i][l]+A[j][l]*B[i][k]
def build_electric(E):
    # C_{n i n j}=E_ij; C tangent = -1/2(h wedge E).
    C={}
    for i,j,k,l in product(range(N),repeat=4):
        z=F(0)
        if all(x != 1 for x in (i,j,k,l)):
            z=-F(1,2)*(g(i,k)*E[j][l]-g(i,l)*E[j][k]-g(j,k)*E[i][l]+g(j,l)*E[i][k])
        elif sum(x==1 for x in (i,j,k,l))==2:
            if i==1 and k==1:z=E[j][l]
            if i==1 and l==1:z=-E[j][k]
            if j==1 and k==1:z=-E[i][l]
            if j==1 and l==1:z=E[i][k]
        C[i,j,k,l]=z
    return C
def Kgen(p,q):
    # K_{p q}=1, K_{q p}=-1, K^i_j=eta^ii K_ij
    K=[[0]*N for _ in range(N)]; K[p][q]=eta[p];K[q][p]=-eta[q]
    return K
def dcomp(C,K,i,j,k,l):
    slots=(i,j,k,l); z=F(0)
    for t in range(4):
        for p in range(N):
            v=list(slots);v[t]=p;z+=K[p][slots[t]]*C[tuple(v)]
    return z
def screenH(C,K=None):
    def component(i,j,k,l):return C[i,j,k,l] if K is None else dcomp(C,K,i,j,k,l)
    return {(a,b):sum(F(eps(a,c,d),4)*component(b,0,c,d)+F(eps(b,c,d),4)*component(a,0,c,d) for c in screen for d in screen) for a in screen for b in screen}
def check(condition, message):
    if not condition: raise AssertionError(message)
def tensor_checks(C):
    for i,j,k,l in product(range(N),repeat=4):
        check(C[i,j,k,l]==-C[j,i,k,l]==-C[i,j,l,k]==C[k,l,i,j], "Weyl pair symmetry")
        check(C[i,j,k,l]+C[i,k,l,j]+C[i,l,j,k]==0, "Weyl cyclic symmetry")
    for j,l in product(range(N),repeat=2): check(sum(eta[i]*C[i,j,i,l] for i in range(N))==0, "Weyl trace")

def generic_weyl(A,B):
    G=[[g(i,j) for j in range(N)] for i in range(N)]
    R={(i,j,k,l):F(wedge(A,B,i,j,k,l)) for i,j,k,l in product(range(N),repeat=4)}
    Ric=[[sum(eta[i]*R[i,j,i,l] for i in range(N)) for l in range(N)] for j in range(N)]
    s=sum(eta[i]*Ric[i][i] for i in range(N))
    return {(i,j,k,l):R[i,j,k,l]-F(1,3)*wedge(G,Ric,i,j,k,l)+F(s,12)*(g(i,k)*g(j,l)-g(i,l)*g(j,k)) for i,j,k,l in product(range(N),repeat=4)}

# A basis of nine symmetric traceless tensors on 4D Lorentzian tangent plane.
basis=[]
for p,q in combinations(tangent,2):
    E=[[F(0)]*N for _ in range(N)]; E[p][q]=E[q][p]=F(1);basis.append(E)
for p in (0,2,3):
    E=[[F(0)]*N for _ in range(N)];E[p][p]=F(1);E[4][4]=-eta[p];basis.append(E)
check(len(basis)==9, "basis dimension")
checked=0
for E in basis:
    C=build_electric(E);tensor_checks(C)
    check(not any(screenH(C).values()), "first screen magnetic tensor")
    for p,q in combinations(range(N),2):
        K=Kgen(p,q)
        check(not any(screenH(C,K).values()), "frame derivative")
        checked+=9
# A generic first Weyl may have H1=0 but respond to an observer change.
A=[[F(0)]*N for _ in range(N)];B=[[F(0)]*N for _ in range(N)]
for i,v in enumerate((0,1,2,4,7)): A[i][i]=F(v)
for i,v in enumerate((0,2,5,3,11)):B[i][i]=F(v)
B[1][2]=B[2][1]=F(2)
B[2][3]=B[3][2]=F(3)
C=generic_weyl(A,B);tensor_checks(C)
check(not any(screenH(C).values()), "first screen magnetic tensor")
found=[]
for p,q in combinations(range(N),2):
    x=screenH(C,Kgen(p,q)); nz={str(k):str(v) for k,v in x.items() if v}
    if nz:found.append({'generator':[p,q], 'nonzero':nz})
check(bool(found), "negative control did not respond")
print(json.dumps({'tangent_symmetric_tracefree_basis_count':len(basis),'lorentz_generators':10,'component_zero_checks':checked,'generic_weyl_H1_zero':True,'generic_counterexample':found},sort_keys=True,indent=2))
