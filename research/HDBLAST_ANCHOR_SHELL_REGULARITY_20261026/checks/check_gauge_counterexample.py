#!/usr/bin/env python3
"""Exact counterexample: H1=0 does not guarantee gauge invariance of H2
when a coordinate-defined frame is rebuilt instead of transported.
All equations concern the fixed point of a Lorentz gauge diffeomorphism.
"""
from fractions import Fraction as F
from itertools import product
import json

N=5
eta=(-1,1,1,1,1)
screen=(2,3,4)
def check(ok,msg):
    if not ok: raise RuntimeError(msg)
def g(i,j): return F(eta[i]) if i==j else F(0)
def wedge(A,B,i,j,k,l):
    return A[i][k]*B[j][l]-A[i][l]*B[j][k]-A[j][k]*B[i][l]+A[j][l]*B[i][k]
def eps(a,b,c):
    if len({a,b,c})!=3: return 0
    p=[screen.index(a),screen.index(b),screen.index(c)]
    return (-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
def W():
    A=[[F(0)]*N for _ in range(N)]
    B=[[F(0)]*N for _ in range(N)]
    for i,v in enumerate((0,1,2,4,7)): A[i][i]=F(v)
    for i,v in enumerate((0,2,5,3,11)): B[i][i]=F(v)
    B[1][2]=B[2][1]=2; B[2][3]=B[3][2]=3
    G=[[g(i,j) for j in range(N)] for i in range(N)]
    R={x:wedge(A,B,*x) for x in product(range(N),repeat=4)}
    Ric=[[sum(eta[i]*R[i,j,i,l] for i in range(N)) for l in range(N)] for j in range(N)]
    scalar=sum(eta[i]*Ric[i][i] for i in range(N))
    return {(i,j,k,l):R[i,j,k,l]-F(1,3)*wedge(G,Ric,i,j,k,l)+F(scalar,12)*(g(i,k)*g(j,l)-g(i,l)*g(j,k)) for i,j,k,l in product(range(N),repeat=4)}
def H(C):
    return {(a,b):sum(F(eps(a,c,d),4)*C[b,0,c,d]+F(eps(b,c,d),4)*C[a,0,c,d] for c in screen for d in screen) for a in screen for b in screen}
def D(C,K):
    return {inds:sum(K[p][inds[s]]*C[inds[:s]+(p,)+inds[s+1:]] for s in range(4) for p in range(N)) for inds in product(range(N),repeat=4)}
def generator(p,q):
    K=[[0]*N for _ in range(N)]; K[p][q]=eta[p]; K[q][p]=-eta[q]
    return K
C=W()
for i,j,k,l in product(range(N),repeat=4):
    check(C[i,j,k,l]==-C[j,i,k,l]==-C[i,j,l,k]==C[k,l,i,j],'curvature symmetries')
    check(C[i,j,k,l]+C[i,k,l,j]+C[i,l,j,k]==0,'Bianchi')
for j,l in product(range(N),repeat=2):
    check(sum(eta[i]*C[i,j,i,l] for i in range(N))==0,'tracefree')
check(not any(H(C).values()),'screen H1 must vanish')
# Local metric q_ab = -C_acbd x^c x^d/3. Its constant second
# derivatives determine its linearized Riemann tensor everywhere.
def q2(a,b,c,d): return -F(1,3)*(C[a,c,b,d]+C[a,d,b,c])
def linear_riemann(a,b,c,d):
    return F(1,2)*(q2(a,d,b,c)+q2(b,c,a,d)-q2(a,c,b,d)-q2(b,d,a,c))
for ijkl in product(range(N),repeat=4):
    check(linear_riemann(*ijkl)==C[ijkl],f'linearized metric realization {ijkl}')
# φ_λ(x)=exp(λK)x fixes x=0 and preserves η. For gλ=η+λq,
# gλ(0)=η and ∂gλ(0)=0, so Cλ(0)=λC, C2(0)=0 exactly.
# Rebuilding the coordinate frame after φ* gives H2'=2 H(D_K C).
# Transporting the frame by φ* instead cancels that D_K term.
K=generator(0,1)
for i,j in product(range(N),repeat=2):
    check(sum(g(i,k)*K[k][j]+g(j,k)*K[k][i] for k in range(N))==0,'K is Lorentz')
delta=H(D(C,K)); rebuilt={k:2*v for k,v in delta.items() if v}
check(bool(rebuilt),'counterexample must have nonzero gauge-dependent H2')
check(rebuilt[(3,4)]==F(-6),'fixed-point witness sign/normalization')
# Independently apply the inverse frame map to each of the four argument slots.
def frame_term(b,c,d):
    slots=[b,0,c,d]
    return -sum(K[p][slots[s]]*C[tuple(slots[:s]+[p]+slots[s+1:])] for s in range(4) for p in range(N))
def transported_derivative(a,b):
    return sum(F(eps(a,c,d),4)*frame_term(b,c,d)+F(eps(b,c,d),4)*frame_term(a,c,d) for c in screen for d in screen)
transported={k:2*(delta[k]+transported_derivative(*k)) for k in delta}
check(not any(transported.values()),'transported frame cancellation')
print(json.dumps({'result':'PASS','curvature_components_checked':N**4,'H1_zero_everywhere_for_quadratic_metric':True,'gauge_generator':'x0-x1 boost fixing origin','coordinate_rebuilt_H2_34':str(rebuilt[(3,4)]),'transported_frame_H2_34':'0','basis':'(-,+,+,+,+), screen (2,3,4)','scope':'generic exact local metric; does not meet special registered first-Weyl structure'},sort_keys=True))
