#!/usr/bin/env python3
"""Exact algebraic test of all 10 Lorentz frame directions; no field equations.

Frame order U,e1,e2,e3,n; signature -++++, C_abcd convention arbitrary.
Reproduces H2 = H(C2) + 2 delta_frame H(C1). Fractions only.
"""
from fractions import Fraction as F
from itertools import product, combinations
import json

N=5
eta=[-1,1,1,1,1]
def metric(a,b): return eta[a] if a==b else 0
def kn(A,B,a,b,c,d):
    return A[a][c]*B[b][d]+B[a][c]*A[b][d]-A[a][d]*B[b][c]-B[a][d]*A[b][c]
def weyl_from_kn(A,B):
    R={(a,b,c,d): F(kn(A,B,a,b,c,d)) for a,b,c,d in product(range(N), repeat=4)}
    Ric={(b,d): sum(eta[a]*R[a,b,a,d] for a in range(N)) for b,d in product(range(N),repeat=2)}
    scalar=sum(eta[a]*Ric[a,a] for a in range(N))
    def C(a,b,c,d):
        return R[a,b,c,d] - F(1,3)*(metric(a,c)*Ric[b,d]-metric(a,d)*Ric[b,c]-metric(b,c)*Ric[a,d]+metric(b,d)*Ric[a,c]) + F(1,12)*scalar*(metric(a,c)*metric(b,d)-metric(a,d)*metric(b,c))
    return C

def eps(a,b,c):
    if len({a,b,c})<3: return 0
    return 1 if (a,b,c) in [(1,2,3),(2,3,1),(3,1,2)] else -1

def H(C,a,b):
    return sum(F(eps(a,c,d),4)*C(b,0,c,d)+F(eps(b,c,d),4)*C(a,0,c,d) for c,d in product(range(1,4),repeat=2))

def delta_component(C,indices,G):
    total=F(0)
    for slot in range(4):
        i=indices[slot]
        for j in range(N):
            trial=list(indices); trial[slot]=j
            total+=G[i][j]*C(*trial)
    return total

def delta_H(C,a,b,G):
    return sum(F(eps(a,c,d),4)*delta_component(C,(b,0,c,d),G)+F(eps(b,c,d),4)*delta_component(C,(a,0,c,d),G) for c,d in product(range(1,4),repeat=2))

def empty(): return [[0]*N for _ in range(N)]
def direction(i,j):
    # E_i -> G_i^j E_j and E_j -> G_j^i E_i, preserving eta.
    G=empty(); G[i][j]=1; G[j][i]=-eta[j]*eta[i]
    return G

def require(condition, msg):
    if not condition: raise RuntimeError(msg)

def count_nonzero(C):
    return sum(C(*idx)!=0 for idx in product(range(N),repeat=4))

def weyl_axioms(C):
    for a,b,c,d in product(range(N),repeat=4):
        require(C(a,b,c,d)==-C(b,a,c,d)==-C(a,b,d,c), 'antisymmetry')
        require(C(a,b,c,d)==C(c,d,a,b), 'pair symmetry')
        require(C(a,b,c,d)+C(a,c,d,b)+C(a,d,b,c)==0, 'Bianchi')
    for b,d in product(range(N),repeat=2):
        require(sum(eta[a]*C(a,b,a,d) for a in range(N))==0, 'trace')

# Build 9 independent trace-free symmetric S over the 4D tangent spacetime.
S_basis=[]
for i in range(1,4):
    S=empty(); S[0][0]=1; S[i][i]=1; S_basis.append(S)
for i,j in combinations(range(4),2):
    S=empty(); S[i][j]=S[j][i]=1; S_basis.append(S)

def structured(S):
    def C(a,b,c,d):
        if 4 not in (a,b,c,d): return F(kn([[metric(i,j) for j in range(N)] for i in range(N)],S,a,b,c,d))
        # Weyl completion enforces C_{i n j n}=-2 S_ij.
        T=[a,b,c,d]
        if T.count(4)!=2: return F(0)
        if a==4 and c==4: return F(-2*S[b][d])
        if a==4 and d==4: return F(2*S[b][c])
        if b==4 and c==4: return F(2*S[a][d])
        if b==4 and d==4: return F(-2*S[a][c])
        return F(0)
    return C

# A generic C2 independently constructed via the 5D Weyl projection.
A=empty(); B=empty()
for i,j,v in [(0,1,1),(1,2,2),(3,4,1),(0,4,-1)]: A[i][j]=A[j][i]=v
for i,j,v in [(0,2,1),(1,3,1),(2,4,2),(0,3,-1)]: B[i][j]=B[j][i]=v
C2=weyl_from_kn(A,B); weyl_axioms(C2)
require(any(H(C2,a,b)!=0 for a,b in product(range(1,4),repeat=2)), 'C2 selected H must be nonzero')

dirs=[direction(i,j) for i,j in combinations(range(5),2)]
checks=0
for S in S_basis:
    C1=structured(S); weyl_axioms(C1)
    require(count_nonzero(C1)>0,'must exercise nonzero C1')
    for a,b in product(range(1,4),repeat=2):
        require(H(C1,a,b)==0,'structured background H1')
        for G in dirs:
            delta=delta_H(C1,a,b,G)
            require(delta==0,'structured frame variation')
            direct_H2=H(C2,a,b)+2*delta
            require(direct_H2==H(C2,a,b),'C2 must survive frame test')
            checks+=1

# Counterexample: generic algebraic Weyl with one-radial, three-tangent
# components, yet all-tangent C_abcd=0 and hence H1=0 in the initial frame.
A=empty(); B=empty(); A[0][4]=A[4][0]=1; B[1][2]=B[2][1]=1
Cbad=weyl_from_kn(A,B); weyl_axioms(Cbad)
require(count_nonzero(Cbad)>0,'generic C1 must be nonzero')
require(all(H(Cbad,a,b)==0 for a,b in product(range(1,4),repeat=2)), 'generic control H1 must vanish')
witness=[]
for k,G in enumerate(dirs):
    for a,b in product(range(1,4),repeat=2):
        v=2*delta_H(Cbad,a,b,G)
        if v: witness.append((list(combinations(range(5),2))[k],a,b,str(v)))
require(witness,'generic H1=0 must not imply frame independence')

# Independently expand the FULL C(lambda) contraction, interpolating its exact
# lambda^2 coefficient from 7 integer samples. Frame orthonormality to first
# order suffices because C0=0; higher frame changes multiply zero background C.
def mul_poly(P,Q):
    out=[F(0)]*(len(P)+len(Q)-1)
    for i,p in enumerate(P):
        for j,q in enumerate(Q): out[i+j]+=p*q
    return out

def coeff2_interpolated(values):
    ans=F(0)
    for k in range(7):
        p=[F(1)]; den=F(1)
        for j in range(7):
            if j!=k:
                p=mul_poly(p,[F(-j),F(1)]); den*=k-j
        ans+=values[k]*p[2]/den
    return ans

def frame_eval(C1,C2,G,t,a,b):
    def frame(i):
        return [(j, F(i==j)+t*G[i][j]) for j in range(N) if F(i==j)+t*G[i][j]]
    def comp(i,j,k,l):
        out=F(0)
        for (aa,pa),(bb,pb),(cc,pc),(dd,pd) in product(frame(i),frame(j),frame(k),frame(l)):
            out+=pa*pb*pc*pd*(t*C1(aa,bb,cc,dd)+t*t*C2(aa,bb,cc,dd)/2)
        return out
    return sum(F(eps(a,c,d),4)*comp(b,0,c,d)+F(eps(b,c,d),4)*comp(a,0,c,d) for c,d in product(range(1,4),repeat=2))

Gfull=empty()
for k,G in enumerate(dirs):
    for i,j in product(range(N),repeat=2): Gfull[i][j]+=(k+1)*G[i][j]
Sfull=empty()
for k,S in enumerate(S_basis):
    for i,j in product(range(N),repeat=2): Sfull[i][j]+=(k+1)*S[i][j]
interpolation_checks=0
for C1,G in [(structured(Sfull),Gfull),(Cbad,direction(1,4))]:
    for a,b in product(range(1,4),repeat=2):
        vals=[frame_eval(C1,C2,G,F(t),a,b) for t in range(7)]
        numerical_second=2*coeff2_interpolated(vals)
        require(numerical_second==H(C2,a,b)+2*delta_H(C1,a,b,G), 'full polynomial interpolation mismatch')
        interpolation_checks+=1
print(json.dumps({'status':'PASS','signature':'-++++','structured_STF_basis':len(S_basis),'independent_Lorentz_generators':len(dirs),'frame_component_checks':checks,'full_contraction_interpolation_checks':interpolation_checks,'nonzero_generic_C2_components':count_nonzero(C2),'nonzero_selected_C2_H_components':sum(H(C2,a,b)!=0 for a,b in product(range(1,4),repeat=2)),'generic_H1_zero_frame_dependence_witnesses':len(witness),'first_witness':witness[0],'scope':'finite-dimensional algebraic test; does not derive metric-to-Weyl structure, gauge invariance, or nonlinear solution'},sort_keys=True))
