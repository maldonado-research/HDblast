#!/usr/bin/env python3
"""Exact algebraic witnesses only; no actual-field or new physical claims."""
from fractions import Fraction as F
from itertools import permutations,product
import json

N=5
I=range(N)

def zero():return [[F(0) for _ in I] for _ in I]
def identity():return [[F(i==j) for j in I] for i in I]
def transpose(a):return list(map(list,zip(*a)))
def matmul(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0)) for j in range(len(b[0]))] for i in range(len(a))]
def rank(rows):
    a=[list(row) for row in rows];r=0
    for c in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][c]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
    return r

def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def main():
    eta=zero();eta[0][0]=-1
    for i in range(1,N):eta[i][i]=1
    # Index 0=time, 1=radial, 2,3,4=screen. Rational rapidity witness.
    stacked=[]
    for axis in (2,3,4):
        b=identity();b[0][0]=b[axis][axis]=F(5,3);b[0][axis]=b[axis][0]=F(4,3)
        if matmul(matmul(transpose(b),eta),b)!=eta:raise ValueError('LORENTZ')
        if b[0][0]<=0:raise ValueError('TIME_ORIENTATION')
        if b[0][0]*b[axis][axis]-b[0][axis]*b[axis][0]!=1:raise ValueError('PROPER_ORIENTATION')
        if [b[i][1] for i in I]!=[F(i==1) for i in I]:raise ValueError('RADIAL_FIXED')
        stacked.extend([[b[i][j]-F(i==j) for j in I] for i in I])
    if rank(stacked)!=4:raise ValueError('COMMON_FIXED_SPACE_DIMENSION')
    if any(row[1] for row in stacked):raise ValueError('RADIAL_FIXED_COLUMN')
    # Seven independent parameters of relative adapted-frame first variation.
    for parameter in range(7):
        d=zero()
        if parameter==0:d[1][0]=d[0][1]=1
        elif parameter<=3:
            a=parameter+1;d[a][0]=d[0][a]=1
        else:
            a=parameter-2;d[a][1]=1;d[1][a]=-1
        left=matmul(transpose(d),eta);right=matmul(eta,d)
        if any(left[i][j]+right[i][j] for i,j in product(I,repeat=2)):raise ValueError('RELATIVE_FRAME_ORTHONORMALITY')
    # Synthetic algebraic Weyl example: C0=0, C1 has zero magnetic part for U0,
    # but a first observer boost has a nonzero second magnetic correction.
    # This does not use the registered actual first field or Einstein equations.
    E={1:F(1),2:F(2),3:F(3),4:F(-6)}
    C={k:F(0) for k in product(I,repeat=4)}
    for i,j in product(range(1,5),repeat=2):
        value=E[i] if i==j else F(0)
        C[0,i,0,j]=C[i,0,j,0]=value
        C[i,0,0,j]=C[0,i,j,0]=-value
    for i,j,k,l in product(range(1,5),repeat=4):
        C[i,j,k,l]=(F(i==k and j==l)*E[j]-F(i==l and j==k)*E[j]-F(j==k and i==l)*E[i]+F(j==l and i==k)*E[i])/2
    for a,b,c,d in product(I,repeat=4):
        if C[a,b,c,d]!=-C[b,a,c,d] or C[a,b,c,d]!=-C[a,b,d,c] or C[a,b,c,d]!=C[c,d,a,b]:raise ValueError('WEYL_SYMMETRIES')
        if C[a,b,c,d]+C[a,c,d,b]+C[a,d,b,c]:raise ValueError('WEYL_BIANCHI')
    for b,d in product(I,repeat=2):
        if sum((eta[a][a]*C[a,b,a,d] for a in I),F(0)):raise ValueError('WEYL_TRACE')
    eps={p:sign(tuple(x-2 for x in p)) for p in permutations((2,3,4))}
    def epsval(a,c,d):return eps.get((a,c,d),0)
    def first(a,b):return sum((F(epsval(a,c,d))*C[b,0,c,d] for c,d in product((2,3,4),repeat=2)),F(0))/2
    if any(first(a,b)+first(b,a) for a,b in product((2,3,4),repeat=2)):raise ValueError('FIRST_H_ZERO')
    # delta U=e2, delta e2=U, other frame changes zero.
    def variation(a,b):
        return sum((F(epsval(a,c,d))*(C[b,2,c,d]+(C[0,0,c,d] if b==2 else 0)+(C[b,0,0,d] if c==2 else 0)+(C[b,0,c,0] if d==2 else 0)) for c,d in product((2,3,4),repeat=2)),F(0))
    correction={(a,b):(variation(a,b)+variation(b,a))/2 for a,b in product((2,3,4),repeat=2)}
    if correction[3,4]!=F(-27,2):raise ValueError('SYNTHETIC_OBSERVER_VARIATION')
    print(json.dumps({'status':'PASS_EXACT_ANCHOR_SYMMETRY_AND_FRAME_WITNESSES','arithmetic':'rational Fraction','boost_pair':['5/3','4/3'],'common_fixed_space':'span(radial spacelike vector)','common_fixed_rank':rank(stacked),'proper_orthochronous_boosts':3,'relative_frame_parameters_checked':7,'synthetic_weyl_symmetry_trace_bianchi':True,'synthetic_H1_zero':True,'synthetic_observer_correction_34':str(correction[3,4]),'actual_field_evaluated':False,'claim':'algebraic controls for standard isotropy obstruction and observer change, not new physical result'},indent=2))

if __name__=='__main__':main()
