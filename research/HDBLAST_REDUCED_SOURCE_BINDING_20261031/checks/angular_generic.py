#!/usr/bin/env python3
"""Coefficient-generic exact polynomial tensor calculus on S^3.

Polynomials are reduced modulo n_3^2=1-n_0^2-n_1^2-n_2^2.
Intrinsic tensors retain four ambient tangent slots.  D projects every slot.
Derived from the frozen parent angular_exact.py except clean() preserves exact
coefficient objects rather than coercing them to Fraction. This permits fully
independent formal radial/time jets.
"""
from fractions import Fraction as Q
from itertools import product, permutations
from functools import lru_cache
from math import factorial

I=range(4)
ZERO=(0,0,0,0)
ONE={ZERO:Q(1)}

def clean(p): return {e:v for e,v in p.items() if v}
def add(*ps):
    out={}
    for p in ps:
        for e,v in p.items(): out[e]=out.get(e,Q(0))+v
    return clean(out)
def scale(c,p): return clean({e:c*v for e,v in p.items()})
@lru_cache(None)
def monomial_reduce(e):
    if e[3]<2: return ((e,Q(1)),)
    out={}; f=list(e); f[3]-=2
    for j,c in [(-1,1),(0,-1),(1,-1),(2,-1)]:
        g=f.copy()
        if j>=0:g[j]+=2
        for h,v in monomial_reduce(tuple(g)):out[h]=out.get(h,Q(0))+c*v
    return tuple(clean(out).items())
def reduce(p):
    out={}
    for e,v in p.items():
        for f,w in monomial_reduce(e):out[f]=out.get(f,Q(0))+v*w
    return clean(out)
def mul(p,q):
    out={}
    for e,v in p.items():
        for f,w in q.items():
            h=tuple(a+b for a,b in zip(e,f));out[h]=out.get(h,Q(0))+v*w
    return reduce(out)
def diff(p,i):
    out={}
    for e,v in p.items():
        if e[i]:
            f=list(e);f[i]-=1;out[tuple(f)]=v*e[i]
    return clean(out)
X=[]
for i in I:
    e=[0]*4;e[i]=1;X.append({tuple(e):Q(1)})
P={(i,j):add(ONE if i==j else {},scale(-1,mul(X[i],X[j]))) for i,j in product(I,repeat=2)}
def zero_tensor(rank):return {k:{} for k in product(I,repeat=rank)}
def tadd(*ts):return {k:add(*(t[k] for t in ts)) for k in ts[0]}
def tscale(c,t):return {k:scale(c,p) for k,p in t.items()}
def tmul(p,t):return {k:mul(p,q) for k,q in t.items()}
def outer(t,u):return {a+b:mul(p,q) for a,p in t.items() for b,q in u.items()}
def covd(t):
    """D_i T_J=P_i^a partial_a T_J+sum_s n_js T_(J[js->i])."""
    out={}
    for j,p in t.items():
        for i in I:
            terms=[mul(P[i,a],diff(p,a)) for a in I]
            for s,b in enumerate(j):
                jj=list(j);jj[s]=i;terms.append(mul(X[b],t[tuple(jj)]))
            out[(i,)+j]=add(*terms)
    return out
def trace(t,first=0,second=1):
    rank=len(next(iter(t)))
    out={}
    for k in product(I,repeat=rank-2):
        terms=[]
        for i in I:
            j=list(k);j.insert(first,i);j.insert(second,i);terms.append(t[tuple(j)])
        out[k]=add(*terms)
    return out
def div(t):return trace(covd(t))
def lap(t):return trace(covd(covd(t)))
def grad(p):return covd({():p})
def hess(p):return covd(grad(p))
def tf(t):return tadd(t,tmul(scale(-Q(1,3),trace(t)[()]),P))
def moment(e):
    if any(k%2 for k in e):return Q(0)
    q=Q(1)
    for k in e:
        for j in range(k//2):q*=2*j+1
    for j in range(sum(e)//2):q/=4+2*j
    return q
def integrate(p):return sum((v*moment(e) for e,v in p.items()),Q(0))
def inner(t,u):return sum((integrate(mul(t[k],u[k])) for k in t),Q(0))
def dot(t,u):return add(*(mul(t[k],u[k]) for k in t))
def equal(t,u):return all(not reduce(add(t[k],scale(-1,u[k]))) for k in t)
def expect(t,u,label):
    if not equal(t,u):raise ValueError(label)
def scalar(p):return {():p}
def project(t,basis,label='PROJECTION'):
    coeff=[inner(t,b)/inner(b,b) for b in basis]
    recon=tadd(*(tscale(c,b) for c,b in zip(coeff,basis)))
    expect(t,recon,label+'_ORTHOGONAL_RECONSTRUCTION')
    return coeff
def sign(p):return (-1)**sum(p[i]>p[j] for i in I for j in range(i+1,4))
EPS={p:sign(p) for p in permutations(I)}
def make_tt():
    d={(0,1):Q(-2,3),(2,3):Q(-2,3),(0,2):Q(1,3),(0,3):Q(1,3),(1,2):Q(1,3),(1,3):Q(1,3)}
    C={k:Q(0) for k in product(I,repeat=4)}
    for (i,j),q in d.items():
        C[i,j,i,j]=C[j,i,j,i]=q;C[i,j,j,i]=C[j,i,i,j]=-q
    star={k:sum((Q(EPS.get((k[0],k[1],a,b),0))*C[a,b,k[2],k[3]]/2 for a,b in product(I,repeat=2)),Q(0)) for k in C}
    tensors=[]
    for s in (1,-1):
        W={k:(C[k]+s*star[k])/2 for k in C}
        tensors.append({(i,j):add(*(scale(W[i,a,j,b],mul(X[a],X[b])) for a,b in product(I,repeat=2))) for i,j in product(I,repeat=2)})
    return tensors
S=add(mul(X[0],X[0]),mul(X[1],X[1]))
Y=add(mul(X[0],X[0]),scale(-1,mul(X[1],X[1])))
Z=add(S,scale(-Q(1,2),ONE))
F=add(mul(Y,Y),scale(-Q(1,2),S),scale(Q(1,12),ONE))
SCALARS=[scalar(p) for p in (ONE,Z,F)]
GRADS=[grad(p) for p in (Z,F)]
HESS=[tf(hess(p)) for p in (Z,F)]
TT=make_tt()
TENSORS=[tmul(p,P) for p in (ONE,Z,F)]+HESS+TT
SCALAR_NAMES=['1','Z','F'];VECTOR_NAMES=['DZ','DF'];TENSOR_NAMES=['kappa','Zkappa','Fkappa','HZ','HF','Tplus','Tminus']

def primitive_matrices(mutation=None):
    result={}
    operations=[('grad',SCALARS,GRADS,lambda t:covd(t)),
                ('hess',SCALARS,TENSORS,lambda t:covd(covd(t))),
                ('div_vector',GRADS,SCALARS,div),
                ('symgrad',GRADS,TENSORS,lambda t:tscale(2,covd(t))),
                ('div_tensor',TENSORS,GRADS,div),
                ('trace',TENSORS,SCALARS,trace),
                ('lap_scalar',SCALARS,SCALARS,lap),
                ('lap_vector',GRADS,GRADS,lap),
                ('lap_tensor',TENSORS,TENSORS,lap)]
    for name,src,dst,fn in operations:
        columns=[project(fn(b),dst,name) for b in src]
        result[name]=[list(row) for row in zip(*columns)]
    if mutation=='rough_tensor_shift':result['lap_tensor'][3][3]-=1
    return result
