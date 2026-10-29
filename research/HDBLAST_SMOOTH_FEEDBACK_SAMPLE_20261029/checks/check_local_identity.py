#!/usr/bin/env python3
"""Independent exact polynomial audit of local unreflected ridge formulas.
A finite synthetic control, not an evaluation of the HDBLAST field.
"""
from fractions import Fraction as Q
from math import comb
import json

def need(x,m):
    if not x:raise RuntimeError(m)
def clean(p):return {k:v for k,v in p.items() if v}
def add(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,Q(0))+v
    return clean(out)
def scale(a,s):return clean({k:v*s for k,v in a.items()})
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
    out={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():out[i+k,j+l]=out.get((i+k,j+l),Q(0))+v*w
    return clean(out)
def dif(a,slot):
    out={}
    for ij,v in a.items():
        n=ij[slot]
        if n:
            kl=list(ij);kl[slot]-=1;out[tuple(kl)]=n*v
    return clean(out)
def val(a,x,r):return sum(v*x**i*r**j for (i,j),v in a.items())
x={(1,0):Q(1)};r={(0,1):Q(1)}
# A z polynomial is a vector of scalar coefficient lists, ascending powers.
A=[[Q(1,2),Q(-2,3),Q(4,5)],[Q(3,7),Q(5,11),Q(-7,13)],[Q(-2,9),Q(4,3),Q(9,8)]]
B=[[Q(2,3),Q(7,5),Q(-1,2)],[Q(-4,9),Q(8,7),Q(2,11)],[Q(6,13),Q(-5,4),Q(3,2)]]
D=[[Q(5,3),Q(-7,8),Q(2,5)],[Q(4,11),Q(3,4),Q(-8,9)],[Q(-5,7),Q(1,6),Q(7,3)]]
u=[[Q(i+1),Q(2-i,3),Q(i+2,5),Q(-3-i,7),Q(4+i,11)] for i in range(3)]
v=[[Q(2-i),Q(3+i,5),Q(-2-i,7),Q(5+i,13)] for i in range(3)]
def pa(a,b):return [sum((p[k] if k<len(p) else Q(0) for p in (a,b)),Q(0)) for k in range(max(len(a),len(b)))]
def ps(a,s):return [q*s for q in a]
def pm(M,p):
    out=[]
    for row in M:
        q=[]
        for s,z in zip(row,p):q=pa(q,ps(z,s))
        out.append(q)
    return out
def pvadd(a,b):return [pa(i,j) for i,j in zip(a,b)]
def pvsub(a,b):return pvadd(a,[ps(j,-1) for j in b])
def pd(a):return [[k*q[k] for k in range(1,len(q))] for q in a]
def pi(a):return [[Q(0)]+[q[k]/(k+1) for k in range(len(q))] for q in a]
def pshift(a):return [[Q(0)]+q for q in a]
def mm(a,b,s=1):return [[q+s*w for q,w in zip(row,rw)] for row,rw in zip(a,b)]
V=pi(v);up=pd(u)
F=[ps(q,Q(1,2)) for q in pvadd(pm(mm(A,B),pvadd(up,v)),pm(D,pvadd(u,V)))]
G=[ps(q,Q(1,2)) for q in pvadd(pm(mm(B,A,-1),pvsub(up,v)),pm(D,pvsub(u,V)))]
def ep(p,sign):
    # z=r+sign*x.
    out={}
    for n,c in enumerate(p):
        for k in range(n+1):out[k,n-k]=out.get((k,n-k),Q(0))+c*comb(n,k)*sign**k
    return clean(out)
def end(p,sign):return [ep(q,sign) for q in p]
def va(a,b):return [add(q,w) for q,w in zip(a,b)]
def vs(a,b):return [sub(q,w) for q,w in zip(a,b)]
def sc(a,s):return [scale(q,s) for q in a]
def xp(a):return [mul(x,q) for q in a]
def rp(a):return [mul(r,q) for q in a]
def mv(M,p):return [sum_polys([scale(q,s) for q,s in zip(p,row)]) for row in M]
def sum_polys(items):
    out={}
    for p in items:out=add(out,p)
    return out
def integ(p):return vs(end(pi(p),1),end(pi(p),-1))
ua,uq=end(u,1),end(u,-1)
Iu,Iv,IV=integ(u),integ(v),integ(V)
weighted_v=vs(integ(pshift(v)),rp(Iv))
weighted_V=vs(integ(pshift(V)),rp(IV))
# Integration-by-parts expression using only local u integral and endpoint data.
local=sc(va(va(va(xp(mv(A,va(ua,uq))),xp(mv(B,vs(ua,uq)))),vs(xp(mv(D,Iu)),mv(A,Iu))),va(va(mv(B,weighted_v),xp(mv(A,Iv))),mv(D,weighted_V))),Q(1,4))
# Original exact unreflected kernel primitive expression.
PF,QF=pi(F),pi(pi(F));PG,QG=pi(G),pi(pi(G))
original=sc(vs(va(vs(xp(end(PF,1)),sc(vs(end(QF,1),end(QF,-1)),Q(1,2))),sc(vs(end(QG,1),end(QG,-1)),Q(1,2))),xp(end(PG,-1))),Q(1,2))
# Local endpoint+integral formulas for first derivatives.
dx=va(sc(xp(va(end(F,1),end(G,-1))),Q(1,2)),sc(va(va(mv(B,vs(ua,uq)),mv(A,Iv)),mv(D,Iu)),Q(1,4)))
dr=va(sc(xp(vs(end(F,1),end(G,-1))),Q(1,2)),sc(va(va(sc(mv(A,vs(ua,uq)),-1),sc(mv(B,Iv),-1)),sc(mv(D,IV),-1)),Q(1,4)))
# Fully compact transcription, eliminating every global velocity primitive.
upa,upq=end(up,1),end(up,-1); va_,vq_=end(v,1),end(v,-1)
Kv=vs([mul(mul(x,x),z) for z in Iv], va(vs(integ(pshift(pshift(v))),sc(rp(integ(pshift(v))),2)),[mul(mul(r,r),z) for z in Iv]))
compact=sc(va(va(va(xp(mv(A,va(ua,uq))),xp(mv(B,vs(ua,uq)))),vs(xp(mv(D,Iu)),mv(A,Iu))),va(va(xp(mv(A,Iv)),mv(B,weighted_v)),sc(mv(D,Kv),Q(1,2)))),Q(1,4))
compact_dx=sc(va(xp(va(va(va(mv(mm(A,B),upa),mv(mm(B,A,-1),upq)),mv(D,va(ua,uq))),va(va(mv(mm(A,B),va_),mv(mm(A,B,-1),vq_)),mv(D,Iv)))),va(va(mv(B,vs(ua,uq)),mv(A,Iv)),mv(D,Iu))),Q(1,4))
compact_dr=sc(va(xp(va(va(vs(mv(mm(A,B),upa),mv(mm(B,A,-1),upq)),mv(D,vs(ua,uq))),va(mv(mm(A,B),va_),mv(mm(B,A,-1),vq_)))),va(va(sc(mv(A,vs(ua,uq)),-1),sc(mv(B,Iv),-1)),mv(D,weighted_v))),Q(1,4))
checks=0
for k in range(3):
    need(original[k]==local[k],'integration by parts local value')
    need(dif(local[k],0)==dx[k],'local time derivative')
    need(dif(local[k],1)==dr[k],'local radial derivative')
    need(local[k]==compact[k],'compact value velocity moments')
    need(dx[k]==compact_dx[k],'compact time derivative velocity moments')
    need(dr[k]==compact_dr[k],'compact radial derivative velocity moments')
    checks+=6
# Nonzero controls: changing factors or signs must change actual fixtures.
pt=Q(2,9),Q(5,3)
need(any(val(original[k],*pt)!=0 for k in range(3)),'nonzero source fixture')
need(any(val(dx[k],*pt)!=val(dr[k],*pt) for k in range(3)),'mixed derivative distinction')
wrong=va(sc(xp(va(end(F,1),end(G,-1))),Q(1,2)),sc(va(va(mv(B,vs(ua,uq)),mv(A,Iv)),mv(D,Iu)),Q(1,2)))
need(any(val(dx[k],*pt)!=val(wrong[k],*pt) for k in range(3)),'reject wrong integral factor')
wrong2=va(sc(xp(vs(end(F,1),end(G,-1))),Q(1,2)),sc(va(va(mv(A,vs(ua,uq)),mv(B,Iv)),mv(D,IV)),Q(1,4)))
need(any(val(dr[k],*pt)!=val(wrong2[k],*pt) for k in range(3)),'reject radial integral sign')
print(json.dumps({'status':'PASS','polynomial_component_identities':checks,'controls':['nonzero generic coupled fixture','wrong integral factor','wrong radial integral sign'],'scope':'synthetic exact constant-matrix unreflected ridge audit, no actual-field evaluation'},sort_keys=True))
