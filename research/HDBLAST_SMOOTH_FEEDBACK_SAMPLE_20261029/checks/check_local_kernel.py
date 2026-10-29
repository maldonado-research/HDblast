#!/usr/bin/env python3
"""Exact polynomial controls for the unreflected local feedback kernel.
Synthetic scalar tests only; no HDBLAST field/residual computation.
"""
from fractions import Fraction as Q
from math import comb
import json

def check(ok,why):
    if not ok: raise RuntimeError(why)
def add(a,b):
    n=max(len(a),len(b)); return [(a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0)) for i in range(n)]
def scale(a,c): return [c*x for x in a]
def diff(a): return [i*a[i] for i in range(1,len(a))] or [Q(0)]
def integ(a): return [Q(0)]+[x/Q(i+1) for i,x in enumerate(a)]
def val(a,z): return sum((x*z**i for i,x in enumerate(a)),Q(0))
def area(a,q,b):
    p=integ(a);return val(p,b)-val(p,q)
def times_z_shift(a,c): return add([Q(0)]+list(a),scale(a,-c))
def triangle(F,G,x,r):
    q,a=r-x,r+x
    total=Q(0)
    # Direct integration of each monomial on s in [0,x], w in [q+s,a-s].
    for n,f in enumerate(F):
        lower=sum((Q(comb(n+1,k))*q**(n+1-k)*2**k*x**(k+1)/Q(k+1) for k in range(n+2)),Q(0))
        total+=f*(x*a**(n+1)-lower)/Q(2*(n+1))
    for n,g in enumerate(G):
        upper=sum((Q(comb(n+1,k))*a**(n+1-k)*(-2)**k*x**(k+1)/Q(k+1) for k in range(n+2)),Q(0))
        total+=g*(upper-x*q**(n+1))/Q(2*(n+1))
    return total

def deriv_exact(f,t,h,n=8):
    weights=[-sum((Q(1,k) for k in range(1,n+1)),Q(0))]+[Q((-1)**(k+1)*comb(n,k),k) for k in range(1,n+1)]
    return sum((weights[k]*f(t+k*h) for k in range(n+1)),Q(0))/h

checks=0
for fixture in range(1,4):
    u=list(map(Q,[1+fixture,-2,3,-fixture,2]));v=list(map(Q,[2,-fixture,1,3]));V=integ(v)
    A,B,D=Q(2*fixture-1,3),Q(-3-fixture,5),Q(7+fixture,11)
    F=scale(add(scale(add(diff(u),v),A+B),scale(add(u,V),D)),Q(1,2))
    G=scale(add(scale(add(diff(u),scale(v,-1)),B-A),scale(add(u,scale(V,-1)),D)),Q(1,2))
    x,r=Q(1,10),Q(3,2);q,a=r-x,r+x
    direct=triangle(F,G,x,r)
    local=Q(1,4)*(area(times_z_shift(F,q),q,a)-area(times_z_shift(G,a),q,a))
    Iu=area(u,q,a)
    parts=Q(1,4)*(A*(x*(val(u,a)+val(u,q))-Iu)+x*B*(val(u,a)-val(u,q))
        +B*area(times_z_shift(v,r),q,a)+x*A*area(v,q,a)+x*D*Iu+D*area(times_z_shift(V,r),q,a))
    check(direct==local==parts,'triangle/local/integration-by-parts disagreement');checks+=2
    dx=x*(val(F,a)+val(G,q))/2+area(add(F,G),q,a)/4
    dr=x*(val(F,a)-val(G,q))/2+area(add(G,scale(F,-1)),q,a)/4
    h=Q(1,1000)
    check(deriv_exact(lambda xx:triangle(F,G,xx,r),x,h)==dx,'x derivative');checks+=1
    check(deriv_exact(lambda rr:triangle(F,G,x,rr),r,h)==dr,'r derivative');checks+=1
    check(add(F,G)==add(add(scale(diff(u),B),scale(v,A)),scale(u,D)),'F+G cancellation');checks+=1
    check(direct != local*2,'wrong kernel factor control');checks+=1
    check(parts != parts-x*B*(val(u,a)-val(u,q))/4,'missing B endpoint control');checks+=1
print(json.dumps({'status':'PASS','fixtures':3,'checks':checks,'scope':'synthetic rational polynomials; local unreflected feedback identity, not actual field or residual','negative_controls':['wrong kernel factor','missing B endpoint term']},sort_keys=True))
