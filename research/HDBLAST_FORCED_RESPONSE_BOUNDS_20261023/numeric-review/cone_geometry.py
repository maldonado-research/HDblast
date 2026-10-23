#!/usr/bin/env python3
"""Independent rational vertical-slice geometry for an exterior backward cone."""
from fractions import Fraction as Q
import json

def rectangle_cone_moments(xlo,xhi,rlo,rhi,T=Q(1,4),R=Q(5,4)):
    """Area, first x moment, first r moment of rectangle intersect backward cone.

    This uses one-dimensional piecewise-linear slice integration, independently
    of polygon-clipping or vertex-triangulation implementations.
    """
    xlo,xhi,rlo,rhi=map(Q,(xlo,xhi,rlo,rhi));T=Q(T);R=Q(R)
    a=max(xlo,Q(0));b=min(xhi,T)
    if a>=b or rlo>=rhi:return (Q(0),Q(0),Q(0))
    left=R-T;right=R+T
    cuts=sorted({a,b,*[q for q in (right-rhi,right-rlo,rlo-left,rhi-left,T) if a<q<b]})
    answer=[Q(0),Q(0),Q(0)]
    def slab(x):
        lo=max(rlo,left+x);hi=min(rhi,right-x)
        return lo,hi
    for u,v in zip(cuts,cuts[1:]):
        mid=(u+v)/2;lm,hm=slab(mid)
        if hm<=lm:continue
        def fun(x):
            lo,hi=slab(x);z=hi-lo
            return z,x*z,(hi*hi-lo*lo)/2
        f0,fm,f1=fun(u),fun(mid),fun(v)
        for k in range(3):answer[k]+=(v-u)*(f0[k]+4*fm[k]+f1[k])/6
    return tuple(answer)

def rectangle_cone_area(*args,**kw):return rectangle_cone_moments(*args,**kw)[0]

def self_check():
    T=Q(1,4);R=Q(5,4)
    expected=(T*T,T**3/3,R*T*T)
    rows=[]
    for n in (1,2,4,8,16):
        sums=[Q(0)]*3;positive=0
        for i in range(n):
            for j in range(2*n):
                m=rectangle_cone_moments(i*T/n,(i+1)*T/n,R-T+j*T/n,R-T+(j+1)*T/n)
                if m[0]>0:positive+=1
                sums=[a+b for a,b in zip(sums,m)]
        if tuple(sums)!=expected:raise ValueError('MOMENT_NORMALIZATION')
        rows.append({'subdivision':n,'positive_cells':positive,'moments':list(map(str,sums))})
    return {'status':'PASS_INDEPENDENT_EXACT_CONE_AREA_AND_FIRST_MOMENTS','expected':list(map(str,expected)),'meshes':rows,
            'method':'Rational vertical slicing and exact degree-two integration; not polygon clipping or source evaluation.'}

if __name__=='__main__':print(json.dumps(self_check(),indent=2))
