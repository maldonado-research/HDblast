#!/usr/bin/env python3
"""Independent rational polynomial audit of half-line reflection and shell defect.

No upstream feedback evaluator is imported. The direct reference integrates
physical source polynomials over split characteristic triangles/segments.
"""
from fractions import Fraction as Q
from math import comb
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def plus(a,b):
    out=dict(a)
    for key,v in b.items():out[key]=out.get(key,Q(0))+v
    return {key:v for key,v in out.items() if v}


def times(a,b):
    out={}
    for (i,j),u in a.items():
        for (k,l),v in b.items():out[(i+k,j+l)]=out.get((i+k,j+l),Q(0))+u*v
    return {key:v for key,v in out.items() if v}


def scaled(a,c):return {key:c*v for key,v in a.items() if c*v}
def derivative(a,axis):
    out={}
    for key,v in a.items():
        if key[axis]:
            z=list(key);z[axis]-=1;out[tuple(z)]=v*key[axis]
    return out


def affine_power(c,d,n):return {k:Q(comb(n,k))*c**(n-k)*d**k for k in range(n+1)}


def line_integral(p,c,d,lo,hi):
    """Integral p(s,c+d*s) ds, evaluated with exact polynomial expansion."""
    out=Q(0)
    for (i,j),v in p.items():
        for k,w in affine_power(c,d,j).items():
            n=i+k+1;out+=v*w*(hi**n-lo**n)/n
    return out


def strip_integral(p,lower,upper,lo,hi):
    """Integral_lo^hi ds integral_lower(s)^upper(s) dy p(s,y)."""
    anti={(i,j+1):v/(j+1) for (i,j),v in p.items()}
    return line_integral(anti,*upper,lo,hi)-line_integral(anti,*lower,lo,hi)


def unipoly(p,z):return sum((c*z**n for n,c in enumerate(p)),Q(0))
def atpoint(p,s,y):return sum((c*s**i*y**j for (i,j),c in p.items()),Q(0))
def uniadd(a,b):return [(a[i] if i<len(a) else Q(0))+(b[i] if i<len(b) else Q(0)) for i in range(max(len(a),len(b)))]
def uni_derivative(p):return [i*p[i] for i in range(1,len(p))]
def primitive(p):return [Q(0)]+[v/Q(i+1) for i,v in enumerate(p)]
def prim(p,k):
    for _ in range(k):p=primitive(p)
    return p


def ridge_poly(p,sign):
    """p(y+sign*s), represented as a physical (s,y) polynomial."""
    out={}
    for n,c in enumerate(p):
        for i in range(n+1):out[(i,n-i)]=c*comb(n,i)*sign**i
    return out


def vad(a,b):return [x+y for x,y in zip(a,b)]
def vsub(a,b):return [x-y for x,y in zip(a,b)]
def vscale(a,c):return [c*x for x in a]
def mv(a,v):return [sum((c*z for c,z in zip(row,v)),Q(0)) for row in a]
def PM(v,sign):return [z if p==sign else Q(0) for z,p in zip(v,P)]
def matpoly(a,v):
    out=[]
    for row in a:
        z={}
        for c,p in zip(row,v):z=plus(z,scaled(p,c))
        out.append(z)
    return out


P=(1,-1,1)
A=[[Q(z) for z in row] for row in ((1,2,-1),(3,-2,4),(2,1,3))]
B=[[Q(z) for z in row] for row in ((-2,3,1),(1,4,-3),(5,-1,2))]
D=[[Q(z) for z in row] for row in ((2,-1,4),(-3,2,1),(1,5,-2))]
M=[[Q(z) for z in row] for row in ((1,2,-1),(0,0,0),(2,-3,4))]
# Global polynomial extensions have the specified parity and compatible origin.
u=[[Q(0)]*4+[Q(1)],[Q(0)]*5+[Q(2)],[Q(0)]*6+[Q(-1)]]
v=[[Q(0)]*6+[Q(-2)],[Q(0)]*7+[Q(1)],[Q(0)]*4+[Q(3)]]
# Build B0 and its derivatives before the physical coupled source.
b0=[]
for U,V in zip(u,v):
    z=plus(ridge_poly(U,1),ridge_poly(U,-1))
    z=plus(z,plus(ridge_poly(primitive(V),1),scaled(ridge_poly(primitive(V),-1),-1)))
    b0.append(scaled(z,Q(1,2)))
f=[]
for fa,fb,fd in zip(matpoly(A,[derivative(z,0) for z in b0]),matpoly(B,[derivative(z,1) for z in b0]),matpoly(D,b0)):
    f.append(plus(plus(fa,fb),fd))


def evalprofiles(profiles,z,k=0):return [unipoly(prim(p,k),z) for p in profiles]

def FF(z,k=0):
    return vscale(vad(mv([[aa+bb for aa,bb in zip(ar,br)] for ar,br in zip(A,B)],vad(evalprofiles([uni_derivative(p) for p in u],z,k),evalprofiles(v,z,k))),mv(D,vad(evalprofiles(u,z,k),evalprofiles(v,z,k+1)))),Q(1,2))

def GG(z,k=0):
    return vscale(vad(mv([[bb-aa for aa,bb in zip(ar,br)] for ar,br in zip(A,B)],vsub(evalprofiles([uni_derivative(p) for p in u],z,k),evalprofiles(v,z,k))),mv(D,vsub(evalprofiles(u,z,k),evalprofiles(v,z,k+1)))),Q(1,2))


def direct(x,r):
    # Integrate full f first; reflected negative-y portion receives output parity.
    if r>=x:
        value=[strip_integral(p,(r-x,1),(r+x,-1),Q(0),x)/2 for p in f]
        pp=[line_integral(p,r+x,-1,Q(0),x) for p in f]
        mm=[line_integral(p,r-x,1,Q(0),x) for p in f]
    else:
        t=x-r
        value=[(strip_integral(p,(0,0),(r+x,-1),0,t)+par*strip_integral(p,(0,0),(t,-1),0,t)+strip_integral(p,(-t,1),(r+x,-1),t,x))/2 for p,par in zip(f,P)]
        pp=[line_integral(p,r+x,-1,0,x) for p in f]
        mm=[par*line_integral(p,t,-1,0,t)+line_integral(p,-t,1,t,x) for p,par in zip(f,P)]
    return {'value':value,'dx':vscale(vad(pp,mm),Q(1,2)),'dr':vscale(vsub(pp,mm),Q(1,2))}


def closed(x,r):
    a,q=r+x,r-x
    pp=vad(vscale(FF(a),x),vscale(vsub(GG(a,1),GG(q,1)),Q(1,2)))
    if q>=0:
        value=vscale(vsub(vad(vsub(vscale(FF(a,1),x),vscale(vsub(FF(a,2),FF(q,2)),Q(1,2))),vscale(vsub(GG(a,2),GG(q,2)),Q(1,2))),vscale(GG(q,1),x)),Q(1,2))
        mm=vad(vscale(vsub(FF(a,1),FF(q,1)),Q(1,2)),vscale(GG(q),x))
    else:
        t=-q
        physical=vsub(vad(vsub(vscale(FF(a,1),x),vscale(vad(FF(a,2),FF(t,2)),Q(1,2))),vscale(vad(GG(a,2),GG(-t,2)),Q(1,2))),vscale(GG(-t,1),r))
        reflected=vad(vsub(vscale(FF(t,1),t),FF(t,2)),vscale(vad(GG(t,2),GG(-t,2)),Q(1,2)))
        value=vscale(vad(physical,[p*z for p,z in zip(P,reflected)]),Q(1,2))
        reflected_m=vad(vscale(FF(t),t),vscale(vsub(GG(t,1),GG(-t,1)),Q(1,2)))
        mm=vad(vad(vscale(vsub(FF(a,1),FF(t,1)),Q(1,2)),vscale(GG(-t),r)),[p*z for p,z in zip(P,reflected_m)])
    return {'value':value,'dx':vscale(vad(pp,mm),Q(1,2)),'dr':vscale(vsub(pp,mm),Q(1,2))}


checks=[]
for x,r in ((Q(0),Q(1,2)),(Q(1,8),Q(3,2)),(Q(1,8),Q(1,8)),(Q(9,8),Q(1,16)),(Q(9,8),Q(0)),(Q(6,5),Q(0))):
    a,b=direct(x,r),closed(x,r)
    for key in ('value','dx','dr'):
        need(a[key]==b[key],f'direct reflected integral {x},{r} {key}')
        checks.append(f'direct_reflected_{x}_{r}_{key}')
# Shell special formulas and full Robin defect; no commutation with P.
for x in (Q(9,8),Q(6,5)):
    H=vad(vsub(vscale(FF(x,1),x),FF(x,2)),vscale(vad(GG(x,2),GG(-x,2)),Q(1,2)))
    pp=vad(vscale(FF(x),x),vscale(vsub(GG(x,1),GG(-x,1)),Q(1,2)))
    K=direct(x,Q(0))
    need(K=={'value':PM(H,1),'dx':PM(pp,1),'dr':PM(pp,-1)},'shell projector identities')
    checks.append(f'shell_projectors_{x}')
    trace0=PM(vad(evalprofiles(u,x),evalprofiles(v,x,1)),1)
    trace0x=PM(vad(evalprofiles([uni_derivative(z) for z in u],x),evalprofiles(v,x)),1)
    g=mv(M,trace0);pg=mv(M,PM(vad(evalprofiles(u,x,1),evalprofiles(v,x,2)),1))
    # Direct polynomial differentiation of the lift -integral_0^(x-r)g,
    # independently verifies both first-derivative signs.
    traceprofiles=[uniadd(U,primitive(V)) if par==1 else [Q(0)] for U,V,par in zip(u,v,P)]
    gprofiles=[]
    for row in M:
        gp=[Q(0)]
        for coeff,profile in zip(row,traceprofiles):gp=uniadd(gp,[coeff*z for z in profile])
        gprofiles.append(gp)
    liftpolys=[scaled(ridge_poly([coeff*(-1)**n for n,coeff in enumerate(primitive(gp))],-1),-1) for gp in gprofiles]
    directlift={key:[atpoint(z if axis is None else derivative(z,axis),x,Q(0)) for z in liftpolys] for key,axis in (('value',None),('dx',0),('dr',1))}
    need(directlift=={'value':vscale(pg,-1),'dx':vscale(g,-1),'dr':g},'direct lift derivative signs')
    checks.append(f'direct_lift_signs_{x}')
    b1={key:vad(K[key],directlift[key]) for key in K}
    Gshell=[x*x,Q(0),-3*x];Gshellx=[2*x,Q(0),Q(-3)]
    total=vad(trace0,b1['value']);totalr=b1['dr']
    defect=vsub(vsub(PM(totalr,1),mv(M,PM(total,1))),Gshell)
    expected=vsub(vscale(mv(M,PM(b1['value'],1)),-1),Gshell)
    need(defect==expected,'frozen affine Robin defect')
    need(PM(total,-1)==[0,0,0] and PM(b1['dx'],-1)==[0,0,0],'homogeneous Dirichlet identities')
    checks.extend((f'frozen_Robin_defect_{x}',f'Dirichlet_{x}'))
    # Boundary derivative defect is differentiated algebra, independent of
    # values: B1_rx=M B0_x, hence derivative residual=-M B1_x-G_x.
    totalrx=mv(M,trace0x);totalx=vad(trace0x,b1['dx'])
    defectx=vsub(vsub(PM(totalrx,1),mv(M,PM(totalx,1))),Gshellx)
    need(defectx==vsub(vscale(mv(M,PM(b1['dx'],1)),-1),Gshellx),'differentiated Robin defect')
    checks.append(f'differentiated_frozen_Robin_defect_{x}')
    # M above is the frozen M0. An actual coefficient can differ; evaluate
    # the target operator directly before comparing its correction identity.
    delta0=[[Q(1,7),Q(0),Q(-2,5)],[Q(0)]*3,[Q(3,11),Q(0),Q(-1,13)]]
    delta=[[x*c for c in row] for row in delta0]
    actualM=[[a+b for a,b in zip(row,change)] for row,change in zip(M,delta)]
    actual_defect=vsub(vsub(PM(totalr,1),mv(actualM,PM(total,1))),Gshell)
    corrected=vsub(expected,mv(delta,PM(total,1)))
    need(actual_defect==corrected,'actual Robin coefficient correction')
    need(actual_defect!=expected,'omitted actual coefficient correction must fail')
    checks.append(f'actual_coefficient_correction_{x}')
    # Here M_actual(x)=M0+x*delta0, so its derivative is delta0.
    actual_defectx=vsub(vsub(vsub(PM(totalrx,1),mv(actualM,PM(totalx,1))),mv(delta0,PM(total,1))),Gshellx)
    correctedx=vsub(vsub(vsub(vscale(mv(M,PM(b1['dx'],1)),-1),mv(delta,PM(totalx,1))),mv(delta0,PM(total,1))),Gshellx)
    need(actual_defectx==correctedx,'actual variable-coefficient derivative correction')
    checks.append(f'actual_variable_coefficient_derivative_{x}')
# A plausible wrong implementation reflects inputs and applies matrices later.
# That uses the globally extended free-wave source on negative y, not Pf(s,-y).
x,r=Q(9,8),Q(1,16)
naive=[strip_integral(p,(r-x,1),(r+x,-1),0,x)/2 for p in f]
need(naive!=direct(x,r)['value'],'wrong input-reflection formula must fail')
need(any(A[i][j] and P[i]!=P[j] for i in range(3) for j in range(3)),'noncommuting coupling control')
checks.extend(('wrong_input_reflection_rejected','noncommuting_coupling_present'))
print(json.dumps({'status':'PASS','named_check_groups':len(checks),'checks':checks,'method':'direct rational physical-source polynomial triangle and characteristic-line integrations; no upstream feedback import','scope':'finite independent algebra controls for exact reflected kernel, shell projectors, Robin lift and candidate defect; no numeric HDBLAST source bounds or actual solution'},sort_keys=True))
