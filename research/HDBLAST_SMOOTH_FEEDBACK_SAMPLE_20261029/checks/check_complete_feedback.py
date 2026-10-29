#!/usr/bin/env python3
"""Interval evaluation of one complete tighter nominal B0+B1 modal sample.
No bulk/shell affine forcing is inserted; no C2-error or Einstein residual.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib.util, json, math

def need(c,m):
    if not c: raise RuntimeError(m)
def pin(path,h):
    raw=path.read_bytes(); need(hashlib.sha256(raw).hexdigest()==h,'INPUT_HASH:'+str(path)); return raw

def box(x):return (F(x),F(x))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def scale(a,c):
    q=(a[0]*c,a[1]*c);return (min(q),max(q))
def widen(a,e):return (a[0]-e,a[1]+e)
def vsum(xs):
    z=box(0)
    for a in xs:z=add(z,a)
    return z
def vadd(a,b):return [add(x,y) for x,y in zip(a,b)]
def vsub(a,b):return [sub(x,y) for x,y in zip(a,b)]
def vscale(a,c):return [scale(x,c) for x in a]
def dot(row,vec):return vsum(scale(v,c) for c,v in zip(row,vec))
def pointbox(a):return (a.lo,a.hi)
def encode(a):return list(map(str,a))

ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--panels',type=int,default=32);args=ap.parse_args()
need(1<=args.panels<=4096,'panel count')
root=args.repo_root.resolve();base=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
init=base/'initial/full23_initial.py'
pin(init,'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1')
report=json.loads(pin(base/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
prof=json.loads(pin(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
response=json.loads(pin(base/'compute/response-normal/RESPONSE_REPORT.json','ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd'))
need(prof['labels']==report['labels']==response['labels'],'label agreement')
need(report['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER', 'initial scope status')
need(response['status']=='PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'prior response status')
need(report['derivative_convention']['du']=='partial_r U' and report['derivative_convention']['v']=='partial_x U; x=t/ell', 'initial derivative conventions')
labels=prof['labels'];weights=list(map(F,prof['gram']));idx=12
need(weights==list(map(F,report['gram']))==list(map(F,response['gram'])) and min(weights)>0, 'positive common Gram norm')
need(labels[idx]=='w:0' and weights[idx]==1,'modal norm')
need(report['derivative_convention']['second_derivatives'][0]=='partial_r^2 U','M2 derivative convention')
M2=list(map(F,report['nominal_second_derivative_component_caps'][0]));need(len(M2)==23 and min(M2)>=0,'component second derivative cap')
M2H=F(report['nominal_second_derivative_H_caps'][0]);need(sum(w*m*m for w,m in zip(weights,M2))<=M2H*M2H,'M2 norm relation')
spec=importlib.util.spec_from_file_location('feedback_nominal_initial',init);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
ctx=mod.Context(root) # authenticates every upstream root before scientific imports
x,r=F(1,1000),F(1499,1000);a,b=r+x,r-x;du=F(1,2**30);dv=F(1,2048)
need(0<x<r and F(66,64)<b-du<a+du<F(7,4),'unreflected branch and analytic initial piece')
need(r>max(F(0),x-F(5,4)) and r<F(7,4)-x,'D+ interior')
au=[m*du*du/18 for m in M2];bu=[m*35*du/128 for m in M2]
def endpoint(z):
    values,_=ctx.profiles(z)
    return ([widen(pointbox(values[0][i]),au[i]) for i in range(23)], [widen(pointbox(values[1][i]),bu[i]) for i in range(23)])
ua,upa=endpoint(a);ub,upb=endpoint(b)
# On each complete local integration panel, the frozen interval evaluator
# encloses u0 at every point. Positive-length interval sums enclose its integral.
# Convolution error is added separately using the global C1,1 bound.
Iu=[box(0) for _ in range(23)];N=args.panels
for k in range(N):
    left=b+(a-b)*k/N;right=b+(a-b)*(k+1)/N
    values,_=ctx.profiles(ctx.I(left,right))
    Iu=vadd(Iu,[scale(pointbox(values[0][i]),right-left) for i in range(23)])
Iu=[widen(Iu[i],(a-b)*au[i]) for i in range(23)]
# Exact compact-kernel convolution of the continuous P1 velocity and its
# first three normalized primitives; the zero collar and affine endpoint
# continuation justify a positive-half-line hinge representation here.
nodes=list(map(F,prof['nodes']));vv=[[F(v) for v in row] for row in prof['profiles']['v']]
need(all(v==0 for row in vv[:2] for v in row),'velocity zero collar')
slopes=[[(vv[j+1][i]-vv[j][i])/(nodes[j+1]-nodes[j]) for i in range(23)] for j in range(len(nodes)-1)]
need(not any(slopes[0]),'initial affine slope zero')
jumps=[(nodes[j],[slopes[j][i]-slopes[j-1][i] for i in range(23)]) for j in range(1,len(slopes))]
def moment(k):
    if k%2:return F(0)
    return F(35,32)*sum(F(2*c,2*j+k+1) for j,c in enumerate((1,-3,3,-1)))
def hinge(z,n):
    if z<=-dv:return F(0)
    if z>=dv:return sum(F(math.comb(n,k),math.factorial(n))*z**(n-k)*(-dv)**k*moment(k) for k in range(n+1))
    s=z/dv;out=F(0)
    for j,c in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            m=2*j+k;out+=c*math.comb(n,k)*s**(n-k)*(-1)**k*(s**(m+1)-(-1)**(m+1))/F(m+1)
    return F(35,32)*dv**n*out/F(math.factorial(n))
def velocity(z,order=0):
    out=[F(0)]*23
    for knot,jump in jumps:
        if z>knot-dv:
            h=hinge(z-knot,order+1)
            out=[v+h*j for v,j in zip(out,jump)]
    return [box(v) for v in out]
va,vb=velocity(a),velocity(b);v1a,v1b=velocity(a,1),velocity(b,1)
v2a,v2b=velocity(a,2),velocity(b,2);v3a,v3b=velocity(a,3),velocity(b,3)
Iv=vsub(v1a,v1b)
Jv=vsub(vscale(vadd(v1a,v1b),x),vsub(v2a,v2b))
Kv=vscale(vsub(vscale(vadd(v2a,v2b),x),vsub(v3a,v3b)),2)
A,B,D=([[F(z) for z in row] for row in prof['matrices'][k]] for k in ('A','B','D'))
Ar,Br,Dr=A[idx],B[idx],D[idx]
# Integration by parts eliminates all global nominal primitives and keeps
# the sole numerical integral local: Iu=integral_b^a u_s.
deltaU=vsub(ua,ub);sumU=vadd(ua,ub)
b1=scale(vsum([scale(dot(Ar,sumU),x),scale(dot(Br,deltaU),x),dot([x*d-aa for d,aa in zip(Dr,Ar)],Iu),scale(dot(Ar,Iv),x),dot(Br,Jv),scale(dot(Dr,Kv),F(1,2))]),F(1,4))
# Endpoint F and G in the exact registered ridge source.
Fa=vscale(vadd([dot([aa+bb for aa,bb in zip(A[i],B[i])],vadd(upa,va)) for i in range(23)],[dot(D[i],vadd(ua,v1a)) for i in range(23)]),F(1,2))
Gb=vscale(vadd([dot([bb-aa for aa,bb in zip(A[i],B[i])],vsub(upb,vb)) for i in range(23)],[dot(D[i],vsub(ub,v1b)) for i in range(23)]),F(1,2))
b1x=add(scale(add(Fa[idx],Gb[idx]),x/2),scale(vsum([dot(Br,deltaU),dot(Ar,Iv),dot(Dr,Iu)]),F(1,4)))
# integral(G-F)=-A deltaU-B Iv-D integral(V), integral(V)=Jv-independent relation.
# Integral V = v2a-v2b, which uses the primitive V=v1.
b1r=add(scale(sub(Fa[idx],Gb[idx]),x/2),scale(neg(vsum([dot(Ar,deltaU),dot(Br,Iv),dot(Dr,vsub(v2a,v2b))])),F(1,4)))
b0=scale(add(add(ua[idx],ub[idx]),Iv[idx]),F(1,2))
b0x=scale(add(sub(upa[idx],upb[idx]),add(va[idx],vb[idx])),F(1,2))
b0r=scale(add(add(upa[idx],upb[idx]),sub(va[idx],vb[idx])),F(1,2))
Vs={k:add(z,w) for k,z,w in [('value',b0,b1),('dx',b0x,b1x),('dr',b0r,b1r)]}
need(all(z[0]<=z[1] for z in list(Vs.values())+[b1,b1x,b1r]),'ordered output intervals')
# Independently recompute the published PR92 analytic transfer ledger for
# THIS exact integral-defined candidate from the same pinned reports. This
# verifies its rational arithmetic; the analytic PDE/transfer proof is inherited.
eU,eUr,eV=map(F,report['P1_H_error_caps'])
need([eU,eUr,eV]==[F(response['ledger'][key]) for key in ('epsilon_u_H','epsilon_du_H','epsilon_v_H')], 'prior P1 interpolation ledger agreement')
Eold=F(response['ledger']['total_X_error_upper'])
nA,nB,nD,nM=(F(response['constant_matrix_H_operator_upper'][key]) for key in ('A','B','D','M'))
blocks=[[i for i,label in enumerate(labels) if label.endswith(':'+degree)] for degree in ('0','2','4','+','-')]
need(sum(map(len,blocks))==23, 'operator block partition')
for key,cap in zip(('A','B','D','M'),(nA,nB,nD,nM)):
    matrix=[[F(v) for v in row] for row in prof['matrices'][key]]
    for i in range(23):
        for j in range(23):
            if not any(i in block and j in block for block in blocks): need(matrix[i][j]==0,'off-block operator entry')
    for block in blocks:
        need(sum(weights[i]*matrix[i][j]**2/weights[j] for i in block for j in block)<=cap**2, 'operator H block cap')
velocity_jump_norm2=max(sum(weights[i]*jump[i]**2 for i in range(23)) for knot,jump in jumps)
ceiling=math.isqrt(velocity_jump_norm2.numerator//velocity_jump_norm2.denominator)
if ceiling**2*velocity_jump_norm2.denominator<velocity_jump_norm2.numerator:ceiling+=1
av=F(35,256)*dv*ceiling
auH=M2H*du**2/18;buH=M2H*35*du/128
free_transfer=max((auH+F(3,2)*av)/3,buH+av)
feedback_transfer=F(3,2)*((nA+nB)*(eUr+buH+av)+nD*(eU+auH+F(3,2)*av))+nM*(eU+auH+av/4)
Ecomputed=Eold+free_transfer+feedback_transfer
E=F(2184*10**25)
need(Ecomputed<E, 'documented tighter candidate transfer bound')
need(x+r-1==F(1,2),'zero-wedge point factor')
U={k:widen(z,E*(x+r-1) if k=='value' else E) for k,z in Vs.items()}
need(U['value'][0]>0,'actual coefficient value lower bound not positive')
need(U['dr'][1]<0,'actual coefficient radial derivative upper bound not negative')
need(F(417*10**26)<U['value'][0] and U['value'][1]<F(636*10**26),'documented strict actual value caps')
need(F(-5288*10**28)<U['dr'][0] and U['dr'][1]<F(-5283*10**28),'documented strict actual radial caps')
# Every 23-mode endpoint and integral is already enclosed, so assemble the
# complete point table without any additional nominal evaluations.
def sqrt_upper(q):
    need(q>=0,'sqrt nonnegative')
    grid=1<<64;k=math.isqrt(q.numerator*grid*grid//q.denominator)
    if k*k*q.denominator<q.numerator*grid*grid:k+=1
    return F(k,grid)
point_table=[]
for i,label in enumerate(labels):
    ar,br,dr=A[i],B[i],D[i]
    z1=scale(vsum([scale(dot(ar,sumU),x),scale(dot(br,deltaU),x),dot([x*d-aa for d,aa in zip(dr,ar)],Iu),scale(dot(ar,Iv),x),dot(br,Jv),scale(dot(dr,Kv),F(1,2))]),F(1,4))
    z1x=add(scale(add(Fa[i],Gb[i]),x/2),scale(vsum([dot(br,deltaU),dot(ar,Iv),dot(dr,Iu)]),F(1,4)))
    z1r=add(scale(sub(Fa[i],Gb[i]),x/2),scale(neg(vsum([dot(ar,deltaU),dot(br,Iv),dot(dr,vsub(v2a,v2b))])),F(1,4)))
    z0=scale(add(add(ua[i],ub[i]),Iv[i]),F(1,2))
    z0x=scale(add(sub(upa[i],upb[i]),add(va[i],vb[i])),F(1,2))
    z0r=scale(add(add(upa[i],upb[i]),sub(va[i],vb[i])),F(1,2))
    candidate={key:add(z,w) for key,z,w in [('value',z0,z1),('dx',z0x,z1x),('dr',z0r,z1r)]}
    factor=sqrt_upper(1/weights[i]);actual={key:widen(z,E*factor*(x+r-1) if key=='value' else E*factor) for key,z in candidate.items()}
    need(all(z[0]<=z[1] for z in list(candidate.values())+list(actual.values())),'full23 point interval order')
    if i==idx:need(candidate==Vs and actual==U,'full23 primary row consistency')
    if i in (21,22):need(all(z==(0,0) for z in candidate.values()),'candidate TT zero')
    point_table.append({'label':label,'B0':{key:encode(z) for key,z in [('value',z0),('dx',z0x),('dr',z0r)]},'B1':{key:encode(z) for key,z in [('value',z1),('dx',z1x),('dr',z1r)]},'candidate':{key:encode(z) for key,z in candidate.items()},'conditional_actual_formal_response':{key:encode(z) for key,z in actual.items()},'error_component_factor_upper':str(factor)})
ctx.reauthenticate()
print(json.dumps({'status':'PASS','point':[str(x),str(r)],'component':'w:0','candidate':'tighter nominal u0 convolution plus smoothed P1 velocity; B0+B1 with frozen constant operators','du':str(du),'dv':str(dv),'local_integral_panels':N,'endpoint_method':'pinned nominal interval evaluator plus proved percomponent convolution remainder','B0':{k:encode(z) for k,z in [('value',b0),('dx',b0x),('dr',b0r)]},'B1':{k:encode(z) for k,z in [('value',b1),('dx',b1x),('dr',b1r)]},'candidate_intervals':{k:encode(z) for k,z in Vs.items()},'conditional_actual_formal_response_intervals':{k:encode(z) for k,z in U.items()},'C1_error_readable_upper':str(E),'C1_error_exact_recomputed_upper':str(Ecomputed),'prior_transfer_checker_reference_sha256':'6d1be89a615280a8e4dd9a5e3b224f0da889d1813190a2ece8e46c26641d195e','zero_wedge_value_error_factor':'1/2','actual_formal_value_positive':True,'actual_formal_dr_negative':True,'strict_primary_value_caps':['417*10^26','636*10^26'],'strict_primary_radial_caps':['-5288*10^28','-5283*10^28'],'full23_point_table':point_table,'scope':'one complete full23 nominal point table, primary scalar-mean sign bounds conditional on inherited certificates and proved candidate transfer; not exact nonlinear spacetime, C2 error, pointwise Einstein residual or physical observable'},sort_keys=True))
