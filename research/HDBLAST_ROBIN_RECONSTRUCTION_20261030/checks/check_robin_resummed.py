#!/usr/bin/env python3
"""Enclose the integral-defined boundary-resummed smooth candidate at one shell point.

The actual constant Robin matrix and actual affine shell source are retained.
The bulk source is still frozen at the free smooth wave, not at the new field.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib.util, json, math


def need(ok, message):
    if not ok: raise RuntimeError(message)
def pin(path, digest):
    data=path.read_bytes()
    need(hashlib.sha256(data).hexdigest()==digest, 'INPUT_HASH:'+str(path))
    return data
def box(z): return F(z),F(z)
def add(a,b): return a[0]+b[0],a[1]+b[1]
def neg(a): return -a[1],-a[0]
def sub(a,b): return add(a,neg(b))
def scale(a,c): return min(a[0]*c,a[1]*c),max(a[0]*c,a[1]*c)
def mul(a,b):
    out=[x*y for x in a for y in b]
    return min(out),max(out)
def widen(a,e): return a[0]-e,a[1]+e
def interval(a):
    out=tuple(map(F,a));need(len(out)==2 and out[0]<=out[1], 'ordered interval')
    return out
def vsum(rows):
    out=box(0)
    for row in rows: out=add(out,row)
    return out
def vadd(a,b): return [add(x,y) for x,y in zip(a,b)]
def vsub(a,b): return [sub(x,y) for x,y in zip(a,b)]
def vscale(a,c): return [scale(z,c) for z in a]
def mv(a,v): return [vsum(scale(z,c) for c,z in zip(row,v)) for row in a]
def imv(a,v): return [vsum(mul(c,z) for c,z in zip(row,v)) for row in a]
def encode(a): return list(map(str,a))
def sqrt_upper(q):
    need(q>=0,'nonnegative norm square')
    grid=1<<64;k=math.isqrt(q.numerator*grid*grid//q.denominator)
    if k*k*q.denominator<q.numerator*grid*grid:k+=1
    need(F(k,grid)**2>=q, 'outward square root')
    return F(k,grid)


ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',type=Path,required=True)
ap.add_argument('--prior-capsule',type=Path)
ap.add_argument('--panels',type=int,default=64)
args=ap.parse_args();need(1<=args.panels<=4096,'panel count')
root=args.repo_root.resolve();base=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
prior=args.prior_capsule or root/'research/HDBLAST_SMOOTH_SHELL_CONTACT_20261030'
old=json.loads(pin(prior/'reference/SHELL.json','9de23a17dfb96361e4d5ad46f31b60cb0da458492d6e21d26c9dd6890e0b1d64'))
coarse=json.loads(pin(prior/'reference/COARSE.json','2b8daad3410478d180ba46cffd07609ff57fce8f8e3fd2f77256ea4a0fe2214c'))
source=json.loads(pin(prior/'reference/AFFINE_SOURCE.json','fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e'))
initial_path=base/'initial/full23_initial.py'
pin(initial_path,'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1')
initial=json.loads(pin(base/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
profile=json.loads(pin(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
response=json.loads(pin(base/'compute/response-normal/RESPONSE_REPORT.json','ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd'))
labels=profile['labels'];weights=list(map(F,profile['gram']));dim=len(labels)
need(dim==23 and labels==initial['labels']==response['labels']==old['labels'], 'registered labels')
need([r['label'] for r in coarse['full23_point_table']]==labels==[r['label'] for r in source['rows']], 'inherited table labels')
need(weights==list(map(F,initial['gram']))==list(map(F,response['gram'])) and min(weights)>0,'registered Gram')
need(initial['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER' and response['status']=='PASS_BOUNDED_FULL23_RESPONSE_X_ERROR', 'inherited scope')
need(all(d['status']=='PASS' and d['point']==['9/8','0'] for d in (old,coarse,source)), 'successful inherited point tables')
P=list(map(int,profile['parity']['u']))
need(P==[-1 if z.startswith(('L:','v:')) else 1 for z in labels], 'mixed parity')
PN=lambda row:[z if s==1 else box(0) for s,z in zip(P,row)]
PD=lambda row:[z if s==-1 else box(0) for s,z in zip(P,row)]
reflect=lambda row:[scale(z,s) for s,z in zip(P,row)]
x=F(9,8);d=x-1;N=3;du=F(1,2**30);dv=F(1,2048)
M2=list(map(F,initial['nominal_second_derivative_component_caps'][0]))
need(initial['derivative_convention']['du']=='partial_r U' and initial['derivative_convention']['v']=='partial_x U; x=t/ell', 'initial derivative convention')
need(initial['derivative_convention']['second_derivatives'][0]=='partial_r^2 U', 'second derivative convention')
need(len(M2)==dim and min(M2)>=0, 'nominal second derivative caps')
alpha=[m*du**2/18 for m in M2];beta=[m*35*du/128 for m in M2]
spec=importlib.util.spec_from_file_location('robin_nominal_initial',initial_path)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
ctx=mod.Context(root)  # authenticates all scientific upstream roots before imports
vals,_=ctx.profiles(x)
u=[widen((z.lo,z.hi),a) for z,a in zip(vals[0],alpha)]
p=[widen((z.lo,z.hi),b) for z,b in zip(vals[1],beta)]
# Uj=J^j u_s. Each nonnegative polynomial weight is integrated exactly over
# whole-panel nominal interval ranges; convolution remainders are added on
# the entire smoothed support, including the tiny part below65/64.
U=[u]+[[box(0) for _ in labels] for j in range(N+2)]
c=F(65,64)-du
for left0,right0 in ((F(65,64),F(66,64)),(F(66,64),x)):
    for k in range(args.panels):
        left=left0+(right0-left0)*k/args.panels
        right=left0+(right0-left0)*(k+1)/args.panels
        vals,_=ctx.profiles(ctx.I(left,right));z=[(v.lo,v.hi) for v in vals[0]]
        for j in range(1,N+3):
            mass=((x-left)**j-(x-right)**j)/math.factorial(j)
            U[j]=vadd(U[j],vscale(z,mass))
for j in range(1,N+3):
    U[j]=[widen(z,a*(x-c)**j/math.factorial(j)) for z,a in zip(U[j],alpha)]
nodes=list(map(F,profile['nodes']));vv=[list(map(F,row)) for row in profile['profiles']['v']]
need(all(z==0 for row in vv[:2] for z in row) and nodes[1]-max(du,dv)>1, 'flat initial collar')
slopes=[[(b-a)/(hi-lo) for a,b in zip(va,vb)] for lo,hi,va,vb in zip(nodes,nodes[1:],vv,vv[1:])]
need(not any(slopes[0]), 'zero first velocity slope')
jumps=[(nodes[k],[b-a for a,b in zip(slopes[k-1],slopes[k])]) for k in range(1,len(slopes))]
def moment(k):
    return F(0) if k%2 else F(35,32)*sum(F(2*z,2*j+k+1) for j,z in enumerate((1,-3,3,-1)))
def hinge(z,n):
    if z<=-dv:return F(0)
    if z>=dv:return sum(F(math.comb(n,k),math.factorial(n))*z**(n-k)*(-dv)**k*moment(k) for k in range(n+1))
    s=z/dv;out=F(0)
    for j,z0 in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            power=2*j+k
            out+=z0*math.comb(n,k)*s**(n-k)*(-1)**k*(s**(power+1)-(-1)**(power+1))/F(power+1)
    return F(35,32)*dv**n*out/math.factorial(n)
def velocity(order):
    out=[F(0)]*dim
    for knot,jump in jumps:
        if x>knot-dv:
            z=hinge(x-knot,order+1);out=[a+z*b for a,b in zip(out,jump)]
    return [box(z) for z in out]
V=[velocity(j) for j in range(N+4)]
A,B,D,M0=([list(map(F,row)) for row in profile['matrices'][key]] for key in ('A','B','D','M'))
need(all(profile['matrices'][key]==response['constant_matrices'][key] for key in ('A','B','D','M')), 'frozen matrices agree')
plus=[[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
minus=[[b-a for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def ridge_F(j):
    deriv=p if j==0 else U[j-1]
    return vscale(vadd(mv(plus,vadd(deriv,V[j])),mv(D,vadd(U[j],V[j+1]))),F(1,2))
def ridge_G(j,negative=False):
    need(j>=1,'only normalized positive-order G primitives used')
    if not negative:return vscale(vadd(mv(minus,vsub(U[j-1],V[j])),mv(D,vsub(U[j],V[j+1]))),F(1,2))
    # Normalized negative primitives preserve their parity inside matrix actions.
    out=vsub(mv(minus,reflect(vadd(U[j-1],V[j]))),mv(D,reflect(vadd(U[j],V[j+1]))))
    return vscale(out,F((-1)**(j-1),2))
H=[];h=[]
for j in range(N+1):
    z=vadd(vsub(vscale(ridge_F(j+1),x),vscale(ridge_F(j+2),j+1)),
           vscale(vadd(ridge_G(j+2),vscale(ridge_G(j+2,True),(-1)**j)),F(1,2)))
    H.append(PN(z));h.append(PN(vadd(vadd(U[j],V[j+1]),z)))
edge=vadd(vscale(ridge_F(0),x),vscale(vsub(ridge_G(1),ridge_G(1,True)),F(1,2)))
hx=PN(vadd(vadd(p,V[0]),edge));base_dr=PD(vadd(vadd(p,V[0]),edge))
M=[[interval(z) for z in row] for row in response['actual_coefficient_boxes']['M']]
need(len(M)==dim and all(len(row)==dim for row in M), 'actual M shape')
err=[list(map(F,row)) for row in response['coefficient_entry_error_caps']['M']]
deltaM=[[sub(z,box(a)) for z,a in zip(row,nom)] for row,nom in zip(M,M0)]
for i in range(dim):
    for j in range(dim):
        need(max(map(abs,deltaM[i][j]))<=err[i][j], 'actual coefficient error entry')
        if P[i]<0 or P[j]<0:need(M[i][j]==(0,0), 'actual Robin homogeneous Dirichlet row and column')
delta=F(response['coefficient_error_H_operator_upper']['M'])
nM0=F(response['constant_matrix_H_operator_upper']['M']);m=nM0+delta
blocks=[[j for j,z in enumerate(labels) if z.endswith(':'+degree)] for degree in ('0','2','4','+','-')]
for i in range(dim):
    for j in range(dim):
        if not any(i in block and j in block for block in blocks):need(M[i][j]==(0,0), 'actual Robin off-block zero')
for block in blocks:
    need(sum(weights[i]*M0[i][j]**2/weights[j] for i in block for j in block)<=nM0**2,'frozen M operator cap')
    need(sum(weights[i]*err[i][j]**2/weights[j] for i in block for j in block)<=delta**2,'coefficient error operator cap')
G=[interval(row['interval']) for row in source['rows']]
need(source['uniform_domain']=='0<=x<=3,r=0', 'source bounds cover full active history')
need(all(G[i]==(0,0) for i,s in enumerate(P) if s<0), 'actual Dirichlet source exactly zero')
Gnorm=F(source['source_H_norm_strict_integer_upper'])
tracecap=F(coarse['local_U0_H_cap'])+d*F(coarse['local_V_H_cap'])
Hcap=tracecap+F(coarse['reflected_bulk_value_H_cap'])
a=m*d;need(0<a<1, 'rational Volterra tail ratio')
htail=Hcap*a**(N+1)/math.factorial(N+1)/(1-a/F(N+2))
# The unresolved actual source is retained via its exact convolution integral.
# The first term -JG is bounded rowwise; the remaining terms use a Gram norm.
gtail=m*d*d*Gnorm/2/(1-a/3)
series=[box(0) for _ in labels]
terms=[]
for j in range(N+1):
    z=h[j]
    for k in range(j):z=imv(M,z)
    z=vscale(z,(-1)**j);terms.append(z);series=vadd(series,z)
factors=[sqrt_upper(1/w) for w in weights]
source_part=[];tail=[]
for i in range(dim):
    zero_row=not any(z!=(0,0) for z in M[i])
    tail.append(F(0) if zero_row else factors[i]*htail)
    source_part.append(widen(scale(G[i],-d),F(0) if zero_row else factors[i]*gtail))
Z=[widen(add(z,g),t) for z,g,t in zip(series,source_part,tail)]
# A cancellation-free comparison against the former first-feedback candidate.
JZ0=PN(vadd(U[1],V[2]))
correction=vsub(vscale(imv(deltaM,JZ0),-1),imv(M,H[1]))
for j in range(2,N+1):correction=vadd(correction,terms[j])
correction=[widen(add(z,g),t) for z,g,t in zip(correction,source_part,tail)]
oldrows=old['full23_point_table'];oldtrace=[interval(row['candidate']['value']) for row in oldrows]
via_old=vadd(oldtrace,correction)
for i in range(dim):
    z=(max(Z[i][0],via_old[i][0]),min(Z[i][1],via_old[i][1]))
    need(z[0]<=z[1], 'direct series and cancellation-free correction overlap')
    Z[i]=z
q=vadd(imv(M,Z),G)
W={'value':Z,'dx':vsub(hx,q),'dr':vadd(base_dr,q)}
for i,s in enumerate(P):
    if s<0:need(W['value'][i]==W['dx'][i]==correction[i]==(0,0), 'exact Dirichlet value/time trace/correction')
    need(all(z[i][0]<=z[i][1] for z in W.values()), 'ordered reconstructed jets')
need(labels[12]=='w:0' and all(W[k][12][0]>0 for k in W), 'positive reconstructed scalar-mean first jet')
if args.panels==64:
    need(F(745*10**18)<correction[12][0]<=correction[12][1]<F(1526*10**18), 'strict positive scalar-mean trace correction')
    for key,lo,hi in [('value',F(829037*10**25),F(829041*10**25)),
                      ('dx',F(197462*10**27),F(197463*10**27)),
                      ('dr',F(215406*10**22),F(215408*10**22))]:
        need(lo<W[key][12][0]<=W[key][12][1]<hi, 'readable reconstructed w0 '+key+' bounds')
need(htail<F(6*10**11) and gtail<F(1043,1000), 'readable Volterra tail bounds')
ctx.reauthenticate()
rows=[]
for i,label in enumerate(labels):
    rows.append({'label':label,'parity':P[i],
                 'trace_moments_h_j':[encode(row[i]) for row in h],
                 'finite_homogeneous_trace_series':encode(series[i]),
                 'homogeneous_series_tail_abs_upper':str(tail[i]),
                 'actual_affine_source_trace_contribution':encode(source_part[i]),
                 'trace_correction_from_old_candidate':encode(correction[i]),
                 'reconstructed_boundary_candidate':{k:encode(row[i]) for k,row in W.items()},
                 'boundary_law_rhs_MZ_plus_G':encode(q[i]),
                 'exact_affine_Robin_residual_by_analytic_identity':['0','0'] if P[i]>0 else None})
print(json.dumps({'status':'PASS','point':[str(x),'0'],'labels':labels,'parity':P,
    'candidate':'integral-defined W=B0+K f0+L(MZ+G); actual constant Robin coefficients and affine shell source retained; bulk f0 frozen at free smooth wave',
    'series_degree':N,'active_history_length':str(d),'local_integral_panels_per_branch':args.panels,
    'actual_Robin_M_H_operator_upper':str(m),'trace_h_H_uniform_cap':str(Hcap),
    'dimensionless_Volterra_ratio':str(a),'homogeneous_trace_series_H_tail_upper':str(htail),
    'affine_source_higher_trace_terms_H_tail_upper':str(gtail),
    'affine_source_method':'actual source convolution retained as an integral-defined function; enclosure uses uniform inherited source bounds, not a point evaluation',
    'source_H_norm_upper':str(Gnorm),'full23_point_table':rows,
    'conditional_boundary_problem':'exact actual affine Robin and homogeneous Dirichlet laws for integral-defined W, principal bulk forcing frozen at f0',
    'zero_residual_semantics':'symbolic identity for the integral-defined W with the same fixed actual M and inherited source G; not an interval subtraction of independently selected boxes and not the residual of the finite series alone',
    'regularity_scope':'C1 trace uses inherited continuous G and smooth h; any physical-side C2/corner claim also requires inherited C1 G and flat zero-corner compatibility, not merely the source magnitude bounds',
    'prior_capsule_reports_sha256':{'SHELL':'9de23a17dfb96361e4d5ad46f31b60cb0da458492d6e21d26c9dd6890e0b1d64','COARSE':'2b8daad3410478d180ba46cffd07609ff57fce8f8e3fd2f77256ea4a0fe2214c','AFFINE_SOURCE':'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e'},
    'scope':'conditional boundary-resummed smooth candidate at one shell sample, not a full second-order bulk solution, numerical C2 error, Einstein/constraint residual, nonlinear spacetime or physical observable'},sort_keys=True))
