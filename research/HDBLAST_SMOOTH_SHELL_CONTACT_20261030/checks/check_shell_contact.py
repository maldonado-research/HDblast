#!/usr/bin/env python3
"""Interval evaluation of the tighter nominal B0+B1 at an active shell point.
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

ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--panels',type=int,default=64);ap.add_argument('--source-report',type=Path,default=Path(__file__).resolve().parents[1]/'reference/AFFINE_SOURCE.json');args=ap.parse_args()
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
x,r=F(9,8),F(0);du=F(1,2**30);dv=F(1,2048)
need(1<x<F(5,4) and x+r<F(7,4),'active shell point in D+')
parity=list(map(int,prof['parity']['u']));need(len(parity)==23 and set(parity)=={-1,1},'parity shape')
need(parity==[-1 if s.startswith(('L:','v:')) else 1 for s in labels],'registered parity')
au=[m*du*du/18 for m in M2];bu=[m*35*du/128 for m in M2]
values,_=ctx.profiles(x)
u=[widen(pointbox(values[0][i]),au[i]) for i in range(23)]
p=[widen(pointbox(values[1][i]),bu[i]) for i in range(23)]
# Local range integration; split exactly at the only nominal analytic join.
# Below c the smoothed initial displacement is identically zero.
c=F(65,64)-du;join=F(66,64)
U1=[box(0) for _ in range(23)];U2=[box(0) for _ in range(23)]
for left0,right0 in [(F(65,64),join),(join,x)]:
    for k in range(args.panels):
        left=left0+(right0-left0)*k/args.panels
        right=left0+(right0-left0)*(k+1)/args.panels
        vals,_=ctx.profiles(ctx.I(left,right))
        U1=vadd(U1,[scale(pointbox(vals[0][i]),right-left) for i in range(23)])
        kernelmass=((x-left)**2-(x-right)**2)/2
        U2=vadd(U2,[scale(pointbox(vals[0][i]),kernelmass) for i in range(23)])
U1=[widen(U1[i],(x-c)*au[i]) for i in range(23)]
U2=[widen(U2[i],(x-c)**2*au[i]/2) for i in range(23)]
# Exact smoothed velocity primitives, with the registered zero collar.
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
    for j,c0 in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            m=2*j+k;out+=c0*math.comb(n,k)*s**(n-k)*(-1)**k*(s**(m+1)-(-1)**(m+1))/F(m+1)
    return F(35,32)*dv**n*out/F(math.factorial(n))
def velocity(z,order=0):
    out=[F(0)]*23
    for knot,jump in jumps:
        if z>knot-dv:
            h=hinge(z-knot,order+1)
            out=[v+h*j for v,j in zip(out,jump)]
    return [box(v) for v in out]
v,V1,V2,V3=[velocity(x,k) for k in range(4)]
A,B,D,M=([list(map(F,row)) for row in prof['matrices'][key]] for key in ('A','B','D','M'))
need(all(prof['matrices'][key]==response['constant_matrices'][key] for key in ('A','B','D','M')),'candidate frozen matrix agreement')
for i in range(23):
    if parity[i]==-1:need(all(c0==0 for c0 in M[i]),'M Dirichlet row')
    for j in range(23):
        if parity[j]==-1:need(M[i][j]==0,'M Dirichlet column')
def mv(m,z):return [dot(row,z) for row in m]
def refl(z):return [scale(v0,s) for s,v0 in zip(parity,z)]
def PN(z):return [v0 if s==1 else box(0) for s,v0 in zip(parity,z)]
def PD(z):return [v0 if s==-1 else box(0) for s,v0 in zip(parity,z)]
plus=[[a+b for a,b in zip(aa,bb)] for aa,bb in zip(A,B)]
minus=[[b-a for a,b in zip(aa,bb)] for aa,bb in zip(A,B)]
F0=vscale(vadd(mv(plus,vadd(p,v)),mv(D,vadd(u,V1))),F(1,2))
F1=vscale(vadd(mv(plus,vadd(u,V1)),mv(D,vadd(U1,V2))),F(1,2))
F2=vscale(vadd(mv(plus,vadd(U1,V2)),mv(D,vadd(U2,V3))),F(1,2))
G1plus=vscale(vadd(mv(minus,vsub(u,V1)),mv(D,vsub(U1,V2))),F(1,2))
G1minus=vscale(vsub(mv(minus,refl(vadd(u,V1))),mv(D,refl(vadd(U1,V2)))),F(1,2))
G2plus=vscale(vadd(mv(minus,vsub(U1,V2)),mv(D,vsub(U2,V3))),F(1,2))
G2minus=vscale(vsub(mv(D,refl(vadd(U2,V3))),mv(minus,refl(vadd(U1,V2)))),F(1,2))
# Complete physical source is coupled before reflection. These are the
# physical shell reduction of the reflected Duhamel kernel, not reflection
# of uncoupled input profiles through matrices.
H=vadd(vsub(vscale(F1,x),F2),vscale(vadd(G2plus,G2minus),F(1,2)))
pp=vadd(vscale(F0,x),vscale(vsub(G1plus,G1minus),F(1,2)))
K={'value':PN(H),'dx':PN(pp),'dr':PD(pp)}
lift_trace=mv(M,PN(vadd(U1,V2)))
lift_derivative=mv(M,PN(vadd(u,V1)))
Lift={'value':vscale(lift_trace,-1),'dx':vscale(lift_derivative,-1),'dr':lift_derivative}
B1={key:vadd(K[key],Lift[key]) for key in K}
B0={'value':PN(vadd(u,V1)),'dx':PN(vadd(p,v)),'dr':PD(vadd(p,v))}
candidate={key:vadd(B0[key],B1[key]) for key in K}
homogeneous_defect=vscale(mv(M,PN(B1['value'])),-1)
# A second arithmetic path computes PN V_r-M PN V using interval matrix
# operations on the assembled candidate. The identity avoids cancellation.
direct_defect=vsub(PN(candidate['dr']),mv(M,PN(candidate['value'])))
need(all(max(a[0],b[0])<=min(a[1],b[1]) for a,b in zip(homogeneous_defect,direct_defect)),'boundary defect identity overlap')
for i in range(23):
    need(all(z[i][0]<=z[i][1] for z in list(candidate.values())+list(B1.values())),'ordered candidate intervals')
    if parity[i]==-1:need(candidate['value'][i]==candidate['dx'][i]==(0,0),'homogeneous Dirichlet value and x derivative')
    else:need(K['dr'][i]==(0,0),'reflected Neumann kernel normal derivative')
    if i in (21,22):need(all(z[i]==(0,0) for z in candidate.values()),'TT candidate zero')

# Recompute the inherited exact C1 transfer ledger.
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
def sqrt_upper(q):
    need(q>=0,'sqrt nonnegative')
    grid=1<<64;k=math.isqrt(q.numerator*grid*grid//q.denominator)
    if k*k*q.denominator<q.numerator*grid*grid:k+=1
    return F(k,grid)
need(x-1==F(1,8),'shell zero-wedge value error factor')
actual={key:[widen(z,E*sqrt_upper(1/weights[i])*(x-1 if key=='value' else 1)) for i,z in enumerate(row)] for key,row in candidate.items()}
source=json.loads(pin(args.source_report,'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e'))
need([row['label'] for row in source['rows']]==labels,'source labels')
G=[tuple(map(F,row['interval'])) for row in source['rows']]
need(all(z[0]<=z[1] for z in G),'source intervals ordered')
frozen_full_defect=vsub(homogeneous_defect,G)
Mactual=[[tuple(map(F,z)) for z in row] for row in response['actual_coefficient_boxes']['M']]
need(len(Mactual)==23 and all(len(row)==23 for row in Mactual),'actual M shape')
need(all(z[0]<=z[1] for row in Mactual for z in row),'actual M ordered intervals')
def mul(a,b):
    z=[x0*y0 for x0 in a for y0 in b];return(min(z),max(z))
def imv(matrix0,z):return [vsum(mul(a,b) for a,b in zip(row,z)) for row in matrix0]
deltaM=[[sub(a,box(b)) for a,b in zip(row,nom)] for row,nom in zip(Mactual,M)]
for i in range(23):
    for j in range(23):
        if parity[i]==-1 or parity[j]==-1:need(Mactual[i][j]==(0,0),'actual M Dirichlet row/column')
M_uncertainty_correction=imv(deltaM,PN(candidate['value']))
full_defect=vsub(frozen_full_defect,M_uncertainty_correction)
# Use the exact inherited shell laws to constrain the actual formal response,
# rather than transferring the candidate normal derivative alone. These laws
# refer to the actual formal boundary problem, not to this candidate.
for i,s in enumerate(parity):
    if s==-1:actual['value'][i]=actual['dx'][i]=box(0)
robin_dr=vadd(imv(Mactual,PN(actual['value'])),G)
actual_dr_boundary=[]
for i,s in enumerate(parity):
    if s==-1:actual_dr_boundary.append(actual['dr'][i]);continue
    a0,b0=actual['dr'][i],robin_dr[i]
    z=(max(a0[0],b0[0]),min(a0[1],b0[1]))
    need(z[0]<=z[1],'C1 and exact Robin radial intervals compatible')
    actual_dr_boundary.append(z)
need(actual['value'][12][0]>0 and actual['dx'][12][0]>0 and actual_dr_boundary[12][0]>0,'conditional actual formal w0 shell first jet positive')
need(full_defect[12][0]>0,'candidate full affine w0 boundary defect positive')
rows=[]
for i,label in enumerate(labels):
    rows.append({'label':label,'parity':parity[i],**{name:{key:encode(row[i]) for key,row in table.items()} for name,table in [('B0',B0),('reflected_bulk_first_feedback',K),('Robin_lift',Lift),('B1',B1),('candidate',candidate),('conditional_actual_formal_response',actual)]},'affine_shell_source':encode(G[i]),'full_affine_boundary_defect':encode(full_defect[i]),'frozen_full_affine_boundary_defect':encode(frozen_full_defect[i]),'actual_minus_frozen_M_times_candidate':encode(M_uncertainty_correction[i]),'conditional_actual_formal_Robin_radial':encode(robin_dr[i]) if parity[i]==1 else None,'conditional_actual_formal_dr_intersection':encode(actual_dr_boundary[i]),'homogeneous_boundary_defect':encode(homogeneous_defect[i]),'direct_homogeneous_boundary_defect':encode(direct_defect[i])})
ctx.reauthenticate()
print(json.dumps({'status':'PASS','point':[str(x),str(r)],'labels':labels,'parity':parity,'candidate':'tighter nominal displacement convolution plus smoothed P1 velocity; reflected B0+B1 including Robin lift with frozen constant operators','displacement_delta':str(du),'velocity_delta':str(dv),'local_integral_panels_per_analytic_branch':args.panels,'local_integral_analytic_branches':[['65/64','66/64'],['66/64','9/8']],'displacement_integral_method':'whole-panel interval ranges times exact nonnegative integral weights plus percomponent convolution remainder over smoothed support','full23_point_table':rows,'C1_error_readable_upper':str(E),'C1_error_exact_recomputed_upper':str(Ecomputed),'zero_wedge_value_error_factor':'1/8','actual_Robin_matrix':'interval boxes from pinned response report; distinct from frozen M0 in candidate','boundary_identity':'PN V_r-M0 PN V=-M0 PN B1 at shell; actual-law affine defect subtracts (M-M0)PN V and G','conditional_actual_formal_w0_value_dx_and_boundary_constrained_dr_positive':True,'candidate_w0_full_affine_boundary_defect_positive':True,'affine_source_report_sha256':'fec70c6c305ceb55b3c1b9917b4cbf5abb25c403e0399770877d2f3dc5b2bf2e','Dirichlet_slots_exact_zero_value_and_dx':True,'boundary_homogeneous_defect_identity_overlap_all23':True,'scope':'complete formal candidate at one active shell point with complete source reflected after matrix action and Robin lift; conditional C1 actual response intervals inherit earlier certificates; affine G is bounded in the candidate defect and exact Robin trace, but not inserted in candidate; no C2 error, pointwise Einstein residual, nonlinear solution or physical observable'},sort_keys=True))
