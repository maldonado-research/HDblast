#!/usr/bin/env python3
"""Independent shell first-feedback enclosure by norm integration, not ridge quadrature."""
from fractions import Fraction as Q
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,subprocess,sys

def need(c,m):
    if not c:raise RuntimeError(m)
def read(p,h):
    b=p.read_bytes();need(hashlib.sha256(b).hexdigest()==h,'INPUT_HASH:'+str(p));return b
def sqrtn(q):
    need(q>=0,'square root nonnegative');n=math.isqrt(q.numerator//q.denominator)
    if n*n*q.denominator<q.numerator:n+=1
    need(n*n>=q,'sqrt outward');return Q(n)
def n2(row,w):return sum((a*a*b for a,b in zip(row,w)),Q(0))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def mul(a,c):return min(a[0]*c,a[1]*c),max(a[0]*c,a[1]*c)
def wide(a,e):return a[0]-e,a[1]+e
def box(a):return a,a
def dot(row,vec):
    out=box(Q(0))
    for a,v in zip(row,vec):out=add(out,mul(v,a))
    return out

a=argparse.ArgumentParser();a.add_argument('--repo-root',type=Path,required=True);a.add_argument('--source-checker',type=Path);args=a.parse_args();root=args.repo_root.resolve()
b=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
p=json.loads(read(b/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
r=json.loads(read(b/'compute/response-normal/RESPONSE_REPORT.json','ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd'))
i=json.loads(read(b/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
init=b/'initial/full23_initial.py';read(init,'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1')
need(r['status']=='PASS_BOUNDED_FULL23_RESPONSE_X_ERROR' and i['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER','inherited status')
labels=p['labels'];w=list(map(Q,p['gram']));N=len(w)
need(N==23 and labels==r['labels']==i['labels'] and w==list(map(Q,r['gram']))==list(map(Q,i['gram'])) and min(w)>0,'space')
parity=[-1 if z.startswith(('L:','v:')) else 1 for z in labels];need(parity.count(-1)==5,'mixed split')
blocks=[[j for j,z in enumerate(labels) if z.endswith(':'+d)] for d in ('0','2','4','+','-')]
ops={k:Q(r['constant_matrix_H_operator_upper'][k]) for k in ('A','B','D','M')}
mat={k:[list(map(Q,row)) for row in p['matrices'][k]] for k in ops}
for name,cap in ops.items():
    M=mat[name]
    for j in range(N):
        for k in range(N):
            need(any(j in block and k in block for block in blocks) or M[j][k]==0,'off-block entry')
    for block in blocks:need(sum(w[j]*M[j][k]**2/w[k] for j in block for k in block)<=cap**2,'H operator cap')
need(all(mat['M'][j][k]==0 for j in range(N) if parity[j]<0 for k in range(N)),'no odd Robin rows')
# Actual shell M is enclosed separately from the candidate's frozen M0.
actualM=r['actual_coefficient_boxes']['M'];entryM=[[Q(z) for z in row] for row in r['coefficient_entry_error_caps']['M']]
deltaM=Q(r['coefficient_error_H_operator_upper']['M'])
for j in range(N):
    for k in range(N):
        lo,hi=map(Q,actualM[j][k]);need(lo<=hi and max(abs(lo-mat['M'][j][k]),abs(hi-mat['M'][j][k]))<=entryM[j][k],'actual minus frozen shell entry')
        need(any(j in block and k in block for block in blocks) or entryM[j][k]==0,'actual shell off-block error')
for block in blocks:need(sum(w[j]*entryM[j][k]**2/w[k] for j in block for k in block)<=deltaM**2,'actual shell error H cap')
need(deltaM>0,'nonzero frozen-vs-actual coefficient error')
nodes=list(map(Q,p['nodes']));tab={name:[list(map(Q,row)) for row in p['profiles'][name]] for name in ('u','du','v')}
need(len(nodes)==len(i['interpolation_cells'])+1,'initial cell roster')
x=Q(9,8);rad=Q(0);du=Q(1,2**30);dv=Q(1,2048);h=x-1
need(1<x<Q(5,4) and x<Q(7,4),'active shell within Dplus')
need(nodes[1]-max(du,dv)>1 and all(v==0 for t in tab.values() for row in t[:2] for v in row),'flat collar')
# Local nominal cap uses P1 convexity PLUS authenticated continuous interpolation
# error on every intersected cell. It is independent of new local ridge integration.
ucap=[Q(0)]*N;pcap=[Q(0)]*N;vcap=[Q(0)]*N
for k,cell in enumerate(i['interpolation_cells']):
    need(list(map(Q,cell['r']))==nodes[k:k+2],'cell geometry')
    if nodes[k]>x+max(du,dv):break
    for j in range(N):
        ucap[j]=max(ucap[j],max(abs(tab['u'][k][j]),abs(tab['u'][k+1][j]))+Q(cell['component_error_caps'][0][j]))
        pcap[j]=max(pcap[j],max(abs(tab['du'][k][j]),abs(tab['du'][k+1][j]))+Q(cell['component_error_caps'][1][j]))
        vcap[j]=max(vcap[j],abs(tab['v'][k][j]),abs(tab['v'][k+1][j]))
U0,U1,Vcap=(sqrtn(n2(c,w)) for c in (ucap,pcap,vcap))
Fcap=(ops['A']+ops['B'])*(U1+Vcap)+ops['D']*(U0+x*Vcap)
# Full coupled source acts on physical z>=0, then output reflection is H-isometric.
# Area of active physical cone: ∫_{alpha=1}^x alpha d alpha=(x²−1)/2.
# Source on the outgoing cone edge can be active throughout all s∈[0,x].
Kvalue=(x*x-1)*Fcap/2;Kder=x*Fcap
tracecap=U0+h*Vcap;Lvalue=h*ops['M']*tracecap;Lder=ops['M']*tracecap
B1value=Kvalue+Lvalue;B1der=Kder+Lder
need(Kvalue==(x*x-1)*Fcap/2 and Kvalue>h*h*Fcap/2,'active cone area control')
# Independent exact positive convolution estimate for free trace. Bound primitive
# mollification by ||v||_infty E|s|, avoiding all new hinge/primitive routines.
spec=importlib.util.spec_from_file_location('coarse_shell_nominal',init);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
ctx=mod.Context(root)
values,_=ctx.profiles(x);M2=list(map(Q,i['nominal_second_derivative_component_caps'][0]))
ua=[wide((values[0][j].lo,values[0][j].hi),M2[j]*du*du/18) for j in range(N)]
pa=[wide((values[1][j].lo,values[1][j].hi),35*M2[j]*du/128) for j in range(N)]
Vint=[Q(0)]*N
for k,(lo,hi) in enumerate(zip(nodes,nodes[1:])):
    if lo>=x:break
    right=min(hi,x);width=right-lo
    for j in range(N):
        slope=(tab['v'][k+1][j]-tab['v'][k][j])/(hi-lo)
        Vint[j]+=width*tab['v'][k][j]+width*width*slope/2
# Global v norm cap over complete P1 plus last affine continuation is unnecessary:
# x+dv is inside mesh, so local convexity cap above covers every convolution point.
vsbox=[];v1box=[]
for j in range(N):
    point=tab['v'][nodes.index(x)][j]
    local_slopes=[]
    for k in range(len(nodes)-1):
        if nodes[k]<x+dv and nodes[k+1]>x-dv:
            local_slopes.append(abs((tab['v'][k+1][j]-tab['v'][k][j])/(nodes[k+1]-nodes[k])))
    vsbox.append(wide(box(point),35*dv*max(local_slopes)/128))
    v1box.append(wide(box(Vint[j]),35*dv*vcap[j]/128))
zero=box(Q(0));b0value=[add(ua[j],v1box[j]) if parity[j]>0 else zero for j in range(N)]
b0der=[add(pa[j],vsbox[j]) for j in range(N)]
# The new wider smoothed-velocity primitive estimate is independent of exact
# polynomial hinge convolution; it still encloses the same integral-defined B0.
# Recompute the same inherited C1 transfer with the exact outward jump norm.
slopes=[[(tab['v'][k+1][j]-tab['v'][k][j])/(nodes[k+1]-nodes[k]) for j in range(N)] for k in range(len(nodes)-1)]
jumps=[[b-a for a,b in zip(a,b)] for a,b in zip(slopes,slopes[1:])]
av=35*dv*sqrtn(max(n2(row,w) for row in jumps))/256
M2H=Q(i['nominal_second_derivative_H_caps'][0]);au=M2H*du**2/18;bu=35*M2H*du/128
EU,EUr,EV=map(Q,i['P1_H_error_caps']);need([EU,EUr,EV]==[Q(r['ledger'][k]) for k in ('epsilon_u_H','epsilon_du_H','epsilon_v_H')],'old interpolation ledger')
E=Q(r['ledger']['total_X_error_upper'])+max((au+Q(3,2)*av)/3,bu+av)+Q(3,2)*((ops['A']+ops['B'])*(EUr+bu+av)+ops['D']*(EU+au+Q(3,2)*av))+ops['M']*(EU+au+av/4)
need(E<Q(2184*10**25),'C1 transfer readable cap')
sourcepath=args.source_checker or Path(__file__).with_name('check_affine_shell_source.py')
sourcepin='7bb0179f76e43b5d9bbaf6ef4df4041aaf232326bf55776cc9566880b8f2a089'
read(sourcepath,sourcepin)
cmd=[sys.executable,'-I','-S','-B']+(['-O'] if sys.flags.optimize else [])+[str(sourcepath),'--repo-root',str(root)]
source=json.loads(subprocess.check_output(cmd,text=True))
need(source['status']=='PASS' and source['point']==[str(x),str(rad)] and [z['label'] for z in source['rows']]==labels,'source scope/labels')
rows=[]
for j,label in enumerate(labels):
    factor=sqrtn(1/w[j]);need(factor*factor*w[j]>=1,'component norm factor')
    if parity[j]>0:
        cand={'value':wide(b0value[j],factor*B1value),'dx':wide(b0der[j],factor*B1der),'dr':dot(mat['M'][j],b0value)}
    else:
        cand={'value':zero,'dx':zero,'dr':wide(b0der[j],factor*B1der)}
    actual={key:wide(v,factor*E*(h if key=='value' else 1)) for key,v in cand.items()}
    if parity[j]<0:actual['value']=zero;actual['dx']=zero # inherited homogeneous Dirichlet trace
    robin_actual=wide(dot(mat['M'][j],b0value),factor*(deltaM*tracecap+(ops['M']+deltaM)*(B1value+h*E))+Q(source['rows'][j]['abs_upper'])) if parity[j]>0 else None
    rows.append({'shell_Robin_actual_radial_derivative':list(map(str,robin_actual)) if robin_actual is not None else None,'label':label,'parity':parity[j],'B0_value':list(map(str,b0value[j])),'B0_active_derivative':list(map(str,b0der[j])),'candidate':{key:list(map(str,v)) for key,v in cand.items()},'conditional_actual_formal_response':{key:list(map(str,v)) for key,v in actual.items()}})
mean=rows[12];need(labels[12]=='w:0' and w[12]==1,'scalar mean')
need(Q(mean['conditional_actual_formal_response']['value'][0])>0,'positive actual scalar mean')
need(Q(mean['conditional_actual_formal_response']['dx'][0])>0,'positive actual scalar mean time derivative')
need(Q(mean['shell_Robin_actual_radial_derivative'][0])>0,'positive actual scalar mean radial derivative from Robin law')
ctx.reauthenticate()
print(json.dumps({'status':'PASS','point':[str(x),str(rad)],'method':'local cone norm bound plus independent nominal endpoint and Lipschitz velocity convolution errors; no new ridge quadrature import','local_U0_H_cap':str(U0),'local_U1_H_cap':str(U1),'local_V_H_cap':str(Vcap),'physical_source_H_cap':str(Fcap),'reflected_bulk_value_H_cap':str(Kvalue),'reflected_bulk_derivative_H_cap':str(Kder),'Robin_lift_value_H_cap':str(Lvalue),'Robin_lift_derivative_H_cap':str(Lder),'B1_value_H_cap':str(B1value),'B1_derivative_H_cap':str(B1der),'C1_error_exact_recomputed_upper':str(E),'zero_wedge_value_error_factor':str(h),'full23_point_table':rows,'scalar_mean_value_positive':True,'scalar_mean_dx_positive':True,'scalar_mean_robin_dr_positive':True,'source_checker_sha256':sourcepin,'actual_shell_M_minus_frozen_M0_H_cap':str(deltaM),'actual_shell_M_H_cap':str(ops['M']+deltaM),'actual_radial_reconstruction_norm_error_H_cap':str(deltaM*tracecap+(ops['M']+deltaM)*(B1value+h*E)),'boundary_identity':'candidate PN V_r=M0 PN B0; actual PN U_r=Mactual PN U+G; PD V=0','scope':'conditional registered formal coefficients and candidate boundary identity; no full nonlinear solution, numerical C2 error or physical observable'},sort_keys=True))
