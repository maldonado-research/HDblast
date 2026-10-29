#!/usr/bin/env python3
"""Independent local cone bound for the tight smooth candidate and target.
Evaluates B0 through a pinned interval checker; bounds B1 by cone norms.
No quantitative C2 error or Einstein residual is computed.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import zip_longest
import argparse, hashlib, json, math, subprocess, sys

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def read(p,pin):
    b=p.read_bytes();need(hashlib.sha256(b).hexdigest()==pin,'INPUT_HASH:'+p.name)
    return json.loads(b)
def sqrtup(q):
    need(q>=0,'negative squared norm')
    n=math.isqrt(q.numerator//q.denominator)
    if n*n*q.denominator<q.numerator:n+=1
    need(F(n)**2>=q,'outward square root')
    return F(n)
def n2(row,w):return sum((a*a*b for a,b in zip(row,w)),F(0))
def box(v):
    a,b=map(F,v);need(a<=b,'box ordering');return (a,b)
def widen(v,e):return (v[0]-e,v[1]+e)

# Nonzero independent cone-integration control.  Integrate the z-primitive
# of f(s,z)=aa+bb*s+cc*z as an exact polynomial in s; separately integrate
# the endpoint sums/differences for the x and r derivative formulas.
def padd(a,b):return [sum(q) for q in zip_longest(a,b,fillvalue=F(0))]
def pscale(a,c):return [c*z for z in a]
def pmul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for j,z in enumerate(a):
        for k,v in enumerate(b):out[j+k]+=z*v
    return out
def integral(a,x):return sum((z*x**(j+1)/F(j+1) for j,z in enumerate(a)),F(0))
tx,tr=F(3,5),F(2);aa,bb,cc=F(2),F(3),F(5)
upper,lower=[tr+tx,F(-1)],[tr-tx,F(1)]
primitive=lambda z:padd(pmul([aa,bb],z),pscale(pmul(z,z),cc/2))
toy=integral(padd(primitive(upper),pscale(primitive(lower),-1)),tx)/2
fu=padd([aa,bb],pscale(upper,cc));fl=padd([aa,bb],pscale(lower,cc))
toy_x=integral(padd(fu,fl),tx)/2;toy_r=integral(padd(fu,pscale(fl,-1)),tx)/2
need(toy==(aa+cc*tr)*tx**2/2+bb*tx**3/6,'cone source integral')
need(toy_x==(aa+cc*tr)*tx+bb*tx**2/2 and toy_r==cc*tx**2/2,'cone endpoint derivatives')
need(toy!=2*toy and toy_r!=-toy_r,'missing half/reversed radial control')

pa=argparse.ArgumentParser();pa.add_argument('--repo-root',required=True,type=Path)
pa.add_argument('--smooth-capsule',type=Path)
a=pa.parse_args();root=a.repo_root.resolve()
cap=(a.smooth_capsule or root/'research/HDBLAST_SMOOTH_CANDIDATES_20261026').resolve()
base=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
p=read(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
r=read(base/'compute/response-normal/RESPONSE_REPORT.json','ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
i=read(base/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
need(r['status']=='PASS_BOUNDED_FULL23_RESPONSE_X_ERROR','old response status')
need(i['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER','old initial status')
labels=p['labels'];w=list(map(F,p['gram']))
need(labels==i['labels']==r['labels'] and w==list(map(F,i['gram']))==list(map(F,r['gram'])),'space mismatch')
need(len(w)==23 and min(w)>0 and labels[12]=='w:0' and w[12]==1,'selected component/Gram')
def weight(label):
    if label.startswith('T:'):return F(1,6)
    ch,l=label.split(':');l=int(l);z={0:F(1),2:F(1,12),4:F(19,720)}[l]
    if ch=='L':z*=2
    if ch=='J':z*=3
    if ch in ('u','v'):z*=2*l*(l+2)
    if ch=='H':z*=F(2,3)*l*(l+2)*(l*(l+2)-3)
    return z
need(w==list(map(weight,labels)),'reconstructed Gram')
blocks=[[j for j,l in enumerate(labels) if l.endswith(':'+s)] for s in ('0','2','4','+','-')]
need(sum(map(len,blocks))==23,'block partition')
ops={k:F(r['constant_matrix_H_operator_upper'][k]) for k in ('A','B','D','M')}
for k,upper in ops.items():
    M=[[F(v) for v in row] for row in p['matrices'][k]]
    for u in range(23):
        for v in range(23):
            if not any(u in b and v in b for b in blocks):need(M[u][v]==0,'off-block coupling')
    for b in blocks:
        need(sum(w[u]*M[u][v]**2/w[v] for u in b for v in b)<=upper**2,'operator cap:'+k)

nodes=list(map(F,p['nodes']));du=F(1,2**30);dv=F(1,2048)
need(nodes[0]==0 and nodes[-1]==F(7,4) and min(b-a for a,b in zip(nodes,nodes[1:]))>=2*dv,'grid')
profiles={k:[list(map(F,row)) for row in p['profiles'][k]] for k in ('u','du','v')}
need(all(all(z==0 for z in row) for tab in profiles.values() for row in tab[:2]),'zero initial collar')
need(nodes[1]-max(du,dv)>1,'smoothed zero wedge')
slopes=[[(b[j]-a[j])/(nodes[k+1]-nodes[k]) for j in range(23)] for k,(a,b) in enumerate(zip(profiles['v'],profiles['v'][1:]))]
jumps=[[b[j]-a[j] for j in range(23)] for a,b in zip(slopes,slopes[1:])]
av=F(35,256)*dv*sqrtup(max(n2(v,w) for v in jumps))
U0,U1,_=map(F,i['nominal_H_caps']);M2=F(i['nominal_second_derivative_H_caps'][0])
au=M2*du**2/18;bu=35*M2*du/128
EU,EUr,EV=map(F,i['P1_H_error_caps'])
need([EU,EUr,EV]==[F(r['ledger'][k]) for k in ('epsilon_u_H','epsilon_du_H','epsilon_v_H')],'interpolation ledger')
old=F(r['ledger']['total_X_error_upper'])
free=max((au+F(3,2)*av)/3,bu+av)
feedback=F(3,2)*((ops['A']+ops['B'])*(EUr+bu+av)+ops['D']*(EU+au+F(3,2)*av))+ops['M']*(EU+au+av/4)
E=old+free+feedback
need(E<F(2184*10**25),'tighter candidate C1 error')

# The cone is unreflected.  Its initial endpoints z+-s lie in [r-x,r+x].
x=F(1,1000);rad=F(1499,1000);lo=rad-x;hi=rad+x
need(0<x<rad and lo>0 and hi==F(3,2),'unreflected cone')
need(F(66,64)<lo-du and hi+du<F(7,4),'nominal convolution footprint')
need(0<lo-dv and hi+dv<F(7,4),'velocity convolution footprint')
need(rad>max(F(0),x-F(5,4)) and rad<F(7,4)-x,'Dplus interior')
Vcap=sqrtup(max(n2(v,w) for v in profiles['v']))
Fcap=(ops['A']+ops['B'])*(U1+Vcap)+ops['D']*(U0+x*Vcap)
B1cap=x*x*Fcap/2;B1der=x*Fcap
need(min(U0,U1,Vcap,Fcap,B1cap,B1der)>=0,'positive cap')
need(B1cap*2==x*x*Fcap and B1der==x*Fcap,'cone factors')
need(Fcap>ops['A']*(U1+Vcap),'nonzero B/D controls')

script=cap/'checks/check_smoothed_freewave.py'
pin='622886d396cb28bffea1e37dafa330460a0159beb199b96cf1678a034e7d4a89'
need(hashlib.sha256(script.read_bytes()).hexdigest()==pin,'B0 checker hash')
cmd=[sys.executable,'-I','-S','-B']
if sys.flags.optimize:cmd.append('-O')
cmd += [str(script),'--repo-root',str(root)]
b0=json.loads(subprocess.check_output(cmd,text=True))
need(b0['status']=='PASS' and b0['point']==[str(x),str(rad)] and b0['component']=='w:0','B0 point identity')
need(b0['delta_u']==str(du) and b0['delta_v']==str(dv),'B0 candidate widths')
source_keys={'value':'B0_value_interval','dx':'B0_x_interval','dr':'B0_r_interval'}
base_boxes={key:box(b0[src]) for key,src in source_keys.items()}
candidate={key:widen(v,B1cap if key=='value' else B1der) for key,v in base_boxes.items()}
# Error is zero at r0=1-x.  This fixed-x radial path stays in D+.
r0=1-x;distance=rad-r0
need(distance==F(1,2) and 0<=r0<rad<F(7,4)-x,'sharp zero-wedge path')
actual={key:widen(v,distance*E if key=='value' else E) for key,v in candidate.items()}
need(actual['value'][0]>0,'actual formal scalar mean positive')
need(actual['dr'][1]<0,'actual formal scalar mean radial derivative negative')
need(actual['value'][0]>F(417*10**26),'readable scalar mean lower bound')
need(actual['dr'][1]<-F(5283*10**28),'readable radial derivative upper bound')
print(json.dumps({'status':'PASS','point':[str(x),str(rad)],'component':'w:0','B0_checker_sha256':pin,'velocity_P1_H_nodal_cap':str(Vcap),'local_bulk_feedback_source_H_cap':str(Fcap),'B1_value_H_cap':str(B1cap),'B1_derivative_H_cap':str(B1der),'global_C1_error':str(E),'local_value_error':str(distance*E),'candidate_component_intervals':{k:list(map(str,v)) for k,v in candidate.items()},'actual_formal_component_intervals':{k:list(map(str,v)) for k,v in actual.items()},'strict_readable_actual_w0_lower':'417*10^26','strict_readable_actual_w0_r_upper':'-5283*10^28','interpretation':'registered second scalar mean coefficient and scaled radial derivative; conditional on inherited full23 formal response and C1 certificates','method':'independent unreflected-cone norm bounds plus pinned nominal B0 interval, not tight B1 quadrature','limits':'no gauge-invariant scalar observable, full response sample, C2 error, pointwise Einstein residual, finite-amplitude solution or cosmological evidence'},sort_keys=True))
