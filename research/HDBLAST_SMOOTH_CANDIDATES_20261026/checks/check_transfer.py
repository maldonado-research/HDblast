#!/usr/bin/env python3
"""Conditional C1 transfer from frozen full23 P1 candidate to a C4 smoothed one.
Authenticates prior reports, reads no scientific project code, no PDE replay.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, math

def need(x,m):
    if not x:raise RuntimeError(m)
def read(path,pin):
    b=path.read_bytes(); need(hashlib.sha256(b).hexdigest()==pin,'INPUT_HASH '+path.name)
    return json.loads(b)
def ceilroot(q):
    need(q>=0,'negative norm2')
    n=math.isqrt(q.numerator//q.denominator)
    return F(n+(n*n*q.denominator<q.numerator))
def vadd(a,b):return [x+y for x,y in zip(a,b)]
def vsub(a,b):return [x-y for x,y in zip(a,b)]
def vscale(a,k):return [x*k for x in a]
def norm2(a,gram):return sum((g*x*x for g,x in zip(gram,a)),F(0))

ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);args=ap.parse_args()
base=args.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
p=read(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc')
r=read(base/'compute/response-normal/RESPONSE_REPORT.json','ce1a07ea7409253e679e84090ae9ffbcfa9a53dbfe8c620600f2b022ae3431fd')
i=read(base/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad')
need(i['status']=='PASS_CONTINUOUS_FULL23_INITIAL_CANDIDATE_AND_ERROR_LEDGER','initial report status')
need(r['status']=='PASS_BOUNDED_FULL23_RESPONSE_X_ERROR','response report status')
nodes=list(map(F,p['nodes']));delta=F(1,2048);gram=list(map(F,p['gram']))
need(len(gram)==23 and min(gram)>0,'positive Gram')
need(len(nodes)>2 and min(b-a for a,b in zip(nodes,nodes[1:]))>=2*delta,'overlapping smoothing bands')
need(nodes[0]==0 and nodes[-1]==F(7,4),'registered footprint')
need(all(nodes[k]>=1+delta for k in range(1,len(nodes)-1)), 'zero wedge margin')
profiles={k:[[F(v) for v in row] for row in p['profiles'][k]] for k in ('u','du','v')}
need(all(len(row)==23 for tab in profiles.values() for row in tab),'shape')
need(all(not any(row) for tab in profiles.values() for row in tab[:2]), 'initial collar')
need(gram==list(map(F,i['gram'])) and gram==list(map(F,r['gram'])),'Gram mismatch')
need(p['labels']==i['labels']==r['labels'],'label mismatch')
slopes={k:[vscale(vsub(tab[j+1],tab[j]),1/(nodes[j+1]-nodes[j])) for j in range(len(nodes)-1)] for k,tab in profiles.items()}
jumps={k:[vsub(v[j],v[j-1]) for j in range(1,len(v))] for k,v in slopes.items()}
maxjump={k:ceilroot(max(norm2(v,gram) for v in tab)) for k,tab in jumps.items()}
# max over P1 cells of ||slope(uP1)-duP1|| occurs at endpoints by convexity.
mismatch=ceilroot(max(norm2(vsub(s,profiles['du'][j+side]),gram) for j,s in enumerate(slopes['u']) for side in (0,1)))
alpha_u=F(35,256)*delta*maxjump['u'];alpha_v=F(35,256)*delta*maxjump['v']
# Within each disjoint smoothing band the CDF differs from its step by <=1/2.
alpha_du=mismatch+maxjump['u']/2
E=F(r['ledger']['total_X_error_upper']);eU,eUr,eV=map(F,i['P1_H_error_caps'])
need([eU,eUr,eV]==[F(r['ledger'][k]) for k in ('epsilon_u_H','epsilon_du_H','epsilon_v_H')],'P1 ledger mismatch')
nA,nB,nD,nM=(F(r['constant_matrix_H_operator_upper'][k]) for k in ('A','B','D','M'))
need(all(v>=0 for v in (E,eU,eUr,eV,nA,nB,nD,nM,alpha_u,alpha_v,alpha_du)), 'negative bound')
free=max((eU+alpha_u+F(3,2)*alpha_v)/3, eUr+alpha_du+alpha_v)
feedback=F(3,2)*((nA+nB)*(alpha_du+alpha_v)+nD*(alpha_u+F(3,2)*alpha_v))+nM*(alpha_u+alpha_v/4)
added=free+feedback; total=E+added
need(total < F(9711*10**28), 'C1 transfer exceeds documented 9.711e31')
# Monotone/outward integers for a portable transparent bound.
ceil=lambda q:(q.numerator+q.denominator-1)//q.denominator
print(json.dumps({'status':'PASS','delta':str(delta),'source_sha256':'00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc','nodes':len(nodes),'interior_hinges':len(nodes)-2,'kernel':'k(z)=35/(32δ)(1−(z/δ)^2)^3 on |z|<=δ, zero outside','derivative_mismatch_bound':str(alpha_du),'u_smoothing_value_error':str(alpha_u),'v_smoothing_value_error':str(alpha_v),'free_transfer_upper':str(free),'feedback_transfer_upper':str(feedback),'old_C1_error':str(E),'new_C1_error_upper':str(total),'new_C1_integer_upper':str(ceil(total)),'added_approx':float(added),'old_approx':float(E),'new_approx':float(total),'scope':'conditional C1 transfer for stated candidate; not a numerical C2 residual or independent proof of inherited bounds'},sort_keys=True))
