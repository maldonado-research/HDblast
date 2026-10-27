#!/usr/bin/env python3
"""Rigorous one-point quadrature for mollified nominal u0; no full C2 error bound."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json

def need(x,s):
    if not x:raise RuntimeError(s)
def read(p,h):
    b=p.read_bytes();need(hashlib.sha256(b).hexdigest()==h,'INPUT_HASH '+str(p));return b
ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);args=ap.parse_args();root=args.repo_root.resolve()
base=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
init=base/'initial/full23_initial.py'
read(init,'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1')
report=json.loads(read(base/'initial/results/INITIAL_REPORT.json','bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'))
labels=report['labels'];gram=list(map(F,report['gram']))
need(labels[12]=='w:0' and gram[12]==1,'H norm component bound')
U1=F(report['nominal_H_caps'][1]);U2=F(report['nominal_second_derivative_H_caps'][0]);need(U1>0 and U2>0,'nominal derivative caps')
spec=importlib.util.spec_from_file_location('verified_full23_initial',init);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
context=module.Context(root) # authenticates all seven upstream roots before project imports
N=32;delta=F(1,2**30);r=F(3,2)
need(F(66,64)<r-delta<r+delta<F(7,4),'quadrature touches cutoff or domain edge')
def primitive(z):return F(35,32)*sum(F(c,2*j+1)*z**(2*j+1) for j,c in enumerate((1,-3,3,-1)))
weights=[];low=[F(0),F(0)];high=[F(0),F(0)]
for k in range(N):
    left=F(-1)+F(2*k,N);right=left+F(2,N);mid=(left+right)/2
    weight=primitive(right)-primitive(left);need(weight>0,'positive exact kernel panel mass');weights.append(weight)
    values,_=context.profiles(r-delta*mid)
    for j in range(2):
        x=values[j][12];low[j]+=weight*x.lo;high[j]+=weight*x.hi
need(sum(weights)==1,'kernel normalization')
# On each panel, |r-delta*z - (r-delta*mid)| <= delta/N.
# The Gram-weighted derivative bounds imply the component bounds since H_12=1.
errors=[U1*delta/N,U2*delta/N]
for j in range(2):low[j]-=errors[j];high[j]+=errors[j];need(low[j]<=high[j],'interval order')
need(high[0]<0,'w:0 candidate value should remain negative in this fixture')
context.reauthenticate()
print(json.dumps({'status':'PASS','profile':'registered nominal displacement w:0 and its r derivative, mollified by normalized C2 kernel','r':str(r),'delta':str(delta),'panels':N,'value_interval':[str(low[0]),str(high[0])],'derivative_interval':[str(low[1]),str(high[1])],'value_quadrature_error_upper':str(errors[0]),'derivative_quadrature_error_upper':str(errors[1]),'inherited_nominal_derivative_caps':'u0_r and u0_rr in H norm; coefficient w:0 has unit Gram weight','comparison':'single point, nominal displacement; no true-field interval, no full-response residual, no nonlinear existence'},sort_keys=True))
