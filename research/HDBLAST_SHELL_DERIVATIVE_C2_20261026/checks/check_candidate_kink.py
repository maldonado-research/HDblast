#!/usr/bin/env python3
"""Exact kink witness for the published P1-based 23-mode candidate, not the exact solution."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, importlib.util, json

def check(ok,msg):
    if not ok: raise RuntimeError(msg)
parser=argparse.ArgumentParser()
parser.add_argument('--repo-root',required=True,type=Path)
a=parser.parse_args()
base=a.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
profile=base/'compute/response-normal/RESPONSE_PROFILES.json'
runner=base/'compute/angular_ridge.py'
initial=base/'initial/full23_initial.py'
def verify(path,pin):
    raw=path.read_bytes(); check(hashlib.sha256(raw).hexdigest()==pin,'SHA256 mismatch: '+path.name); return raw
p=json.loads(verify(profile,'00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
verify(runner,'d6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f')
initial_source=verify(initial,'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1').decode('utf-8')
check('mean_u = [zero, zero, -u*u/3' in initial_source and "for name in ('K', 'L', 'N', 'J', 'w')" in initial_source, 'K:0 nominal displacement provenance')
ns=list(map(F,p['nodes'])); k=ns.index(F(3,2)); x=F(1,1000); r=F(1499,1000); q=r-x; h=F(1,10000)
check(ns[k-1]==F(383,256) and ns[k+1]==F(385,256),'unexpected knot neighbors')
check(0<x<r and x+r==ns[k] and 0<r<F(7,4)-x and r>max(F(0),x-F(5,4)) and x+r>1,'domain/characteristic')
check(all(q+j*h not in ns for j in range(-3,4)) and ns[k-1]<q-3*h<q+3*h<ns[k], 'second characteristic also crosses a knot')
check(ns[k-1]<ns[k]-3*h and ns[k]+3*h<ns[k+1], 'stencil crosses other knots')
def jump(key):
    v=p['profiles'][key]
    return [(F(v[k+1][i])-F(v[k][i]))/(ns[k+1]-ns[k])-(F(v[k][i])-F(v[k-1][i]))/(ns[k]-ns[k-1]) for i in range(23)]
d={key:jump(key) for key in ('u','du','v')}
A,B,D=([[F(e) for e in row] for row in p['matrices'][key]] for key in ('A','B','D'))
check(p['labels'][0]=='K:0','component mismatch')
fp=sum(((A[0][i]+B[0][i])*(d['du'][i]+d['v'][i])+D[0][i]*d['u'][i] for i in range(23)), F(0))/2
free=d['v'][0]/2
feedback=x*fp/2
claimed=free+feedback
check(claimed!=0, 'candidate jump vanished')
# Second, independent calculation differentiates the frozen direct B1_x oracle.
# Its x-dependence is piecewise cubic on the stencil, so the four-point
# one-sided derivative formula is exact rather than a finite-difference limit.
spec=importlib.util.spec_from_file_location('authenticated_angular_ridge',runner)
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
obj=mod.from_data(p); M=p['matrices']['M']
f=lambda z:obj.feedback(z,r,M)['dx'][0]
f0=f(x)
right=(-11*f0+18*f(x+h)-9*f(x+2*h)+2*f(x+3*h))/(6*h)
left=(11*f0-18*f(x-h)+9*f(x-2*h)-2*f(x-3*h))/(6*h)
check(right-left==feedback,'direct B1 oracle versus analytic knot jump')
expected=F(122949678051127779469013307270796835622755434515356118845205572291510230097770149650546701,336999333339382997433337688587745383420464305281757156013795128115200)
check(claimed==expected,'regression mismatch')
print(json.dumps({'status':'PASS','component':'K:0','point':[str(x),str(r)],'characteristic':'r+x=3/2','jump_right_minus_left_V_xx':str(claimed),'free_B0_jump':str(free),'feedback_B1_jump':str(feedback),'approximate_total_jump':float(claimed),'direct_pinned_B1_oracle':'four-point exact one-sided derivatives agree','scope':'existing C1 numerical candidate only; not exact field or Einstein residual'},sort_keys=True))
