#!/usr/bin/env python3
"""Rigorous one-point quadrature for mollified nominal u0; no full C2 error bound."""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,importlib.util,json,math

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
# These are the frozen candidate's former kink coordinates, with both
# characteristic endpoints strictly in one analytic nominal branch.
N=32;delta_u=F(1,2**30);delta_v=F(1,2048);x=F(1,1000);r=F(1499,1000)
endpoints=[r+x,r-x]
need(all(F(66,64)<z-delta_u<z+delta_u<F(7,4) for z in endpoints),'quadrature touches cutoff or domain edge')
def primitive(z):return F(35,32)*sum(F(c,2*j+1)*z**(2*j+1) for j,c in enumerate((1,-3,3,-1)))
weights=[]
for k in range(N):
    left=F(-1)+F(2*k,N);right=left+F(2,N)
    weight=primitive(right)-primitive(left);need(weight>0,'positive exact kernel panel mass');weights.append(weight)
need(sum(weights)==1,'kernel normalization')

def quad(endpoint):
    low=[F(0),F(0)];high=[F(0),F(0)]
    for k,w in enumerate(weights):
        mid=F(-1)+F(2*k+1,N)
        values,_=context.profiles(endpoint-delta_u*mid)
        for j in range(2):
            a=values[j][12];low[j]+=w*a.lo;high[j]+=w*a.hi
    errors=[U1*delta_u/N,U2*delta_u/N]
    return [(low[j]-errors[j],high[j]+errors[j]) for j in range(2)]

qa,qb=map(quad,endpoints)
# A distinct, exact local P1 convolution for w:0 velocity and primitive.  An
# affine continuation is used past the last grid point; the zero collar makes
# the positive-side hinge representation valid without negative-side terms.
prof=json.loads(read(base/'compute/response-normal/RESPONSE_PROFILES.json','00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc'))
need(prof['labels'][12]=='w:0' and F(prof['gram'][12])==1,'velocity label/Gram')
nodes=list(map(F,prof['nodes']));v=[F(a[12]) for a in prof['profiles']['v']]
need(v[0]==v[1]==0 and nodes[0]==0,'velocity zero collar')
slopes=[(v[j+1]-v[j])/(nodes[j+1]-nodes[j]) for j in range(len(nodes)-1)]
need(slopes[0]==0,'velocity first cell')
jumps=[(nodes[j],slopes[j]-slopes[j-1]) for j in range(1,len(slopes))]
def moment(k):
    if k%2:return F(0)
    return F(35,32)*sum(F(2*c,2*j+k+1) for j,c in enumerate((1,-3,3,-1)))
def hinge(z,n):
    if z<=-delta_v:return F(0)
    if z>=delta_v:
        return sum(F(math.comb(n,k),math.factorial(n))*z**(n-k)*(-delta_v)**k*moment(k) for k in range(n+1))
    s=z/delta_v; out=F(0)
    for j,c in enumerate((1,-3,3,-1)):
        for k in range(n+1):
            m=2*j+k
            out+=c*math.comb(n,k)*s**(n-k)*(-1)**k*(s**(m+1)-(-1)**(m+1))/F(m+1)
    return F(35,32)*delta_v**n*out/F(math.factorial(n))
def vs(z):return sum((jump*hinge(z-a,1) for a,jump in jumps if z>a-delta_v),F(0))
def vsprimitive(z):return sum((jump*hinge(z-a,2) for a,jump in jumps if z>a-delta_v),F(0))
va,vb=map(vs,endpoints);ia,ib=map(vsprimitive,endpoints)
# Independent nonzero control rejects omitting velocity or a characteristic end.
need(va!=0 and vb!=0 and ia!=ib,'nonzero velocity controls')
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[1],a[1]-b[0])
def half(a):return (a[0]/2,a[1]/2)
b0=half(add(add(qa[0],qb[0]),(ia-ib,ia-ib)))
b0x=half(add(sub(qa[1],qb[1]),(va+vb,va+vb)))
b0r=half(add(add(qa[1],qb[1]),(va-vb,va-vb)))
need(all(a<=b for a,b in (b0,b0x,b0r)),'interval ordering')
context.reauthenticate()
print(json.dumps({'status':'PASS','point':[str(x),str(r)],'component':'w:0','candidate':'free-wave B0 only; mollified nominal displacement plus mollified P1 velocity','delta_u':str(delta_u),'delta_v':str(delta_v),'panels_per_endpoint':N,'B0_value_interval':list(map(str,b0)),'B0_x_interval':list(map(str,b0x)),'B0_r_interval':list(map(str,b0r)),'source':'nominal initial profiles; exact kernel masses and H-derived midpoint error','scope':'not B1, total candidate, true field, Einstein residual or nonlinear solution'},sort_keys=True))
