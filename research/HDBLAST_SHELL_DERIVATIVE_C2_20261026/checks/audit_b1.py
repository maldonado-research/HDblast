#!/usr/bin/env python3
"""Independent exact direct-ridge evaluation for one full23 candidate kink."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json, bisect

def check(x,m):
    if not x:raise RuntimeError(m)
parser=argparse.ArgumentParser()
parser.add_argument('--repo-root',required=True,type=Path)
args=parser.parse_args()
base=args.repo_root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026'
pins={'compute/angular_ridge.py':'d6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f','compute/full23_response.py':'efa6a6c2a9c7db596fe024246684f760b807e727fd2ddac8e8126182ba941a73','initial/full23_initial.py':'2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1'}
for relative, expected_hash in pins.items():
    check(hashlib.sha256((base/relative).read_bytes()).hexdigest()==expected_hash,'provenance source hash: '+relative)
pth=base/'compute/response-normal/RESPONSE_PROFILES.json'
raw=pth.read_bytes()
check(hashlib.sha256(raw).hexdigest()=='00165ede35939b768150a9cb9b04cd0b01f0ec6be6d54c8dbcf1b37b90959fdc','profile hash')
p=json.loads(raw);check(p['schema']=='HDBLAST_FULL23_NATIVE_RESPONSE_PROFILES_V1','profile schema');nodes=list(map(F,p['nodes']));vs={k:[[F(a) for a in row] for row in p['profiles'][k]] for k in ('u','du','v')}
A,B,D=([[F(a) for a in row] for row in p['matrices'][k]] for k in ('A','B','D'))
check(p['labels'][0]=='K:0','label');plus=[A[0][i]+B[0][i] for i in range(23)];minus=[B[0][i]-A[0][i] for i in range(23)];drow=D[0]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def combine(a,b,s=1):return [x+s*y for x,y in zip(a,b)]
def profile(name,t):
    vals=vs[name]; j=bisect.bisect_right(nodes,t)-1
    check(0<=j<len(nodes)-1,'target beyond P1 grid')
    area=[F(0)]*23; moment=[F(0)]*23
    for k in range(j+1):
        l=nodes[k]; end=min(t,nodes[k+1]); width=end-l
        if width<=0:continue
        slopes=[(vals[k+1][i]-vals[k][i])/(nodes[k+1]-l) for i in range(23)]
        for i in range(23):
            area[i]+=vals[k][i]*width+slopes[i]*width**2/2
            moment[i]+=l*vals[k][i]*width+(vals[k][i]+l*slopes[i])*width**2/2+slopes[i]*width**3/3
    slope=[(vals[j+1][i]-vals[j][i])/(nodes[j+1]-nodes[j]) for i in range(23)]
    val=[vals[j][i]+slope[i]*(t-nodes[j]) for i in range(23)]
    twice=[t*area[i]-moment[i] for i in range(23)]
    return val,area,twice

def coefficient(name,t,order):
    u=profile('u',t);d=profile('du',t);v=profile('v',t)
    if name=='F':return (dot(plus,combine(d[order],v[order]))+dot(drow,combine(u[order],v[order+1])))/2
    return (dot(minus,combine(d[order],v[order],-1))+dot(drow,combine(u[order],v[order+1],-1)))/2

def b1x(x,r):
    a,q=r+x,r-x
    check(0<x<r and q>=0,'only unreflected branch supported')
    return (x*coefficient('F',a,0)+(coefficient('G',a,1)-coefficient('G',q,1))/2+(coefficient('F',a,1)-coefficient('F',q,1))/2+x*coefficient('G',q,0))/2

x=F(1,1000);r=F(1499,1000);h=F(1,10000)
k=nodes.index(F(3,2)); q=r-x
check(nodes[k-1]==F(383,256) and nodes[k+1]==F(385,256),'unexpected adjacent P1 knots')
check(0<x<r and x+r==nodes[k] and r<F(7,4)-x and r>max(F(0),x-F(5,4)) and x+r>1,'witness not in retained interior')
check(nodes[k-1]<q-3*h<q+3*h<nodes[k], 'incoming characteristic crosses an unintended knot')
check(nodes[k-1]<nodes[k]-3*h and nodes[k]+3*h<nodes[k+1], 'outgoing stencil crosses another knot')
f0=b1x(x,r)
right=(-11*f0+18*b1x(x+h,r)-9*b1x(x+2*h,r)+2*b1x(x+3*h,r))/(6*h)
left=(11*f0-18*b1x(x-h,r)+9*b1x(x-2*h,r)-2*b1x(x-3*h,r))/(6*h)
def jump(name):return [(vs[name][k+1][i]-vs[name][k][i])/(nodes[k+1]-nodes[k])-(vs[name][k][i]-vs[name][k-1][i])/(nodes[k]-nodes[k-1]) for i in range(23)]
du,dv,uu=jump('du'),jump('v'),jump('u')
dF=(dot(plus,combine(du,dv))+dot(drow,uu))/2
feedback=x*dF/2;free=dv[0]/2;total=free+feedback
check(right-left==feedback,'independent direct B1 jump')
check(total!=0,'zero total')
expected=F(122949678051127779469013307270796835622755434515356118845205572291510230097770149650546701,336999333339382997433337688587745383420464305281757156013795128115200)
check(total==expected,'wrong full jump')
print(json.dumps({'status':'PASS','point':[str(x),str(r)],'a':str(x+r),'q':str(r-x),'free_jump':str(free),'feedback_jump':str(feedback),'total_jump':str(total),'method':'independent P1 segment integrals of p,m analytic formulas; no project code imported','scope':'kink of one nominal candidate component, not exact field'},sort_keys=True))
