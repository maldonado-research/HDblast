#!/usr/bin/env python3
"""Independent rational reconstruction from retained raw boxes; no oracle/imports."""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import argparse
import hashlib
import json

import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from portable_io import arguments, authenticate, finish
HERE=Path(__file__).resolve().parent
ROOT=None
NUM=None
ell=F(1,10**6);T=F(5,4);delta=T-1
def need(ok,msg):
    if not ok:raise ValueError(msg)
@dataclass(frozen=True)
class Box:
    lo:F
    hi:F
    def __post_init__(self):need(self.lo<=self.hi,'invalid box')
    def __add__(self,b):
        b=point(b);return Box(self.lo+b.lo,self.hi+b.hi)
    __radd__=__add__
    def __neg__(self):return Box(-self.hi,-self.lo)
    def __sub__(self,b):return self+-point(b)
    def __rsub__(self,b):return point(b)+-self
    def __mul__(self,b):
        b=point(b);v=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return Box(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=point(b);need(b.hi<0 or b.lo>0,'division through zero')
        return self*Box(1/b.hi,1/b.lo)
    def __rtruediv__(self,b):return point(b)/self
    def __pow__(self,n):
        need(n>=0 and int(n)==n,'power')
        if n==0:return point(1)
        if n%2:return self*(self**(n-1))
        return Box(F(0) if self.lo<=0<=self.hi else min(self.lo**n,self.hi**n),max(self.lo**n,self.hi**n))
    def mag(self):return max(abs(self.lo),abs(self.hi))
def point(v):return v if isinstance(v,Box) else Box(F(v),F(v))
def readbox(row):return Box(*map(F,row['exact']))
def symmetric(v):return Box(-v,v)
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rootup(v):
    scale=1<<80
    scaled=v*scale*scale
    n=-((-scaled.numerator)//scaled.denominator)
    return F(isqrt(n)+1,scale)
def gramcap(v):return rootup(sum(F(w)*x*x for w,x in zip((1,2,1,3,1),v)))

def evaluate(report_path, initial_path, conditioned_path, repo_root, mutant=None):
    global ROOT
    ROOT=repo_root/"research"
    report=load(report_path)
    original_path=initial_path
    raw=load(original_path);conditioned=load(conditioned_path)
    source=ROOT/'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    shell_path=source/'numeric/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'
    bg_path=source/'numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json'
    bounds_path=ROOT/'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023/golden/BOUNDS.json'
    pins={}
    py=readbox(load(bg_path)['y_jets']['p'][1])

    def data(row):
        p={k:readbox(v) for k,v in row['profiles'].items()}
        rr,aa,pp,uu,du,uy,ff,fy,xi,xiy,qq=[p[k] for k in ('rho','a','p','u','ut','uy','f','fy','xi','xiy','QM')]
        initial=[ff,point(0),-(uu*uu)/3,ff,-pp*xi]
        radial=[-ell*rr*fy,point(0),F(2,3)*ell*rr*uu*uy,-ell*rr*fy,ell*rr*(py*xi+pp*xiy)]
        velocity=[ell*(F(2,3)*rr*qq/aa+F(32,3)*uu*du),
                  ell*rr*(-fy+aa*ff-F(13,3)*uu*uy-aa*uu*uu/3),
                  ell*(-F(2,3)*rr*qq/aa+F(4,3)*uu*du),-F(4,3)*ell*uu*du,point(0)]
        return (initial,radial,velocity)
    need(len(raw['profile_cells'])==994,'expected frozen coverage')
    prior=F(1);count=0
    for original,adapted in zip(raw['profile_cells'],conditioned['profile_cells']):
        lo,hi=map(F,original['r']);need(lo==prior and hi>lo,'partition gap/overlap');prior=hi
        need(original['r']==adapted['r'],'conditioned partition')
        for key,vec in zip(('U0','U0_r','V0'),data(original)):
            for slot,val in zip(('K0','L0','N0','J0','w0'),vec):
                actual=readbox(adapted['conditioned'][key][slot])
                need(actual.lo<=val.lo<=val.hi<=actual.hi,'conditioned enclosure inclusion')
                count+=1
    need(prior==2,'packet endpoint')
    used=[r for r in raw['profile_cells'] if F(r['r'][1])<=T]
    need(F(used[-1]['r'][1])==T,'target initial footprint')
    caps=[[F(0)]*5 for _ in range(3)]
    moments=[point(0),point(0)];absolute=[F(0),F(0)]
    for row in used:
        lo,hi=map(F,row['r']);vecs=data(row)
        for cs,vs in zip(caps,vecs):
            for j,v in enumerate(vs):cs[j]=max(cs[j],v.mag())
        for n in range(2):
            weight=((T-lo)**(n+1)-(T-hi)**(n+1))/(n+1)
            moments[n]+=weight*vecs[0][4]
            absolute[n]+=weight*vecs[0][4].mag()
    d0,dr,dv=map(gramcap,caps)
    free=max((d0+T*dv)/3,dr+dv)
    q=T*(F(2,10000)+delta/F(10**6))+delta/2000
    bounds=load(bounds_path)['exact']
    source_cap=T*F(bounds['F_H_upper'])+F(bounds['full_G_H_upper'])
    err=(q*free+source_cap)/(1-q)

    first=report['first']['raw_actual_shell_2jet']
    need(first['x']==['1','5/4'] and first['r']==['0','0'],'first source domain')
    m={k:readbox(v) for k,v in first['m'].items()}
    s={k:readbox(v) for k,v in first['s'].items()}
    rho=readbox(first['background']['rho'])
    X=(m['D']/2).mag();Z=(m['L']/2).mag();d=(m['M']/2).mag();shift=m['E'].mag()
    clock=rho.hi*ell*delta*X
    Hcap=ell*T/rho.lo;Hpcap=1/rho.lo**2
    Ecap=Z+Hcap*clock
    Epcap=(m['Lt']/(2*rho)).mag()+Hpcap*clock+Hcap*X
    Dpcap=(m['Mt']/(2*rho)).mag()+shift/rho.lo+clock/rho.lo**2
    psi=s['s'].mag();psip=(s['st']/rho).mag()
    correct=clock*psip/3+(Ecap+F(8,3)*d)*psi
    if mutant=='omit_measure':correct=clock*psip/3
    shell=load(shell_path)['boxes']
    sigma,sp,spp,Vp,Vpp=map(lambda k:readbox(shell[k]),('sigma_b','sigma_prime_b','sigma_second_b','U_prime_b','U_second_b'))
    Wp=2*Vp-sp*spp/2+F(2,3)*sigma*sp
    Wpp=2*Vpp-(spp*spp+4*sp)/2+F(2,3)*(sp*sp+sigma*spp)
    remainder=(Wpp.mag()*psi**2/12+psip**2/6+F(4,3)*psi**2/rho.lo**2
        +F(5,3)*Ecap**2/rho.lo**2+Epcap**2+F(16,3)*Epcap*Dpcap+F(8,3)*Dpcap**2)
    target=next(r for r in raw['targets'] if F(r['r'])==T)
    target_scalar=data(target)[0][4]
    if mutant=='halve_scalar':target_scalar=target_scalar/2
    M2=target_scalar+symmetric(delta*err+correct)
    history_errors=[delta**2*err/2+delta*correct,delta**3*err/6+delta**2*correct/2]
    if mutant=='omit_history_error':history_errors=[F(0),F(0)]
    green=[]
    for n in range(2):
        e=history_errors[n]
        green.append(moments[n]+symmetric(e+F(1,10**11)*(absolute[n]+e)))
    volume=(rho*ell)**2*(Wp*green[1]/2+symmetric(remainder)*delta**2/2)
    theta=rho*ell*(Wp*green[0]/2+symmetric(remainder)*delta*(1+F(1,10**11)))
    recomputed={'q':q,'free_X_cap':free,'full_resolvent_minus_free_X_cap':err,
        'proper_M2_target':M2,'proper_v2':volume,'proper_Theta2':theta}
    matches={}
    for key,value in recomputed.items():
        wanted=readbox(report[key]) if isinstance(value,Box) else F(report[key]['exact'])
        need(value==wanted,'independent reconstructed mismatch: '+key)
        matches[key]=True
    need(M2.lo>0 and volume.hi<0 and theta.hi<0,'strict signs')
    for name,cap in [('M2_correction',correct),('Q1',remainder),('clock',clock),('E',Ecap),('E_s',Epcap),('D_s',Dpcap)]:
        need(cap==F(report['first'][name]['exact']),'first bound mismatch: '+name)
        matches[name]=True
    for p,digest in pins.items():need(sha(Path(p))==digest,'input changed during independent reconstruction')
    def encode(x):
        if isinstance(x,Box):return {'exact':[str(x.lo),str(x.hi)],'display_only':[float(x.lo),float(x.hi)]}
        return {'exact':str(x),'display_only':float(x)}
    return {'status':'PASS_INDEPENDENT_RATIONAL_RESPONSE_RECONSTRUCTION',
        'conditioned_component_inclusions':count,'exact_report_matches':matches,
        'independent_results':{k:encode(v) for k,v in recomputed.items()},
        'scope':'No scientific producer import or oracle call; retained raw input enclosures and analytic theorems remain premises.'}

if __name__=='__main__':
    args=arguments(mutants=['omit_measure','halve_scalar','omit_history_error'])
    authenticate(args)
    try:
        result=evaluate(args.response_report,args.initial_report,args.conditioned_report,args.repo_root,args.mutant)
        finish(args,result)
    except ValueError as error:
        import sys
        print(json.dumps({'status':'REJECT','stage':'SCIENCE','reason':str(error)}))
        raise SystemExit(2)
