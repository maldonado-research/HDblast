#!/usr/bin/env python3
"""Independent raw-derivative/binomial enclosure of the common-gauge source.

Unlike the producer's normalized Taylor tables, all jets here store actual
derivatives and every product explicitly uses binomial coefficients. Interval
rounding is outward on a distinct 10^50 grid. No producer module is imported.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb,factorial
from pathlib import Path
import json
import sys

GRID=10**50
ELL=F(1,10**6)
PINS={'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json':'4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53',
      'HIGHER_GAUGE_NORMS_NORMAL.json':'f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d',
      'ACTUAL_GAUGE_SOLUTION_NORMAL.json':'dd16312b2d5490259d39ec784b45c3a401f1378676e8ddbff97d54d85f2371f6',
      'HIGHER_JETS_NORMAL.json':'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
      'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0'}
POLY=list(map(F,['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27']))


def need(ok,msg):
    if not ok:raise ValueError(msg)


def up(x):return F(-(-F(x)*GRID).__floor__(),GRID)


class Iv:
    def __init__(self,a,b=None):
        a,b=F(a),F(a if b is None else b)
        need(a<=b,'interval order');self.a=F((a*GRID).__floor__(),GRID);self.b=up(b)
    @staticmethod
    def cast(x):return x if isinstance(x,Iv) else Iv(x)
    def __add__(self,o):
        o=self.cast(o);return Iv(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return Iv(-self.b,-self.a)
    def __sub__(self,o):return self+-self.cast(o)
    def __rsub__(self,o):return self.cast(o)+-self
    def __mul__(self,o):
        o=self.cast(o);q=[a*b for a in (self.a,self.b) for b in (o.a,o.b)];return Iv(min(q),max(q))
    __rmul__=__mul__
    def reciprocal(self):
        need(not self.a<=0<=self.b,'nonzero denominator');return Iv(1/self.b,1/self.a)
    def mag(self):return max(abs(self.a),abs(self.b))


class Jet:
    """Actual univariate derivatives, never divided by factorials."""
    def __init__(self,d):self.d=[Iv.cast(x) for x in d]
    def cast(self,x):return x if isinstance(x,Jet) else Jet([x]+[0]*(len(self.d)-1))
    def __add__(self,o):
        o=self.cast(o);return Jet([a+b for a,b in zip(self.d,o.d)])
    __radd__=__add__
    def __neg__(self):return Jet([-a for a in self.d])
    def __sub__(self,o):return self+-self.cast(o)
    def __mul__(self,o):
        o=self.cast(o);return Jet([sum((comb(n,k)*self.d[k]*o.d[n-k] for k in range(n+1)),Iv(0)) for n in range(min(len(self.d),len(o.d)))])
    __rmul__=__mul__
    def reciprocal(self):
        r=[self.d[0].reciprocal()]
        for n in range(1,len(self.d)):
            r.append(-r[0]*sum((comb(n,k)*self.d[k]*r[n-k] for k in range(1,n+1)),Iv(0)))
        return Jet(r)
    def __truediv__(self,o):return self*self.cast(o).reciprocal()
    def __rtruediv__(self,o):return self.cast(o)*self.reciprocal()
    def __pow__(self,n):
        if n<0:return self.reciprocal()**(-n)
        r=self.cast(1)
        for _ in range(n):r=r*self
        return r


class Bound:
    """Actual mixed derivatives with explicit Leibniz binomials."""
    def __init__(self,data=None,order=3):
        self.n=order;self.d={(i,j):up((data or {}).get((i,j),0)) for i in range(order+1) for j in range(order+1-i)}
    def cast(self,x):return x if isinstance(x,Bound) else Bound({(0,0):abs(F(x))},self.n)
    def __add__(self,o):
        o=self.cast(o);n=min(self.n,o.n);return Bound({k:self.d[k]+o.d[k] for k in self.d if sum(k)<=n},n)
    __radd__=__add__
    def __mul__(self,o):
        o=self.cast(o);n=min(self.n,o.n)
        return Bound({(i,j):sum((comb(i,a)*comb(j,b)*self.d[a,b]*o.d[i-a,j-b] for a in range(i+1) for b in range(j+1)),F(0))
                     for i in range(n+1) for j in range(n+1-i)},n)
    __rmul__=__mul__
    def derivative(self,axis):
        return Bound({(i,j):self.d[i+(axis==0),j+(axis==1)] for i in range(self.n) for j in range(self.n-i)},self.n-1)
    @staticmethod
    def lift(jet,axis):
        return Bound({(k,0) if axis==0 else (0,k):v.mag() for k,v in enumerate(jet.d)},len(jet.d)-1)


def polyder(p):return [n*p[n] for n in range(1,len(p))]
def polyval(p,x):
    y=x.cast(0)
    for a in reversed(p):y=y*x+a
    return y


def evaluate(inputs,candidate,injection):
    docs={}
    for name,pin in PINS.items():
        raw=(inputs/name).read_bytes();need(sha256(raw).hexdigest()==pin,'input hash:'+name);docs[name]=json.loads(raw)
    c=docs['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'];g=docs['HIGHER_GAUGE_NORMS_NORMAL.json']
    solution=docs['ACTUAL_GAUGE_SOLUTION_NORMAL.json']['time_window_certificates']['3']
    master=docs['HIGHER_JETS_NORMAL.json']
    need(docs['ACTUAL_POTENTIAL_JET_BOUNDS.json']['actual_model']['U_coefficients']==list(map(str,POLY)),
         'actual potential polynomial')
    jdecode=lambda rows:Jet([factorial(n)*Iv(*row['exact']) for n,row in enumerate(rows)])
    static={k:jdecode(v) for k,v in c['static_z_jets'].items()}
    rho,eta,p,a,alpha=[static[k] for k in ('rho','eta','p','a','alpha')]
    ri=1/rho;U=[polyval(POLY,eta)];pol=POLY
    for _ in range(3):pol=polyder(pol);U.append(polyval(pol,eta))
    rec=alpha/rho**2;B=3*alpha/(rho**3*p);k=rho*(F(5,2)*a-U[1]/p)
    H=jdecode(c['temporal_t_jets']['tanh_t']);qt=jdecode(c['temporal_t_jets']['cosh_minus_three_halves'])
    cv=[Iv(1,1/(1-(3*ELL)**2/2))]
    for n in range(3):cv.append(sum((comb(n,j)*H.d[j]*cv[n-j] for j in range(n+1)),Iv(0)))
    cosh=Jet(cv);sech=cosh.reciprocal()
    lift=Bound.lift
    riM,aM,pM,kM=[lift(v,1) for v in (ri,a,p,k)]
    cM,ciM,HM,qtM=[lift(v,0) for v in (cosh,sech,H,qt)]
    wcaps=list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
    W=Bound({(i,j):wcaps[i+j]/ELL**(i+j) for i in range(4) for j in range(4-i)})
    raw=[]
    for channel in range(3):
        data={}
        for key,v in g['unscaled_coordinate_mixed_partial_caps'].items():
            i,j=map(int,key.split(','));cap=F(v[channel])
            if i+j==0:cap=min(cap,F(solution['actual_amplitude_value_by_channel'][channel]['exact_upper']))
            if i+j==1:cap=min(cap,F(solution['actual_each_scaled_first_derivative_by_channel'][channel]['exact_upper'])/ELL)
            data[i,j]=cap
        raw.append(Bound(data))
    f,b,d=raw
    u=lift(rec,1)*qtM*W
    v=lift(B,1)*qtM*(W.derivative(1)+kM*W)
    scalar=v+pM*f
    terms={'yy':[(4*u+2*riM*f.derivative(1),(1,2,4))],
           'y0':[(riM*f.derivative(0)+b.derivative(1),(1,2,4))],
           'ya':[(riM*ciM*f+cM*d.derivative(1),(2,4,8))],
           '00':[(2*(u+aM*f+b.derivative(0)),(1,2,4))],
           '0a':[(cM*d.derivative(0)+ciM*b,(2,4,8))],
           'ab':[(2*(u+aM*f+HM*b),(1,2,4)),(2*d,(4,8,16))]}
    ri0=ri.d[0].mag();a0=a.d[0].mag();H0time=H.d[0].mag();s=ri0*sech.d[0].mag()
    def baseline(tt):
        value=F(0);grad=[F(0)]*3;second=[F(0)]*9
        for amp,ang in tt:
            F0,Ft,Fz,Ftt,Ftz,Fzz=[amp.d[key] for key in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))]
            t0,t1,t2=ang
            value+=F0*t0
            yy=ri0*ri0*Fzz+a0*ri0*Fz
            gr=[ri0*Fz*t0,ri0*Ft*t0,s*F0*t1]
            sec=[yy*t0,(ri0*ri0*Ftz+ri0*a0*Ft)*t0,ri0*ri0*Ftz*t0,ri0*ri0*Ftt*t0,
                 s*(ri0*Fz+a0*F0)*t1,s*ri0*Fz*t1,
                 s*ri0*(Ft+H0time*F0)*t1,s*ri0*Ft*t1,s*s*F0*t2]
            grad=[a+b for a,b in zip(grad,gr)];second=[a+b for a,b in zip(second,sec)]
        return value,max(grad),max(second)
    baselines={k:baseline(v) for k,v in terms.items()}
    h0=max(v[0] for v in baselines.values());dh=max(v[1] for v in baselines.values());ddh=max(v[2] for v in baselines.values())
    ay=static['a_y'].d[0].mag();beta=H0time*ri0
    K0=a0+beta;K1=max(ay+a0*beta,ri0*ri0)
    h1=dh+2*K0*h0;h2=ddh+5*K0*dh+(2*K1+6*K0*K0)*h0
    e0,e1,de2=baseline([(scalar,(1,2,4))]);e2=de2+K0*e1
    p0=p.d[0].mag();py=(U[1]-4*a*p).d[0].mag();phi2=max(py,a0*p0)
    u1,u2,u3=[U[i].d[0].mag() for i in (1,2,3)]
    P=F(3,2)*h1;DP=F(3,2)*h2;DJ=10*(h1*P+h0*DP)
    QE=10*DJ+100*P*P+2*e1*e1+F(2,3)*(u2*e0*e0+2*u1*e0*h0)
    QC=50*h0*P;DQC=50*(h1*P+h0*DP)
    QS=50*(h0*h0*phi2+h0*e2+h0*P*p0)+u3*e0*e0
    fields={'metric':[h0,h1,h2],'scalar':[e0,e1,e2]}
    sources={'QE':QE,'QEred':QE+DQC,'Qphi':QS,'QC':QC,'nabla_QC':DQC}
    for name,values in fields.items():
        caps=list(map(F,candidate['strict_readable_field_covariant_caps'][name]))
        if injection==name+'_cap':caps=[v/1000 for v in caps]
        need(all(x<y for x,y in zip(values,caps)),'INDEPENDENT_'+name.upper()+'_CAP')
    for name,value in sources.items():
        cap=F(candidate['actual_source_bounds'][name]['strict_simple_upper'])
        if injection=='source_cap':cap/=1000
        need(value<cap,'INDEPENDENT_SOURCE_CAP:'+name)
    need(candidate['coverage']['Einstein_component_count']==15 and candidate['coverage']['nabla_QC_components']==25,'full tensor coverage')
    return {'status':'PASS_INDEPENDENT_RAW_DERIVATIVE_COMMON_GAUGE_BOUNDS',
            'method':'Actual derivative interval jets with explicit binomial products, outward 10^50 grid; no producer imports',
            'field_caps_verified':candidate['strict_readable_field_covariant_caps'],
            'source_caps_verified':{k:v['strict_simple_upper'] for k,v in candidate['actual_source_bounds'].items()},
            'input_sha256':PINS,'metric_components':15,'gauge_derivative_components':25,
            'scope':'Continuous actual common-gauge C2 field and C0 quadratic-source caps, conditional on frozen actual master/gauge/branch theorems; no evolution or source derivative claim.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',type=Path,required=True)
    p.add_argument('--candidate',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--inject',choices=['metric_cap','scalar_cap','source_cap'])
    a=p.parse_args();need(not a.output.exists(),'output already exists')
    raw=a.candidate.read_bytes();candidate=json.loads(raw)
    need(candidate['status']=='PASS_ACTUAL_FULL_COMMON_GAUGE_BULK_HESSIAN_ENCLOSURE','candidate status')
    r=evaluate(a.inputs,candidate,a.inject)
    r['candidate_sha256']=sha256(raw).hexdigest();r['checker_sha256']=sha256(Path(__file__).read_bytes()).hexdigest()
    with a.output.open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    print(r['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError) as e:print(str(e),file=sys.stderr);sys.exit(2)
