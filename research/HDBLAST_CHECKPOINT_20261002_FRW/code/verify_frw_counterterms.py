"""Exact rational Taylor-jet checks, not a universal symbolic proof or mode evolution."""
from fractions import Fraction as F
from math import isqrt
import json, hashlib
from pathlib import Path
N=6
class Jet:
    def __init__(self,x=0):
        if isinstance(x,Jet):self.c=x.c[:]
        elif isinstance(x,(list,tuple)):self.c=[F(v) for v in (list(x)+[0]*(N+1))[:N+1]]
        else:self.c=[F(x)]+[F(0)]*N
    @staticmethod
    def of(x):return x if isinstance(x,Jet) else Jet(x)
    def __add__(self,b):
        b=Jet.of(b);return Jet([x+y for x,y in zip(self.c,b.c)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.c])
    def __sub__(self,b):return self+(-Jet.of(b))
    def __rsub__(self,b):return Jet.of(b)+(-self)
    def __mul__(self,b):
        b=Jet.of(b);return Jet([sum((self.c[i]*b.c[n-i] for i in range(n+1)),F(0)) for n in range(N+1)])
    __rmul__=__mul__
    def inv(self):
        out=[1/self.c[0]]+[F(0)]*N
        for n in range(1,N+1):out[n]=-sum((self.c[i]*out[n-i] for i in range(1,n+1)),F(0))/self.c[0]
        return Jet(out)
    def __truediv__(self,b):return self*Jet.of(b).inv()
    def __rtruediv__(self,b):return Jet.of(b)*self.inv()
    def __pow__(self,p):
        assert isinstance(p,int)
        if p<0:return self.inv()**(-p)
        out,base=Jet(1),self
        while p:
            if p&1:out=out*base
            base=base*base;p//=2
        return out
    def d(self):return Jet([(i+1)*self.c[i+1] for i in range(N)]+[0])
    def sqrt(self):
        z=self.c[0];assert z>0
        root=F(isqrt(z.numerator),isqrt(z.denominator));assert root*root==z
        out=[root]+[F(0)]*N
        for n in range(1,N+1):out[n]=(self.c[n]-sum((out[i]*out[n-i] for i in range(1,n)),F(0)))/(2*root)
        return Jet(out)
def counterterms(a,h,k,r):
    w=(k*k+a*a*r).sqrt();L=a.d()/a;D=a*a*h
    B=-k*k/F(3)-a*a*r;C=1-B/w**2
    u=(D-a.d().d()/a)/(2*w)-w.d().d()/(4*w**2)+3*w.d()**2/(8*w**3)
    v=(-u**2/(2*w)-u.d().d()/(4*w**2)+w.d().d()*u/(4*w**3)
       +3*w.d()*u.d()/(4*w**3)-3*w.d()**2*u/(4*w**4))
    shared2=L*w.d()/w**2+w.d()**2/(4*w**3)
    shared4=(L*u.d()/w**2-2*L*w.d()*u/w**3
       +w.d()*u.d()/(2*w**3)-3*w.d()**2*u/(4*w**4))
    return dict(w=w,L=L,u=u,v=v,rho0=w/(2*a**4),
       rho2=((D+L**2)/w+shared2)/(4*a**4),
       rho4=(u**2/w-(D+L**2)*u/w**2+shared4)/(4*a**4),
       p0=k*k/(6*a**4*w),p2=(C*u+(L**2-D)/w+shared2)/(4*a**4),
       p4=(C*v+B*u**2/w**3-(L**2-D)*u/w**2+shared4)/(4*a**4),
       Q0=1/(2*a**2*w),Q2=-u/(2*a**2*w**2))
def zero(j,label):assert j.c[0]==0,(label,str(j.c[0]))
def main():
    bad_p=False;bad_Q=False
    for i in range(24):
        a0=F(4+i%5,4);k=F(1+i%3,2);w0=F(5+i%4,2);r=(w0*w0-k*k)/(a0*a0);assert r>0
        a=Jet([a0,F(i+2,19),F(-i-1,37),F(i+3,83),F(-i-2,173),F(i+1,349),F(-i-1,701)])
        h0=-r if i%3==0 else (-r/2 if i%3==1 else r/3)
        h1=F(0) if i%3==0 else F(i+1,23)
        h=Jet([h0,h1,F(i+2,43),F(-i-1,89),F(i+1,181),F(-i-1,367),F(i+2,739)])
        z=counterterms(a,h,k,r);L=z["L"]
        R1=z["rho0"].d()+3*L*(z["rho0"]+z["p0"])
        R3=z["rho2"].d()+3*L*(z["rho2"]+z["p2"])-h.d()*z["Q0"]/2
        R5=z["rho4"].d()+3*L*(z["rho4"]+z["p4"])-h.d()*z["Q2"]/2
        for grade,R in [(1,R1),(3,R3),(5,R5)]:zero(R,("Ward",i,grade))
        bad_p|=(R5-3*L*z["p4"]).c[0]!=0;bad_Q|=(R5+h.d()*z["Q2"]/2).c[0]!=0
    assert bad_p and bad_Q
    for i in range(12):
        a=Jet(1);k=F(1+i%3,2);w0=F(5+i%4,2);r=w0*w0-k*k
        h=Jet([-r if i%2==0 else r/3,F(i+1,17),F(i+2,31),F(-i-1,67)])
        z=counterterms(a,h,k,r);w=z["w"]
        er=w/2+h/(4*w)-h**2/(16*w**3)
        ep=k*k/(6*w)-h*k*k/(12*w**3)+h**2*k*k/(16*w**5)-h.d().d()*(2*k*k/F(3)+r)/(16*w**5)
        eq=1/(2*w)-h/(4*w**3)
        zero(z["rho0"]+z["rho2"]+z["rho4"]-er,("flat rho",i))
        zero(z["p0"]+z["p2"]+z["p4"]-ep,("flat p",i));zero(z["Q0"]+z["Q2"]-eq,("flat Q",i))
    out=dict(status="PASS",method="exact rational Taylor jets; selected local jets, not universal symbolic proof",
       general_jet_cases=24,graded_Ward_assertions=72,flat_matching_assertions=36,
       deliberate_pressure_omission_detected=bad_p,deliberate_current_omission_detected=bad_Q,
       field_evolution_executed=False,absolute_FRW_source_claimed=False,jet_order=N,
       source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path("outputs").mkdir(exist_ok=True)
    Path("outputs/counterterm_checks.json").write_text(json.dumps(out,indent=2)+"\n")
    print("HDBLAST_FRW_COUNTERTERMS_JSON="+json.dumps(out,separators=(",",":")))
if __name__=="__main__":main()

