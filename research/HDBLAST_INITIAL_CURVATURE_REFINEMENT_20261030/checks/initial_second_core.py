"""Order-five extension of authenticated paired initial-profile algebra.

Same analytic data and three telescoping channels as the retained producer.
The Taylor remainder of each background derivative is enclosed separately;
an unknown remainder function is never differentiated.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

PRODUCER_PIN = '2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1'
INITIAL_PIN = 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'
PROFILES_PIN = '555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def pin(path, digest):
    raw=Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest()==digest,'INPUT_HASH:'+str(path))
    return raw


class HigherInitial:
    def __init__(self, root):
        self.root=Path(root).resolve()
        base=self.root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial'
        self.producer=base/'full23_initial.py'
        pin(self.producer,PRODUCER_PIN)
        spec=importlib.util.spec_from_file_location('higher_initial_authenticated',self.producer)
        self.mod=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.mod)
        self.ctx=self.mod.Context(self.root)
        self.I,self.J,self.P=self.ctx.I,self.ctx.J,self.ctx.P
        self.ell=self.ctx.ell
        self.order=5
        need(self.ctx.data['bump']['maximum_derivative_order']>=5 and
             len(self.ctx.data['bump']['polynomial_coefficients_ascending'])>=6,
             'authenticated fifth bump derivative exists')

    def constant(self,z):
        return self.ctx.constant(z,self.order)

    def background(self,r):
        c=self.ctx;z=self.ell*self.J([r,1/self.ell]+[0]*4)
        out={}
        for key in ('rho','eta','p','a','alpha'):
            nominal=c.polynomial(c.cs[key],z)
            coefficients=[self.I(*v) for v in c.frozen['shell_Taylor_coefficients'][key]]
            fifth=self.I(*c.frozen['uniform_fifth_normalized_ODE_derivative'][key])
            branch=[];remainder=[]
            for n in range(6):
                branch.append(sum((math.comb(j,n)*(coefficients[j]-c.cs[key][j])*z.c[0]**(j-n)
                    for j in range(n,5)),self.I(0)))
                remainder.append(math.comb(5,n)*fifth*z.c[0]**(5-n))
            differences=[self.J(branch),self.constant(0),self.J(remainder)]
            out[key]=self.P(nominal+differences[0]+differences[2],nominal,differences)
        return out

    def kernel(self,t):
        c=self.ctx;b=1+t/c.power;rational=c.square40(1/b)
        t0,b0=t.val(),1+t.val()/c.power;R=rational.val()
        need(t0.lo>=0,'nonnegative exponential argument')
        differences=[];rising=F(1)
        for k in range(6):
            if k:rising*=1+F(k-1,c.power)
            cs=[F(math.comb(k,j),c.power**j) for j in range(k+1)]
            cs[0]-=rising
            pdiff=c.init.poly(cs,t0)
            D0=self.I(-(R*t0**2/(2*c.power)).hi,0)
            differences.append((-1)**k*(D0+R*pdiff/b0**k))
        dt=self.J([self.I(0)]+t.c[1:]);difference=self.constant(0)
        actual=self.constant(0);power_dt=self.constant(1);E=c.init.expbox(t0)
        for k in range(6):
            difference+=power_dt*differences[k]/math.factorial(k)
            actual+=power_dt*((-1)**k*E)/math.factorial(k)
            power_dt*=dt
        return self.P(actual,rational,[self.constant(0),difference,self.constant(0)])

    def bump(self,r,tail=False):
        c=self.ctx
        polys=[list(map(F,row)) for row in c.data['bump']['polynomial_coefficients_ascending'][:6]]
        if tail:
            edge=F(1,64);smin=1/(edge*(1-edge));em=c.init.expminus(smin-4);values=[]
            for n in range(6):
                need(smin>=2*n,'flat-tail monotonicity at fifth derivative')
                cap=smin**(2*n)*em.hi
                values.append(self.I(0,cap) if n==0 else self.I(-cap,cap))
        else:
            s=1/(F(1,4)-(r-F(3,2))**2);v=3-2*r;ex=c.init.expbox(s-4)
            values=[c.init.centered_poly(polys[n],s)*ex*(v if n%2 else 1) for n in range(6)]
        return self.J([values[n]/(self.ell**n*math.factorial(n)) for n in range(6)])

    def packet(self,r):
        c=self.ctx;rr=self.J([r,1/self.ell]+[0]*4)
        if r.hi<=F(65,64):
            actual=self.bump(r,True)
            return self.P(actual,self.constant(0),[self.constant(0),actual,self.constant(0)])
        need(r.lo>=F(65,64) and not(r.lo<F(66,64)<r.hi),'packet branch crossing')
        dr=rr-F(3,2);den=F(1,4)-dr*dr
        den.c[0]=F(1,4)-(r-F(3,2))**2
        q=self.kernel(1/den-4)
        if r.lo<F(66,64):
            tt=64*rr-65
            S=126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9
            q=q*self.P(self.constant(1),S,[self.constant(0),self.constant(1)-S,self.constant(0)])
        oracle=self.bump(r)
        q.a=self.J([self.I(max(a.lo,b.lo),min(a.hi,b.hi)) for a,b in zip(q.a.c,oracle.c)])
        return q

    def paired(self,r):
        r=self.I.cv(r);need(1<=r.lo<=r.hi<=F(7,4),'retained active footprint')
        bg=self.background(r)
        E=self.kernel(self.ell*self.J([r,1/self.ell]+[0]*4))
        u,v=self.ctx._registered(bg,E,self.packet(r))
        need(all(z.n>=2 for z in u) and all(z.n>=1 for z in v),'required output jet orders')
        return u,v

    def nominal(self,r,branch=None):
        c=self.ctx;r=self.I.cv(r);order=6
        if branch=='zero' or (branch is None and r.hi<=F(65,64)):
            return [c.constant(0,order)]*23,[c.constant(0,order)]*23
        bg,rr,z=c.nominal_background(r,order)
        dr=rr-F(3,2);den=F(1,4)-dr*dr;den.c[0]=F(1,4)-(r-F(3,2))**2
        packet=c.square40(1/(1+(1/den-4)/c.power))
        if branch=='cutoff' or (branch is None and r.lo<F(66,64)):
            tt=64*rr-65
            packet*=126*tt**5-420*tt**6+540*tt**7-315*tt**8+70*tt**9
        E=c.square40(1/(1+z/c.power))
        return c._registered(bg,E,packet)

    def authenticate(self):
        self.ctx.reauthenticate();pin(self.producer,PRODUCER_PIN)
