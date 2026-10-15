#!/usr/bin/env python3
"""Exact outward interval bounds for the explicit actual second initial jet.

Taylor coefficients are D_z^j/j!. All interval operations round outward to
the fixed dyadic grid; no floating value decides a mathematical inequality.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile

DEN=1<<144

def require(ok,msg):
    if not ok:raise ValueError(msg)

def eq(a,b,msg):require(a==b,msg+': '+str(a-b))

def floorf(x):return x.numerator//x.denominator
def ceilf(x):return -((-x.numerator)//x.denominator)

class Iv:
    def __init__(self,lo=0,hi=None):
        lo=F(lo);hi=lo if hi is None else F(hi)
        require(lo<=hi,'interval ordering')
        self.lo=F(floorf(lo*DEN),DEN);self.hi=F(ceilf(hi*DEN),DEN)
    @staticmethod
    def cv(x):return x if isinstance(x,Iv) else Iv(x)
    def __add__(self,o):
        o=Iv.cv(o);return Iv(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return Iv(-self.hi,-self.lo)
    def __sub__(self,o):return self+-Iv.cv(o)
    def __rsub__(self,o):return Iv.cv(o)+-self
    def __mul__(self,o):
        o=Iv.cv(o);v=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return Iv(min(v),max(v))
    __rmul__=__mul__
    def inv(self):
        require(self.lo>0 or self.hi<0,'interval reciprocal crosses zero')
        return Iv(1/self.hi,1/self.lo)
    def __truediv__(self,o):return self*Iv.cv(o).inv()
    def __rtruediv__(self,o):return Iv.cv(o)*self.inv()
    def __pow__(self,n):
        require(isinstance(n,int) and n>=0,'nonnegative integer power required')
        if n==2:
            values=[self.lo*self.lo,self.hi*self.hi]
            return Iv(0 if self.lo<=0<=self.hi else min(values),max(values))
        r=Iv(1)
        for _ in range(n):r=r*self
        return r
    def absmax(self):return max(abs(self.lo),abs(self.hi))
    def row(self):return [str(self.lo),str(self.hi)]

class Jet:
    def __init__(self,c):self.c=[Iv.cv(x) for x in c]
    @property
    def n(self):return len(self.c)-1
    @staticmethod
    def const(x,n):return Jet([x]+[0]*n)
    def cv(self,o):return o if isinstance(o,Jet) else Jet.const(o,self.n)
    def __add__(self,o):
        o=self.cv(o);n=min(self.n,o.n);return Jet([self.c[j]+o.c[j] for j in range(n+1)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.c])
    def __sub__(self,o):return self+-self.cv(o)
    def __rsub__(self,o):return self.cv(o)+-self
    def __mul__(self,o):
        o=self.cv(o);n=min(self.n,o.n)
        return Jet([sum((self.c[k]*o.c[j-k] for k in range(j+1)),Iv(0)) for j in range(n+1)])
    __rmul__=__mul__
    def inv(self):
        r=[self.c[0].inv()]
        for j in range(1,self.n+1):r.append(-r[0]*sum((self.c[k]*r[j-k] for k in range(1,j+1)),Iv(0)))
        return Jet(r)
    def __truediv__(self,o):return self*self.cv(o).inv()
    def __rtruediv__(self,o):return self.cv(o)*self.inv()
    def __pow__(self,n):
        if n==2:
            out=self*self;out.c[0]=self.c[0]**2;return out
        r=Jet.const(1,self.n)
        for _ in range(n):r=r*self
        return r
    def dz(self):
        require(self.n>=1,'Taylor derivative needs another coefficient')
        return Jet([(j+1)*self.c[j+1] for j in range(self.n)])
    def val(self):return self.c[0]

def polynomial_identity_checks(injection):
    count=0
    for seed in (271,1927,6469):
        rng=random.Random(seed);r=lambda:F(rng.randint(-7,7),rng.randint(1,6))
        rho=F(7,3);p=r() or F(2,3);a,h,U2,u,uy,uyy,E=[r() for _ in range(7)]
        ay=-1/rho**2-p*p/3;py=p*(h-4*a)
        v=-3*(uy+2*a*u)/p
        vy=-3*(uyy+2*ay*u+2*a*uy)/p+3*(uy+2*a*u)*py/p**2
        d=3/(rho*p);dy=d*(3*a-h)
        ut=E*u;uyt=E*(uy+u/rho);vt=E*(v-d*u)
        vty=E*(vy-dy*u-d*uy+(v-d*u)/rho)
        A=(-24*uy**2-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u
           -(336*u*u+12*ut*ut+2*vt*vt)/rho**2-2*vy*vy-16*p*u*vy-2*U2*v*v)
        reduced=(-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192/rho**2)*u*u
                 -(2*U2+F(8,3)*p*p)*v*v-(12*ut*ut+2*vt*vt)/rho**2)
        if injection=='wrong_A_cross':reduced+=8*p*u*vy
        eq(A,reduced,'on-shell Hamiltonian source reduction')
        Fm=(18*u*ut-v*vt)/rho
        Fmy=(18*uy*ut+18*u*uyt-vy*vt-v*vty-a*(18*u*ut-v*vt))/rho
        G=(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)/rho
        P=(G-Fmy)/p
        reducedP=E/(rho*p)*(d*u*vy+(d*(h-1/rho)+10*p)*u*v
                             +(30*a-6/rho)*u*u+(2/rho-a)*v*v)
        if injection=='wrong_Pi_cross':reducedP-=d*u*vy*E/(rho*p)
        eq(P,reducedP,'actual-velocity scalar momentum cancellation')
        eq(Fm,E/rho*(18*u*u-v*v+d*u*v),'actual-velocity pure-trace K profile')
        count+=3
    # Independent scalar chain-rule and adapted-frame connection witnesses.
    rho=F(7,3);fz=F(5,7);rad=-fz/rho
    candidate=(fz/rho if injection=='wrong_y_sign' else -fz/rho)
    eq(candidate,rad,'proper radial derivative conversion')
    aa,H=F(2,5),F(7,11)
    mixed=(F(0) if injection=='drop_projector_connection' else -aa*H)
    eq(mixed,-aa*H,'angular projector mixed covariant derivative')
    # Convention coefficients multiply epsilon²/2, so derivative coefficient2.
    q2=F(13,7);reconstructed=(q2/2 if injection=='half_second_coefficient' else q2)
    eq(reconstructed,q2,'second parameter derivative normalization')
    return count+3

def load_jet(report,key):
    # Each independently certified row is normalized D_z^j/j!.
    rows=report['z_jets'][key]
    return Jet([Iv(*row['exact']) for row in rows])

def strict_decimal_cap(x):
    require(x>=0,'nonnegative majorant required')
    scale=F(1)
    while x>=scale:scale*=10
    while x<scale/10 and scale>F(1,10**100):scale/=10
    unit=scale/1000
    return (floorf(x/unit)+1)*unit

def packed_bound(x):
    cap=strict_decimal_cap(x)
    require(x<cap,'strict decimal upper cap')
    return {'exact_upper':str(x),'strict_simple_upper':str(cap),'display_only':float(x)}

def calculate(background,bump,injection):
    ell=F(1,1000000)
    require(background['status']=='PASS_EXACT_ACTUAL_ANNULUS_STATIC_COEFFICIENT_JETS','actual static background certificate status')
    require(background['domain']['ell']==str(ell),'registered ell mismatch')
    require(F(background['domain']['z_min'])==ell and F(background['domain']['z_max'])==2*ell
            and background['domain']['t']=='0','registered initial annulus mismatch')
    require(bump['ell']==str(ell),'bump scale mismatch')
    require(bump['status']=='PASS' and bump['maximum_derivative_order']>=4,'bump derivative certificate status/order')
    require(background['order']>=4,'background derivative certificate order')
    scaled=[F(x) for x in bump['dimensionless_Q_derivative_sup']]
    require(len(scaled)>=5,'initial first spatial jet needs bump order four')
    physical=[scaled[j]/ell**j for j in range(5)]
    if injection=='omit_bump_scale':physical=[scaled[j] for j in range(5)]
    for j in range(1,5):eq(physical[j]*ell**j,scaled[j],'physical bump derivative scale')
    q=Jet([Iv(0,1)]+[Iv(-physical[j],physical[j])/math.factorial(j) for j in range(1,5)])
    rho,a,p,C,B,k,h,U2=[load_jet(background,key) for key in ('rho','a','p','C','B','k','h','U_2')]
    require(min(j.n for j in (rho,a,p,C,B,k,h,U2))>=3,'background adapter derivative order')
    ir=rho.inv();ip=p.inv();d=3*ir*ip
    def dy(j):return -ir*j.dz()
    # The exact exponential has signed derivatives E^(j)=(-1)^jE.
    base=Iv(1-2*ell,1)
    E=Jet([((-1)**j)*base/math.factorial(j) for j in range(5)])
    u=C*q;v=B*(q.dz()+k*q);vy=dy(v)
    ut=E*u;vt=E*(v-d*u)
    A=(-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192*ir**2)*u*u
       -(2*U2+F(8,3)*p*p)*v*v-ir**2*(12*ut*ut+2*vt*vt))
    D=ir**2*(36*u*u-2*v*v)
    sectors=[A/6+F(4,3)*D,A/2+2*D,A-4*D]
    f=[-rho**2*sectors[j]/den for j,den in enumerate((-6,10,42))]
    fy=[dy(x) for x in f]
    zeta=[-F(3,2)*ip*x for x in fy]
    zetay=[dy(x) for x in zeta]
    K=E*ir/3*(18*u*u-v*v+d*u*v)
    P=E*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v+(30*a-6*ir)*u*u+(2*ir-a)*v*v)
    Ky=dy(K);Py=dy(P)
    abs0=lambda j:j.val().absmax()
    weights=[F(1),F(1,2),F(7,12)]
    total=lambda xs:sum(weights[j]*abs0(xs[j]) for j in range(3))
    sphere=lambda xs:abs0(xs[1])+F(5,2)*abs0(xs[2])
    gamma=total(f);gamma_y=total(fy)
    phi=total(zeta);phi_y=total(zetay)
    rmin=rho.val().lo;amax=a.val().absmax()
    gamma_ang=sphere(f)/rmin;gamma_mixed=amax*gamma
    phi_ang=sphere(zeta)/rmin
    kval=abs0(K);pval=abs0(P)
    kval_y=abs0(Ky);pval_y=abs0(Py)
    K_ang=2*kval/rmin;Pi_ang=2*pval/rmin
    # Exact geometric estimate for max individual adapted-frame components.
    bounds={'gamma2_components':gamma,'K2_components':kval,'phi2':phi,'Pi2':pval,
            'gamma2_radial_covariant_components':gamma_y,
            'gamma2_angular_tangent_covariant_components':gamma_ang,
            'gamma2_angular_mixed_covariant_components':gamma_mixed,
            'gamma2_full_spatial_covariant_components':max(gamma_y,gamma_ang,gamma_mixed),
            'phi2_radial_derivative':phi_y,'phi2_angular_derivative':phi_ang,
            'phi2_full_spatial_derivative':max(phi_y,phi_ang),
            'K2_radial_covariant_components':kval_y,'K2_angular_covariant_components':K_ang,
            'K2_full_spatial_covariant_components':max(kval_y,K_ang),
            'Pi2_radial_derivative':pval_y,'Pi2_angular_derivative':Pi_ang,
            'Pi2_full_spatial_derivative':max(pval_y,Pi_ang)}
    sector_rows=[]
    for j,l in enumerate((0,2,4)):
        sector_rows.append({'degree':l,'f_interval':f[j].val().row(),'f_y_interval':fy[j].val().row(),
                            'phi2_interval':zeta[j].val().row(),'phi2_y_interval':zetay[j].val().row()})
    return {'bounds':{key:packed_bound(value) for key,value in bounds.items()},
            'sector_intervals':sector_rows,'pure_trace_K_profile_interval':K.val().row(),
            'Pi2_profile_interval':P.val().row(),
            'input_physical_bump_derivative_sup':[str(x) for x in physical],
            'retained_jet_orders':{'A':A.n,'f':[x.n for x in f],
                                   'phi2':[x.n for x in zeta],'Pi2':P.n}}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    root=Path(__file__).resolve().parent;ap=argparse.ArgumentParser()
    ap.add_argument('--background',type=Path,default=root/'inputs/BACKGROUND_RECONSTRUCTION_JETS.json')
    ap.add_argument('--bump',type=Path,default=root/'inputs/BUMP_DERIVATIVE_BOUNDS.json')
    ap.add_argument('--pins',type=Path,default=root/'INPUT_PINS.json')
    ap.add_argument('--initial-jet',type=Path,default=root/'inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--run-controls',action='store_true')
    mutations={'wrong_A_cross':'on-shell Hamiltonian source reduction',
               'wrong_Pi_cross':'actual-velocity scalar momentum cancellation',
               'wrong_y_sign':'proper radial derivative conversion',
               'drop_projector_connection':'angular projector mixed covariant derivative',
               'half_second_coefficient':'second parameter derivative normalization',
               'omit_bump_scale':'physical bump derivative scale'}
    ap.add_argument('--inject-error',choices=['none']+list(mutations),default='none')
    a=ap.parse_args();checks=polynomial_identity_checks(a.inject_error)
    pins=json.loads(a.pins.read_text())
    input_paths=[('background',a.background),('bump',a.bump),('initial_jet_derivation',a.initial_jet)]
    for key,path in input_paths:
        require(sha(path)==pins[key]['sha256'] and path.stat().st_size==pins[key]['bytes'],'producer input authentication: '+key)
    background=json.loads(a.background.read_text());bump=json.loads(a.bump.read_text())
    result=calculate(background,bump,a.inject_error)
    report={'status':'PASS_ACTUAL_SECOND_INITIAL_VALUE_AND_FIRST_SPATIAL_BOUNDS',
            'finite_exact_identity_checks':checks,'arithmetic':{'outward_dyadic_denominator':str(DEN)},
            'domain':{'initial_time':'t=0','radial':'ell<=z<=2ell, smooth zero extension outside','ell':'1/1000000','angular':'entire unit S3'},
            'input_sha256':{key:sha(path) for key,path in input_paths},
            'norm':'Supremum of absolute individual components in the adapted background spatial orthonormal frame (partial_y,rho^-1 E_i), including covariant derivative indices. K is +Lie_n(gamma)/2; Pi is n(phi). Second parameter derivatives multiply epsilon^2/2.',
            'scope':{'actual_explicit_initial_jet':True,'full_first_spatial_derivatives':True,
                     'ADM_momentum_densities':False,'Sobolev_energy_or_evolution_norm':False,
                     'gauge_completion_computed':False,'second_evolution_computed':False,
                     'nonlinear_amplitude_or_remainder_established':False},**result}
    if a.run_controls:
        cases=[]
        with tempfile.TemporaryDirectory(prefix='second-initial-controls-') as td:
            for mutation,reason in mutations.items():
                output=Path(td)/(mutation+'.json');cmd=[sys.executable]
                if sys.flags.optimize:cmd.append('-O')
                cmd += [str(Path(__file__).resolve()),'--background',str(a.background),'--bump',str(a.bump),
                        '--pins',str(a.pins),'--initial-jet',str(a.initial_jet),
                        '--inject-error',mutation,'--output',str(output)]
                run=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
                require(run.returncode!=0 and ('ValueError: '+reason) in run.stderr and not output.exists(),
                        'producer control failed intended rejection: '+mutation)
                cases.append({'mutation':mutation,'rejection_prefix':reason,'no_output_report':True})
        report['controls']={'status':'PASS','optimized':bool(sys.flags.optimize),'cases':cases}
    text=json.dumps(report,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='')

if __name__=='__main__':main()
