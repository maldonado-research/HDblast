#!/usr/bin/env python3
"""Exact Laurent-polynomial incoming-state transfer checks; no physical data."""
import argparse
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

# Exponents are (k,L,c,p,q). Negative k,L exponents are permitted. All
# calculations below use exact rational arithmetic, not sampled identities.
N = 5
class Poly:
    def __init__(self, terms=None):
        self.terms = {m:F(v) for m,v in (terms or {}).items() if v}
    @staticmethod
    def coerce(v):
        return v if isinstance(v,Poly) else Poly({(0,)*N:F(v)})
    def __add__(self,o):
        d=dict(self.terms)
        for m,v in self.coerce(o).terms.items(): d[m]=d.get(m,F(0))+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self): return Poly({m:-v for m,v in self.terms.items()})
    def __sub__(self,o): return self+-self.coerce(o)
    def __rsub__(self,o): return self.coerce(o)+-self
    def __mul__(self,o):
        d={}
        for m,v in self.terms.items():
            for n,w in self.coerce(o).terms.items():
                e=tuple(a+b for a,b in zip(m,n));d[e]=d.get(e,F(0))+v*w
        return Poly(d)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            if len(self.terms)!=1: raise ValueError('Only monomial negative powers supported')
            m,v=next(iter(self.terms.items()))
            return Poly({tuple(n*e for e in m):v**n})
        out=self.coerce(1)
        for _ in range(n):out=out*self
        return out
    def diff(self,i):
        out={}
        for m,v in self.terms.items():
            if m[i]:
                e=list(m);e[i]-=1;out[tuple(e)]=v*m[i]
        return Poly(out)
    def subst_zero(self,i):
        return Poly({m:v for m,v in self.terms.items() if not m[i]})
    def value(self,values):
        return sum((v*__import__('functools').reduce(lambda a,z:a*z,
            (F(x)**e for x,e in zip(values,m)),F(1)) for m,v in self.terms.items()),F(0))

def var(i):
    e=[0]*N;e[i]=1;return Poly({tuple(e):1})
k,L,c,p,q=(var(i) for i in range(N))
invk=k**-1
def D(e): return L**2*e.diff(1)-2*k*q*e.diff(3)+2*k*p*e.diff(4)

def require(ok,label):
    if not ok: raise RuntimeError(label)

def rational(v): return f'{v.numerator}/{v.denominator}'
def upper_decimal(v):
    with localcontext() as ctx:
        ctx.prec=12;ctx.rounding=ROUND_CEILING
        return str(Decimal(v.numerator)/Decimal(v.denominator))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    require(not a.output.exists(),'Fresh output required')
    checks=[];mutations=[]
    def eq(name,x,y=0):
        require(not (x-y).terms,name);checks.append(name)
    def reject(name,e):
        require(bool(e.terms),'Mutation survived: '+name);mutations.append(name)
    X=c+p;Z=-2*k*q;T=2*k*p
    R=((2*k**2+3*L**2)*X-k*T-L*Z)*invk*F(1,2)
    P=((F(2,3)*k**2-L**2)*X-k*T-L*Z)*invk*F(1,2)
    Rc=(k+F(3,2)*L**2*invk)*c
    Pc=(F(1,3)*k-F(1,2)*L**2*invk)*c
    Ra=F(3,2)*L**2*invk*p+L*q
    Pa=-(F(2,3)*k+F(1,2)*L**2*invk)*p+L*q
    eq('Direct density canonical and phase decomposition',R,Rc+Ra)
    eq('Direct pressure canonical and phase decomposition',P,Pc+Pa)
    eq('Real homogeneous U flow',D(X),Z)
    eq('Real homogeneous W flow',D(Z),-2*k*T)
    eq('Imaginary homogeneous W flow',D(T),2*k*Z)
    eq('Canonical normalization coefficient',X-T*invk*F(1,2),c)
    eq('Canonical coefficient is constant',D(c))
    eq('Homogeneous phase magnitude invariant',D(p**2+q**2))
    eq('Direct incoming-state Ward identity',D(R),L*(R-3*P))
    eq('Canonical sector direct Ward identity',D(Rc),L*(Rc-3*Pc))
    eq('Phase sector direct Ward identity',D(Ra),L*(Ra-3*Pa))
    eq('Direct pressure time primitive',D(-F(1,3)*R*L**-1),P)
    eq('Canonical direct pressure time primitive',D(-F(1,3)*Rc*L**-1),Pc)
    eq('Phase direct pressure time primitive',D(-F(1,3)*Ra*L**-1),Pa)
    eq('Squared density phase coefficient',
       (F(3,2)*L**2*invk)**2+L**2,
       L**2*(1+F(9,4)*L**2*invk**2))
    eq('Squared pressure phase coefficient',
       (F(2,3)*k+F(1,2)*L**2*invk)**2+L**2,
       F(4,9)*k**2+F(5,3)*L**2+F(1,4)*L**4*invk**2)
    # The following exact nonnegative differences prove the pointwise norm
    # bounds because all candidate coefficients are positive for k,L>0.
    rbound=L+F(9,8)*L**3*invk**2
    pbound=F(2,3)*k+F(5,4)*L**2*invk
    eq('Density coefficient majorant square remainder',
       rbound**2-((F(3,2)*L**2*invk)**2+L**2),
       F(81,64)*L**6*invk**4)
    eq('Pressure coefficient majorant square remainder',
       pbound**2-((F(2,3)*k+F(1,2)*L**2*invk)**2+L**2),
       F(21,16)*L**4*invk**2)
    # Pi^2 is factored out. Differentiate the finite-K polynomial bounds;
    # k now denotes the upper cutoff variable for these primitive checks.
    CRa=L*k**3*F(1,6)+L**3*k*F(9,16)
    CPa=k**4*F(1,12)+L**2*k**2*F(5,16)
    CRc=k**4*F(1,8)+L**2*k**2*F(3,8)
    CPc=k**4*F(1,24)+L**2*k**2*F(1,8)
    for name,C,majorant in [
        ('phase density',CRa,rbound),('phase pressure',CPa,pbound),
        ('canonical density',CRc,k+F(3,2)*L**2*invk),
        ('canonical pressure',CPc,F(1,3)*k+F(1,2)*L**2*invk)]:
        eq(name+' bound momentum primitive',C.diff(0),k**2*F(1,2)*majorant)
        require(all(e[0]>0 for e in C.terms),name+' primitive has zero lower limit')
        checks.append(name+' primitive zero lower limit')
    reject('Drop density phase term',R-(Rc+F(3,2)*L**2*invk*p))
    reject('Drop pressure phase term',P-(Pc-(F(2,3)*k+F(1,2)*L**2*invk)*p))
    reject('Drop canonical density state term',R-Ra)
    reject('Drop canonical pressure state term',P-Pa)
    reject('Wrong density phase sign',R-(Rc+F(3,2)*L**2*invk*p-L*q))
    reject('Wrong pressure leading coefficient',P-(Pc-(F(1,3)*k+F(1,2)*L**2*invk)*p+L*q))
    reject('Wrong pressure primitive sign',D(F(1,3)*R*L**-1)-P)
    reject('Freeze geometry in state work',D(R)-(-2*k*q*R.diff(3)+2*k*p*R.diff(4)))
    reject('Drop geometric canonical work',D(Rc))
    reject('Omit density norm geometry correction',rbound-L)
    reject('Omit pressure norm geometry correction',pbound-F(2,3)*k)
    eq('Canonical density leading k term has zero time derivative',D(k*c))
    La,Lb=F(2,9),F(2,7);pi2lower=F(9)
    sigma_example=F(1,10**16);chi_example=F(1,10**16)
    rows=[]
    for K in (64,128,256):
        def evaluate(expr,l): return expr.value((K,l,0,0,0))/pi2lower
        r_a=evaluate(CRa,Lb);p_a=evaluate(CPa,Lb)
        r_c=evaluate(CRc,Lb);p_c=evaluate(CPc,Lb)
        work_a=evaluate(CRa,La)+r_a
        work_c=F(3,8)*(Lb**2-La**2)*K**2/pi2lower
        intp_a=(evaluate(CRa,La)/La+r_a/Lb)/3
        intp_c=F(K**4,24)/pi2lower+(Lb-La)*K**2/(8*pi2lower)
        vals={'normalized_phase_density':r_a,'normalized_phase_pressure':p_a,
              'canonical_density':r_c,'canonical_pressure':p_c,
              'normalized_phase_continuity_work':work_a,
              'canonical_continuity_work':work_c,
              'normalized_phase_time_integrated_pressure':intp_a,
              'canonical_time_integrated_pressure':intp_c}
        require(r_a==F(K**3,189)+F(K,686),'Density all-time rational bound')
        require(p_a==F(K**4,108)+F(5*K**2,1764),'Pressure all-time rational bound')
        checks.extend([f'K={K} exact density upper coefficient',f'K={K} exact pressure upper coefficient'])
        rows.append({'K':K,'uniform_norm_coefficients_upper':{n:rational(v) for n,v in vals.items()},
                     'illustration_only_sigma_and_chi':'1/10000000000000000',
                     'phase_error_example_upper_decimal':{n:upper_decimal(v*sigma_example) for n,v in vals.items() if n.startswith('normalized_phase')},
                     'canonical_error_example_upper_decimal':{n:upper_decimal(v*chi_example) for n,v in vals.items() if n.startswith('canonical')}})
    out={'status':'PASS_EXACT_INCOMING_STATE_TRANSFER',
         'check_count':len(checks),'checks':checks,'mutation_count':len(mutations),'mutations_rejected':mutations,
         'physical_source_evaluations':0,'saved_array_decodes':0,'physical_trajectory_evaluations':0,
         'full_twelve_case_pressure_contact_certificate':'UNRESOLVED',
         'metric_calibration':'FAIL','higher_dimensional_Big_Bang':'NOT_ESTABLISHED',
         'external_novelty':'NOT_ASSESSED','actual_incoming_state_error':'NOT_ENCLOSED',
         'python_optimization':sys.flags.optimize,'python_version':sys.version.split()[0],
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'interval':['-9/2','-7/2'],'fixed_measure':'k^2 dk/(2 Pi^2), same declared Pi everywhere',
         'rational_bound_premise':'Pi >= 3; exact fixed geometry; identical real forcing and complete contacts',
         'scope':'Incoming-state homogeneous contribution only, conditional on independently established continuum state error envelopes. Uniform example norms are hypothetical sensitivities, not achieved errors, acceptance gates or evidence of physical accuracy.',
         'sensitivity_rows':rows}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({n:out[n] for n in ['status','check_count','mutation_count','physical_source_evaluations','saved_array_decodes']}))
if __name__=='__main__':main()
