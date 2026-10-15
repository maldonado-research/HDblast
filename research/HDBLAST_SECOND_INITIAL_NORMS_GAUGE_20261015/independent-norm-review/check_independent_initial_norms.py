#!/usr/bin/env python3
"""Independent physical second-Cauchy norm enclosure using derivative arrays.

All array entries here are UNNORMALIZED derivatives D_z^j, unlike the
producer's normalized Taylor arithmetic. No project program is imported.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb, factorial
from pathlib import Path

GRID=1<<96
INPUT_HASHES={
    'background':'3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22',
    'bump':'add927ba57117a967da7d6725821d95e894ef954405d503731b3055fc64377bd',
    'initial_jet_derivation':'2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed'}
ANCESTORS={
    'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md':'0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c'}
CAPS={
    'gamma2_components':'1.06e24','K2_components':'7.96e7','phi2':'5.13e33','Pi2':'2e11',
    'gamma2_radial_covariant_components':'4.25e29',
    'gamma2_angular_tangent_covariant_components':'1.34e22',
    'gamma2_angular_mixed_covariant_components':'1.18e23',
    'gamma2_full_spatial_covariant_components':'4.25e29',
    'phi2_radial_derivative':'3.06e39','phi2_angular_derivative':'6.48e31',
    'phi2_full_spatial_derivative':'3.06e39',
    'K2_radial_covariant_components':'1.96e13','K2_angular_covariant_components':'1.23e6',
    'K2_full_spatial_covariant_components':'1.96e13',
    'Pi2_radial_derivative':'5.02e16','Pi2_angular_derivative':'3.07e9',
    'Pi2_full_spatial_derivative':'5.02e16'}


def need(condition,message):
    if not condition:raise ValueError(message)


class Box:
    def __init__(self,lo=0,hi=None):
        lo=Q(lo);hi=lo if hi is None else Q(hi)
        need(lo<=hi,'ordered interval')
        self.lo=Q((lo*GRID).numerator//(lo*GRID).denominator,GRID)
        self.hi=Q(-((-(hi*GRID).numerator)//(hi*GRID).denominator),GRID)
    @staticmethod
    def of(x):return x if isinstance(x,Box) else Box(x)
    def __add__(self,other):
        other=Box.of(other);return Box(self.lo+other.lo,self.hi+other.hi)
    __radd__=__add__
    def __neg__(self):return Box(-self.hi,-self.lo)
    def __sub__(self,other):return self+-Box.of(other)
    def __rsub__(self,other):return Box.of(other)+-self
    def __mul__(self,other):
        other=Box.of(other)
        endpoints=[a*b for a in (self.lo,self.hi) for b in (other.lo,other.hi)]
        return Box(min(endpoints),max(endpoints))
    __rmul__=__mul__
    def reciprocal(self):
        need(self.lo>0 or self.hi<0,'reciprocal stays off zero')
        return Box(1/self.hi,1/self.lo)
    def square(self):
        return Box(0 if self.lo<=0<=self.hi else min(self.lo*self.lo,self.hi*self.hi),
                   max(self.lo*self.lo,self.hi*self.hi))
    def absolute(self):return max(abs(self.lo),abs(self.hi))
    def contains(self,value):return self.lo<=value<=self.hi


class Derivatives:
    def __init__(self,rows):self.d=[Box.of(x) for x in rows]
    @property
    def order(self):return len(self.d)-1
    def other(self,x):return x if isinstance(x,Derivatives) else Derivatives([x]+[0]*self.order)
    def __add__(self,x):
        x=self.other(x)
        return Derivatives([self.d[j]+x.d[j] for j in range(min(self.order,x.order)+1)])
    __radd__=__add__
    def __neg__(self):return Derivatives([-x for x in self.d])
    def __sub__(self,x):return self+-self.other(x)
    def __rsub__(self,x):return self.other(x)+-self
    def __mul__(self,x):
        x=self.other(x)
        return Derivatives([sum((comb(n,j)*self.d[j]*x.d[n-j] for j in range(n+1)),Box(0))
                            for n in range(min(self.order,x.order)+1)])
    __rmul__=__mul__
    def square(self):
        result=self*self;result.d[0]=self.d[0].square();return result
    def reciprocal(self):
        result=[self.d[0].reciprocal()]
        for n in range(1,self.order+1):
            result.append(-result[0]*sum((comb(n,j)*self.d[j]*result[n-j] for j in range(1,n+1)),Box(0)))
        return Derivatives(result)
    def __truediv__(self,x):return self*self.other(x).reciprocal()
    def differentiate(self):
        need(self.order>0,'derivative needs one more input order')
        return Derivatives(self.d[1:])
    def absolute_value(self):return self.d[0].absolute()


def arithmetic_checks():
    # Exact polynomial and reciprocal derivatives exercise the unnormalized
    # binomial arithmetic independently of any physical expression.
    x=Derivatives([1,1,0,0,0])
    power=(x*x)*(x*x)
    expected=[1,4,12,24,24]
    need(all(power.d[j].contains(Q(v)) for j,v in enumerate(expected)),
         'unnormalized product/binomial derivative identity')
    inverse=x.reciprocal()
    need(all(inverse.d[j].contains(Q((-1)**j*factorial(j))) for j in range(5)),
         'unnormalized reciprocal derivative identity')
    need(all(power.differentiate().d[j].contains(expected[j+1]) for j in range(4)),
         'derivative-array shift has no extra Taylor factor')
    for value in (Q(1,3),Q(-2,7),Q(12345,17)):
        need(Box(value).contains(value),'outward rational input enclosure')
    return ['unnormalized binomial products','unnormalized reciprocal derivatives',
            'derivative-array shift','outward rational endpoint enclosures']


def reduction_checks():
    # Formal Laurent expansion, including p^-1, verifies every coefficient
    # of the two reductions. This is not a seeded point-jet check.
    variables=('u','v','j','a','p','h','b','E','U2')
    count=len(variables)
    class Polynomial:
        def __init__(self,terms=0):
            self.t=({k:Q(v) for k,v in terms.items() if v} if isinstance(terms,dict)
                    else {(0,)*count:Q(terms)} if terms else {})
        def __add__(self,x):
            x=cast(x);terms=self.t.copy()
            for k,v in x.t.items():terms[k]=terms.get(k,Q(0))+v
            return Polynomial(terms)
        __radd__=__add__
        def __neg__(self):return Polynomial({k:-v for k,v in self.t.items()})
        def __sub__(self,x):return self+-cast(x)
        def __rsub__(self,x):return cast(x)+-self
        def __mul__(self,x):
            x=cast(x);terms={}
            for i,a in self.t.items():
                for j,b in x.t.items():
                    key=tuple(v+w for v,w in zip(i,j))
                    terms[key]=terms.get(key,Q(0))+a*b
            return Polynomial(terms)
        __rmul__=__mul__
        def __pow__(self,n):
            if n<0:
                need(len(self.t)==1,'Laurent negative power is a monomial')
                k,v=next(iter(self.t.items()))
                return Polynomial({tuple(n*j for j in k):v**n})
            result=Polynomial(1)
            for _ in range(n):result=result*self
            return result
        def __truediv__(self,x):return self*cast(x)**-1
        def __rtruediv__(self,x):return cast(x)*self**-1
    def cast(x):return x if isinstance(x,Polynomial) else Polynomial(x)
    def variable(i):
        powers=[0]*count;powers[i]=1;return Polynomial({tuple(powers):1})
    u,v,j,a,p,h,b,E,U2=[variable(i) for i in range(count)]
    ay=-b*b-p*p/3
    uy=-p*v/3-2*a*u
    uyy=-p*j/3+p*(2*a-h/3)*v+(4*a*a-2*ay)*u
    d=3*b/p
    ut=E*u;vt=E*(v-d*u)
    original_A=(-24*uy*uy-288*a*u*uy-24*u*uyy-(192*ay+384*a*a+32*p*p)*u*u
                -b*b*(336*u*u+12*ut*ut+2*vt*vt)-2*j*j-16*p*u*j-2*U2*v*v)
    reduced_A=(-2*j*j-8*p*u*j+8*p*(h+2*a)*u*v+(16*p*p-192*b*b)*u*u
               -(2*U2+Q(8,3)*p*p)*v*v-b*b*(12*ut*ut+2*vt*vt))
    need(not (original_A-reduced_A).t,'all-point reduced Hamiltonian Laurent residual')
    uyt=E*(uy+b*u)
    dt=d*(3*a-h)
    vty=E*(j-dt*u-d*uy+b*(v-d*u))
    F=(18*u*ut-v*vt)*b
    Fy=(18*uy*ut+18*u*uyt-j*vt-v*vty)*b-a*F
    G=(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*j)*b
    reduced_P=E*b/p*(d*u*j+(d*(h-b)+10*p)*u*v+(30*a-6*b)*u*u+(2*b-a)*v*v)
    need(not ((G-Fy)/p-reduced_P).t,'all-point reduced scalar-velocity Laurent residual')
    need(not (F/3-E*b*(18*u*u-v*v+d*u*v)/3).t,'all-point pure-trace K Laurent residual')
    return ['Hamiltonian source reduction','scalar velocity cancellation','pure-trace K reduction']


def authenticate(paths):
    objects={};digests={}
    for name,path in paths.items():
        raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
        need(digest==INPUT_HASHES[name],'authenticated actual dependency: '+name)
        digests[name]=digest
        if name!='initial_jet_derivation':objects[name]=json.loads(raw)
    bg,bump=objects['background'],objects['bump']
    need(bg['status']=='PASS_EXACT_ACTUAL_ANNULUS_STATIC_COEFFICIENT_JETS','actual static background status')
    need(bg['domain']['t']=='0' and Q(bg['domain']['z_min'])==Q(1,10**6) and
         Q(bg['domain']['z_max'])==Q(1,500000) and bg['domain']['ell']=='1/1000000',
         'whole original initial annulus')
    need(bg['input_sha256']==ANCESTORS,'static background dependency chain')
    model=bg['actual_model']
    need(model['eta']=='phi-1' and model['delta']=='1/1000' and
         model['c']=='96714321425692/161839258930731' and
         model['U_coefficients']==['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27'],
         'registered actual model and polynomial potential')
    need(bg['normalization']=='coefficient j = D_z^j(value)/j! or D_y^j(value)/j!, as named',
         'normalized input coefficient convention')
    need(bg['coefficient_definitions']['h']=='U_1/p' and bg['coefficient_definitions']['B']=='3alpha/(rho^3 p)',
         'actual reconstruction coefficient convention')
    need(bump['status']=='PASS' and bump['ell']=='1/1000000' and bump['maximum_derivative_order']>=4,
         'exact dimensionless bump fourth derivative input')
    return bg,bump,digests


def load_derivatives(background,key):
    rows=background['z_jets'][key]
    need(len(rows)>=5,'fourth-order background derivative supply: '+key)
    return Derivatives([Box(*row['exact'])*factorial(j) for j,row in enumerate(rows)])


def independent_enclosures(bg,bump):
    ell=Q(1,10**6)
    norms=list(map(Q,bump['dimensionless_Q_derivative_sup']))
    physical=[norms[j]/ell**j for j in range(5)]
    q=Derivatives([Box(0,1)]+[Box(-x,x) for x in physical[1:]])
    rho,a,p,C,B,k,h,U2,b,g,E=[load_derivatives(bg,name) for name in
        ('rho','a','p','C','B','k','h','U_2','rho_inv','p_inv','exp_minus_z')]
    need(rho.d[0].lo>0 and p.d[0].hi<0,'actual radius and scalar-gradient gaps')
    # Use the independently certified reciprocal jets directly. The producer
    # recomputes these from rho/p. Both enclose the same actual functions.
    d=3*b*g
    def y(x):return -(b*x.differentiate())
    u=C*q
    v=B*(q.differentiate()+k*q)
    j=y(v)
    ut=E*u
    vt=E*(v-d*u)
    A=(-2*j.square()-8*p*u*j+8*p*(h+2*a)*u*v+(16*p.square()-192*b.square())*u.square()
       -(2*U2+Q(8,3)*p.square())*v.square()-b.square()*(12*ut.square()+2*vt.square()))
    D=b.square()*(36*u.square()-2*v.square())
    sources=[A/6+Q(4,3)*D,A/2+2*D,A-4*D]
    f=[-rho.square()*source/den for source,den in zip(sources,(-6,10,42))]
    fy=[y(x) for x in f]
    zeta=[-Q(3,2)*g*x for x in fy]
    zetay=[y(x) for x in zeta]
    K=E*b*(18*u.square()-v.square()+d*u*v)/3
    P=E*b*g*(d*u*j+(d*(h-b)+10*p)*u*v+(30*a-6*b)*u.square()+(2*b-a)*v.square())
    Ky,Py=y(K),y(P)
    weights=(Q(1),Q(1,2),Q(7,12))
    def sup(xs):return sum((weights[n]*xs[n].absolute_value() for n in range(3)),Q(0))
    def angular(xs):return xs[1].absolute_value()+Q(5,2)*xs[2].absolute_value()
    rmin=rho.d[0].lo
    gamma,radial=sup(f),sup(fy)
    gamma_tan=angular(f)/rmin
    gamma_mixed=a.d[0].absolute()*gamma
    phi,phi_radial,phi_ang=sup(zeta),sup(zetay),angular(zeta)/rmin
    Kval,Pval=K.absolute_value(),P.absolute_value()
    Kr,Pr=Ky.absolute_value(),Py.absolute_value()
    Ka,Pa=2*Kval/rmin,2*Pval/rmin
    bounds={'gamma2_components':gamma,'K2_components':Kval,'phi2':phi,'Pi2':Pval,
        'gamma2_radial_covariant_components':radial,
        'gamma2_angular_tangent_covariant_components':gamma_tan,
        'gamma2_angular_mixed_covariant_components':gamma_mixed,
        'gamma2_full_spatial_covariant_components':max(radial,gamma_tan,gamma_mixed),
        'phi2_radial_derivative':phi_radial,'phi2_angular_derivative':phi_ang,
        'phi2_full_spatial_derivative':max(phi_radial,phi_ang),
        'K2_radial_covariant_components':Kr,'K2_angular_covariant_components':Ka,
        'K2_full_spatial_covariant_components':max(Kr,Ka),
        'Pi2_radial_derivative':Pr,'Pi2_angular_derivative':Pa,
        'Pi2_full_spatial_derivative':max(Pr,Pa)}
    need(set(bounds)==set(CAPS),'complete list of published bounds')
    for key,value in bounds.items():need(value<Q(CAPS[key]),'independent enclosure below published strict cap: '+key)
    need(A.order>=2 and all(x.order>=2 for x in f) and all(x.order>=1 for x in zeta),
         'all required second-source and first-field derivatives retained')
    return bounds,{'A':A.order,'f':[x.order for x in f],'phi2':[x.order for x in zeta],'Pi2':P.order}


def compare_report(path):
    raw=path.read_bytes();report=json.loads(raw)
    need(report['status']=='PASS_ACTUAL_SECOND_INITIAL_VALUE_AND_FIRST_SPATIAL_BOUNDS','producer physical initial bound status')
    need(report['input_sha256']==INPUT_HASHES,'producer uses precisely the reviewed actual inputs')
    need(set(report['bounds'])==set(CAPS),'all published component/direction caps covered')
    for key,value in CAPS.items():
        row=report['bounds'][key]
        need(Q(row['strict_simple_upper'])==Q(value),'reviewed rounded cap: '+key)
        need(Q(row['exact_upper'])<Q(value),'producer exact enclosure is strictly below its rounded cap: '+key)
    scope=report['scope']
    need(scope['full_first_spatial_derivatives'] and not scope['Sobolev_energy_or_evolution_norm'] and
         not scope['ADM_momentum_densities'] and not scope['second_evolution_computed'],
         'precise physical first-spatial component scope')
    return hashlib.sha256(raw).hexdigest()


def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument('--background',type=Path,default=root/'producer/inputs/BACKGROUND_RECONSTRUCTION_JETS.json')
    ap.add_argument('--bump',type=Path,default=root/'producer/inputs/BUMP_DERIVATIVE_BOUNDS.json')
    ap.add_argument('--initial-jet',type=Path,default=root/'producer/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md')
    ap.add_argument('--producer-report','--producer',dest='producer_report',type=Path,
                    default=root/'producer/SECOND_INITIAL_NORMS_NORMAL.json')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    paths={'background':args.background,'bump':args.bump,'initial_jet_derivation':args.initial_jet}
    bg,bump,pins=authenticate(paths)
    checks=arithmetic_checks()
    reductions=reduction_checks()
    producer_hash=compare_report(args.producer_report)
    bounds,orders=independent_enclosures(bg,bump)
    result={'status':'PASS_INDEPENDENT_PHYSICAL_SECOND_INITIAL_NORM_REVIEW','optimized':not __debug__,
        'arithmetic':'unnormalized derivative arrays with binomial products; outward rational grid 2^96',
        'norm':'Supremum of absolute individual tensor components, with each first spatial covariant derivative index in the adapted background frame (partial_y,rho^-1 E_i). This is not a summed tensor norm.',
        'input_sha256':pins,'input_ancestor_sha256':ANCESTORS,
        'reviewed_producer_report_sha256':producer_hash,'arithmetic_identity_checks':checks,
        'exact_zero_Laurent_reductions':reductions,
        'independently_certified_rounded_caps':len(bounds),'retained_derivative_orders':orders,
        'bounds':{key:{'independent_upper_exact':str(value),'published_strict_cap':str(Q(CAPS[key])),
                       'display_only':float(value)} for key,value in bounds.items()},
        'scope':{'all_individual_adapted_frame_components':True,'all_first_spatial_covariant_directions':True,
            'physical_Cauchy_variables_not_momentum_densities':True,'coordinate_component_norm':False,
            'Sobolev_or_energy_norm':False,'gauge_completion_certified_by_this_review':False,
            'second_evolution_certified':False,'nonlinear_amplitude_certified':False}}
    text=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')


if __name__=='__main__':main()
