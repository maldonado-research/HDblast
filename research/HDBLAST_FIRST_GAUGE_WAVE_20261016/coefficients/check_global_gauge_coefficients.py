#!/usr/bin/env python3
"""Portable exact/outward-rational whole-collar gauge coefficient certificate.

Python standard library only. No annulus jet report, inherited program,
floating point computation, numerical trajectory, or point sampling is used.
Normalized Taylor jets enclose derivatives at every possible collar point.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
PINS = {
 'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
 'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
 'HIGHER_JETS_NORMAL.json':'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
 'REFEREE_REVIEW.md':'1341cb2488f9fe5338d8355305b55cbcccf305597974914d6d1fdc5ace35f74e',
 'REFEREE_ACCEPTANCE_NORMAL.json':'94ea961e9dc2990cd5c0b40a61ff7824c63c8b1630e2905d367fbe320f222143',
 'FIRST_GAUGE_VECTOR_SYSTEM_UNPUBLISHED.md':'1fec12bc9f729ef7ca93fd4ec06fdce8fc1f61a2d295b6b04e3925ba609b479c',
}
ELL, TMAX, ZMAX = F(1,10**6), F(3,10**6), F(8,10**6)
ORDER = 3
U = [F(-2,27),F(0),F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)]


def need(ok, message):
    if not ok:
        raise ValueError(message)


class I:
    grid = 0
    def __init__(self, lo, hi=None):
        lo, hi = F(lo), F(lo if hi is None else hi)
        need(lo <= hi, 'reversed interval')
        if self.grid:
            q = self.grid
            lo = F((lo*q).__floor__(), q)
            hi = F(-(-hi*q).__floor__(), q)
        self.lo, self.hi = lo, hi
    @staticmethod
    def box(x):
        return x if isinstance(x,I) else I(x)
    def __add__(self, y):
        y=self.box(y); return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        y=self.box(y);v=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
        return I(min(v),max(v))
    __rmul__=__mul__
    def recip(self):
        need(not self.lo<=0<=self.hi,'zero in interval denominator')
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,y):return self*self.box(y).recip()
    def __rtruediv__(self,y):return self.box(y)*self.recip()
    def __pow__(self,n):
        need(isinstance(n,int),'integer power required')
        if n<0:return self.recip()**(-n)
        if n==0:return I(1)
        if n%2==0:
            x=self**(n//2)
            return I(0 if x.lo<=0<=x.hi else min(x.lo*x.lo,x.hi*x.hi),
                     max(x.lo*x.lo,x.hi*x.hi))
        return self*self**(n-1)
    def mag(self):return max(abs(self.lo),abs(self.hi))
    def record(self):return {'exact':[str(self.lo),str(self.hi)],'abs_upper':str(self.mag())}


class J:
    """Univariate coefficients D^n/n!, for intervals or formal polynomials."""
    def __init__(self,c):
        need(all(not isinstance(x,float) for x in c),'floating Taylor coefficient forbidden')
        self.c=[F(x) if isinstance(x,(int,F)) else x for x in c]
    @property
    def n(self):return len(self.c)
    def box(self,y):return y if isinstance(y,J) else J([y]+[0]*(self.n-1))
    def __add__(self,y):
        y=self.box(y);return J([self.c[i]+y.c[i] for i in range(min(self.n,y.n))])
    __radd__=__add__
    def __neg__(self):return J([-v for v in self.c])
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        y=self.box(y);return J([sum((self.c[j]*y.c[i-j] for j in range(i+1)),0)
                               for i in range(min(self.n,y.n))])
    __rmul__=__mul__
    def recip(self):
        out=[1/self.c[0]]
        for n in range(1,self.n):
            out.append(-out[0]*sum((self.c[j]*out[n-j] for j in range(1,n+1)),0))
        return J(out)
    def __truediv__(self,y):return self*self.box(y).recip()
    def __rtruediv__(self,y):return self.box(y)*self.recip()
    def __pow__(self,n):
        if n<0:return self.recip()**(-n)
        out=self.box(1)
        for _ in range(n):out=out*self
        return out
    def dz(self):return J([(n+1)*self.c[n+1] for n in range(self.n-1)])


class P:
    """Independent Laurent total derivatives on rho,eta,p,a,alpha."""
    n=5
    def __init__(self,d=None):self.d={m:F(c) for m,c in (d or {}).items() if c}
    @classmethod
    def box(cls,x):return x if isinstance(x,cls) else cls({(0,)*cls.n:F(x)})
    @classmethod
    def var(cls,i):
        e=[0]*cls.n;e[i]=1;return cls({tuple(e):F(1)})
    def __add__(self,y):
        d=self.d.copy()
        for m,c in self.box(y).d.items():d[m]=d.get(m,F(0))+c
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({m:-c for m,c in self.d.items()})
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        d={}
        for m,c in self.d.items():
            for k,b in self.box(y).d.items():
                e=tuple(a+v for a,v in zip(m,k));d[e]=d.get(e,F(0))+c*b
        return P(d)
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:
            need(len(self.d)==1,'formal inverse requires monomial')
            (m,c),=self.d.items();return P({tuple(a*n for a in m):c**n})
        result=P.box(1)
        for _ in range(n):result=result*self
        return result
    def __truediv__(self,y):return self*self.box(y)**-1
    def __rtruediv__(self,y):return self.box(y)*self**-1
    def diff(self,i):
        d={}
        for m,c in self.d.items():
            if m[i]:
                e=list(m);e[i]-=1;d[tuple(e)]=c*m[i]
        return P(d)


def derivative(poly):return [i*poly[i] for i in range(1,len(poly))]
UP=[U]
for _ in range(4):UP.append(derivative(UP[-1]))


def evaluate(poly,x):
    ans=0
    for c in reversed(poly):ans=ans*x+c
    return ans


def build_static(initial, injection='none'):
    state={name:J([value]+[0]*ORDER) for name,value in initial.items()}
    for n in range(ORDER):
        r,e,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')]
        h=evaluate(UP[1],e)/p
        rhs={'rho':-r*r*a,'eta':-r*p,'p':-r*(evaluate(UP[1],e)-4*a*p),
             'a':1/r+r*p*p/3,'alpha':r*(F(5,2)*a-h)*alpha}
        divisor=1 if injection=='wrong_jet_factorial' else n+1
        for name in state:state[name].c[n+1]=rhs[name].c[n]/divisor
    r,e,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')]
    out=dict(state)
    out.update({'a_y':-1/r**2-p*p/3,'rho_inv':1/r,'h':evaluate(UP[1],e)/p,
                'C':alpha/r**2,'k':r*(F(5,2)*a-evaluate(UP[1],e)/p)})
    return out


def formal_checks(injection):
    r,e,p,a,alpha=[P.var(i) for i in range(5)]
    initial=dict(zip(('rho','eta','p','a','alpha'),(r,e,p,a,alpha)))
    jets=build_static(initial,injection)
    velocity=(-r*r*a,-r*p,-r*(evaluate(UP[1],e)-4*a*p),
              1/r+r*p*p/3,r*(F(5,2)*a-evaluate(UP[1],e)/p)*alpha)
    def D(poly):return sum((poly.diff(i)*v for i,v in enumerate(velocity)),P())
    for name,series in jets.items():
        value=series.c[0]
        for n in range(1,ORDER+1):
            value=D(value)
            need(not (factorial(n)*series.c[n]-value).d,
                 'formal normalized jet identity: '+name+' order '+str(n))
    C=jets['C'].c[0];Cz=D(C)
    factor=-6 if injection=='wrong_source_sign' else 6
    ps=factor*ELL*r*C
    p0=6*ELL*ELL*r*(Cz-4*r*a*C)
    raw_ps=-6*ELL*ELL*r*r*(-C/r/ELL)
    raw_p0=-6*ELL*ELL*r*r*(-Cz/r+4*a*C)
    need(not (ps-raw_ps).d and not (p0-raw_p0).d,'proper-y/scaled-r source identity')
    return {'normalized_static_jets':'PASS_DZ_TOTAL_DERIVATIVES_1_TO_3',
            'source_chain_rule':'PASS_FRACTION_LAURENT_POLYNOMIAL_IDENTITY'}


class BJ:
    """Normalized mixed scaled-coordinate jets through total order two."""
    indices=[(i,j) for n in range(3) for i in range(n+1) for j in [n-i]]
    def __init__(self,d=None):self.d={k:I.box((d or {}).get(k,0)) for k in self.indices}
    @classmethod
    def box(cls,x):return x if isinstance(x,cls) else cls({(0,0):x})
    @classmethod
    def lift(cls,jets,axis):
        return cls({(n,0) if axis=='x' else (0,n):value*ELL**n
                    for n,value in enumerate(jets.c[:3])})
    def __add__(self,y):
        y=self.box(y);return BJ({k:self.d[k]+y.d[k] for k in self.indices})
    __radd__=__add__
    def __neg__(self):return BJ({k:-v for k,v in self.d.items()})
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        y=self.box(y);out={}
        for i,j in self.indices:
            out[i,j]=sum((self.d[a,b]*y.d[i-a,j-b]
                          for a in range(i+1) for b in range(j+1)),I(0))
        return BJ(out)
    __rmul__=__mul__
    def __pow__(self,n):
        need(n>=0,'mixed nonnegative integer powers only')
        result=self.box(1)
        for _ in range(n):result=result*self
        return result
    def deriv(self,i,j):return factorial(i)*factorial(j)*self.d[i,j]
    def record(self):
        return {f'{i},{j}':self.deriv(i,j).record() for i,j in self.indices}


def time_jets():
    H=J([I(0,TMAX)]+[0]*ORDER)
    q0=I(1-F(3,4)*TMAX*TMAX,1)
    q=J([q0]+[0]*ORDER)
    for n in range(ORDER):
        dh=1-H*H;dq=-F(3,2)*H*q
        H.c[n+1]=dh.c[n]/(n+1)
        q.c[n+1]=dq.c[n]/(n+1)
    return H,q,1-H*H


def inverse3(matrix):
    aug=[list(row)+[F(int(i==j)) for j in range(3)] for i,row in enumerate(matrix)]
    for k in range(3):
        need(aug[k][k]!=0,'zero inverse pivot')
        pivot=aug[k][k];aug[k]=[x/pivot for x in aug[k]]
        for i in range(3):
            if i!=k:
                factor=aug[i][k];aug[i]=[x-factor*y for x,y in zip(aug[i],aug[k])]
    return [row[3:] for row in aug]


def run(injection='none'):
    I.grid=0
    for name,expected in PINS.items():
        need(sha256((ROOT/'inputs'/name).read_bytes()).hexdigest()==expected,
             'input authentication: '+name)
    bg=json.loads((ROOT/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
    shell=json.loads((ROOT/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    master=json.loads((ROOT/'inputs/HIGHER_JETS_NORMAL.json').read_text())
    referee=json.loads((ROOT/'inputs/REFEREE_ACCEPTANCE_NORMAL.json').read_text())
    need(bg['actual_model']['U_coefficients']==list(map(str,U)),'actual polynomial model')
    need(bg['domain']['Z']==str(ZMAX) and bg['domain']['T']==str(TMAX),'global source report domain')
    need(shell['parameters_exact']=={'c':bg['actual_model']['c'],'delta':bg['actual_model']['delta']},
         'same retained shell branch')
    outer={k:I(*bg['background_intervals'][k]['exact']) for k in ('rho','eta','p','a')}
    r,e,p,a=[outer[k] for k in ('rho','eta','p','a')]
    rhs={'rho':-r*r*a,'eta':-r*p,'p':-r*(evaluate(UP[1],e)-4*a*p),'a':1/r+r*p*p/3}
    tube={};firstexit={}
    for name in outer:
        need(rhs[name].record()['exact']==bg['background_intervals'][name+'_z']['exact'],
             'authenticated global ODE image: '+name)
        x0=I(*shell['boxes'][name+'_b']['exact'])
        tube[name]=x0+I(0,ZMAX)*rhs[name]
        margins=[tube[name].lo-outer[name].lo,outer[name].hi-tube[name].hi]
        need(min(margins)>0,'strict whole-collar first exit: '+name)
        need(list(map(str,margins))==referee['strict_first_exit_margins'][name],
             'fresh independent referee first-exit agreement: '+name)
        firstexit[name]={'outer':outer[name].record(),'shell':x0.record(),
                        'outer_ODE_image':rhs[name].record(),'whole_collar_tube':tube[name].record(),
                        'strict_face_margins':list(map(str,margins))}
    need(131*ZMAX<9,'proper collar continuation endpoint')
    r,e,p,a=[tube[k] for k in ('rho','eta','p','a')]
    kglobal=r*(F(5,2)*a-evaluate(UP[1],e)/p)
    need(kglobal.hi<0,'negative global alpha logarithmic derivative')
    alpha_length=2*ELL if injection=='annulus_alpha_reuse' else ZMAX
    need(alpha_length==ZMAX,'alpha integral must cover the whole collar from zero')
    alpha=I(1-kglobal.mag()*alpha_length,1)
    if injection=='wrong_alpha_origin':alpha=I(alpha.lo,F(999,1000))
    need(alpha.lo>0 and alpha.lo<=1<=alpha.hi,'alpha(0)=1 lies in global enclosure')
    # Pointwise jets only need values in the accepted whole-collar tube.
    # Every subsequent interval operation is rounded outward on a rational grid.
    I.grid=10**36
    static=build_static(dict(zip(('rho','eta','p','a','alpha'),(r,e,p,a,alpha))),injection)
    formal=formal_checks(injection)
    static['C_z']=static['C'].dz()
    H,q,s2=time_jets()
    rr,aa,ay,CC=[BJ.lift(static[k],'r') for k in ('rho','a','a_y','C')]
    Cz=BJ.lift(static['C_z'],'r')
    hh,qq,ss=[BJ.lift(j,'x') for j in (H,q,s2)]
    zz=BJ.box(0)
    a_factor=1 if injection=='omit_ell_A' else ELL
    A=[[a_factor*x for x in row] for row in
       [[-3*hh,-2*aa*rr**2,zz],[-2*aa,-3*hh,zz],[zz,zz,-5*hh]]]
    B=[[zz for _ in range(3)] for _ in range(3)]
    for j,factor in enumerate((3,5,5)):B[j][j]=-ELL*factor*aa*rr
    D=[[ELL**2*x for x in row] for row in [
       [-8*ss-rr**2*(4*ay+8*aa**2),-6*aa*rr**2*hh,16*aa*rr**2],
       [zz,3+3*hh**2-8*ss,-16*hh],
       [2*aa*ss,2*hh*ss,-4*ss]]]
    Ps=6*ELL*rr*CC*qq
    P0=6*ELL**2*rr*(Cz-4*rr*aa*CC)*qq
    exact_master=list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
    simple_master=list(map(F,['1.000193','16.000488','256.002081','6665.267480']))
    need(all(x<y for x,y in zip(exact_master,simple_master)),'strict certified master C3 envelopes')
    need(master['potential_report_sha256']==PINS['ACTUAL_POTENTIAL_JET_BOUNDS.json'],
         'master source report uses this actual potential')
    need(master['source_bound_adapter']['domain']=='0<=x<=3, r>=0, x+r<=8',
         'known master domain')
    need(F(2)+2*F(3)<=8,'finite-support triangle covers rectangle nonzero support')
    source_caps={}
    for i,j in BJ.indices:
        val=F(0)
        for di in range(i+1):
            for dj in range(j+1):
                coefficient=comb(i,di)*comb(j,dj)
                val+=coefficient*(Ps.deriv(di,dj).mag()*simple_master[i+j-di-dj+1]
                                  +P0.deriv(di,dj).mag()*simple_master[i+j-di-dj])
        source_caps[f'{i},{j}']=str(val)
    matrices={'A':A,'B':B,'D':D}
    elementwise={name:[[x.deriv(0,0).mag() for x in row] for row in mat]
                 for name,mat in matrices.items()}
    rowcaps={name:[sum(row,F(0)) for row in mat] for name,mat in elementwise.items()}
    rownorm={name:max(row) for name,row in rowcaps.items()}
    mixed_rownorm={name:{f'{i},{j}':str(max(sum((entry.deriv(i,j).mag() for entry in row),F(0))
                                        for row in mat)) for i,j in BJ.indices}
                   for name,mat in matrices.items()}
    M=[[3*(elementwise['A'][i][j]+elementwise['B'][i][j]
             +F(3,2)*elementwise['D'][i][j]) for j in range(3)] for i in range(3)]
    contraction=max(map(sum,M))
    need(contraction<1,'weighted C1 Volterra matrix contraction')
    inverse=inverse3([[F(int(i==j))-M[i][j] for j in range(3)] for i in range(3)])
    need(all(x>=0 for row in inverse for x in row),'nonnegative contraction resolvent')
    source=F(source_caps['0,0'])
    channel_bounds=[3*inverse[i][0]*source for i in range(3)]
    freeze=Ps.deriv(0,0)-6*ELL/F(130)
    simple={'A_maximum_absolute_row_sum':F('0.003761'),
            'B_maximum_absolute_row_sum':F('0.00007236'),
            'D_maximum_absolute_row_sum':F('0.00000003177'),
            'source_absolute':F('0.0000007392'),
            'source_all_first_scaled_partials':F('0.00001183'),
            'source_all_second_scaled_partials':F('0.000308'),
            'P_s_minus_6ell_over130':F('0.00000000007062'),
            'P_0_absolute':F('0.00000000001280'),
            'weighted_C1_contraction':F('0.011412')}
    values={'A_maximum_absolute_row_sum':rownorm['A'],
            'B_maximum_absolute_row_sum':rownorm['B'],
            'D_maximum_absolute_row_sum':rownorm['D'],
            'source_absolute':source,
            'source_all_first_scaled_partials':max(F(source_caps[k]) for k in ('0,1','1,0')),
            'source_all_second_scaled_partials':max(F(source_caps[k]) for k in ('0,2','1,1','2,0')),
            'P_s_minus_6ell_over130':freeze.mag(),'P_0_absolute':P0.deriv(0,0).mag(),
            'weighted_C1_contraction':contraction}
    need(all(values[k]<cap for k,cap in simple.items()),'strict readable summary caps')
    return {
      'status':'PASS_WHOLE_COLLAR_GAUGE_COEFFICIENT_AND_SOURCE_CERTIFICATE',
      'input_sha256':PINS,
      'domain':{'ell':str(ELL),'t':[str(F(0)),str(TMAX)],'z':[str(F(0)),str(ZMAX)],
                'scaled_x':[0,3],'scaled_r':[0,8],
                'coverage':'Entire closed rectangle; background by first exit, master by finite-support zero extension.'},
      'arithmetic':{'method':'fractions.Fraction only; outward rational interval grid after exact first-exit check',
                    'grid_denominator':str(I.grid),'floating_arithmetic_used':False,
                    'normalization':'univariate coefficient n is D^n/n!; reported mixed derivatives are actual D_x^i D_r^j'},
      'first_exit':firstexit,
      'alpha_global':{'value':alpha.record(),'log_derivative':kglobal.record(),
                      'integral_domain':[str(F(0)),str(ZMAX)],
                      'proof':'alpha(0)=1, k<0; exp(-K Z)<=alpha<=1 and exp(-K Z)>=1-K Z.'},
      'static_z_jets':{name:[I.box(x).record() for x in series.c] for name,series in static.items()},
      'static_actual_derivative_caps':{name:[str(I.box(x).mag()*factorial(n)) for n,x in enumerate(series.c)]
                                       for name,series in static.items()},
      'temporal_t_jets':{name:[I.box(x).record() for x in j.c] for name,j in
                         {'tanh_t':H,'cosh_minus_three_halves':q,'sech_squared':s2}.items()},
      'temporal_proof':'0<=tanh(t)<=T; 1-3T^2/4<=cosh(t)^(-3/2)<=1; differential recurrences H_t=1-H^2, q_t=-3Hq/2.',
      'matrix_convention':'U_xx-U_rr=A U_x+B U_r+D U+(S,0,0), U=(xi^y/Y,xi^t/Y,angular_gradient_amplitude); raw amplitudes, no rho rescaling.',
      'matrix_mixed_derivative_intervals':{name:[[x.record() for x in row] for row in mat]
                                            for name,mat in matrices.items()},
      'matrix_elementwise_absolute_caps':{name:[[str(x) for x in row] for row in mat]
                                           for name,mat in elementwise.items()},
      'matrix_absolute_row_sums':{name:list(map(str,v)) for name,v in rowcaps.items()},
      'matrix_maximum_absolute_row_sum':{name:str(v) for name,v in rownorm.items()},
      'matrix_maximum_absolute_row_sum_by_mixed_derivative':mixed_rownorm,
      'source_definition':'S=P_s W_r+P_0 W; P_s=6ell rho C q; P_0=6ell^2 rho(C_z-4rho a C)q; C=alpha/rho^2,q=cosh(t)^(-3/2),W=w(ell x,ell r).',
      'source_profile_mixed_derivative_intervals':{'P_s':Ps.record(),'P_0':P0.record()},
      'source_strict_absolute_mixed_derivative_caps':source_caps,
      'master_scaled_order_strict_caps':list(map(str,simple_master)),
      'master_domain_extension':'Inherited smooth finite propagation gives supp(w) subset z<=2ell+t. On 0<=t<=3ell, z+t<=2ell+2t<=8ell on that support. Outside it every derivative through order3 is zero; the triangle certificate therefore covers the whole rectangle after smooth zero extension.',
      'weighted_C1_matrix':[[str(x) for x in row] for row in M],
      'weighted_C1_contraction_upper':str(contraction),
      'weighted_C1_resolvent':[[str(x) for x in row] for row in inverse],
      'conditional_weighted_C1_channel_bounds':list(map(str,channel_bounds)),
      'weighted_C1_norm':'N_i=max((2/3)||U_i||_infinity,||partial_x U_i||_infinity,||partial_r U_i||_infinity); kernel bound 3, value bound <=3N_i/2; conditional on associated D/N Volterra theorem.',
      'freezing_reference_P_s':str(6*ELL/F(130)),
      'freezing_P_s_minus_reference':freeze.record(),
      'simple_strict_caps':{k:str(v) for k,v in simple.items()},
      'formal_checks':formal,
      'scope':'Conditional on authenticated actual branch/shell box and retained smooth first-master/contact theorem. This supplies coefficient/source inequalities; the vector solution construction and its regularity are separate.'
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--inject-error',choices=['none','wrong_jet_factorial','wrong_source_sign',
                         'annulus_alpha_reuse','wrong_alpha_origin','omit_ell_A'],default='none')
    args=parser.parse_args()
    destination=args.output.resolve()
    need(destination!=ROOT and ROOT not in destination.parents,'output must be outside coefficient input tree')
    need(not destination.exists(),'refuse to overwrite output')
    report=run(args.inject_error)
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],
                     'matrix_maximum_absolute_row_sum':report['matrix_maximum_absolute_row_sum'],
                     'weighted_C1_contraction_upper':report['weighted_C1_contraction_upper'],
                     'output':str(args.output)},indent=2))


if __name__=='__main__':main()
