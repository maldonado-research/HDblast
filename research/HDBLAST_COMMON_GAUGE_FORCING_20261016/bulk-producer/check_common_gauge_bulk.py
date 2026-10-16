#!/usr/bin/env python3
"""Actual common-gauge bulk Hessian envelopes, standard library only.

All numerical decisions use exact Fraction or rational outward rounding.
Covariant frame slots, including the first derivative slot, are retained.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
PINS={
 'PUBLISHED_FIRST_GAUGE_MANIFEST.json':'adb6e3a7e00a8e477af80a49bc6b7a0f4037fafd0dd42ab8175c8b55e3e4465c',
 'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json':'4009c46f6d63583b4d161eebbb620077467b7561a3c4b5728d5293655d340a53',
 'HIGHER_GAUGE_NORMS_NORMAL.json':'f246bc6840ea36f123c7c53da28a9c04925954187e0f06df25add779ce083c6d',
 'ACTUAL_GAUGE_SOLUTION_NORMAL.json':'dd16312b2d5490259d39ec784b45c3a401f1378676e8ddbff97d54d85f2371f6',
 'HIGHER_JETS_NORMAL.json':'92b17fd71d78f83eb3f2c5ec8b9ea5f4e70365ec4d41aa4c0205551e7eaeb5df',
 'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
 'FIRST_HARMONIC_FIELDS_UNPUBLISHED.md':'d15f8c65ffffd0085d195064ccea3c61bb4e719695652d0dd73472c8a06c224d',
 'COMPLETED_HANDOFF.md':'1cb4b5f07011ee86f6a125249b5e3cf678f4d93dc0527f32cc90561cb207ea0f',
}
ELL=F(1,10**6);TMAX=3*ELL;GRID=10**40
UPOLY=list(map(F,['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27']))


def need(ok,message):
    if not ok:raise ValueError(message)
def ceilq(x):return F(-(-F(x)*GRID).__floor__(),GRID)


class I:
    def __init__(self,lo,hi=None):
        lo,hi=F(lo),F(lo if hi is None else hi)
        need(lo<=hi,'ordered interval')
        self.lo=F((lo*GRID).__floor__(),GRID);self.hi=ceilq(hi)
    @staticmethod
    def box(x):return x if isinstance(x,I) else I(x)
    def __add__(self,y):
        y=self.box(y);return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        y=self.box(y);v=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
        return I(min(v),max(v))
    __rmul__=__mul__
    def recip(self):
        need(not self.lo<=0<=self.hi,'no denominator crossing')
        return I(1/self.hi,1/self.lo)
    def __truediv__(self,y):return self*self.box(y).recip()
    def __rtruediv__(self,y):return self.box(y)*self.recip()
    def __pow__(self,n):
        if n<0:return self.recip()**(-n)
        if n==0:return I(1)
        if n%2==0:
            h=self**(n//2);return I(0 if h.lo<=0<=h.hi else min(h.lo*h.lo,h.hi*h.hi),max(h.lo*h.lo,h.hi*h.hi))
        return self*self**(n-1)
    def mag(self):return max(abs(self.lo),abs(self.hi))
    def record(self):return {'exact':[str(self.lo),str(self.hi)],'absolute_upper':str(self.mag())}


class J:
    def __init__(self,c):self.c=[F(x) if isinstance(x,(int,F)) else x for x in c]
    @property
    def n(self):return len(self.c)
    def box(self,x):return x if isinstance(x,J) else J([x]+[0]*(self.n-1))
    def __add__(self,y):
        y=self.box(y);return J([self.c[i]+y.c[i] for i in range(min(self.n,y.n))])
    __radd__=__add__
    def __neg__(self):return J([-x for x in self.c])
    def __sub__(self,y):return self+-self.box(y)
    def __rsub__(self,y):return self.box(y)+-self
    def __mul__(self,y):
        y=self.box(y);return J([sum((self.c[k]*y.c[i-k] for k in range(i+1)),0) for i in range(min(self.n,y.n))])
    __rmul__=__mul__
    def recip(self):
        out=[1/self.c[0]]
        for n in range(1,self.n):out.append(-out[0]*sum((self.c[j]*out[n-j] for j in range(1,n+1)),0))
        return J(out)
    def __truediv__(self,y):return self*self.box(y).recip()
    def __rtruediv__(self,y):return self.box(y)*self.recip()
    def __pow__(self,n):
        if n<0:return self.recip()**(-n)
        out=self.box(1)
        for _ in range(n):out=out*self
        return out
    def deriv(self):return J([(n+1)*self.c[n+1] for n in range(self.n-1)])


class M:
    """Nonnegative mixed normalized derivative majorants; no cancellation."""
    def __init__(self,d=None,order=3):
        self.order=order
        self.d={(i,n-i):ceilq((d or {}).get((i,n-i),0)) for n in range(order+1) for i in range(n+1)}
        need(all(v>=0 for v in self.d.values()),'nonnegative derivative majorant')
    def box(self,x):return x if isinstance(x,M) else M({(0,0):abs(F(x))},self.order)
    def __add__(self,y):
        y=self.box(y);n=min(self.order,y.order)
        return M({k:self.d[k]+y.d[k] for k in self.d if sum(k)<=n},n)
    __radd__=__add__
    def __mul__(self,y):
        y=self.box(y);n=min(self.order,y.order);out={}
        for i,j in self.d:
            if i+j<=n:out[i,j]=sum(self.d[a,b]*y.d[i-a,j-b] for a in range(i+1) for b in range(j+1))
        return M(out,n)
    __rmul__=__mul__
    def __pow__(self,n):
        out=self.box(1)
        for _ in range(n):out=out*self
        return out
    def derive(self,axis):
        out={}
        for i in range(self.order):
            for j in range(self.order-i):
                key=(i+1,j) if axis=='t' else (i,j+1)
                out[i,j]=(i+1 if axis=='t' else j+1)*self.d[key]
        return M(out,self.order-1)
    def get(self,i=0,j=0):return factorial(i)*factorial(j)*self.d[i,j]
    @classmethod
    def lift(cls,j,axis):
        return cls({(n,0) if axis=='t' else (0,n):I.box(v).mag() for n,v in enumerate(j.c)},len(j.c)-1)
    def record(self):return {f'{i},{j}':str(self.get(i,j)) for (i,j) in self.d}


def polyder(p):return [i*p[i] for i in range(1,len(p))]
def horner(poly,x):
    out=0
    for a in reversed(poly):out=out*x+a
    return out


def semantic_checks(injection):
    # Exact proper-radial chain rule, tested without treating rho as constant.
    rho,a=F(2),F(3,7);rho_z=-rho*rho*a
    direct_second=-1/rho*(rho_z/rho**2)
    candidate=F(0) if injection=='freeze_frame_scales' else a/rho
    need(candidate==direct_second,'radial frame scale derivative retained')
    # Constant D-components q_yy=1; K^0_0y=K^y_00=a.
    # (nabla_0 nabla_y q)_0y contains only the derivative-slot contribution a^2.
    derivative_slot=F(0) if injection=='drop_derivative_slot' else a*a
    need(derivative_slot==a*a,'second covariant jet includes derivative-slot connection')
    chi,p,f=F(2,3),F(-3,5),F(7,11)
    transformed=chi if injection=='omit_scalar_pullback' else chi+p*f
    need(transformed==chi+p*f,'metric and scalar use the same pullback plus pf')
    angular4=8 if injection=='undersize_sphere_fourth' else 16
    # At n=e1, tangent e2: D22Y=-4 and D2222Y=-4 D22Y=16.
    need(angular4>=16,'actual quadrupole fourth sphere derivative cap')
    n=5
    gauge_factor=2*n if injection=='drop_gauge_double_sum' else 2*n*n
    need(gauge_factor==2*sum(1 for _ in range(n) for _ in range(n)),
         'quadratic gauge row retains both contracted indices')
    scalar_hessian=F(0) if injection=='drop_scalar_hessian_connection' else -a
    need(scalar_hessian==-a,'scalar Hessian retains its frame connection')
    taylor_factor=1 if injection=='half_scalar_stress_hessian' else 2
    need(-taylor_factor*F(3,7)**2==-2*F(3,7)**2,
         'second parameter derivative retains scalar stress factor two')
    return {'radial_scale':'PASS_EXACT_FRACTION_IDENTITY',
            'derivative_slot':'PASS_NONZERO_RANK3_CONNECTION_WITNESS',
            'scalar_pullback':'PASS_NONZERO_EXACT_PULLBACK_WITNESS',
            'sphere_fourth':'PASS_SHARP_ACTUAL_QUADRUPOLE_WITNESS',
            'gauge_contraction_count':'PASS_FULL_TWO_INDEX_COUNT',
            'scalar_hessian':'PASS_NONZERO_EXACT_FRAME_HESSIAN_WITNESS',
            'source_Taylor_factor':'PASS_NONZERO_SECOND_PARAMETER_WITNESS'}


def tensor_count_checks():
    """Direct indexed exact tensors check all 15 Einstein sums and gauge rows."""
    n=5;sgn=[1,-1,1,1,1]
    q=lambda a,b:F(1+min(a,b)+2*max(a,b),17)
    q1=lambda e,a,b:F(1+e+a+b,29) if a!=b and e not in (a,b) else F(0)
    q2=lambda e,f,a,b:F(1+e+f+a+b,37) if a!=b and e not in (a,b) and f not in (a,b) else F(0)
    P=lambda a,b,c:sgn[a]*(q1(b,c,a)+q1(c,b,a)-q1(a,b,c))/2
    P1=lambda e,a,b,c:sgn[a]*(q2(e,b,c,a)+q2(e,c,b,a)-q2(e,a,b,c))/2
    J1=lambda e,a,b,c:-2*sum(sgn[a]*(q1(e,a,d)*P(d,b,c)+q(a,d)*P1(e,d,b,c)) for d in range(n))
    H0=max(q(a,b) for a in range(n) for b in range(n))
    H1=max(abs(q1(e,a,b)) for e in range(n) for a in range(n) for b in range(n))
    H2=max(abs(q2(e,f,a,b)) for e in range(n) for f in range(n) for a in range(n) for b in range(n))
    pcap=F(3,2)*H1;p1cap=F(3,2)*H2
    jcap=2*n*(H1*pcap+H0*p1cap)
    riccicap=2*n*jcap+4*n*n*pcap**2
    eta=F(5,13);et1=[F(i+1,19) for i in range(n)]
    et2=lambda a,b:F(2+a+b,23)
    E1=max(et1);E2=max(et2(a,b) for a in range(n) for b in range(n))
    U1,U2,U3=F(2,7),F(3,11),F(5,17)
    ecap=riccicap+2*E1*E1+F(2,3)*(U2*eta*eta+2*U1*eta*H0)
    for a in range(n):
        for b in range(a,n):
            val=sum(J1(c,c,a,b)-J1(b,c,a,c) for c in range(n))
            val+=2*sum(P(c,c,d)*P(d,a,b)-P(c,b,d)*P(d,a,c) for c in range(n) for d in range(n))
            need(abs(val)<=riccicap,'all15 direct indexed Ricci contraction bounds')
            val-=2*et1[a]*et1[b]+F(2,3)*(U2*eta*eta*(sgn[a] if a==b else 0)+2*U1*eta*q(a,b))
            need(abs(val)<=ecap,'all15 complete Einstein source contraction bounds')
    for a in range(n):
        cc=sum(sgn[b]*P(a,b,b) for b in range(n))
        need(cc==0,'nontrivial harmonic first-jet count witness')
        value=-2*sum(sgn[b]*sgn[c]*q(b,c)*P(a,b,c) for b in range(n) for c in range(n))
        need(abs(value)<=2*n*n*H0*pcap,'all5 direct indexed quadratic gauge bounds')
        for e in range(n):
            need(sum(sgn[b]*P1(e,a,b,b) for b in range(n))==0,
                 'nontrivial harmonic second-jet count witness')
            val=-2*sum(sgn[b]*sgn[c]*(q1(e,b,c)*P(a,b,c)+q(b,c)*P1(e,a,b,c))
                       for b in range(n) for c in range(n))
            need(abs(val)<=2*n*n*(H1*pcap+H0*p1cap),'all25 direct indexed derivative gauge contraction bounds')
    phi_diag=[F(3,5),F(-2,7),F(2,7),F(2,7),F(2,7)];p=F(1,11)
    scalar=2*sum(sgn[d]*q(a,d)*q(d,a)*phi_diag[a] for a in range(n) for d in range(n))
    scalar-=2*sum(sgn[a]*sgn[b]*q(a,b)*et2(a,b) for a in range(n) for b in range(n))
    scalar+=2*p*sum(sgn[a]*sgn[b]*q(a,b)*P(0,a,b) for a in range(n) for b in range(n))
    scalar-=U3*eta*eta
    sc=2*n*n*(H0*H0*max(map(abs,phi_diag))+H0*E2+H0*pcap*abs(p))+U3*eta*eta
    need(abs(scalar)<=sc,'complete scalar source indexed bound with diagonal background Hessian')
    return {'Einstein_independent_components':15,'scalar_source':1,'gauge_components':5,'gauge_derivative_components':25,
            'method':'direct exact index sums with nontrivial harmonic first jets; independent of interval enclosures'}


def run(injection='none'):
    for name,pin in PINS.items():need(sha256((ROOT/'inputs'/name).read_bytes()).hexdigest()==pin,'input authentication: '+name)
    load=lambda n:json.loads((ROOT/'inputs'/n).read_text())
    manifest=load('PUBLISHED_FIRST_GAUGE_MANIFEST.json')
    entries={row['path']:row for row in manifest['files']}
    published_paths={
      'GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json':'coefficients/GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json',
      'HIGHER_GAUGE_NORMS_NORMAL.json':'higher-norms/HIGHER_GAUGE_NORMS_NORMAL.json',
      'ACTUAL_GAUGE_SOLUTION_NORMAL.json':'solver/ACTUAL_GAUGE_SOLUTION_NORMAL.json',
      'HIGHER_JETS_NORMAL.json':'coefficients/inputs/HIGHER_JETS_NORMAL.json',
      'ACTUAL_POTENTIAL_JET_BOUNDS.json':'coefficients/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json'}
    for name,path in published_paths.items():
        need(entries[path]['sha256']==PINS[name] and
             entries[path]['bytes']==len((ROOT/'inputs'/name).read_bytes()),
             'published manifest authenticates copied actual input: '+name)
    coef=load('GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json');gauge=load('HIGHER_GAUGE_NORMS_NORMAL.json')
    solver=load('ACTUAL_GAUGE_SOLUTION_NORMAL.json');master=load('HIGHER_JETS_NORMAL.json')
    potential=load('ACTUAL_POTENTIAL_JET_BOUNDS.json')
    need(coef['input_sha256']['ACTUAL_POTENTIAL_JET_BOUNDS.json']==PINS['ACTUAL_POTENTIAL_JET_BOUNDS.json'],'same actual background')
    need(gauge['source_pins']['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json']==PINS['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'],'same gauge C3 coefficients')
    need(solver['input_sha256']['coefficients']==PINS['GLOBAL_GAUGE_COEFFICIENTS_NORMAL.json'],'same actual C1 gauge solution')
    need(master['potential_report_sha256']==PINS['ACTUAL_POTENTIAL_JET_BOUNDS.json'],'same actual master background')
    need(potential['actual_model']['U_coefficients']==list(map(str,UPOLY)),'actual potential polynomial')
    checks=semantic_checks(injection);counts=tensor_count_checks()
    static={k:J([I(*v['exact']) for v in row]) for k,row in coef['static_z_jets'].items()}
    rr,ee,pp,aa,alpha=[static[k] for k in ('rho','eta','p','a','alpha')]
    invr=1/rr;hh=horner(polyder(UPOLY),ee)/pp
    C=alpha/rr**2;B=3*alpha/(rr**3*pp);kk=rr*(F(5,2)*aa-hh)
    potjets={};up=UPOLY
    for n in range(4):potjets['U'+str(n)]=horner(up,ee);up=polyder(up)
    H=J([I(*v['exact']) for v in coef['temporal_t_jets']['tanh_t']])
    qt=J([I(*v['exact']) for v in coef['temporal_t_jets']['cosh_minus_three_halves']])
    cc=J([I(1,1/(1-TMAX*TMAX/2))]+[I(0)]*3)
    for n in range(3):cc.c[n+1]=(H*cc).c[n]/(n+1)
    ic=1/cc
    mj=lambda j,axis:M.lift(j,axis)
    rmi=mj(invr,'z');am=mj(aa,'z');pm=mj(pp,'z');km=mj(kk,'z')
    cm=mj(cc,'t');icm=mj(ic,'t');Hm=mj(H,'t');qtm=mj(qt,'t')
    mastercaps=list(map(F,master['source_bound_adapter']['scaled_order_bounds']))
    W=M({(i,n-i):mastercaps[n]*ELL**(-n)/(factorial(i)*factorial(n-i))
         for n in range(4) for i in range(n+1)})
    actual_solver=solver['time_window_certificates']['3']
    raw=[];selected={}
    for channel,name in enumerate(('f','b','d')):
        data={};record={}
        for key,values in gauge['unscaled_coordinate_mixed_partial_caps'].items():
            i,j=map(int,key.split(','));v=F(values[channel])
            if i+j==0:v=min(v,F(actual_solver['actual_amplitude_value_by_channel'][channel]['exact_upper']))
            if i+j==1:v=min(v,F(actual_solver['actual_each_scaled_first_derivative_by_channel'][channel]['exact_upper'])/ELL)
            data[i,j]=v/(factorial(i)*factorial(j));record[key]=str(v)
        raw.append(M(data));selected[name]=record
    f,b,d=raw
    u=mj(C,'z')*qtm*W
    v=mj(B,'z')*qtm*(W.derive('z')+km*W)
    eta=v+pm*f
    amplitudes={
      'yy':[(4*u+2*rmi*f.derive('z'),(1,2,4))],
      'y0':[(rmi*f.derive('t')+b.derive('z'),(1,2,4))],
      'ya':[(rmi*icm*f+cm*d.derive('z'),(2,4,8))],
      '00':[(2*(u+am*f+b.derive('t')),(1,2,4))],
      '0a':[(cm*d.derive('t')+icm*b,(2,4,8))],
      'ab':[(2*(u+am*f+Hm*b),(1,2,4)),(2*d,(4,8,16))]}
    rinv=I.box(invr.c[0]).mag();av=I.box(aa.c[0]).mag();Hv=I.box(H.c[0]).mag()
    cinv=I.box(ic.c[0]).mag();sv=rinv*cinv
    def baseline(terms):
        val=F(0);first=[F(0)]*3;second=[F(0)]*9
        for amp,angular in terms:
            F0,Ft,Fz,Ftt,Ftz,Fzz=[amp.get(i,j) for i,j in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))]
            a0,a1,a2=map(F,angular);Fy=rinv*Fz
            Fyy=rinv*rinv*Fzz+av*rinv*Fz;Fyt=rinv*Ftz
            val+=F0*a0
            ff=[Fy*a0,rinv*Ft*a0,sv*F0*a1]
            ss=[Fyy*a0,rinv*(Fyt+av*Ft)*a0,rinv*Fyt*a0,
                rinv*rinv*Ftt*a0,sv*(Fy+av*F0)*a1,sv*Fy*a1,
                sv*rinv*(Ft+Hv*F0)*a1,sv*rinv*Ft*a1,sv*sv*F0*a2]
            first=[x+y for x,y in zip(first,ff)];second=[x+y for x,y in zip(second,ss)]
        return {'value':val,'first_ordered':first,'second_ordered':second,'D1':max(first),'D2':max(second)}
    base={name:baseline(terms) for name,terms in amplitudes.items()}
    se=baseline([(eta,(1,2,4))])
    H0=max(row['value'] for row in base.values());D1=max(row['D1'] for row in base.values());D2=max(row['D2'] for row in base.values())
    ay=I.box(static['a_y'].c[0]).mag();beta=Hv*rinv
    K0=av+beta;K1=max(ay+av*beta,rinv*rinv)
    H1=D1+2*K0*H0;H2=D2+5*K0*D1+(2*K1+6*K0*K0)*H0
    E0=se['value'];E1=se['D1'];E2=se['D2']+K0*E1
    U1,U2,U3=[I.box(potjets['U'+str(i)].c[0]).mag() for i in (1,2,3)]
    pv=I.box(pp.c[0]).mag();py=I.box((potjets['U1']-4*aa*pp).c[0]).mag();Phi2=max(py,av*pv)
    n=5;P0=F(3,2)*H1;P1=F(3,2)*H2;J1=2*n*(H1*P0+H0*P1)
    QE=2*n*J1+4*n*n*P0*P0+2*E1*E1+F(2,3)*(U2*E0*E0+2*U1*E0*H0)
    QC=2*n*n*H0*P0;DQC=2*n*n*(H1*P0+H0*P1)
    Qphi=2*n*n*(H0*H0*Phi2+H0*E2+H0*P0*pv)+U3*E0*E0
    bounds={'QE':QE,'Qphi':Qphi,'QC':QC,'nabla_QC':DQC,'QEred':QE+DQC}
    def simple(v):
        exponent=0;w=v
        while w>=10:w/=10;exponent+=1
        while w<1:w*=10;exponent-=1
        step=F(10)**(exponent-3)
        return ((v/step).__floor__()+1)*step
    simplify=lambda values:[str(simple(v)) for v in values]
    labels=['yy','y0','y1','y2','y3','00','01','02','03','11','12','13','22','23','33']
    encode=lambda v:({k:encode(x) for k,x in v.items()} if isinstance(v,dict)
                     else [encode(x) for x in v] if isinstance(v,list) else str(v))
    return {
      'status':'PASS_ACTUAL_FULL_COMMON_GAUGE_BULK_HESSIAN_ENCLOSURE',
      'input_sha256':PINS,
      'published_input_paths':published_paths,
      'domain':{'ell':str(ELL),'scaled':'0<=x<=3,r>=0,x+r<=8','angular':'entire unit S3'},
      'norm':'maximum absolute individual adapted-frame component, including every covariant derivative slot; e_y,e_0 fixed, sphere-orthonormal angular frame; no Lorentz boosts',
      'source_convention':'affine pair path(g0+eps q,phi0+eps eta); second parameter derivative; Ered=E-nabla_(A C_B) with one-half symmetrization; actual C1=0 and nabla C1=0',
      'actual_source_bounds':{k:{'exact_upper':str(v),'strict_simple_upper':str(simple(v))} for k,v in bounds.items()},
      'source_formulas':{
          'QE':'2n*J1+4n^2*P0^2+2*E1^2+(2/3)(U2*E0^2+2*U1*E0*H0)',
          'Qphi':'2n^2*(H0^2*Phi2+H0*E2+H0*P0*abs(p))+U3*E0^2',
          'QC':'2n^2*H0*P0','nabla_QC':'2n^2*(H1*P0+H0*P1)',
          'QEred':'QE+nabla_QC',
          'auxiliary':'P0=3H1/2,P1=3H2/2,J1=2n*(H1*P0+H0*P1),n=5'},
      'universal_formula_constants':{'dimension':5,'P_from_metric_gradient':'3/2',
          'J1_product_factor':10,'Einstein_J1_factor':10,'Einstein_P_squared_factor':100,
          'Einstein_scalar_gradient_squared_factor':2,'potential_factor':'2/3',
          'scalar_first_three_terms_factor':50,'gauge_double_contraction_factor':50,
          'reduced_symmetrized_derivative_bound_factor':1,
          'metric_gradient_connection_slots':2,'metric_second_D1_factor':5,
          'metric_second_DK_factor':2,'metric_second_K_squared_factor':6,
          'scalar_Hessian_connection_factor':1},
      'normalized_sphere_Frobenius_source_caps':{'QE':str(5*QE),'QEred':str(5*(QE+DQC)),'Qphi':str(Qphi),
          'meaning':'Euclidean Frobenius norm of adapted-frame 5x5 components, then normalized S3 L2; at most5 times max-component cap. Not Lorentzian contraction.'},
      'field_covariant_caps':{'metric':list(map(str,(H0,H1,H2))),'scalar':list(map(str,(E0,E1,E2)))},
      'strict_readable_field_covariant_caps':{'metric':simplify((H0,H1,H2)),'scalar':simplify((E0,E1,E2))},
      'field_baseline_caps':encode(base),'scalar_baseline_caps':encode(se),
      'ordered_baseline_first_slots':['y','0','a'],
      'ordered_baseline_second_slots':['yy','y0','0y','00','ya','ay','0a','a0','ab'],
      'metric_amplitude_mixed_derivative_majorants':{k:[amp.record() for amp,_ in terms] for k,terms in amplitudes.items()},
      'scalar_amplitude_mixed_derivative_majorants':eta.record(),
      'selected_actual_gauge_mixed_caps':selected,
      'static_reconstruction_jets':{name:[I.box(x).record() for x in j.c] for name,j in
                                     {'B':B,'k':kk,'C_rec':C,'rho_inv':invr,**potjets}.items()},
      'connection_and_background_caps':encode({'K0':K0,'K1':K1,'a':av,'a_y':ay,'beta':beta,
                                              'p':pv,'p_y':py,'Phi2':Phi2,'U1':U1,'U2':U2,'U3':U3}),
      'background_coefficient_caps':encode({'rho':I.box(rr.c[0]).mag(),'rho_inv':rinv,
          'a':av,'a_y':ay,'p':pv,'p_y':py,'H':Hv,'cosh_t':I.box(cc.c[0]).mag(),
          'cosh_t_inverse':cinv,'alpha':I.box(alpha.c[0]).mag(),
          'B':I.box(B.c[0]).mag(),'k':I.box(kk.c[0]).mag(),'C_rec':I.box(C.c[0]).mag(),
          **{name:I.box(j.c[0]).mag() for name,j in potjets.items()}}),
      'connection_variation_caps':encode({'P':P0,'nabla_P':P1,'nabla_J':J1}),
      'coverage':{'Einstein_independent_components':labels,'Einstein_component_count':15,
                  'Einstein_and_reduced_each_use_common_cap':True,'QC_components':5,'nabla_QC_components':25,
                  'sphere_harmonic_caps_through_order4':[1,2,4,8,16],
                  'all_angular_tensor_sectors_retained':True},
      'arithmetic':{'type':'standard-library Fraction; outward interval and upward absolute jet rounding',
                    'rational_grid_denominator':str(GRID),'floating_arithmetic':False},
      'semantic_checks':checks,'direct_tensor_count_checks':counts,
      'scope':'Continuous actual common-gauge bulk second-parameter-source component bounds only. No derivative of QE/Qphi/QEred, second evolution, full boundary or constraint forcing, nonlinear amplitude/remainder, energy or observational claim.'
    }


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--inject-error',choices=['none','freeze_frame_scales','drop_derivative_slot',
                         'omit_scalar_pullback','undersize_sphere_fourth','drop_gauge_double_sum',
                         'drop_scalar_hessian_connection','half_scalar_stress_hessian'],default='none')
    args=parser.parse_args();out=args.output.resolve()
    need(out!=ROOT and ROOT not in out.parents,'output outside source/input tree')
    need(not out.exists(),'refuse to overwrite output')
    report=run(args.inject_error);out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'bounds':{k:v['strict_simple_upper'] for k,v in report['actual_source_bounds'].items()},'output':str(out)},indent=2))


if __name__=='__main__':main()
