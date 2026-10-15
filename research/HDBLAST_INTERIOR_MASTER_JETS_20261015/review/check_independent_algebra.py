#!/usr/bin/env python3
"""Independent standard-library interval/jet and reflected-kernel checks.

No scientific implementation is imported. This uses Taylor coefficient
propagation of the background ODE rather than the source's handwritten Wzz
formula. Further acceptance of actual solution regularity is analytical.
"""
from fractions import Fraction as F
from csv import DictReader
from itertools import product
from math import comb
from pathlib import Path
import argparse
import hashlib
import json

INPUT_PINS = {
    'BACKGROUND_DERIVATIVE_BOUNDS.json':'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
    'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md':'0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
    'OCT11_JACOBI_ENDPOINT_MONOTONICITY.md':'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8',
    'OCT11_RESCALED_BOUND_CERTIFICATE.json':'22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'BACKGROUND_DERIVATIVE_AUDIT.md':'c85130916eb36a4a81c512cf594db8f9d4327dba926233aef126da95d45b4fff',
}


def need(test, why):
    if not test:
        raise ValueError(why)


class Box:
    def __init__(self, low, high=None):
        self.low, self.high = F(low), F(low if high is None else high)
        need(self.low <= self.high, 'box order')

    @classmethod
    def asbox(cls, x):
        return x if isinstance(x, cls) else cls(x)

    def __add__(self, x):
        x = self.asbox(x)
        return Box(self.low+x.low, self.high+x.high)
    __radd__ = __add__

    def __neg__(self):
        return Box(-self.high, -self.low)

    def __sub__(self, x):
        return self+-self.asbox(x)

    def __rsub__(self, x):
        return self.asbox(x)+-self

    def __mul__(self, x):
        x = self.asbox(x)
        products = [a*b for a in (self.low,self.high) for b in (x.low,x.high)]
        return Box(min(products),max(products))
    __rmul__ = __mul__

    def inverse(self):
        need(self.low*self.high>0, 'zero denominator')
        return Box(1/self.high,1/self.low)

    def __truediv__(self, x):
        return self*self.asbox(x).inverse()

    def __rtruediv__(self, x):
        return self.asbox(x)*self.inverse()

    def absmax(self):
        return max(abs(self.low),abs(self.high))

    def record(self):
        return [str(self.low),str(self.high)]


class Jet:
    """Three Taylor coefficients, including factorial denominators."""
    def __init__(self, coeff):
        self.c = [Box.asbox(x) for x in coeff]
        need(len(self.c)==3,'jet degree')

    @classmethod
    def lift(cls, x):
        return x if isinstance(x,cls) else cls([x,0,0])

    def __add__(self,x):
        x=self.lift(x)
        return Jet([a+b for a,b in zip(self.c,x.c)])
    __radd__=__add__

    def __neg__(self):
        return Jet([-a for a in self.c])

    def __sub__(self,x):
        return self+-self.lift(x)

    def __rsub__(self,x):
        return self.lift(x)+-self

    def __mul__(self,x):
        x=self.lift(x)
        return Jet([sum((self.c[j]*x.c[i-j] for j in range(i+1)),Box(0)) for i in range(3)])
    __rmul__=__mul__

    def inverse(self):
        coeff=[self.c[0].inverse()]
        for i in (1,2):
            coeff.append(-sum((self.c[j]*coeff[i-j] for j in range(1,i+1)),Box(0))/self.c[0])
        return Jet(coeff)

    def __truediv__(self,x):
        return self*self.lift(x).inverse()

    def __rtruediv__(self,x):
        return self.lift(x)*self.inverse()


def polynomial(coeff,x):
    out=Jet.lift(0)
    for a in reversed(coeff):
        out=out*x+a
    return out


def derivative(coeff):
    return [i*coeff[i] for i in range(1,len(coeff))]


def ode(state,potential):
    rho,eta,p,a=state
    U1=polynomial(derivative(potential),eta)
    return [-rho*rho*a,-rho*p,-rho*(U1-4*a*p),1/rho+rho*p*p/3]


def state_jets(boxes,potential):
    values=[Jet([b,0,0]) for b in boxes]
    for degree in (1,2):
        rhs=ode(values,potential)
        for j in range(4):
            values[j].c[degree]=rhs[j].c[degree-1]/degree
    return values


def independent_potential(root):
    report=json.loads((root/'potential/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
    bg=json.loads((root/'potential/inputs/BACKGROUND_DERIVATIVE_BOUNDS.json').read_text())
    for name,pin in INPUT_PINS.items():
        fresh=(root/'potential/inputs'/name).read_bytes()
        need(hashlib.sha256(fresh).hexdigest()==pin,'independent retained input pin '+name)
    # Derive U directly from the explicit shifted cubic superpotential.
    superpot=[F(1,3),F(0),F(1),F(1,3)]
    def conv(a,b):
        return [sum((a[j]*b[i-j] for j in range(len(a)) if 0<=i-j<len(b)),F(0)) for i in range(len(a)+len(b)-1)]
    square=conv(superpot,superpot)
    dsquare=conv(derivative(superpot),derivative(superpot))
    potential=[(dsquare[i]/2 if i<len(dsquare) else 0)-F(2,3)*square[i] for i in range(len(square))]
    need(potential==list(map(F,report['actual_model']['U_coefficients'])),'independent potential coefficients')
    boxes=[Box(*bg['first_exit']['component_margins'][k]['outer']['exact']) for k in ('rho','eta','p','a')]
    rho,eta,p,a=state_jets(boxes,potential)
    h=polynomial(derivative(potential),eta)/p
    U2=polynomial(derivative(derivative(potential)),eta)
    Q=2*h*h-10*a*h+F(35,4)*a*a-U2-p*p/6
    V=rho*rho*Q
    # Uncancelled interval automatic differentiation is intentionally broader.
    global_bounds={name:V.c[i]*factor for i,(name,factor) in enumerate((('W_spatial',1),('W_z',1),('W_zz',2)))}
    need(global_bounds['W_spatial'].absmax()<F(102000),'automatic potential bound')
    # The unsimplified AD box is wider than the source's cancellation-aware
    # expression. Cartesian interval subdivisions prove a separate whole-box
    # bound; every subrectangle is covered, with no numerical trajectory.
    independent_bound=None
    for panels in (2,4):
        ranges=[]
        for indices in product(range(panels),repeat=4):
            small=[Box(b.low+(b.high-b.low)*F(j,panels),b.low+(b.high-b.low)*F(j+1,panels)) for b,j in zip(boxes,indices)]
            rr,ee,pp,aa=state_jets(small,potential)
            hh=polynomial(derivative(potential),ee)/pp
            qq=2*hh*hh-10*aa*hh+F(35,4)*aa*aa-polynomial(derivative(derivative(potential)),ee)-pp*pp/6
            vv=rr*rr*qq
            ranges.append(vv.c[2]*2)
        independent_bound=Box(min(b.low for b in ranges),max(b.high for b in ranges))
        if independent_bound.absmax()<F(30000000000):break
    need(independent_bound.absmax()<F(30000000000),'subdivided independent ODE Taylor Wzz bound')
    # At exact rational states AD and the separate reported handwritten formula
    # have identical scalar jets. This checks independent algebra ordering.
    witnesses=0
    for n in range(9):
        point=[Box(b.low+(b.high-b.low)*F(n,8)) for b in boxes]
        rj,ej,pj,aj=state_jets(point,potential)
        hh=polynomial(derivative(potential),ej)/pj
        qq=2*hh*hh-10*aj*hh+F(35,4)*aj*aj-polynomial(derivative(derivative(potential)),ej)-pj*pj/6
        vv=rj*rj*qq
        need(all(z.low==z.high for z in vv.c),'point jets exact')
        # Every exact witness lies in the claimed true enclosure.
        for i,name in enumerate(('W','W_z','W_zz')):
            lo,hi=map(F,report['potential_intervals'][name]['exact'])
            value=vv.c[i].low*(2 if i==2 else 1)
            need(lo<=value<=hi,'reported W jet contains independent ODE witness '+name)
            witnesses+=1
    return {'copied_inputs_authenticated':len(INPUT_PINS),
            'independent_polynomial_degree':len(potential)-1,
            'ode_taylor_point_jet_witnesses':witnesses,
            'uncancelled_global_automatic_jet_boxes':{k:v.record() for k,v in global_bounds.items()},
            'uncancelled_global_Wzz_bound_below_claim':global_bounds['W_zz'].absmax()<F(30000000000),
            'rigorous_cartesian_subrectangles':panels**4,
            'independent_subdivided_Wzz_box':independent_bound.record(),
            'independent_subdivided_Wzz_bound_below_claim':True}


def integrate_poly(m,n,A,B,low,high):
    return sum((comb(n,j)*A**(n-j)*B**j*(high**(m+j+1)-low**(m+j+1))/F(m+j+1) for j in range(n+1)),F(0))


def kernel(m,n,x,r):
    # Integral of the spatial primitive of sign(q)*s^m*|q|^n.
    split=max(F(0),min(x,x-r))
    top=integrate_poly(m,n+1,r+x,-1,0,x)
    below=integrate_poly(m,n+1,x-r,-1,0,split)
    above=integrate_poly(m,n+1,r-x,1,split,x)
    return (top-below-above)/F(2*(n+1))


def finite_derivative(function,degree,order,step):
    # Solve exact moment equations for one-sided interpolation derivatives.
    count=degree+1
    matrix=[[F(j)**k for j in range(count)]+[F(1) if k==order else F(0)] for k in range(count)]
    # Desired moment is order!; use falling-product after solve.
    fac=1
    for j in range(2,order+1):fac*=j
    for i in range(count):
        pivot=next(k for k in range(i,count) if matrix[k][i])
        matrix[i],matrix[pivot]=matrix[pivot],matrix[i]
        scale=matrix[i][i];matrix[i]=[a/scale for a in matrix[i]]
        for k in range(count):
            if k!=i:
                scale=matrix[k][i];matrix[k]=[a-scale*b for a,b in zip(matrix[k],matrix[i])]
    return sum((matrix[i][-1]*function(i*step) for i in range(count)),F(0))*fac/step**order


def reflected_ray(m,n,x,r,which):
    split=max(F(0),min(x,x-r))
    # g=f_r is even, h=f_rr is odd. Explicit shell ray terms retained.
    d=1 if which in ('rr','xr') else 2
    if n<d:upper=lower=F(0)
    else:
        factor=n if d==1 else n*(n-1)
        upper=factor*integrate_poly(m,n-d,r+x,-1,0,x)
        lower=factor*((-1 if d==2 else 1)*integrate_poly(m,n-d,x-r,-1,0,split)+integrate_poly(m,n-d,r-x,1,split,x))
    shell=(x-r)**m if x>r and n==0 else F(0)
    shellprime=m*(x-r)**(m-1) if x>r and n==0 and m else F(0)
    if which=='rr':return (upper-lower)/2-shell
    if which=='xr':return (upper+lower)/2+shell
    if which=='rrr':return (upper-lower)/2+shellprime
    if which=='xrr':return (upper+lower)/2-shellprime
    raise ValueError(which)


def reflection_checks():
    count=0
    for m in range(4):
        for n in range(4):
            degree=m+n+2
            for x,r in ((F(2),F(1)),(F(1),F(2)),(F(2),F(0))):
                for which,dx,dr in (('rr',0,2),('xr',1,1),('rrr',0,3),('xrr',1,2)):
                    step=F(1,30*(degree+1))
                    if dx:
                        value=finite_derivative(lambda a:finite_derivative(lambda b:kernel(m,n,x+a,r+b),degree,dr,step),degree,dx,step)
                    else:
                        value=finite_derivative(lambda b:kernel(m,n,x,r+b),degree,dr,step)
                    need(value==reflected_ray(m,n,x,r,which),'independent reflected derivative '+which)
                    count+=1
    # Contact characteristic: choose m>=2 for C3 corner flatness, compare both
    # one-sided third spatial derivatives at a characteristic point.
    joins=0
    for m in (2,3,4):
        for n in (0,1,2,3):
            degree=m+n+2
            values=[finite_derivative(lambda d:kernel(m,n,F(1),F(1)+d),degree,3,sgn*F(1,30*(degree+1))) for sgn in (-1,1)]
            need(values[0]==values[1],'C3 characteristic join for flat forcing corner')
            joins+=1
    # Nonflat constant forcing must fail C2 joining, and x forcing must fail C3.
    bad2=[finite_derivative(lambda d:kernel(0,0,F(1),F(1)+d),2,2,sgn*F(1,100)) for sgn in (-1,1)]
    bad3=[finite_derivative(lambda d:kernel(1,0,F(1),F(1)+d),3,3,sgn*F(1,100)) for sgn in (-1,1)]
    need(bad2[0]!=bad2[1] and bad3[0]!=bad3[1],'detect incompatible forcing corner jumps')
    return {'C2_C3_reflected_exact_witnesses':count,'C3_characteristic_join_witnesses':joins,'incompatible_corner_counterexamples':2}


def incoming_checks():
    # Manufacture U=x^m r^n, define f=Urr-Uxx, and integrate the exact
    # incoming characteristic. No source implementation is imported.
    def monomial(m,n,x,r):
        return x**m*r**n
    ell,gamma=F(1,10),F(3)
    count=0
    for m,n in product(range(5),repeat=2):
        for x in (F(0),F(1,2),F(2)):
            I=F(0);Ir=F(0);fend=F(0)
            if n>=2:
                I+=n*(n-1)*integrate_poly(m,n-2,x,-1,0,x)
                if n>=3:Ir+=n*(n-1)*(n-2)*integrate_poly(m,n-3,x,-1,0,x)
                fend+=n*(n-1)*monomial(m,n-2,x,F(0))
            if m>=2:
                I-=m*(m-1)*integrate_poly(m-2,n,x,-1,0,x)
                if n:Ir-=m*(m-1)*n*integrate_poly(m-2,n-1,x,-1,0,x)
                fend-=m*(m-1)*monomial(m-2,n,x,F(0))
            y=monomial(m,n,x,F(0))
            yp=m*monomial(m-1,n,x,F(0)) if m else F(0)
            ypp=m*(m-1)*monomial(m-2,n,x,F(0)) if m>=2 else F(0)
            J0=(x**n if m==1 else F(0))+(n*x**(n-1) if m==0 and n else F(0))
            J0p=(n*x**(n-1) if m==1 and n else F(0))+(n*(n-1)*x**(n-2) if m==0 and n>=2 else F(0))
            H=(n*monomial(m,n-1,x,F(0)) if n else F(0))-ell*gamma*y
            Hp=(m*n*monomial(m-1,n-1,x,F(0)) if m and n else F(0))-ell*gamma*yp
            need(H==J0-yp-ell*gamma*y-I,'incoming H identity')
            need(Hp==J0p-ypp-ell*gamma*yp-fend-Ir,'incoming Hprime endpoint identity')
            count+=2
    # Polynomial positivity after u=s-4 provides the global bump envelopes.
    # Q''=e^(4-s)s²(s²-6s+6); Q'''=v e^(4-s)s³(s³-10s²+24s-12).
    positivity_certificates={
        's2_minus_Qsecond_factor':[18,6],
        's2_plus_Qsecond_factor':[14,10,2],
        's3_minus_Qthird_factor':[76,56,10],
        's3_plus_Qthird_factor':[52,40,14,2],
    }
    need(all(c>=0 for v in positivity_certificates.values() for c in v),'positive-coefficient bump certificates')
    # Differentiate e^(4-s)s^k: maximum at s=max(4,k), proving Q1<=16,
    # Q2<=256. For k=6 use e²>sum_(j=0)^4 2^j/j!=7.
    factorial,total=1,F(1)
    for j in range(1,5):
        factorial*=j;total+=F(2**j,factorial)
    need(total==7 and F(6**6,7)<6666,'strict third-bump rational envelope')
    return {'incoming_H_Hprime_exact_witnesses':count,
            'bump_positive_polynomial_certificates':positivity_certificates,
            'bump_derivative_bounds':['1','16','256',str(F(6**6,7))],
            'bump_scope':'Analytic envelope over s>=4 and flat zero extension; no sampled exponential claims.'}


def independent_jet_bounds(root,inputs_root):
    pins={
        'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv':'2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7',
        'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json':'d2ffbcec5413580ff5703da58b062a9b83024431748fa40726555c09c2a968d8',
        'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
        'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json':'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
        'inputs/reference/PHYSICAL_SHELL_RESPONSE.md':'a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b',
    }
    for name,pin in pins.items():
        need(hashlib.sha256((inputs_root/name).read_bytes()).hexdigest()==pin,'independent waveform input '+name)
    old=json.loads((inputs_root/'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json').read_text())
    jets=json.loads((root/'jets/HIGHER_JETS_NORMAL.json').read_text())
    ell=F(1,10**6);X=F(3)
    Hmax=Gmax=Lmax=F(0)
    # G=int H and L=int G; propagate integrated endpoint boxes, avoiding the
    # primary program's p0*x-p1 representation.
    G=[F(0),F(0)];L=[F(0),F(0)];last=F(0);cells=0
    with (inputs_root/'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv').open() as stream:
        for row in DictReader(stream):
            a=F(row['left_grid'])/F(old['cell_data']['time_grid_denominator'])
            b=F(row['right_grid'])/F(old['cell_data']['time_grid_denominator'])
            need(a==last and a<b,'independent complete H partition')
            H=[F(row[k])/F(old['cell_data']['value_denominator']) for k in ('H_lo','H_hi')]
            need(H[0]<=H[1],'independent H box order')
            dt=b-a
            for k in (0,1):
                times=[F(0),dt]
                if H[k] and 0<-G[k]/H[k]<dt:times.append(-G[k]/H[k])
                Lmax=max(Lmax,*[abs(L[k]+G[k]*t+H[k]*t*t/2) for t in times])
                L[k]+=G[k]*dt+H[k]*dt*dt/2
                Gmax=max(Gmax,abs(G[k]),abs(G[k]+H[k]*dt))
                G[k]+=H[k]*dt
            Hmax=max(Hmax,*map(abs,H));last=b;cells+=1
    need(last==3 and cells==4353,'independent H complete coverage')
    gain=ell/min(map(F,old['shell_parameter_intervals']['mu']))
    gamma=max(map(F,old['shell_parameter_intervals']['gamma']))
    A=F(5,2)*ell*ell;Ap=F(35,2)*ell**3*X*ell
    y=gain*Lmax/(1-A*X*X/2)
    yp=gain*Gmax+X*A*y
    ypp=gain*Hmax+A*y
    M=F(102000);Mz=F(30000000);Mt=F(21,400000);Mzz=F(30000000000)
    U=(1+X*ell+y)/(1-M*(X*ell)**2/2)
    first=16+ell+yp+X*ell*ell*M*U
    f=ell*ell*M*U;fr=ell**3*Mz*U+ell*ell*M*first
    fx=ell**3*Mt*U+ell*ell*M*first
    f0=ell*ell*M*y;fx0=ell**3*Mt*y+ell*ell*M*yp
    Hp=256+ell*16+ell*ell+ypp+ell*gamma*yp+f0+X*fr
    yppp=gain*Hp+Ap*y+A*yp
    free2=256+ell*(16+ell)+ypp
    rr=free2+X*fr+f0;second=rr+f
    frr=ell**4*Mzz*U+2*ell**3*Mz*first+ell*ell*M*rr
    third=F(46656,7)+ell*(256+32*ell+ell*ell)+yppp+X*frr+fx0+max(fr,fx)
    compare={'y':y,'yprime':yp,'ysecond':ypp,'Hprime':Hp,'ythird':yppp,
             'U':U,'first':first}
    for name,value in compare.items():
        field=jets['exact']['boundary'] if name in ('y','yprime','ysecond') else jets['exact']['interior']
        need(F(field[name])==value,'independent exact jet scalar '+name)
    adapter=[U,first,second,third]
    need(adapter==list(map(F,jets['source_bound_adapter']['scaled_order_bounds'])),'independent all-order source adapter')
    return {'whole_H_cells_reconstructed':cells,'retained_input_pins':len(pins),
            'independent_scalar_equalities':len(compare),'independent_order_norm_equalities':4,
            'exact_scaled_order_bounds':list(map(str,adapter))}


def independent_bulk_bounds(root,adapter):
    bg=json.loads((root/'potential/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
    source=json.loads((root/'bulk/ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS.json').read_text())
    def magnitude(section,name):
        return max(map(abs,map(F,bg[section][name]['exact'])))
    R=min(map(F,bg['background_intervals']['rho']['exact']))
    A=magnitude('background_intervals','a')
    Ay=magnitude('background_intervals','a_z')/R
    P=[magnitude('background_intervals','U_'+str(j)) for j in (1,2,3)]
    C,C1,C2=[magnitude('reconstruction_coefficients','u_w_coefficient'+suffix) for suffix in ('','_z','_zz')]
    B,B1=[magnitude('reconstruction_coefficients','v_wz_coefficient'+suffix) for suffix in ('','_z')]
    k,k1=[magnitude('reconstruction_coefficients',name) for name in ('alpha_log_z','alpha_log_zz')]
    ell=F(1,10**6);T=3*ell
    need(T*T<F(4,5),'clock second logarithmic quotient bound')
    W=[F(x)/ell**j for j,x in enumerate(adapter)]
    for j in range(4):need(W[j]*ell**j==F(adapter[j]),'independent physical derivative conversion')
    q1=3*T/2;q2=F(3,2)
    u=C*W[0]
    ut=C*W[1]+q1*u
    utt=C*W[2]+2*q1*C*W[1]+q2*u
    uz=C1*W[0]+C*W[1]
    uy=uz/R
    uyy=(C2*W[0]+2*C1*W[1]+C*W[2])/R**2+A*uy
    uyt=(C1*W[1]+C*W[2]+q1*uz)/R
    v=B*(W[1]+k*W[0])
    vt=B*(W[2]+k*W[1])+q1*v
    vtt=B*(W[3]+k*W[2])+2*q1*B*(W[2]+k*W[1])+q2*v
    vy=(B1*(W[1]+k*W[0])+B*(W[2]+k*W[1]+k1*W[0]))/R
    dpsi=max(ut,2*u)/R;dchi=max(vt,2*v)/R
    dpsiy=max(uyt,2*uy)/R
    hess=max(utt,2*ut+2*T*u,4*u+T*ut)/R**2
    wave=(utt+3*T*ut+8*u)/R**2
    contraction=(ut**2+4*u**2)/R**2
    radial=8*uy**2+104*A*u*uy+8*u*uyy+48*(Ay+4*A**2)*u**2
    Eyy=sum([8*uy**2,32*A*u*uy,16*u*uyy,16*contraction,8*u*wave,
             2*vy**2,F(2,3)*P[1]*v**2,F(16,3)*P[0]*u*v])
    Eya=12*u*dpsiy+48*A*u*dpsi+2*vy*dchi
    Eab=24*u*hess+12*dpsi**2+2*dchi**2+4*u*wave+4*contraction+radial+F(2,3)*P[1]*v**2+F(8,3)*P[0]*u*v
    scalar=12*u*(vtt+3*T*vt+8*v)/R**2+8*P[1]*u*v+12*uy*vy+P[2]*v**2
    values={'u':u,'u_t':ut,'u_tt':utt,'u_y':uy,'u_yy':uyy,'u_yt':uyt,
            'v':v,'v_t':vt,'v_tt':vtt,'v_y':vy,'Einstein_nn':Eyy,
            'Einstein_na_each':Eya,'Einstein_ab_each':Eab,'scalar':scalar}
    for key,value in values.items():
        need(value==F(source['bounds'][key]['exact_upper']),'independent bulk arithmetic '+key)
        need(value<F(source['bounds'][key]['strict_simple_upper']),'independent bulk cap '+key)
    return {'independent_reconstruction_equalities':10,'independent_full_source_envelope_equalities':4,
            'strict_Einstein_component_bound':'1200000000000000000000',
            'strict_scalar_bound':'10000000000000',
            'source_gauge':'first longitudinal representative; adapted background orthonormal frame',
            'second_order_solution_certified':False,'relative_harmonic_source_certified':False}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--root',type=Path,required=True)
    ap.add_argument('--inputs-root',type=Path,help='optional retained input root; defaults to capsule root')
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    actual_jets=independent_jet_bounds(args.root,args.inputs_root or args.root)
    report={'status':'PASS','optimized':not __debug__,'potential':independent_potential(args.root),'reflection':reflection_checks(),'incoming_and_bump':incoming_checks(),
            'actual_jets':actual_jets,'actual_bulk':independent_bulk_bounds(args.root,actual_jets['exact_scaled_order_bounds']),
            'scope':'Independent interval Taylor/ODE and exact polynomial witnesses; analytical smoothness and retained hypotheses reviewed separately.'}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','potential':report['potential']['ode_taylor_point_jet_witnesses'],'reflection':report['reflection']}))


if __name__=='__main__':main()
