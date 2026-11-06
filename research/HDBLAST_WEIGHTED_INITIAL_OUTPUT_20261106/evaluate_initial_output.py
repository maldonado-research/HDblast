#!/usr/bin/env python3
"""Pure exact initial-output arithmetic; analytic and PDE premises are prose.

Only Python's standard library is used. No source producer is imported, no
floating-point operation is performed, and evaluate(root) writes no files.
"""
from fractions import Fraction as Q
import hashlib
import json
from math import factorial
from pathlib import Path


def need(condition, message):
    if not condition:
        raise ValueError(message)


PINS = {
    'EXACT_INITIAL_THEOREM.md': '95ccaa550309b72226f8d03d3a950e4e10760d0d98810f0e818d19057a8cfd0b',
    'REGISTERED_EXACT_FAMILY_BINDING.md': '69bfe846cc4fda8ec6a981969fbdd12fff7862b5ca7cf31a19d696b0a839165a',
    'CANONICAL_CONSTRAINTS.md': '4e25e993f9653b530e30e2a2d3fefcb0a006f04050534190b99a98126ba04f63',
    'EXPLICIT_COMPACT_SECOND_INITIAL_JET.md': '2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed',
    'ACTUAL_POTENTIAL_JET_BOUNDS.json': 'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
}

# Rehashed arithmetic corruption control: changing this must fail mathematics.
MOMENTUM_GAIN = 36


def load_json(raw):
    def unique(items):
        out = {}
        for key, value in items:
            need(key not in out, 'Duplicate JSON key: ' + key)
            out[key] = value
        return out
    def nonfinite(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    # Display-only decimal numbers in historical reports are never floats.
    return json.loads(raw.decode('utf-8'), object_pairs_hook=unique,
                      parse_float=str, parse_constant=nonfinite)


class I:
    """Closed rational interval, with sign-aware integer powers."""
    def __init__(self, lo, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        need(self.lo <= self.hi, 'Ordered rational interval required')
    @staticmethod
    def cast(x):
        return x if isinstance(x, I) else I(x)
    def __add__(self, x):
        x = I.cast(x)
        return I(self.lo+x.lo, self.hi+x.hi)
    __radd__ = __add__
    def __neg__(self):
        return I(-self.hi, -self.lo)
    def __sub__(self, x):
        return self+-I.cast(x)
    def __rsub__(self, x):
        return I.cast(x)+-self
    def __mul__(self, x):
        x = I.cast(x)
        vals = [a*b for a in (self.lo, self.hi) for b in (x.lo, x.hi)]
        return I(min(vals), max(vals))
    __rmul__ = __mul__
    def reciprocal(self):
        need(not self.lo <= 0 <= self.hi, 'Interval inverse away from zero')
        return I(1/self.hi, 1/self.lo)
    def __truediv__(self, x):
        return self*I.cast(x).reciprocal()
    def __pow__(self, n):
        need(type(n) is int, 'Integer interval power required')
        if n < 0:
            return self.reciprocal()**(-n)
        if n == 0:
            return I(1)
        lo = 0 if n % 2 == 0 and self.lo <= 0 <= self.hi else min(self.lo**n, self.hi**n)
        return I(lo, max(self.lo**n, self.hi**n))
    def absmax(self):
        return max(abs(self.lo), abs(self.hi))


VARS = ('rho', 'a', 'p', 'V1', 'V2', 'V3', 'f', 'f1', 'f2', 'f3',
        'n', 'n1', 'n2', 'n3', 'u', 'u1', 'u2', 'v', 'v1', 'w', 'vt',
        't', 'tp', 'chi', 'Y', 'G', 'F', 'J', 'H0', 'H2', 'H4', 'r', 'g')
ZERO = (0,)*len(VARS)


class P:
    """Small sparse Laurent polynomials with exact Fraction coefficients.

    Negative exponents are used only for the nonzero background rho,p.
    Equality is coefficient equality, never finite sampling.
    """
    def __init__(self, terms=None):
        self.terms = {k: Q(v) for k, v in (terms or {}).items() if v}
    @staticmethod
    def cast(x):
        return x if isinstance(x, P) else P({ZERO: Q(x)})
    @staticmethod
    def var(name):
        key = list(ZERO)
        key[VARS.index(name)] = 1
        return P({tuple(key): Q(1)})
    def __add__(self, x):
        out = self.terms.copy()
        for key, value in P.cast(x).terms.items():
            out[key] = out.get(key, Q(0))+value
        return P(out)
    __radd__ = __add__
    def __neg__(self):
        return P({k: -v for k, v in self.terms.items()})
    def __sub__(self, x):
        return self+-P.cast(x)
    def __rsub__(self, x):
        return P.cast(x)+-self
    def __mul__(self, x):
        out = {}
        for ka, va in self.terms.items():
            for kb, vb in P.cast(x).terms.items():
                key = tuple(a+b for a, b in zip(ka, kb))
                out[key] = out.get(key, Q(0))+va*vb
        return P(out)
    __rmul__ = __mul__
    def __pow__(self, n):
        need(type(n) is int, 'Integer Laurent power required')
        if n < 0:
            need(len(self.terms) == 1, 'Only Laurent monomials invertible')
            key, value = next(iter(self.terms.items()))
            return P({tuple(n*j for j in key): value**n})
        out = P.cast(1)
        for _ in range(n):
            out = out*self
        return out
    def __truediv__(self, x):
        return self*P.cast(x)**-1
    def __rtruediv__(self, x):
        return P.cast(x)*self**-1
    def partial(self, name):
        idx, out = VARS.index(name), {}
        for key, value in self.terms.items():
            if key[idx]:
                new = list(key)
                new[idx] -= 1
                out[tuple(new)] = value*key[idx]
        return P(out)


def algebra_checks(check):
    x = {key: P.var(key) for key in VARS}
    rho, a, p, V1, V2 = (x[k] for k in ('rho', 'a', 'p', 'V1', 'V2'))
    f, n, n1, n2 = (x[k] for k in ('f', 'n', 'n1', 'n2'))
    velocity = {'rho': a*rho, 'a': -rho**-2-p*p/3, 'p': V1-4*a*p,
                'V1': V2*p, 'V2': x['V3']*p,
                'f': x['f1'], 'f1': x['f2'], 'f2': x['f3'],
                'n': n1, 'n1': n2, 'n2': x['n3'],
                'u': x['u1'], 'u1': x['u2']}
    def D(expr):
        return sum((expr.partial(key)*value for key, value in velocity.items()), P())
    def zero(expr, name):
        check(not expr.terms, name)
    for lam in (0, 8, 24):
        Rict = 2/rho**2-D(a)-3*a*a
        A = (-3*a*n1+(3*a*a+(lam-3)/rho**2)*n)/2
        B = (-n2-2*a*n1+(Q(2*lam, 3)/rho**2-Rict)*n)/2
        S = (n1-a*n)/2
        divy = D(A)+3*a*(A-B)-lam*S/rho**2
        diva = D(S)+3*a*S+B-n*(lam-3)/(3*rho**2)
        zero(divy-p*p*(n1-a*n)/2, 'adjoint radial divergence lambda'+str(lam))
        zero(diva, 'adjoint angular divergence lambda'+str(lam))
        c = 2*rho*rho/(lam-3)
        b = 6*a*rho*rho/((lam-3)*p*p)
        zero(c*A+b*divy-n, 'differential scalar left inverse lambda'+str(lam))
        Ar, Br, Cr = c*f-D(b*f), -b*f/2, -a*b*f
        V = p*p/2-6*a*a+6/rho**2
        Pg = (3*a*D(Ar)/2+((lam+6)/(2*rho**2)-V)*Ar
              -3*D(D(Cr))/2-6*a*D(Cr)+(lam-3)*Cr/rho**2
              -lam*(D(Br)+2*a*Br)/rho**2)
        zero(Pg-f, 'compact scalar differential right inverse lambda'+str(lam))
        zero(Ar*Ar+3*Cr*Cr+2*lam*Br*Br/rho**2
             -((c*f-D(b*f))**2+3*a*a*b*b*f*f+lam*b*b*f*f/(2*rho**2)),
             'compact scalar tensor norm lambda'+str(lam))
        C = x['r']/(6*a)
        zero(6*a*C-x['r'], 'compact momentum radial equation lambda'+str(lam))
        if lam:
            dl = Q(2*(lam-3), 3)
            Dt = (C+x['g']/2)/dl
            zero(-2*C+2*dl*Dt-x['g'], 'compact momentum angular equation lambda'+str(lam))
            Dmetric = 3*rho*rho*f/(lam*(lam-3))
            zero(lam*(lam-3)*Dmetric/(3*rho*rho)-f,
                 'angular scalar right inverse lambda'+str(lam))
            zero(Dmetric*Dmetric*Q(2*lam*(lam-3), 3)
                 -6*rho**4*f*f/(lam*(lam-3)),
                 'angular scalar tensor norm lambda'+str(lam))
    u, u1, u2, v, v1, w, vt = (x[k] for k in ('u', 'u1', 'u2', 'v', 'v1', 'w', 'vt'))
    ay = D(a)
    FF = (-12*u*u2-12*u1*u1-144*a*u*u1
          -(192*a*a+96*ay+168/rho**2)*u*u
          -(v1+4*u*p)**2-V2*v*v-(6*w*w+vt*vt)/rho**2)
    JJ = (18*u*u-v*v)/rho**2
    oldA = (-24*u1*u1-288*a*u*u1-24*u*u2
            -(192*ay+384*a*a+32*p*p)*u*u
            -(336*u*u+12*w*w+2*vt*vt)/rho**2
            -2*v1*v1-16*p*u*v1-2*V2*v*v)
    oldD = (36*u*u-2*v*v)/rho**2
    zero(FF-oldA/2, 'canonical affine Hamiltonian equals half physical Hessian')
    zero(JJ-oldD/2, 'canonical affine angular source equals half physical Hessian')
    t, chi = -3*w/(2*rho), vt/rho
    zero(4*(t*t/3-t*t)-chi*chi+(6*w*w+vt*vt)/rho**2,
         'canonical density momentum Hamiltonian normalization')
    # Density divergence gives the first expression; use tp+3at=chi*p/2.
    Rfull = 16*u*(chi*p/2)+8*u1*t+2*chi*v1
    zero(Rfull-(8*u1*t+2*chi*(v1+4*u*p)), 'canonical momentum Hessian radial reduction')
    Fp, Jp, H0, H2, H4 = (x[k] for k in ('F', 'J', 'H0', 'H2', 'H4'))
    ysquare, grad = H0/6+H2+H4, 4*H0/3+4*H2-4*H4
    zero((-Fp/6-4*Jp/3)*H0+(-Fp-4*Jp)*H2+(-Fp+4*Jp)*H4+Fp*ysquare+Jp*grad,
         'negative actual Hamiltonian harmonic source decomposition')
    zero((Fp/6+4*Jp/3)**2+(Fp+4*Jp)**2/48+Q(19, 720)*(Fp-4*Jp)**2
         -(Q(3, 40)*Fp*Fp+Q(2, 5)*Fp*Jp+Q(38, 15)*Jp*Jp),
         'exact Hamiltonian angular source norm polynomial')
    check(Q(1, 36)+Q(1, 48)+Q(19, 720) == Q(3, 40), 'quadratic harmonic probability normalization')
    check(Q(8, 4*48)+Q(24*19, 4*720) == Q(1, 5), 'momentum gradient angular normalization')


def input_joins(root, check):
    pins = load_json((root/'SOURCE_PINS.json').read_bytes())
    check(isinstance(pins, dict) and isinstance(pins.get('files'), list), 'upstream pin inventory shape')
    locals_seen = set()
    for row in pins['files']:
        name = row.get('local')
        check(isinstance(name, str) and name.startswith('inputs/')
              and name[7:] in PINS and name not in locals_seen, 'upstream pin local member '+str(name))
        locals_seen.add(name)
        raw = (root/name).read_bytes()
        check(hashlib.sha256(raw).hexdigest() == PINS[name[7:]] == row.get('sha256')
              and len(raw) == row.get('bytes'), 'retained upstream bytes '+name[7:])
        check(row.get('ref') == '2355be3e7a6c3ecbdb9d7f41ab50aa05b1f3ea49', 'retained immutable upstream ref '+name[7:])
    check(locals_seen == {'inputs/'+name for name in PINS}, 'five retained upstream sources')
    actual = load_json((root/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_bytes())
    check(actual['status'] == 'PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ', 'retained conditional background status')
    check(actual['domain']['ell'] == '1/1000000', 'same actual packet length')
    check(actual['actual_model']['U_coefficients'] == ['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27'],
          'same exact model potential polynomial')
    expected = {'rho': (Q(129),Q(131)), 'a': (Q(111,1000),Q(14,125)),
                'p': (Q(-7,50000),Q(-3,25000)), 'eta': (Q(-9,100000),Q(-1,12500))}
    for key, pair in expected.items():
        check(tuple(map(Q,actual['background_intervals'][key]['exact'])) == pair,
              'actual wider collar background '+key)
    return expected


def source_bounds(check):
    eta = I(Q(-9,100000),Q(-8,100000))
    V1 = Q(28,9)*eta+Q(50,9)*eta**2-Q(2,3)*eta**3-Q(20,9)*eta**4-Q(4,9)*eta**5
    V2 = Q(28,9)+Q(100,9)*eta-2*eta**2-Q(80,9)*eta**3-Q(20,9)*eta**4
    V3 = Q(100,9)-4*eta-Q(80,3)*eta**2-Q(80,9)*eta**3
    check(V1.absmax()<Q(3,10000) and 3<V2.lo<=V2.hi<4 and V3.absmax()<12,
          'exact model potential first three derivative boxes')
    check(Q(129,131)**3*Q(12,14)**2>Q(4,5)**2
          and Q(131,129)**3*Q(14,12)**2<Q(6,5)**2, 'actual reconstruction alpha box')
    ell=Q(1,10**6); rlo=Q(129); rhi=Q(131)
    plo=Q(12,100000); phi=Q(14,100000); aa=Q(112,1000)
    ap=1/rlo**2+phi**2/3; hh=Q(3,10000)/plo
    K=hh+Q(9,2)*aa; Kp=4+hh**2+4*aa*hh+Q(9,2)*ap
    check(K<Q(301,100) and Kp<12, 'actual radial coefficient first derivative boxes')
    B=Q(6,5)/rlo**2; Bmin=Q(4,5)/rhi**2
    Qy=16/(ell*rlo); Qyy=736/(ell**2*rlo**2)+aa*Qy
    check(B<Q(1,10000) and B*(Qy+K)<10
          and B*(Qyy+2*K*Qy+K*K+Kp)<4000000, 'actual bump metric through second derivative')
    U=Q(1,10000); U1=Q(10); U2=Q(4000000); ratio=hh+4*aa
    vb=3*(U1+2*aa*U)/plo
    vpb=3*(U2+2*aa*U1+2*ap*U+(U1+2*aa*U)*ratio)/plo
    wy=U1+U/rlo
    vtb=3*(wy+2*aa*U)/plo
    check(vb<260000 and vpb<110000000000 and vtb<270000,
          'actual scalar and velocity profile bounds')
    VB=Q(260000); VP=Q(110000000000); VT=Q(270000)
    Fmax=(12*U*U2+12*U1**2+144*aa*U*U1+(192*aa**2+96*ap+168/rlo**2)*U**2
          +(VP+4*U*phi)**2+4*VB**2+(6*U**2+VT**2)/rlo**2)
    Jmax=(18*U**2+VB**2)/rlo**2
    Rmax=8*U1*(3*U/(2*rlo))+2*(VT/rlo)*(VP+4*U*phi)
    Lmax=8*U*(3*U/(2*rlo))+2*(VT/rlo)*VB
    check(Fmax<122*10**20 and Jmax<4100000, 'actual affine Hamiltonian source ceilings')
    check(Rmax<470000000000000 and Lmax<1100000000, 'actual affine momentum source ceilings')
    # b>=4 and exp(-1)<1/2 give the degree-five/six bump maxima.
    check(sum(Q(1,factorial(j)) for j in range(4))>2, 'elementary exponential one greater than two')
    check(Q(6**6,4)+Q(10*5**5,2)+24*4**4+12*4**3<35000,
          'global third bump derivative polynomial majorant')
    hp=4+hh**2+4*aa*hh
    app=2*aa/rlo**2+Q(2,3)*phi**2*ratio
    hpp=12*phi+(2*hh+4*aa)*hp+4*ap*hh
    Kpp=hpp+Q(9,2)*app
    check(Kpp<64, 'actual radial coefficient second derivative box')
    Qyyy=35000/(rlo**3*ell**3)+3*aa*736/(rlo**2*ell**2)+(ap+aa*aa)*Qy
    U3max=B*(Qyyy+3*K*Qyy+3*(K*K+Kp)*Qy+K**3+3*K*Kp+Kpp)
    check(U3max<1300000000000, 'actual metric third proper derivative box')
    U3=Q(1300000000000)
    G1=U1+2*aa*U; G1p=U2+2*aa*U1+2*ap*U
    G1pp=U3+2*aa*U2+4*ap*U1+2*app*U
    vpp=3*(G1pp+2*ratio*G1p+(ratio**2+hp+4*ap)*G1)/plo
    wyy=U2+2*U1/rlo+U*(1/rlo**2+aa/rlo)
    vtp=3*(wyy+2*aa*wy+2*ap*U+(wy+2*aa*U)*ratio)/plo
    check(vpp<33000000000000000 and vtp<110000000000,
          'actual scalar second and velocity first proper derivatives')
    cc=192*aa**2+96*ap+168/rlo**2
    cp=384*aa*ap+96*app+336*aa/rlo**2
    Fprime=(12*(U1*U2+U*U3)+24*U1*U2
            +144*(ap*U*U1+aa*U1**2+aa*U*U2)+cp*U**2+2*cc*U*U1
            +2*(VP+4*U*phi)*(33000000000000000+4*U1*phi+4*U*phi*ratio)
            +12*phi*VB**2+8*VB*VP+(12*U*wy+2*VT*110000000000)/rlo**2
            +2*aa*(6*U**2+VT**2)/rlo**2)
    Jprime=(36*U*U1+2*VB*VP)/rlo**2+2*aa*4100000
    check(Fprime<73*10**26 and Jprime<3500000000000, 'actual affine Hamiltonian source derivative ceilings')
    centeru=-Bmin*(32/(rhi**2*ell**2)-(K*K+Kp))
    centererr=2*ap*U+2*aa*U*K+(U*K+2*aa*U)*ratio
    centerv=3*(-86000+centererr)/phi
    check(centeru<-86000 and centerv+4*U*phi<-1800000000,
          'actual nonzero bump-center shifted scalar derivative')
    remainder=Fmax-(VP+4*U*phi)**2
    check(-1800000000**2+remainder<-3*10**18,
          'actual Hamiltonian center source strictly negative')
    check(Q(273,1000)**2<Q(3,40) and Q(1592,1000)**2>Q(38,15)
          and Q(273,1000)*3*10**18-Q(1592,1000)*4100000>8*10**17,
          'actual pointwise angular source norm nonzero witness')
    return {'F':Fmax,'J':Jmax,'R':Rmax,'L':Lmax,'F_y':Fprime,'J_y':Jprime,
            'u_yyy':U3max,'v_yy':vpp,'vt_y':vtp}


def budget_checks(check):
    amin=Q(111,1000);amax=Q(112,1000);rmin=Q(129);rmax=Q(131)
    pmin=Q(3,25000);pmax=Q(7,50000);ell=Q(1,10**6)
    matrices={}
    for lam in (8,24):
        dl=Q(2*(lam-3),3)
        m11=(3+Q(lam)/dl)/(36*amin**2)
        m22=rmax**2/(4*dl)
        m12sq=lam*rmax**2/(144*amin**2*dl**2)
        gap=(MOMENTUM_GAIN**2-m11)*(MOMENTUM_GAIN**2-m22)-m12sq
        check(m11<MOMENTUM_GAIN**2 and m22<MOMENTUM_GAIN**2 and gap>0,
              'momentum positive matrix gain36 lambda'+str(lam))
        matrices[str(lam)]={'M11':m11,'M22':m22,'M12_squared':m12sq,'det_gap':gap}
    check(1/(12*amin**2)<MOMENTUM_GAIN**2, 'momentum monopole gain36')
    exp_lower=sum(Q(1,3**j*factorial(j)) for j in range(9))
    exp_upper=exp_lower+Q(1,3**9*factorial(9))/(1-Q(1,30))
    check(exp_lower<exp_upper and Q(1,30)<1, 'rational exp one third Taylor tail enclosure')
    contrastH=Q(4,3)**4*exp_upper;contrastM=Q(4,3)**2*exp_upper
    check(contrastH<Q(442,100) and contrastM<Q(249,100), 'distinguished weight contrasts below4.42 and2.49')
    bmax=2*amax*rmax**2/pmin**2
    bpmax=2*rmax**2/pmin**2*(1/rmin**2+pmax**2/3+10*amax**2+2*amax*Q(28,100000)/pmin)
    cmax=2*rmax**2/3+bpmax+Q(7,4)*amax*bmax
    check(Q(7,4)**2>3 and bmax<267*10**9 and cmax<16*10**11,
          'explicit compact monopole comparison operator coefficients')
    check(Q(442,100)*bmax<1181*10**9 and Q(442,100)*cmax<708*10**10,
          'distinguished monopole output operator coefficients')
    check(6*rmax**4/(8*5)<6650**2 and 6*rmax**4/(24*21)<1873**2,
          'explicit angular scalar comparison gains')
    check(Q(442,100)*6650<29400 and Q(442,100)*1873<8280,
          'distinguished nonconstant scalar output operator gains')
    FF=Q(122)*10**20;JJ=Q(4100000);Fy=Q(73)*10**26;Jy=Q(3500000000000)
    B0=FF/6+4*JJ/3;B0y=Fy/6+4*Jy/3;B24=FF+4*JJ
    volume=rmax**4*ell;vhalf=Q(17161,1000)
    check(vhalf*vhalf==volume, 'exact packet support probability volume')
    H0=vhalf*(16*10**11*B0+267*10**9*B0y)
    bound2=H0**2+(vhalf*6650*B24)**2/48+Q(19,720)*(vhalf*1873*B24)**2
    check(bound2<(558*10**37)**2, 'compact metric comparison L2 below5.58e39')
    check(Q(442,100)**2*bound2<(25*10**39)**2, 'distinguished metric probability L2 below2.5e40')
    check(Q(99,70)**2>2 and Q(22,7)*Q(99,70)<Q(445,100), 'geometric sphere volume square root below4.45')
    check(Q(445,100)*(25*10**39)<112*10**39, 'distinguished metric geometric L2 below1.12e41')
    momentum2=Q(3,40)*470000000000000**2+Q(1100000000)**2/(5*rmin**2)
    check(momentum2<129000000000000**2, 'actual angular momentum source norm below1.29e14')
    check(vhalf*129000000000000<222*10**13, 'actual probability momentum source L2 below2.22e15')
    check(MOMENTUM_GAIN*Q(249,100)*(222*10**13)<2*10**17, 'distinguished momentum probability L2 below2e17')
    check(20*volume<77**2 and 77*129000000000000<10**16,
          'actual geometric momentum source L2 below1e16')
    check(MOMENTUM_GAIN*Q(249,100)*10**16<9*10**17, 'distinguished momentum geometric L2 below9e17')
    kappa=ell/100;qdelta=4*kappa*kappa/ell
    check(qdelta==ell/2500 and 4*qdelta<1 and 4*qdelta<rmin,
          'distinguished defining function fourfold ball inclusion')
    eps=8*kappa*kappa/(rmin*ell)
    check(eps<Q(1,10**11) and eps/(1-eps)<Q(1,10**11),
          'distinguished defining function and exponential weight comparability')
    return {'exp_one_third_lower':exp_lower,'exp_one_third_upper':exp_upper,
            'metric_contrast_upper':contrastH,'momentum_contrast_upper':contrastM,
            'momentum_PSD_matrices':matrices,'monopole_b_abs_upper':bmax,
            'monopole_value_coefficient_upper':cmax,'source_probability_volume_upper':volume,
            'comparison_metric_L2_squared_upper':bound2,'momentum_source_angular_squared_upper':momentum2,
            'metric_probability_L2_upper':25*10**39,'metric_geometric_L2_upper':112*10**39,
            'momentum_probability_L2_upper':2*10**17,'momentum_geometric_L2_upper':9*10**17,
            'cover_relative_r_variation_upper':eps}


def exact_strings(value):
    if isinstance(value,Q):
        return str(value)
    if isinstance(value,dict):
        return {key:exact_strings(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):
        return [exact_strings(item) for item in value]
    need(type(value) in (int,str,bool) or value is None, 'Exact report value required')
    return value


def evaluate(root):
    root=Path(root)
    names=[]
    def check(condition,name):
        need(condition,name)
        need(name not in names,'Duplicate mathematical check name: '+name)
        names.append(name)
    input_joins(root,check)
    algebra_checks(check)
    sources=source_bounds(check)
    budgets=budget_checks(check)
    return {'schema':'hdblast-weighted-initial-output-v1',
            'status':'PASS','named_checks':names,'named_check_count':len(names),
            'exact':exact_strings({'source_upper_bounds':sources,**budgets}),
            'scope':'Finite rational interval and Laurent-polynomial identities only. '
                    'Retained background and weighted-domain existence are premises. '
                    'Absolute L2 ceilings apply to the explicitly distinguished member, '
                    'not an unidentified historical weight. No H14, evolved second-field, '
                    'nonlinear radius or common target cylinder is certified.'}
