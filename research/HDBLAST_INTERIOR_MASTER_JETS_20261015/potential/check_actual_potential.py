#!/usr/bin/env python3
"""Exact retained-model potential and reconstruction coefficient audit.

Standard library only. No historical producer is imported or executed.
All inequalities use Fraction endpoints; float output is display only.
"""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json

ROOT = Path(__file__).resolve().parent
PINS = {
    'BACKGROUND_DERIVATIVE_BOUNDS.json': 'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
    'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md': '0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
    'OCT11_JACOBI_ENDPOINT_MONOTONICITY.md': 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8',
    'OCT11_RESCALED_BOUND_CERTIFICATE.json': '22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'BACKGROUND_DERIVATIVE_AUDIT.md': 'c85130916eb36a4a81c512cf594db8f9d4327dba926233aef126da95d45b4fff',
}


def require(ok, text):
    if not ok:
        raise ValueError(text)


class I:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        require(self.lo <= self.hi, 'reversed interval')

    @staticmethod
    def box(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, y):
        y = self.box(y)
        return I(self.lo+y.lo, self.hi+y.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, y):
        return self + -self.box(y)

    def __rsub__(self, y):
        return self.box(y) + -self

    def __mul__(self, y):
        y = self.box(y)
        corners = [a*b for a in (self.lo, self.hi) for b in (y.lo, y.hi)]
        return I(min(corners), max(corners))

    __rmul__ = __mul__

    def recip(self):
        require(not self.lo <= 0 <= self.hi, 'division across zero')
        return I(1/self.hi, 1/self.lo)

    def __truediv__(self, y):
        return self*self.box(y).recip()

    def __rtruediv__(self, y):
        return self.box(y)*self.recip()

    def __pow__(self, n):
        require(isinstance(n, int), 'integer power required')
        if n < 0:
            return self.recip()**(-n)
        if n == 0:
            return I(1)
        if n % 2 == 0:
            h = self**(n//2)
            return I(0 if h.lo <= 0 <= h.hi else min(h.lo*h.lo, h.hi*h.hi),
                     max(h.lo*h.lo, h.hi*h.hi))
        return self*(self**(n-1))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def record(self):
        return {'exact': [str(self.lo), str(self.hi)],
                'display_only': [float(self.lo), float(self.hi)]}


def deriv(poly):
    return [i*poly[i] for i in range(1, len(poly))]


def add(a, b):
    return [(a[i] if i < len(a) else F(0))+(b[i] if i < len(b) else F(0))
            for i in range(max(len(a), len(b)))]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def scale(a, x):
    return [F(x)*v for v in a]


def evaluate(poly, x):
    result = 0
    for c in reversed(poly):
        result = result*x+c
    return result


def same_polynomial(a, b):
    return all(v == 0 for v in add(a, scale(b, -1)))


class P:
    """Sparse Laurent polynomials, used as an independent formal derivative."""
    n = 7  # rho,a,p,h,U2,U3,U4

    def __init__(self, terms=None):
        self.d = {k:F(v) for k,v in (terms or {}).items() if v}

    @classmethod
    def box(cls, x):
        return x if isinstance(x, cls) else cls({(0,)*cls.n:F(x)})

    @classmethod
    def var(cls, i):
        e = [0]*cls.n
        e[i] = 1
        return cls({tuple(e):F(1)})

    def __add__(self, other):
        out = self.d.copy()
        for m,c in self.box(other).d.items():
            out[m] = out.get(m, F(0))+c
        return P(out)

    __radd__ = __add__

    def __neg__(self):
        return P({m:-c for m,c in self.d.items()})

    def __sub__(self, other):
        return self+-self.box(other)

    def __rsub__(self, other):
        return self.box(other)+-self

    def __mul__(self, other):
        out = {}
        for a,x in self.d.items():
            for b,y in self.box(other).d.items():
                m = tuple(i+j for i,j in zip(a,b))
                out[m] = out.get(m,F(0))+x*y
        return P(out)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self*F(1, scalar)

    def __pow__(self, n):
        if n < 0:
            require(len(self.d) == 1, 'negative polynomial power needs monomial')
            (m,c), = self.d.items()
            return P({tuple(i*n for i in m):c**n})
        result = P.box(1)
        for _ in range(n):
            result = result*self
        return result

    def diff(self, i):
        out = {}
        for m,c in self.d.items():
            if m[i]:
                mm = list(m)
                mm[i] -= 1
                out[tuple(mm)] = c*m[i]
        return P(out)


def formal_identities(corrupt):
    r,a,p,h,u2,u3,u4 = [P.var(i) for i in range(P.n)]
    rz, az, pz = -r*r*a, r**-1+r*p*p/3, -r*p*(h-4*a)
    hz = -r*(u2-h*h+4*a*h)
    velocity = (rz, az, pz, hz, -r*p*u3, -r*p*u4, P())

    def D(poly):
        return sum((poly.diff(i)*v for i,v in enumerate(velocity)), P())

    hzz = -rz*(u2-h*h+4*a*h)-r*(-r*p*u3+(-2*h+4*a)*hz+4*az*h)
    azz = -rz*r**-2+rz*p*p/3+2*r*p*pz/3
    pzz = -rz*p*(h-4*a)-r*(pz*(h-4*a)+p*(hz-4*az))
    Q = 2*h*h-10*a*h+F(35,4)*a*a-u2-p*p/6
    Qz = (4*h-10*a)*hz+(-10*h+F(35,2)*a)*az+r*p*u3-p*pz/3
    Qzz = 4*hz*hz+(4*h-10*a)*hzz-20*az*hz+F(35,2)*az*az+(-10*h+F(35,2)*a)*azz+(rz*p+r*pz)*u3-r*r*p*p*u4-(pz*pz+p*pzz)/3
    Wz = -2*r**3*a*Q+r*r*Qz
    cross = -2 if corrupt == 'drop_Wzz_cross_factor' else -4
    Wzz = (6*r**4*a*a-2*r**3*az)*Q+cross*r**3*a*Qz+r*r*Qzz
    checks = {'hzz_total_derivative':D(hz)-hzz,
              'azz_total_derivative':D(az)-azz,
              'pzz_total_derivative':D(pz)-pzz,
              'Qz_total_derivative':D(Q)-Qz,
              'Qzz_total_derivative':D(Qz)-Qzz,
              'Wz_total_derivative':D(r*r*Q)-Wz,
              'Wzz_total_derivative':D(Wz)-Wzz}
    if corrupt == 'drop_Qzz_reflection_cross':
        checks['Qzz_total_derivative'] = D(Qz)-(Qzz+20*az*hz)
    for name, defect in checks.items():
        require(not defect.d, 'formal identity failed: '+name)
    return {name:'PASS_ZERO_LAURENT_POLYNOMIAL' for name in checks}


def run(corrupt='none'):
    for name,pin in PINS.items():
        data = (ROOT/'inputs'/name).read_bytes()
        if corrupt == 'input_identity' and name == 'BACKGROUND_DERIVATIVE_BOUNDS.json':
            data += b' '
        require(sha256(data).hexdigest() == pin, 'input authentication failed: '+name)
    prior = json.loads((ROOT/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    bg = json.loads((ROOT/'inputs/BACKGROUND_DERIVATIVE_BOUNDS.json').read_text())
    cert = json.loads((ROOT/'inputs/OCT11_RESCALED_BOUND_CERTIFICATE.json').read_text())
    c, delta = F(96714321425692,161839258930731), F(1,1000)
    ell, T, Z = F(1,10**6), F(3,10**6), F(8,10**6)
    require(prior['parameters_exact'] == {'c':str(c),'delta':str(delta)}, 'retained shell parameters differ')
    require(bg['candidate_state_exact'] == {'delta':str(delta),'c':str(c),'ell':str(ell),'Z':str(Z),'T':str(T)}, 'wider audit parameters differ')
    for name in ('OCT11_JACOBI_ENDPOINT_MONOTONICITY.md','OCT11_RESCALED_BOUND_CERTIFICATE.json','PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'):
        require(bg['input_pins'][name] == PINS[name], 'wider audit input identity differs')
    W = [F(1,3),F(0),F(1),F(1,3)]
    U = add(scale(mul(deriv(W),deriv(W)),F(1,2)),scale(mul(W,W),F(-2,3)))
    sigma = add(scale(W,2), [delta*(1+c),delta*c])
    f = [1+c,c]
    Apoly = add(scale(mul(W,f),F(1,9)),scale(deriv(W),-c/12))
    Bpoly = add(scale(mul(f,f),F(1,36)),[-c*c/48])
    direct = add(add(scale(mul(sigma,sigma),F(1,36)),scale(mul(deriv(sigma),deriv(sigma)),F(-1,48))),scale(U,F(1,6)))
    require(same_polynomial(direct,add(scale(Apoly,delta),scale(Bpoly,delta*delta))), 'Hamiltonian shell cancellation failed')
    constants = cert['primary_uniform_certificate']['exact_constants']
    require(same_polynomial(Apoly,list(map(F,constants['A']))) and same_polynomial(Bpoly,list(map(F,constants['B_fun']))), 'registered polynomial certificate differs')
    outer = {k:I(*entry['outer']['exact']) for k,entry in bg['first_exit']['component_margins'].items()}
    r,e,p,a = [outer[k] for k in ('rho','eta','p','a')]
    Up = []
    poly = U
    for _ in range(5):
        Up.append(evaluate(poly,e))
        poly = deriv(poly)
    u0,u1,u2,u3,u4 = Up
    rz, ez, pz, az = -r*r*a, -r*p, -r*(u1-4*a*p), 1/r+r*p**2/3
    rhs = {'rho':rz,'eta':ez,'p':pz,'a':az}
    initial_names = {'rho':'rho_b','eta':'eta_b','p':'p_b','a':'a_b'}
    for k,derivative in rhs.items():
        require(derivative.record()['exact'] == bg['first_exit']['component_margins'][k]['rhs']['exact'], 'retained ODE image differs: '+k)
        init = I(*prior['boxes'][initial_names[k]]['exact'])
        tube = init+I(0,Z)*derivative
        require(outer[k].lo < tube.lo <= tube.hi < outer[k].hi, 'first-exit inclusion failed: '+k)
    require(131*Z < 9, 'proper collar endpoint margin failed')
    h = u1/p
    L = u2-h**2+4*a*h
    hz = -r*L
    hzz = -rz*L-r*(u3*ez+(-2*h+4*a)*hz+4*az*h)
    azz = -rz/r**2+rz*p**2/3+2*r*p*pz/3
    pzz = -rz*(u1-4*a*p)-r*(u2*ez-4*az*p-4*a*pz)
    rzz = 2*r**3*a**2-r**2*az
    Q = 2*h**2-10*a*h+F(35,4)*a**2-u2-p**2/6
    Qz = (4*h-10*a)*hz+(-10*h+F(35,2)*a)*az+r*p*u3-p*pz/3
    Qzz = 4*hz**2+(4*h-10*a)*hzz-20*az*hz+F(35,2)*az**2+(-10*h+F(35,2)*a)*azz+(rz*p+r*pz)*u3+r*p*u4*ez-(pz**2+p*pzz)/3
    Wz = -2*r**3*a*Q+r**2*Qz
    Wzz = (6*r**4*a**2-2*r**3*az)*Q-4*r**3*a*Qz+r**2*Qzz
    for name,value in {'h':h,'h_z':hz,'Q':Q,'Q_z':Qz,'W_z':Wz}.items():
        require(value.record()['exact'] == bg['outer_derived'][name]['exact'], 'retained derived interval differs: '+name)
    potential = r**2*Q-I(0,F(35,4)*T*T)
    bounds = {'absolute_W':F(102000),'absolute_W_z':F(30000000),
              'absolute_W_t':F(21,400000),'absolute_W_zz':F(30000000000)}
    if corrupt == 'undersized_Wzz':
        bounds['absolute_W_zz'] = F(1)
    for name,value in [('absolute_W',potential),('absolute_W_z',Wz),('absolute_W_zz',Wzz)]:
        require(value.mag() < bounds[name], 'claimed strict bound failed: '+name)
    require(F(35,2)*T == bounds['absolute_W_t'], 'temporal bound failed')
    symbolic = formal_identities(corrupt)
    # alpha(0)=1 and alpha_z/alpha=k; exp(KZ)<=1/(1-KZ).
    k = r*(F(5,2)*a-h)
    kz = rz*(F(5,2)*a-h)+r*(F(5,2)*az-hz)
    kzz = rzz*(F(5,2)*a-h)+2*rz*(F(5,2)*az-hz)+r*(F(5,2)*azz-hzz)
    require(k.mag()*Z < 1, 'alpha exponential majorant denominator failed')
    alpha_mag = 1/(1-k.mag()*Z)
    alpha = I(1/alpha_mag,alpha_mag)
    C = alpha/r**2
    l = r*(F(9,2)*a-h)
    lz = rz*(F(9,2)*a-h)+r*(F(9,2)*az-hz)
    Cz,Czz = l*C,(lz+l**2)*C
    B = 3*alpha/(r**3*p)
    m = F(3,2)*r*a
    mz = F(3,2)*(rz*a+r*az)
    Bz,Bzz = m*B,(mz+m**2)*B
    reconstruction = {'alpha':alpha,'alpha_log_z':k,'alpha_log_zz':kz,'alpha_log_zzz':kzz,
                      'u_w_coefficient':C,'u_w_coefficient_z':Cz,'u_w_coefficient_zz':Czz,
                      'v_wz_coefficient':B,'v_wz_coefficient_z':Bz,'v_wz_coefficient_zz':Bzz}
    report = {
        'status':'PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ',
        'scope':'Conditional on authenticated retained regular branch theorem; exact collar bounds, not a numerical trajectory or nonlinear solution.',
        'input_sha256':PINS,
        'domain':{'ell':str(ell),'T':str(T),'Z':str(Z),'triangle':'0<=t<=T, z>=0, z+t<=Z','verified_enlargement':'[0,T] x [0,Z]'},
        'actual_model':{'eta':'phi-1','c':str(c),'delta':str(delta),'W_super_coefficients':list(map(str,W)),'U_coefficients':list(map(str,U))},
        'retained_identity':'Authenticated exact h,h_z,Q,Q_z,W_z intervals and all four vector-field interval images equal retained BACKGROUND_DERIVATIVE_BOUNDS.json; same registered Hamiltonian polynomial coefficients.',
        'simple_bounds':{k:str(v) for k,v in bounds.items()},
        'potential_intervals':{'W':potential.record(),'W_z':Wz.record(),'W_zz':Wzz.record()},
        'background_intervals':{name:value.record() for name,value in {'rho':r,'eta':e,'p':p,'a':a,'rho_z':rz,'eta_z':ez,'p_z':pz,'a_z':az,'rho_zz':rzz,'p_zz':pzz,'a_zz':azz,'h':h,'h_z':hz,'h_zz':hzz,'Q':Q,'Q_z':Qz,'Q_zz':Qzz,'U':u0,'U_1':u1,'U_2':u2,'U_3':u3,'U_4':u4}.items()},
        'reconstruction_coefficients':{name:value.record() for name,value in reconstruction.items()},
        'formal_identities':symbolic,
        'analytic_temporal_bound':'|W_t|=(35/2)|tanh(t)|sech(t)^2 <= (35/2)T = 21/400000',
    }
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--corrupt',choices=['none','input_identity','drop_Wzz_cross_factor','drop_Qzz_reflection_cross','undersized_Wzz'],default='none')
    parser.add_argument('--output',type=Path,default=ROOT/'ACTUAL_POTENTIAL_JET_BOUNDS.json')
    args=parser.parse_args()
    report=run(args.corrupt)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'simple_bounds':report['simple_bounds'],'W_zz_display_only':report['potential_intervals']['W_zz']['display_only'],'output':str(args.output)},indent=2))


if __name__=='__main__':
    main()
