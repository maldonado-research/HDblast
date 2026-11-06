#!/usr/bin/env python3
"""Actual affine initial source jets: rational, global, no bump sampling.

All series coefficients bound |d_t^n h|/n!, t=(y-y0)/ell, uniformly
over y0 in the closed seed support. Coefficient rounding is outward and
performed using integers/Fractions only. No producer is imported; all finite majorants are recomputed here.
"""
from fractions import Fraction as F
from math import factorial

N = 16
SIG = 14
ELL = F(1, 10**6)
KAPPA = ELL/100
VCOEFF = [F(-2,27),0,F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)]


def need(c, message):
    if not c:
        raise ValueError(message)


def pow10(n):
    return F(10**n) if n >= 0 else F(1, 10**(-n))


def up(x):
    """Ceiling to SIG significant decimal digits, exactly."""
    x = F(x)
    need(x >= 0, 'nonnegative majorant coefficient')
    if not x:
        return F(0)
    e = len(str(x.numerator))-len(str(x.denominator))
    while x < pow10(e):
        e -= 1
    while x >= pow10(e+1):
        e += 1
    unit = pow10(e-SIG+1)
    z = x/unit
    value = ((z.numerator+z.denominator-1)//z.denominator)*unit
    need(value >= x, 'outward rational rounding')
    return value


class J:
    """Nonnegative normalized Taylor majorant."""
    def __init__(self, seq=()):
        self.c = [up(v) for v in list(seq)[:N+1]]
        self.c += [F(0)]*(N+1-len(self.c))
    @staticmethod
    def cast(x):
        return x if isinstance(x, J) else J([abs(F(x))])
    def __add__(self, other):
        other = J.cast(other)
        return J([a+b for a,b in zip(self.c,other.c)])
    __radd__ = __add__
    def __mul__(self, other):
        other = J.cast(other)
        return J([sum((self.c[k]*other.c[n-k] for k in range(n+1)),F(0))
                  for n in range(N+1)])
    __rmul__ = __mul__
    def __truediv__(self, value):
        value=F(value)
        need(value!=0,'nonzero scalar majorant divisor')
        return self*(1/abs(value))
    def __pow__(self, n):
        need(type(n) is int and n >= 0, 'nonnegative integral majorant power')
        out = J([1])
        for _ in range(n):
            out = out*self
        return out
    def derivative(self):
        # Only coefficients up to N-1 are claimed after differentiation.
        return J([(n+1)*self.c[n+1] for n in range(N)])
    def dy(self):
        return (1/ELL)*self.derivative()
    def inverse(self, value_lower):
        """Reciprocal of a signed jet with |h(0)| >= value_lower > 0."""
        lower = F(value_lower)
        need(lower > 0, 'reciprocal constant away from zero')
        out = [up(1/lower)]
        for n in range(1,N+1):
            out.append(up(sum((self.c[k]*out[n-k] for k in range(1,n+1)),F(0))/lower))
        return J(out)
    def exponential(self, value_upper):
        out = [up(value_upper)]
        for n in range(1,N+1):
            out.append(up(sum((k*self.c[k]*out[n-k] for k in range(1,n+1)),F(0))/n))
        return J(out)


def polynomial(coeff, x):
    out = J()
    for c in reversed(coeff):
        out = out*x+abs(c)
    return out


def v_derivative(n):
    return [VCOEFF[j]*F(factorial(j),factorial(j-n)) for j in range(n,len(VCOEFF))]


def background():
    """Simultaneous coefficient induction in a polynomial ODE system."""
    initial = {'rho':131,'a':F(112,1000),'p':F(14,100000),
               'eta':F(9,100000),'b':F(1,129),'c':F(100000,12),
               'B':F(6,5*129**2),'H':1500,'E':1}
    need(131**3 < 1500**2, 'rho to the three halves initial ceiling')
    cols = {key:[up(value)] for key,value in initial.items()}
    for n in range(N):
        j = {key:J(value) for key,value in cols.items()}
        rho,a,p,eta,b,c,B,H,E = (j[k] for k in initial)
        V1 = polynomial(v_derivative(1),eta)
        rhs = {'rho':a*rho,'a':b*b+p*p/F(3),
               'p':V1+4*a*p,'eta':p,'b':a*b,
               'c':V1*c*c+4*a*c,'B':(V1*c+F(9,2)*a)*B,
               'H':F(3,2)*a*H,'E':b*E}
        for key in cols:
            cols[key].append(up(ELL*rhs[key].c[n]/(n+1)))
    return {key:J(value) for key,value in cols.items()}


def bump_bounds():
    """Exact Bell-polynomial bound for sup |Q^(n)(x)|/n!, n<=16.

    With B=1/[x(1-x)]>=4, |B^(k)|/k! <= B^(k+1).
    P_n(B) is [t^n]exp(sum_{k>=1} B^(k+1)t^k).
    """
    polys = [{0:F(1)}]
    for n in range(1,N+1):
        out = {}
        for k in range(1,n+1):
            for degree,value in polys[n-k].items():
                d = degree+k+1
                out[d] = out.get(d,F(0))+F(k,n)*value
        polys.append(out)
    def maximum(m):
        # e>2 by 1+1+1/2>2. For m>=4, maximum at B=m.
        return F(4**m) if m <= 4 else F(m**m,2**(m-4))
    bounds = [sum((value*maximum(m) for m,value in p.items()),F(0)) for p in polys]
    need(bounds[0] == 1 and bounds[1] == 16, 'global bump constants zero and one')
    need(all(all(v>=0 for v in p.values()) for p in polys), 'positive Bell-polynomial recursion')
    return bounds,polys


def compose_bump(bounds, b):
    # x_t=-1/rho=-b; x(y0) in [0,1]. h is its zero-constant increment.
    h = J([0]+[b.c[n-1]/n for n in range(1,N+1)])
    out = J()
    power = J([1])
    for n in range(N+1):
        out = out+bounds[n]*power
        power = power*h
    return out,h


def source_jets(j,Q):
    rho,a,p,eta,b,c,B,H,E = (j[k] for k in j)
    u = B*Q
    uy,uyy = u.dy(),u.dy().dy()
    w = E*u
    v = 3*(uy+2*a*u)*c
    vt = 3*(w.dy()+2*a*w)*c
    v_y = v.dy()
    t = F(3,2)*w*b
    chi = vt*b
    V2 = polynomial(v_derivative(2),eta)
    Fjet = (12*u*uyy+12*uy*uy+144*a*u*uy
            +(192*a*a+96*a.dy()+168*b*b)*u*u
            +(v_y+4*u*p)**2+V2*v*v+(6*w*w+vt*vt)*b*b)
    Jjet = (18*u*u+v*v)*b*b
    R = 8*uy*t+2*chi*(v_y+4*u*p)
    L = 8*u*t+2*chi*v
    f0 = Fjet*F(1,6)+Jjet*F(4,3)
    f24 = Fjet+4*Jjet
    return {'u':u,'v':v,'w_velocity':w,'vt':vt,
            'F':Fjet,'J':Jjet,'R':R,'L':L,
            'f0':f0,'f2':f24,'f4':f24,'r0':R*F(1,6),
            'r2':R,'r4':R,'g2':L*F(1,2),'g4':L*F(1,2)}


def ds_values(jet, q_over_ell, max_order):
    out=[]
    for order in range(max_order+1):
        out.append(jet.c[0])
        jet = q_over_ell*jet.derivative()
    return out


def decade_ceiling(value):
    value = F(value)
    if value == 0:
        return 0
    e = len(str(value.numerator))-len(str(value.denominator))
    while value > pow10(e):
        e += 1
    while value <= pow10(e-1):
        e -= 1
    return e


def sqrt_ceiling(value):
    """Coarse integer power of ten >= sqrt(value), exactly checked."""
    e = decade_ceiling(value)
    out = pow10((e+1)//2)
    need(out*out >= value, 'rational square-root upper bound')
    return out


def stringify(x):
    if isinstance(x,F):
        return str(x)
    if isinstance(x,dict):
        return {k:stringify(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [stringify(v) for v in x]
    return x


def derive():
    need(N>=16,'sixteen stored Taylor orders for source and coefficient derivative margin')
    need(ELL==F(1,10**6) and KAPPA==F(1,10**8),'same explicitly distinguished weight member')
    need(F(129,131)**3*F(12,14)**2 > F(4,5)**2
         and F(131,129)**3*F(14,12)**2 < F(6,5)**2, 'alpha reconstruction bound')
    j=background()
    bump,bump_polys=bump_bounds()
    Q,h=compose_bump(bump,j['b'])
    src=source_jets(j,Q)
    # Known canonical-to-physical second-data conversion profiles. These use
    # only first seed fields and are distinct from the background B=alpha/rho².
    adm_profiles={'u_w_velocity_over_rho':src['u']*src['w_velocity']*j['b'],
                  'u_vt_over_rho':src['u']*src['vt']*j['b']}
    adm_dy={key:[factorial(n)*value.c[n]/ELL**n for n in range(14)]
            for key,value in adm_profiles.items()}
    adm_ceiling={'u_w_velocity_over_rho':F(10**73),'u_vt_over_rho':F(10**84)}
    for key,values in adm_dy.items():
        need(all(value<=adm_ceiling[key] for value in values),
             'known ADM conversion profile Dy derivatives through13 '+key)
    # r=1-(x-1/2)^2, x in [0,1]: r0<=1, |1-2x0|<=1.
    r=J([1])+h+h*h
    q=KAPPA*KAPPA*r*r
    q_over_ell=q*(1/ELL)
    rinv=r.inverse(F(3,4))
    # Prove exp(4/3)<4 by Taylor plus geometric upper tail.
    m=16; z=F(4,3)
    exp_upper=sum((z**k/factorial(k) for k in range(m+1)),F(0))
    exp_upper+=(z**(m+1)/factorial(m+1))/(1-z/F(m+2))
    need(exp_upper<4,'rational exp four-thirds upper bound')
    erinv=rinv.exponential(4)
    factor=KAPPA*r*j['H']*erinv
    transformed={key:factor*src[key] for key in ('f0','f2','f4','r0','r2','r4')}
    transformed['momentum_gradient_24']=factor*5*src['g2']*j['b']
    sup={key:ds_values(value,q_over_ell,12) for key,value in transformed.items()}
    # |K|_dy <= rho_max ell; q >= kappa² (3/4)², so |K|_ds <= this.
    support_ds=131*ELL/(KAPPA*KAPPA*F(3,4)**2)
    sqrt_support_ds=sqrt_ceiling(support_ds)
    l2={key:[sqrt_support_ds*v for v in values] for key,values in sup.items()}
    scalar={key:sum(l2[key],F(0)) for key in ('f0','f2','f4')}
    momentum={'0':sum(l2['r0'],F(0)),
              '2':sum((l2['r2'][k]+l2['momentum_gradient_24'][k] for k in range(13)),F(0)),
              '4':sum((l2['r4'][k]+l2['momentum_gradient_24'][k] for k in range(13)),F(0))}
    # C matrix entry l1 majorant; sqrt2<=2, sqrt12<=4, sqrt(lambda)<=5,
    # sqrt(d_lambda)<=4, uniformly lambda=8,24. Entry l1 controls op norm.
    # For coefficient estimates on the WHOLE annulus x in [-1/2,3/2],
    # |1-2x0|<=2. The source factor above only uses the seed x in [0,1].
    rann=J([1])+2*h+h*h
    qann=KAPPA*KAPPA*rann*rann
    qann_over_ell=qann*(1/ELL)
    b0=KAPPA*KAPPA*((1+rann)*rann.dy()+F(3,2)*j['a']*rann*rann)
    k=5*j['b']
    entries=[2*b0,2*qann*k,2*(b0+qann*j['a']),
             4*3*qann*j['a'],F(4,3)*qann*k,8*qann*j['b']]
    Cmaj=sum(entries,J())
    Cds=ds_values(Cmaj,qann_over_ell,14)
    need(all(x<F(1,10**10) for x in Cds),'whole-annulus momentum C Ds derivatives through14 below1e-10')
    geometric_coefficients={'q_over_rho':qann*j['b'], 'q_times_a':qann*j['a'],
                            'q_y':qann.dy(),'alpha_Ds_log_U':b0,
                            'alpha_physical_Ds_log_m':b0}
    geometric_ds={key:ds_values(value,qann_over_ell,14)
                  for key,value in geometric_coefficients.items()}
    need(all(v<=1 for vals in geometric_ds.values() for v in vals),
         'whole-annulus geometric coefficient Ds derivatives through14 below1')
    # Make consumer ceilings unambiguous; these assertions survive python -O.
    scalar_ceiling={key:pow10(decade_ceiling(v)) for key,v in scalar.items()}
    momentum_ceiling={key:pow10(decade_ceiling(v)) for key,v in momentum.items()}
    for key,value in scalar.items():
        need(value<=scalar_ceiling[key], 'scalar H12 norm ceiling '+key)
    for key,value in momentum.items():
        need(value<=momentum_ceiling[key], 'momentum H12 norm ceiling '+key)
    report={'schema':'hdblast-actual-source-jets-v1','status':'PASS',
            'arithmetic':'Python stdlib Fraction only; outward 14-significant-digit rational majorants',
            'orders':{'background':16,'global_bump':16,'source_Dy':12,'source_Ds':12,
                      'momentum_C_Ds':14,'known_ADM_conversion_profiles_Dy':13},
            'conventions':{'t':'(y-y0)/ell','Ds':'q Dy','q':'kappa^2 r^2',
                           'scalar_G':'kappa r rho^(3/2) exp(1/r) f_lambda',
                           'momentum_G':'kappa r rho^(3/2) exp(1/r) (r_lambda,sqrt(lambda) g_lambda/rho)',
                           'source_basis':'public H0,H2,H4, unnormalized; unit-harmonic conversion multiplies by ||H_l||_prob <=1',
                           'H12':'sum_{j=0}^{12} ||Ds^j G||_{L2(ds)}'},
            'scalar_H12_upper':scalar_ceiling,'momentum_H12_upper':momentum_ceiling,
            'source_Dy_derivative_upper':{key:[factorial(n)*src[key].c[n]/ELL**n for n in range(13)]
                                          for key in ('F','J','R','L','f0','f2','f4','r0','r2','r4','g2','g4')},
            'quadratic_ADM_conversion_profile_Dy_upper':adm_dy,
            'quadratic_ADM_conversion_profile_Dy13_uniform_upper':adm_ceiling,
            'transformed_source_Ds_sup_upper':sup,'transformed_source_Ds_L2_upper':l2,
            'scalar_H12_rational_majorant':scalar,'momentum_H12_rational_majorant':momentum,
            'momentum_C_Ds_opnorm_upper':Cds,'support_ds_upper':support_ds,
            'whole_annulus_geometric_coefficient_Ds_upper':geometric_ds,
            'sqrt_support_ds_upper':sqrt_support_ds,
            'background_normalized_taylor_upper':{key:value.c for key,value in j.items()},
            'bump_normalized_x_derivative_upper':bump,
            'bump_majorant_polynomials':[{str(k):v for k,v in p.items()} for p in bump_polys],
            'r_normalized_taylor_upper':r.c,'q_over_ell_normalized_taylor_upper':q_over_ell.c,
            'scope':'Conditional on retained smooth actual branch and all-point box. '
                    'Finite global derivative majorants for actual canonical affine source only. '
                    'No inverse regularity theorem, evolved field, amplitude radius, nonlinear remainder or physical calibration certified.'}
    return stringify(report)

