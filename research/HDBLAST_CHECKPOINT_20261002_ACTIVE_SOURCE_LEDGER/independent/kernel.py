"""Independent direct-stress and dense-source kernels; importing evaluates no source."""
from __future__ import annotations
import math
import numpy as np
from numpy.polynomial.legendre import leggauss

LD = np.longdouble
CD = np.clongdouble


class Pair:
    """First directional derivative, with no finite differencing."""
    __slots__ = ('v', 'd')
    __array_priority__ = 10000
    def __init__(self, value, direction=0):
        self.v, self.d = value, direction
    @staticmethod
    def of(value):
        return value if isinstance(value, Pair) else Pair(value)
    def __add__(self, other):
        o = Pair.of(other)
        return Pair(self.v + o.v, self.d + o.d)
    __radd__ = __add__
    def __neg__(self):
        return Pair(-self.v, -self.d)
    def __sub__(self, other):
        return self + (-Pair.of(other))
    def __rsub__(self, other):
        return Pair.of(other) + (-self)
    def __mul__(self, other):
        o = Pair.of(other)
        return Pair(self.v*o.v, self.d*o.v+self.v*o.d)
    __rmul__ = __mul__
    def __truediv__(self, other):
        o = Pair.of(other)
        value = self.v/o.v
        return Pair(value, (self.d-value*o.d)/o.v)
    def __rtruediv__(self, other):
        return Pair.of(other)/self
    def __pow__(self, power):
        return Pair(self.v**power, power*self.v**(power-1)*self.d)


def polynomials(order=5):
    p = [1]
    out = [p]
    for n in range(order):
        # Ascending coefficient recurrence (1-x^2)^2 P_n' +
        # 4n*x*(1-x^2)*P_n -2x*P_n.
        nxt = [0]*(len(p)+4)
        for j,c in enumerate(p):
            if j:
                nxt[j-1] += j*c
                nxt[j+1] -= 2*j*c
                nxt[j+3] += j*c
            nxt[j+1] += (4*n-2)*c
            nxt[j+3] -= 4*n*c
        while len(nxt)>1 and nxt[-1] == 0:
            nxt.pop()
        p = nxt
        out.append(p)
    return tuple(tuple(p) for p in out)

BUMP_POLYNOMIALS = polynomials()


def source_jet(eta, source):
    if source not in ('positive_B','signed_uB'):
        raise ValueError('Unregistered source')
    original_shape=np.shape(eta)
    z = np.atleast_1d(np.asarray(eta,dtype=LD))+LD(4)
    inside = np.abs(z)<1
    selected = z[inside]
    d = (1-selected)*(1+selected)
    bump = np.exp(-selected*selected/d)
    out = np.zeros((6,)+z.shape,dtype=LD)
    for n,coefficients in enumerate(BUMP_POLYNOMIALS):
        p = np.zeros_like(selected)
        for c in reversed(coefficients):
            p = p*selected+c
        out[n][inside] = bump*p/d**(2*n)
    if source == 'signed_uB':
        out = np.stack([z*out[n]+(n*out[n-1] if n else 0) for n in range(6)])
    return out.reshape((6,)+original_shape)


def forcing(eta, jet):
    L = -LD(1)/np.asarray(eta,dtype=LD)
    return 4*L*L*jet[0]-2*L*jet[1]-jet[2]


def propagator(k, distance):
    phase = CD(2j)*k*distance
    increment = np.expm1(phase)
    return 1+increment, increment/(CD(2j)*k)


def frequency_pairs(k, L, jet):
    """Differentiate omega^2=k^2+2L^2 and its four time derivatives."""
    M = 2*L*L
    mn = [math.factorial(n+1)*M*L**n for n in range(5)]
    dm = [2*sum(math.comb(n,j)*mn[j]*jet[n-j] for j in range(n+1)) for n in range(5)]
    value = [np.sqrt(k*k+M)]
    direction = [dm[0]/(2*value[0])]
    for n in range(1,5):
        vn = (mn[n]-sum(math.comb(n,j)*value[j]*value[n-j] for j in range(1,n)))/(2*value[0])
        value.append(vn)
        dn = dm[n]-sum(math.comb(n,j)*(direction[j]*value[n-j]+value[j]*direction[n-j])
                       for j in range(1,n))-2*direction[0]*vn
        direction.append(dn/(2*value[0]))
    return [Pair(v,d) for v,d in zip(value,direction)]


def full_subtractions(k, L, jet):
    """Published W0/W2/W4 inventory, applied with independent pair algebra."""
    w,w1,w2,w3,w4 = frequency_pairs(k,L,jet)
    C = Pair(2*L*L,jet[2]+2*L*jet[1])
    C1 = Pair(4*L**3,jet[3]+2*L*jet[2]+2*L*L*jet[1])
    C2 = Pair(12*L**4,jet[4]+2*L*jet[3]+4*L*L*jet[2]+4*L**3*jet[1])
    lp = Pair(L,jet[1])
    M = Pair(2*L*L,4*L*L*jet[0])
    U = -C/(2*w)-w2/(4*w**2)+3*w1**2/(8*w**3)
    U1 = -C1/(2*w)+C*w1/(2*w**2)-w3/(4*w**2)+5*w1*w2/(4*w**3)-9*w1**3/(8*w**4)
    U2 = (-C2/(2*w)+C1*w1/w**2+C*w2/(2*w**2)-C*w1**2/w**3-w4/(4*w**2)
          +7*w1*w3/(4*w**3)+5*w2**2/(4*w**3)-57*w1**2*w2/(8*w**4)+9*w1**4/(2*w**5))
    V = -U**2/(2*w)-U2/(4*w**2)+w2*U/(4*w**3)+3*w1*U1/(4*w**3)-3*w1**2*U/(4*w**4)
    b = -k*k/3-M
    c = 1-b/w**2
    J2 = lp*w1/w**2+w1**2/(4*w**3)
    J4 = lp*U1/w**2-2*lp*w1*U/w**3+w1*U1/(2*w**3)-3*w1**2*U/(4*w**4)
    R = w/2+(lp**2/w+J2)/4+(U**2/w-lp**2*U/w**2+J4)/4
    P = k*k/(6*w)+(c*U+lp**2/w+J2)/4+(c*V+b*U**2/w**3-lp**2*U/w**2+J4)/4
    return R.d,P.d


def baseline(k,L):
    """Stable published bare-minus-full-subtraction rational identities."""
    omega = np.sqrt(k*k+2*L*L)
    v = k/omega
    d = 2*L*L/(omega+k)
    def horner(coefficients):
        value = np.zeros_like(v)+LD(coefficients[-1])
        for c in reversed(coefficients[:-1]):
            value = value*v+c
        return value
    d2=d*d
    d4=d2*d2
    r0 = d4*horner((384,493,-204,-1188,-1940,-1990,-868,308,420,105))/(1024*k*omega**2)
    p0 = -d4*horner((384,1536,2944,3039,252,-5740,-15260,-19250,-8652,3612,4620,1155))/(3072*k*omega**2)
    return r0,p0


def direct_stresses(k,eta,u,w,jet,epsilon):
    L = -LD(1)/eta
    U,W = u/epsilon,w/epsilon
    common = -L*W.real-k*W.imag+L*jet[1]
    bare_r = ((2*k*k+3*L*L)*U.real+common+2*L*L*jet[0])/(2*k)
    bare_p = ((2*k*k/3-L*L)*U.real+common-2*L*L*jet[0])/(2*k)
    from specialized_contacts import directional_subtractions
    dr,dp = directional_subtractions(k,L,jet)
    r0,p0 = baseline(k,L)
    R,P = bare_r-dr-4*jet[0]*r0, bare_p-dp-4*jet[0]*p0
    work = -3*jet[1]*(r0+p0)
    F = L*(R-3*P)+work
    return R,P,F,work


class DenseKernel:
    """Fixed local source and dense time quadrature; no endpoint-ledger identity."""
    def __init__(self,k,dt,source_order,ledger_order):
        self.k,self.dt = k,dt
        x,v = leggauss(source_order)
        y,z = leggauss(ledger_order)
        self.s = (x.astype(LD)+1)/2
        self.sw = v.astype(LD)/2
        self.t = (y.astype(LD)+1)*dt/2
        self.tw = z.astype(LD)*dt/2
        self.E,self.D = propagator(k[:,None],self.t[None,:])
        distance = self.t[:,None]*(1-self.s[None,:])
        nested_E,nested_D = propagator(k[:,None,None],distance[None,:,:])
        nested_weights = self.t[:,None]*self.sw[None,:]
        self.kw = nested_E*nested_weights
        self.ku = nested_D*nested_weights
        self.Eend,self.Dend = propagator(k,dt)
        end_E,end_D = propagator(k[:,None],dt*(1-self.s[None,:]))
        self.end_kw = end_E*(dt*self.sw[None,:])
        self.end_ku = end_D*(dt*self.sw[None,:])
    def dense(self,u,w,g_nested,epsilon):
        wu = np.einsum('kqr,qr->kq',self.ku,g_nested,optimize=False)
        ww = np.einsum('kqr,qr->kq',self.kw,g_nested,optimize=False)
        return u[:,None]+self.D*w[:,None]-epsilon*wu,self.E*w[:,None]-epsilon*ww
    def step(self,u,w,g_end,epsilon):
        a = np.einsum('kr,r->k',self.end_ku,g_end,optimize=False)
        b = np.einsum('kr,r->k',self.end_kw,g_end,optimize=False)
        return u+self.Dend*w-epsilon*a,self.Eend*w-epsilon*b
