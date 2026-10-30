"""Exact sparse Laurent differential algebra; standard library only.

Every expression is a collected polynomial over Q in the displayed generators,
allowing negative exponents only where division by a monomial is requested.
The derivation is the authenticated actual background/packet z derivation.
"""
from fractions import Fraction as F

NAMES = ('rho','p','a','eta','alpha','E')+tuple('Q'+str(j) for j in range(9))
INDEX = {name:i for i,name in enumerate(NAMES)}
ZERO_POWER = (0,)*len(NAMES)

class Expr:
    def __init__(self, terms=None):
        self.terms={tuple(k):F(v) for k,v in (terms or {}).items() if v}
    @classmethod
    def cv(cls,value):
        return value if isinstance(value,cls) else cls({ZERO_POWER:F(value)})
    @classmethod
    def symbol(cls,name):
        key=list(ZERO_POWER);key[INDEX[name]]=1
        return cls({tuple(key):F(1)})
    def __add__(self,other):
        out=dict(self.terms)
        for k,v in self.cv(other).terms.items():out[k]=out.get(k,F(0))+v
        return Expr(out)
    __radd__=__add__
    def __neg__(self):return Expr({k:-v for k,v in self.terms.items()})
    def __sub__(self,other):return self+-self.cv(other)
    def __rsub__(self,other):return self.cv(other)+-self
    def __mul__(self,other):
        out={}
        for k,v in self.terms.items():
            for l,w in self.cv(other).terms.items():
                key=tuple(a+b for a,b in zip(k,l));out[key]=out.get(key,F(0))+v*w
        return Expr(out)
    __rmul__=__mul__
    def inv(self):
        if len(self.terms)!=1:raise ValueError('division requires a nonzero monomial')
        (k,v),=self.terms.items()
        return Expr({tuple(-a for a in k):1/v})
    def __truediv__(self,other):return self*self.cv(other).inv()
    def __rtruediv__(self,other):return self.cv(other)*self.inv()
    def __pow__(self,n):
        if not isinstance(n,int):raise ValueError('integer exponent required')
        if n<0:return self.inv()**(-n)
        out=Expr.cv(1);base=self
        while n:
            if n%2:out=out*base
            n//=2
            if n:base=base*base
        return out
    def __eq__(self,other):return self.terms==self.cv(other).terms
    def __bool__(self):return bool(self.terms)
    def partial(self,name):
        i=INDEX[name];out={}
        for k,v in self.terms.items():
            if k[i]:
                key=list(k);key[i]-=1;out[tuple(key)]=v*k[i]
        return Expr(out)
    def dz(self):
        rhs=derivative_table();out=Expr.cv(0)
        for name in NAMES:
            if name=='Q8':
                if self.partial(name):raise ValueError('packet derivative exceeds formal Q8')
            else:out+=self.partial(name)*rhs[name]
        return out
    def encode(self):
        return [[list(key),str(value)] for key,value in sorted(self.terms.items())]
    def __repr__(self):return 'Expr('+str(len(self.terms))+' terms)'

def symbols():return {name:Expr.symbol(name) for name in NAMES}
def poly(coefficients,x):
    out=Expr.cv(0)
    for c in reversed(coefficients):out=out*x+c
    return out
POTENTIAL = tuple(map(F,('-2/27','0','14/9','50/27','-1/6','-4/9','-2/27')))
UP=[POTENTIAL]
for _ in range(5):UP.append(tuple(k*c for k,c in enumerate(UP[-1]) if k))
_DERIVATIVE_TABLE=None
def derivative_table():
    global _DERIVATIVE_TABLE
    if _DERIVATIVE_TABLE is None:
        s=symbols();rho,p,a,eta,alpha,E=[s[k] for k in NAMES[:6]]
        vp=poly(UP[1],eta)
        _DERIVATIVE_TABLE={'rho':-rho*rho*a,'eta':-rho*p,
          'p':-rho*(vp-4*a*p),'a':1/rho+rho*p*p/3,
          'alpha':rho*(F(5,2)*a-vp/p)*alpha,'E':-E}
        _DERIVATIVE_TABLE.update({'Q'+str(j):s['Q'+str(j+1)] for j in range(8)})
    return _DERIVATIVE_TABLE
def dy(value):return -Expr.cv(value).dz()/Expr.symbol('rho')

def physical_table():
    """Independent equations (3),(6),(7),(9),(10),(12), physical-Cauchy affine."""
    s=symbols();rho,p,a,eta,alpha,E=[s[k] for k in NAMES[:6]]
    vp,vpp=poly(UP[1],eta),poly(UP[2],eta)
    u=alpha/rho**2*s['Q0']
    ut=E*u
    # First coupled reconstruction, rather than the registered v formula.
    uy=dy(u);v=-3*(uy+2*a*u)/p
    uyt=dy(ut);vt=-3*(uyt+2*a*ut)/p
    vy=dy(v);uyy=dy(uy)
    A=(-24*uy**2-288*a*u*uy-24*u*uyy
       -(192*dy(a)+384*a*a+32*p*p)*u*u
       -(336*u*u+12*ut*ut+2*vt*vt)/rho**2
       -2*vy*vy-16*p*u*vy-2*vpp*v*v)
    D=(36*u*u-2*v*v)/rho**2
    Fmom=(18*u*ut-v*vt)/rho
    Gmom=(12*u*uyt-48*a*u*ut-6*uy*ut-2*vt*vy)/rho
    sectors=(A/6+F(4,3)*D,A/2+2*D,A-4*D)
    fs=tuple(-rho*rho*source/den for source,den in zip(sectors,(-6,10,42)))
    fy=tuple(dy(f) for f in fs)
    return dict(s,u=u,v=v,ut=ut,vt=vt,uy=uy,vy=vy,uyt=uyt,uyy=uyy,
                A=A,D=D,Fmom=Fmom,Gmom=Gmom,sectors=sectors,fs=fs,fy=fy,
                K=Fmom/3,Pi=(Gmom-dy(Fmom))/p,vp=vp,vpp=vpp)
