"""Small exact Laurent polynomial ring used only for verification.

Division is permitted by a single nonzero monomial, which is all that the
frame rescalings require. There is no simplification by sampled values.
"""
from fractions import Fraction as Q
NAMES=('ell','rho','a','h','c','Ht','p','py','b','d','V','V1','V2','sigma1','sigma2')
ZERO=(0,)*len(NAMES)
class P:
    def __init__(self,value=0):
        if isinstance(value,P):self.terms=value.terms.copy()
        elif isinstance(value,dict):self.terms={e:Q(v) for e,v in value.items() if v}
        else:self.terms={} if value==0 else {ZERO:Q(value)}
    def __add__(self,o):
        o=P(o);d=self.terms.copy()
        for e,v in o.terms.items():d[e]=d.get(e,Q(0))+v
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({e:-v for e,v in self.terms.items()})
    def __sub__(self,o):return self+-P(o)
    def __rsub__(self,o):return P(o)+-self
    def __mul__(self,o):
        o=P(o);d={}
        for e,v in self.terms.items():
            for f,w in o.terms.items():
                g=tuple(a+b for a,b in zip(e,f));d[g]=d.get(g,Q(0))+v*w
        return P(d)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=P(o)
        if len(o.terms)!=1:raise ValueError('EXACT_DIVISOR_MUST_BE_MONOMIAL')
        f,w=next(iter(o.terms.items()))
        return P({tuple(a-b for a,b in zip(e,f)):v/w for e,v in self.terms.items()})
    def __rtruediv__(self,o):return P(o)/self
    def __pow__(self,n):
        if not isinstance(n,int):raise ValueError('INTEGER_EXPONENT_REQUIRED')
        if n<0:return P(1)/(self**(-n))
        out=P(1)
        for _ in range(n):out=out*self
        return out
    def __eq__(self,o):return self.terms==P(o).terms
    def __bool__(self):return bool(self.terms)
    def __repr__(self):return repr(self.terms)
def variables():
    out=[]
    for i in range(len(NAMES)):
        e=[0]*len(NAMES);e[i]=1;out.append(P({tuple(e):Q(1)}))
    return out
