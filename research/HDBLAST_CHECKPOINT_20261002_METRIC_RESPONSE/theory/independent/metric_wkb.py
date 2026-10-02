"""Independent full metric directional WKB evaluator.

No field evolution, quadrature, source selection, or acceptance gate is defined.
The evaluator differentiates the inherited fixed-r, minimal x=r subtraction
at fixed conformal time and fixed comoving k. Taylor coefficients carry n!.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import comb, factorial
import numpy as np

ORDER = 4
HALF = 0.5
ONE = 1

@dataclass
class Jet:
    __array_priority__ = 1000
    coefficients: tuple

    @classmethod
    def constant(cls, value):
        zero = np.zeros_like(value)
        return cls((value, zero, zero, zero, zero))

    @classmethod
    def from_derivatives(cls, derivatives):
        if len(derivatives) != 5:
            raise ValueError('Exactly five derivatives, orders zero through four, are required')
        return cls(tuple(value/factorial(n) for n, value in enumerate(derivatives)))

    def __add__(self, other):
        if not isinstance(other, Jet): other = Jet.constant(other)
        return Jet(tuple(a+b for a,b in zip(self.coefficients, other.coefficients)))
    __radd__ = __add__
    def __neg__(self): return Jet(tuple(-x for x in self.coefficients))
    def __sub__(self, other): return self + (-other if isinstance(other,Jet) else -other)
    def __rsub__(self, other): return -self + other
    def __mul__(self, other):
        if not isinstance(other, Jet): other = Jet.constant(other)
        return Jet(tuple(sum(self.coefficients[j]*other.coefficients[n-j] for j in range(n+1)) for n in range(5)))
    __rmul__ = __mul__
    def inverse(self):
        a = self.coefficients
        if np.any(a[0] == 0): raise ValueError('Zero Taylor-jet denominator')
        out = [ONE/a[0]]
        for n in range(1,5):
            out.append(-sum(a[j]*out[n-j] for j in range(1,n+1))/a[0])
        return Jet(tuple(out))
    def __truediv__(self, other):
        return self*(other.inverse() if isinstance(other,Jet) else ONE/other)
    def __rtruediv__(self, other): return self.inverse()*other
    def __pow__(self, exponent):
        if isinstance(exponent,int):
            if exponent < 0: return self.inverse()**(-exponent)
            out = Jet.constant(np.ones_like(self.coefficients[0]))
            for _ in range(exponent): out = out*self
            return out
        a0=self.coefficients[0]
        if np.any(a0 <= 0): raise ValueError('Nonpositive reference frequency squared')
        q=self/a0-1
        out=Jet.constant(np.ones_like(a0)); power=out; coefficient=1
        for n in range(1,5):
            power=power*q; coefficient=coefficient*(exponent-n+1)/n
            out=out+coefficient*power
        return a0**exponent*out
    def derivative(self):
        c=self.coefficients
        return Jet(tuple((n+1)*c[n+1] for n in range(4))+(np.zeros_like(c[0]),))
    def actual_derivative(self, n=0): return factorial(n)*self.coefficients[n]

@dataclass
class DualJet:
    __array_priority__ = 1000
    value: Jet
    delta: Jet
    @classmethod
    def constant(cls,value): return cls(Jet.constant(value),Jet.constant(np.zeros_like(value)))
    def __add__(self, other):
        if not isinstance(other,DualJet): other=DualJet.constant(other)
        return DualJet(self.value+other.value,self.delta+other.delta)
    __radd__=__add__
    def __neg__(self): return DualJet(-self.value,-self.delta)
    def __sub__(self,other): return self+(-other if isinstance(other,DualJet) else -other)
    def __rsub__(self,other): return -self+other
    def __mul__(self,other):
        if not isinstance(other,DualJet): other=DualJet.constant(other)
        return DualJet(self.value*other.value,self.delta*other.value+self.value*other.delta)
    __rmul__=__mul__
    def __truediv__(self,other):
        if not isinstance(other,DualJet): other=DualJet.constant(other)
        inverse=other.value.inverse()
        return DualJet(self.value*inverse,self.delta*inverse-self.value*other.delta*inverse**2)
    def __rtruediv__(self,other): return DualJet.constant(other)/self
    def __pow__(self,exponent): return DualJet(self.value**exponent,exponent*self.value**(exponent-1)*self.delta)
    def derivative(self): return DualJet(self.value.derivative(),self.delta.derivative())


def generic_subtractions(k2, a_derivatives, delta_a_derivatives, r=2):
    """Return all canonical baseline and metric directional subtractions.

    Inputs are exact a^(n), delta a^(n), n=0..4; r and k2 are held fixed.
    R/P outputs require division by a^4; S outputs by a^2. Their physical
    directional variations also have -4h R/P and -2h S prefactor contacts.
    Only S, S0, S2, omega, and W2 have first/second derivative outputs: R/P higher
    derivatives would need additional background jets and are not exposed.
    """
    a=DualJet(Jet.from_derivatives(a_derivatives),Jet.from_derivatives(delta_a_derivatives))
    L=a.derivative()/a
    C=a.derivative().derivative()/a
    M=r*a*a
    w=(k2+M)**HALF
    w1=w.derivative();w2=w1.derivative()
    U=-C/(2*w)-w2/(4*w**2)+3*w1**2/(8*w**3)
    U1=U.derivative();U2=U1.derivative()
    V=-U**2/(2*w)-U2/(4*w**2)+w2*U/(4*w**3)+3*w1*U1/(4*w**3)-3*w1**2*U/(4*w**4)
    b=-k2/3-M;c=1-b/w**2
    J2=L*w1/w**2+w1**2/(4*w**3)
    J4=L*U1/w**2-2*L*w1*U/w**3+w1*U1/(2*w**3)-3*w1**2*U/(4*w**4)
    R0=w/2;R2=(L**2/w+J2)/4;R4=(U**2/w-L**2*U/w**2+J4)/4
    P0=k2/(6*w);P2=(c*U+L**2/w+J2)/4;P4=(c*V+b*U**2/w**3-L**2*U/w**2+J4)/4
    S0=1/(2*w);S2=-U/(2*w**2)
    expressions={'R0':R0,'R2':R2,'R4':R4,'P0':P0,'P2':P2,'P4':P4,'S0':S0,'S2':S2,
                 'R':R0+R2+R4,'P':P0+P2+P4,'S':S0+S2,'W2':U,'W4':V,'J2':J2,'J4':J4,'omega':w}
    out={}
    for name,expr in expressions.items():
        out[name]={'value':expr.value.actual_derivative(),'delta':expr.delta.actual_derivative()}
        if name in ('S','S0','S2','omega','W2'):
            for n,label in [(1,'first'),(2,'second')]:
                out[name][label]=expr.value.actual_derivative(n)
                out[name]['delta_'+label]=expr.delta.actual_derivative(n)
    return out


def metric_subtractions(eta, k, h_jet):
    """Specialize to H=1, a0=-1/eta, x=r=2 and h derivatives through four."""
    if len(h_jet) < 5: raise ValueError('h derivatives through fourth order are required')
    a=-1/eta
    base=tuple(factorial(n)*a**(n+1) for n in range(5))
    variation=tuple(sum(comb(n,j)*base[n-j]*h_jet[j] for j in range(n+1)) for n in range(5))
    return generic_subtractions(k*k,base,variation,r=2)
