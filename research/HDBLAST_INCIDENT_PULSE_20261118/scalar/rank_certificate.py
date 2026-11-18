"""Exact rational box enclosure for the compact-control endpoint rank test."""
import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from math import ceil, floor
from pathlib import Path

PIN = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"
BOX = {"E":["61.43","61.45"], "Q":["2.740","2.744"],
       "phi":["0.00060","0.00064"], "r":["7.600","7.610"]}


@dataclass(frozen=True)
class Interval:
    lo: Fraction
    hi: Fraction

    def __neg__(self):
        return Interval(-self.hi,-self.lo)

    def __add__(self, other):
        other = interval(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-interval(other))

    def __rsub__(self, other):
        return interval(other) + (-self)

    def __mul__(self, other):
        other = interval(other)
        products = (self.lo*other.lo,self.lo*other.hi,self.hi*other.lo,self.hi*other.hi)
        return Interval(min(products),max(products))

    __rmul__ = __mul__

    def __truediv__(self, denominator):
        if not isinstance(denominator,int) or denominator <= 0:
            raise ValueError("Only positive integer denominators occur here")
        return self*Fraction(1,denominator)


def interval(value):
    if isinstance(value,Interval):
        return value
    return Interval(Fraction(value),Fraction(value))


def jet_columns(E,Q,phi,r):
    """Return K, (d/dZ-L)K and (d/dZ-L)^2K using the reference radial ODE."""
    w = 1-phi+phi*phi*phi/3
    w1, w2, w3 = phi*phi-1, 2*phi, 2
    U = w1*w1/2 - 2*w*w/3
    U1 = w1*w2 - 4*w*w1/3
    U2 = w2*w2+w1*w3-4*(w1*w1+w*w2)/3
    U3 = 3*w2*w3-4*(3*w1*w2+w*w3)/3
    E1 = 2*Q*E
    Q1 = 1-Q*Q-r*r/6-E*U/3
    r1 = E*U1-3*Q*r
    E2 = 2*(Q1*E+Q*E1)
    Q2 = -2*Q*Q1-r*r1/3-(E1*U+E*U1*r)/3
    r2 = E1*U1+E*U2*r-3*(Q1*r+Q*r1)
    r3 = E2*U1+2*E1*U2*r+E*(U3*r*r+U2*r1)-3*(Q2*r+2*Q1*r1+Q*r2)

    def L(vec):
        b,q,v = vec
        return q, -2*E*U*b/3-2*Q*q+2*v, 0

    def dL(vec):
        b,q,_ = vec
        return 0,-2*(E1*U+E*U1*r)*b/3-2*Q1*q,0

    K = (-r/3,(2-Q)*r/3,-r1/3)
    K1 = (-r1/3,((2-Q)*r1-Q1*r)/3,-r2/3)
    K2 = (-r2/3,((2-Q)*r2-2*Q1*r1-Q2*r)/3,-r3/3)
    NK = tuple(a-b for a,b in zip(K1,L(K)))
    NN = tuple(a-b-2*c+d for a,b,c,d in zip(K2,dL(K),L(K1),L(L(K))))
    return K,NK,NN


def determinant(cols):
    (a,d,g),(b,e,h),(c,f,i) = cols
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)


def certificate(raw):
    if hashlib.sha256(raw).hexdigest()!=PIN:
        raise ValueError("Unexpected B1 summary")
    summary = json.loads(raw,parse_float=Fraction)
    phi, rho = summary["phi_b_static"],summary["rho_b"]
    p = summary["params"]
    c,d,det = p["c"],p["d"],p["delta"]
    sigma = 2*(1-phi+phi**3/3)+det*(1+c*phi+d*phi**2/2)
    sigma1 = 2*(phi**2-1)+det*(c+d*phi)
    point = dict(E=rho*rho,Q=rho*sigma/6,phi=phi,r=-rho*sigma1/2)
    box = {k:Interval(Fraction(a),Fraction(b)) for k,(a,b) in BOX.items()}
    if any(not box[k].lo <= v <= box[k].hi for k,v in point.items()):
        raise ValueError("Wall value outside reference box")
    # The static first integral gives an exact normalized wall radius for the
    # specified decimal scalar and tension data. Compare its squares to avoid
    # approximating the positive square root.
    w = 1-phi+phi**3/3
    U = (phi**2-1)**2/2-2*w*w/3
    exact_E = 1/(sigma*sigma/36-sigma1*sigma1/48+U/6)
    if not box["E"].lo <= exact_E <= box["E"].hi or sigma <= 0 or sigma1 >= 0:
        raise ValueError("Static wall radius or orientation outside the certified box")
    for key,square in [("Q",exact_E*sigma*sigma/36),("r",exact_E*sigma1*sigma1/4)]:
        if not 0 < box[key].lo**2 <= square <= box[key].hi**2:
            raise ValueError("Static gradient outside the certified box")
    bounds = determinant(jet_columns(**box))
    value = determinant(jet_columns(**{k:interval(v) for k,v in point.items()}))
    if value.lo != value.hi or not bounds.lo <= value.lo <= bounds.hi < 0:
        raise ValueError("Wall-jet determinant did not certify rank three")
    return {
        "wall_data_box_decimal": BOX,
        "determinant_integer_enclosure": [floor(bounds.lo),ceil(bounds.hi)],
        "point_determinant_rounded": float(value.lo),
        "source_wall_data_inside_box": True,
        "exact_static_normalization_inside_box": True,
        "source_summary_sha256": PIN,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("summary",type=Path)
    arg = p.parse_args()
    print(json.dumps(certificate(arg.summary.read_bytes()),indent=2))
