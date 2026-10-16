"""Integer implementation of the same outward 2^112 interval operations.

Endpoint integers count grid ticks. Every operation has exactly the same
rounding order as interval_core.I. This avoids repeated Fraction gcd work
inside the large algebraically collected source maps. Input/output endpoints
remain exact Fraction values; numerical decisions use no binary floats.
"""
from fractions import Fraction as F
GRID=1<<112

def ceil_ratio(n,d):return -((-n)//d)

class D:
 __slots__=('l','h')
 def __init__(self,lo=0,hi=None):
  if hi is None and isinstance(lo,D):self.l,self.h=lo.l,lo.h;return
  if hi is None and hasattr(lo,'lo') and hasattr(lo,'hi'):lo,hi=lo.lo,lo.hi
  lo=F(lo);hi=lo if hi is None else F(hi)
  if lo>hi:raise ValueError('dyadic interval order')
  self.l=(lo*GRID).__floor__();self.h=(hi*GRID).__ceil__()
 @classmethod
 def ticks(cls,l,h):
  out=object.__new__(cls);out.l=l;out.h=h;return out
 @staticmethod
 def cv(x):return x if isinstance(x,D) else D(x)
 @property
 def lo(self):return F(self.l,GRID)
 @property
 def hi(self):return F(self.h,GRID)
 def __add__(self,b):b=D.cv(b);return D.ticks(self.l+b.l,self.h+b.h)
 __radd__=__add__
 def __neg__(self):return D.ticks(-self.h,-self.l)
 def __sub__(self,b):return self+-D.cv(b)
 def __rsub__(self,b):return D.cv(b)+-self
 def __mul__(self,b):
  b=D.cv(b);v=(self.l*b.l,self.l*b.h,self.h*b.l,self.h*b.h)
  return D.ticks(min(v)//GRID,ceil_ratio(max(v),GRID))
 __rmul__=__mul__
 def inv(self):
  if not (self.l>0 or self.h<0):raise ValueError('zero denominator')
  return D.ticks((GRID*GRID)//self.h,ceil_ratio(GRID*GRID,self.l))
 def __truediv__(self,b):return self*D.cv(b).inv()
 def __rtruediv__(self,b):return D.cv(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  if n==0:return D(1)
  if n%2==0:
   a=self**(n//2);lower=0 if a.l<=0<=a.h else min(a.l*a.l,a.h*a.h)
   return D.ticks(lower//GRID,ceil_ratio(max(a.l*a.l,a.h*a.h),GRID))
  return self*self**(n-1)
 def mag(self):return F(max(abs(self.l),abs(self.h)),GRID)
 def row(self):return [str(self.lo),str(self.hi)]

def verify_equivalence(reference):
 """Deterministic rational and signed-endpoint differential controls."""
 vals=[(-7,-3),(-3,0),(-2,5),(0,0),(0,7),(3,8),(F(1,3),F(7,11)),(-F(7,11),-F(1,3))]
 checks=0
 for lo,hi in vals:
  a=reference(lo,hi);b=D(lo,hi)
  if a.row()!=b.row():raise ValueError('dyadic constructor equivalence')
  for n in range(7):
   if (a**n).row()!=(b**n).row():raise ValueError('dyadic integer-power equivalence')
   checks+=1
  if lo>0 or hi<0:
   if a.inv().row()!=b.inv().row():raise ValueError('dyadic reciprocal equivalence')
   checks+=1
  for lo2,hi2 in vals:
   c=reference(lo2,hi2);d=D(lo2,hi2)
   for expected,actual in ((a+c,b+d),(a-c,b-d),(a*c,b*d)):
    if expected.row()!=actual.row():raise ValueError('dyadic arithmetic equivalence')
    checks+=1
   if lo2>0 or hi2<0:
    if (a/c).row()!=(b/d).row():raise ValueError('dyadic quotient equivalence')
    checks+=1
 return checks
