#!/usr/bin/env python3
"""Independent exact-coefficient identities and directed kernel spot controls."""
from fractions import Fraction as F
import argparse,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'DIFFERENCE_CONTROLS.json');args=ap.parse_args()
def require(ok,msg):
 if not ok:raise ValueError(msg)
N=5
def add(x,y):return [a+b for a,b in zip(x,y)]
def neg(x):return [-a for a in x]
def sub(x,y):return add(x,neg(y))
def mul(x,y):return [sum((x[k]*y[j-k] for k in range(j+1)),F(0)) for j in range(N)]
def inv(x):
 out=[1/x[0]]
 for j in range(1,N):out.append(-out[0]*sum(x[k]*out[j-k] for k in range(1,j+1)))
 return out
x=list(map(F,['7/3','1/5','-2/9','5/21','11/29']));xh=list(map(F,['7/3','1/5','-2/9','5/21','11/29']));xh[0]-=F(1,10**17);xh[3]+=F(3,10**18)
y=list(map(F,['2/7','-3/4','5/9','13/23','-7/11']));yh=y[:];yh[2]+=F(7,10**19)
dx,dy=sub(x,xh),sub(y,yh)
require(sub(mul(x,y),mul(xh,yh))==add(mul(dx,yh),mul(x,dy)),'paired convolution')
require(sub(inv(x),inv(xh))==neg(mul(mul(inv(x),dx),inv(xh))),'paired reciprocal')
require(sub(inv(x),inv(x))==[0]*N,'zero reciprocal difference')
dt=list(map(F,['0','2/3','-7/11','5/13','-3/17']));df=list(map(F,['7/19','-11/23','13/29','-17/31','19/37']))
composed=[F(0)]*N;pw=[F(1)]+[F(0)]*4
for k in range(5):
 composed=add(composed,[df[k]*v/math.factorial(k) for v in pw]);pw=mul(pw,dt)
t1,t2,t3,t4=dt[1:]
explicit=[df[0],df[1]*t1,df[1]*t2+df[2]*t1*t1/2,df[1]*t3+df[2]*t1*t2+df[3]*t1**3/6,df[1]*t4+df[2]*(t1*t3+t2*t2/2)+df[3]*t1*t1*t2/2+df[4]*t1**4/24]
require(composed==explicit,'normalized Bell composition through fourth')
# An independent 512-bit directed scalar interval class makes the spot
# comparisons accurate enough that direct subtraction can serve as a control.
G=1<<512
class B:
 def __init__(self,a,b=None):
  a=F(a);b=a if b is None else F(b);self.a=F((a*G).__floor__(),G);self.b=F((b*G).__ceil__(),G)
 @staticmethod
 def cv(q):return q if isinstance(q,B) else B(q)
 def __add__(self,q):q=B.cv(q);return B(self.a+q.a,self.b+q.b)
 __radd__=__add__
 def __neg__(self):return B(-self.b,-self.a)
 def __sub__(self,q):return self+-B.cv(q)
 def __rsub__(self,q):return B.cv(q)+-self
 def __mul__(self,q):q=B.cv(q);v=[a*b for a in (self.a,self.b) for b in (q.a,q.b)];return B(min(v),max(v))
 __rmul__=__mul__
 def inv(self):require(self.a>0,'positive reciprocal');return B(1/self.b,1/self.a)
 def __truediv__(self,q):return self*B.cv(q).inv()
 def __rtruediv__(self,q):return B.cv(q)*self.inv()
 def __pow__(self,n):
  out=B(1)
  for _ in range(n):out=out*self
  return out
def expref(t):
 h=0
 while t>1:t/=2;h+=1
 upper=sum(((-t)**j/math.factorial(j) for j in range(81)),F(0));lower=upper+(-t)**81/math.factorial(81);r=B(lower,upper)
 for _ in range(h):r=r*r
 return r
def kernelref(t,n):
 r=1/B(1+t/n)
 for _ in range(40):r=r*r
 return r
rows=[];n=1<<40
for t in map(F,['0','1/1000000','1','7','60']):
 R=kernelref(t,n);E=expref(t);A=F(1)
 for k in range(5):
  if k:A*=1+F(k-1,n)
  actual=(-1)**k*E;nominal=(-1)**k*A*R/(1+t/n)**k;direct=actual-nominal
  # Expanded difference, evaluated by exact scalar fractions plus a
  # separately bounded nonpositive D0.
  p=sum((F(math.comb(k,j),n**j)*t**j for j in range(k+1)),F(0))-A
  d0=B(-(R*t*t/(2*n)).b,0);bound=(-1)**k*(d0+R*p/(1+t/n)**k)
  # At exact algebraic zeros the reference enclosure width can straddle
  # zero; the identity itself proves that point without reference roundoff.
  if t==0 and k<2:require(bound.a==bound.b==0,'exact zero kernel difference')
  else:require(bound.a<=direct.a<=direct.b<=bound.b,'kernel derivative containment')
  rows.append({'t':str(t),'k':k,'direct_interval':[str(direct.a),str(direct.b)],'proved_difference_bound':[str(bound.a),str(bound.b)]})
report={'status':'PASS_EXACT_PAIRED_ALGEBRA_AND_INDEPENDENT_KERNEL_CONTROLS','paired_convolution_through_order':4,'paired_reciprocal_through_order':4,'independent_Bell_composition_through_order':4,'scalar_kernel_cases':len(rows),'control_grid_bits':512,'scope':'Independent spot controls supplement the analytic all-cell proof; not a replacement for it.','rows':rows}
args.output.write_text(json.dumps(report,indent=2)+'\n');print(report['status'])
