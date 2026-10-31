from fractions import Fraction as F

class P:
 def __init__(self, t=0):
  if isinstance(t,P):self.t=dict(t.t)
  elif isinstance(t,dict):self.t={k:F(v) for k,v in t.items() if v}
  else:self.t={():F(t)} if t else {}
 @classmethod
 def var(cls,n):return cls({((n,1),):F(1)})
 @classmethod
 def cv(cls,x):return x if isinstance(x,cls) else cls(x)
 def __add__(self,b):
  t=dict(self.t)
  for k,v in self.cv(b).t.items():t[k]=t.get(k,F(0))+v
  return P(t)
 __radd__=__add__
 def __neg__(self):return P({k:-v for k,v in self.t.items()})
 def __sub__(self,b):return self+-self.cv(b)
 def __rsub__(self,b):return self.cv(b)+-self
 def __mul__(self,b):
  t={}
  for k,v in self.t.items():
   for l,w in self.cv(b).t.items():
    z=dict(k)
    for n,e in l:z[n]=z.get(n,0)+e
    z=tuple(sorted((n,e) for n,e in z.items() if e))
    t[z]=t.get(z,F(0))+v*w
  return P(t)
 __rmul__=__mul__
 def __pow__(self,n):
  if n<0:
   if len(self.t)!=1:raise ValueError('NONMONOMIAL_DENOMINATOR')
   (k,v),=self.t.items();return P({tuple((s,e*n) for s,e in k):v**n})
  out=P(1)
  while n:
   if n&1:out=out*self
   self=self*self;n//=2
  return out
 def __truediv__(self,b):return self*self.cv(b)**-1
 def __rtruediv__(self,b):return self.cv(b)*self**-1
 def __bool__(self):return bool(self.t)
 def __eq__(self,b):return self.t==self.cv(b).t
 def diff(self,n):
  t={}
  for k,v in self.t.items():
   z=dict(k);e=z.get(n,0)
   if e:
    z[n]=e-1;t[tuple(sorted((s,p) for s,p in z.items() if p))]=v*e
  return P(t)
 def names(self):return sorted({n for k in self.t for n,e in k})
 def encode(self):return [[[[n,e] for n,e in k],str(v)] for k,v in sorted(self.t.items())]
 def __repr__(self):return str(self.encode())

def dt(x):
 out=P()
 rules={s:P.var(s+'t') for s in ('h','X','V','L','d','v')}
 rules.update({'c':P.var('c')*P.var('theta'),'theta':1-P.var('theta')**2})
 for n in P.cv(x).names():
  out+=P.cv(x).diff(n)*rules.get(n,P())
 return out

def encoded(terms,values):
 out=P()
 for row in terms:
  q=P(F(row['coefficient']))
  for n,p in row['powers'].items():q=q*values[n]**p
  out+=q
 return out

def shell_values(report,normal_scalar=False):
 vals={n:P.var(n) for n in report['input_symbols']}
 g={n:P.var(n) for n in ('h','X','V','L','d','v','rho','c','theta','a','p','sigma2')}
 h,X,V,L,d,v,rho,c,theta,a,p,sigma2=[g[n] for n in g]
 sig1=-2*p;beta=a*h+sig1*v/3
 vals.update(A=h,B=P(),C=P(),D=-2*X,E=V,L=2*L,M=2*d,s=v)
 for key,x in [('A',h),('D',-2*X),('E',V),('L',2*L),('M',2*d),('s',v)]:
  vals[key+'t']=dt(x)
 # B=C zero as physical shell traces: only tangential derivatives vanish.
 for key in ('B','C'):
  vals[key+'t']=P();vals[key+'tt']=P()
 vals.update(Dy=-beta,Ly=beta,Ey=P(),My=P())
 for k in ('D','L','E','M'):vals[k+'yt']=dt(vals[k+'y'])
 if normal_scalar:vals['sy']=p*h/2-sigma2*v/2
 return vals,g,beta
