#!/usr/bin/env python3
"""Exact rational static jets on the retained original packet annulus.

This is interval automatic differentiation at every possible annulus point,
not a sampled trajectory or a convergent Taylor approximation claim.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
PINS={
 'ACTUAL_POTENTIAL_JET_BOUNDS.json':'fd93a0de3d321ff5182115138eed2725f3c327083cf5148e3efc1aca3659d5c0',
 'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
 'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md':'0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
}
ORDER=4


def require(ok,text):
 if not ok:raise ValueError(text)


class I:
 def __init__(self,lo,hi=None):
  self.lo,self.hi=F(lo),F(lo if hi is None else hi)
  require(self.lo<=self.hi,'reversed rational interval')
 @staticmethod
 def box(x):return x if isinstance(x,I) else I(x)
 def __add__(self,y):
  y=self.box(y);return I(self.lo+y.lo,self.hi+y.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,y):return self+-self.box(y)
 def __rsub__(self,y):return self.box(y)+-self
 def __mul__(self,y):
  y=self.box(y);v=[a*b for a in (self.lo,self.hi) for b in (y.lo,y.hi)]
  return I(min(v),max(v))
 __rmul__=__mul__
 def recip(self):
  require(not self.lo<=0<=self.hi,'unsafe interval denominator')
  return I(1/self.hi,1/self.lo)
 def __truediv__(self,y):return self*self.box(y).recip()
 def __rtruediv__(self,y):return self.box(y)*self.recip()
 def __pow__(self,n):
  require(isinstance(n,int),'integer interval exponent required')
  if n<0:return self.recip()**(-n)
  if n==0:return I(1)
  if n%2==0:
   h=self**(n//2)
   return I(0 if h.lo<=0<=h.hi else min(h.lo*h.lo,h.hi*h.hi),max(h.lo*h.lo,h.hi*h.hi))
  return self*self**(n-1)
 def intersect(self,y):return I(max(self.lo,y.lo),min(self.hi,y.hi))
 def mag(self):return max(abs(self.lo),abs(self.hi))
 def record(self):return {'exact':[str(self.lo),str(self.hi)],'display_only':[float(self.lo),float(self.hi)]}


class P:
 """Independent formal total-derivative algebra on rho,eta,p,a,alpha."""
 n=5
 def __init__(self,d=None):self.d={m:F(c) for m,c in (d or {}).items() if c}
 @classmethod
 def box(cls,x):return x if isinstance(x,cls) else cls({(0,)*cls.n:F(x)})
 @classmethod
 def var(cls,i):
  e=[0]*cls.n;e[i]=1;return cls({tuple(e):F(1)})
 def __add__(self,y):
  d=self.d.copy()
  for m,c in self.box(y).d.items():d[m]=d.get(m,F(0))+c
  return P(d)
 __radd__=__add__
 def __neg__(self):return P({m:-c for m,c in self.d.items()})
 def __sub__(self,y):return self+-self.box(y)
 def __rsub__(self,y):return self.box(y)+-self
 def __mul__(self,y):
  d={}
  for m,c in self.d.items():
   for n,b in self.box(y).d.items():
    v=tuple(i+j for i,j in zip(m,n));d[v]=d.get(v,F(0))+c*b
  return P(d)
 __rmul__=__mul__
 def __pow__(self,n):
  if n<0:
   require(len(self.d)==1,'formal inverse requires a monomial')
   (m,c),=self.d.items();return P({tuple(i*n for i in m):c**n})
  x=P.box(1)
  for _ in range(n):x=x*self
  return x
 def __truediv__(self,y):return self*self.box(y)**-1
 def __rtruediv__(self,y):return self.box(y)*self**-1
 def diff(self,i):
  d={}
  for m,c in self.d.items():
   if m[i]:
    e=list(m);e[i]-=1;d[tuple(e)]=c*m[i]
  return P(d)


class J:
 """Normalized finite Taylor jets; coefficients may be intervals or polynomials."""
 def __init__(self,c):
  require(all(not isinstance(v,float) for v in c),'floating coefficient forbidden')
  self.c=[F(v) if isinstance(v,(int,F)) else v for v in c]
 @property
 def n(self):return len(self.c)
 def box(self,x):return x if isinstance(x,J) else J([x]+[0]*(self.n-1))
 def __add__(self,y):
  y=self.box(y);return J([self.c[i]+y.c[i] for i in range(min(self.n,y.n))])
 __radd__=__add__
 def __neg__(self):return J([-v for v in self.c])
 def __sub__(self,y):return self+-self.box(y)
 def __rsub__(self,y):return self.box(y)+-self
 def __mul__(self,y):
  y=self.box(y);return J([sum((self.c[j]*y.c[i-j] for j in range(i+1)),0) for i in range(min(self.n,y.n))])
 __rmul__=__mul__
 def recip(self):
  out=[1/self.c[0]]
  for n in range(1,self.n):out.append(-out[0]*sum((self.c[j]*out[n-j] for j in range(1,n+1)),0))
  return J(out)
 def __truediv__(self,y):return self*self.box(y).recip()
 def __rtruediv__(self,y):return self.box(y)*self.recip()
 def __pow__(self,n):
  if n<0:return self.recip()**(-n)
  out=self.box(1)
  for _ in range(n):out=out*self
  return out
 def dz(self):return J([(i+1)*self.c[i+1] for i in range(self.n-1)])


U=[F(-2,27),F(0),F(14,9),F(50,27),F(-1,6),F(-4,9),F(-2,27)]


def derivative(poly):return [i*poly[i] for i in range(1,len(poly))]


def evaluate(poly,x):
 v=0
 for c in reversed(poly):v=v*x+c
 return v


UP=[U]
for _ in range(4):UP.append(derivative(UP[-1]))


def seed(r,e,p,a,alpha):
 h=evaluate(UP[1],e)/p
 return {'rho':r,'eta':e,'p':p,'a':a,'rho_inv':1/r,'h':h,'alpha':alpha,
         'C':alpha/r**2,'B':3*alpha/(r**3*p)}


def build(initial,order,injection='none'):
 state={k:J([v]+[0]*order) for k,v in initial.items()}
 for n in range(order):
  r,e,p,a,b,h,alpha,C,B=[state[k] for k in ('rho','eta','p','a','rho_inv','h','alpha','C','B')]
  k=r*(F(5,2)*a-h);L=r*(F(9,2)*a-h)
  rhoprime=-r**2*a
  if injection=='flip_rho_ODE':rhoprime=-rhoprime
  Bfactor=F(3,4) if injection=='halve_B_log_derivative' else F(3,2)
  rhs={'rho':rhoprime,'eta':-r*p,'p':-r*(evaluate(UP[1],e)-4*a*p),
       'a':b+r*p**2/3,'rho_inv':a,'h':-r*(evaluate(UP[2],e)-h**2+4*a*h),
       'alpha':k*alpha,'C':L*C,'B':Bfactor*r*a*B}
  divisor=1 if injection=='wrong_jet_normalization' else n+1
  for name in state:state[name].c[n+1]=rhs[name].c[n]/divisor
 r,e,p,a,h=[state[k] for k in ('rho','eta','p','a','h')]
 out=dict(state)
 out.update({'p_inv':1/p,'k':r*(F(5,2)*a-h)})
 for j,poly in enumerate(UP):out['U_'+str(j)]=evaluate(poly,e)
 return out


def formal_checks(injection):
 r,e,p,a,alpha=[P.var(i) for i in range(5)]
 initial=seed(r,e,p,a,alpha)
 jets=build(initial,ORDER,injection)
 velocity=(-r*r*a,-r*p,-r*(evaluate(UP[1],e)-4*a*p),1/r+r*p*p/3,
           r*(F(5,2)*a-evaluate(UP[1],e)/p)*alpha)
 def D(f):return sum((f.diff(i)*v for i,v in enumerate(velocity)),P())
 records={}
 for name,series in jets.items():
  value=series.c[0]
  for n in range(1,ORDER+1):
   value=D(value)
   require(not (factorial(n)*series.c[n]-value).d,'formal jet identity failed: '+name+' order '+str(n))
  records[name]='PASS_ORDERS_1_TO_'+str(ORDER)
 # Product-rule conversion is exact in the same independent algebra.
 f=alpha/r**2;b=1/r
 require(not ((-b*D(-b*D(f)))-(b*b*D(D(f))+a*b*D(f))).d,'formal y conversion identity failed')
 return records


def load_inputs(injection):
 for name,pin in PINS.items():
  raw=(ROOT/'inputs'/name).read_bytes()
  if injection=='input_tamper' and name=='ACTUAL_POTENTIAL_JET_BOUNDS.json':raw+=b' '
  require(sha256(raw).hexdigest()==pin,'actual retained input pin mismatch: '+name)
 bg=json.loads((ROOT/'inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json').read_text())
 shell=json.loads((ROOT/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
 require(bg['actual_model']['U_coefficients']==list(map(str,U)),'actual potential polynomial differs')
 require(shell['parameters_exact']=={'c':bg['actual_model']['c'],'delta':bg['actual_model']['delta']},'actual shell branch differs')
 return bg,shell


def run(injection='none'):
 bg,shell=load_inputs(injection)
 ell=F(1,1000000);left,right=ell,2*ell
 outer={k:I(*bg['background_intervals'][k]['exact']) for k in ('rho','eta','p','a')}
 r,e,p,a=[outer[k] for k in ('rho','eta','p','a')]
 rhs={'rho':-r**2*a,'eta':-r*p,'p':-r*(evaluate(UP[1],e)-4*a*p),'a':1/r+r*p**2/3}
 annulus={}
 refined={}
 initial_names={'rho':'rho_b','eta':'eta_b','p':'p_b','a':'a_b'}
 for name in outer:
  require(rhs[name].record()['exact']==bg['background_intervals'][name+'_z']['exact'],'actual outer ODE image differs: '+name)
  x0=I(*shell['boxes'][initial_names[name]]['exact'])
  tube=x0+I(left,right)*rhs[name]
  annulus[name]=outer[name].intersect(tube)
  refined[name]={'outer':outer[name].record(),'initial_shell':x0.record(),
                 'outer_ODE_image':rhs[name].record(),'annulus_integral_image':tube.record(),
                 'accepted_intersection':annulus[name].record()}
 # This global coefficient is needed along the integral from z=0.
 hglobal=evaluate(UP[1],e)/p;kglobal=r*(F(5,2)*a-hglobal)
 require(kglobal.record()['exact']==bg['reconstruction_coefficients']['alpha_log_z']['exact'],'retained global k identity differs')
 require(kglobal.hi<0,'alpha sign certificate failed')
 K,m=-kglobal.lo,-kglobal.hi
 alpha=I(1-K*right,1/(1+m*left))
 require(0<alpha.lo<alpha.hi<1,'annulus alpha enclosure failed')
 r,e,p,a=[annulus[k] for k in ('rho','eta','p','a')]
 if injection=='unsafe_p_denominator':p=I(-1,1)
 require(r.lo>0 and p.hi<0,'actual annulus denominator signs failed')
 jets=build(seed(r,e,p,a,alpha),ORDER,injection)
 E=I(1-right,1/(1+left))
 jets['exp_minus_z']=J([(-1)**n*E/factorial(n) for n in range(ORDER+1)])
 formal=formal_checks(injection)
 yjets={}
 for name,series in jets.items():
  current=series;values=[current.c[0]]
  for n in range(1,ORDER+1):
   current=-jets['rho_inv']*current.dz()
   values.append(current.c[0]/factorial(n))
  yjets[name]=values
 def caps(values):
  return [str((value.mag()*factorial(n)).numerator//(value.mag()*factorial(n)).denominator+1)
          for n,value in enumerate(values)]
 report={
  'status':'PASS_EXACT_ACTUAL_ANNULUS_STATIC_COEFFICIENT_JETS',
  'domain':{'t':'0','ell':str(ell),'z_min':str(left),'z_max':str(right),'coverage':'Every point of the original packet annulus; no samples.'},
  'order':ORDER,'normalization':'coefficient j = D_z^j(value)/j! or D_y^j(value)/j!, as named',
  'input_sha256':PINS,'actual_model':bg['actual_model'],
  'annulus_refinement':refined,
  'alpha_enclosure':{'exact':alpha.record()['exact'],'global_log_derivative':kglobal.record(),
                     'proof':'alpha(0)=1, k in [-K,-m]<0 on [0,2ell]; 1-2ell K <= alpha(z) <= 1/(1+ell m).'},
  'z_jets':{name:[I.box(v).record() for v in series.c] for name,series in jets.items()},
  'y_jets':{name:[I.box(v).record() for v in values] for name,values in yjets.items()},
  'strict_abs_z_derivative_caps':{name:caps(series.c) for name,series in jets.items()},
  'strict_abs_y_derivative_caps':{name:caps(values) for name,values in yjets.items()},
  'formal_total_derivative_checks':formal,
  'formal_y_conversion':'PASS: D_y^2 f = rho^-2 D_z^2 f + (a/rho)D_z f',
  'coefficient_definitions':{'C':'alpha/rho^2','B':'3alpha/(rho^3 p)','h':'U_1/p','k':'rho(5a/2-h)',
                             'alpha':'(rho/rho_b)^(3/2) p/p_b','exp_minus_z':'exp(-z)'},
  'reconstruction':'u=C q; u_t=C exp(-z) q; v=B(q_z+kq); v_t=B exp(-z)(q_z+(k-1)q).',
  'scope':'Conditional actual retained static background coefficient derivative enclosures only. No second initial norms, gauge completion, evolution, exact nonlinear family or publication.'}
 return report


def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--output',type=Path,default=ROOT/'ACTUAL_ANNULUS_STATIC_JETS.json')
 ap.add_argument('--inject-error',choices=['none','input_tamper','wrong_jet_normalization','flip_rho_ODE','halve_B_log_derivative','unsafe_p_denominator'],default='none')
 args=ap.parse_args()
 report=run(args.inject_error)
 args.output.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({'status':report['status'],'order':ORDER,'output':str(args.output),
                   'formal_jet_families':len(report['formal_total_derivative_checks'])},indent=2))


if __name__=='__main__':main()
