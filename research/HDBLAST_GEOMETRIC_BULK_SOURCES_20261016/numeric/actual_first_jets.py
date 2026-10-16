#!/usr/bin/env python3
"""Certified actual first common-gauge frame-sector 2-jets.

All decisions are rational.  The intervals enclose derivatives of the actual
solutions, using authenticated differentiated C3 error certificates; they are
not derivatives of independently selected value intervals.  ``actual_jet`` is
valid on every rectangular cell contained in 0<=x<=3, r>=0, x+r<=8.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction as F
from pathlib import Path

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('interval_core',ROOT/'interval_core.py')
core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
I=core.I; ELL=core.ELL; need=core.need
N=3
KEYS=tuple((i,j) for i in range(N+1) for j in range(N+1-i))

def hull(*xs):
 xs=[I.cv(x) for x in xs]
 return I(min(x.lo for x in xs),max(x.hi for x in xs))

class T:
 """Bivariate normalized Taylor jet; only coefficients through total 3."""
 def __init__(self,c=None):
  c={} if c is None else c
  self.c={k:I.cv(c.get(k,0)) for k in KEYS}
 @staticmethod
 def cv(x):return x if isinstance(x,T) else T({(0,0):x})
 def __add__(self,b):
  b=T.cv(b);return T({k:self.c[k]+b.c[k] for k in KEYS})
 __radd__=__add__
 def __neg__(self):return T({k:-self.c[k] for k in KEYS})
 def __sub__(self,b):return self+-T.cv(b)
 def __rsub__(self,b):return T.cv(b)+-self
 def __mul__(self,b):
  b=T.cv(b)
  return T({(i,j):sum((self.c[a,b0]*b.c[i-a,j-b0] for a in range(i+1) for b0 in range(j+1)),I(0)) for i,j in KEYS})
 __rmul__=__mul__
 def inv(self):
  out={(0,0):self.c[0,0].inv()}
  for n in range(1,N+1):
   for i in range(n+1):
    j=n-i
    out[i,j]=-out[0,0]*sum((self.c[a,b]*out[i-a,j-b] for a in range(i+1) for b in range(j+1) if a+b),I(0))
  return T(out)
 def __truediv__(self,b):return self*T.cv(b).inv()
 def __rtruediv__(self,b):return T.cv(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  out=T.cv(1)
  for _ in range(n):out=out*self
  return out
 def dx(self):return T({(i,j):(i+1)*self.c.get((i+1,j),I(0)) for i,j in KEYS})
 def dr(self):return T({(i,j):(j+1)*self.c.get((i,j+1),I(0)) for i,j in KEYS})
 def val(self):return self.c[0,0]
 def derivative(self,i,j):return math.factorial(i)*math.factorial(j)*self.c[i,j]

_CACHE=None
_EXTERIOR=None
def inputs():
 global _CACHE
 if _CACHE is not None:return _CACHE
 pins=json.loads((ROOT/'PREDECESSOR_PINS.json').read_text())
 for path,pin in pins.items():
  raw=(ROOT/path).read_bytes()
  need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'predecessor authentication: '+path)
 data,_=core.load()
 errors=json.loads((ROOT/'inputs/APPROXIMATION_ERRORS_NORMAL.json').read_text())
 review=json.loads((ROOT/'inputs/C3_ERROR_FINAL_NORMAL.json').read_text())
 need(errors['status']=='PASS_ACTUAL_MASTER_AND_GAUGE_CONTINUOUS_SCALED_C3_APPROXIMATION_ERRORS','actual differentiated error status')
 need(review['status']=='PASS_INDEPENDENT_C3_ERROR_REVIEW','independent inherited error review')
 need(review['producer_report_sha256']==pins['inputs/APPROXIMATION_ERRORS_NORMAL.json']['sha256'],'review binds inherited error certificate')
 need(review['producer_source_sha256']==pins['check_inherited_errors.py']['sha256'],'review binds inherited error producer')
 need(errors['input_sha256']==review['input_sha256'],'review coefficient/master dependency')
 _CACHE=(data,errors,pins)
 return _CACHE

def exterior_caps():
 global _EXTERIOR
 if _EXTERIOR is None:
  spec=importlib.util.spec_from_file_location('exterior_errors',ROOT/'exterior_error_caps.py')
  module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
  # The same implementation with the retained lift restored must exactly
  # reproduce all inherited derivative errors before any local use.
  certificate=module.run()
  _EXTERIOR=([F(v) for v in certificate['master_actual_minus_comparison_by_total_order']],
             {k:list(map(F,v)) for k,v in certificate['gauge_actual_minus_comparison_mixed'].items()})
 return _EXTERIOR

def qpositive(r,order,data):
 """Q derivatives on a positive interval, including analytic flat tails."""
 out=[];edge=F(1,64)
 if r.lo<=1 or r.hi>=2:out.append([I(0)]*(order+1))
 for lo,hi in ((max(r.lo,F(1)),min(r.hi,F(1)+edge)),(max(r.lo,F(2)-edge),min(r.hi,F(2)))):
  if lo<=hi:
   # |Q^(n)| <= s^(2n) exp(4-s), s>=4096/63>2n.
   # The common tail threshold avoids losing the tiny exponential to the
   # fixed output grid before multiplying a very large inverse distance.
   near=hi if hi<=1+edge else lo
   if near==1 or near==2:out.append([I(0)]*(order+1));continue
   s0=1/(edge*(1-edge));ex=core.expminus(s0-4)
   need(s0>=2*order,'flat-tail monotonicity')
   out.append([I(0,s0**(2*n)*ex.hi) if n==0 else I(-s0**(2*n)*ex.hi,s0**(2*n)*ex.hi) for n in range(order+1)])
 lo=max(r.lo,F(1)+edge);hi=min(r.hi,F(2)-edge)
 if lo<=hi:
  rr=I(lo,hi);s=1/(F(1,4)-(rr-F(3,2))**2);v=3-2*rr;ex=core.expbox(s-4)
  polys=data['bump']['polynomial_coefficients_ascending']
  out.append([core.centered_poly(list(map(F,polys[n])),s)*ex*(v if n%2 else 1) for n in range(order+1)])
 if not out:return [I(0)]*(order+1)
 bounds=list(map(F,data['bump']['dimensionless_Q_derivative_sup']))
 result=[]
 for n in range(order+1):
  value=hull(*(row[n] for row in out));cap=bounds[n]
  result.append(I(max(value.lo,F(0) if n==0 else -cap),min(value.hi,cap)))
 return result

def qodd(s,order,data):
 """Derivatives of the smooth odd compact extension Qo."""
 out=[]
 if s.lo<=0<=s.hi:out.append([I(0)]*(order+1))
 if s.hi>=0:out.append(qpositive(I(max(s.lo,F(0)),s.hi),order,data))
 if s.lo<=0:
  z=qpositive(I(max(-s.hi,F(0)),-s.lo),order,data)
  out.append([(-1)**(n+1)*z[n] for n in range(order+1)])
 return [hull(*(row[n] for row in out)) for n in range(order+1)]

def comparisons(x,r,data):
 """Explicit Wbar,P jets, with Pe-rQo cancellations done first."""
 qp=qodd(r+x,N,data);qm=qodd(r-x,N,data);c=6*ELL/130
 w={};p={};sm=r-x
 # Negative and positive characteristic branches have an identically zero
 # neighbourhood at sm=0.  Taking their hull across that gap is sound.
 branches=[]
 if sm.hi>=0:branches.append((False,qodd(I(max(sm.lo,F(0)),sm.hi),N,data)))
 if sm.lo<=0:branches.append((True,qodd(I(sm.lo,min(sm.hi,F(0))),N,data)))
 for i,j in KEYS:
  n=i+j;den=math.factorial(i)*math.factorial(j)
  w[i,j]=(qp[n]+(-1)**i*qm[n])/(2*den)
  plus=x*qp[n]+(i*qp[n-1] if i else 0)
  pieces=[]
  for negative,q in branches:
   lin=x-2*r if negative else -x
   lx=1 if negative else -1;lr=-2 if negative else 0
   v=lin*(-1)**i*q[n]
   if i:v+=i*lx*(-1)**(i-1)*q[n-1]
   if j:v+=j*lr*(-1)**i*q[n-1]
   pieces.append(v)
  p[i,j]=c*(plus+hull(*pieces))/(4*den)
 return T(w),T(p)

def temporal(x):
 """c=cosh(ell*x), H=tanh(ell*x), v=c^-3/2, scaled jets."""
 t=ELL*x;en=core.expbox(t);ep=1/en
 cc=(ep+en)/2;ss=(ep-en)/2;hh=ss/cc;jj=1-hh**2
 c=T({(n,0):ELL**n*(cc if n%2==0 else ss)/math.factorial(n) for n in range(4)})
 h=T({(0,0):hh,(1,0):ELL*jj,(2,0):-ELL**2*hh*jj,(3,0):-ELL**3*jj*(1-3*hh**2)/3})
 # c^(-3/2) uses monotonic rational square-root brackets; derivatives follow
 # v_t=-(3/2)H v, independently of any value-interval differentiation.
 rootlo=F(math.isqrt((cc.lo*core.GRID**2).__floor__()),core.GRID)
 roothi=F(math.isqrt((cc.hi*core.GRID**2).__ceil__())+1,core.GRID)
 need(rootlo**2<=cc.lo and roothi**2>=cc.hi,'rational cosh square-root brackets')
 v0=1/(cc*I(rootlo,roothi));v=T({(0,0):v0})
 rhs=-F(3,2)*ELL*h
 for n in range(3):v.c[n+1,0]=(rhs*v).c[n,0]/(n+1)
 return c,h,v

def radial(j):return T({(0,n):ELL**n*j.c[n] for n in range(4)})

def physical_derivatives(j,rho,a):
 L=1/(ELL*rho.val())
 return {'':j.val(),'y':-L*j.derivative(0,1),'t':j.derivative(1,0)/ELL,
         'yy':L**2*j.derivative(0,2)+a.val()*L*j.derivative(0,1),
         'yt':-L*j.derivative(1,1)/ELL,'tt':j.derivative(2,0)/ELL**2}

def frame_amplitudes(W,U,rho,p,a,R,B,k,c,H,vartheta,injection='none'):
 """Registered Eq.(7) and scalar pullback, kept as one correlated algebra."""
 f,b,d=U;u=R*vartheta*W;v=B*vartheta*(W.dr()/ELL+k*W)
 amps={'A':-4*u-2*f.dr()/(ELL*rho),
       'B':f.dx()/(ELL*rho)+b.dr()/ELL,
       'C':f/(rho*c)-c*d.dr()/ELL,
       'D':-2*(u+a*f+b.dx()/ELL),
       'E':c*d.dx()/ELL-b/c,
       'L':2*(u+a*f+H*b),'M':2*d}
 if injection=='drop_radial_rho':amps['A']=-4*u-2*f.dr()/ELL
 scalar=v-p*f if injection=='wrong_scalar_shift' else v+p*f
 return amps,scalar

def actual_jet(x,r,mode='actual',injection='none',clip_causal=False):
 """Return (m,s,b), intervals in ordinary physical y,t derivatives.

 m keys are A,B,C,D,E,L,M plus suffix '',y,t,yy,yt,tt.  They multiply
 Y,Y,DY,Y,DY,(Y*kappa),HessY in registered orthonormal frame order.
 s keys are s,sy,st,syy,syt,stt.  Mode comparator is diagnostic only.
 With clip_causal=True the assertion is only on the nonempty intersection
 of the supplied rectangular box and x+r<=8.  Evaluating its comparison
 function on the bounding box only widens the valid enclosure on that
 intersection; no actual solution outside the causal triangle is assumed.
 """
 x=I(x.lo,x.hi) if hasattr(x,'lo') and hasattr(x,'hi') else I.cv(x)
 r=I(r.lo,r.hi) if hasattr(r,'lo') and hasattr(r,'hi') else I.cv(r)
 need(0<=x.lo<=x.hi<=3 and 0<=r.lo<=r.hi<=8 and x.lo+r.lo<=8,'nonempty cell in controlled time/collar')
 need(clip_causal or x.hi+r.hi<=8,'whole cell inside causal triangle')
 need(mode in ('actual','comparator'),'jet mode')
 data,errors,_=inputs();W,P=comparisons(x,r,data);U=[P,T(),T()]
 if mode=='actual':
  ew=[F(v['exact_upper']) for v in errors['master_error_budget_by_scaled_total_order']['actual_minus_comparison']]
  eu_all={k:[F(v['exact_upper']) for v in values] for k,values in errors['actual_minus_gauge_comparison_mixed_scaled_derivative_caps'].items()}
  if r.lo>=x.hi:ew,eu_all=exterior_caps()
  for i,j in KEYS:
   den=math.factorial(i)*math.factorial(j);W.c[i,j]+=I(-ew[i+j],ew[i+j])/den
   eu=eu_all[f'{i},{j}']
   for channel in range(3):
    cap=eu[channel];U[channel].c[i,j]+=I(-cap,cap)/den
  # Tangential derivatives of the physical mixed rows, as identities of
  # actual smooth solutions; no odd/even extension of coefficients is used.
  if r.lo==r.hi==0:
   for i in range(N+1):U[0].c[i,0]=I(0)
   for i in range(N):
    U[1].c[i,1]=I(0);U[2].c[i,1]=I(0)
 rr,eta,pp,aa,hh,RR,BB,kk,_,_=core.background(r,data,'none')
 rho,p,a,R,B,k=map(radial,(rr,pp,aa,RR,BB,kk))
 c,H,vartheta=temporal(x)
 amps,scalar=frame_amplitudes(W,U,rho,p,a,R,B,k,c,H,vartheta,injection)
 m={key+suffix:value for key,amp in amps.items() for suffix,value in physical_derivatives(amp,rho,a).items()}
 s={'s'+suffix:value for suffix,value in physical_derivatives(scalar,rho,a).items()}
 ay=-1/rho.val()**2-p.val()**2/3;py=p.val()*(radial(hh).val()-4*a.val())
 if injection=='drop_radial_conversion':
  for key,amp in amps.items():m[key+'yy']=(1/(ELL*rho.val()))**2*amp.derivative(0,2)
  s['syy']=(1/(ELL*rho.val()))**2*scalar.derivative(0,2)
 bg={'rho':rho.val(),'c':c.val(),'theta':H.val(),'a':a.val(),'ay':ay,'p':p.val(),'py':py,
     'Vp':core.poly(core.UP[1],eta.val()),'Vpp':core.poly(core.UP[2],eta.val()),'Vppp':core.poly(core.UP[3],eta.val()),'ell':ELL}
 return m,s,bg

def serialize(j):return {k:v.row() if isinstance(v,I) else str(v) for k,v in j.items()}

def report(cells):
 data,errors,pins=inputs();rows=[]
 for label,x,r in cells:
  m,s,b=actual_jet(x,r);mc,sc,bc=actual_jet(x,r,mode='comparator')
  rows.append({'label':label,'scaled_x':I.cv(x).row(),'scaled_r':I.cv(r).row(),
               'actual_m':serialize(m),'actual_s':serialize(s),'background':serialize(b),
               'comparison_m':serialize(mc),'comparison_s':serialize(sc)})
 return {'status':'PASS_ACTUAL_FIRST_COMMON_GAUGE_ORDINARY_2JET_CELL_ENCLOSURES',
         'domain':'0<=x<=3,r>=0,x+r<=8; t=ell*x,z=ell*r,dy=-rho dz',
         'ell':str(ELL),'arithmetic':{'outward_grid':str(core.GRID),'exact_rational_decisions':True},
         'input_sha256':{k:v['sha256'] for k,v in pins.items()},
         'cells':rows,'scope':{'actual_first_solution_jets':True,'numerical_picard_refinement':False,
         'predecessor_smooth_PDE_existence_is_conditional_premise':True,'second_response_evaluated':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 need(not a.output.exists(),'refuse overwrite')
 cells=[('evolved_packet',I(F(1,4),F(1,4)+F(1,4096)),I(F(5,4),F(5,4)+F(1,4096))),
        ('evolved_reflected_packet',I(F(3,2),F(3,2)+F(1,4096)),I(F(1,4),F(1,4)+F(1,4096)))]
 r=report(cells);a.output.write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
if __name__=='__main__':main()
