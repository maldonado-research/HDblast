#!/usr/bin/env python3
"""Continuous whole-cell rational evaluation of the actual second initial data."""
import argparse,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent;ELL=F(1,10**6);GRID=1<<112

def need(ok,msg):
 if not ok:raise ValueError(msg)
class I:
 def __init__(self,lo=0,hi=None):
  lo=F(lo);hi=lo if hi is None else F(hi);need(lo<=hi,'interval order')
  self.lo=F((lo*GRID).__floor__(),GRID);self.hi=F(-(-hi*GRID).__floor__(),GRID)
 @staticmethod
 def cv(x):return x if isinstance(x,I) else I(x)
 def __add__(self,b):b=I.cv(b);return I(self.lo+b.lo,self.hi+b.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,b):return self+-I.cv(b)
 def __rsub__(self,b):return I.cv(b)+-self
 def __mul__(self,b):
  b=I.cv(b);v=[a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi)];return I(min(v),max(v))
 __rmul__=__mul__
 def inv(self):need(self.lo>0 or self.hi<0,'zero denominator');return I(1/self.hi,1/self.lo)
 def __truediv__(self,b):return self*I.cv(b).inv()
 def __rtruediv__(self,b):return I.cv(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  if n==0:return I(1)
  if n%2==0:
   a=self**(n//2);return I(0 if a.lo<=0<=a.hi else min(a.lo*a.lo,a.hi*a.hi),max(a.lo*a.lo,a.hi*a.hi))
  return self*self**(n-1)
 def mag(self):return max(abs(self.lo),abs(self.hi))
 def row(self):return [str(self.lo),str(self.hi)]
class J:
 def __init__(self,c):self.c=list(map(I.cv,c))
 @property
 def n(self):return len(self.c)-1
 def cv(self,b):return b if isinstance(b,J) else J([b]+[0]*self.n)
 def __add__(self,b):b=self.cv(b);return J([self.c[k]+b.c[k] for k in range(min(self.n,b.n)+1)])
 __radd__=__add__
 def __neg__(self):return J([-x for x in self.c])
 def __sub__(self,b):return self+-self.cv(b)
 def __rsub__(self,b):return self.cv(b)+-self
 def __mul__(self,b):b=self.cv(b);return J([sum((self.c[k]*b.c[j-k] for k in range(j+1)),I(0)) for j in range(min(self.n,b.n)+1)])
 __rmul__=__mul__
 def inv(self):
  out=[self.c[0].inv()]
  for j in range(1,self.n+1):out.append(-out[0]*sum((self.c[k]*out[j-k] for k in range(1,j+1)),I(0)))
  return J(out)
 def __truediv__(self,b):return self*self.cv(b).inv()
 def __rtruediv__(self,b):return self.cv(b)*self.inv()
 def __pow__(self,n):
  if n<0:return self.inv()**(-n)
  out=self.cv(1)
  for _ in range(n):out=out*self
  return out
 def dz(self):return J([(k+1)*self.c[k+1] for k in range(self.n)])
 def val(self):return self.c[0]
U=list(map(F,['-2/27','0','14/9','50/27','-1/6','-4/9','-2/27']))
UP=[U]
for _ in range(4):UP.append([k*c for k,c in enumerate(UP[-1]) if k])
def poly(cs,x):
 out=0
 for c in reversed(cs):out=out*x+c
 return out

def expminus(q):
 """Enclose exp(-q), rational Taylor at reduced argument, then squaring."""
 need(q>=0,'nonnegative exponential argument')
 if not q:return I(1)
 halvings=0;t=q
 while t>1:t/=2;halvings+=1
 degree=28;upper=sum(((-t)**k/math.factorial(k) for k in range(degree+1)),F(0));lower=upper+(-t)**(degree+1)/math.factorial(degree+1)
 need(0<lower<=upper,'exponential alternating enclosure');out=I(lower,upper)
 for _ in range(halvings):out=out*out
 return out

def expbox(q):return I(expminus(q.hi).lo,expminus(q.lo).hi)
def centered_poly(cs,s):
 mid=(s.lo+s.hi)/2;v=s-mid;coeff=[sum((cs[k]*math.comb(k,j)*mid**(k-j) for k in range(j,len(cs))),F(0)) for j in range(len(cs))]
 return poly(coeff,v)

def load():
 pins=json.loads((ROOT/'INPUT_PINS.json').read_text());data={}
 for k,pin in pins.items():
  raw=(ROOT/pin['path']).read_bytes();need(len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256'],'input authentication: '+k);data[k]=json.loads(raw)
 return data,pins

def background(r,data,injection):
 co=data['coefficients'];z=ELL*r;shell=data['shell'];vals={}
 for name in ('rho','eta','p','a'):
  entry=co['first_exit'][name];base=I(*shell['boxes'][name+'_b']['exact']);rhs=I(*entry['outer_ODE_image']['exact']);vals[name]=base+z*rhs
 need(vals['rho'].lo>0 and vals['p'].hi<0,'actual branch signs')
 kg=I(*co['alpha_global']['log_derivative']['exact']);need(kg.hi<=0 and 1-kg.mag()*z.hi>0,'negative alpha logarithmic derivative and positive enclosure');vals['alpha']=expbox(-kg*z)
 state={k:J([v]+[0]*4) for k,v in vals.items()}
 for n in range(4):
  rho,e,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')];h=poly(UP[1],e)/p
  rhs={'rho':-rho*rho*a,'eta':-rho*p,'p':-rho*(poly(UP[1],e)-4*a*p),'a':1/rho+rho*p*p/3,'alpha':rho*(F(5,2)*a-h)*alpha}
  for k in state:state[k].c[n+1]=rhs[k].c[n]/(n+1)
 rho,e,p,a,alpha=[state[k] for k in ('rho','eta','p','a','alpha')]
 h=poly(UP[1],e)/p;k=rho*(F(5,2)*a-h);C=alpha/rho**2;B=3*alpha/(rho**3*p)
 E0=expbox(z);E=J([(-1)**j*E0/math.factorial(j) for j in range(5)])
 return rho,e,p,a,h,C,B,k,poly(UP[2],e),E

def bump(r,data,tail=False):
 polys=[list(map(F,row)) for row in data['bump']['polynomial_coefficients_ascending'][:5]]
 if tail:
  edge=F(1,64);smin=1/(edge*(1-edge));em=expminus(smin-4);out=[]
  for n in range(5):
   cap=smin**(2*n)*em.hi
   need(smin>=2*n,'flat-tail exponential-polynomial monotonicity')
   out.append(I(0,cap) if n==0 else I(-cap,cap))
 else:
  s=1/(F(1,4)-(r-F(3,2))**2);v=3-2*r;ex=expbox(s-4)
  out=[centered_poly(polys[n],s)*ex*(v if n%2 else 1) for n in range(5)]
 return J([out[n]/(ELL**n*math.factorial(n)) for n in range(5)])

def calculate(r,data,tail=False,injection='none'):
 rho,e,p,a,h,C,B,k,U2,E=background(r,data,injection);q=bump(r,data,tail)
 ir=1/rho;ip=1/p;dy=lambda j:-ir*j.dz();u=C*q;v=B*(q.dz()+k*q);vy=dy(v);ut=E*u;d=3*ir*ip;vt=E*(v-d*u)
 A=-2*vy**2-8*p*u*vy+8*p*(h+2*a)*u*v+(16*p*p-192*ir**2)*u*u-(2*U2+F(8,3)*p*p)*v*v-ir**2*(12*ut*ut+2*vt*vt)
 D=ir**2*(36*u*u-2*v*v)
 sectors=[A/6+F(4,3)*D,A/2+2*D,A-4*D]
 fs=[-rho**2*sectors[j]/den for j,den in enumerate((-6,10,42))];fy=[dy(x) for x in fs]
 zeta=[-F(3,2)*ip*x for x in fy];zetay=[dy(x) for x in zeta]
 K=E*ir/3*(18*u*u-v*v+d*u*v)
 Pi=E*ir*ip*(d*u*vy+(d*(h-ir)+10*p)*u*v+(30*a-6*ir)*u*u+(2*ir-a)*v*v)
 weights=[F(1),F(1,2),F(7,12)]
 abs0=lambda j:j.val().mag();total=lambda seq:sum((w*abs0(j) for w,j in zip(weights,seq)),F(0));sphere=lambda seq:abs0(seq[1])+F(5,2)*abs0(seq[2])
 gv=total(fs);gyr=total(fy);pv=total(zeta);pyr=total(zetay);rmin=rho.val().lo;amax=a.val().mag();ustar=abs0(u);utstar=abs0(ut);uystar=abs0(dy(u));Kstar=abs0(K);Pstar=abs0(Pi)
 normal=[8*Kstar+44*ustar*utstar/rmin,26*ustar*uystar+F(3,2)*gyr+3*amax*gv,44*ustar**2/rmin+sphere(fs)/(2*rmin),2*Kstar+8*ustar*utstar/rmin,2*Kstar+4*ustar*utstar/rmin]
 bounds={'metric_value':max(gv,2*ustar**2),'metric_spatial_covariant_derivative':max(gyr,sphere(fs)/rmin,amax*gv,4*ustar*uystar,8*ustar**2/rmin),'metric_normal_covariant_derivative':max(normal),'scalar_value':pv,'scalar_spatial_derivative':max(pyr,sphere(zeta)/rmin),'scalar_normal_derivative':Pstar+2*ustar*abs0(vt)/rmin}
 if injection=='omit_scalar_initial':bounds['scalar_value']=F(0)
 need(bounds['scalar_value']>0,'nonzero second scalar initial data retained')
 # Full frame-component partial time velocity: omega0 mixes 0/y indices.
 # k_0y=0, k_yy=0, k_00=-2u²: only normal 0y changes by -2a u².
 timepartial=max(normal[0],normal[1]+2*amax*ustar**2,normal[2],normal[3],normal[4])
 bounds.update({'metric_scaled_radial_partial':ELL*rho.val().hi*bounds['metric_spatial_covariant_derivative'],'scalar_scaled_radial_partial':ELL*rho.val().hi*bounds['scalar_spatial_derivative'],'metric_scaled_time_partial':ELL*rho.val().hi*timepartial,'scalar_scaled_time_partial':ELL*rho.val().hi*bounds['scalar_normal_derivative']})
 return bounds,{'r':r.row(),'scalar_mean':zeta[0].val().row(),'scalar_degree2_amplitude':zeta[1].val().row(),'scalar_degree4_amplitude':zeta[2].val().row(),'scalar_mean_scaled_radial':(ELL*zeta[0].dz()).val().row(),'metric_angular_mean':fs[0].val().row(),'lapse_square_amplitude':(-2*u*u).val().row(),'pure_trace_K':K.val().row(),'scalar_Pi':Pi.val().row()}

def run(cells,injection):
 data,pins=load();need(cells>=64 and cells%64==0,'partition divisible by64');caps={};rows=[]
 # Complete cover; compact data outside [1,2] are exactly zero.
 pieces=[(F(1),F(65,64),True)]+[(F(k,cells),F(k+1,cells),False) for k in range(cells+cells//64,2*cells-cells//64)]+[(F(127,64),F(2),True)]
 for lo,hi,tail in pieces:
  bounds,row=calculate(I(lo,hi),data,tail,injection);row['analytic_flat_endpoint_tail']=tail;rows.append(row)
  for k,v in bounds.items():caps[k]=max(caps.get(k,F(0)),v)
 need(pieces[0][0]==1 and pieces[-1][1]==2 and all(pieces[j][1]==pieces[j+1][0] for j in range(len(pieces)-1)),'continuous whole-annulus partition cover')
 return {'status':'PASS_CONTINUOUS_ACTUAL_SECOND_INITIAL_PROFILE_ENCLOSURES','domain':{'ell':str(ELL),'scaled_r':['1','2'],'time':'0','whole_sphere':True,'smooth_zero_extension_outside_packet':True},'conventions':{'coefficients':'Second parameter derivatives; epsilon^2/2','metric_frame':'e0=rho^-1 partial_t, ey=partial_y, ea=rho^-1 Ea at t=0; parallel along y','radial':'partial_r=-ell*rho partial_y','scaled_time':'partial_x=ell*rho frame-component partial_0; frame connection included','scalar_mean':'Probability-normalized sphere mean, W0=1; zeta0'},'input_sha256':{k:v['sha256'] for k,v in pins.items()},'arithmetic':{'only_exact_rational_decisions':True,'outward_grid_denominator':str(GRID),'exponential':'alternating degree28 Taylor after halving, squared with outward intervals','bump':'actual derivative polynomials through4, not independent global derivative boxes','background':'local actual shell box plus ODE integral image; ODE differentiated recurrences through4; alpha=exp(integral k) with actual k box. Shared functional expressions retained. Shell uncertainty not reduced.'},'partition':{'interior_step':str(F(1,cells)),'whole_cells':len(rows),'endpoint_tail_width':'1/64','endpoint_tail_proof':'s>=4096/63>2k; |Q^(k)|<=s^(2k)exp(4-s), k<=4; monotone tail bound, flat endpoint zero'},'component_caps':{k:{'exact_upper':str(v),'display_only':float(v)} for k,v in caps.items()},'profile_cells':rows,'scope':{'actual_initial_profiles_enclosed':True,'continuous_cell_and_tail_coverage':True,'initial_lapse_square_retained':True,'mixed_shift_velocity_retained':True,'scalar_initial_and_velocity_retained':True,'actual_background_uncertainty_resolved':False,'initial_interpolant_approximation_error_certified':False,'second_response_evaluated':False}}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--cells',type=int,default=512);p.add_argument('--inject-error',choices=['none','omit_scalar_initial'],default='none');a=p.parse_args();need(not a.output.exists(),'refuse overwrite');r=run(a.cells,a.inject_error);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'status':r['status'],'cells':r['partition']['whole_cells'],'caps':{k:v['display_only'] for k,v in r['component_caps'].items()}},indent=2))
if __name__=='__main__':main()
