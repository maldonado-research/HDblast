#!/usr/bin/env python3
"""Exact rational analytic comparison for the retained existent regular branch."""
import json, hashlib, math, shutil
from pathlib import Path
from fractions import Fraction as F
from authenticate import authenticate
HERE=Path(__file__).resolve().parent

def require(b,msg):
 if not b: raise ValueError(msg)
def clean(p):return {k:v for k,v in p.items() if v}
def add(*ps):
 r={}
 for p in ps:
  for k,v in p.items():r[k]=r.get(k,F(0))+v
 return clean(r)
def scale(p,c):return clean({k:v*c for k,v in p.items()})
def mul(p,q):
 r={}
 for (i,j),v in p.items():
  for (a,b),w in q.items():r[i+a,j+b]=r.get((i+a,j+b),F(0))+v*w
 return clean(r)
def power(p,n):
 r={(0,0):F(1)}
 for _ in range(n):r=mul(r,p)
 return r
def dx(p):return clean({(i-1,j):i*v for (i,j),v in p.items() if i})
def dh(p):return clean({(i,j-1):j*v for (i,j),v in p.items() if j})
def val(p,x,h):return sum((v*x**i*h**j for (i,j),v in p.items()),F(0))
def major(p,x,h):return sum((abs(v)*x**i*h**j for (i,j),v in p.items()),F(0))
def coarse(p,x,h):
 c=p.get((0,0),F(0));r=major(add(p,{(0,0):-c}),x,h)
 return c-r,c+r
class I:
 def __init__(self,a,b=None):self.a=F(a);self.b=F(a if b is None else b);require(self.a<=self.b,'interval order')
 def __add__(self,q):
  if not isinstance(q,I):q=I(q)
  return I(self.a+q.a,self.b+q.b)
 __radd__=__add__
 def __neg__(self):return I(-self.b,-self.a)
 def __sub__(self,q):return self+-q if isinstance(q,I) else self+(-q)
 def __rsub__(self,q):return -self+q
 def __mul__(self,q):
  if not isinstance(q,I):q=I(q)
  v=[self.a*q.a,self.a*q.b,self.b*q.a,self.b*q.b];return I(min(v),max(v))
 __rmul__=__mul__
 def __pow__(self,n):
  if n==0:return I(1)
  if n%2==0:
   return I(0 if self.a<=0<=self.b else min(self.a**n,self.b**n),max(self.a**n,self.b**n))
  return I(self.a**n,self.b**n)
 def out(self):return [str(self.a),str(self.b)]
def iv(p,x,h):return sum((v*(x**i)*(h**j) for (i,j),v in p.items()),I(0))
def sqrt_out(x,digits=32):
 n=10**digits;s=math.isqrt((x.numerator*n*n)//x.denominator);lo=F(s,n);hi=F(s+1,n)
 require(lo*lo<=x<=hi*hi,'sqrt rational bound');return lo,hi

def run():
 pre=authenticate()
 graph_source=HERE.parent/'graph-algebra/graph_degree_Z4_P5.json'
 graph_local=HERE/'GRAPH_POLYNOMIALS.json'
 if graph_source.exists() and not graph_local.exists():shutil.copyfile(graph_source,graph_local)
 raw=graph_local.read_bytes();data=json.loads(raw)
 ps={key:{(int(r['i']),int(r['j'])):F(int(r['num']),int(r['den'])) for r in data[key]} for key in ['U','M','Z','A','B','R','D']}
 U,M,Z,A,B,R,D=(ps[k] for k in ['U','M','Z','A','B','R','D'])
 x={(1,0):F(1)};h={(0,1):F(1)};one={(0,0):F(1)}
 W={(0,0):F(1,3),(2,0):F(1),(3,0):F(1,3)}
 require(U==add(scale(power(dx(W),2),F(1,2)),scale(power(W,2),F(-2,3))),'exact U model')
 require(mul(x,M)==dx(U),'M identity')
 require(B==add(scale(mul(h,dh(Z)),2),scale(Z,-4)),'B identity')
 require(R==add(M,scale(power(Z,2),-1),scale(mul(mul(x,dx(Z)),Z),-1),mul(A,B)),'R identity')
 require(D==add(h,scale(mul(power(x,2),power(Z,2)),F(1,12)),scale(U,F(-1,6)),scale(power(A,2),-1)),'D identity')
 require(all(i+j>=5 for i,j in R),'R degree cancellation')
 require(all(i+j>=5 for i,j in D),'D degree cancellation')
 eps=F(9,100000);h0=F(1,100);amin=F(1,9)-F(1,1000000);amax=F(3,20);kap=F(37,20)
 zlo,zhi=coarse(Z,eps,h0);alo,ahi=coarse(A,eps,h0);blo,bhi=coarse(B,eps,h0)
 xzx=major(mul(x,dx(Z)),eps,h0)
 require(zlo>0 and zhi<2,'0<Z<2');require(alo>0,'A>0');require(bhi<0,'B<0')
 require(4*amin+zlo-xzx>=kap,'secant damping')
 require(kap-12*amax>0,'memory damping')
 Mlo,Mhi=coarse(M,eps,0)
 require(0<Mlo and Mhi<4,'0<M<4 gives 0<z<2')
 Udiff_poly=add(U,{(0,0):F(2,27)})
 Udiff_over_x2={(i-2,j):v for (i,j),v in Udiff_poly.items()}
 require(coarse(Udiff_over_x2,eps,0)[0]>0,'U>=U0')
 Udiff=major(Udiff_poly,eps,0)
 require(-F(2,27)+Udiff<0,'U<0')
 require(F(1,81)-Udiff/6>amin**2,'constraint a lower bound')
 require(h0+F(1,81)+eps**2/3<amax**2,'constraint a upper bound')
 # Exact coefficientwise majorant for |R+B*D/(a+A)|; all x powers use eps.
 Babs=max(abs(blo),abs(bhi));fac=Babs/(amin+alo)
 C=[F(0) for _ in range(6)]
 for poly,factor in [(R,F(1)),(D,fac)]:
  for (i,j),v in poly.items():C[min(j,5)]+=factor*abs(v)*eps**i*h0**max(j-5,0)
 coeff=[C[j]/(kap-2*j*amax) for j in range(6)]
 def E(H):return 2*(H/h0)**6+sum((coeff[j]*H**j for j in range(6)),F(0))
 # Old theorem selection, used only as branch existence and initial enclosure.
 c=F(96714321425692,161839258930731);delta=F(1,1000)
 beta0=-F(9,64)*c*delta;oldrad=F(799,250)*delta**2
 old=I(beta0-oldrad,beta0+oldrad)
 require(old.b<0 and max(abs(old.a),abs(old.b))<eps,'retained negative beta region')
 f=add(scale(one,1+c),scale(x,c))
 H=add(scale(add(scale(mul(W,f),F(1,9)),scale(dx(W),-c/F(12))),delta),scale(add(scale(power(f,2),F(1,36)),scale(one,-c*c/F(48))),delta**2))
 p=add(scale(dx(W),-1),scale(one,-delta*c/2))
 a=add(scale(W,F(1,3)),scale(f,delta/6))
 Hbox=iv(H,old,I(0));require(0<Hbox.a<Hbox.b<h0,'shell inside comparison entry H')
 Hp=dx(H)
 def g(t):return t*val(Z,t,val(H,t,0))+val(dx(W),t,0)+delta*c/2
 def gp(t):
  ht=val(H,t,0)
  return val(Z,t,ht)+t*(val(dx(Z),t,ht)+val(dh(Z),t,ht)*val(Hp,t,0))+val(dx(dx(W)),t,0)
 seed=beta0;seed_log=[];scale24=10**24
 for n in range(4):
  proposal=seed-g(seed)/gp(seed)
  seed=F((proposal.numerator*scale24)//proposal.denominator,scale24)
  seed_log.append({'iteration':n+1,'beta_seed':str(seed),'polynomial_residual':str(g(seed))})
 gpbox=iv(Z,old,Hbox)+old*(iv(dx(Z),old,Hbox)+iv(dh(Z),old,Hbox)*iv(Hp,old,I(0)))+iv(dx(dx(W)),old,I(0))
 require(gpbox.a>0,'shell polynomial derivative positive')
 error=eps*E(Hbox.b)
 rad=(abs(g(seed))+error)/gpbox.a
 require(old.a<=seed-rad and seed+rad<=old.b,'contracted interval inside retained root region')
 # Widen outward to terminating rational endpoints, independently check.
 scale20=10**20
 beta_lo=F(((seed-rad).numerator*scale20)//(seed-rad).denominator,scale20)
 beta_hi=F(-((-(seed+rad)).numerator*scale20//(-(seed+rad)).denominator),scale20)
 require(beta_lo<=seed-rad and beta_hi>=seed+rad,'outward rational beta endpoints')
 beta=I(beta_lo,beta_hi);newH=iv(H,beta,I(0));newp=iv(p,beta,I(0));newa=iv(a,beta,I(0))
 rholo=sqrt_out(1/newH.b)[0];rhohi=sqrt_out(1/newH.a)[1]
 post=authenticate();require(pre==post,'old closure changed during arithmetic')
 def pack(v):return {'exact':str(v),'display':float(v)}
 result={'status':'PASS_RETAINED_BRANCH_ANALYTIC_COMPARISON_ENCLOSURE','route':'degree4 normalized-slope graph with exact polynomial defect and H-power comparison; retained existence used','scope':'Necessary enclosure of the same existing unique regular branch in its retained tube. No independent finite endpoint-map solve, no total U evaluation.','graph_sha256':hashlib.sha256(raw).hexdigest(),'dependency_manifest_sha256':pre['manifest_sha256'],'finite_cap':{'graph_degree':4,'comparison_H_domains':1,'radial_panels':0,'rational_seed_steps':4,'validated_contractions':1,'interval_newton_steps':0},'domain':{'eta_abs_max':str(eps),'H_max':str(h0),'actual_z':[0,2],'a_min':str(amin),'a_max':str(amax)},'exact_checks':{'U_model':True,'M_identity':True,'B_identity':True,'R_identity':True,'D_identity':True,'R_D_vanish_through_degree4':True,'M_bounds':[pack(Mlo),pack(Mhi)],'U_minus_U0_nonnegative':True,'U_negative':True,'Z_bounds':[pack(zlo),pack(zhi)],'A_bounds':[pack(alo),pack(ahi)],'B_bounds':[pack(blo),pack(bhi)],'eta_Zeta_abs':pack(xzx),'secant_damping_lower':pack(4*amin+zlo-xzx),'selected_damping':str(kap),'memory_damping_margin':str(kap-12*amax)},'residual_majorant_Cj':[str(v) for v in C],'barrier_Ej':[str(v) for v in coeff],'old_beta':old.out(),'old_H2':Hbox.out(),'endpoint_error_z':pack(E(Hbox.b)),'endpoint_error_scalar_residual':pack(error),'scalar_polynomial_derivative':gpbox.out(),'seed_steps':seed_log,'pre_rounding_radius':pack(rad),'beta':beta.out(),'beta_midpoint':str((beta_lo+beta_hi)/2),'beta_radius':pack((beta_hi-beta_lo)/2),'provisional_beta_precision_target_passed':(beta_hi-beta_lo)/2<=F(1,10**11),'shell_correlated_definition':{'parameter':'beta in certified beta interval','H2':'delta*(W*f/9-W_prime*c/12)+delta^2*(f^2/36-c^2/48)','p':'-W_prime-delta*c/2','a':'W/3+delta*f/6','rho':'1/sqrt(H2)','W':'1/3+beta^2+beta^3/3','f':'1+c+c*beta','c':str(c),'delta':str(delta)},'shell_outer_projections':{'H2':newH.out(),'p':newp.out(),'a':newa.out(),'rho':[str(rholo),str(rhohi)]},'post_authentication_passed':True}
 (HERE/'BRANCH_COMPARISON_REPORT.json').write_text(json.dumps(result,indent=2)+'\n')
 print(result['status']);print('beta',float(beta.a),float(beta.b),'radius',float((beta.b-beta.a)/2));print('zerr',float(E(Hbox.b)),'damping',float(4*amin+zlo-xzx),'gp_lower',float(gpbox.a))
 return result
if __name__=='__main__':run()
