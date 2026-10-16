#!/usr/bin/env python3
"""Exact algebra and provenance controls for the actual first-jet adapter."""
import argparse,importlib.util,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('actual_first_jets',ROOT/'actual_first_jets.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
I=a.I;need=a.need

def contains(v,x):return v.lo<=x<=v.hi
def overlap(v,w):return max(v.lo,w.lo)<=min(v.hi,w.hi)
def falling(n,k):return F(math.factorial(n),math.factorial(n-k)) if n>=k else F(0)

def run(injection='none'):
 data,errors,pins=a.inputs();checks=[]
 spec=importlib.util.spec_from_file_location('retained_error_replay',ROOT/'check_inherited_errors.py')
 old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 need(old.run()==errors,'inherited full actual derivative certificate exact replay');checks.append('authenticated inherited complete C3 certificate exact replay')
 # An independently differentiated product, with Pe(s)=s^5 and signed
 # Qo(s)=sign(s)s^4, checks every comparator mixed derivative on both
 # characteristic branches.  The direct product expression is not the
 # branch-simplified expression used in actual_first_jets.py.
 witnesses=0
 for x,r in ((F(1,2),F(5,4)),(F(3,2),F(1,4))):
  sp,sm=r+x,r-x;sg=1 if sm>0 else -1
  for i,j in a.KEYS:
   n=i+j
   direct=falling(5,n)*(sp**(5-n)+(-1)**i*sm**(5-n))
   direct-=r*falling(4,n)*(sp**(4-n)+sg*(-1)**i*sm**(4-n))
   if j:direct-=j*falling(4,n-1)*(sp**(5-n)+sg*(-1)**i*sm**(5-n))
   plus=x*falling(4,n)*sp**(4-n)
   if i:plus+=i*falling(4,n-1)*sp**(5-n)
   negative=sm<0;lin=x-2*r if negative else -x;lx=1 if negative else -1;lr=-2 if negative else 0
   minus=lin*(-1)**i*sg*falling(4,n)*sm**(4-n)
   if i:minus+=i*lx*(-1)**(i-1)*sg*falling(4,n-1)*sm**(5-n)
   if j and injection!='drop_gauge_product':minus+=j*lr*(-1)**i*sg*falling(4,n-1)*sm**(5-n)
   need(direct==plus+minus,'characteristic gauge product derivative');witnesses+=1
 checks.append('20 independent characteristic product derivatives, including reflected branch')
 # Test the normalized two-variable Taylor arithmetic on a nontrivial
 # rational function, with products, inverses, and mixed derivatives.
 x=a.T({(0,0):F(2,3),(1,0):1});r=a.T({(0,0):F(5,4),(0,1):1})
 fun=(x*x+3*x*r+r*r*r)/(2+r)
 residual=fun*(2+r)-(x*x+3*x*r+r*r*r)
 need(all(contains(v,0) for v in residual.c.values()),'bivariate rational Taylor identity')
 checks.append('bivariate normalized Taylor product/inverse identity through3')
 c,h,v=a.temporal(I(F(1,4),F(1,4)+F(1,4096)))
 identity=v*v*c*c*c
 need(contains(identity.val(),1) and all(contains(value,0) for key,value in identity.c.items() if key!=(0,0)),'temporal cosh power differentiated identity')
 checks.append('temporal cosh^-3/2 differentiated identity through3')
 q=a.qodd(I(F(3,2)),3,data)
 need(all(contains(v,x) for v,x in zip(q,(1,0,-32,0))),'exact midpoint packet derivatives')
 for s in (I(F(5,4),F(5,4)+F(1,4096)),I(F(65,64),F(66,64))):
  qp=a.qodd(s,3,data);qm=a.qodd(-s,3,data)
  need(all(overlap(qm[n],(-1)**(n+1)*qp[n]) for n in range(4)),'signed packet jet parity')
 for s in (I(-3,-2),I(-1,1),I(2,8)):
  need(all(v.lo==v.hi==0 for v in a.qodd(s,3,data)),'exact compact zero collar')
 checks.append('packet midpoint, parity, exact zero collar and flat-tail branch')
 # Derive physical y^2 by applying -1/(ell rho) d_r twice.  rho_r is
 # fixed by the actual background ODE; this is independent of Eq.(8).
 rho0=F(13);a0=F(2,7);rho=a.T({(0,0):rho0,(0,1):-a.ELL*a0*rho0**2})
 radial=a.T({(0,0):F(5,4),(0,1):1});fun=radial**3
 dy=lambda j:-(1/(a.ELL*rho))*j.dr()
 direct=dy(dy(fun)).val();converted=a.physical_derivatives(fun,rho,a.T.cv(a0))['yy']
 if injection=='drop_radial_conversion':converted=(1/(a.ELL*rho0))**2*fun.derivative(0,2)
 need(overlap(direct,converted),'variable-rho physical second derivative conversion')
 checks.append('variable-rho ordinary second radial derivative')
 # Reconstruct one independently manufactured raw vector and longitudinal
 # pair.  Pulling back phi_0 gives +f*p directly from xi^y d_y phi_0;
 # the yy Lie derivative gives 2 f_y=-2 f_r/(ell*rho).
 W=a.T({(0,0):2,(0,1):3});f=a.T({(0,0):31,(0,1):43});b=a.T.cv(37);d=a.T.cv(41)
 values=list(map(a.T.cv,(13,-11,17,3,5,7,19,23,29)))
 manufactured,scalar=a.frame_amplitudes(W,[f,b,d],*values,injection=injection)
 expected_scalar=5*29*(3/a.ELL+7*2)-11*31
 need(contains(scalar.val(),expected_scalar),'registered scalar pullback shift is +pf')
 expected_yy=-4*3*29*2-2*43/(a.ELL*13)
 need(contains(manufactured['A'].val(),expected_yy),'registered raw gauge radial frame factor includes rho')
 checks.append('independent manufactured scalar pullback and metric radial-frame factor')
 shell_m,shell_s,shell_b=a.actual_jet(I(F(5,4),F(5,4)+F(1,4096)),I(0))
 need(all(shell_m[key].lo==shell_m[key].hi==0 for key in ('B','Bt','Btt','C','Ct','Ctt')),'physical mixed shell rows and tangential derivatives')
 checks.append('exact physical mixed shell rows and their tangential jets')
 ext=a.exterior_caps()
 need(ext[0][3]<F(1,2000),'exterior causal-lift restriction improves actual third derivative')
 checks.append('backward-closed exterior derivative certificate exact baseline replay')
 # Distinct true and comparison solutions remain explicitly distinguished.
 actual,scalar,bg=a.actual_jet(I(F(1,4),F(1,4)+F(1,4096)),I(F(5,4),F(5,4)+F(1,4096)),injection=injection)
 comparison,_,_=a.actual_jet(I(F(1,4),F(1,4)+F(1,4096)),I(F(5,4),F(5,4)+F(1,4096)),mode='comparator')
 need(any(actual[k].hi-actual[k].lo>comparison[k].hi-comparison[k].lo for k in actual),'actual derivative uncertainty retained')
 checks.append('actual-minus-comparison derivative uncertainty retained')
 return {'status':'PASS_ACTUAL_FIRST_JET_ADAPTER_CONTROLS','checks':checks,'manufactured_characteristic_derivative_witnesses':witnesses,
         'input_sha256':{k:v['sha256'] for k,v in pins.items()},
         'scope':'Internal computational check; retained smooth actual master/gauge/background PDE constructions are conditional premises. No finite-difference derivative or C0-error differentiation is used.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--inject-error',default='none',choices=['none','drop_gauge_product','drop_radial_conversion','wrong_scalar_shift','drop_radial_rho']);args=p.parse_args()
 need(not args.output.exists(),'refuse overwrite');report=run(args.inject_error);args.output.write_text(json.dumps(report,indent=2)+'\n');print(report['status'])
if __name__=='__main__':main()
