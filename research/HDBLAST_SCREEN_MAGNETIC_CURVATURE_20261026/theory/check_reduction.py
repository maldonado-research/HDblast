#!/usr/bin/env python3
"""Exact, all-sphere angular reductions for scalar screen magnetic H''."""
import argparse
import sys
sys.dont_write_bytecode=True
import hashlib
import types
import json
from fractions import Fraction as Q
from itertools import product
from pathlib import Path

HERE=Path(__file__).resolve().parent
DEFAULT=HERE.parent/'angular'/'inputs'/'angular_exact.py'
ANGULAR_SHA256='2c0bafec777581c79c3fd6bd9879d0d3cb8114bb899722bd33882086751ded60'


def validate_output(path):
 """Reject overwrite and source/capsule writes before scientific work."""
 path=Path(path)
 if path.exists() or path.is_symlink():raise ValueError('FRESH_NONSYMLINK_OUTPUT_REQUIRED')
 target=path.resolve()
 protected=[HERE]+[a for a in HERE.parents if a.parent.name=='research']
 if any(target.is_relative_to(root) for root in protected):
  raise ValueError('OUTPUT_INSIDE_SOURCE_OR_RESEARCH_CAPSULE')
 if not target.parent.is_dir():raise ValueError('OUTPUT_PARENT_DIRECTORY_REQUIRED')
 return target


def authenticated_angular_bytes(path):
 """Authenticate the exact bytes subsequently compiled; no path reread."""
 raw=Path(path).read_bytes()
 if hashlib.sha256(raw).hexdigest()!=ANGULAR_SHA256:
  raise ValueError('ANGULAR_SOURCE_SHA256_MISMATCH')
 return raw


def load_angular(path):
 raw=authenticated_angular_bytes(path)
 r=types.ModuleType('frozen_angular')
 r.__file__=str(Path(path).resolve())
 exec(compile(raw,r.__file__,'exec'),r.__dict__)
 return r


def run(path):
 r=load_angular(path)
 I=range(4);Y=r.Y;P=r.P;v=r.grad(Y);h=r.hess(Y);dh=r.covd(h)
 eps={(a,c,d):r.add(*(r.scale(r.EPS.get((j,a,c,d),0),r.X[j]) for j in I)) for a,c,d in product(I,repeat=3)}
 def sym(t):return {(a,b):r.scale(Q(1,2),r.add(t[a,b],t[b,a])) for a,b in product(I,repeat=2)}
 def dual(j):return {(a,b):r.add(*(r.mul(eps[a,c,d],j[c,b,d]) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)}
 def curl(t):return sym(dual(r.covd(t)))
 def test(a,b,name):
  if not r.equal(a,b):raise ValueError(name)
  checks.append(name)
 checks=[]
 zero=r.zero_tensor(2)
 third={(c,a,b):r.scale(-1,r.add(r.mul(P[c,a],v[b,]),r.mul(P[c,b],v[a,]),r.scale(2,r.mul(P[a,b],v[c,])))) for c,a,b in product(I,repeat=3)}
 test(dh,third,'third_derivative_quadrupole_all_point')
 A=sym({(a,b):r.add(*(r.mul(eps[a,c,d],r.mul(v[c,],h[b,d])) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)})
 coeff=[r.inner(A,t)/r.inner(t,t) for t in r.TT]
 test(A,r.tadd(*(r.tscale(c,t) for c,t in zip(coeff,r.TT))),'axial_tensor_complete_TT_reconstruction')
 test(curl(r.tmul(Y,P)),zero,'scalar_trace_curl_zero')
 test(curl(h),zero,'scalar_hessian_curl_zero')
 test(curl(r.tmul(Y,h)),A,'Y_hessian_curl')
 test(curl(r.outer(v,v)),r.tscale(-1,A),'gradient_square_curl')
 # Screen connection variations for p=Y*kappa and p=HessY.
 gamL={(e,c,b):r.scale(Q(1,2),r.add(r.mul(P[e,b],v[c,]),r.mul(P[e,c],v[b,]),r.scale(-1,r.mul(P[c,b],v[e,])))) for e,c,b in product(I,repeat=3)}
 gamM={(e,c,b):r.scale(-1,r.add(r.mul(P[e,c],v[b,]),r.mul(P[e,b],v[c,]))) for e,c,b in product(I,repeat=3)}
 for name,p,g in [('L',r.tmul(Y,P),gamL),('M',h,gamM)]:
  dp=r.covd(p)
  direct={(e,c,b):r.scale(Q(1,2),r.add(dp[c,b,e],dp[b,c,e],r.scale(-1,dp[e,c,b]))) for e,c,b in product(I,repeat=3)}
  test(direct,g,'connection_'+name)
  for tag,S,target in [('trace',r.tmul(Y,P),zero),('hessian',h,r.tscale(-Q(1,2) if name=='L' else 1,A))]:
   correction=sym({(a,b):r.scale(-1,r.add(*(r.mul(eps[a,c,d],r.mul(g[e,c,b],S[e,d])) for c,d,e in product(I,repeat=3)))) for a,b in product(I,repeat=2)})
   test(correction,target,'connection_curl_'+name+'_'+tag)
  # e1_a^cd=p_a^e eps_e^cd - tr(p)/2 eps_a^cd.
  tr=r.trace(p)[()]
  e1={(a,c,d):r.add(*(r.mul(p[a,e],eps[e,c,d]) for e in I),r.scale(-Q(1,2),r.mul(tr,eps[a,c,d]))) for a,c,d in product(I,repeat=3)}
  for tag,S,co in [('trace',r.tmul(Y,P),-1),('hessian',h,1)]:
   dS=r.covd(S)
   correction=sym({(a,b):r.add(*(r.mul(e1[a,c,d],dS[c,b,d]) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)})
   test(correction,zero if name=='L' else r.tscale(co,A),'volume_raise_correction_'+name+'_'+tag)
 test(sym({(a,b):r.add(*(r.mul(eps[a,c,d],r.mul(h[b,c],v[d,])) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)}),r.tscale(-1,A),'radial_shape_gradient_curl')
 if coeff != [Q(2),Q(-2)]:raise ValueError('unexpected axial TT convention '+str(coeff))
 return {'status':'PASS_EXACT_SCALAR_SCREEN_MAGNETIC_ANGULAR_REDUCTION','orientation':'epsilon_abc=sum_j EPS[j,a,b,c]*X[j], EPS[0,1,2,3]=+1','checks':checks,'axial_TT_coefficients':[str(x) for x in coeff],'axial_TT_norms':[str(r.inner(t,t)) for t in r.TT],'axial_norm_squared':str(r.inner(A,A)),'all_point_polynomial_reduction':True,'angular_source_sha256':ANGULAR_SHA256}

def main():
 p=argparse.ArgumentParser();p.add_argument('--angular',type=Path,default=DEFAULT);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 validate_output(args.output)
 res=run(args.angular)
 target=validate_output(args.output)
 with target.open('x',encoding='utf-8') as stream:stream.write(json.dumps(res,indent=2)+'\n')
 print(json.dumps(res))
if __name__=='__main__':main()
