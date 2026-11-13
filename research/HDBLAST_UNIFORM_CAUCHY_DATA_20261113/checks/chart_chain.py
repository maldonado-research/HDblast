#!/usr/bin/env python3
"""Independent exact nonzero-state chart/time-velocity chain witnesses.

A genuine square-root chart is composed with a non-diagonal cubic polynomial
map. Literal one/two-parameter composition is differentiated independently of
the proposed chain formula. No actual flow-jet or Sobolev bound is inferred.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import itertools
import json
from pathlib import Path
import sys
import sympy as S
import mpmath

FACTORS=("velocity_second_P_second_Q_factor","velocity_second_P_first_Q_first_factor",
 "velocity_second_third_differential_factor","polarized_mixed_cross_factor","composed_root_third_differential_factor")
F={k:S.Symbol("factor_"+k) for k in FACTORS}
MUTATIONS=(("velocity_second_P_second_Q_factor","0"),("velocity_second_P_first_Q_first_factor","1"),
           ("velocity_second_third_differential_factor","0"),("polarized_mixed_cross_factor","0"),
           ("composed_root_third_differential_factor","0"))


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def clean(v):return S.cancel(S.expand(v))


def build_checks():
    u=S.symbols("u0:3");r=S.symbols("r0:3");n=3
    root=S.sqrt(1+u[0])-1
    rv=[root,u[1],u[2]]
    # All three levels are non-diagonal and nonseparable. Symmetric numerical
    # tensors yield exact rational cubic maps; these are auxiliary witnesses,
    # not substitutes for the retained boundary-flow derivatives.
    linear=S.Matrix(n,n,lambda i,j:S.Rational(2*i+3*j+1,7))
    def b(i,j,k):return S.Rational((i+1)*(j+k+3)+j*k+1,11)
    def c(i,j,k,l):return S.Rational((i+2)*(j+k+l+2)+j*k+k*l+j*l+1,13)
    poly=[sum(linear[i,j]*r[j] for j in range(n))+
          sum(b(i,j,k)*r[j]*r[k]/2 for j,k in itertools.product(range(n),repeat=2))+
          sum(c(i,j,k,l)*r[j]*r[k]*r[l]/6 for j,k,l in itertools.product(range(n),repeat=3)) for i in range(n)]
    T=[v.subs(dict(zip(r,rv)),simultaneous=True) for v in poly]
    p0=S.Matrix([S.Rational(20001,10**8),S.Rational(1,5000),-S.Rational(1,6000)])
    q0=S.Matrix([S.Rational(2,7),-S.Rational(1,11),S.Rational(3,13)])
    at=dict(zip(u,p0))
    d1=S.Matrix(n,n,lambda i,j:S.diff(T[i],u[j]).subs(at))
    d2=[[[S.diff(T[i],u[j],u[k]).subs(at) for k in range(n)] for j in range(n)] for i in range(n)]
    d3=[[[[S.diff(T[i],u[j],u[k],u[l]).subs(at) for l in range(n)] for k in range(n)] for j in range(n)] for i in range(n)]
    vectors=[S.Matrix(S.symbols(name+"0:3")) for name in ("a","b","c","d","m","n")]
    a,bv,cv,dv,mv,nv=vectors
    def hessian(left,right):return S.Matrix([sum(d2[i][j][k]*left[j]*right[k] for j,k in itertools.product(range(n),repeat=2)) for i in range(n)])
    rat=dict(zip(r,[v.subs(at) for v in rv]))
    flow_first_root=S.Matrix([S.diff(v,r[0]).subs(rat) for v in poly])
    root3=S.diff(root,u[0],3).subs(at)
    def third(left,middle,right):
        exact=S.Matrix([sum(d3[i][j][k][l]*left[j]*middle[k]*right[l]
                      for j,k,l in itertools.product(range(n),repeat=3)) for i in range(n)])
        return exact+(F["composed_root_third_differential_factor"]-1)*flow_first_root*root3*left[0]*middle[0]*right[0]
    e,f=S.symbols("e f")
    curve=p0+e*a+e*e*bv/2
    velocity=q0+e*cv+e*e*dv/2
    derivative=S.Matrix(n,n,lambda i,j:S.diff(T[i],u[j]))
    literal=derivative.subs(dict(zip(u,curve)),simultaneous=True)*velocity
    one_oracle=S.Matrix([clean(S.diff(v,e,2).subs(e,0)) for v in literal])
    one_production=(d1*dv+F["velocity_second_P_second_Q_factor"]*hessian(bv,q0)+
       F["velocity_second_P_first_Q_first_factor"]*hessian(a,cv)+
       F["velocity_second_third_differential_factor"]*third(a,a,q0))
    polarized_curve=p0+e*a+f*bv+e*f*mv
    polarized_velocity=q0+e*cv+f*dv+e*f*nv
    literal_polarized=derivative.subs(dict(zip(u,polarized_curve)),simultaneous=True)*polarized_velocity
    two_oracle=S.Matrix([clean(S.diff(v,e,f).subs({e:0,f:0})) for v in literal_polarized])
    two_production=(d1*nv+F["velocity_second_P_second_Q_factor"]*hessian(mv,q0)+
        F["polarized_mixed_cross_factor"]*(hessian(a,dv)+hessian(bv,cv))+
        F["velocity_second_third_differential_factor"]*third(a,bv,q0))
    # Both potentially omitted terms are nonzero at the auxiliary nonzero base.
    if all(clean(v)==0 for v in hessian(bv,q0)):raise ValueError("P second/Q base witness degenerate")
    if all(clean(v)==0 for v in third(a,a,q0).subs({F["composed_root_third_differential_factor"]:1})):raise ValueError("third-differential witness degenerate")
    return [("nonzero_state_second_velocity_chain",list(one_production-one_oracle)),
            ("polarized_second_velocity_chain",list(two_production-two_oracle))],{
      "auxiliary_state_dimension":n,"root_state":str(p0[0]),"base_velocity_nonzero":True,
      "root_first_second_third":[str(S.diff(root,u[0],k).subs(at)) for k in (1,2,3)],
      "formula":"DT Q'' + D2T[P'',Q] + 2 D2T[P',Q'] + D3T[P',P',Q]",
      "scope":"Exact chain identities for a square root composed with an auxiliary non-diagonal cubic map at nonzero state/velocity. No actual boundary-flow derivative, analytic Nemytskii constant, or uniform Cauchy bound is proved."}


def evaluate(claims,check):
    if S.__version__!="1.14.0" or mpmath.__version__!="1.3.0":raise ValueError("requires pinned exact engine")
    if claims.get("schema")!="hdblast-uniform-cauchy-chart-chain-claims-v1" or set(claims)!={"schema",*FACTORS}:raise ValueError("bad chart-chain claims")
    checks,meta=build_checks();sub={F[k]:S.Rational(claims[k]) for k in FACTORS};records=[];bad=[]
    for name,residuals in checks:
        failures=[clean(v.subs(sub)) for v in residuals];failures=[v for v in failures if v!=0]
        check("chart_chain/"+name,not failures)
        records.append({"name":name,"components":len(residuals),"status":"fail" if failures else "pass"})
        if failures:bad.append({"name":name,"nonzero_exact_residual":str(failures[0])})
    if bad:raise ValueError("chart-chain oracle rejects claims: "+json.dumps(bad,sort_keys=True))
    return {"schema":"hdblast-uniform-cauchy-chart-chain-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,"metadata":meta}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
          "rejected_by":"fresh literal nonzero-state chart-composition oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving chart-chain fault: "+name)
    return {"schema":"hdblast-uniform-cauchy-chart-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("chart_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        c=json.loads(args.claims.read_text());r=evaluate(c,lambda name,ok:None);payload=canonical(r)
        if args.report:args.report.write_bytes(payload)
        if args.controls:args.controls.write_bytes(canonical(controls(c)))
        sys.stdout.buffer.write(payload);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1
