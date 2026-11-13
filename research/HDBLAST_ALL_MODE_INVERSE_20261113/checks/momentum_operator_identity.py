#!/usr/bin/env python3
"""Exact all-sector momentum normal-operator decomposition, no inverse claim."""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as s

t=s.Symbol("s",real=True)


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def normal(P):return {j:s.expand(v) for j,v in P.items() if s.expand(v)!=0}
def add(*operators):
    out={}
    for P in operators:
        for j,v in P.items():out[j]=out.get(j,0)+v
    return normal(out)
def scale(P,k):return normal({j:k*v for j,v in P.items()})
def compose(P,Q):
    out={}
    for i,a in P.items():
        for j,b in Q.items():
            for k in range(i+1):
                order=i-k+j
                out[order]=out.get(order,0)+a*math.comb(i,k)*s.diff(b,t,k)
    return normal(out)
def adjoint(P):
    out={}
    for i,a in P.items():
        for k in range(i+1):out[k]=out.get(k,0)+(-1)**i*math.comb(i,k)*s.diff(a,t,i-k)
    return normal(out)


def evaluate(claims,check):
    if s.__version__!="1.14.0":raise ValueError("requires sympy==1.14.0")
    if claims["schema"]!="hdblast-all-mode-momentum-normal-identity-claims-v1":raise ValueError("wrong momentum identity schema")
    f={k:s.Rational(v) for k,v in claims.items() if k!="schema"}
    b,j,c,zeta=[s.Function(name)(t) for name in ("b","j","c","zeta")]
    h=s.diff(zeta,t)/zeta
    # Actual first three normalized sym-gradient rows. The fourth row's
    # squared contribution is represented directly as multiplication, removing
    # an unnecessary square root and preserving its trace-free correction.
    rows=[[{1:2,0:-2*b},{}],
          [{0:s.sqrt(2)*zeta},{1:s.sqrt(2),0:-s.sqrt(2)*(b+j)}],
          [{0:s.sqrt(12)*j},{0:-s.sqrt(12)*zeta/3}]]
    oracle=[[add(*(compose(adjoint(row[a]),row[d]) for row in rows)) for d in range(2)] for a in range(2)]
    oracle[1][1]=add(oracle[1][1],{0:s.Rational(8,3)*zeta*zeta-8*c*c})
    P={2:-s.Integer(1),0:zeta*zeta}
    G=[{1:s.Integer(1)},{0:-zeta}]
    v1=f["v1_b_squared_factor"]*b*b+f["v1_j_squared_factor"]*j*j+f["v1_b_derivative_factor"]*s.diff(b,t)
    v2=f["v2_combined_squared_factor"]*(b+j)**2+f["v2_c_squared_factor"]*c*c+f["v2_combined_derivative_factor"]*s.diff(b+j,t)
    potential=[v1,v2]
    production=[[scale(compose(adjoint(G[a]),G[d]),f["principal_G_gram_factor"]) for d in range(2)] for a in range(2)]
    for a in range(2):production[a][a]=add(production[a][a],scale(P,f["principal_P_factor"]),{0:potential[a]})
    for a,d in ((0,1),(1,0)):
        production[a][d]=add(production[a][d],{0:f["offdiagonal_factor"]*(h+b+3*j)*zeta})
    residuals=[]
    for a in range(2):
        for d in range(2):
            delta=add(production[a][d],scale(oracle[a][d],-1))
            for order in range(3):residuals.append(("poloidal_"+str(a)+str(d)+"_D"+str(order),delta.get(order,0)))
    # Coexact mode: its angular tensor term is 2(mu-4)c^2 =2zeta^2-8c^2.
    coexactrow={1:s.sqrt(2),0:-s.sqrt(2)*(b+j)}
    coexact_oracle=add(compose(adjoint(coexactrow),coexactrow),{0:2*zeta*zeta-8*c*c})
    coexact_production=add(scale(P,2),{0:v2})
    mono_oracle=add(compose(adjoint(rows[0][0]),rows[0][0]),{0:12*j*j})
    mono_production=add({2:-s.Integer(4)},{0:v1})
    for name,delta in (("coexact",add(coexact_production,scale(coexact_oracle,-1))),
                       ("monopole",add(mono_production,scale(mono_oracle,-1)))):
        for order in range(3):residuals.append((name+"_D"+str(order),delta.get(order,0)))
    # The highest commutative symbols cancel because A^6 is scalar. The
    # single-commutation symbol then has order thirteen, not fourteen. This
    # checks the actual mixed matrix principal symbol, including zeta D.
    dv,zv,hv=s.symbols("D_symbol zeta_symbol relative_derivative",real=True)
    Asymbol=(-dv*dv+zv*zv)**6
    Lsymbol=s.Matrix([[-4*dv*dv+2*zv*zv,2*zv*dv],
                      [-2*zv*dv,-2*dv*dv+4*zv*zv]])
    top14=Lsymbol*(Asymbol*s.eye(2))-(Asymbol*s.eye(2))*Lsymbol
    for index,value in enumerate(top14):
        residuals.append(("commutator_order14_symbol_component_"+str(index),s.expand(value)))
    top13=Lsymbol.applyfunc(lambda value:s.expand(s.diff(value,dv)*hv*zv*s.diff(Asymbol,zv)-
                                                  s.diff(Asymbol,dv)*hv*zv*s.diff(value,zv)))
    actual_top=max(s.Poly(value,dv,zv).total_degree() for value in top13 if value!=0)
    residuals.append(("commutator_exact_top_weighted_order",actual_top-f["commutator_top_weighted_order"]))
    residuals.append(("commutator_top_symbol_nonzero_witness",top13[0,1].subs({dv:0,zv:1,hv:1})-24))
    records=[];failures=[]
    for name,value in residuals:
        value=s.expand(value);ok=value==0;check("momentum_identity/"+name,ok)
        records.append({"name":name,"status":"pass" if ok else "fail"})
        if not ok:failures.append({"name":name,"exact_residual":str(value)})
    if failures:raise ValueError("momentum normal decomposition fails: "+json.dumps(failures,sort_keys=True))
    return {"schema":"hdblast-all-mode-momentum-normal-identity-v1","status":"pass",
        "claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
        "scope":"Exact generic real-coefficient operator identities for poloidal, both coexact families, and monopole. No graph estimate, momentum H14 inverse bound, or analytic domain assertion follows."}


def controls(claims):
    records=[]
    for name in (k for k in claims if k!="schema"):
        bad=copy.deepcopy(claims);bad[name]=str(s.Rational(bad[name])+1);payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name+"_plus_one","self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
            "rejected_by":"fresh exact momentum operator oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving momentum identity fault: "+name)
    return {"schema":"hdblast-all-mode-momentum-normal-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("momentum_identity_claims.json"))
    parser.add_argument("--report",type=Path);parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text());report=evaluate(claims,lambda name,ok:None)
        payload=canonical(report)
        if args.report:args.report.write_bytes(payload)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(payload);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":raise SystemExit(main())
