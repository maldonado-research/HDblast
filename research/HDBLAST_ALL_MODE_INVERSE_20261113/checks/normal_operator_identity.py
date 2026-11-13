#!/usr/bin/env python3
"""Exact formal differential-operator identity underlying scalar all-mode work.

Rows and formal adjoints are composed independently of the proposed square
plus remainder formula.  Coefficients are arbitrary real smooth functions.
This is an algebra identity, not a regularity estimate or inverse bound.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as s

x=s.Symbol("s",real=True)


def canonical(value):
    return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()


def simplify(P):
    return {j:s.expand(a) for j,a in P.items() if s.expand(a)!=0}


def add(*operators):
    answer={}
    for P in operators:
        for j,a in P.items():answer[j]=answer.get(j,0)+a
    return simplify(answer)


def scale(P,a):return simplify({j:a*b for j,b in P.items()})


def compose(P,Q):
    """Leibniz normal order for P(Q u), coefficients on the left."""
    answer={}
    for i,a in P.items():
        for j,b in Q.items():
            for k in range(i+1):
                order=i-k+j
                answer[order]=answer.get(order,0)+a*math.comb(i,k)*s.diff(b,x,k)
    return simplify(answer)


def adjoint(P):
    """Flat-ds formal adjoint for real coefficients, including derivative terms."""
    answer={}
    for i,a in P.items():
        for k in range(i+1):
            answer[k]=answer.get(k,0)+(-1)**i*math.comb(i,k)*s.diff(a,x,i-k)
    return simplify(answer)


def evaluate(claims,check):
    if s.__version__!="1.14.0":raise ValueError("requires pinned sympy==1.14.0")
    if claims["schema"]!="hdblast-all-mode-scalar-normal-identity-claims-v1":raise ValueError("wrong identity schema")
    a0,a1,b0,b1,m,e,h=[s.Function(name)(x) for name in ("alpha0","alpha1","beta0","beta1","m","e","h")]
    alpha={0:a0,1:a1};beta={0:b0,1:b1}
    Bbase={2:-s.Rational(1,2),0:m/3}
    Arow=add({0:m/2},alpha)
    Brow=add(Bbase,beta)
    Cbase={1:s.Integer(1),0:-e}
    # Oracle starts with the energy rows of the exact normal operator.
    oracle=add(compose(adjoint(Arow),Arow),scale(compose(adjoint(Brow),Brow),3),
               scale(compose(compose(adjoint(Cbase),{0:m}),Cbase),s.Rational(1,2)),
               {0:m*m/6-h*h*m/2})
    f={name:s.Rational(value) for name,value in claims.items() if name!="schema"}
    model={2:-s.Integer(1),0:m}
    alpha_cross=add(compose({0:m},alpha),compose(adjoint(alpha),{0:m}))
    beta_cross=add(compose(Bbase,beta),compose(adjoint(beta),Bbase),compose(adjoint(beta),beta))
    production=add(scale(compose(model,model),f["model_square_factor"]),
        {0:f["m_second_derivative_factor"]*s.diff(m,x,2)},
        scale(alpha_cross,f["alpha_cross_factor"]),compose(adjoint(alpha),alpha),
        scale(beta_cross,f["beta_cross_factor"]),
        {0:f["e_cross_factor"]*(s.diff(m*e,x)+e*e*m)+f["h_squared_factor"]*h*h*m})
    residual=add(production,scale(oracle,-1))
    records=[]
    for j in range(5):
        ok=s.expand(residual.get(j,0))==0
        check("normal_identity/coefficient_D"+str(j),ok)
        records.append({"coefficient":"D^"+str(j),"status":"pass" if ok else "fail"})
    if residual:raise ValueError("scalar normal identity fails: "+str(residual))
    return {"schema":"hdblast-all-mode-scalar-normal-identity-v1","status":"pass",
        "claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
        "oracle":"Arow* Arow + 3 Brow* Brow + (D-e)* m (D-e)/2 + m^2/6 - h^2 m/2",
        "production":"3/4(-D^2+m)^2 + m''/4 + (m alpha+alpha* m)/2+alpha*alpha +3(Bbase beta+beta* Bbase+beta*beta)+((m e)'+e^2 m-h^2 m)/2",
        "scope":"Exact real-coefficient formal differential-operator algebra only; no norm estimate, all-mode inverse constant, or nonlinear claim."}


def controls(claims):
    records=[]
    for name in (k for k in claims if k!="schema"):
        bad=copy.deepcopy(claims);bad[name]=str(s.Rational(bad[name])+1)
        payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:
            records.append({"fault":name+"_plus_one","self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
                "rejected_by":"fresh generic differential-operator oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving normal-identity factor mutation: "+name)
    return {"schema":"hdblast-all-mode-scalar-normal-identity-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("normal_identity_claims.json"))
    parser.add_argument("--report",type=Path)
    parser.add_argument("--controls",type=Path)
    args=parser.parse_args(argv)
    try:
        claims=json.loads(args.claims.read_text());report=evaluate(claims,lambda name,ok:None)
        encoded=canonical(report)
        if args.report:args.report.write_bytes(encoded)
        if args.controls:args.controls.write_bytes(canonical(controls(claims)))
        sys.stdout.buffer.write(encoded);return 0
    except (ValueError,KeyError,TypeError,OSError) as exc:
        sys.stderr.write(str(exc)+"\n");return 1


if __name__=="__main__":raise SystemExit(main())
