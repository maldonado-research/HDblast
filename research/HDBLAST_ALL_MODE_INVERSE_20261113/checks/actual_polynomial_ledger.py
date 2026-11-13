#!/usr/bin/env python3
"""Actual polynomial normal-ordering ledger for the scalar all-mode proof.

R is expanded exactly.  A^6 R uses a certified coefficient-norm recurrence
instead of expanding huge multivariate polynomials.  This supports coefficient
accounting only; the analytic graph and norm estimates remain separate.
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

variables=s.symbols("x b a p eta")
x,b,a,p,eta=variables


def canonical(value):return (json.dumps(value,sort_keys=True,indent=2)+"\n").encode()


def normal(P):
    answer={}
    for key,value in P.items():
        expanded=s.expand(value)
        if expanded!=0:answer[key]=expanded
    return answer


def add(*ops):
    answer={}
    for op in ops:
        for key,value in op.items():answer[key]=answer.get(key,0)+value
    return normal(answer)


def scale(op,factor):return normal({key:factor*value for key,value in op.items()})


def norm_degree(poly):
    P=s.Poly(poly,*variables)
    return sum((abs(c) for c in P.coeffs()),s.Rational(0)),max((sum(k) for k,_ in P.terms()),default=0)


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        value=bool(ok);check("actual_polynomials/"+name,value)
        records.append({"name":name,"status":"pass" if value else "fail"})
        if not value:raise ValueError("actual polynomial ledger rejects "+name)
    need("schema",claims["schema"]=="hdblast-all-mode-scalar-coefficients-claims-v1")
    need("sympy_pin",s.__version__=="1.14.0")
    need("review_status",claims["review_status"]=="accepted_conditional_analytic_inputs")
    retained=[s.Rational(z) for z in ("-2/27","0","14/9","50/27","-1/6","-4/9","-2/27")]
    U=s.Add(*(s.Rational(z)*eta**j for j,z in enumerate(claims["U_coefficients"])))
    need("retained_actual_potential",s.expand(U-s.Add(*(z*eta**j for j,z in enumerate(retained))))==0)
    q=s.Rational(1,10**16)*(1-x*x)**2
    vf=(-s.Rational(1,10**10)*(1-x*x)**2*b,-q*a*b,
        q*(-b*b-p*p/3),q*(s.diff(U,eta)-4*a*p),q*p)
    vf_norm=max(norm_degree(v)[0] for v in vf)
    vf_degree=max(norm_degree(v)[1] for v in vf)
    V=s.Rational(claims["vectorfield_norm_ceiling"])
    need("vectorfield_coefficient_norm",vf_norm<=V)
    need("vectorfield_degree",vf_degree<=claims["vectorfield_degree_ceiling"])
    need("required_vectorfield_ceiling",V==s.Rational(1,10**8) and claims["vectorfield_degree_ceiling"]==9)
    def D(poly):return s.expand(sum(v*s.diff(poly,z) for z,v in zip(variables,vf)))
    d=4*s.Rational(1,10**10)*x*(1-x*x)*b
    c=d/2+s.Rational(3,2)*q*a+2*s.Rational(1,10**10)*x*b
    relative=d-q*a
    h=q*b
    Rt=3*b*b+p*p/3-3*a*a
    alpha={(0,1):s.Rational(claims["alpha1_factor"])*a*q,
           (0,0):s.Rational(3,2)*a*q*c+s.Rational(3,2)*q*q*(a*a-b*b)}
    beta={(0,1):c+d/2-a*q,
          (0,0):-(c*c+d*c-D(c))/2+a*q*c-q*q*Rt/2}
    e=c+a*q
    basic_coefficients={"unitary_drift":c,"q_y":d,"q_a":q*a,"q_over_rho":h,"relative_angular_drift":relative}
    basic_records={}
    for name,poly in basic_coefficients.items():
        coefficient_norm,degree=norm_degree(poly)
        envelope=coefficient_norm
        for k in range(14):envelope*=V*(degree+8*k)
        need("basic_coefficient_jet/"+name,coefficient_norm<1 and envelope<1 and
             all(V*(degree+8*k)<1 for k in range(14)))
        basic_records[name]={"coefficient_norm":str(coefficient_norm),"degree":degree,"D14_envelope":str(envelope)}
    alpha_norm=sum(norm_degree(z)[0] for z in alpha.values())
    beta_norm=sum(norm_degree(z)[0] for z in beta.values())
    alpha_degree=max(norm_degree(z)[1] for z in alpha.values())
    beta_degree=max(norm_degree(z)[1] for z in beta.values())
    need("alpha_actual_coefficient_norm",alpha_norm<=s.Rational(claims["row_coefficient_norm_ceiling"]))
    need("beta_actual_coefficient_norm",beta_norm<=s.Rational(claims["row_coefficient_norm_ceiling"]))
    need("alpha_actual_degree",alpha_degree<=claims["alpha_degree_ceiling"])
    need("beta_actual_degree",beta_degree<=claims["beta_degree_ceiling"])
    def relative_derivative(poly,power):return s.expand(D(poly)+power*relative*poly)
    def differentiated(poly,power,count):
        for _ in range(count):poly=relative_derivative(poly,power)
        return poly
    def compose(P,Q):
        answer={}
        for (bp,ap),u in P.items():
            for (bq,aq),v in Q.items():
                for k in range(ap+1):
                    key=(bp+bq,ap-k+aq)
                    answer[key]=answer.get(key,0)+u*math.comb(ap,k)*differentiated(v,bq,k)
        return normal(answer)
    def adjoint(P):
        answer={}
        for (power,order),u in P.items():
            for k in range(order+1):
                key=(power,k)
                answer[key]=answer.get(key,0)+(-1)**order*math.comb(order,k)*differentiated(u,power,order-k)
        return normal(answer)
    m={(2,0):s.Integer(1)}
    Bbase={(0,2):-s.Rational(1,2),(2,0):s.Rational(1,3)}
    # Oracle radial row comes directly from the retained full row before the
    # proposed alpha split. In particular it does not consume alpha1_factor.
    Arow={(2,0):s.Rational(1,2),(0,1):-s.Rational(3,2)*a*q,
          (0,0):s.Rational(3,2)*a*q*c+s.Rational(3,2)*q*q*(a*a-b*b)}
    Brow=add(Bbase,beta)
    Cbase={(0,1):s.Integer(1),(0,0):-e}
    direct=add(compose(adjoint(Arow),Arow),scale(compose(adjoint(Brow),Brow),3),
               scale(compose(compose(adjoint(Cbase),m),Cbase),s.Rational(1,2)),
               {(4,0):s.Rational(1,6),(2,0):-h*h/2})
    P={(0,2):-s.Integer(1),(2,0):s.Integer(1)}
    oracle_R=add(direct,scale(compose(P,P),-s.Rational(3,4)))
    # Independent expansion of the proposed remainder, including derivatives
    # of zeta^2 and of e*zeta^2 through the actual polynomial vector field.
    production_R=add({(2,0):differentiated(s.Integer(1),2,2)/4},
        scale(add(compose(m,alpha),compose(adjoint(alpha),m)),s.Rational(1,2)),
        compose(adjoint(alpha),alpha),
        scale(add(compose(Bbase,beta),compose(adjoint(beta),Bbase),compose(adjoint(beta),beta)),3),
        {(2,0):(relative_derivative(e,2)+e*e-h*h)/2})
    need("actual_normal_ordered_remainder_identity",not add(production_R,scale(oracle_R,-1)))
    coefficient_records=[]
    Rnorm=s.Rational(0);Rdegree=0;Rorder=0
    for (power,order),poly in sorted(oracle_R.items()):
        norm,degree=norm_degree(poly)
        Rnorm+=norm;Rdegree=max(Rdegree,degree);Rorder=max(Rorder,power+order)
        coefficient_records.append({"zeta_power":power,"D_order":order,
                                    "coefficient_norm_exact":str(norm),"degree":degree,
                                    "polynomial_monomial_count":len(s.Poly(poly,*variables).terms())})
    need("R_actual_degree",Rdegree<=claims["R_degree_ceiling"])
    need("R_actual_coefficient_sum",Rnorm<s.Rational(claims["R_norm_ceiling"]))
    need("R_actual_weighted_order",Rorder<=claims["R_weighted_order_ceiling"])
    need("R_target_ceiling",claims["R_degree_ceiling"]==50 and s.Rational(claims["R_norm_ceiling"])==10**6 and claims["R_weighted_order_ceiling"]==3)
    # Certify the left-A recurrence for every possible degree/power used,
    # rather than asserting a universal coefficient multiplier.  If H=C*zeta^b,
    # then (D H)/zeta^b=D C+b*relative*C, and repeat once for D^2 H.
    relnorm,reldegree=norm_degree(relative)
    for degree in range(201):
        for power in range(17):
            first=V*degree+power*relnorm
            # D C increases degree by8; multiplication by relative increases
            # degree by<=5, so a common degree bound is degree+8.
            second=first*(V*(degree+8)+power*relnorm)
            need_condition=first<=1 and second<=1
            if not need_condition:raise ValueError("left-A derivative coefficient bound failed")
    need("left_A_first_two_relative_derivatives_contract",True)
    # A(C zeta^b D^a) has groups C, -C D^2, -2C' D,
    # -C'', and C zeta^2. Their norms sum to at most6||C||.
    multiplier=claims["A_left_coefficient_factor"]
    increment=claims["A_left_degree_increment"]
    need("left_A_group_weight",multiplier==1+1+2+1+1)
    need("left_A_degree_increment",increment==16)
    A6budget=multiplier**6
    A6Rbudget=A6budget*s.Rational(claims["R_norm_ceiling"])
    need("A6_coefficient_budget",A6budget<s.Rational(claims["A6_coefficient_ceiling"]))
    need("A6R_coefficient_budget",A6Rbudget<s.Rational(claims["A6R_coefficient_ceiling"]))
    need("A6R_degree_budget",claims["R_degree_ceiling"]+6*increment==claims["A6R_degree_ceiling"]==146)
    need("A6R_weighted_order",claims["R_weighted_order_ceiling"]+12==claims["A6R_weighted_order_ceiling"]==15)
    need("relative_jet_factor",all(s.Rational(1,10**8)*(8*j+power)<s.Rational(4,10**6)
                                  for j in range(32) for power in range(33)))
    return {"schema":"hdblast-all-mode-actual-polynomial-ledger-v1","status":"pass",
        "review_status":claims["review_status"],"claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),
        "checks":records,"exact":{
            "vectorfield_norm":str(vf_norm),"vectorfield_degree":vf_degree,
            "alpha_coefficient_norm":str(alpha_norm),"alpha_degree":alpha_degree,
            "beta_coefficient_norm":str(beta_norm),"beta_degree":beta_degree,
            "R_coefficient_norm_sum":str(Rnorm),"R_degree":Rdegree,"R_weighted_order":Rorder,
            "relative_coefficient_norm":str(relnorm),"relative_coefficient_degree":reldegree,
            "A6_coefficient_budget":str(A6budget),"A6R_coefficient_budget":str(A6Rbudget),
            "A6R_degree_budget":146,"A6R_weighted_order_budget":15},
        "R_coefficients":coefficient_records,"basic_coefficient_jets":basic_records,
        "scope":"Actual exact polynomial/coefficient ledger plus certified norm-degree recurrences. Accepted analytic graph, domain, all-mode coercivity, and norm-bridge proofs remain separate pinned inputs."}


def controls(claims):
    records=[]
    for name,value in (("alpha1_factor","-1"),("vectorfield_norm_ceiling","1/1000000000000"),
                       ("R_degree_ceiling",1),("R_norm_ceiling","1/100000000000000000000"),
                       ("R_weighted_order_ceiling",1),("A_left_coefficient_factor",2)):
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),
            "rejected_by":"fresh actual polynomial or exact coefficient-budget oracle","diagnostic":str(exc)})
        else:raise ValueError("surviving polynomial-ledger control: "+name)
    return {"schema":"hdblast-all-mode-polynomial-controls-v1","status":"pass","controls":records}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--claims",type=Path,default=Path(__file__).with_name("coefficient_claims.json"))
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
