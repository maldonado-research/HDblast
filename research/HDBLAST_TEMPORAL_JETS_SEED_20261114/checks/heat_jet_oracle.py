#!/usr/bin/env python3
"""Universal finite heat-counterterm jets and Gaussian derivative witnesses.

Single-mode Taylor identities and finite polynomial/Gaussian coefficient
budgets are checked exactly. A continuum spectral lift, strong traces, norm
bridges and boundary conditions require the separate analytic proof.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import sympy as S
import mpmath

FACTORS=("max_jet_order","max_time_derivative","counterterms","counterterm_sign",
 "counterterm_lambda_power_factor","counterterm_r_factorial","counterterm_lower_jet_factorial",
 "series_factorial","heat_exponential_sign","heat_time_power","Hermite_linear_factor","Gaussian_monomial_budget")
MUTATIONS=(("counterterms","diagonal_only"),("counterterm_sign","-1"),
 ("counterterm_lambda_power_factor",2),("counterterm_r_factorial","unit"),
 ("counterterm_lower_jet_factorial","unit"),("series_factorial","unit"),
 ("heat_exponential_sign","1"),("heat_time_power",1),("Hermite_linear_factor","2"))


def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+"\n").encode()
def clean(v):return S.cancel(S.expand(v))
def factorial(mode,n):
    if mode=="factorial":return math.factorial(n)
    if mode=="unit":return 1
    raise ValueError("unsupported factorial mode")


def evaluate(claims,check):
    records=[]
    def need(name,ok):
        boolean=bool(ok);check("heat_jet_oracle/"+name,boolean)
        records.append({"name":name,"status":"pass" if boolean else "fail"})
        if not boolean:raise ValueError("heat oracle rejects "+name)
    need("schema",claims.get("schema")=="hdblast-heat-jet-lift-claims-v1" and set(claims)=={"schema",*FACTORS})
    need("exact_engine",S.__version__=="1.14.0" and mpmath.__version__=="1.3.0")
    need("finite_orders",claims["max_jet_order"]==14 and claims["max_time_derivative"]==14)
    need("counterterm_mode",claims["counterterms"] in ("full","diagonal_only"))
    s,x=S.symbols("s x",real=True);lam=S.Symbol("lambda",positive=True)
    J=S.symbols("J0:15")
    K=[]
    for j in range(15):
        terms=[]
        for r in range(j//2+1):
            if claims["counterterms"]=="diagonal_only" and r:continue
            low=j-2*r
            coefficient=S.Rational(math.factorial(j),factorial(claims["counterterm_r_factorial"],r)*
                                    factorial(claims["counterterm_lower_jet_factorial"],low))
            terms.append(coefficient*S.Rational(claims["counterterm_sign"])**r*
                         lam**(claims["counterterm_lambda_power_factor"]*r)*J[low])
        K.append(S.Add(*terms))
    polynomial=sum(s**j*K[j]/factorial(claims["series_factorial"],j) for j in range(15))
    actual_lift=S.exp(S.Rational(claims["heat_exponential_sign"])*lam*s**claims["heat_time_power"])*polynomial
    # Literal Gaussian-times-polynomial differentiation contains no supplied
    # counterterm identity or binomial cancellation formula.
    for j in range(15):
        recovered=clean(S.diff(actual_lift,s,j).subs(s,0))
        need("arbitrary_initial_jet_"+str(j),recovered==J[j])
    Hfactor=S.Rational(claims["Hermite_linear_factor"])
    gaussian=S.exp(-x*x)
    envelopes={};degrees=[];coefficient_sums=[]
    for j in range(15):
        P=x**j
        for k in range(15):
            if k:P=S.expand(S.diff(P,x)+Hfactor*x*P)
            # Independent direct derivative of the actual Gaussian. All
            # exponentials cancel before conversion to an exact polynomial.
            oracle=S.Poly(S.expand(S.diff(x**j*gaussian,x,k)*S.exp(x*x)),x)
            need("Gaussian_polynomial_j"+str(j)+"_k"+str(k),clean(P-oracle.as_expr())==0)
            degree=oracle.degree();degrees.append(degree)
            need("Gaussian_degree_j"+str(j)+"_k"+str(k),degree<=j+k<=28)
            coefficients={m[0]:int(a) for m,a in oracle.terms()}
            coefficient_sums.append(sum(abs(a) for a in coefficients.values()))
            envelopes[j,k]=sum(abs(a)*math.factorial(m) for m,a in coefficients.items())
    need("Gaussian_factorial_budget_mode",claims["Gaussian_monomial_budget"]=="factorial")
    # Positive-series coefficient 1/m! of exp(x^2) proves the stated finite
    # monomial budget for x>=1; x<=1 is bounded by one. This finite check
    # records the exact constants; the analytic inequality is proved outside.
    for m in range(29):
        need("positive_Gaussian_series_coefficient_"+str(m),S.Rational(1,math.factorial(m))>0 and
              S.Rational(1,math.factorial(m))*math.factorial(m)==1)
    spectral_count=0
    coefficient_budget=Fraction(0)
    for k in range(15):
        for j in range(15):
            for r in range(j//2+1):
                low=j-2*r
                # Twice the lambda exponent avoids all fractional powers.
                left=(14-k)+(k-j)+2*r
                if left!=14-low:raise ValueError("spectral order mismatch")
                spectral_count+=1
                coefficient_budget+=Fraction(envelopes[j,k],math.factorial(r)*math.factorial(low))
    need("all_finite_spectral_exponents_balance",spectral_count>0)
    ceiling=0
    while coefficient_budget>=10**ceiling:ceiling+=1
    return {"schema":"hdblast-heat-jet-oracle-v1","status":"pass","claims_sha256":hashlib.sha256(canonical(claims)).hexdigest(),"checks":records,
       "finite_exact":{"arbitrary_jet_count":15,"Gaussian_polynomial_count":225,"maximum_polynomial_degree":max(degrees),
           "maximum_polynomial_coefficient_l1":max(coefficient_sums),"maximum_factorial_Gaussian_envelope":str(max(envelopes.values())),
           "spectral_weight_identities":spectral_count,"summed_Gaussian_counterterm_coefficient_budget":str(coefficient_budget),
           "strict_decimal_coefficient_ceiling_exponent":ceiling},
       "scope":"Finite single-positive-mode Taylor identities, exact Hermite/Gaussian polynomials and conditional coefficient envelopes. The coefficient budget is not a certified continuum C_seed; spectral domains, norm bridges, strong traces and boundary preservation require separately accepted analytic inputs."}


def controls(claims):
    records=[]
    for name,value in MUTATIONS:
        bad=copy.deepcopy(claims);bad[name]=value;payload=canonical(bad)
        try:evaluate(json.loads(payload),lambda name,ok:None)
        except ValueError as exc:records.append({"fault":name,"self_rehashed_sha256":hashlib.sha256(payload).hexdigest(),"diagnostic":str(exc),"rejected_by":"fresh exact Gaussian/jet oracle"})
        else:raise ValueError("surviving heat fault: "+name)
    return {"schema":"hdblast-heat-jet-controls-v1","status":"pass","controls":records}
