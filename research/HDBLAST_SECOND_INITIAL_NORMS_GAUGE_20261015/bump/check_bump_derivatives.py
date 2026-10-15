#!/usr/bin/env python3
"""Exact rational smooth-bump derivative envelopes, with no external inputs.

The all-point proof is in BUMP_DERIVATIVES_AND_CAUCHY_NORMS.md. Supremum
bounds use polynomial positivity and exponential monotonicity, never sampling.
"""
import argparse
from fractions import Fraction as F
import json
from math import comb, factorial
from pathlib import Path
import subprocess
import sys

ELL = F(1, 10**6)
CONTROLS = {
    "drop_v_derivative": "odd-to-even recurrence retains vprime=-2",
    "replace_v_squared_by_one": "odd-to-even recurrence retains v squared equals one minus four over s",
    "wrong_sprime_sign": "even-to-odd recurrence uses s squared times P minus Pprime",
    "maximize_degree_eight_at_four": "degree-eight exponential envelope is maximized at s=8",
    "false_fourth_derivative_bound": "actual center fourth derivative exceeds 256",
    "omit_fourth_physical_ell": "physical fourth derivative carries ell to the minus fourth power",
}


class Failure(Exception):
    pass


def need(ok, label):
    if not ok:
        raise Failure(label)


def add_to(poly, degree, value):
    poly[degree] = poly.get(degree, F(0)) + value
    if not poly[degree]:
        del poly[degree]


def derivative(poly):
    return {k-1:k*v for k,v in poly.items() if k}


def next_polynomial(poly, order):
    difference = dict(poly)
    for k,value in derivative(poly).items():
        add_to(difference,k,-value)
    result = {}
    for k,value in difference.items():
        add_to(result,k+2,value)
        if order % 2:
            add_to(result,k+1,-4*value)
    if order % 2:
        for k,value in poly.items():
            add_to(result,k,-2*value)
    return result


def polynomials(max_order):
    result = [{0:F(1)}]
    for n in range(max_order):
        result.append(next_polynomial(result[-1],n))
    return result


def value(poly, x):
    return sum((coefficient*x**k for k,coefficient in poly.items()),F(0))


def shifted(poly, shift, degree):
    return [sum((coefficient*comb(k,j)*shift**(k-j)
                 for k,coefficient in poly.items() if k>=j),F(0))
            for j in range(degree+1)]


def envelope_multiplier(poly, order):
    # This algorithm works at arbitrary finite order: (u+4)^(2n) has
    # strictly positive coefficients, so a finite rational multiplier exists.
    degree = 2*order
    coefficients = shifted(poly,F(4),degree)
    ratios = [abs(coefficient)/F(comb(degree,j)*4**(degree-j))
              for j,coefficient in enumerate(coefficients)]
    largest = max([F(1)]+ratios)
    multiplier = F(-(-largest.numerator//largest.denominator))
    for sign in (-1,1):
        difference = {k:sign*v for k,v in poly.items()}
        add_to(difference,degree,multiplier)
        need(all(c>=0 for c in shifted(difference,F(4),degree)),
             "nonnegative shifted coefficients certify polynomial envelope")
    return multiplier


def reciprocal(coefficients):
    need(coefficients[0] != 0,"formal reciprocal constant is nonzero")
    result=[1/coefficients[0]]
    for n in range(1,len(coefficients)):
        result.append(-sum((coefficients[k]*result[n-k] for k in range(1,n+1)),F(0))/coefficients[0])
    return result


def exponential_zero_constant(coefficients):
    need(coefficients[0] == 0,"formal exponential constant removed")
    result=[F(1)]
    for n in range(1,len(coefficients)):
        result.append(sum((k*coefficients[k]*result[n-k] for k in range(1,n+1)),F(0))/n)
    return result


def independent_derivative_ratios(r, max_order):
    # Generate q(r+d)/q(r) by directly expanding the quadratic h(r+d).
    # This route does not use the P_n recurrence or any supremum bound.
    h0=(r-1)*(2-r)
    v=3-2*r
    h=[h0,v,F(-1)]+[F(0)]*max(0,max_order-2)
    inverse=reciprocal(h[:max_order+1])
    exponent=[-a for a in inverse]
    exponent[0]+=1/h0
    expansion=exponential_zero_constant(exponent)
    return [factorial(n)*coefficient for n,coefficient in enumerate(expansion)]


def identity_checks(polys):
    witnesses=0
    for r in (F(9,8),F(5,4),F(3,2),F(7,4),F(15,8)):
        ratios=independent_derivative_ratios(r,len(polys)-1)
        s=1/((r-1)*(2-r))
        v=3-2*r
        for n,poly in enumerate(polys):
            expected=value(poly,s)*(v if n%2 else 1)
            need(ratios[n] == expected,"independent direct-composition derivative identity")
            witnesses+=1
    need(value(polys[2],F(4))==-32,"actual center second derivative")
    need(value(polys[4],F(4))==1536,"actual center fourth derivative")
    return witnesses


def exp_negative_interval(q, degree=40):
    # Taylor's theorem gives alternating remainder signs for exp(-q), q>0.
    need(q>0 and degree%2==0,"even exponential upper partial sum")
    upper=sum(((-q)**j/factorial(j) for j in range(degree+1)),F(0))
    lower=upper+(-q)**(degree+1)/factorial(degree+1)
    need(0<lower<upper,"positive rational exponential witness enclosure")
    return lower,upper


def certificate(max_order):
    need(max_order>=4,"certificate includes the required fourth derivative")
    polys=polynomials(max_order)
    need(all(max(p)==2*n for n,p in enumerate(polys)),"recurrence polynomial degree is twice derivative order")
    multipliers=[envelope_multiplier(p,n) for n,p in enumerate(polys)]
    # exp(2)>this positive partial sum, strictly because all remaining terms>0.
    need(sum((F(2)**j/factorial(j) for j in range(5)),F(0))==7,
         "positive exponential series lower bound exp(2)>7")
    bounds=[]
    for n,multiplier in enumerate(multipliers):
        if n==0:
            bound=F(1)
        elif n==1:
            bound=16*multiplier
        else:
            # d/ds of s^(2n)exp(4-s) has sign 2n-s; its maximum
            # is at 2n (at the endpoint when n=2).
            bound=multiplier*F((2*n)**(2*n),7**(n-2))
        bounds.append(bound)
    retained=[F(1),F(16),F(256),F(46656,7),F(16777216,49),F(10000000000,343),F(8916100448256,2401)]
    need(bounds[:min(len(bounds),7)]==retained[:min(len(bounds),7)],"named derivative constants through order six")
    witnesses=identity_checks(polys)
    return {"schema":"hdblast-bump-derivative-envelopes-v1","status":"PASS",
        "optimized":not __debug__,"ell":str(ELL),"maximum_derivative_order":max_order,
        "definition":"Q(r)=exp(4-1/((r-1)(2-r))) for 1<r<2; Q=0 elsewhere",
        "dimensionless_Q_derivative_sup":[str(b) for b in bounds],
        "physical_q_derivative_sup":[str(b/ELL**n) for n,b in enumerate(bounds)],
        "polynomial_coefficients_ascending":[[str(p.get(k,F(0))) for k in range(2*n+1)] for n,p in enumerate(polys)],
        "polynomial_envelope_multipliers":[str(c) for c in multipliers],
        "formal_composition_identity_witnesses":witnesses,
        "named_checks":["exact parity recurrence","independent formal composition derivative identities",
            "all-point shifted polynomial positivity","exponential envelope maxima and exp(2)>7",
            "named rational constants","physical ell scaling"],
        "scope":{"supremum_bounds_do_not_use_sampled_maxima":True,"flat_smooth_zero_extension":True,
            "all_orders_flatness_uses_analytic_recurrence":True,"fourth_order_suffices_for_first_Cauchy_spatial_derivative":True,
            "physical_Cauchy_norm_evaluated":False,"U2_evolution_evaluated":False}}


def inject(name):
    reason=CONTROLS[name]
    polys=polynomials(4)
    if name=="drop_v_derivative":
        wrong=dict(polys[2])
        for k,coefficient in polys[1].items():add_to(wrong,k,2*coefficient)
        need(value(wrong,F(4))==independent_derivative_ratios(F(3,2),2)[2],reason)
    elif name=="replace_v_squared_by_one":
        wrong=dict(polys[2])
        difference=dict(polys[1])
        for k,coefficient in derivative(polys[1]).items():add_to(difference,k,-coefficient)
        for k,coefficient in difference.items():add_to(wrong,k+1,4*coefficient)
        need(value(wrong,F(4))==independent_derivative_ratios(F(3,2),2)[2],reason)
    elif name=="wrong_sprime_sign":
        wrong=dict(polys[3])
        for k,coefficient in derivative(polys[2]).items():add_to(wrong,k+2,2*coefficient)
        r=F(5,4);s=1/((r-1)*(2-r));v=3-2*r
        need(v*value(wrong,s)==independent_derivative_ratios(r,3)[3],reason)
    elif name=="maximize_degree_eight_at_four":
        lower,_=exp_negative_interval(F(4))
        need(F(4)**8>=F(8)**8*lower,reason)
    elif name=="false_fourth_derivative_bound":
        need(F(256)>=abs(independent_derivative_ratios(F(3,2),4)[4]),reason)
    elif name=="omit_fourth_physical_ell":
        wrong=F(16777216,49)/ELL**3
        actual_center=abs(independent_derivative_ratios(F(3,2),4)[4])/ELL**4
        need(wrong>=actual_center,reason)
    raise Failure("deliberate corruption unexpectedly passed")


def controls(output):
    cases=[]
    for optimized in (False,True):
        for name,reason in CONTROLS.items():
            cmd=[sys.executable,"-I"]+(["-O"] if optimized else [])+["-B",str(Path(__file__).resolve()),"--inject",name]
            run=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
            need(run.returncode==1 and run.stdout=="" and run.stderr=="FAIL: "+reason+"\n",
                 "intended explicit negative-control rejection: "+name)
            cases.append({"name":name,"optimized":optimized,"status":"EXPECTED_EXPLICIT_FAILURE"})
    result={"status":"PASS","control_runs":len(cases),"cases":cases}
    text=json.dumps(result,indent=2)+"\n"
    if output:output.write_text(text)
    print(json.dumps({"status":"PASS","control_runs":len(cases)}))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--order',type=int,default=6)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--inject',choices=list(CONTROLS))
    ap.add_argument('--run-controls',action='store_true')
    args=ap.parse_args()
    try:
        if args.inject:inject(args.inject)
        if args.run_controls:
            controls(args.output)
            return 0
        result=certificate(args.order)
        text=json.dumps(result,indent=2)+"\n"
        if args.output:args.output.write_text(text)
        else:print(text,end='')
    except (Failure,OSError,ValueError) as error:
        print("FAIL: "+str(error),file=sys.stderr)
        return 1
    return 0


if __name__=='__main__':sys.exit(main())
