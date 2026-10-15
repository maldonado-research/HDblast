#!/usr/bin/env python3
"""Exact planning checks for one S3 quadrupole; outside the frozen replay."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import sympy as sp

x=sp.symbols("x1 x2 x3 x4", real=True)
s=sum(v*v for v in x)
a=x[0]**2+x[1]**2
b=x[0]**2-x[1]**2
p4=b*b
z2=a-s/2
h4=p4-s*a/2+s*s/12

def lap(p):
    return sp.expand(sum(sp.diff(p,v,2) for v in x))

def euler(p):
    return sp.expand(sum(v*sp.diff(p,v) for v in x))

def sphere_lap(p):
    # Unit-S3 restriction of r^2 Delta_R4 - Euler^2 - 2 Euler.
    return sp.expand(s*lap(p)-euler(euler(p))-2*euler(p))

def norm_grad(p):
    # Homogeneous tangential numerator, restricted to r=1 at the end.
    return sp.expand(s*sum(sp.diff(p,v)**2 for v in x)-euler(p)**2)

def average_s3(p):
    # Normalized S3 even moments: product (1/2)_(alpha_i)/(2)_(sum alpha_i).
    total=sp.Rational(0)
    for powers,coefficient in sp.Poly(sp.expand(p),*x).terms():
        if any(k%2 for k in powers):
            continue
        half=[k//2 for k in powers]
        moment=sp.prod(sp.rf(sp.Rational(1,2),k) for k in half)/sp.rf(2,sum(half))
        total+=coefficient*moment
    return sp.simplify(total)

def stringify(expression):
    return sp.sstr(sp.factor(expression))

def checked_group(name,expressions):
    residuals={key:sp.expand(value) for key,value in expressions.items()}
    failures={key:stringify(value) for key,value in residuals.items() if value!=0}
    if failures:
        raise ArithmeticError(f"{name}: nonzero residuals {failures}")
    return {"name":name,"passed":True,"exact_residuals":{key:stringify(value) for key,value in residuals.items()}}

def rejected_control(name,residual):
    residual=sp.expand(residual)
    if residual==0:
        raise ArithmeticError(f"wrong-formula control unexpectedly passed: {name}")
    return {"name":name,"wrong_formula_rejected":True,"exact_nonzero_residual":stringify(residual)}

def build_receipt():
    gradient=norm_grad(b)
    decomposed_gradient=sp.Rational(4,3)*s*s+2*s*z2-4*h4
    groups=[]
    groups.append(checked_group("ambient harmonicity and homogeneity",{
        "Delta B":lap(b),"Delta Z2":lap(z2),"Delta H4":lap(h4),
        "Euler B minus 2 B":euler(b)-2*b,
        "Euler Z2 minus 2 Z2":euler(z2)-2*z2,
        "Euler H4 minus 4 H4":euler(h4)-4*h4,
    }))
    groups.append(checked_group("ambient Fischer decomposition",{
        "P4 minus H4 minus s Z2/2 minus s^2/6":p4-h4-s*z2/2-s*s/6,
        "Delta P4 minus 8 A":lap(p4)-8*a,
        "Delta sA minus 16 A minus 4 s":lap(s*a)-16*a-4*s,
        "Delta s^2 minus 24 s":lap(s*s)-24*s,
    }))
    groups.append(checked_group("unit S3 Laplace eigenvalues",{
        "DeltaS3 B plus 8 B":sphere_lap(b)+8*b,
        "DeltaS3 Z2 plus 8 Z2":sphere_lap(z2)+8*z2,
        "DeltaS3 H4 plus 24 H4":sphere_lap(h4)+24*h4,
    }))
    groups.append(checked_group("direct tangential-gradient projection",{
        "projected norm minus 4 s A plus 4 P4":gradient-4*s*a+4*p4,
        "projected norm minus harmonic decomposition":gradient-decomposed_gradient,
    }))
    product_route=(sphere_lap(p4)+16*p4)/2
    groups.append(checked_group("independent Laplace product rule",{
        "DeltaS3 P4 plus 16 P4 minus 2 gradient":sphere_lap(p4)+16*p4-2*gradient,
        "product-rule gradient minus direct projection":product_route-gradient,
        "product-rule gradient minus harmonic decomposition":product_route-decomposed_gradient,
    }))
    groups.append(checked_group("normalized S3 moments and Dirichlet identity",{
        "mean Y^2 minus 1/6":average_s3(p4)-sp.Rational(1,6),
        "mean Z2":average_s3(z2),"mean H4":average_s3(h4),
        "mean gradient norm minus 4/3":average_s3(gradient)-sp.Rational(4,3),
        "mean gradient norm minus 8 mean Y^2":average_s3(gradient)-8*average_s3(p4),
    }))
    wrong_h4=p4-s*a/2+s*s/24
    eps=sp.Symbol("epsilon")
    # F(u,v)=v+u^2, U0=(0,0), U1=(1,0), expansion eps^2 U2/2.
    correct_second=sp.expand(-2*eps**2/2+eps**2)
    if correct_second!=0:
        raise ArithmeticError("correct Taylor normalization failed")
    controls=[
        rejected_control("H4 with r^4/24 instead of r^4/12",lap(wrong_h4)),
        rejected_control("omit degree-two term in Y^2 decomposition",p4-h4-s*s/6),
        rejected_control("use S2 eigenvalue -6 for Y on S3",sphere_lap(b)+6*b),
        rejected_control("use S2 eigenvalue -20 for H4 on S3",sphere_lap(h4)+20*h4),
        rejected_control("use plus 4 H4 in gradient identity",gradient-(sp.Rational(4,3)*s*s+2*s*z2+4*h4)),
        rejected_control("extra factor one-half in second-order forcing",sp.expand(-eps**2/2+eps**2)),
    ]
    script=Path(__file__).resolve()
    note=script.with_name("NEXT_SECOND_ORDER_TARGET.md")
    return {
        "schema":"hdblast-second-order-harmonic-planning-v1",
        "classification":"public planning identities only",
        "scope":"selected unnormalized Y=n1^2-n2^2 on unit S3; exact angular algebra and derivative-convention controls",
        "runtime":{"python":platform.python_version(),"sympy":sp.__version__,"optimization":sys.flags.optimize},
        "arithmetic":"exact symbolic rational polynomial arithmetic; no floating-point decisions",
        "check_groups":groups,"check_group_count":len(groups),
        "wrong_formula_controls":controls,"wrong_formula_control_count":len(controls),
        "normalization":{"mean_Y_squared":"1/6","mean_gradient_squared":"4/3","sphere_volume":"2*pi^2"},
        "checker_sha256":hashlib.sha256(script.read_bytes()).hexdigest(),
        "planning_note_sha256":hashlib.sha256(note.read_bytes()).hexdigest(),
        "uses_assert_for_scientific_checks":False,
        "frozen_waveform_replay_modified":False,
        "counted_in_frozen_replay_roster":False,
        "quadratic_bulk_or_junction_sources_verified":False,
        "second_order_equations_solved_by_this_check":False,
        "prior_conditional_complete_coefficient_construction_retained":True,
        "prior_conditional_exact_family_hypotheses_retained":True,
        "nonlinear_existence_or_remainder_established":False,
        "novelty_claimed":False,
    }

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).resolve().with_name("SECOND_ORDER_HARMONIC_IDENTITY_RECEIPT.json"))
    args=parser.parse_args()
    receipt=build_receipt()
    args.output.write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(f"Passed {receipt['check_group_count']} exact planning groups; rejected {receipt['wrong_formula_control_count']} wrong-formula controls. Receipt: {args.output}")

if __name__=="__main__":
    main()
