#!/usr/bin/env python3
"""Pure rational global-envelope certificate; no sampled pulse evaluation.

Bounds extrema algebraically using rational root intervals and exponential
Taylor inequalities. No floating arithmetic, quadrature, mode/response solve,
mpmath interval evaluation, or call to the pulse evaluator is performed.
"""
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

ROOT_STEPS=40
SQRT_STEPS=40
EXP_ORDER=80


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def add(a,b):
    return a[0]+b[0],a[1]+b[1]


def multiply(a,b):
    values=(a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1])
    return min(values),max(values)


def negate(a):
    return -a[1],-a[0]


def absolute(a):
    return (F(0) if a[0]<=0<=a[1] else min(abs(a[0]),abs(a[1])),max(abs(a[0]),abs(a[1])))


def poly(coefficients,x):
    out=(F(0),F(0))
    for coefficient in coefficients:
        out=add(multiply(out,x),(F(coefficient),F(coefficient)))
    return out


def poly_exact(coefficients,x):
    return poly(coefficients,(x,x))[0]


def root_enclosure(coefficients,numerator):
    lo,hi=F(numerator,256),F(numerator+1,256)
    fl=poly_exact(coefficients,lo)
    require(fl*poly_exact(coefficients,hi)<0,'Root bracket lacks sign change')
    for _ in range(ROOT_STEPS):
        mid=(lo+hi)/2
        fm=poly_exact(coefficients,mid)
        if fm==0:
            return mid,mid
        if fl*fm<0:
            hi=mid
        else:
            lo,fl=mid,fm
    return lo,hi


def square_root_enclosure(value):
    lo,hi=F(0),F(1)
    require(0<=value<=1,'Root argument outside unit interval')
    for _ in range(SQRT_STEPS):
        mid=(lo+hi)/2
        if mid*mid<=value:
            lo=mid
        else:
            hi=mid
    return lo,hi


def exp_positive_bounds(t):
    require(0<=t<EXP_ORDER+2,'Geometric exponential remainder does not apply')
    term=total=F(1)
    for n in range(1,EXP_ORDER+1):
        term*=t/n
        total+=term
    remainder=term*t/(EXP_ORDER+1)/(1-t/(EXP_ORDER+2))
    return total,total+remainder


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),'Refusing to overwrite output')
    path=Path(__file__).with_name('pulse_derivative_budgets.py')
    spec=importlib.util.spec_from_file_location('pulse_envelope_inputs',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    results={}
    for source in ('positive_B','signed_uB'):
        rows=[]
        for m in (2,3,4):
            positive=[]
            for numerator in module.CRITICAL_BRACKET_NUMERATORS[source][m]:
                zl,zr=root_enclosure(module.CRITICAL_POLYNOMIALS[source][m],numerator)
                ul=square_root_enclosure(zl)[0]
                ur=square_root_enclosure(zr)[1]
                numerator_interval=poly(module.DERIVATIVE_POLYNOMIALS[source][m],(ul,ur))
                tl,tr=zl/(1-zl),zr/(1-zr)
                lower_exp_at_left=exp_positive_bounds(tl)[0]
                upper_exp_at_right=exp_positive_bounds(tr)[1]
                bump=(1/upper_exp_at_right,1/lower_exp_at_left)
                reciprocal_denominator=((1-zl)**(-2*m),(1-zr)**(-2*m))
                value=multiply(multiply(numerator_interval,bump),reciprocal_denominator)
                positive.append(value)
            odd=(m%2==1) if source=='positive_B' else (m%2==0)
            critical=[negate(x) if odd else x for x in reversed(positive)]
            if module.ZERO_IS_CRITICAL[source][m]:
                atzero=F(module.DERIVATIVE_POLYNOMIALS[source][m][-1])
                critical.append((atzero,atzero))
            critical+=positive
            previous=(F(0),F(0))
            variation=(F(0),F(0))
            for value in critical+[(F(0),F(0))]:
                variation=add(variation,absolute(add(value,negate(previous))))
                previous=value
            upper=variation[1]
            ceiling=(upper.numerator+upper.denominator-1)//upper.denominator
            require(upper<=ceiling,'Integer upper envelope failed exact inequality')
            rows.append({'j':m-2,'global_N_upper_integer':ceiling,
                         'exact_rational_upper_le_integer_verified':True,
                         'number_of_interior_critical_points':len(critical)})
        results[source]=rows
    report={'status':'PASS','scope':'Exact rational analytic global envelopes, not pulse/source sampling or response evaluation',
            'optimization_level':sys.flags.optimize,'root_refinement_bisections':ROOT_STEPS,
            'square_root_bisections':SQRT_STEPS,'positive_exponential_Taylor_order':EXP_ORDER,
            'proof':'N_j is nondecreasing and its global maximum is the full total variation of f^(j+2). Algebraic critical intervals and rational exponential inequalities enclose this variation.',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'pulse_module_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'envelopes':results}
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
