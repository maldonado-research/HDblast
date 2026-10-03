#!/usr/bin/env python3
"""Exact derivative/root certificate only: no pulse or interval evaluation."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import traceback
import sympy as S


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module',type=Path,default=Path(__file__).with_name('pulse_derivative_budgets.py'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output.exists(),f'Refusing to overwrite {args.output}')
    report={'status':'RUNNING','scope':'Exact polynomial identities and rational root isolation; no source/response/mode/interval evaluation',
            'python':platform.python_version(),'sympy':S.__version__,'optimization_level':sys.flags.optimize,
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'module_sha256':hashlib.sha256(args.module.read_bytes()).hexdigest(),
            'exact_identities':[],'root_certificates':[],'detected_mutations':[]}
    error=None
    try:
        spec=importlib.util.spec_from_file_location('pulse_budget_exact_inputs',args.module)
        require(spec is not None and spec.loader is not None,'Cannot load adjacent pulse module')
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        require(module.ROOT_BISECTIONS==256 and module.INTERVAL_DECIMAL_PRECISION==60,'Frozen precision changed')
        u,z=S.symbols('u z')
        def polynomial(coefficients,var):
            return S.Poly.from_list(list(coefficients),var).as_expr()
        def equal(name,left,right):
            require(S.expand(left-right)==0,f'{name}: exact polynomial defect')
            report['exact_identities'].append(name)
        for source,initial in [('positive_B',S.Integer(1)),('signed_uB',u)]:
            p=initial
            for n in range(6):
                actual=polynomial(module.DERIVATIVE_POLYNOMIALS[source][n],u)
                equal(f'{source}_derivative_n{n}',actual,p)
                p=S.expand((1-u*u)**2*S.diff(p,u)+(4*n*u*(1-u*u)-2*u)*p)
            for order in (2,3,4):
                derivative=polynomial(module.DERIVATIVE_POLYNOMIALS[source][order+1],u)
                has_zero=module.ZERO_IS_CRITICAL[source][order]
                parity=int(has_zero)
                terms=S.Poly(derivative/u**parity,u).terms()
                require(all(monomial[0]%2==0 for monomial,_ in terms),'Critical-root parity mismatch')
                reduced=S.Poly(sum(coefficient*z**(monomial[0]//2) for monomial,coefficient in terms),z).primitive()[1]
                if reduced.LC()<0:
                    reduced=-reduced
                declared=S.Poly(polynomial(module.CRITICAL_POLYNOMIALS[source][order],z),z)
                equal(f'{source}_critical_polynomial_m{order}',declared.as_expr(),reduced.as_expr())
                require(S.gcd(declared,declared.diff()).degree()==0,'Multiple critical root requires separate treatment')
                require(declared.eval(0)!=0 and declared.eval(1)!=0,'Critical root at support endpoint')
                if has_zero:
                    require(S.diff(derivative,u).subs(u,0)!=0,'Zero critical root is not simple')
                else:
                    require(derivative.subs(u,0)!=0,'Missing zero critical root')
                numerators=module.CRITICAL_BRACKET_NUMERATORS[source][order]
                total=declared.count_roots(0,1)
                require(total==len(numerators),'Critical-root bracket list is incomplete')
                brackets=[]
                previous=S.Integer(0)
                for numerator in numerators:
                    left,right=S.Rational(numerator,256),S.Rational(numerator+1,256)
                    require(0<left<right<1 and previous<left,'Root brackets are not disjoint and ordered')
                    require(declared.eval(left)*declared.eval(right)<0,'Root bracket has no sign change')
                    require(declared.count_roots(left,right)==1,'Root bracket does not isolate one critical root')
                    # All six inherited observations can be ordered before
                    # any interval/source values are evaluated.
                    for eta in module.INHERITED_OBSERVATIONS:
                        endpoint=eta-module.CENTER
                        if -1<endpoint<1 and endpoint!=0:
                            endpoint2=S.Rational(endpoint.numerator**2,endpoint.denominator**2)
                            require(endpoint2<left or endpoint2>right,'Observation/root order needs refinement')
                    brackets.append([str(left),str(right)])
                    previous=right
                report['root_certificates'].append({'source':source,'primitive_order':order,
                    'critical_z_polynomial':str(declared.as_expr()),'roots_in_open_unit_interval':int(total),
                    'simple_zero_u_root':bool(has_zero),'isolating_brackets_z':brackets})
                require(total!=len(numerators)-1,'Omitted-root mutation was not detected')
                report['detected_mutations'].append(f'{source}_omit_one_critical_root_m{order}')
        # Recurrence detects a representative missing product-rule term.
        correct=polynomial(module.DERIVATIVE_POLYNOMIALS['signed_uB'][4],u)
        base=polynomial(module.DERIVATIVE_POLYNOMIALS['signed_uB'][3],u)
        wrong=(1-u*u)**2*S.diff(base,u)-2*u*base
        require(S.expand(correct-wrong)!=0,'Missing denominator-derivative mutation was not detected')
        report['detected_mutations'].append('omit_denominator_derivative_in_recurrence')
        report['status']='PASS'
    except Exception as exc:
        error=exc
        report['status']='FAIL'
        report['failure']=traceback.format_exc()
    report['exact_identity_count']=len(report['exact_identities'])
    report['critical_root_certificate_count']=len(report['root_certificates'])
    report['detected_mutation_count']=len(report['detected_mutations'])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(report,stream,indent=2)
        stream.write('\n')
    print(json.dumps(report,indent=2))
    if error is not None:
        raise RuntimeError('Exact pulse proof failed; original exception preserved in JSON') from error


if __name__=='__main__':
    main()
