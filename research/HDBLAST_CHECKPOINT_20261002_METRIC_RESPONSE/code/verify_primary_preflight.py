#!/usr/bin/env python3
"""Pure symbolic source bridge and import/authorization guards; no sampling."""
from __future__ import annotations
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as S
from source_jet import POLYNOMIALS


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--algebra',type=Path,required=True)
    parser.add_argument('--independent-inventory',type=Path)
    parser.add_argument('--experiment',type=Path)
    args=parser.parse_args()
    if args.output.exists(): raise RuntimeError('Refusing to overwrite preflight evidence')
    checks=[]
    def eq(name,lhs,rhs):
        if S.simplify(lhs-rhs)!=0: raise RuntimeError('Exact mismatch '+name)
        checks.append({'name':name,'passed':True})
    u=S.symbols('u',real=True)
    d=1-u*u
    B=S.exp(-u*u/d)
    polys=[sum(S.Integer(c)*u**(len(coeffs)-j-1) for j,c in enumerate(coeffs)) for coeffs in POLYNOMIALS]
    for n,p in enumerate(polys):
        eq('exact_bump_derivative_'+str(n),B*p/d**(2*n),S.diff(B,u,n))
        eq('exact_signed_derivative_'+str(n),u*B*p/d**(2*n)+(n*B*polys[n-1]/d**(2*n-2) if n else 0),S.diff(u*B,u,n))
    L=S.symbols('L',positive=True)
    hs=S.symbols('h0:7',real=True)
    D=lambda f:L*L*S.diff(f,L)+sum(S.diff(f,hs[j])*hs[j+1] for j in range(6))
    g=4*L*L*hs[0]-2*L*hs[1]-hs[2]
    for n in range(4):
        candidate=4*sum(S.binomial(n,j)*S.factorial(n-j+1)*L**(n-j+2)*hs[j] for j in range(n+1))-2*sum(S.binomial(n,j)*S.factorial(n-j)*L**(n-j+1)*hs[j+1] for j in range(n+1))-hs[n+2]
        eq('exact_canonical_forcing_derivative_'+str(n),candidate,g)
        g=D(g)
    algebra=json.loads(args.algebra.read_text())
    generated=Path(__file__).with_name('metric_contact_coefficients.py')
    tree=ast.parse(generated.read_text())
    function=next(x for x in tree.body if isinstance(x,ast.FunctionDef))
    returned=next(x for x in function.body if isinstance(x,ast.Return)).value
    for key,value in zip(returned.keys,returned.values):
        name=key.value
        expr=ast.unparse(value).replace('math.pi','pi')
        eq('generated_runtime_expression_'+name,S.sympify(expr),S.sympify(algebra['closed_expressions'][name]))
    if args.independent_inventory:
        independent=json.loads(args.independent_inventory.read_text())['contacts']
        for name,othername in (('q','Q'),('rho','rho'),('p','p')):
            ours=algebra['coefficient_inventory'][name]
            theirs=independent[othername]['integrand_coefficients_of_omega_minus_n']
            if set(ours)!=set(theirs): raise RuntimeError('Independent inventory powers differ')
            for power,x in ours.items():eq('independent_contact_'+name+'_'+power,S.sympify(x),S.sympify(theirs[power]).subs(S.Symbol('h'),S.Symbol('h0')))
    path=Path(__file__).with_name('metric_primary.py')
    # Any source evaluation during import would raise, before the real producer
    # is created. No synthetic physical source is passed to the source function.
    source=sys.modules['source_jet']
    original=source.source_jet_over_epsilon
    def forbidden(*unused):raise RuntimeError('Physical source sampled during import')
    source.source_jet_over_epsilon=forbidden
    try:
        spec=importlib.util.spec_from_file_location('metric_preflight_producer',path)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        if module.AUTHORIZED:raise RuntimeError('Import authorized physical evaluation')
        checks.append({'name':'import_has_no_physical_evaluation','passed':True})
        if args.experiment:
            module.primary_settings(json.loads(args.experiment.read_text()))
            checks.append({'name':'prospective_experiment_source_contract','passed':True})
        for name,call in [('source',lambda:module.metric_jet('positive_B',-4)),('quadrature',lambda:module.integrate(forbidden,-5,-3,{})),('response',lambda:module.direct_response(-4,[0]*6,[0]*4,[0]*3,[0]*3)),('run',lambda:module.run({},forbidden))]:
            try:call()
            except RuntimeError as exc:
                if 'verified prospective registration' not in str(exc):raise
            else:raise RuntimeError('Physical guard failed '+name)
            checks.append({'name':'registration_guard_'+name,'passed':True})
    finally:source.source_jet_over_epsilon=original
    result={'passed':True,'check_count':len(checks),'checks':checks,'scope':'Pure symbols, AST checks and blocked calls only; no physical source, contact, integral or mode evaluated','python_optimization':sys.flags.optimize,'sympy':S.__version__,'provenance':{'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'algebra_sha256':hashlib.sha256(args.algebra.read_bytes()).hexdigest()}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('Primary source/import preflight passed:',len(checks),'checks')


if __name__=='__main__':main()
