#!/usr/bin/env python3
"""Exact endpoint algebra, arithmetic/import guards and synthetic receipts."""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sympy as S
import mpmath as mp
import continuum_metric_mp as hp


def audit():
    checks=[]
    def eq(name,lhs,rhs):
        if S.simplify(lhs-rhs)!=0:raise RuntimeError('High-precision exact proof mismatch '+name)
        checks.append({'name':name,'passed':True})
    ell,y,z=S.symbols('ell y z',positive=True)
    f,f0,f1,f2=S.symbols('f f0 f1 f2',real=True)
    eq('endpoint_subtraction_decomposition',ell*f*(S.log(ell)+S.log(y)),ell*(f-f0)*(S.log(ell)+S.log(y))+ell*f0*(S.log(ell)+S.log(y)))
    eq('endpoint_constant_integral',S.integrate(S.log(ell)+S.log(y),(y,0,1)),S.log(ell)-1)
    eq('squared_substitution_jacobian',S.diff(z*z,z),2*z)
    eq('squared_substitution_logarithm',S.expand_log(S.log(z*z),force=True),2*S.log(z))
    formal=f0-f1*ell*y+f2*ell**2*y*y
    original=ell*S.integrate(formal*(S.log(ell)+S.log(y)),(y,0,1))
    transformed=ell*f0*(S.log(ell)-1)+2*ell*S.integrate(z*(formal.subs(y,z*z)-f0)*(S.log(ell)+2*S.log(z)),(z,0,1))
    eq('formal_polynomial_transform_identity',original,transformed)
    integrand=2*ell*z*(formal.subs(y,z*z)-f0)*(S.log(ell)+2*S.log(z))
    for n in range(3):eq('continuous_endpoint_derivative_'+str(n),S.limit(S.diff(integrand,z,n),z,0,dir='+'),0)
    L=S.symbols('L',positive=True);h=S.symbols('h0:7',real=True)
    D=lambda x:L*L*S.diff(x,L)+sum(S.diff(x,h[j])*h[j+1] for j in range(6))
    implemented=hp.canonical_forcing_jet_mp(list(h[:6]),-1/L)
    g=4*L*L*h[0]-2*L*h[1]-h[2]
    for n in range(4):eq('actual_implemented_exact_forcing_'+str(n),implemented[n],g);g=D(g)
    tree=ast.parse(Path(hp.__file__).read_text())
    functions={node.name:node for node in tree.body if isinstance(node,ast.FunctionDef)}
    for name in ('source_jet_mp','forcing_jet_mp','log_history_mp','one_precision'):
        node=functions[name]
        if any(isinstance(x,ast.Call) and isinstance(x.func,ast.Name) and x.func.id=='float' for x in ast.walk(node)):raise RuntimeError('Premature binary64 conversion '+name)
        if any(isinstance(x,ast.Name) and x.id=='source_jet_over_epsilon' for x in ast.walk(node)):raise RuntimeError('Binary64 source wrapper '+name)
        checks.append({'name':'genuine_mp_arithmetic_'+name,'passed':True})
    if not any(isinstance(x,ast.With) and any(isinstance(item.context_expr,ast.Call) and isinstance(item.context_expr.func,ast.Attribute) and item.context_expr.func.attr=='workdps' for item in x.items) for x in ast.walk(functions['one_precision'])):raise RuntimeError('Missing independent precision context')
    checks.append({'name':'independent_full_precision_contexts','passed':True})
    # Fabricated error vectors, not physical source or quadrature values.
    with mp.workdps(70):
        fabricated=[mp.mpf(0),mp.mpf('7e-13'),mp.mpf(1)/3,mp.mpf('1e-60'),mp.mpf('1e-300')]
        raw,normalized=hp.serialize_error_allowances(fabricated)
        for n,(allowance,r,e) in enumerate(zip(fabricated,raw,normalized)):
            if e!=r/(8*hp.math.pi**2) or mp.mpf(e)<allowance:raise RuntimeError('Synthetic upward serialization mismatch')
            checks.append({'name':'synthetic_error_serialization_'+str(n),'passed':True})
    return {'passed':True,'check_count':len(checks),'checks':checks,'physical_evaluations':0,'scope':'Exact symbols, AST guards and fabricated error vectors only; no profile values, physical modes, memory or response quadratures','helper_sha256':hashlib.sha256(Path(hp.__file__).read_bytes()).hexdigest()}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Refusing to overwrite exact algorithm proof')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result=audit();args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print('High-precision pure algorithm audit passed:',result['check_count'],'checks')

if __name__=='__main__':main()
