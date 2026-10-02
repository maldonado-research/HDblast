#!/usr/bin/env python3
"""Exact baseline rationalization audit: symbolic variables, zero physical samples."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import sympy as s
from raw_metric_audit import BASELINE_COEFFICIENTS

ROOT=Path(__file__).resolve().parents[1]
path=ROOT/'theory/metric_wkb.py'
name='metric_raw_baseline_symbolic_wkb'
spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
k,a,w=s.symbols('k a w',positive=True)
base=tuple(s.factorial(n)*a**(n+1) for n in range(5))
out=module.generic_subtractions(k*k,base,(s.S.Zero,)*5,r=2)
def exact(expr): return s.nsimplify(expr,rational=True)
ratio=k/w;difference=2*a*a/(w+k)
poly={name:sum(s.Integer(c)*ratio**n for n,c in enumerate(coefficients)) for name,coefficients in BASELINE_COEFFICIENTS.items()}
stable={'Q':difference**3*poly['Q']/(32*k*w**3),'R':difference**4*poly['R']/(1024*k*w**2),'P':-difference**4*poly['P']/(3072*k*w**2)}
bare={'Q':1/(2*k),'R':k/2+3*a*a/(4*k),'P':k/6-a*a/(4*k)}
checks=[];mutations=[]
for name,subname in [('Q','S'),('R','R'),('P','P')]:
    subtraction=exact(out[subname]['value']).subs(s.sqrt(k*k+2*a*a),w)
    residual=s.factor(s.together(bare[name]-subtraction-stable[name]).subs(a*a,(w*w-k*k)/2))
    if residual!=0: raise RuntimeError(name+' exact baseline mismatch: '+str(residual))
    checks.append(name+' exact generic W0/W2/W4 baseline minus rationalized polynomial')
    for change in [stable[name]+difference**4/(k*w*w),2*stable[name]]:
        residual=s.factor(s.together(bare[name]-subtraction-change).subs(a*a,(w*w-k*k)/2))
        if residual==0: raise RuntimeError(name+' baseline mutation survived')
        mutations.append(name+' coefficient/normalization mutation')
report={'status':'PASS','identities':checks,'mutations':mutations,'physical_evaluations':0,'python_optimization':sys.flags.optimize,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'generic_wkb_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'raw_audit_sha256':hashlib.sha256((Path(__file__).parent/'raw_metric_audit.py').read_bytes()).hexdigest(),'scope':'Exact symbolic polynomial identities and six mutation controls only; no modes, sources, response integrals or baseline quadrature evaluated.'}
print(json.dumps(report,indent=2))
