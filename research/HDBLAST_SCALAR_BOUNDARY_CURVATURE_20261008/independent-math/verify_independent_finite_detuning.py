#!/usr/bin/env python3
"""Exact algebra checks for independent weak-coupling theorem; no numerical BVP proof."""
import hashlib
import json
from pathlib import Path
import sympy as s
k,w,d,f0,g1,D,A,c,b,v,g2=s.symbols('k w delta f0 g1 D A c b v g2', real=True)
eta=s.symbols('eta', real=True)
W=3*k+w*eta**2/2+v*eta**3/6
f=f0+c*(g1*eta+g2*eta**2/2)
U=s.diff(W,eta)**2/2-s.Rational(2,3)*W**2
sigma=2*W+d*f
exact=s.expand(sigma**2/36-s.diff(sigma,eta)**2/48+U/6)
reduced=d*(W*f/9-s.diff(W,eta)*s.diff(f,eta)/12)+d**2*(f**2/36-s.diff(f,eta)**2/48)
checks={}
def check(name,expr):
    value=s.factor(s.cancel(expr))
    if value!=0:
        raise RuntimeError(f'{name}: {value}')
    checks[name]=True
check('boundary_identity',exact-reduced)
check('scalar_mass',s.diff(U,eta,2).subs(eta,0)-w*(w-4*k))
check('base_curvature',exact.subs({eta:0,c:0})-((k+d*f0/6)**2-k**2))
expanded=s.series(exact.subs(eta,c*b),c,0,3).removeO().expand()
check('linear_curvature_zero',expanded.coeff(c,1))
raw=expanded.coeff(c,2)
expected=d*g1*b*(4*k-w)/12+d*w*f0*b**2/18+d**2*(f0*g1*b/18-g1**2/48)
check('quadratic_before_matching',raw-expected)
matched=s.factor(raw.subs(b,-d*g1/(2*(D+w))))
Dp=w*(w-4*k)-D**2-4*A*D
candidate=-d**2*g1**2*(4*A*(D+w)-Dp)/(48*(D+w)**2)
check('finite_detuning_coefficient',matched-candidate.subs(A,k+d*f0/6))
# Differentiate E along A'=-(A²-k²), D'=m²-D²-4AD.
E=4*A*(D+w)-Dp
Eprime=s.diff(E,A)*(-(A**2-k**2))+s.diff(E,D)*Dp
# Independently differentiate E=4A(D+w)-D' and eliminate D''.
Eprime_alternative=-4*(A**2-k**2)*(2*D+w)+(2*D+8*A)*Dp
check('barrier_derivative',Eprime-Eprime_alternative)
barrier=Eprime_alternative.subs(Dp,4*A*(D+w))
positive=4*(2*A*D**2+2*A*D*w+6*A**2*D+7*A**2*w+k**2*(2*D+w))
# Explicit replace is used because SymPy substitution in expanded expressions is not reliable.
barrier=-4*(A**2-k**2)*(2*D+w)+(2*D+8*A)*4*A*(D+w)
check('strict_positive_barrier_polynomial',barrier-positive)
lim=candidate.subs({A:k,D:w-4*k})/d**2
check('small_detuning_limit',lim+k*g1**2/(24*(w-2*k)))
check('threshold_massless_coefficient',candidate.subs({w:4*k,D:0})+d**2*g1**2*A/(48*k))
report={'status':'PASS_EXACT_ALGEBRA_ONLY','checks':checks,'count':len(checks),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'sympy_version':s.__version__,'limitations':['Symbolic identities do not validate the nonlinear branch radius or numerical errors.','This file does not run or import the attached numerical benchmark.']}
print(json.dumps(report,indent=2,sort_keys=True))
