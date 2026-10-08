"""Portable exact algebra check of the conditional local response formula.

No archived inputs or numerical output are needed. SymPy is the only dependency.
This verifies algebra under branch assumptions, not existence or external novelty.
"""
import argparse, hashlib, json, sys
from pathlib import Path
parser=argparse.ArgumentParser()
parser.add_argument('--library-path',help='Optional existing directory containing SymPy')
args=parser.parse_args()
if args.library_path: sys.path.insert(0,args.library_path)
import sympy as s
k,w,v,u,f0,f1,f2,f3,e,d,b=s.symbols('k w v u f0 f1 f2 f3 eta delta b')
W=3*k+w*e**2/2+v*e**3/6+u*e**4/24
f=f0+f1*e+f2*e**2/2+f3*e**3/6
Wp=s.diff(W,e); fp=s.diff(f,e)
U=Wp**2/2-2*W**2/3
sigma=2*W+d*f
H=sigma**2/36-s.diff(sigma,e)**2/48+U/6
A=W*f/9-Wp*fp/12
B=f*f/36-fp*fp/48
lam=w-4*k
a=-f1/(2*(lam+w))
coef1=A.subs(e,0)
coef2=s.factor(a*s.diff(A,e).subs(e,0)+B.subs(e,0))
expected=f0*f0/36-k*f1*f1/(24*(w-2*k))
checks={
 'exact_boundary_reduction':s.expand(H-d*A-d*d*B)==0,
 'stationary_bulk':s.diff(U,e).subs(e,0)==0,
 'background_curvature':s.simplify(U.subs(e,0)+6*k*k)==0,
 'bulk_scalar_mass':s.simplify(s.diff(U,e,2).subs(e,0)-w*(w-4*k))==0,
 'growing_indicial_root':s.expand(lam*lam+4*k*lam-w*(w-4*k))==0,
 'leading_scalar_junction':s.simplify((lam+w)*a+f1/2)==0,
 'leading_curvature':s.simplify(coef1-k*f0/3)==0,
 'second_order_curvature_with_higher_derivatives':s.simplify(coef2-expected)==0,
 'higher_derivative_independence':all(s.diff(coef2,q)==0 for q in (v,u,f2,f3)),
}
result=dict(pass_all=all(checks.values()),checks=checks,sympy=s.__version__,
            scalar_shift=str(s.factor(a)),second_order_curvature=str(coef2),
            assumptions=['k>0,w>4k,f0>0,delta approaches zero from positive side',
              'specified regular growing-mode branch exists, with eta_b=O(delta),H2=O(delta)',
              'boundary scalar derivative=lambda*eta_b+O(delta**2); smooth local expansion'],
            limitations=['no global branch existence or remainder proof','no stability or physical detection',
                         'no external priority claim'],
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('portable-algebra-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
if not result['pass_all']:raise SystemExit(1)
