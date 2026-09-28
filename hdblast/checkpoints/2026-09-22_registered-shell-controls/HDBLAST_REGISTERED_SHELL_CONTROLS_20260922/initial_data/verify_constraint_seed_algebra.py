"""Exact symbolic checks of the initial-data construction and its obstruction."""
from pathlib import Path
import json
import sympy as S

r,v,p,s,U,Up,eps,b=S.symbols('rho v phi s U Up epsilon bump', real=True)
sig,sig1=S.symbols('sigma sigma1',real=True)
ryy=(1-v*v)/r-r*s*s/6-r*U/3
syy=Up-4*v*s/r+eps*b
D=1-v*v+r*r*(s*s/12-U/6)
Dy=S.diff(D,r)*v+S.diff(D,v)*ryy+S.diff(D,s)*syy+S.diff(D,U)*Up*s
Aacc=r*ryy-3+3*v*v+2*r*r*U/3
Bacc=r*ryy+3-3*v*v+r*r*s*s/2-r*r*U/3
Facc=r*r*syy+4*r*v*s-r*r*Up
H=-2*U*r*r+6-12*v*v+6*v*v-6*r*ryy-r*r*s*s
checks={}
def check(name,expr):
    result=S.simplify(expr)==0;checks[name]=bool(result)
    if not result:raise RuntimeError(name+': '+str(S.simplify(expr)))
check('Hamiltonian identity',H)
check('Mass defect transport',Dy+2*v*D/r-eps*r*r*s*b/6)
check('A acceleration',Aacc+2*D)
check('B acceleration',Bacc-4*D)
check('Scalar acceleration',Facc-eps*r*r*b)
Db=S.symbols('Db',real=True)
# Near the brane b vanishes identically; D_y=-2 v D/r.
firstA=0;firstB=0;firstF=0
gA2=r*sig*4*Db/6;gF2=-r*sig1*4*Db/2
a2z=4*v*Db;b2z=-8*v*Db;f2z=0
check('A second-time corner',(a2z-gA2).subs(sig,6*v/r))
check('B second-time obstruction',(b2z-gA2+12*v*Db).subs(sig,6*v/r))
check('Scalar second-time obstruction',(f2z-gF2+4*r*s*Db).subs(sig1,-2*s))
out={'status':'PASS_EXACT_ALGEBRA_AND_OBSTRUCTION','checks':checks,'count':len(checks),
     'scope':'Identities only; generic positive-bump seeds fail second-time corner compatibility.'}
Path(__file__).with_name('CONSTRAINT_SEED_ALGEBRA.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
