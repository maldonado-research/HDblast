#!/usr/bin/env python3
"""Pure exact review of a high-precision logarithmic-memory repair proposal.

No B/uB values, source samples, physical response, modes or quadratures run.
Only formal symbolic coefficients and synthetic polynomials are integrated.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,platform
from pathlib import Path
import sympy as s
p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
here=Path(__file__).resolve().parent
checks=[];mutations=[]
def eq(label,left,right=0):
    residual=s.simplify(left-right)
    if residual!=0:raise RuntimeError(label+': '+str(residual))
    checks.append(label)
def reject(label,residual):
    if s.simplify(residual)==0:raise RuntimeError('Invisible mutation: '+label)
    mutations.append(label)
ell,y,z=s.symbols('ell y z',positive=True)
eta=s.symbols('eta',real=True)
f0,f1,f2=s.symbols('f0 f1 f2',real=True)
f=s.symbols('f',real=True)
eq('endpoint substitution orientation',s.diff(eta-ell*y,y),-ell)
eq('positive logarithmic scale split',s.expand_log(s.log(ell*y),force=True),s.log(ell)+s.log(y))
eq('exact endpoint subtraction decomposition',f*(s.log(ell)+s.log(y)),(f-f0)*(s.log(ell)+s.log(y))+f0*(s.log(ell)+s.log(y)))
eq('analytic endpoint constant integral',s.integrate(s.log(ell)+s.log(y),(y,0,1)),s.log(ell)-1)
eq('square substitution Jacobian',s.diff(z*z,z),2*z)
eq('square substitution logarithm',s.expand_log(s.log(z*z),force=True),2*s.log(z))
eq('square substitution analytic constant',s.integrate(2*z*(s.log(ell)+2*s.log(z)),(z,0,1)),s.log(ell)-1)
# These are arbitrary formal polynomial coefficients, not any physical source.
synthetic=f0-f1*ell*y+f2*ell**2*y**2
original=ell*s.integrate(synthetic*(s.log(ell)+s.log(y)),(y,0,1))
subtracted=ell*f0*(s.log(ell)-1)+ell*s.integrate((synthetic-f0)*(s.log(ell)+s.log(y)),(y,0,1))
squared=ell*f0*(s.log(ell)-1)+2*ell*s.integrate(z*(synthetic.subs(y,z*z)-f0)*(s.log(ell)+2*s.log(z)),(z,0,1))
eq('synthetic formal polynomial subtraction agreement',subtracted,original)
eq('synthetic formal polynomial square-transform agreement',squared,original)
regularized=2*ell*z*(-f1*ell*z*z+f2*ell**2*z**4)*(s.log(ell)+2*s.log(z))
for n in (0,1,2):eq('continuous square endpoint derivative '+str(n),s.limit(s.diff(regularized,z,n),z,0,dir='+'))
for n in (1,2,3,4):eq('flat endpoint dominates logarithm power '+str(n),s.limit(z**n*s.log(z),z,0,dir='+'))
eq('post-source substitution Jacobian',s.diff(-5+2*z,z),2)
eq('post-source logarithm argument',(eta-(-5+2*z)),eta+5-2*z)
# Independently verify the exact integer bump-derivative polynomial recurrence.
path=here/'reference_inputs/source_jet.py'
spec=importlib.util.spec_from_file_location('reference_source_jet',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
u=s.symbols('u',real=True)
P=[sum(s.Integer(c)*u**(len(coefficients)-j-1) for j,c in enumerate(coefficients)) for coefficients in module.POLYNOMIALS]
eq('bump logarithmic derivative',s.diff(-u*u/(1-u*u),u),-2*u/(1-u*u)**2)
for n in range(5):eq('integer bump derivative polynomial recurrence '+str(n),(1-u*u)**2*s.diff(P[n],u)+4*n*u*(1-u*u)*P[n]-2*u*P[n],P[n+1])
# Exact metric forcing derivatives use h0 through h5 and no finite differences.
L=s.symbols('L',positive=True);h=s.symbols('h0:6',real=True)
def D(expression):return s.expand(L*L*s.diff(expression,L)+sum(h[n+1]*s.diff(expression,h[n]) for n in range(5)))
g=4*L*L*h[0]-2*L*h[1]-h[2]
gn=g
for n in range(4):
    if n:gn=D(gn)
    recurrence=4*sum(s.binomial(n,j)*s.factorial(n-j+1)*L**(n-j+2)*h[j] for j in range(n+1))-2*sum(s.binomial(n,j)*s.factorial(n-j)*L**(n-j+1)*h[j+1] for j in range(n+1))-h[n+2]
    eq('canonical forcing exact derivative recurrence '+str(n),gn,recurrence)
eq('metric forcing first derivative',D(g),8*L**3*h[0]+2*L**2*h[1]-2*L*h[2]-h[3])
eq('metric forcing second derivative',D(D(g)),24*L**4*h[0]+12*L**3*h[1]-2*L*h[3]-h[4])
eq('metric forcing third derivative',D(D(D(g))),96*L**5*h[0]+60*L**4*h[1]+12*L**3*h[2]-2*L**2*h[3]-2*L*h[4]-h[5])
# The differentiable convolution has Fn'=F(n+1)+L*g(n-1).
F1,F2,F3=s.symbols('F1 F2 F3',real=True);g0,g1,g2=s.symbols('g0 g1 g2',real=True)
def DM(expression):return s.expand(s.diff(expression,L)*L*L+s.diff(expression,F1)*(F2+L*g0)+s.diff(expression,F2)*(F3+L*g1)+s.diff(expression,g0)*g1+s.diff(expression,g1)*g2)
eq('subtraction scale logarithmic derivative',D(s.log(s.sqrt(2)*L)),L)
eq('first memory derivative contact',DM(-F1),-F2-L*g0)
eq('second memory derivative contacts',DM(DM(-F1)),-F3-2*L*g1-L*L*g0)
reject('omit endpoint constant analytic term',ell*f0*(s.log(ell)-1))
reject('omit square Jacobian',2*z-1)
reject('lose factor two in square logarithm',2*s.log(z)-s.log(z))
reject('wrong memory second derivative metric drift',2*L*g1+L*L*g0)
reject('truncate h fifth derivative from g third',h[5])
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
receipt={'status':'PASS','scope':'Pure symbolic identities and formal synthetic polynomials; no physical source samples, response values, mode solutions, or physical quadrature','identity_count':len(checks),'identities':checks,'mutation_count':len(mutations),'mutations':mutations,'python':platform.python_version(),'sympy':s.__version__,'source_sha256':sha(Path(__file__)),'input_sha256':{str(q.relative_to(here)):sha(q) for q in [path,here/'reference_inputs/original_metric_primary.py']},'physical_evaluations':0}
Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ('status','identity_count','mutation_count','physical_evaluations')}))
