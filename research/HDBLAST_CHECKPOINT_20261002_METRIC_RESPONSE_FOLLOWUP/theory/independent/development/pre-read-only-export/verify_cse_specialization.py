#!/usr/bin/env python3
"""Independent exact audit of the generated numeric route's CSE evaluator."""
from __future__ import annotations
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sympy as s
import metric_wkb as generic
p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
here=Path(__file__).resolve().parent
candidate=here/'reference_inputs/CSE_metric_wkb.py'
spec=importlib.util.spec_from_file_location('cse_candidate',candidate);cse=importlib.util.module_from_spec(spec);spec.loader.exec_module(cse)
cse.LD=s.Integer;generic.ONE=s.Integer(1);generic.HALF=s.Rational(1,2)
L,w=s.symbols('L w',positive=True);hj=s.symbols('hj0:5',real=True)
ad=tuple(s.factorial(n)*L**(n+1) for n in range(5))
dad=tuple(sum(s.binomial(n,j)*ad[n-j]*hj[j] for j in range(n+1)) for n in range(5))
base=generic.generic_subtractions(w*w-2*L*L,ad,dad,r=2)
ns={'__file__':str(here/'explore_contacts.py')}
with contextlib.redirect_stdout(io.StringIO()):exec(compile((here/'explore_contacts.py').read_text(),ns['__file__'],'exec'),ns)
xx={ns['L']:L,ns['w']:w,ns['M']:2*L*L,**dict(zip([ns[n] for n in ('h','h1','h2','h3','h4')],hj))}
D=ns['D'];M=ns['M'];ww=ns['w'];z=M*ns['h']/ww
omega=[];deltaomega=[]
for n in range(5):
    if n==0:ob=ww;od=z
    else:ob=D(ob);od=D(od)
    omega.append(s.expand(ob.subs(xx)));deltaomega.append(s.expand(od.subs(xx)))
C=[2*L*L,4*L**3,12*L**4]
dC=[hj[2]+2*L*hj[1],hj[3]+2*L*L*hj[1]+2*L*hj[2],hj[4]+4*L**3*hj[1]+4*L*L*hj[2]+2*L*hj[3]]
out=cse.evaluate_wkb(w*w-2*L*L,L,2*L*L,omega,deltaomega,C,dC,hj[1],4*L*L*hj[0])
checks=[]
def eq(label,left,right):
    residual=s.expand(left-right)
    if residual!=0:raise RuntimeError(label+': '+str(s.factor(residual)))
    checks.append(label)
for name in ('R0','R2','R4','P0','P2','P4','S0','S2','R','P','S','W2','W4','J2','J4'):
    eq('CSE baseline '+name,out[name],base[name]['value'])
    eq('CSE directional '+name,out['delta'+name],base[name]['delta'])
for name,label in [('S_prime','first'),('S_second','second')]:
    eq('CSE baseline '+name,out[name],base['S'][label])
    eq('CSE directional '+name,out['delta'+name],base['S']['delta_'+label])
for name,reference in [('W2_prime',ns['Up']),('W2_second',ns['Upp'])]:
    eq('CSE baseline '+name,out[name],s.expand(reference.subs(xx)))
for name,reference in [('deltaW2_prime',ns['u1']),('deltaW2_second',ns['u2'])]:
    eq('CSE directional '+name,out[name],s.expand(reference.subs(xx)))
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
receipt={'status':'PASS','scope':'Exact symbolic comparison; zero physical evaluations','identity_count':len(checks),'identities':checks,'candidate_sha256':sha(candidate),'generic_dual_jet_sha256':sha(here/'metric_wkb.py'),'expanded_chain_rule_sha256':sha(here/'explore_contacts.py'),'source_sha256':sha(Path(__file__)),'sympy':s.__version__,'physical_results_computed':False}
Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ('status','identity_count','candidate_sha256','physical_results_computed')}))
