#!/usr/bin/env python3
"""Exact AST/symbolic review of high-precision code, without source sampling."""
import argparse,ast,hashlib,json,math
from pathlib import Path
from types import SimpleNamespace
import sympy as s
p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
here=Path(__file__).resolve().parent;candidate=here/'reference_inputs/continuum_metric_mp.py'
tree=ast.parse(candidate.read_text());functions={node.name:node for node in tree.body if isinstance(node,ast.FunctionDef)}
checks=[]
def require(label,condition):
 if not condition:raise RuntimeError(label)
 checks.append(label)
def eq(label,left,right):require(label,s.simplify(left-right)==0)
def expression(node,namespace):return eval(compile(ast.Expression(node),'pure_symbolic_AST','eval'),namespace)
contract_node=next(node.value for node in tree.body if isinstance(node,ast.Assign) and any(isinstance(target,ast.Name) and target.id=='CONTRACT' for target in node.targets))
contract=ast.literal_eval(contract_node)
expected={'method':'mpmath_tanh_sinh_endpoint_subtraction_squared','decimal_precisions':[50,70],'maxdegree':10,'panels':['0','0.5','1'],'reported_precision':70,'error_policy':'max_precision_gap_quad_estimates_and_float_conversion','source_arithmetic':'analytic_integer_polynomials_mpmath_only'}
require('exact frozen high-precision algorithm contract',contract==expected)
for name in ('source_jet_mp','log_history_mp','one_precision','continuum_q_jet_mp'):
 node=functions[name];first=node.body[0]
 require('physical guard before '+name,isinstance(first,ast.Expr) and isinstance(first.value,ast.Call) and isinstance(first.value.func,ast.Name) and first.value.func.id=='guard')
# Execute only the pure integer forcing recurrence with formal SymPy inputs.
namespace={'math':math};pure=ast.Module(body=[functions['canonical_forcing_jet_mp']],type_ignores=[])
exec(compile(pure,'pure_symbolic_forcing','exec'),namespace)
L=s.symbols('L',positive=True);h=s.symbols('h0:6',real=True)
actual=namespace['canonical_forcing_jet_mp'](h,-1/L)
g=4*L*L*h[0]-2*L*h[1]-h[2]
D=lambda expr:s.expand(L*L*s.diff(expr,L)+sum(h[n+1]*s.diff(expr,h[n]) for n in range(5)))
for n in range(4):
 if n:g=D(g)
 eq('actual high-precision recurrence symbolic order '+str(n),actual[n],g)
# Compare the actual integrand return expressions with formal symbolic data.
log_function=functions['log_history_mp']
integrands=[node for node in ast.walk(log_function) if isinstance(node,ast.FunctionDef) and node.name=='integrand']
require('exactly two registered continuum branch integrands',len(integrands)==2)
ell,z,logell,difference=s.symbols('ell z logell difference',positive=True)
eta,t,F=s.symbols('eta t F',real=True)
mp=SimpleNamespace(log=s.log,sqrt=s.sqrt,euler=s.Symbol('EulerGamma'),mpf=s.sympify)
first_return=next(node.value for node in reversed(integrands[0].body) if isinstance(node,ast.Return))
eq('actual squared endpoint integrand algebra',expression(first_return,{'ell':ell,'z':z,'difference':difference,'logell':logell,'mp':mp}),2*ell*z*difference*(logell+2*s.log(z)))
second_return=next(node.value for node in reversed(integrands[1].body) if isinstance(node,ast.Return))
eq('actual post-source integrand algebra',expression(second_return,{'z':z,'eta':eta,'t':t,'source':None,'guard':None,'order':1,'forcing_jet_mp':lambda *x:[F]*4,'mp':mp}),2*F*s.log(eta-t))
constant_node=next(node.value for node in log_function.body if isinstance(node,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='constant' for x in node.targets))
eq('actual inherited scale and Euler contact',expression(constant_node,{'mp':mp,'eta':eta}),s.log(s.sqrt(2)/(-eta))+mp.euler+1)
context=next(node for node in functions['one_precision'].body if isinstance(node,ast.With))
qjet_node=next(node.value for node in context.body if isinstance(node,ast.Assign) and any(isinstance(x,ast.Name) and x.id=='qjet' for x in node.targets))
F1,F2,F3,g0,g1,pref=s.symbols('F1 F2 F3 g0 g1 pref',real=True)
qjet=expression(qjet_node,{'pref':pref,'L':L,'memories':[(F1,0),(F2,0),(F3,0)],'endpoint':[g0,g1,0,0]})
for n,target in enumerate((-pref*F1,-pref*(F2+L*g0),-pref*(F3+2*L*g1+L*L*g0))):eq('actual q derivative contact '+str(n),qjet[n],target)
# Configuration and structural checks do not invoke any physical entry point.
text=candidate.read_text()
require('high precision before decimal serialization',"mp.nstr(x,dps)" in text)
require('higher context for precision comparison',"with mp.workdps(config['reported_precision'])" in text)
require('error scope explicitly empirical',"not interval error certificates" in text)
require('no double source evaluator import','source_jet_over_epsilon' not in text)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
receipt={'status':'PASS','scope':'Static AST and exact symbolic source review; no B/uB values, physical responses or numerical quadratures','identity_count':len(checks),'identities':checks,'candidate_sha256':sha(candidate),'source_sha256':sha(Path(__file__)),'physical_evaluations':0}
Path(args.output).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ('status','identity_count','candidate_sha256','physical_evaluations')}))
