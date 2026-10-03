"""Pure symbolic regrouping: linear subtraction inventory to discrete weighted moments."""
import json,hashlib
from pathlib import Path
import sympy as S
from sympy.printing.pycode import PythonCodePrinter
here=Path(__file__).resolve().parent
L,W=S.symbols('L W',positive=True)
H=S.symbols('h0:6')
inventory=json.loads((here/'SYMBOLIC_SPECIALIZATION.json').read_text())
expr=[S.sympify(e,locals={'L':L,'W':W,**{str(x):x for x in H}}) for e in inventory['expressions']]
mu={n:S.Symbol('m'+str(n)) for n in (1,3,5,7,9,11,13,15)}
reduced=[]
for e in expr:
 total=0
 for term in S.Add.make_args(S.expand(e)):
  n=-term.as_powers_dict().get(W,0)
  if n not in mu:raise ValueError('Unexpected inverse omega grade: '+str(n))
  total+=S.cancel(term*W**n)*mu[n]
 reduced.append(S.expand(total))
checks=[S.expand(a-b.xreplace({m:W**(-n) for n,m in mu.items()})) for a,b in zip(expr,reduced)]
if any(c!=0 for c in checks):raise RuntimeError('Exact moment regrouping failed')
rep,red=S.cse(reduced,symbols=S.numbered_symbols('b'),order='canonical')
class Printer(PythonCodePrinter):
 def _print_Rational(self,e):return f'(LD({e.p})/LD({e.q}))'
 def _print_Pow(self,e):
  if e.base==L and e.exp>1 and e.exp.is_Integer:return 'l'+str(e.exp)
  return super()._print_Pow(e)
p=Printer()
lines=['"""Exact moment regrouping of independent Pair inventory; no physical evaluations."""','import numpy as np','LD=np.longdouble','',
       'def weighted_subtractions(L,moments,jet):','    h0,h1,h2,h3,h4,h5=jet','    m1,m3,m5,m7,m9,m11,m13,m15=moments','    l2=L*L']
lines+=['    l'+str(n)+'=l'+str(n-1)+'*L' for n in range(3,17)]
lines+=['    '+str(v)+'='+p.doprint(e) for v,e in rep]
lines+=['    return '+','.join(p.doprint(e) for e in red),'']
text='\n'.join(lines)
(here/'moment_contacts.py').write_text(text)
(here/'MOMENT_REGROUPING_PROOF.json').write_text(json.dumps({'status':'exact_symbolic_PASS','physical_evaluations':0,'expressions':[str(x) for x in reduced],
 'checks':['density_exact_zero','pressure_exact_zero'],'cse_terms':len(rep),'generated_sha256':hashlib.sha256(text.encode()).hexdigest()},indent=2)+'\n')
print(json.dumps({'status':'exact_symbolic_PASS','physical_evaluations':0,'cse_terms':len(rep)}))
