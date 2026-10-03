"""Pure symbolic specialization of this route's directional pair inventory."""
import hashlib,json,math
from pathlib import Path
import sympy as S
from sympy.printing.pycode import PythonCodePrinter
import kernel as K

L,W=S.symbols('L W',positive=True)
H=S.symbols('h0:6')
mn=[2*S.factorial(n+1)*L**(n+2) for n in range(5)]
dm=[2*sum(S.binomial(n,j)*mn[j]*H[n-j] for j in range(n+1)) for n in range(5)]
a=[W]; d=[dm[0]/(2*W)]
for n in range(1,5):
 a.append(S.expand((mn[n]-sum(S.binomial(n,j)*a[j]*a[n-j] for j in range(1,n)))/(2*W)))
 d.append(S.expand((dm[n]-sum(S.binomial(n,j)*(d[j]*a[n-j]+a[j]*d[n-j]) for j in range(1,n))-2*d[0]*a[n])/(2*W)))
K.frequency_pairs=lambda k,lp,jet:[K.Pair(v,dv) for v,dv in zip(a,d)]
k=S.sqrt(W**2-2*L**2)
r,p=K.full_subtractions(k,L,H)
expressions=[S.expand(r),S.expand(p)]
rep,red=S.cse(expressions,symbols=S.numbered_symbols('a'),order='canonical')
class Printer(PythonCodePrinter):
 def _print_Rational(self,e):return f'(LD({e.p})/LD({e.q}))'
 def _print_Pow(self,e):
  if e.base==W and e.exp<0 and e.exp.is_Integer:return 'iw'+str(-e.exp)
  if e.base==L and e.exp>1 and e.exp.is_Integer:return 'l'+str(e.exp)
  return super()._print_Pow(e)
p=Printer()
lines=['"""Pure-symbolic specialized independent Pair inventory; no physical evaluations."""','import numpy as np','LD=np.longdouble','',
       'def directional_subtractions(k,L,jet):','    W=np.sqrt(k*k+2*L*L)','    h0,h1,h2,h3,h4,h5=jet','    iw1=1/W','    iw2=iw1*iw1']
lines+=['    iw'+str(n)+'=iw'+str(n-2)+'*iw2' for n in range(3,16,2)]
lines+=['    l2=L*L']+['    l'+str(n)+'=l'+str(n-1)+'*L' for n in range(3,17)]
lines+=['    '+str(v)+'='+p.doprint(e) for v,e in rep]
lines+=['    return '+','.join(p.doprint(e) for e in red),'']
text='\n'.join(lines)
here=Path(__file__).resolve().parent
(here/'specialized_contacts.py').write_text(text)
(here/'SYMBOLIC_SPECIALIZATION.json').write_text(json.dumps({'status':'pure_symbolic_specialization','physical_evaluations':0,
 'expressions': [str(x) for x in expressions],'cse_terms':len(rep),'generated_sha256':hashlib.sha256(text.encode()).hexdigest()},indent=2)+'\n')
print(json.dumps({'status':'pure_symbolic_specialization','physical_evaluations':0,'cse_terms':len(rep)}))
