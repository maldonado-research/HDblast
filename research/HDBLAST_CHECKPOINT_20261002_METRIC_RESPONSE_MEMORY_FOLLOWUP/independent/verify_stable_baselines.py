#!/usr/bin/env python3
"""Pure exact proof of the rationalized baseline integrands; zero physical samples."""
from pathlib import Path
import hashlib
import json
import sympy as s
import metric_wkb as wk
import stable_baselines as sb

wk.LD=s.Integer
L,w,k=s.symbols('L w k',positive=True)
D=lambda expr:s.diff(expr,L)*L**2+s.diff(expr,w)*2*L**3/w
omega=[w]
for n in range(4):
    omega.append(s.expand(D(omega[-1])))
out=wk.evaluate_wkb(k*k,L,2*L*L,omega,[s.Integer(0)]*5,
                    [2*L*L,4*L**3,12*L**4],[s.Integer(0)]*3,0,0)
# SymPy integer coercion and a scalar zero preserve the exact runtime Horner
# expression while avoiding any numerical point or physical source.
sb.LD=s.Integer
saved_zeros=sb.np.zeros_like
sb.np.zeros_like=lambda value:s.Integer(0)
stable=sb.baseline_differences(k,w,L)
sb.np.zeros_like=saved_zeros
bare=(1/(2*k),k/2+3*L*L/(4*k),k/6-L*L/(4*k))
checks=[]
for name,barevalue,subname,combined in zip(('Q','rho','p'),bare,('S','R','P'),stable):
    residual=s.factor(s.together(barevalue-out[subname]-combined).subs(L**2,(w*w-k*k)/2))
    if residual!=0:
        raise RuntimeError(name+' baseline identity failed: '+str(residual))
    checks.append(name+' exact bare minus complete subtraction rationalization')
here=Path(__file__).resolve().parent
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
proof={'schema_version':1,'status':'PASS','scope':'pure exact algebra; zero physical evaluations',
       'identities':checks,'baseline_module_sha256':sha(here/'stable_baselines.py'),
       'wkb_module_sha256':sha(here/'metric_wkb.py'),'source_sha256':sha(Path(__file__)),
       'literal_bare_and_subtracted_observation_pieces_retained':True,
       'integrated_quantity':'exact combined directional/baseline operators; no Ward reconstruction'}
(here/'BASELINE_RATIONALIZATION_PROOF.json').write_text(json.dumps(proof,indent=2)+'\n')
print(json.dumps(proof,indent=2))
