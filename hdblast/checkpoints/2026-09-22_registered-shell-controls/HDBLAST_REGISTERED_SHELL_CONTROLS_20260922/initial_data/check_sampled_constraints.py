"""Differentiate saved profiles independently of their defining ODE RHS."""
import json
from pathlib import Path
import numpy as np

p=Path(__file__).resolve().parent
q=np.load(p/'BALANCED_CONSTRAINT_SEED_PROFILES.npz');rows=[]
for key in q:
    y,r,ry,psi,sy,z,I=q[key];h=y[1]-y[0];sl=slice(2,-2)
    def d1(v):return (v[:-4]-8*v[1:-3]+8*v[3:-1]-v[4:])/(12*h)
    w=5/3-psi**2+psi**3/3;wp=psi*(psi-2);U=.5*wp*wp-2*w*w/3
    H=6-6*ry[sl]**2-6*r[sl]*d1(ry)-r[sl]**2*sy[sl]**2-2*r[sl]**2*U[sl]
    scale=1+6*ry[sl]**2+(r[sl]*sy[sl])**2;mask=y[sl]>.05
    row={'epsilon':float(key),
      'finite_difference_Hamiltonian_normalized_max':float(max(abs(H[mask])/scale[mask])),
      'finite_difference_rho_first_derivative_error':float(max(abs(d1(r)[mask]-ry[sl][mask]))),
      'finite_difference_phi_first_derivative_error':float(max(abs(d1(psi)[mask]-sy[sl][mask])))}
    if row['finite_difference_Hamiltonian_normalized_max']>1e-8:raise RuntimeError('Hamiltonian profile check failed')
    if row['finite_difference_rho_first_derivative_error']>1e-7:raise RuntimeError('rho derivative profile check failed')
    if row['finite_difference_phi_first_derivative_error']>1e-8:raise RuntimeError('scalar derivative profile check failed')
    rows.append(row)
result={'status':'PASS_15_SAMPLED_CONTROLS','method':'independent fourth-order first differences of saved dense profiles, excluding y<.05 and endpoint stencils','rows':rows}
(p/'SAMPLED_CONSTRAINT_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
