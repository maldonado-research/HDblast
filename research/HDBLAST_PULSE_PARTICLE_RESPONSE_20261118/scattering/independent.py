"""Independent DOP853 control for the convex-combination mass and outgoing coefficient."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
root=Path(__file__).resolve().parent
d=json.loads((root/'PULSES.json').read_text())
checks={}
for name in ('g8_full_positive','g8_full_negative','g8_half_positive'):
    c=np.array(d['protocols'][name],float)
    knots=np.array(d['centers_in_h'])*.0025
    expected=dict(json.loads((root/'SCATTER.json').read_text())['protocols'][name]['spectrum'])
    vals=[]
    for k in (5.,10.,18.):
        w=np.hypot(k,100*(.5-float(d['phi0'])))
        def rhs(t,y):
            z=np.abs((t-knots)/.0025)
            B=(np.maximum(2-z,0)**3-4*np.maximum(1-z,0)**3)/6
            return [y[1],-(w*w+c@B)*y[0]]
        sol=solve_ivp(rhs,(.1,.15),[1.+0j,-1j*w],rtol=1e-12,atol=1e-14,method='DOP853',max_step=.00006)
        if not sol.success:raise ValueError(sol.message)
        q,p=sol.y[:,-1]
        n=abs(q-1j*p/w)**2/4
        if abs(n-expected[k])>1.1e-11:raise ValueError((n,k,expected[k]))
        vals.append((k,n))
    checks[name]=vals
print(json.dumps(checks,indent=2))
