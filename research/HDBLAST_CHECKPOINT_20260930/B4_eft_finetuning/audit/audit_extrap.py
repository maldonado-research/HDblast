#!/usr/bin/env python3
"""Audit: alternative extrapolations of Lambda_res(d) to d* from the converged static-branch points (audit rerun JSON).
Inverse interpolation d(H2) and fits on different subsets/degrees. Output audit_extrap.json"""
import json, numpy as np
S=json.load(open('B4_STATIC_LAMBDA.json')); D=-3.106933495673783
rows=[r for r in S['branch'] if r['ok']]
d=np.array([r['d'] for r in rows]); L=np.array([r['H2_over_H0sq'] for r in rows]); q=np.array([r['q'] for r in rows])
out={'n_ok':len(rows),'q_ok_max':float(q.max())}
for lo in (0.8,0.85,0.9,0.94):
    for deg in (2,3,4):
        m=q>=lo-1e-9
        if m.sum()<=deg+1: continue
        c=np.polyfit(d[m],L[m],deg); ci=np.polyfit(L[m],d[m],deg)
        out['lo%.2f_deg%d'%(lo,deg)]=dict(npts=int(m.sum()),Lambda_at_dstar=float(np.polyval(c,D)),dstar_exact_inverse=float(np.polyval(ci,0.0)),
            Lambda_at_1p05=float(np.polyval(c,1.05*D)))
json.dump(out,open('audit_extrap.json','w'),indent=1); print(json.dumps(out,indent=1))
