#!/usr/bin/env python3
"""B3: stiffness of the friction boundary term.  Semi-discrete model problem phi_TT = phi_ZZ on a uniform grid (spacing dz),
shell condition phi_Z(0) = -(Y/2) phi_T imposed through the same degree-5 Hermite ghost as the solver (evolve_a1.build_operators,
order 4).  The most negative eigenvalue lambda scales as Y/dz; explicit RK4 is stable on the negative real axis only for
|lambda| dt <= 2.785.  Writes D2_FRICTION_STIFFNESS.json.  Numerical (eigenvalues of a 800 x 800 matrix)."""
import json, sys, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import evolve_a1 as E
RK4_REAL_LIMIT = 2.785293563405282   # RK4 stability boundary on the negative real axis
out = dict(status='numerical', rk4_real_axis_limit=RK4_REAL_LIMIT, cases=[])
for dz in [1e-3, 5e-4]:
    for Y in [0.5, 1, 2, 3, 4, 5]:
        n = 400; zp = np.full(n, dz); zpp = np.zeros(n)
        D1, D2, g1, g2 = E.build_operators(zp, zpp, 4)
        M = np.zeros((2*n, 2*n)); M[:n, n:] = np.eye(n); M[n:, :n] = D2.toarray(); M[n:, 2*n - 1] += g2*(-Y/2)
        M[0, :] = 0; M[n, :] = 0
        ev = np.linalg.eigvals(M)
        lam = float(np.abs(ev).max())
        for cfl in [0.5, 0.25]:
            out['cases'].append(dict(dz=dz, Y=Y, cfl=cfl, max_abs_lambda_times_dz=lam*dz, max_real=float(ev.real.max()*dz),
                                     abs_lambda_dt=lam*cfl*dz, rk4_stable=bool(lam*cfl*dz <= RK4_REAL_LIMIT)))
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'D2_FRICTION_STIFFNESS.json').write_text(json.dumps(out, indent=1) + '\n')
for c in out['cases']: print(c)
