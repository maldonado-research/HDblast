"""Audit V3b: continue the c<0 (phi_h < -1) initial-shell branch found in V3 from its last converged row toward
c* = -0.99307 (delta = 0.1). Same solver as V3 (solve_signed, sgn = -1). Time box 600 s.
Output: V3B_CONTINUE_TO_CSTAR.json"""
import json, math, time, numpy as np
from pathlib import Path
import v3_c_family_through_zero as V
import rolloff5d_v1 as R_
HERE = Path(__file__).resolve().parent
prev = json.load(open(HERE/'V3_C_FAMILY_THROUGH_ZERO.json'))['negative_side_rows'][-1]
td = 0.1; x = np.array([math.log10(-1 - prev['phi_h']), prev['y_b']]); c = prev['c']; step = 0.05
rows = []; t0 = time.time(); stop = None
while c > -0.99307 + 1e-12:
    if time.time() - t0 > 600: stop = 'time box'; break
    cn = max(c - step, -0.99307)
    try:
        xn, rn = V.solve_signed(R_.Tension(td, cn, 0), x, -1)
        if not np.all(np.isfinite(rn)) or np.max(np.abs(rn)) > 1e-9: raise RuntimeError('residual %g' % np.max(np.abs(rn)))
        x, c = xn, cn; r = V.row(td, c, xn, rn, -1); rows.append(r)
        print('c=%.5f phi_h=%.6e phi_b=%.5f H0^2/d=%.5f fH=%.4f res=%.1e [%.0fs]' % (c, r['phi_h'], r['phi_b'], r['H0sq_over_delta'], r['fH_series'], r['max_junction_residual'], time.time()-t0), flush=True)
    except Exception as e:
        step /= 2
        if step < 1e-4: stop = 'step underflow at c=%g: %s' % (cn, str(e)[:80]); break
(HERE/'V3B_CONTINUE_TO_CSTAR.json').write_text(json.dumps(dict(start=prev, rows=rows, stop=stop,
    reached_c_star=bool(rows and rows[-1]['c'] <= -0.99307 + 1e-12), runtime_s=time.time()-t0), indent=1) + '\n')
