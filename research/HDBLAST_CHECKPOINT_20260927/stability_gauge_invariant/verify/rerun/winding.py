"""Argument-principle count of complex scalar eigenvalues.

M(mu2) (regular solution normalised X ~ y^s, s = -3/2 + sqrt(9/4 - mu2), principal branch) is analytic
in mu2 off the continuum cut [9/4, inf).  The number of zeros inside a closed contour equals the
winding number of M around 0.  Contours: rectangles Re mu2 in [-60, 2] x Im in [-30, 30], and Re in [-400, 2] x Im in [-200, 200].
Calibration: the original shell must give exactly 1 (its real tachyon at -7.718).
Adaptive refinement: every segment is bisected until |d arg M| < 0.3 rad.
Output: WINDING_RESULTS.json
"""
import sys, json, time, math, hashlib, os
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np
from multiprocessing import Pool
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gi_core import Background

RECTS = [(-60.0, 2.0, -30.0, 30.0), (-400.0, 2.0, -200.0, 200.0)]

def contour_count(task):
    branch, delta, e_h, y_b, RECT = task
    t0 = time.time()
    bg = Background(+1 if branch == 'plus' else -1, delta, e_h=e_h, y_b=y_b, polish=False)
    a, b, c, d = RECT
    corners = [complex(a, c), complex(b, c), complex(b, d), complex(a, d), complex(a, c)]
    cache = {}
    def M(z):
        if z not in cache:
            cache[z] = complex(bg.scalar_M(complex(z), full=True)['M'])
        return cache[z]
    total = 0.0; nseg = 0; maxstep = 0.0
    for k in range(4):
        pts = list(np.linspace(0, 1, 25))
        z0, z1 = corners[k], corners[k + 1]
        stack = [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)][::-1]
        while stack:
            s0, s1 = stack.pop()
            za, zb = z0 + (z1 - z0)*s0, z0 + (z1 - z0)*s1
            dphi = np.angle(M(zb)/M(za))
            if abs(dphi) > 0.3 and (s1 - s0) > 1e-7:
                sm = 0.5*(s0 + s1); stack.append((sm, s1)); stack.append((s0, sm))
            else:
                total += dphi; nseg += 1; maxstep = max(maxstep, abs(dphi))
    w = total/(2*math.pi)
    res = dict(branch=branch, delta=delta, rectangle=list(RECT), winding=w, winding_rounded=int(round(w)), segments=nseg,
               evaluations=len(cache), max_abs_dphase=maxstep, min_abs_M_on_contour=float(min(abs(v) for v in cache.values())),
               runtime_s=time.time() - t0)
    print(json.dumps(res), flush=True)
    return res

if __name__ == '__main__':
    t0 = time.time()
    B = json.loads((HERE/'BACKGROUNDS.json').read_text())
    tasks = [('original', r['delta'], r['gi_core']['e_h'], r['gi_core']['y_b'], R) for R in RECTS for r in B['original'] if r['delta'] == 0.001]
    tasks += [('plus', r['delta'], r['gi_core']['e_h'], r['gi_core']['y_b'], R) for R in RECTS for r in B['plus']]
    with Pool(2) as pool:
        results = pool.map(contour_count, tasks, chunksize=1)
    orig = [r for r in results if r['branch'] == 'original']; plus = [r for r in results if r['branch'] == 'plus']
    out = dict(status='NUMERICAL', contour_rectangles=RECTS, results=results,
               calibration_pass=bool(all(r['winding_rounded'] == 1 and abs(r['winding'] - 1) < 1e-6 for r in orig)),
               plus_all_zero=bool(all(r['winding_rounded'] == 0 and abs(r['winding']) < 1e-6 for r in plus)),
               runtime_s=time.time() - t0, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (HERE/'WINDING_RESULTS.json').write_text(json.dumps(out, indent=2) + '\n')
    print('calibration', out['calibration_pass'], 'plus zero', out['plus_all_zero'], out['runtime_s'])
