"""Auditor: re-run the audited winding.contour_count for a subset (calibration delta=0.001 and +1 delta=0.001 and 0.1,
small rectangle) single-process, and compare with the audited WINDING_RESULTS.json.  Output: winding_subset.json"""
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
from winding import contour_count, RECTS
B = json.load(open('BACKGROUNDS.json'))
g = {(k, r['delta']): r['gi_core'] for k in ('plus', 'original') for r in B[k]}
tasks = [('original', 0.001), ('plus', 0.001), ('plus', 0.1)]
res = [contour_count((b, d, g[(b, d)]['e_h'], g[(b, d)]['y_b'], RECTS[0])) for b, d in tasks]
aud = json.load(open('../../WINDING_RESULTS.json'))['results']
cmp = []
for r in res:
    a = [x for x in aud if x['branch'] == r['branch'] and x['delta'] == r['delta'] and x['rectangle'] == r['rectangle']][0]
    cmp.append(dict(branch=r['branch'], delta=r['delta'], rerun=r['winding'], audited=a['winding'], min_abs_M=r['min_abs_M_on_contour']))
json.dump(dict(results=res, comparison=cmp), open('winding_subset.json', 'w'), indent=2)
print(json.dumps(cmp, indent=1))
