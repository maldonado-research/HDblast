"""Compare re-run JSON outputs (verify/rerun) with the workstream's saved JSON (read-only)."""
import json, sys, math
from pathlib import Path
M = Path(__file__).resolve().parent.parent; V = Path(__file__).resolve().parent/'rerun'
SKIP = ('runtime', 'runtime_s')
def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if k in SKIP: continue
            if k not in a or k not in b: out['missing'].append(path+'/'+k); continue
            walk(a[k], b[k], path+'/'+k, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out['len_mismatch'].append(path); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, path+'[%d]' % i, out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        d = abs(a-b)/max(abs(a), abs(b), 1e-300); out['n'] += 1
        if d > out['max_rel'][0]: out['max_rel'] = [d, path, a, b]
    else:
        out['n'] += 1
        if a != b: out['other_diff'].append([path, str(a)[:80], str(b)[:80]])
res = {}
for name in sys.argv[1:]:
    a = json.load(open(M/name)); b = json.load(open(V/name))
    out = dict(n=0, max_rel=[0.0, None, None, None], missing=[], len_mismatch=[], other_diff=[])
    walk(a, b, '', out); res[name] = out
    print(name, 'n=%d max_rel=%.3e at %s' % (out['n'], out['max_rel'][0], out['max_rel'][1]), 'other_diff=%d missing=%d' % (len(out['other_diff']), len(out['missing'])))
    for d in out['other_diff'][:5]: print('   ', d)
json.dump(res, open(Path(__file__).resolve().parent/('COMPARE_'+'_'.join(n.split('_')[0] for n in sys.argv[1:])+'.json'), 'w'), indent=1)
