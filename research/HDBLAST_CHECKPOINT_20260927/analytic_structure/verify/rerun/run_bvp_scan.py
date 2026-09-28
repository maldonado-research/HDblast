"""Run the high-precision BVP at several detunings / values of c.  Writes BVP_SCAN_<label>.json.
usage: python3 run_bvp_scan.py <label> <c: 'reg' or decimal> <dps> <order> <h> <u_cone> <delta1,delta2,...>
"""
import sys, json, time, hashlib
from pathlib import Path
import mpmath as mp
import sympy as sp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mp_bvp import Problem

label, cstr, dps, order, h, ucone, dl = sys.argv[1:8]
dps = int(dps); order = int(order)
mp.mp.dps = dps + 10
if cstr == 'reg':
    c = mp.mpf(2) / mp.mpf('1.0357712571566784') - mp.mpf(4) / 3
else:
    c = mp.mpf(cstr)
deltas = [mp.mpf(x) for x in dl.split(',')]
ser = json.loads((HERE / 'SERIES_COEFFICIENTS.json').read_text())
cs = sp.Symbol('c')
etah_coef = [sp.sympify(x) for x in ser['eta_h_over_delta8']]
rows = []
for d in deltas:
    t0 = time.time()
    guess = sum(mp.mpf(str(sp.N(e.subs(cs, sp.Float(str(c), 30)), 30))) * d ** i for i, e in enumerate(etah_coef[:3])) * d ** 8
    P = Problem(d, c, dps=dps, order=order, h=mp.mpf(h), u_cone=mp.mpf(ucone))
    out = P.solve(guess if c != 0 else mp.mpf(0))
    row = {k: mp.nstr(v, dps) for k, v in out.items()}
    row['runtime_s'] = time.time() - t0
    rows.append(row)
    print(json.dumps({k: row[k] for k in ['delta', 'eta_b', 'H2', 'J2', 'first_integral_residual', 'runtime_s']}), flush=True)
res = dict(status='NUMERICAL high-precision shooting (Taylor-series integrator, mpmath)', label=label, c=mp.nstr(c, dps),
           dps=dps, taylor_order=order, step_h=h, cone_series_start_u=ucone, cone_series_order=90,
           mp_version=mp.__version__, rows=rows,
           solver_sha256=hashlib.sha256((HERE / 'mp_bvp.py').read_bytes()).hexdigest())
(HERE / 'runs').mkdir(exist_ok=True)
(HERE / 'runs' / f'BVP_SCAN_{label}.json').write_text(json.dumps(res, indent=2) + '\n')
