"""Audit V5: attempt the workstream's failed c* pilot (delta=0.1, c=-0.9930694751920546, Y=0, dc=1e-2, L=17) using
the SAME pilot code (rolloff5d_matter.py, unmodified copy) but with the static-shell root finder allowed to use
phi_h < -1 (sign-flipped parameterisation, audit V3) and the phi_b stop window widened to [-3, 1.3].
Question: which way does the continued c* shell roll?  Time-boxed; exploratory only."""
import math, sys, numpy as np
sys.path.insert(0, '..')
import rolloff5d_v1 as R_
import rolloff5d_matter as P
from v3_c_family_through_zero import solve_signed
def solve_shell_signed(ten, guess=None):
    # guess from V3B at c*: phi_h = -1.006946, y_b from JSON
    import json
    g = json.load(open('../V3B_CONTINUE_TO_CSTAR.json'))['rows'][-1]
    x, r = solve_signed(ten, (math.log10(-1 - g['phi_h']), g['y_b']), -1)
    print('shell root c=%.6f phi_h=%.6f y_b=%.5f residual=%.1e' % (ten.c, -1 - 10**x[0], x[1], max(abs(r))), flush=True)
    return -1 - 10**x[0], x[1], r
R_.solve_shell = solve_shell_signed
import types
_m = types.SimpleNamespace(**{k: getattr(math, k) for k in dir(math) if not k.startswith('_')})
_m.log10 = lambda v: math.log10(abs(v)) if v != 0 else -300.0   # the pilot only uses log10(ph_h+1) as an (ignored) guess
P.math = _m
P.R_.solve_shell = solve_shell_signed
import types
_m = types.SimpleNamespace(**{k: getattr(math, k) for k in dir(math) if not k.startswith('_')})
_m.log10 = lambda v: math.log10(abs(v)) if v != 0 else -300.0   # the pilot only uses log10(ph_h+1) as an (ignored) guess
P.math = _m
tf = float(sys.argv[1]) if len(sys.argv) > 1 else 15.0
P.run(Y=0.0, dc=1e-2, L=17.0, t_final=tf, t_det=0.1, tag='v5_cstar_Y0_signed', stop_phi=(-3.0, 1.3),
      log=lambda s: print(s, flush=True), c_override=-0.9930694751920546)
