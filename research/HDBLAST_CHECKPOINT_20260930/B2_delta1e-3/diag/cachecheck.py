import sys, math; sys.path.insert(0,'.')
sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
import numpy as np, static_w as S, rolloff5d_v1 as R_
delta=1e-3; c=S.C_REG
g0 = (math.log10(0.2207*delta**1.8), 8.2282 + 0.9*math.log(1e-3/delta))
ph_h, y_b, _ = R_.solve_shell(R_.Tension(delta, c, 0), guess=g0)
sh = S.StaticShell(S.Tension(delta, c, 0), phh_guess=ph_h, table_dx=2.5e-5)
ch = S.Chart(math.inf)
z = -np.linspace(0, 10, 2000)
r0 = S.reference(sh, ch, np.zeros_like(z), z); r1 = S.reference(sh, ch, np.full_like(z, 1.234), z)
for k in r0:
    d = r1[k] - r0[k] - (1.234 if k == 'A' else 0)
    i = int(np.argmax(np.abs(d))); print(k, '%.3e at z=%.4f  (|val| %.3e)' % (abs(d[i]), z[i], abs(r0[k][i])))
