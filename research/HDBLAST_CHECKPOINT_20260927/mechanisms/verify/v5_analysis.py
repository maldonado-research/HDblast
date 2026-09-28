"""Audit V5 analysis: reliable window (m6 rule: first t>1 with d(H0 tau)/dt < 0.2) of the signed c* pilot run.
Output: V5_CSTAR_PILOT_SIGNED.json"""
import json, numpy as np
from pathlib import Path
H = Path(__file__).resolve().parent
s = json.load(open(H/'pilot/v5_cstar_Y0_signed_summary.json')); d = np.load(H/'pilot/v5_cstar_Y0_signed_timeseries.npz')['rec']
t, ss = d[:, 0], d[:, 8]; l = np.gradient(ss, t); bad = np.where((l < 0.2) & (t > 1))[0]; ie = int(bad[0]) - 1 if len(bad) else len(t)-1
rb = s['rho_b']; c = s['c']; W = lambda p: 1 - p + p**3/3
p = d[ie, 1]; Hv2_series_over_H02 = (0.1*(1+c)/27 + 0.01*((1+c)**2/36 - c*c/384))*rb*rb
out = dict(summary=s, window_end=dict(t=float(t[ie]), H0tau=float(ss[ie]), phi_b=float(p), H_over_H0=float(d[ie, 2]),
           Wy_over_H0sq=float(d[ie, 7]), Wy_over_H2=float(d[ie, 7]/d[ie, 2]**2)),
           phi_b_min=float(d[:, 1].min()), phi_b_max=float(d[:, 1].max()), crossed_phi_minus1=bool((d[:, 1] > -1).any() and d[0, 1] < -1),
           Hvac2_series_over_H02_at_cstar=float(Hv2_series_over_H02),
           samples=[[float(t[i]), float(ss[i]), float(d[i, 1]), float(d[i, 2]), float(d[i, 7]), float(l[i])] for i in range(0, len(t), 50)],
           caveats='single grid (dz_fine=1e-3, L=17), single seed sign (dc=+1e-2), Y=0; values after the window end are not physical by the m6 rule')
(H/'V5_CSTAR_PILOT_SIGNED.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps({k: v for k, v in out.items() if k != 'samples' and k != 'summary'}, indent=1))
