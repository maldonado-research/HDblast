#!/usr/bin/env python3
"""B2 diagnosis D1: does the static reference (cubic-Hermite tables in x = ln|w|) satisfy the shell junctions accurately
enough for the delta = 1e-3 growth calibration (seed dc = 1e-8)?

For the registered model (d = 0), delta = 1e-3, tension c (background) and c + dc (seed), evaluate the reference at the shell
(T = 0, Z = 0, identity chart) and compute the junction residuals
    eA   = e^{-B} A_Z - sigma(phi)/6,     ephi = e^{-B} phi_Z + sigma'(phi)/2
for table spacings dx = 2.5e-4 (A1 default) and finer.  The physical mismatch that drives the calibration run is the tension
difference: e^{-B} phi_Z(seed) + sigma'_c(phi)/2 = -delta*dc/2 (= -5e-12 at dc = 1e-8).  A table error ephi is equivalent to a
spurious tension shift dc_eff = -2 ephi / delta.  Output: diag/D1_TABLE_KICK.json.  Numerical (floating point)."""
import json, math, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
import numpy as np
import static_w as S
import rolloff5d_v1 as R_

delta, c = 1e-3, S.C_REG
out = dict(status='numerical', purpose=__doc__, delta=delta, c=c, cases=[])
g0 = (math.log10(0.2207*delta**1.8), 8.2282 + 0.9*math.log(1e-3/delta))
ph_h, y_b, _ = R_.solve_shell(R_.Tension(delta, c, 0.0), guess=g0)
for dc in [0.0, 1e-8]:
    for dx in [2.5e-4, 1e-4, 5e-5, 2.5e-5]:
        t0 = time.time()
        ten = S.Tension(delta, c + dc, 0.0)
        sh = S.StaticShell(ten, phh_guess=ph_h, table_dx=dx)
        ch = S.Chart(math.inf)
        rf = S.reference(sh, ch, np.array([0.0]), np.array([0.0]))
        B = rf['B'][0]; ph = rf['phi'][0]; eB = math.exp(-B)
        ten_c = S.Tension(delta, c, 0.0)
        eA = eB*rf['A_Z'][0] - ten.s(ph)/6
        ephi = eB*rf['phi_Z'][0] + ten.s1(ph)/2
        ephi_vs_c = eB*rf['phi_Z'][0] + ten_c.s1(ph)/2
        # table node offset of the shell point
        xb = sh.xb; tb = sh.tab[1]; frac = ((xb - tb['x0'])/tb['dx']) % 1.0
        row = dict(dc=dc, table_dx=dx, rho_b=sh.rhob, phi_b=sh.phb, x_b=xb, shell_node_fraction=frac,
                   dense_junction_residuals=list(sh.junction_residuals),
                   table_eA=eA, table_ephi=ephi, dc_eff_from_table=-2*ephi/delta,
                   ephi_against_background_tension=ephi_vs_c, physical_mismatch_expected=-delta*dc/2,
                   phi_b_table=ph, phi_b_dense=sh.phb, dphi_b_table=ph - sh.phb, seconds=round(time.time() - t0, 1))
        out['cases'].append(row); print(json.dumps(row), flush=True)
json.dump(out, open(HERE/'diag/D1_TABLE_KICK.json', 'w'), indent=1, default=float)
