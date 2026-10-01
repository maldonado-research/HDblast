#!/usr/bin/env python3
"""B2 diagnosis D5: why the main runs fail after the switch to full-field mode, and what resolution would be needed.

At T_switch the solver stops evolving the deviation from the exact (analytic) seed and evolves the full fields with 4th-order
finite differences.  The finite-difference truncation error of the steep static background then acts as a source.  Here the
EXACT seed solution (static shell, tension c+dc, in the bounded chart) is put on each grid in full-field form at a given chart
time T and we evaluate, without evolving:
  * the relative Hamiltonian / momentum constraint residuals (as recorded by the solver; near shell z > -1 and domain),
  * the relative residual of the evolution equations |accel_FD(S) - S_TT| / (|S_TT| + |terms|) (max over z > -1).
For an exact solution all of these vanish in the continuum; on the grid they are the truncation error that the full-field run
starts from.  Compared with: the same quantities at T = 0 (identity-chart part of the slice) and the A1 delta = 0.1 main-run
configuration (where full-field mode ran stably to a plateau).  Output diag/D5_FULLMODE_TRUNCATION.json.  Numerical."""
import json, math, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d')
import evolve_a1 as E
import static_w as S
M8 = json.load(open('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json'))

def probe(run, T):
    rf = run.ref(T)
    run.mode = 'full'
    Yv = np.array([rf['A'], rf['B'], rf['phi'], rf['A_T'], rf['B_T'], rf['phi_T']])
    F = run.fields(T, Yv, 0.0)
    ra, rb, rfv = run.accel(F)
    eBb = math.exp(F['B'][-1]); RT = run.Y*F['phi_T'][-1]**2/eBb
    H, M, sH, sM = run.constraints(F, Yv[3], 0.0, RT, rf['phi_TT'][-1])
    ii = run.inner; near = (run.z > -1.0)[ii]
    hr = (np.abs(H)/sH)[ii]; mr = (np.abs(M)/sM)[ii]
    sc = np.abs(rf['A_TT']) + np.abs(F['A_ZZ']) + 3*F['A_T']**2 + 3*F['A_Z']**2 + 1e-300
    ea = (np.abs(ra - rf['A_TT'])/sc)[ii]
    scp = np.abs(rf['phi_TT']) + np.abs(F['phi_ZZ']) + 1e-300
    ep = (np.abs(rfv - rf['phi_TT'])/scp)[ii]
    # shell Neumann data vs exact Z-derivatives (what the record's Weyl estimate sees)
    dAz = float(abs(F['A_Z'][-1] - rf['A_Z'][-1])/abs(rf['A_Z'][-1]))
    zi = run.z[ii]
    return dict(T=T, lapse_shell=float(math.exp(rf['B'][-1])/run.rb), H_near=float(hr[near].max()), M_near=float(mr[near].max()),
                z_H_near=float(zi[near][np.argmax(hr[near])]), H_dom=float(hr.max()), M_dom=float(mr.max()),
                accelA_near=float(ea[near].max()), accelphi_near=float(ep[near].max()), z_accel=float(zi[near][np.argmax(ea[near])]),
                rel_err_A_Z_shell=dAz)

out = dict(status='numerical', purpose=__doc__, cases=[])
cases = []
for (lab, dzf, dzc) in [('S1', 3e-4, 3e-2), ('S2', 1.5e-4, 1.5e-2), ('S3', 7.5e-5, 7.5e-3), ('S4', 3.75e-5, 3.75e-3)]:
    cases.append(dict(label='delta1e-3 ' + lab, delta=1e-3, dc=1e-2, dzf=dzf, dzc=dzc, L=24.0, grid='rho', xc=11.3, Ts=[0.0, 8.0, 10.3], table_dx=2.5e-5))
cases.append(dict(label='delta0.1 A1 main dzf5e-4', delta=0.1, dc=1e-2, dzf=5e-4, dzc=2e-3, L=16.0, grid='tanh', xc=3.1, Ts=[0.0, 2.1], table_dx=2.5e-4))
cases.append(dict(label='delta0.1 A1 main dzf1e-3', delta=0.1, dc=1e-2, dzf=1e-3, dzc=4e-3, L=16.0, grid='tanh', xc=3.1, Ts=[0.0, 2.1], table_dx=2.5e-4))
shells = {}
for c in cases:
    t0 = time.time()
    key = (c['delta'], c['dc'], c['table_dx'])
    d = M8['d_star'][{0.1: '0.1', 0.001: '0.001'}[c['delta']]]
    run = E.Run(delta=c['delta'], d=d, dc=c['dc'], Y=1.0, dzf=c['dzf'], dzc=c['dzc'], L=c['L'], xc=c['xc'], kappa=10.0, grid_kind=c['grid'],
                table_dx=c['table_dx'], shells=shells.get(key), log=lambda s: None)
    shells[key] = (run.bg, run.seed)
    res = dict(label=c['label'], n_points=run.n, points_per_wall=float((1/run.rb)/run.zp[-1]), rho_b=run.rb, dzf=c['dzf'], xc=c['xc'],
               T_switch=run.T_switch, probes=[probe(run, T) for T in c['Ts']], seconds=round(time.time() - t0, 1))
    out['cases'].append(res)
    print(json.dumps(res), flush=True)
(HERE/'diag/D5_FULLMODE_TRUNCATION.json').write_text(json.dumps(out, indent=1) + '\n')
