#!/usr/bin/env python3
"""B4 step 4 (REGISTRATION.md sections 3.3-3.4): fine-tuning measure and recollapse estimate.

Late sector (exact once the shell scalar is frozen; v = 0):  H^2/H0^2 = lam + rad_s e^{-4N} (+ Weyl_s e^{-4N}),
N = ln(a/a_s), lam = Lambda_res/H0^2.  Radiation-dominated (|Omega_vac| <= 0.1) for N e-folds after a_s requires
|lam| <= 0.1 rad_s e^{-4N}  (Weyl neglected; it only adds to the a^-4 term).
For lam < 0 the exact solution is a^2 = a_*^2 sin(2 sqrt|lam| (tau - tau_0)), so H = 0 at tau_0 + pi/(4 sqrt|lam|) and the
crunch is at tau_0 + pi/(2 sqrt|lam|) (checked symbolically and against an ODE below).
Tolerance on d: |d - d*_exact| <= |lam_tol| / (dLambda/dd), slope from the exact static branch (B4_STATIC_LAMBDA.json) at
delta = 0.1 and from the series delta/54 with H0^2 -> 0.16091 delta for delta -> 0.
Inputs labelled: rad_s from the A1 5D runs (numerical, read-only) and from the reduced model (conditional).
Output: B4_FINETUNING.json."""
import json, math, hashlib, sys
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum')
ST = json.loads((HERE/'B4_STATIC_LAMBDA.json').read_text())
PR = json.loads((HERE/'B4_PREDICTIONS.json').read_text())
C = 2/1.0357712571566784 - 4/3

def a1_stop(tag):
    d = np.load(A1/'runs'/'main'/(tag + '_timeseries.npz')); cols = list(d['cols']); R = d['rec']
    g = lambda k: R[:, cols.index(k)]
    i = int(np.flatnonzero((np.abs(g('v_over_H0')) < 0.01) & (g('phi_b') > 1.0))[0])
    return dict(H0tau=float(g('H0tau')[i]), rad_s=float(g('rad')[i]), Weyl_s=float(g('Wy')[i]), vac_s=float(g('vac')[i]), H_s=float(g('H_over_H0')[i]))

def main():
    out = dict(status='exact algebra (late-sector formulas) + numerical inputs; conditional where the reduced model is used')
    # ---- exact algebra: Lambda < 0 radiation-Lambda solution
    t, L, K = sp.symbols('t L K', positive=True)
    a2 = sp.sqrt(K/L)*sp.sin(2*sp.sqrt(L)*t)              # u = a^2 ;  (u'/2)^2 = K - L u^2  <=>  H^2 = K/a^4 - L
    chk = sp.simplify((sp.diff(a2, t)/2)**2 - (K - L*a2**2))
    out['exact_check_rad_negLambda'] = dict(residual=str(chk), H_zero_at='pi/(4 sqrt|Lambda|)', crunch_at='pi/(2 sqrt|Lambda|)', PASS=bool(chk == 0))
    # ODE cross-check
    lam = -3.4e-4; Ks = 0.02
    sol = solve_ivp(lambda tt, y: [y[1], -2*Ks*math.exp(-4*y[0])], (0, 200), [0.0, math.sqrt(Ks + lam)], rtol=1e-11, atol=1e-13,
                    events=lambda tt, y: y[1], dense_output=True)
    # time from a = 0 to a_s: tau_s = asin(sqrt(|lam|/Ks))/(2 sqrt|lam|)
    tz_pred = (math.pi/2 - math.asin(math.sqrt(-lam/Ks)))/(2*math.sqrt(-lam))
    out['ode_check'] = dict(lam=lam, rad_s=Ks, H_zero_ode=float(sol.t_events[0][0]), H_zero_closed=tz_pred, rel=float(abs(sol.t_events[0][0] - tz_pred)/tz_pred))
    # ---- inputs
    f3 = ST['fits']['deg3']; slope = f3['slope_at_dstar']; dex = f3['d_star_exact']
    lam_dstarM8 = f3['at']['1.00']['Lambda_over_H0sq']
    slope0 = 1/(54*0.16091)            # delta -> 0: dLambda/dd = delta/54, H0^2 = 0.16091 delta (EFT / M2 limit)
    dstar0 = -2*(1 + C)
    stops = {tag: a1_stop(tag) for tag in ('main_dstar_Y1_dc1e-2_dzf5e-4', 'main_dstar_Y0.3_dc1e-2_dzf5e-4')}
    out['inputs'] = dict(slope_delta0p1=slope, d_star_exact_delta0p1=dex, d_star_M8=-3.106933495673783, lam_at_dstar_M8=lam_dstarM8,
                         d_star_rel_offset=(-3.106933495673783 - dex)/dex, slope_delta_to_0=slope0, d_star_delta_to_0=dstar0, A1_stop=stops,
                         reduced_rad_match={Y: v['radonly']['rad_match'] for Y, v in PR['Y_scan'].items()})
    print(out['inputs'])
    # ---- (4) A1 d* run: radiation-dominated e-folds and recollapse estimate
    rs = stops['main_dstar_Y1_dc1e-2_dzf5e-4']['rad_s']; ws = stops['main_dstar_Y1_dc1e-2_dzf5e-4']['Weyl_s']; ts = stops['main_dstar_Y1_dc1e-2_dzf5e-4']['H0tau']
    lamA = lam_dstarM8
    N_RD = 0.25*math.log(0.1*rs/abs(lamA))
    Kt = rs + ws
    t_since_origin_s = math.asin(math.sqrt(abs(lamA)/Kt))/(2*math.sqrt(abs(lamA)))
    out['A1_dstar_Y1'] = dict(lam=lamA, rad_s=rs, Weyl_s=ws, tau_s=ts, N_RD_e_folds=N_RD,
                              H_zero_H0tau_estimate=ts - t_since_origin_s + math.pi/(4*math.sqrt(abs(lamA))),
                              crunch_H0tau_estimate=ts - t_since_origin_s + math.pi/(2*math.sqrt(abs(lamA))),
                              Omega_vac_when_rad_diluted_by_e4=lamA/(lamA + Kt*math.exp(-4)))
    for qq, key in ((1.01, '1.01'), (1.05, '1.05')):
        l2 = f3['at'][key]['Lambda_over_H0sq']; Kt2 = Kt
        if abs(l2) < Kt2:
            t0 = math.asin(math.sqrt(abs(l2)/Kt2))/(2*math.sqrt(abs(l2)))
            out['A1_dstar_Y1']['H_zero_estimate_q%s' % key] = ts - t0 + math.pi/(4*math.sqrt(abs(l2)))
    print(out['A1_dstar_Y1'])
    # ---- (3a) one radiation-dominated e-fold, (3b) radiation era to BBN
    def tol(rad_s, N=None, rho_ratio=None, eps=0.1):
        lam_tol = eps*rad_s*(math.exp(-4*N) if N is not None else rho_ratio)
        return dict(lam_tol=lam_tol, dd_tol_delta0p1=lam_tol/slope, rel_dd_tol_delta0p1=lam_tol/slope/abs(dex),
                    dd_tol_delta_to_0=lam_tol/slope0, rel_dd_tol_delta_to_0=lam_tol/slope0/abs(dstar0))
    out['tolerance_one_efold'] = {'A1_5D_Y1': tol(rs, N=1), 'A1_5D_Y0.3': tol(stops['main_dstar_Y0.3_dc1e-2_dzf5e-4']['rad_s'], N=1),
                                  **{'reduced_Y%s' % Y: tol(v, N=1) for Y, v in out['inputs']['reduced_rad_match'].items()}}
    TB = 1e-3; gB = 10.75                               # BBN at T = 1 MeV, g* = 10.75 (SM)
    rho_L = 0.685*8.0992e-47*0.674**2                  # GeV^4, h = 0.674, Omega_Lambda = 0.685 (as in M2/M5)
    rows = []
    for TRH, g in ((1e-2, 10.75), (1.0, 61.75), (1e3, 106.75), (1e10, 106.75), (1e16, 106.75)):
        rr = (gB*TB**4)/(g*TRH**4)                      # rho_r(BBN)/rho_r(T_RH), entropy and g* changes folded into rho ~ g* T^4
        N = math.log(TRH/TB) + math.log((g/gB)**(1/3))
        rho_RH = math.pi**2/30*g*TRH**4
        obs = rho_L/rho_RH                              # observed vacuum relative to radiation at T_RH
        rows.append(dict(T_RH_GeV=TRH, g_star=g, efolds_RH_to_BBN=N, rho_BBN_over_rho_RH=rr, **tol(rs, rho_ratio=rr),
                         observed_Lambda_over_rho_r_RH=obs, rel_dd_for_observed_Lambda_delta0p1=obs*rs/slope/abs(dex)))
    out['tolerance_BBN'] = rows
    out['observed_vacuum_GeV4'] = rho_L
    out['planck_units_statement'] = 'about 1e-122 in Planck units (LITERATURE_2022_2026.md, dark-dimension entry; snippet level)'
    for r in rows: print('T_RH=%.0e N=%.1f lam_tol=%.2e rel_dd=%.2e obs=%.2e rel_dd_obs=%.2e' % (r['T_RH_GeV'], r['efolds_RH_to_BBN'], r['lam_tol'], r['rel_dd_tol_delta0p1'], r['observed_Lambda_over_rho_r_RH'], r['rel_dd_for_observed_Lambda_delta0p1']))
    # ---- numerical check of the one-e-fold tolerance with the reduced model (radiation-only variant), Y = 1
    import eft_reduced as M
    from b4_predictions import one
    e = M.EFT(); qex = dex/M.DSTAR
    chk = []
    for dq in (-4e-4, -2e-4, -1e-4, 0.0, 1e-4, 2e-4, 4e-4):
        o = one(e, 1.0, qex + dq, 'radonly'); chk.append(dict(q=qex + dq, dq=dq, efolds_RD=o['efolds_RD'], classification=o['classification']))
        print('q=%.6f efolds_RD=%.3f %s' % (qex + dq, o['efolds_RD'], o['classification']))
    out['reduced_check_RD_efolds_near_dstar_exact'] = dict(q_exact=qex, rows=chk,
        formula_dq_for_1_efold=0.1*out['inputs']['reduced_rad_match']['1']*math.exp(-4)/slope/abs(M.DSTAR))
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'B4_FINETUNING.json').write_text(json.dumps(out, indent=1) + '\n')

if __name__ == '__main__':
    main()
