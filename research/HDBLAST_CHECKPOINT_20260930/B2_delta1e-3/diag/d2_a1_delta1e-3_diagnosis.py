#!/usr/bin/env python3
"""B2 diagnosis D2: why the A1 delta = 1e-3 runs did not complete.  Reads only the archived A1 outputs (read-only):
  runs/pre/pre_d1e-3_dstar_Y1_dc1e-2*, runs/main/d3_*, runs/cal/cal_growth_d0.001*   (A1 folder, 28 Sept 2026)
and the B2 table diagnosis diag/D1_TABLE_KICK.json.  Output: diag/D2_A1_DELTA1E-3_DIAGNOSIS.json.  Numerical (floating point).

Questions: (1) why the pre-run never met d b_b/dT < -0.9 (so xc = inf); (2) why the two spacings stopped at H0 tau 3.1 and 11.0;
(3) why the growth-rate calibration C2 drifted; (4) the growth rate of the tuned Y = 1 runs at delta = 1e-3 (sets the run length)."""
import json, math
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum')

def load(stem):
    s = json.load(open(str(stem) + '_summary.json'))
    d = np.load(str(stem) + '_timeseries.npz')
    ts = {c: d['rec'][:, i] for i, c in enumerate(list(d['cols']))}
    return s, ts, d

def first(mask, arr):
    i = np.flatnonzero(mask)
    return (float(arr[i[0]]) if len(i) else None)

def logslope(t, y, a, b):
    m = (t >= a) & (t <= b) & np.isfinite(y) & (y > 0)
    return float(np.polyfit(t[m], np.log(y[m]), 1)[0]) if m.sum() > 10 else None

out = dict(status='numerical', purpose=__doc__, inputs={}, pre_run={}, main_runs={}, growth_calibration={}, tuned_growth_rate={}, conclusions={})

# ---------------------------------------------------------------- (1) pre-run
s, ts, _ = load(A1/'runs/pre/pre_d1e-3_dstar_Y1_dc1e-2')
T, Bb, lapse = ts['T'], ts['B_b'], ts['lapse']
dbdt = np.gradient(Bb, T)
out['inputs']['pre'] = dict(params={k: s['params'][k] for k in ['dzf', 'dzc', 'zfine', 'L', 'kappa', 'project', 'dc', 'Y']}, n_points=s['n_points'],
                            points_per_wall=(1/s['rho_b'])/s['params']['dzf'], stop_reason=s['stop_reason'], T_end=s['T_end'], H0tau_end=s['H0tau_end'])
hn = ts['Hmax_near']; mn = ts['Mmax_near']; ok = np.isfinite(hn)
out['pre_run'] = dict(
    min_dbdt=float(np.nanmin(dbdt[:-3])), max_dbdt=float(np.nanmax(dbdt[:-3])),
    criterion_dbdt_lt_minus_0p9_met=bool(np.nanmin(dbdt[:-3]) < -0.9),
    lapse_at_T={('%.1f' % t0): float(np.interp(t0, T, lapse)) for t0 in [0, 1, 2, 3, 4, 5, 5.5, 5.6]},
    lapse_monotone_increase_from_T=first((np.gradient(lapse, T) > 0) & (T > 1.0), T),
    first_T_lapse_gt_1p1=first(lapse > 1.1, T), first_T_lapse_gt_2=first(lapse > 2.0, T),
    near_shell_H_rel_residual={('%.2f' % t): float(h) for t, h in zip(T[ok], hn[ok])},
    near_shell_M_rel_residual={('%.2f' % t): float(h) for t, h in zip(T[ok], mn[ok])},
    first_T_near_shell_residual_gt_0p05=first(ok & ((hn > 0.05) | (mn > 0.05)), T),
    phi_b_at_stop=float(ts['phi_b'][-1]), H_over_H0_at_stop=float(ts['H_over_H0'][-1]))

# ---------------------------------------------------------------- (2) the two main runs
for tag in ['d3_dstar_Y1_dc1e-2_dzf2e-4', 'd3_dstar_Y1_dc1e-2_dzf1.5e-4']:
    s, ts, dd = load(A1/'runs/main'/tag)
    T = ts['T']; ok = np.isfinite(ts['Hmax_near'])
    hn, mn = ts['Hmax_near'], ts['Mmax_near']
    pl = s['projection_log']
    res = pl[:-1] if pl else []
    mono = all(res[i + 1]['max_rel_H_window'] <= res[i]['max_rel_H_window'] for i in range(len(res) - 1)) if res else None
    # location of the largest near-shell residual (z of the domain maximum when it is near the shell)
    zmax = ts['Hmax_z']
    late = ok & (T > 4)
    r = dict(params={k: s['params'][k] for k in ['dzf', 'dzc', 'zfine', 'width', 'L', 'kappa', 'project', 'damp_off']}, n_points=s['n_points'],
             points_per_wall=(1/s['rho_b'])/s['params']['dzf'], stop_reason=s['stop_reason'], T_end=s['T_end'], H0tau_end=s['H0tau_end'],
             projection=dict(initial_rel_H=res[0]['max_rel_H_window'] if res else None, final_rel_H=pl[-1]['max_rel_H_window'] if pl else None,
                             iterations=[x['max_rel_H_window'] for x in res], monotone_decrease=mono,
                             max_abs_warp_shift=pl[-1].get('max_abs_shift') if pl else None),
             at_T0=dict(H_over_H0=float(ts['H_over_H0'][0]), Wy=float(ts['Wy'][0]), lapse=float(ts['lapse'][0]), Hmax_domain=float(ts['Hmax'][0])),
             lapse_at_T={('%.1f' % t0): float(np.interp(t0, T, ts['lapse'])) for t0 in np.arange(0, T[-1], 1.0)},
             first_T_near_shell_residual_gt_0p05=first(ok & ((hn > 0.05) | (mn > 0.05)), T),
             first_H0tau_near_shell_residual_gt_0p05=first(ok & ((hn > 0.05) | (mn > 0.05)), ts['H0tau']),
             z_of_domain_max_H_residual_T_gt_4=dict(median=float(np.median(zmax[late])) if late.any() else None,
                                                    min=float(np.min(zmax[late])) if late.any() else None, max=float(np.max(zmax[late])) if late.any() else None),
             near_shell_M_rel_residual={('%.1f' % t): float(m) for t, m in zip(T[ok][::3], mn[ok][::3])},
             Wy_range_last_2_T=[float(np.min(ts['Wy'][T > T[-1] - 2])), float(np.max(ts['Wy'][T > T[-1] - 2]))],
             H_range_last_2_T=[float(np.min(ts['H_over_H0'][T > T[-1] - 2])), float(np.max(ts['H_over_H0'][T > T[-1] - 2]))])
    # grid spacing at the location of the residual (reconstruct the A1 tanh grid)
    import importlib.util
    spec = importlib.util.spec_from_file_location('ev', str(A1/'evolve_a1.py')); ev = importlib.util.module_from_spec(spec); spec.loader.exec_module(ev)
    z, zp, _ = ev.build_grid(s['params']['L'], s['params']['dzf'], s['params']['dzc'], s['params']['zfine'], s['params']['width'])
    r['grid_dz_at'] = {str(zz): float(np.interp(zz, z, zp)) for zz in [0.0, -0.02, -0.04, -0.05, -0.06, -0.08, -0.1, -0.2]}
    r['grid_ratio_coarse_to_fine'] = s['params']['dzc']/s['params']['dzf']
    # static profile scales from the T = 0 snapshot (every 4th node): B slope and phi tail near z = -0.05
    zz = dd['z']; B0 = dd['snap_0.000_B']; ph0 = dd['snap_0.000_phi']
    for z0 in [-0.01, -0.05, -0.1, -0.5, -1.0]:
        i = int(np.argmin(np.abs(zz - z0)))
        Bz = (B0[i + 1] - B0[i - 1])/(zz[i + 1] - zz[i - 1])
        r.setdefault('static_profile', {})[str(z0)] = dict(B=float(B0[i]), e_B=float(math.exp(B0[i])), B_Z=float(Bz), phi=float(ph0[i]),
                                                          local_scale_1_over_eB=float(math.exp(-B0[i])), wall_width_scaled=float((1/s['rho_b'])*s['rho_b']/math.exp(B0[i])),
                                                          grid_points_per_scaled_wall=float((1/s['rho_b'])*s['rho_b']/math.exp(B0[i])/np.interp(z0, z, zp)))
    out['main_runs'][tag] = r

# ---------------------------------------------------------------- (4) growth rate of the tuned Y = 1 runs at delta = 1e-3 (friction-limited)
for tag, stem in [('pre (dzf 4e-4)', A1/'runs/pre/pre_d1e-3_dstar_Y1_dc1e-2'), ('main dzf 2e-4', A1/'runs/main/d3_dstar_Y1_dc1e-2_dzf2e-4')]:
    s, ts, _ = load(stem)
    dev = np.abs(ts['phi_b'] - s['phi_b_static'])
    rates = {('%.0f-%.0f' % (a, a + 2)): logslope(ts['T'], dev, a, a + 2) for a in [1, 2, 3, 4, 5, 6, 7, 8] if a + 2 <= ts['T'][-1] - 0.5}
    out['tuned_growth_rate'][tag] = dict(sliding_2T_windows_in_chart_time=rates, phi_b_static=s['phi_b_static'], seed_phi_b=s['seed_phi_b'])
ra = [v for v in out['tuned_growth_rate']['main dzf 2e-4']['sliding_2T_windows_in_chart_time'].values() if v is not None]
mu = float(np.median(ra[:5]))
seed_dev_1e2 = out['tuned_growth_rate']['main dzf 2e-4']['seed_phi_b'] - out['tuned_growth_rate']['main dzf 2e-4']['phi_b_static']
out['tuned_growth_rate']['median_rate_T1_to_7'] = mu
out['tuned_growth_rate']['e_folds_from_dc1e-2_seed_to_phi_b_0p5'] = math.log(0.5/seed_dev_1e2)
out['tuned_growth_rate']['T_estimate_to_phi_b_0p5_dc1e-2'] = math.log(0.5/seed_dev_1e2)/mu
out['tuned_growth_rate']['extra_T_for_dc1e-4_seed_(linear)'] = math.log(100)/mu
out['tuned_growth_rate']['note'] = ('seed deviation scales linearly with dc (1e-4 seed ~100x smaller); rates are log-slopes of |phi_b - phi_b(static)| '
                                    'in chart time (lapse within 0.3% of 1 in these windows)')

# ---------------------------------------------------------------- (3) growth calibration C2
s, ts, _ = load(A1/'runs/cal/cal_growth_d0.001')
tau = ts['H0tau']; dev = np.abs(ts['phi_b'] - s['phi_b_static'])
sl = {('%.1f-%.1f' % (a, a + 1)): logslope(tau, np.where(dev < 1e-3, dev, np.nan), a, a + 1) for a in np.arange(1.0, 5.01, 0.5)}
D1 = json.load(open(HERE/'diag/D1_TABLE_KICK.json'))
kick = [c for c in D1['cases'] if c['dc'] == 1e-8 and c['table_dx'] == 2.5e-4][0]
out['growth_calibration'] = dict(
    params={k: s['params'][k] for k in ['dzf', 'dzc', 'zfine', 'L', 'dc', 'kappa', 'project']}, points_per_wall=(1/s['rho_b'])/s['params']['dzf'],
    A1_fit=dict(rate=1.6978203383432817, window=[1.0095, 5.2053], target=1.65719),
    sliding_1_unit_windows_dev_lt_1em3=sl,
    dev_initial=float(dev[0]), dev_at_tau={('%.0f' % t0): float(np.interp(t0, tau, dev)) for t0 in range(0, 8)},
    table_kick=dict(table_dx=2.5e-4, junction_error_phi=kick['table_ephi'], junction_error_A=kick['table_eA'],
                    physical_mismatch_phi=kick['physical_mismatch_expected'], ratio=kick['table_ephi']/kick['physical_mismatch_expected'],
                    dc_eff=kick['dc_eff_from_table']),
    chat14_reference=dict(windows='2-3 ... 4-5 (1-unit windows, |dev| < 1e-3)', rate_75pts=[1.65714, 1.65717], rate_134pts=[1.65718, 1.65719], grid='identical tanh grid (dzf 1.69e-4, dzc 8e-3, zfine 0.06, L 10, 1392 points)'))

out['conclusions'] = dict(
    pre_run=('The coarse old-chart pre-run (32 points per wall, dzc/dzf = 20, no damping) went numerically unstable before the physical '
             'light-cone approach: the shell lapse grew monotonically from T ~ 1 (1.03 at T = 3.4, 2.2 at T = 5.5) instead of decaying, the '
             'near-shell residual grew to O(1), and the run stopped at |B_b| > 40 (T = 5.68).  d b_b/dT stayed > -0.9 throughout, so the rule '
             'returned xc = inf.  In the finer main run the lapse decays (0.78 at T = 10), i.e. the growing lapse is a numerical artefact.'),
    spacings=('dzf = 1.5e-4: Remedy A (Newton projection of the initial warp on Z >= -9.9) diverged (residual 7e-7 -> 1.6e-2 -> 2.2e-6, '
              'warp shift 9e-3 instead of 1e-5), so the run started from constraint-violating data (H/H0 = 0.991, Weyl/H0^2 = -1.3e-2 at T = 0) '
              'and its lapse grew from T ~ 0.3 until |B_b| > 40 at H0 tau = 3.05.  dzf = 2e-4: the projection converged (slowly, linearly); the '
              'run was clean early, but a near-shell residual localised at z ~ -0.05 (the fine-to-coarse transition, dz 2e-4 -> 4e-3 over width '
              '0.02, where the damping also switches on) grew at a rate ~1 per unit T, exceeded 0.05 at T ~ 7.2 and made H and the Weyl term '
              'swing wildly after T ~ 10.5; the stop at H0 tau = 11.0 is that breakdown.  The two stops are two different numerical failures, not physics.'),
    grid=('At delta = 1e-3 the near-shell scale is 1/rho_b = 0.0127 and the local scale grows as 1/e^B(z) into the bulk (x3 at z = -0.05, '
          'x5 at z = -0.1).  The A1 tanh grid keeps 64-85 points per wall only for z > -0.04 and then drops to ~3-6 points per local scale '
          'at z ~ -0.1 (dzc/dzf = 20 vs 4 at delta = 0.1), which is where the residual grew.'),
    growth_calibration=('The A1 C2 run used static-reference tables with dx = 2.5e-4; at the shell they violate the phi junction by '
                        '-1.8e-10, 35x the physical seed mismatch -5e-12 (dc = 1e-8), plus an A-junction error -4.8e-11 with no physical counterpart.  '
                        'This extra, differently shaped kick produced a transient (sliding rate 1.82 -> 1.66 over H0 tau 1-5) and moved the '
                        'onset of nonlinearity earlier, so the fitted window [1.0, 5.2] mixed transient and mode.  The grid is not the cause: it is '
                        'the grid Chat 14 used for its 1.65714-1.65717 result.'),
    tuned_rate=('With the Y = 1 friction closure the tuned shell (d*) at delta = 1e-3 rolls slowly: |phi_b - phi_b(static)| grows at ~%.3f per '
                'unit T (the same model with Y = 0 grows at ~3.7, diag/D3_GRID_DEV.json: the roll is friction-limited, Y rho_b = 79 in Hubble '
                'units vs 7.8 at delta = 0.1).  From dc = 1e-2 (seed deviation 2.3e-3) phi_b ~ 0.5 needs ~%.0f T at that rate; the dc = 1e-4 '
                'seed needs ~%.0f T more.' % (mu, math.log(0.5/seed_dev_1e2)/mu, math.log(100)/mu)))

def clean(v):
    if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)): return [clean(x) for x in v]
    if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, np.bool_): return bool(v)
    return v
json.dump(clean(out), open(HERE/'diag/D2_A1_DELTA1E-3_DIAGNOSIS.json', 'w'), indent=1)
print(json.dumps(clean(out['conclusions']), indent=1))
print(json.dumps(clean(out['tuned_growth_rate']), indent=1))
print(json.dumps(clean(out['growth_calibration']['sliding_1_unit_windows_dev_lt_1em3'])))
for k, v in out['main_runs'].items(): print(k, v['projection'], v['first_T_near_shell_residual_gt_0p05'], v['z_of_domain_max_H_residual_T_gt_4'], v['grid_dz_at'])
