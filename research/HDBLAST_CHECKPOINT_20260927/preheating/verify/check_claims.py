#!/usr/bin/env python3
"""Audit script 1: recompute every headline number of the preheating README from the saved JSON files
(PREHEATING_RESULTS.json, CONTROLS.json, CORRECTION_FIT.json, SUMMARY.json) and from the checkpoint inputs,
independently of build_summary.py.  Read-only on the audited folder; writes verify/CHECK_CLAIMS.json.
Optional argument: a directory holding re-run JSON files (default: the audited folder)."""
import json, math, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
SRC = HERE.parent if len(sys.argv) < 2 else Path(sys.argv[1])
CKPT = Path('/home/user/unified-theory-maldonado/new-files/latest-work/'
            'HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922')

R = json.loads((SRC/'PREHEATING_RESULTS.json').read_text())
C = json.loads((SRC/'CONTROLS.json').read_text())
F = json.loads((SRC/'CORRECTION_FIT.json').read_text())
S = json.loads((SRC/'SUMMARY.json').read_text())
out = {'source_dir': str(SRC)}

# ---- constants from the checkpoint, recomputed independently ------------------------------------------
br = json.loads((CKPT/'static_branch/PLUS_BRANCH_RESULTS.json').read_text())
row = next(r for r in br['rows'] if abs(r['delta'] - 1e-3) < 1e-14)['refined']
c_reg = br['c'] if isinstance(br['c'], float) else 0.5975949350280132
W = lambda p: 1 - p + p**3/3
sig = lambda p: 2*W(p) + 1e-3*(1 + 0.5975949350280132*p)
rho_b0 = R['constants']['rho_b0']
sf = sig(row['phi_b'])*rho_b0; hv = math.sqrt(row['H2'])*rho_b0
crit = math.sqrt(sf**2 + 36*hv**2) - sf
out['constants'] = {'c_in_branch_json': br['c'], 'c_used_by_ph_lib': 2/1.0357712571566784 - 4/3,
                    'sigma_f_over_H0': sf, 'h_vac': hv, 'crit_hat': crit, 'crit_hat_reported': R['constants']['crit_hat'],
                    'crit_rel_diff': crit/R['constants']['crit_hat'] - 1,
                    'crit_lowenergy_3H2_over_kappa4^2 (=18 h^2/sigma)': 18*hv**2/sf,
                    'rho_b0_reported': rho_b0, 'rho_b_branch': row['rho_b'],
                    'H_vac_over_H0_from_rho_b_ratio': rho_b0/row['rho_b']}

# ---- screen thresholds ------------------------------------------------------------------------------------
scr = json.loads((CKPT/'matter/CROSSING_ELIGIBILITY_AND_TOY_BUDGET.json').read_text())
gmin = {float(k): v['Gmin_by_grid']['v2_t1e3_plus_c4_finecoarse'] for k, v in scr['comparisons'].items()}
out['Gmin_finecoarse'] = gmin

# ---- scan table recomputed ----------------------------------------------------------------------------------
rows = []
for key, c in R['primary'].items():
    if 'failed' in c:
        rows.append({'key': key, 'failed': True}); continue
    p = c['post']; ps = c['phistar']; G = c['G']
    Lf, Lp = p['Lambda_hat_full_interval'], p['Lambda_hat_post_crossing']
    Rf = p['rho_hat_end']/crit/Lf**3
    Rp = p['rho_hat_end']/crit/Lp**3
    lam = Lf*(crit/p['rho_hat_end'])**(1/3)
    lamp = Lp*(crit/p['rho_hat_end'])**(1/3)
    # max over the stored post-window profile (every 40th point + end) of b*rho/crit at the full cutoff
    Rprof = max(q['rho'] for q in p['profile'])/crit/Lf**3
    rows.append({'key': key, 'phistar': ps, 'G': G, 'screen': G >= gmin.get(ps, 1e99), 'q': c['q'],
                 'R_full': Rf, 'R_post': Rp, 'R_full_reported': p['ratio_to_crit_at_cutoff_full'],
                 'lambda_c_N0': lam, 'lambda_c_N1': lam*math.exp(4/3), 'lambda_c_post_N0': lamp,
                 'R_cross': p['ratio_to_crossing_vacuum_at_cutoff_full'],
                 'tension_ratio': p['tension_ratio_at_cutoff_full_end'],
                 'Rj_cut': p['Rj_max_at_cutoff_full'], 'Rj_eq': p.get('Rj_max_at_equality_b'),
                 'cw': p['cw_ratio_to_sigma_p_at_cutoff_full_end'],
                 'ward': p.get('ward_identity_max_rel_residual'),
                 'decay_max': {y: p['decay'][y]['max_ratio_to_crit_at_cutoff_full'] for y in p['decay']},
                 'decay_y1_conv_by_end': p['decay']['1.0']['rho_r_hat_at_data_end']/p['rho_hat_end'],
                 'decay_during_prod_y1': p['decay']['1.0']['decay_during_production_estimate'],
                 'R_max_over_stored_profile_at_full_cutoff': Rprof,
                 'wronskian': c['wronskian_dev'], 'order_used': c['order_used'],
                 'A_right_edge': c['window']['A_right_edge'], 'right_truncated': c['window']['right_truncated'],
                 'rel_dev_instant': c['analytic']['rel_diff_number'], 'D_meas': c['analytic']['rel_diff_number']*c['q']})
ok = [r for r in rows if not r.get('failed')]
scr_rows = [r for r in ok if r['screen']]
def arg(rs, k, f=max):
    r = f(rs, key=lambda x: x[k]); return {'value': r[k], 'at': r['key']}
out['aggregates'] = {
    'n_rows': len(rows), 'n_failed': len(rows) - len(ok), 'n_screen_passing': len(scr_rows),
    'max_R_full_screen': arg(scr_rows, 'R_full'), 'max_R_full_all': arg(ok, 'R_full'),
    'max_R_post_screen': arg(scr_rows, 'R_post'), 'min_lambda_N0_screen': arg(scr_rows, 'lambda_c_N0', min),
    'min_lambda_post_N0_screen': arg(scr_rows, 'lambda_c_post_N0', min),
    'max_R_cross_screen': arg(scr_rows, 'R_cross'), 'max_R_cross_all': arg(ok, 'R_cross'),
    'max_tension_ratio_all': arg(ok, 'tension_ratio'),
    'max_Rj_cut_screen': arg(scr_rows, 'Rj_cut'),
    'Rj_eq_range_screen': [min(r['Rj_eq'] for r in scr_rows), max(r['Rj_eq'] for r in scr_rows)],
    'max_decay_any_y_screen': max((v, r['key'], y) for r in scr_rows for y, v in r['decay_max'].items()),
    'max_decay_y1_screen': arg([dict(r, d1=r['decay_max']['1.0']) for r in scr_rows], 'd1'),
    'max_ward_all': arg([r for r in ok if r['ward'] is not None], 'ward'),
    'max_wronskian_G>=100': arg([r for r in ok if r['G'] >= 100], 'wronskian'),
    'max_wronskian_all': arg(ok, 'wronskian'),
    'max_R_over_profile_screen': arg(scr_rows, 'R_max_over_stored_profile_at_full_cutoff'),
    'cw_phi0.5_G300_1000': [r['cw'] for r in ok if r['phistar'] == 0.5 and r['G'] in (300.0, 1000.0)],
    'cw_phi0.5_G1e6': [r['cw'] for r in ok if r['phistar'] == 0.5 and r['G'] == 1e6],
    'cw_ratio_screen_range': [min(r['cw'] for r in scr_rows), max(r['cw'] for r in scr_rows)],
    'n_screen_rows_with_cw_gt_1': sum(r['cw'] > 1 for r in scr_rows),
    'decay_y1_conversion_phi0.5_0.75_G300_1000': {r['key']: r['decay_y1_conv_by_end'] for r in ok
                                                  if r['phistar'] in (0.5, 0.75) and r['G'] in (300.0, 1000.0)},
    'rows_order0_fallback': [r['key'] for r in ok if r['order_used'] == 0],
    'rows_right_truncated_A_gt_1e-2': [r['key'] for r in ok if r['right_truncated'] and r['A_right_edge'] > 1e-2],
}
# SUMMARY.json consistency
out['summary_json_consistency'] = {
    k: {'summary': S[k], 'recomputed': v} for k, v in [
        ('max_ratio_at_cutoff_full_screen_passing', out['aggregates']['max_R_full_screen']['value']),
        ('max_ratio_at_cutoff_full_all_scanned', out['aggregates']['max_R_full_all']['value']),
        ('max_ratio_at_cutoff_post_screen_passing', out['aggregates']['max_R_post_screen']['value']),
        ('min_lambda_c_needed_full_N0_screen_passing', out['aggregates']['min_lambda_N0_screen']['value']),
        ('min_lambda_c_needed_post_N0_screen_passing', out['aggregates']['min_lambda_post_N0_screen']['value']),
        ('max_Rj_at_cutoff_full_screen_passing', out['aggregates']['max_Rj_cut_screen']['value']),
        ('max_decay_y1_ratio_at_cutoff_full_screen_passing', out['aggregates']['max_decay_y1_screen']['value'])]}

# ---- README table (decay column) vs JSON ------------------------------------------------------------------------
readme_decay = {'0.25|100': 5.9e-06, '0.25|1000': 5.2e-06, '0.25|1e+06': 1.7e-06, '0.5|10': 9.0e-05, '0.5|100': 7.4e-05,
                '0.5|300': 7.5e-05, '0.5|1000': 6.2e-05, '0.5|10000': 1.1e-04, '0.5|1e+06': 1.3e-05, '0.75|100': 1.6e-05,
                '0.75|300': 2.0e-05, '0.75|1000': 2.0e-05, '0.75|10000': 2.2e-05, '0.75|1e+06': 2.4e-06,
                '0.9|3000': 2.8e-06, '0.9|1e+06': 3.0e-07}
out['readme_decay_column_vs_json'] = {k: {'readme': v, 'json': next(r['decay_max']['1.0'] for r in ok if r['key'] == k)}
                                      for k, v in readme_decay.items()}
out['readme_decay_column_mismatch_count'] = sum(abs(d['json']/d['readme'] - 1) > 0.1 for d in out['readme_decay_column_vs_json'].values())

# ---- K/sqrt(G) law, recomputed ---------------------------------------------------------------------------
law = {}
for ps in (0.25, 0.5, 0.75, 0.9):
    c = R['primary']['%s|1e+06' % ps]; p = c['post']
    De = p['m_end']/c['G']; Dmax = p['Lambda_hat_full_interval']/c['G']
    K = abs(c['v_star'])**1.5*De/(8*math.pi**3*p['a_end_over_a_star']**3*Dmax**3*crit)
    law[str(ps)] = {'K': K, 'R_sqrtG_G1e6': p['ratio_to_crit_at_cutoff_full']*math.sqrt(c['G']),
                    'rel_dev': p['ratio_to_crit_at_cutoff_full']*math.sqrt(c['G'])/K - 1}
out['K_law'] = law

# ---- D values: measured (primary scan, largest reliable G) vs closed form -----------------------------------------
pred = F['trajectory_prediction']
Dm = {}
for ps in ('0.25', '0.5', '0.75', '0.9'):
    meas = {G: R['primary']['%s|%s' % (ps, G)]['analytic']['rel_diff_number']*R['primary']['%s|%s' % (ps, G)]['q']
            for G in ('10000', '100000', '1e+06')}
    cf = pred[ps]['D_predicted_by_closed_form']
    Dm[ps] = {'D_measured': meas, 'D_closed_form': cf,
              'rel_dev_closed_vs_meas': {G: cf/v - 1 for G, v in meas.items()}}
out['D_trajectory'] = Dm
asym = C['asymptotic_scaling']
out['G_times_rel_spread'] = {ps: {'min': min(a['G_times_rel_diff']), 'max': max(a['G_times_rel_diff']),
                                  'rel_spread_1e3_to_1e6': max(a['G_times_rel_diff'])/min(a['G_times_rel_diff']) - 1,
                                  'rel_spread_1e4_to_1e6': max(a['G_times_rel_diff'][2:])/min(a['G_times_rel_diff'][2:]) - 1,
                                  'slope': a['loglog_slope']} for ps, a in asym.items()}
# calibration numbers
cal = C['calibration_flat_linear']
out['calibration'] = {'max_rel_nk': max(v['max_rel_diff_nk_where_exact_n>1e-6'] for v in cal.values()),
                      'max_abs_relN': max(abs(v['rel_diff_number_correct_formula']) for v in cal.values()),
                      'wrong_half_q_mu0=0': [v['rel_diff_number_wrong_half_q'] for v in cal.values() if v['mu0'] == 0],
                      'wrong_double_mu0=0': [v['rel_diff_number_wrong_double_exponent'] for v in cal.values() if v['mu0'] == 0],
                      'analytic_wrong_half_q': 2**-1.5 - 1, 'analytic_wrong_double': 2**1.5 - 1,
                      'note': 'the wrong-formula numbers are fixed by Gaussian algebra: N ratios 2^{3/2} and 2^{-3/2}'}
out['fit'] = {'closed_resid': F['closed_form_max_abs_residual_on_toys'], 'wrong_resid': F['wrong_exponent_only_max_abs_residual_on_toys'],
              'fit_resid': F['fit_max_abs_residual'], 'closed_minus_fit': F['closed_form_minus_fit'],
              'n_toys_in_fit': len(F['runs']) - 1}
(HERE/'CHECK_CLAIMS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps({k: out[k] for k in ('constants', 'aggregates', 'summary_json_consistency', 'readme_decay_column_mismatch_count',
                                      'K_law', 'D_trajectory', 'G_times_rel_spread', 'calibration', 'fit')}, indent=1, default=float))
