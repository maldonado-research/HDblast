#!/usr/bin/env python3
"""Collect the key numbers from PREHEATING_RESULTS.json and CONTROLS.json into SUMMARY.json and print the README tables.
Also tests the large-G closed form for the cutoff-limited radiation ratio:
    R_cut(G) ~ K / sqrt(G),  K = |v*|^{3/2} |phi_e - phi*| / (8 pi^3 a_e^3 Dmax^3 crit_hat)
which follows from n_k = exp(-pi k^2/q), nonrelativistic particles of mass G|phi_e-phi*| at the data end, and the
EFT condition max m_chi <= M5 (b <= (G Dmax)^-3)."""
import json, math
from pathlib import Path
HERE = Path(__file__).resolve().parent
import sys; sys.path.insert(0, str(HERE))
import ph_lib as L

R = json.loads((HERE/'PREHEATING_RESULTS.json').read_text())
C = json.loads((HERE/'CONTROLS.json').read_text())
screen = json.loads((L.CKPT/'matter/CROSSING_ELIGIBILITY_AND_TOY_BUDGET.json').read_text())
gmin = {}
for ps, v in screen['comparisons'].items():
    gmin[float(ps)] = v['Gmin_by_grid']['v2_t1e3_plus_c4_finecoarse']

crit = R['constants']['crit_hat']
rows = []; law = {}
for key, c in R['primary'].items():
    if 'failed' in c:
        rows.append({'key': key, 'failed': c['failed']}); continue
    p = c['post']; ps = c['phistar']; G = c['G']
    screened = G >= gmin.get(ps, float('inf'))
    row = {'phistar': ps, 'G': G, 'q': c['q'], 'screen_Gmin': gmin.get(ps), 'passes_local_screen': screened,
           'order_used': c['order_used'], 'window_right_truncated': c['window']['right_truncated'],
           'A_right_edge': c['window']['A_right_edge'],
           'rel_dev_from_instant_number': c['analytic']['rel_diff_number'],
           'N_comoving_at_astar': c['analytic']['N_comoving_numeric'],
           'rho_hat_end': p['rho_hat_end'], 'min_m_over_h_post': p['min_m_over_h_post_window'],
           'ratio_at_cutoff_full': p['ratio_to_crit_at_cutoff_full'], 'ratio_at_cutoff_post': p['ratio_to_crit_at_cutoff_post'],
           'lambda_c_needed_full_N0': p['lambda_c_needed_full']['0.0'], 'lambda_c_needed_full_N1': p['lambda_c_needed_full']['1.0'],
           'lambda_c_needed_post_N0': p['lambda_c_needed_post']['0.0'],
           'b_required_for_equality': p['b_required_for_equality'],
           'Rj_max_at_cutoff_full': p['Rj_max_at_cutoff_full'], 'Rj_max_at_equality_b': p.get('Rj_max_at_equality_b'),
           'tension_ratio_at_cutoff_full_end': p['tension_ratio_at_cutoff_full_end'],
           'cw_ratio_at_cutoff_full_end': p['cw_ratio_to_sigma_p_at_cutoff_full_end'],
           'ward_identity_max_rel_residual': p.get('ward_identity_max_rel_residual'),
           'decay_y1_max_ratio_at_cutoff_full': p['decay']['1.0']['max_ratio_to_crit_at_cutoff_full'],
           'decay_y0.1_max_ratio_at_cutoff_full': p['decay']['0.1']['max_ratio_to_crit_at_cutoff_full'],
           'decay_y0.01_max_ratio_at_cutoff_full': p['decay']['0.01']['max_ratio_to_crit_at_cutoff_full'],
           'decay_y1_fraction_converted_by_data_end': p['decay']['1.0']['rho_r_hat_at_data_end']/p['rho_hat_end'],
           'ratio_to_crossing_vacuum_at_cutoff_full': p['ratio_to_crossing_vacuum_at_cutoff_full']}
    # closed-form large-G law
    ve = abs(c['v_star']); ae = p['a_end_over_a_star']; me = p['m_end']; De = me/G
    Dmax = p['Lambda_hat_full_interval']/G
    K = ve**1.5*De/(8*math.pi**3*ae**3*Dmax**3*crit)
    row['K_closed_form'] = K; row['ratio_times_sqrtG'] = p['ratio_to_crit_at_cutoff_full']*math.sqrt(G)
    row['closed_form_rel_dev'] = row['ratio_times_sqrtG']/K - 1
    rows.append(row)

valid = [r for r in rows if 'failed' not in r and r['passes_local_screen']]
allrows = [r for r in rows if 'failed' not in r]
summ = {
    'status': 'NUMERICAL_CONDITIONAL',
    'crit_hat': crit, 'constants': R['constants'],
    'max_ratio_at_cutoff_full_all_scanned': max(r['ratio_at_cutoff_full'] for r in allrows),
    'max_ratio_at_cutoff_full_screen_passing': max(r['ratio_at_cutoff_full'] for r in valid),
    'max_ratio_at_cutoff_post_screen_passing': max(r['ratio_at_cutoff_post'] for r in valid),
    'min_lambda_c_needed_full_N0_screen_passing': min(r['lambda_c_needed_full_N0'] for r in valid),
    'min_lambda_c_needed_post_N0_screen_passing': min(r['lambda_c_needed_post_N0'] for r in valid),
    'max_Rj_at_cutoff_full_screen_passing': max(r['Rj_max_at_cutoff_full'] for r in valid),
    'range_Rj_at_equality_b_screen_passing': [min(r['Rj_max_at_equality_b'] for r in valid), max(r['Rj_max_at_equality_b'] for r in valid)],
    'max_decay_y1_ratio_at_cutoff_full_screen_passing': max(r['decay_y1_max_ratio_at_cutoff_full'] for r in valid),
    'closed_form_rel_dev_G_ge_1e4': {str(r['phistar']): r['closed_form_rel_dev'] for r in allrows if r['G'] == 1e6},
    'max_ward_residual': max(r['ward_identity_max_rel_residual'] or 0 for r in allrows),
    'calibration_pass': C['calibration_pass'],
    'grid_spread': R['grid_spread'],
    'rows': rows}
(HERE/'SUMMARY.json').write_text(json.dumps(summ, indent=1, allow_nan=False) + '\n')

def f(x, n=3):
    return '%.*g' % (n, x) if isinstance(x, (int, float)) else str(x)
print('| phi* | G | q/H0^2 | screen | dev. from instant N | rho_hat(end) | R at M5 cutoff | lambda_c for R=1 | lambda_c for N=1 | R_j at cutoff | R_j if R=1 | decay y=1: R_r at cutoff |')
print('|---|---|---|---|---|---|---|---|---|---|---|---|')
for r in allrows:
    print('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        r['phistar'], f(r['G']), f(r['q']), 'pass' if r['passes_local_screen'] else 'fail', f(r['rel_dev_from_instant_number'], 2),
        f(r['rho_hat_end']), f(r['ratio_at_cutoff_full'], 2), f(r['lambda_c_needed_full_N0']), f(r['lambda_c_needed_full_N1']),
        f(r['Rj_max_at_cutoff_full'], 2), f(r['Rj_max_at_equality_b'], 2), f(r['decay_y1_max_ratio_at_cutoff_full'], 2)))
print(json.dumps({k: v for k, v in summ.items() if k not in ('rows', 'constants', 'grid_spread')}, indent=1))
