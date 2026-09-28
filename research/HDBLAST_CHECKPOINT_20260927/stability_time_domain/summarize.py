"""Key results and pass/fail checks, read from results/ANALYSIS.json and the spectra JSONs.
Writes results/KEY_RESULTS.json.  No new PDE solves."""
import json, glob
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent; RES = HERE/'results'
A = json.loads((RES/'ANALYSIS.json').read_text())
GN = A['references']['GN_rate']; PRIOR = A['references']['prior_finest_spectral']
spec = A['spectra']; td = A['time_domain']


def shift_modes(name, maxres=1e-6):
    d = json.loads((RES/'spectra'/(name + '.json')).read_text()); out = {}
    for blk in d['shift_invert']:
        for m in blk.get('modes', []):
            if m['eigen_residual'] > maxres or m['imag'] < -1e-9: continue
            key = (round(m['real'], 6), round(m['imag'], 4))
            out.setdefault(key, m)
    return sorted(out.values(), key=lambda m: -m['real'])


key = {}; checks = {}
# ---------------- calibration
cal = [c for c in A['calibration_time_domain']]
key['calibration_time_domain'] = cal
orig_shift = {h: [m for m in shift_modes('shift_original_h' + h) if abs(m['real'] - 1.6572) < 1e-3][0]['real']
              for h in ['2e-4', '1e-4', '5e-5']}
key['calibration_semi_discrete_eigenvalue'] = orig_shift
fin = [c for c in cal if c['hmin'] == 5e-5 and c['linear']][0]
checks['calibration_td_pencil_finest_within_1e-6_of_GN'] = abs(fin['diff_pencil_vs_GN']) < 1e-6
checks['calibration_td_window_finest_within_2e-6_of_GN'] = abs(fin['diff_window_vs_GN']) < 2e-6
checks['calibration_eigen_finest_within_1e-7_of_prior'] = abs(orig_shift['5e-5'] - PRIOR) < 1e-7
checks['calibration_td_error_decreases_2e-4_to_1e-4'] = abs([c for c in cal if c['hmin'] == 1e-4 and c['linear']][0]['diff_pencil_vs_GN']) < \
    abs([c for c in cal if c['hmin'] == 2e-4 and c['linear']][0]['diff_pencil_vs_GN'])
nlcal = [c for c in cal if c['hmin'] == 5e-5 and not c['linear']][0]
checks['calibration_nonlinear_small_amplitude_agrees_1e-6'] = abs(nlcal['rate_pencil'] - fin['rate_pencil']) < 1e-6

# ---------------- +1 branch spectra
plus_dense = {k: v for k, v in spec.items() if k.startswith('plus_h') and 'max_real_semi_discrete' in v}
key['plus_dense'] = {k: dict(nodes=v['nodes'], L=v['L'], stretch=v['stretch'], pot_mode=v['pot_mode'],
                             max_real_semi_discrete=v['max_real_semi_discrete'], rk4_max_effective_rate=v['rk4_max_effective_rate'],
                             real_axis=[(r['real'], r['label'], r['shell_phi_rel'], r['far_weight']) for r in v['real_axis']],
                             n_on_dS_line=v['n_on_dS_line'], max_dev_from_line=v['max_dev_from_line'],
                             n_shell_supported_above_line=len(v['shell_supported_above_line']),
                             shell_supported_above_line_max_imag_min=min([abs(r['imag']) for r in v['shell_supported_above_line']] or [None]) if v['shell_supported_above_line'] else None,
                             shell_supported_above_line_max_shell_phi=max([r['shell_phi_rel'] for r in v['shell_supported_above_line']] or [0]))
                         for k, v in plus_dense.items()}
# real-axis eigenvalues above the line with non-negligible shell scalar: none expected
bad = []
for k, v in plus_dense.items():
    for r in v['real_axis']:
        if r['real'] > -1.5 and r['shell_phi_rel'] > 1e-6: bad.append((k, r['real'], r['shell_phi_rel']))
    for r in v['shell_supported_above_line']:
        # any shell-supported mode above the line that is NOT a grid-scale oscillation (|Im| > 1000)
        if abs(r['imag']) < 1000 and r['shell_phi_rel'] > 1e-6: bad.append((k, r['real'], r['imag'], r['shell_phi_rel']))
key['plus_dense_shell_supported_modes_above_line_excluding_grid_scale'] = bad
checks['plus_dense_no_shell_supported_mode_above_Re_-1.5_with_Im_lt_1000'] = len(bad) == 0
checks['plus_dense_rk4_fully_discrete_max_rate_is_gauge_1'] = all(abs(v['rk4_max_effective_rate'] - 1) < 1e-4 for v in plus_dense.values())
checks['plus_dense_on_line_to_1e-7'] = all(v['max_dev_from_line'] < 1e-7 for v in plus_dense.values())
plus_shift = {h: shift_modes('shift_plus_h' + h) for h in ['2e-4', '1e-4', '5e-5']}
key['plus_shift_invert_modes_Re_gt_-3'] = {h: [(m['real'], m['imag'], m['label'], m['shell_phi_rel'], m['far_weight'], m['gauge_template_residual'])
                                              for m in ms if m['real'] > -3] for h, ms in plus_shift.items()}
checks['plus_shift_all_modes_above_line_have_shell_phi_lt_1e-11'] = all(
    m['shell_phi_rel'] < 1e-11 for ms in plus_shift.values() for m in ms if m['real'] > -1.5 + 1e-6)
g = {h: [m for m in ms if abs(m['real'] - 1) < 1e-3][0] for h, ms in plus_shift.items()}
key['plus_gauge_mode_lambda_1'] = {h: dict(real=m['real'], gauge_template_residual=m['gauge_template_residual'], shell_phi_rel=m['shell_phi_rel'])
                                   for h, m in g.items()}
checks['plus_lambda1_matches_analytic_gauge_template_1e-5'] = all(m['gauge_template_residual'] < 1e-5 for m in g.values())
fb_modes = {k: [r for r in v['real_axis'] if 0 < r['real'] < 0.5][0]['real'] for k, v in spec.items()
            if k in ('plus_h4_L6', 'plus_h4_L8', 'plus_h4_L10')}
key['far_boundary_mode_vs_L'] = fb_modes
checks['far_boundary_mode_decreases_with_L'] = fb_modes['plus_h4_L6'] > fb_modes['plus_h4_L8'] > fb_modes['plus_h4_L10']
# Precision control
pp = [(a['real'], b['real']) for a, b in zip(spec['plus_h4_L6_phipoly']['real_axis'], spec['plus_h4_L6']['real_axis'])]
key['precision_control_real_axis_phi_poly_vs_eta_poly'] = pp
# tolerance 1e-7: dense-eigensolver eigen-residuals are ~2e-8 on this grid (spectra JSON)
checks['phi_poly_vs_eta_poly_real_eigs_above_line_agree_1e-7'] = all(abs(a - b) < 1e-7 for a, b in pp if b > -1.5)
key['precision_control_note'] = ('Real eigenvalues above the continuum line agree to <1e-7; the far-boundary partner modes '
                                 'below -3 differ by up to %.1e because the registered phi-polynomial loses U-derivative accuracy '
                                 'where |phi-1| < 1e-12 near the truncated outer boundary.' % max(abs(a - b) for a, b in pp))
lam1 = [m['real'] for m in shift_modes('shift_plus_h1e-4_L8') if abs(m['real'] - 1) < 1e-3][0]
key['gauge_lambda1_L6_vs_L8_h1e-4'] = dict(L6=g['1e-4']['real'], L8=lam1)
checks['gauge_lambda1_deviation_shrinks_with_L'] = abs(lam1 - 1) < abs(g['1e-4']['real'] - 1)

# ---------------- +1 branch time domain
rates = {}
for k, v in td.items():
    if 'scalar_energy' in v and k.startswith('plus'):
        rates[k] = v['scalar_energy']['amplitude_rate_from_lnE_1_to_tf']
key['plus_td_amplitude_rate_from_scalar_energy'] = rates
seq = [rates['plus_lin_h%s' % h] for h in ['2e-4', '1e-4', '5e-5']]
r1, r2 = seq[1] - seq[0], seq[2] - seq[1]
p = np.log2(abs(r1/r2)); extrap = seq[2] + r2/(2**p - 1)
key['plus_td_rate_richardson'] = dict(sequence_h_2e4_1e4_5e5=seq, observed_order=float(p), extrapolated=float(extrap))
seqB = [rates['plus_lin_bumpB_h%s' % h] for h in ['2e-4', '1e-4', '5e-5']]
r1, r2 = seqB[1] - seqB[0], seqB[2] - seqB[1]; pB = np.log2(abs(r1/r2)); extrapB = seqB[2] + r2/(2**pB - 1)
key['plus_td_rate_richardson_bumpB'] = dict(sequence=seqB, observed_order=float(pB), extrapolated=float(extrapB))
checks['plus_td_rate_converges_toward_-1.5'] = abs(seq[2] + 1.5) < abs(seq[1] + 1.5) < abs(seq[0] + 1.5) and abs(seq[2] + 1.5) < 2e-3
checks['plus_td_rate_bumpB_converges_toward_-1.5'] = abs(seqB[2] + 1.5) < abs(seqB[1] + 1.5) < abs(seqB[0] + 1.5) and abs(seqB[2] + 1.5) < 1e-3
checks['plus_td_no_growth_all_runs'] = all(r < -1.49 for r in rates.values())
checks['plus_td_L10_equals_L6_to_1e-3'] = abs(rates['plus_lin_L10_h1e-4'] - rates['plus_lin_h1e-4']) < 1e-3
fbl = {h: td['plus_lin_h' + h]['f_b_rms_last'] for h in ['4e-4', '2e-4', '1e-4', '5e-5']}
key['plus_td_late_shell_scalar_rms_last_window'] = fbl
key['plus_td_shell_scalar_peak_rms'] = td['plus_lin_h5e-5']['f_b_rms_peak']
checks['plus_td_late_shell_scalar_converges_to_zero'] = fbl['5e-5'] < fbl['1e-4'] < fbl['2e-4']
key['plus_td_EX_drift'] = {k: v['scalar_energy']['E_X_max_relative_drift'] for k, v in td.items() if k.startswith('plus_lin_h')}
key['plus_td_near_shell_constraints_final'] = {k: v['final_state_constraints_by_region']['near_shell_z_gt_-0.2'] for k, v in td.items() if k.startswith('plus_lin_h')}
nsH = [key['plus_td_near_shell_constraints_final']['plus_lin_h%s' % h]['H_rel_l2'] for h in ['4e-4', '2e-4', '1e-4', '5e-5']]
checks['plus_td_near_shell_H_constraint_converges'] = all(nsH[i + 1] < nsH[i] for i in range(3))
key['plus_td_far_region_constraints_final'] = {k: v['final_state_constraints_by_region']['far_-0.8L_to_-1'] for k, v in td.items() if k.startswith('plus_lin_h')}
nl = A['nonlinear_vs_linear']
key['nonlinear_vs_linear'] = nl
checks['plus_nonlinear_deviation_scales_linearly_with_eps'] = abs(nl['plus_nl_eps2e-6_h1e-4']['relative']/nl['plus_nl_eps1e-6_h1e-4']['relative'] - 2) < 0.05

# ---------------- controls
cs = {}
for h in ['1e-4', '5e-5']:
    cs['shift_h' + h] = [m['real'] for m in shift_modes('control_shift_plus_s20flip_h' + h) if m['real'] > 100][0]
cs['dense_h4e-4'] = spec['control_plus_s20flip_h4']['max_real_semi_discrete']
for h in ['2e-4', '1e-4']:
    cs['td_pencil_h' + h] = td['control_plus_s20flip_lin_h' + h]['matrix_pencil_order8']['dominant_real_pole']['real']
key['control_sigma2_flipped_growth_rate'] = cs
checks['control_s20flip_detected_unstable_by_both_methods'] = all(v > 100 for v in cs.values()) and abs(cs['td_pencil_h1e-4'] - cs['shift_h1e-4']) < 1e-4
pr = [r['real'] for r in spec['control_plus_s20flip_h4']['real_axis']]
checks['control_s20flip_dS_pairing_lambda_plus_partner_eq_-3'] = abs(max(pr) + min(pr) + 3) < 1e-6
checks['control_s20zero_stable'] = spec['control_plus_s20zero_h4']['max_real_semi_discrete'] < 1.0 + 1e-4 and td['control_plus_s20zero_lin_h1e-4']['scalar_energy']['amplitude_rate_from_lnE_1_to_tf'] < -1.49
key['perturbed_parameter_tdet_3e-3'] = dict(td_rate=rates.get('plus_tdet3e-3_lin_h1e-4'),
                                           modes_above_line=[(m['real'], m['imag'], m['label'], m['shell_phi_rel']) for m in shift_modes('shift_plus_tdet3e-3_h1e-4') if m['real'] > -1.5 + 1e-6])
key['perturbed_parameter_tdet_1e-2'] = dict(modes_above_line=[(m['real'], m['imag'], m['label'], m['shell_phi_rel']) for m in shift_modes('shift_plus_tdet1e-2_h1e-4') if m['real'] > -1.5 + 1e-6])
checks['perturbed_detunings_no_shell_supported_mode_above_line'] = all(x[3] < 1e-11 for kk in ['perturbed_parameter_tdet_3e-3', 'perturbed_parameter_tdet_1e-2'] for x in key[kk]['modes_above_line'])

out = dict(passed=all(checks.values()), n_checks=len(checks), checks=checks, key_results=key)
(RES/'KEY_RESULTS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps(dict(passed=out['passed'], failures=[k for k, v in checks.items() if not v]), indent=1))
print(json.dumps(dict(cal=[(c['run'], c['rate_pencil'], c['diff_pencil_vs_GN']) for c in cal], eig=orig_shift,
                      rich=key['plus_td_rate_richardson'], richB=key['plus_td_rate_richardson_bumpB'], ctrl=cs, fb=fb_modes,
                      gauge=key['plus_gauge_mode_lambda_1'], fbl=fbl, nsH=nsH), indent=1, default=float))
