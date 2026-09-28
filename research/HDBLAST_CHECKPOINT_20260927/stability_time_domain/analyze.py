"""Aggregate spectra and time-domain runs into results/ANALYSIS.json (no new PDE solves).

Time-domain extraction:
  * original shell: least-squares slope of ln|f_b| over late windows (growth rate), 3+ resolutions;
  * +1 branch: (i) envelope rate from sliding-window RMS of f_b, (ii) matrix-pencil (Hua-Sarkar)
    complex exponents of f_b(t), (iii) compensated signal e^{1.5 t} f_b to expose any component
    decaying more slowly than the de Sitter continuum edge (a bound state) or growing.
Fully-discrete check: RK4 amplification |R(lambda dt)| for every dense eigenvalue with dt=0.4*hmin.
"""
import json, glob, math, os
from pathlib import Path
import numpy as np
import td_core as T

HERE = Path(__file__).resolve().parent
RES = HERE/'results'
GN_RATE = 1.6571936312453648          # Chat 9 independent Gaussian-normal shooting (frozen_input)
PRIOR_SPECTRAL = 1.657193663767       # registered-shell controls, finest eigenvalue (frozen_input)


def rk4_R(z):
    return 1 + z + z*z/2 + z**3/6 + z**4/24


def spectra_summary():
    out = {}
    for f in sorted(glob.glob(str(RES/'spectra'/'*.json'))):
        d = json.loads(Path(f).read_text()); name = Path(f).stem
        if d['mode'] != 'dense':
            out[name] = d; continue
        w = np.load(f.replace('.json', '_eigs.npz'))['w']
        dt = .4*d['hmin']*1.0
        amp = np.abs(rk4_R(w*dt)); rate = np.log(amp)/dt
        k = int(np.argmax(rate))
        top_non_gauge = [r for r in d['top_by_real_part'] if r['label'] != 'gauge']
        out[name] = dict(
            bg=d['bg'], hmin=d['hmin'], L=d['L'], stretch=d['stretch'], nodes=d['nodes'], s20_scale=d['s20_scale'],
            pot_mode=d.get('pot_mode'), matrix_vs_rhs=d['matrix_vs_complex_step_rhs_rel_error'],
            max_real_semi_discrete=d['max_real_all'],
            max_eigen_residual=float(np.max(np.load(f.replace('.json', '_eigs.npz'))['eigen_residual'])),
            rk4_max_effective_rate=float(rate[k]), rk4_argmax_eigenvalue=[float(w[k].real), float(w[k].imag)],
            n_on_dS_line=d['n_on_dS_line_upper_half'], max_dev_from_line=d['max_abs_real_plus_1p5_on_line'],
            real_axis=[dict(real=r['real'], label=r['label'], shell_phi_rel=r['shell_phi_rel'], far_weight=r['far_weight'],
                            gauge_template_residual=r['gauge_template_residual'], cH=r['constraint_H_l2_rel'],
                            cM=r['constraint_M_l2_rel']) for r in d['real_axis_eigenvalues']],
            shell_supported_above_line=[dict(real=r['real'], imag=r['imag'], label=r['label'], shell_phi_rel=r['shell_phi_rel'],
                                             cH=r['constraint_H_l2_rel']) for r in d['shell_supported_above_line']],
            lowest_shell_supported_line_mode=(d['lowest_frequency_shell_supported_line_modes'][0]
                                              if d['lowest_frequency_shell_supported_line_modes'] else None),
            top5_by_real=[dict(real=r['real'], imag=r['imag'], label=r['label'], shell_phi_rel=r['shell_phi_rel'])
                          for r in d['top_by_real_part'][:5]])
    return out


def growth_fit(t, y, t0, t1):
    m = (t >= t0) & (t <= t1)
    c = np.polyfit(t[m], np.log(np.abs(y[m])), 1)
    resid = np.log(np.abs(y[m])) - np.polyval(c, t[m])
    return float(c[0]), float(np.max(np.abs(resid)))


def matrix_pencil(y, dt, order, Lfrac=1/3):
    N = len(y); Lp = int(N*Lfrac)
    Y = np.array([y[i:i + Lp + 1] for i in range(N - Lp)])
    U, S, Vh = np.linalg.svd(Y, full_matrices=False)
    V = Vh.conj().T[:, :order]
    V1, V2 = V[:-1], V[1:]
    z = np.linalg.eigvals(np.linalg.pinv(V1)@V2)
    lam = np.log(z.astype(complex))/dt
    # amplitudes by least squares
    Z = np.vander(z, N, increasing=True).T
    amp = np.linalg.lstsq(Z, y.astype(complex), rcond=None)[0]
    rec = (Z@amp).real
    err = float(np.linalg.norm(rec - y)/np.linalg.norm(y))
    o = np.argsort(-np.abs(amp))
    return [dict(real=float(lam[k].real), imag=float(lam[k].imag), abs_amplitude=float(abs(amp[k]))) for k in o], \
        err, (S[:order + 5]/S[0]).tolist()


def window_rms(t, y, width):
    edges = np.arange(t[0], t[-1] - width + 1e-12, width)
    tc, r = [], []
    for e in edges:
        m = (t >= e) & (t < e + width); tc.append(e + width/2); r.append(np.sqrt(np.mean(y[m]**2)))
    return np.array(tc), np.array(r)


def slope(t, y, a0, a1):
    m = (t >= a0) & (t <= a1)
    return float(np.polyfit(t[m], np.log(y[m]), 1)[0])


def td_summary():
    out = {}
    for f in sorted(glob.glob(str(RES/'td'/'*.json'))):
        d = json.loads(Path(f).read_text()); name = Path(f).stem
        z = np.load(f.replace('.json', '.npz')); t = z['t']; fb = z['f_b']; hb = z['h_b']
        A_ = d['args']
        row = dict(bg=A_['bg'], hmin=A_['hmin'], L=A_['L'], linear=A_['linear'], eps=A_['eps'], tdet=A_['tdet'],
                   s20_scale=A_['s20_scale'], bumps=[A_['z1'], A_['w1'], A_['z2'], A_['w2']], tf=A_['tf'],
                   nodes=d['nodes'], dt=d['dt'], stop_reason=d['stop_reason'], runtime_seconds=d['runtime_seconds'],
                   initial_data=d['initial_data'], pot_mode=d.get('pot_mode'))
        diag = d['diagnostics']
        late = [r for r in diag if r['t'] >= .5] or diag
        row['constraints'] = dict(
            note=('linearised constraints, relative to the sum of |terms|' if A_['linear'] else
                  'full nonlinear constraints, relative to the sum of |perturbative terms|') + '; z>-0.8L, nodes>5',
            H_rel_l2_max_t_ge_0p5=max(r['H_rel_l2'] for r in late), M_rel_l2_max_t_ge_0p5=max(r['M_rel_l2'] for r in late),
            H_rel_l2_final=diag[-1]['H_rel_l2'], M_rel_l2_final=diag[-1]['M_rel_l2'],
            registered_normalised_H_max_final=diag[-1]['registered_H_max'],
            registered_normalised_M_max_final=diag[-1]['registered_M_max'])
        row['max_abs_a_final'] = diag[-1]['max_abs_a']; row['max_abs_f_final'] = diag[-1]['max_abs_f']
        te = t[-1]
        if A_['bg'] == 'original' or A_['s20_scale'] < 0:
            if A_['s20_scale'] < 0: wins = [(0.4*te, 0.7*te), (0.6*te, te)]; tp = 0.3*te
            else: wins = [(te - 1.5, te), (te - 1.0, te), (te - 2.0, te - 0.5)]; tp = 2.0
            row['growth_fits_ln_abs_f_b'] = []
            for a0, a1 in wins:
                r, res = growth_fit(t, fb, a0, a1)
                row['growth_fits_ln_abs_f_b'].append(dict(window=[a0, a1], rate=r, max_log_residual=res))
            m = t >= tp
            poles, err, sv = matrix_pencil(fb[m], t[1] - t[0], 8)
            real_poles = [p for p in poles if abs(p['imag']) < 1e-6]
            row['matrix_pencil_order8'] = dict(t_start=tp, relative_reconstruction_error=err, poles_by_amplitude=poles[:6],
                                               dominant_real_pole=real_poles[0] if real_poles else None)
            row['energy_rate_half_slope_lnE'] = slope(t, z['E'], wins[0][0], wins[0][1])/2
        else:
            E = z['E']; EN = z['E_near']; EX = z['E_X']
            row['scalar_energy'] = dict(
                amplitude_rate_from_lnE_1_to_tf=slope(t, E, 1.0, te)/2,
                amplitude_rate_from_lnE_0p5_to_3=slope(t, E, 0.5, 3.0)/2,
                amplitude_rate_from_lnE_3_to_tf=slope(t, E, 3.0, te)/2 if te > 3.5 else None,
                compensated_e3t_E_end_over_start=float(np.exp(3*te)*E[-1]/E[0]),
                E_X_max_relative_drift=float(np.max(np.abs(EX/EX[0] - 1))),
                E_X_final_relative_drift=float(EX[-1]/EX[0] - 1),
                E_near_shell_final_over_initial=float(EN[-1]/EN[0]) if EN[0] > 0 else None,
                E_near_shell_max_over_initial_after_t1=float(np.max(EN[t >= 1])/max(np.max(EN), 1e-300)))
            tc, r = window_rms(t, fb, .1)
            row['f_b_rms_windows_0p1'] = dict(t_center=tc.tolist(), rms=r.tolist())
            row['f_b_rms_peak'] = float(np.max(r)); row['f_b_rms_t_ge_1_max'] = float(np.max(r[tc >= 1]))
            row['f_b_rms_last'] = float(r[-1])
            tc3, r3 = window_rms(t, hb, .5)
            row['shell_hubble_rms_last_over_peak'] = float(r3[-1]/np.max(r3))
        row['final_state_constraints_by_region'] = region_constraints(d, z['final'])
        out[name] = row
    return out


_SOLVERS = {}


def region_constraints(d, v):
    A_ = d['args']; key = (A_['bg'], A_['hmin'], A_['L'], A_['stretch'], A_['tdet'], A_['s20_scale'])
    if key not in _SOLVERS:
        bg = T.plus_background(A_['tdet']) if A_['bg'] == 'plus' else T.original_background(A_['tdet'])
        s = T.GeneralSolver(bg, tdet=A_['tdet'], hmin=A_['hmin'], L=A_['L'], stretch=A_['stretch'])
        if A_['s20_scale'] != 1.: s.s20 *= A_['s20_scale']
        _SOLVERS[key] = s
    s = _SOLVERS[key]
    H, M, Hs, Ms = (T.linear_constraints if A_['linear'] else T.nonlinear_constraints)(s, v)
    idx = np.arange(s.n) > 5; res = {}
    gH = np.max(Hs[idx]); gM = np.max(Ms[idx])
    for nm, m in [('near_shell_z_gt_-0.2', s.z > -.2), ('mid_-1_to_-0.2', (s.z <= -.2) & (s.z > -1)),
                  ('far_-0.8L_to_-1', (s.z <= -1) & (s.z > -.8*s.L))]:
        m = m & idx
        res[nm] = dict(H_rel_l2=float(np.linalg.norm(H[m])/max(np.linalg.norm(Hs[m]), 1e-300)),
                       M_rel_l2=float(np.linalg.norm(M[m])/max(np.linalg.norm(Ms[m]), 1e-300)),
                       H_max_over_global_term_max=float(np.max(np.abs(H[m]))/gH),
                       M_max_over_global_term_max=float(np.max(np.abs(M[m]))/gM))
    return res


def pair_nonlinear(td):
    """Nonlinear-vs-linear differences in the shell scalar (normalised by eps) and eps-scaling."""
    res = {}
    for k, v in td.items():
        if v['linear'] or k.startswith('control'): continue
        lin = [kk for kk, vv in td.items() if vv['linear'] and vv['bg'] == v['bg'] and vv['hmin'] == v['hmin']
               and vv['L'] == v['L'] and vv['s20_scale'] == v['s20_scale'] and vv['bumps'] == v['bumps'] and vv['tdet'] == v['tdet']]
        if not lin: continue
        zn = np.load(RES/'td'/(k + '.npz')); zl = np.load(RES/'td'/(lin[0] + '.npz'))
        n = min(len(zn['t']), len(zl['t']))
        a = zn['f_b'][:n]/v['eps']; b = zl['f_b'][:n]/td[lin[0]]['eps']
        res[k] = dict(linear_partner=lin[0], max_abs_diff_fb_over_eps=float(np.max(np.abs(a - b))),
                      max_abs_fb_over_eps=float(np.max(np.abs(b))),
                      relative=float(np.max(np.abs(a - b))/np.max(np.abs(b))))
    return res


def main():
    spec = spectra_summary(); td = td_summary(); nl = pair_nonlinear(td)
    calib = []
    for k, v in sorted(td.items(), key=lambda kv: -kv[1]['hmin']):
        if v['bg'] != 'original': continue
        calib.append(dict(run=k, hmin=v['hmin'], linear=v['linear'], eps=v['eps'],
                          rate_window_last1p5=v['growth_fits_ln_abs_f_b'][0]['rate'],
                          rate_pencil=(v['matrix_pencil_order8']['dominant_real_pole'] or {}).get('real'),
                          diff_window_vs_GN=v['growth_fits_ln_abs_f_b'][0]['rate'] - GN_RATE,
                          diff_pencil_vs_GN=((v['matrix_pencil_order8']['dominant_real_pole'] or {}).get('real', np.nan) - GN_RATE)))
    out = dict(spectra=spec, time_domain=td, nonlinear_vs_linear=nl, calibration_time_domain=calib,
               references=dict(GN_rate=GN_RATE, prior_finest_spectral=PRIOR_SPECTRAL,
                               sources='frozen_input/CHAT9_INDEPENDENT_GN_RESULT.json; frozen_input/spectrum_reg_h05.json'))
    T.dump(RES/'ANALYSIS.json', out)
    for c in calib: print('CAL', c['run'], '%.9f %.9f  %.2e %.2e' % (c['rate_window_last1p5'], c['rate_pencil'] or np.nan, c['diff_window_vs_GN'], c['diff_pencil_vs_GN']))
    for k, v in spec.items():
        if 'max_real_semi_discrete' in v:
            print('SPEC', k, 'maxRe %.4g rk4max %.3g' % (v['max_real_semi_discrete'], v['rk4_max_effective_rate']),
                  'real:', [(round(r['real'], 6), r['label'][:5]) for r in v['real_axis']][:6])
    for k, v in td.items():
        if 'scalar_energy' in v:
            e = v['scalar_energy']
            print('TD', k, 'rate(E) %.5f [%.5f, %s] EXdrift %.1e fb peak %.2e late %.2e H %.1e' % (
                e['amplitude_rate_from_lnE_1_to_tf'], e['amplitude_rate_from_lnE_0p5_to_3'], e['amplitude_rate_from_lnE_3_to_tf'],
                e['E_X_max_relative_drift'], v['f_b_rms_peak'], v['f_b_rms_last'], v['constraints']['H_rel_l2_final']))
        else:
            print('TD', k, [round(g['rate'], 6) for g in v['growth_fits_ln_abs_f_b']], v['matrix_pencil_order8']['dominant_real_pole'])
    for k, v in nl.items(): print('NL', k, v)


if __name__ == '__main__':
    main()
