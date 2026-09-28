#!/usr/bin/env python3
"""Explicit chi mode-function production on the archived Chat14 shell trajectories, with energy comparison,
one-shot backreaction estimate and an optional perturbative decay ledger.  Writes PREHEATING_RESULTS.json.

All rates in units of H0 = 1/rho_b0 (initial static shell of the archived runs).  Gravitational strength enters only
through b = kappa5^2 H0^3 = (H0/M5)^3 with M5 = kappa5^(-2/3); every density comparison is linear in b.
The bulk-coupled feedback (all three modified junctions) is NOT solved here: the shell trajectory is held fixed.
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ph_lib as L

PRIMARY = 'v2_t1e3_plus_c4_finecoarse'
GRID_CONTROLS = ['v2_t1e3_plus_c4', 'v2_t1e3_plus_c4_fine']
MINUS = 'v2_t1e3_minus_c4_finecoarse'
PHISTARS = [0.25, 0.5, 0.75, 0.9]
GS = [10.0, 30.0, 100.0, 300.0, 1e3, 3e3, 1e4, 1e5, 1e6]
YUK = [0.01, 0.1, 1.0]
N_EFOLD_TARGET = [0.0, 1.0, 5.0]


def constants(traj):
    rho_b0 = traj.meta['rho_b']
    br = L.load_branch()
    sig_f = L.sigma(br['phi_b'])*rho_b0            # sigma_f / H0
    h_vac = math.sqrt(br['H2'])*rho_b0             # H_vac / H0
    crit = math.sqrt(sig_f**2 + 36*h_vac**2) - sig_f   # kappa5^2 rho_crit / H0   (exact benchmark threshold)
    crit_dec = (math.sqrt(sig_f**2 + 108*h_vac**2) - sig_f)/3
    return dict(rho_b0=rho_b0, sigma_f_over_H0=sig_f, h_vac=h_vac, crit_hat=crit, crit_decel_hat=crit_dec,
                branch_phi_b=br['phi_b'], branch_H2=br['H2'])


def post_production(traj, case, C, n_grid=400, extrapolate=True):
    """Densities after the nonadiabatic window, energy comparisons, backreaction and decay ledger."""
    G, ps, mu0 = case['G'], case['phistar'], case['mu0']
    kaps, wts, nk, lna_star = case['_kaps'], case['_wts'], case['_n'], case['_lna_star']
    # densities are evaluated with the final occupation from the point where the kappa=0 adiabaticity drops below 0.05
    # (between that point and the tol=1e-3 window end the occupation is not yet frozen to the quoted accuracy)
    s_b = min(case['window']['s_particle_A0.05'], case['window']['s_b']); s_end = traj.s_max
    s_grid = np.linspace(s_b, s_end, n_grid) if s_end > s_b + 1e-9 else np.array([s_end])
    rows = [L.densities(traj, s, kaps, wts, nk, G, ps, mu0, lna_star) for s in s_grid]
    for r in rows:
        ph, hh = traj.bg(r['s']); r['h'] = float(hh); r['v'] = float(traj.bg(r['s'], 1)[0])
        r['sigma_hat'] = L.sigma(ph)*C['rho_b0']; r['sigma_p_hat'] = L.sigma_p(ph)*C['rho_b0']
        r['m_over_h'] = r['m']/abs(hh)
    end = rows[-1]
    # EFT scale: largest chi mass over the whole analysed interval (pre-crossing included) and sqrt(q)
    s_all = np.linspace(traj.s_min, s_end, 4000)
    phis = traj.bg(s_all)[:, 0]
    m_all = np.sqrt(mu0**2 + G**2*(phis - ps)**2)
    post = s_all >= case['s_star']
    Lam_full = float(max(np.max(m_all), math.sqrt(case['q'])))
    Lam_post = float(max(np.max(m_all[post]), math.sqrt(case['q'])))
    b_cut_full = Lam_full**-3; b_cut_post = Lam_post**-3
    crit = C['crit_hat']
    rho_end = end['rho']
    # crossing-time local comparison: exact benchmark with sigma(phi*) and H(s*)
    sig_star = L.sigma(ps)*C['rho_b0']; hs = case['h_star']
    crit_star = math.sqrt(sig_star**2 + 36*hs**2) - sig_star
    out = {
        'rho_hat_end': rho_end, 'n_hat_end': end['n'], 'a_end_over_a_star': end['a'], 'm_end': end['m'],
        'p_over_rho_end': end['p']/rho_end if rho_end > 0 else None,
        'min_m_over_h_post_window': float(min(r['m_over_h'] for r in rows)),
        'Lambda_hat_full_interval': Lam_full, 'Lambda_hat_post_crossing': Lam_post,
        'b_cutoff_full': b_cut_full, 'b_cutoff_post': b_cut_post,
        # ratio to the residual-vacuum radiation threshold if ALL chi energy at s_end became radiation instantly
        'coef_rho_over_crit_per_b': rho_end/crit,
        'ratio_to_crit_at_cutoff_full': b_cut_full*rho_end/crit,
        'ratio_to_crit_at_cutoff_post': b_cut_post*rho_end/crit,
        'lambda_c_needed_full': {str(N): Lam_full*(math.exp(4*N)*crit/rho_end)**(1/3) for N in N_EFOLD_TARGET} if rho_end > 0 else None,
        'lambda_c_needed_post': {str(N): Lam_post*(math.exp(4*N)*crit/rho_end)**(1/3) for N in N_EFOLD_TARGET} if rho_end > 0 else None,
        'b_required_for_equality': crit/rho_end if rho_end > 0 else None,
        'crossing_local_threshold_hat': crit_star,
        'rho_hat_just_after_window': rows[0]['rho'],
        'ratio_to_crossing_vacuum_at_cutoff_full': b_cut_full*rows[0]['rho']/crit_star,
        'tension_ratio_at_cutoff_full_end': b_cut_full*rho_end/end['sigma_hat'],
    }
    # one-shot backreaction on the scalar junction: R_j = kappa5^2 j / |sigma'| along the post-window interval
    Rj = np.array([abs(r['j'])/max(abs(r['sigma_p_hat']), 1e-300) for r in rows])     # per unit b
    out['Rj_per_b_max'] = float(np.max(Rj)); out['Rj_per_b_end'] = float(Rj[-1])
    out['Rj_max_at_cutoff_full'] = float(np.max(Rj))*b_cut_full
    if rho_end > 0:
        out['Rj_max_at_equality_b'] = float(np.max(Rj))*crit/rho_end
    # work done on chi by the rolling scalar, int j v ds, versus comoving energy (per a_*^3); checks the Ward identity
    sg = np.array([r['s'] for r in rows]); a3 = np.array([r['a']**3 for r in rows])
    jv = np.array([r['j']*r['v'] for r in rows]); rh = np.array([r['rho'] for r in rows]); pp = np.array([r['p'] for r in rows])
    hh = np.array([r['h'] for r in rows])
    if len(sg) > 5:
        drho = np.gradient(rh, sg, edge_order=2)
        ward = drho + 3*hh*(rh + pp) - jv
        scale = np.abs(drho) + 3*np.abs(hh)*(rh + pp) + np.abs(jv)
        out['ward_identity_max_rel_residual'] = float(np.max(np.abs(ward[2:-2])/scale[2:-2]))
    # Coleman-Weinberg size (not included in j; assumed absorbed into the renormalized sigma(phi))
    d_end = end['phi'] - ps
    out['cw_j_hat_end_nolog'] = G**2*abs(d_end)*end['m']**2/(16*math.pi**2)
    out['cw_ratio_to_sigma_p_at_cutoff_full_end'] = b_cut_full*out['cw_j_hat_end_nolog']/abs(end['sigma_p_hat'])
    # decay ledger chi -> psi psibar, Gamma_rest = y^2 m/(8 pi); per-mode dilated rate y^2 m^2/(8 pi omega)
    dec = {}
    for y in YUK:
        S = np.ones_like(kaps); rad_a4 = 0.0; best = 0.0; best_s = None
        prev = None
        for r in rows:
            a = r['a']; m2 = r['m']**2; om = np.sqrt((kaps/a)**2 + m2)
            pref = 1/(2*math.pi**2*a**3)
            Q = pref*np.sum(wts*kaps**2*nk*S*y*y*m2/(8*math.pi))
            if prev is not None:
                ds = r['s'] - prev['s']
                rad_a4 += 0.5*ds*(Q*a**4 + prev['Q']*prev['a']**4)
                gam = y*y*m2/(8*math.pi*om); gam0 = prev['gam']
                S = S*np.exp(-0.5*ds*(gam + gam0))
            r_rad = rad_a4/a**4
            if r_rad > best: best, best_s = r_rad, r['s']
            prev = {'s': r['s'], 'Q': Q, 'a': a, 'gam': y*y*m2/(8*math.pi*om)}
        rho_chi_left = float(np.sum(wts*kaps**2*nk*S*np.sqrt((kaps/end['a'])**2 + end['m']**2))/(2*math.pi**2*end['a']**3))
        rad_end = rad_a4/end['a']**4
        # conditional extrapolation beyond the archived data: phi frozen at phi_end, h = h_vac constant
        best_ext, s_ext = best, best_s
        if extrapolate:
            hv = C['h_vac']; m = end['m']; gamma_r = y*y*m/(8*math.pi)
            Nchi = rho_chi_left/m*end['a']**3      # comoving number density (nonrelativistic by then)
            ds_max = min(60.0/max(gamma_r, 1e-30), 40.0/hv)
            xs = np.linspace(0, ds_max, 20001)
            aext = end['a']*np.exp(hv*xs)
            src = Nchi*np.exp(-gamma_r*xs)*gamma_r*m*aext      # d(rho_r a^4)/ds = Gamma rho_chi a^4
            ra4 = rad_a4 + np.concatenate([[0.0], np.cumsum(0.5*np.diff(xs)*(src[1:] + src[:-1]))])
            rr = ra4/aext**4
            i = int(np.argmax(rr))
            if rr[i] > best_ext: best_ext, s_ext = float(rr[i]), float(s_end + xs[i])
        dec[str(y)] = {'Gamma_rest_over_H0_end': y*y*end['m']/(8*math.pi),
                       'decay_during_production_estimate': y*y/(8*math.pi),
                       'rho_r_hat_at_data_end': rad_end, 'rho_chi_hat_remaining_at_data_end': rho_chi_left,
                       'rho_r_hat_max_incl_extrapolation': best_ext, 's_of_max': s_ext,
                       'max_ratio_to_crit_at_cutoff_full': b_cut_full*best_ext/crit,
                       'Nmax_at_cutoff_full': 0.25*math.log(b_cut_full*best_ext/crit) if best_ext > 0 else None}
    out['decay'] = dec
    out['profile'] = [{k: r[k] for k in ('s', 'a', 'phi', 'm', 'n', 'rho', 'p', 'chi2', 'j', 'h', 'm_over_h')}
                      for r in rows[::40]] + [{k: end[k] for k in ('s', 'a', 'phi', 'm', 'n', 'rho', 'p', 'chi2', 'j', 'h', 'm_over_h')}]
    return out


def strip(case):
    return {k: v for k, v in case.items() if not k.startswith('_')}


def main():
    t0 = time.time()
    res = {'status': 'NUMERICAL_CONDITIONAL', 'archive_sha256': L.sha256_file(L.ARCHIVE),
           'branch_sha256': L.sha256_file(L.BRANCH),
           'units': 'rates in H0 = 1/rho_b0; densities in H0^4 (number H0^3); b = kappa5^2 H0^3 = (H0/M5)^3',
           'scope': 'Fixed archived shell trajectory; chi on the shell in the adiabatic vacuum before crossing; '
                    'no bulk feedback, no nonlinear chi dynamics, no rescattering, no thermalisation',
           'settings': {'tol_adiabatic_window': 1e-3, 'wkb_order': 2, 'rtol': 1e-10, 'nk': 48, 'mu0': 0.0,
                        't_cut_coordinate': L.T_CAUSAL}, 'primary': {}, 'grid_controls': {}, 'minus_branch': {}}
    tr = L.Trajectory(PRIMARY)
    C = constants(tr); res['constants'] = C
    res['trajectory_checks'] = {PRIMARY: {'s_range': [tr.s_min, tr.s_max], 'lna_crosscheck_max_abs': tr.lna_crosscheck,
                                          'spline_node_residual': tr.phi_fit_residual}}
    for ps in PHISTARS:
        for G in GS:
            t1 = time.time()
            try:
              case = L.run_case(tr, G, ps, history_kappa=[0.3*math.sqrt(G*0.5/math.pi)] if (G == 1e4 and ps in (0.5, 0.75)) else None,
                              n_hist=120)
              case_out = strip(case)
              case_out['post'] = post_production(tr, case, C)
            except Exception as e:
              res['primary']['%s|%g' % (ps, G)] = {'failed': repr(e)}
              print('FAILED', ps, G, repr(e), flush=True)
              continue
            case_out['runtime_s'] = round(time.time() - t1, 2)
            res['primary']['%s|%g' % (ps, G)] = case_out
            print('primary', ps, G, 'q=%.3g' % case['q'], 'relN=%.3e' % case['analytic']['rel_diff_number'],
                  'rho_end=%.4g' % case_out['post']['rho_hat_end'],
                  'ratio@cut=%.3e' % case_out['post']['ratio_to_crit_at_cutoff_full'], '%.1fs' % (time.time() - t1), flush=True)
    for tag in GRID_CONTROLS:
        tg = L.Trajectory(tag)
        res['trajectory_checks'][tag] = {'s_range': [tg.s_min, tg.s_max], 'lna_crosscheck_max_abs': tg.lna_crosscheck}
        for ps in (0.5, 0.75):
            for G in (100.0, 1e3, 1e4, 1e5):
                case = L.run_case(tg, G, ps)
                co = strip(case); co['post'] = post_production(tg, case, constants(tg), extrapolate=False)
                for k in ('kappa', 'weights', 'n_k_order0_basis'): co.pop(k)
                res['grid_controls']['%s|%s|%g' % (tag, ps, G)] = co
                print('grid', tag, ps, G, flush=True)
    # grid spread of the key quantities (data end differs slightly between runs: compare at the common crossing quantities
    spread = {}
    for ps in (0.5, 0.75):
        for G in (100.0, 1e3, 1e4, 1e5):
            vals = [res['primary']['%s|%g' % (ps, G)]] + [res['grid_controls']['%s|%s|%g' % (t, ps, G)] for t in GRID_CONTROLS]
            Nn = [v['analytic']['N_comoving_numeric'] for v in vals]
            qq = [v['q'] for v in vals]
            spread['%s|%g' % (ps, G)] = {'N_comoving_rel_spread': (max(Nn) - min(Nn))/min(Nn),
                                          'q_rel_spread': (max(qq) - min(qq))/min(qq)}
    res['grid_spread'] = spread
    # minus-branch secondary: crossing phi* = -0.5 before the reversal; production only
    tm0 = L.Trajectory(MINUS)
    ra = tm0.rec_all
    t_rev = float(ra[np.argmax(ra[:, 2] < 0.3), 0])      # stop before the reversal/collapse (h < 0.3 H0)
    tm = L.Trajectory(MINUS, t_cut=t_rev)
    res['minus_branch_t_cut'] = t_rev
    res['trajectory_checks'][MINUS] = {'s_range': [tm.s_min, tm.s_max], 'lna_crosscheck_max_abs': tm.lna_crosscheck}
    for G in (100.0, 1e3, 1e4):
        case = L.run_case(tm, G, -0.5)
        co = strip(case)
        for k in ('kappa', 'weights', 'n_k_order0_basis'): co.pop(k)
        res['minus_branch']['-0.5|%g' % G] = co
        print('minus', G, co['window'], co['analytic']['rel_diff_number'], flush=True)
    res['runtime_s'] = round(time.time() - t0, 1)
    (HERE/'PREHEATING_RESULTS.json').write_text(json.dumps(res, indent=1, allow_nan=False, default=float) + '\n')
    print('done', res['runtime_s'])


if __name__ == '__main__':
    main()
