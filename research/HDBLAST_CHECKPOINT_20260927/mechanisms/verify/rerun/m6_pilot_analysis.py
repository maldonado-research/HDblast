"""Mechanism screen, part 6: analysis of the PILOT coupled 5D runs (pilot_5d/runs, delta = 0.1).

Compares, for each friction coupling Y, the shell radiation left at the end of the recorded roll with
(i) the 4D effective energy-budget bound  rho_r/rho_vac <= (1-f_E)/f_E  from M2 at delta = 0.1, and
(ii) the radiation-domination threshold R_crit = sqrt(sigma_f^2 + 36 H_vac^2) - sigma_f.
Control: Y = 0 must reproduce the archived Chat14 run v2_t01_plus (same numerics, same parameters).
"""
import io, json, math, zipfile, glob
import numpy as np
from common import *

def main():
    m2 = json.load(open(HERE/'M2_ENERGY_BUDGET_AND_SCALES.json'))
    row = [r for r in m2['rows'] if abs(r['delta']-0.1) < 1e-12][0]
    sf, Hv2, fE = row['final']['sigma'], row['final']['H2'], row['budget']['f_E']
    Rcrit = math.sqrt(sf*sf + 36*Hv2) - sf
    bound = (1-fE)/fE
    with zipfile.ZipFile(CHAT14_ZIP) as z:
        ref = json.loads(z.read(CHAT14_PREFIX+'runs/v2_t01_plus_summary.json'))
    out = dict(delta=0.1, sigma_f=sf, H_vac2=Hv2, R_crit=Rcrit, budget_bound_rho_r_over_rho_vac=bound,
               N_max_4D_budget=row['budget']['N_max_E'], runs={})
    for f in sorted(glob.glob(str(HERE/'pilot_5d/runs/*_summary.json'))):
        s = json.load(open(f)); tag = Path(f).name.replace('_summary.json', '')
        d = np.load(f.replace('_summary.json', '_timeseries.npz'))['rec']
        rb = s['rho_b']; R = d[:, 6]; h = d[:, 2]; ss = d[:, 8]
        lna = np.concatenate(([0.], np.cumsum(0.5*(h[1:]+h[:-1])*np.diff(ss))))
        # reliable window: stop where the conformal chart freezes (d(H0 tau)/dt < 0.2, lapse collapse as in Chat14)
        lapse = np.gradient(ss, d[:, 0])
        bad = np.where((lapse < 0.2) & (d[:, 0] > 1.0))[0]
        i_end = int(bad[0]) - 1 if len(bad) else len(d)-1
        # dilution-invariant comparison: R a^4 extrapolated to the moment the roll ends is what is measured;
        ratio_end = R[i_end]/Rcrit
        v = d[:, 5]/rb; tau = ss*rb; Y = s['Y']
        work = float(np.sum((0.5*(Y*v[1:]**2 + Y*v[:-1]**2)*np.diff(tau))[:i_end]))
        from common import sigma as sig_
        cc = s.get('c', C); dq = s.get('d_quad', 0.0)
        sgf = lambda p: 2*W(p) + 0.1*(1 + cc*p + dq*p*p/2)
        dsig = float(sgf(d[0, 1]) - sgf(d[i_end, 1]))
        dphi = np.abs(d[:, 1] - s['phi_b_static']); t = d[:, 0]
        mfit = (t > 3) & (dphi < 0.05) & (dphi > 0)
        gfit = float(np.polyfit(t[mfit], np.log(dphi[mfit]), 1)[0]) if mfit.sum() > 20 else None
        out['runs'][tag] = dict(summary=s, R_end=float(R[i_end]), R_max=float(R.max()), H0tau_at_Rmax=float(ss[np.argmax(R)]),
                                R_end_over_Rcrit=float(ratio_end),
                                N_rad_after_end=float(0.25*math.log(ratio_end)) if ratio_end > 1 else 0.0,
                                comoving_radiation_R_a4_end=float(R[i_end]*math.exp(4*lna[i_end])),
                                within_4D_budget=bool(ratio_end <= bound),
                                reliable_window_end_t=float(d[i_end, 0]), reliable_window_end_H0tau=float(ss[i_end]),
                                phi_b_at_window_end=float(d[i_end, 1]), H_over_H0_at_window_end=float(h[i_end]),
                                R_max_over_Rcrit=float(R.max()/Rcrit),
                                seed_artefact_flag=bool(s['dc'] > 1e-3),
                                Wy_over_H0sq_end=float(d[i_end, 7]),
                                matter_work_int_Yv2_dtau=work, tension_drop=dsig,
                                local_matter_fraction_of_tension_drop=(work/dsig if dsig > 1e-6 else None),
                                growth_rate_per_H0_fit_small_amplitude=gfit, phi_b_end=s['phi_b_end'], H_over_H0_end=s['H_over_H0_end'])
    y0 = out['runs'].get('pilot_Y0')
    if y0:
        out['control_Y0_vs_archived_v2_t01_plus'] = dict(
            archived_phi_b_end=ref['phi_b_end'], pilot_phi_b_end=y0['phi_b_end'],
            archived_H_end=ref['HJ_over_H0_end'], pilot_H_end=y0['H_over_H0_end'],
            abs_diff_phi=abs(ref['phi_b_end']-y0['phi_b_end']), abs_diff_H=abs(ref['HJ_over_H0_end']-y0['H_over_H0_end']),
            passes=abs(ref['phi_b_end']-y0['phi_b_end']) < 1e-8 and abs(ref['HJ_over_H0_end']-y0['H_over_H0_end']) < 1e-8)
    out['notes'] = ('All quantities are evaluated at the end of the reliable window (before the conformal chart freezes, '
                    'd(H0 tau)/dt < 0.2). Values recorded after the freeze (e.g. the dc=1e-2 Y=0 summary H/H0=0.294) are not physical. '
                    'Seed-size control: Y=0.3 with dc=1e-2 vs 1e-4 agrees at window end (R/Rcrit, Wy, capture fraction) to ~2%.')
    a, b = out['runs'].get('pilot_Y0.3_L21'), out['runs'].get('pilot_Y0.3_dc1e-2')
    if a and b:
        out['seed_size_control_Y0.3'] = dict(R_over_Rcrit=[a['R_end_over_Rcrit'], b['R_end_over_Rcrit']],
            capture_fraction=[a['local_matter_fraction_of_tension_drop'], b['local_matter_fraction_of_tension_drop']],
            Wy=[a['Wy_over_H0sq_end'], b['Wy_over_H0sq_end']],
            max_rel_diff=max(abs(a[k]-b[k])/abs(a[k]) for k in ['R_end_over_Rcrit', 'local_matter_fraction_of_tension_drop', 'Wy_over_H0sq_end']))
    # c-tuned pilot: residual vacuum ~0 by construction (series c*); report the radiation share of H^2 at window end
    ct = {}
    for tag, r in out['runs'].items():
        if 'cstar' not in tag and 'dstar' not in tag: continue
        f = HERE/'pilot_5d'/'runs'/(tag + '_timeseries.npz'); d = np.load(f)['rec']; s_ = r['summary']
        lapse = np.gradient(d[:, 8], d[:, 0]); bad = np.where((lapse < 0.2) & (d[:, 0] > 1.0))[0]
        i = int(bad[0]) - 1 if len(bad) else len(d)-1
        rb = s_['rho_b']; phi = d[i, 1]; Rv = d[i, 6]; h = d[i, 2]
        dq = s_.get('d_quad', 0.0)
        sgm = 2*W(phi) + 0.1*(1 + s_['c']*phi + dq*phi*phi/2)
        rad_terms = (sgm*Rv/18 + Rv*Rv/36)*rb*rb           # in H0^2 units
        ct[tag] = dict(c=s_['c'], d_quad=dq, H0tau_of_R_share_samples=None, window_end_t=float(d[i, 0]), window_end_H0tau=float(d[i, 8]), phi_b=float(phi),
                       H_over_H0=float(h), R=float(Rv), radiation_H2_share=float(rad_terms/h**2),
                       Wy_over_H0sq=float(d[i, 7]), Wy_over_radiation_terms=float(d[i, 7]/rad_terms) if rad_terms > 0 else None,
                       H_over_H0_min_in_window=float(d[:i+1, 2].min()))
    for tag in ct:
        d = np.load(HERE/'pilot_5d'/'runs'/(tag + '_timeseries.npz'))['rec']; rb = out['runs'][tag]['summary']['rho_b']; dq = ct[tag]['d_quad']
        sg_ = 2*W(d[:, 1]) + 0.1*(1 + ct[tag]['c']*d[:, 1] + dq*d[:, 1]**2/2)
        share = (sg_*d[:, 6]/18 + d[:, 6]**2/36)*rb*rb/np.maximum(d[:, 2]**2, 1e-30)
        idx = np.linspace(0, len(d)-1, 25).astype(int)
        ct[tag]['H0tau_of_R_share_samples'] = [[float(d[i, 8]), float(d[i, 2]), float(share[i]), float(d[i, 7])] for i in idx]
        ss = d[:, 8]; hh = d[:, 2]
        lna = np.concatenate(([0.], np.cumsum(0.5*(hh[1:]+hh[:-1])*np.diff(ss))))
        lapse = np.gradient(ss, d[:, 0]); bad = np.where((lapse < 0.2) & (d[:, 0] > 1.0))[0]
        ie = int(bad[0]) - 1 if len(bad) else len(d)-1
        i0 = int(np.searchsorted(ss, ss[ie] - 0.2))
        wa4 = d[:, 7]*np.exp(4*lna); ra4 = d[:, 6]*np.exp(4*lna)
        ct[tag]['last_0.2_H0tau'] = dict(Wy_a4_rel_change=float((wa4[ie]-wa4[i0])/wa4[i0]),
                                          R_a4_rel_change=float((ra4[ie]-ra4[i0])/ra4[i0]) if ra4[i0] > 0 else None,
                                          ln_a_advance=float(lna[ie]-lna[i0]))
    out['tuned_vacuum_pilot'] = ct
    out['c_tuned_pilot'] = dict(status='FAILED: the initial unstable shell could not be constructed at c = c*(0.1) = -0.99307; '
        'see pilot_5d/runs/pilot_cstar_*.log and M7_INITIAL_SHELL_VS_C.json (the registered initial-shell family ends near c -> 0+, '
        'where phi_b -> -1). The vacuum is instead tuned with the quadratic tension term d (see tuned_vacuum_pilot).')
    dump('M6_PILOT_COUPLED_5D.json', out)
    for k, v in ct.items(): print('TUNED', k, {kk: vv for kk, vv in v.items() if kk != 'H0tau_of_R_share_samples'})
    for k, r in out['runs'].items():
        print(k, 'Y=%g dc=%g win_end_t=%.2f H0tau=%.3f phi=%.4f h=%.4f R=%.3e R/Rc=%.2e Rmax/Rc=%.2e frac=%s g=%s Wy=%.4f' % (
            r['summary']['Y'], r['summary']['dc'], r['reliable_window_end_t'], r['reliable_window_end_H0tau'], r['phi_b_at_window_end'],
            r['H_over_H0_at_window_end'], r['R_end'], r['R_end_over_Rcrit'], r['R_max_over_Rcrit'],
            r['local_matter_fraction_of_tension_drop'], r['growth_rate_per_H0_fit_small_amplitude'], r['Wy_over_H0sq_end']))
    print('Rcrit', Rcrit, 'bound', bound, out.get('control_Y0_vs_archived_v2_t01_plus', {}).get('passes'))
    return


if __name__ == '__main__':
    main()
