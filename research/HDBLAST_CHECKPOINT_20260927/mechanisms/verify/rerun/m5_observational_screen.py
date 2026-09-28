"""Mechanism screen, part 5: observational mapping and order-of-magnitude screens.

Inputs: M1 (Weyl term from Chat14), M2 (energy budget, scales). Outputs M5_OBSERVATIONAL_SCREEN.json.
  (a) Delta N_eff -> allowed ratio of Weyl ('dark') radiation to Standard-Model radiation, at BBN and at
      a production epoch with g*=106.75 (entropy-conserving evolution in between).
  (a') applied to the best-resolved Chat14 '+1' trajectory (Weyl term at chart freeze).
  (c) collapse branch: fraction of the resolved contraction with w_eff > 1 (ekpyrotic smoothing needs
      w >> 1 sustained); growth of an anisotropic shear fraction if w_eff stays at its median value;
      the Hubble rate a non-singular bounce would need to reach for a BBN-compatible radiation era.
  (e) gravitational particle production: rho ~ C_g g H^4 vs the released vacuum energy (C_g <= 1e-2 is an
      order-of-magnitude assumption, not a computed Bogoliubov spectrum).
"""
import io, json, math, zipfile
import numpy as np
from common import *

def dneff_to_ratio(dn, g_bbn=10.75, g_prod=106.75):
    r_bbn = (7/4)*dn/g_bbn                       # rho_DR/rho_SM at T ~ 1 MeV (neutrinos still coupled)
    r_prod = r_bbn*(g_prod/g_bbn)**(1/3)         # rho_SM a^4 ∝ g^{-1/3}
    return dict(delta_N_eff=dn, rho_DR_over_rho_SM_at_BBN=r_bbn, rho_DR_over_rho_SM_at_production=r_prod)

def collapse_stats():
    out = {}
    with zipfile.ZipFile(CHAT14_ZIP) as z:
        for r in ['v2_t1e3_minus_c4', 'v2_t1e3_minus_c4_fine', 'v2_t1e3_minus_c4_finecoarse']:
            rec = np.load(io.BytesIO(z.read(CHAT14_PREFIX+'runs/'+r+'_timeseries.npz')))['rec']
            s, h = rec[:, 9], rec[:, 2]
            keep = np.concatenate(([True], np.diff(s) > 1e-12)); s, h = s[keep], h[keep]
            eps = -np.gradient(h, s, edge_order=2)/h**2; w = 2*eps/3 - 1
            sel = (h < -0.5) & (h > -10)        # range where the two wall resolutions agree to <= 8%
            dlna = np.abs(h[sel][:-1]*np.diff(s[sel]))
            wmid = 0.5*(w[sel][:-1] + w[sel][1:])
            out[r] = dict(ln_a_contraction_resolved=float(dlna.sum()),
                          fraction_of_ln_a_with_w_gt_1=float(dlna[wmid > 1].sum()/dlna.sum()),
                          ln_a_weighted_mean_w=float(np.sum(wmid*dlna)/dlna.sum()),
                          H_range=[-0.5, -10])
    return out

def main():
    m1 = json.load(open(HERE/'M1_WEYL_ALONG_CHAT14.json')); m2 = json.load(open(HERE/'M2_ENERGY_BUDGET_AND_SCALES.json'))
    reg = [r for r in m2['rows'] if r['delta'] == 1e-3][0]; bud = reg['budget']
    bounds = {'conservative_BBN_level_0.3': dneff_to_ratio(0.3), 'combined_CMB_BBN_BAO_0.09': dneff_to_ratio(0.09)}
    # (a') Chat14 best-resolved plus run
    a = m1['runs']['v2_t1e3_plus_c4_finecoarse']['analysis']
    Wy_end = a['Wy_over_H0sq_end']; h_end = a['H_over_H0_end']
    fH = bud['H_vac2_over_H02']
    weyl = dict(run='v2_t1e3_plus_c4_finecoarse', Wy_over_H0sq_at_freeze=Wy_end, H_over_H0_at_freeze=h_end,
                Wy_over_Hvac2=Wy_end/fH, peak_Wy_over_H0sq=a['Wy_over_H0sq_max'], peak_at_H0tau=a['s_at_Wy_max'],
                same_quantity_other_runs={k: v['analysis']['Wy_over_H0sq_end'] for k, v in m1['runs'].items() if 'plus' in k},
                required_radiation_H2_contribution_over_H0sq={k: Wy_end/b['rho_DR_over_rho_SM_at_production'] for k, b in bounds.items()},
                max_radiation_H2_contribution_over_H0sq_from_4D_budget=bud['radiation_to_vacuum_max_E']*fH,
                caveat='value at chart freeze; default-far-grid runs give 0.0004-0.0018 and the ext run turns negative, '
                       'so the late Weyl amplitude is NOT converged (inconclusive); the peak ~0.081 H0^2 is grid-robust')
    # (c) collapse / bounce
    cs = collapse_stats()
    wmed = np.median([v['ln_a_weighted_mean_w'] for v in cs.values()])
    target = m2['observational_target']['rho_r_over_rho_Lambda_at_T5MeV']
    Hb_over_Hvac = math.sqrt(target)                     # radiation-dominated right after a bounce: rho_r = 3 M^2 H_b^2
    Hb_over_H0 = Hb_over_Hvac*math.sqrt(fH)
    # shear/background ∝ a^{-6}/a^{-3(1+w)} = a^{-3(1-w)} ; |H| ∝ a^{-3(1+w)/2}  => shear fraction ∝ |H|^{2(1-w)/(1+w)}
    p = 2*(1-wmed)/(1+wmed)
    collapse = dict(runs=cs, lna_weighted_mean_w_median=float(wmed),
                    bounce_H_required_over_H_vac=Hb_over_Hvac, bounce_H_required_over_H0=Hb_over_H0,
                    bounce_H_required_eV_if_Hvac_is_HLambda=Hb_over_Hvac*m2['physical_scales_registered']['H_Lambda_eV'],
                    shear_fraction_growth_exponent_in_H=p,
                    shear_fraction_growth_from_H0_to_bounce=float(Hb_over_H0**p),
                    note='ekpyrotic smoothing requires w >> 1 sustained; resolved collapse has mean w well below 1')
    # (e) gravitational particle production
    phys = m2['physical_scales_registered']; Cg, gdof = 1e-2, 100
    gpp = {}
    for label, H0 in [('registered_c_Hvac_eq_HLambda', phys['H0_blast_eV']),
                      ('c_tuned_Tmax_5MeV_minimal', m2['reheating_temperature_window_if_c_tuned'][0]['H0_at_delta_min_eV'])]:
        released = (1-bud['f_E'])/bud['f_E']*3*MPL_RED_EV**2*fH*H0**2     # (1-f)/f * V_vac (upper bound, brane units)
        gpp[label] = dict(H0_eV=H0, rho_gpp_eV4=Cg*gdof*H0**4, ratio_to_released_energy=Cg*gdof*H0**4/released,
                          H_needed_for_order_one_over_Mpl=math.sqrt(3/(Cg*gdof)))
    out = dict(status='numerical/conditional order-of-magnitude screen', delta_N_eff_bounds=bounds,
               weyl_from_blast=weyl, collapse_branch=collapse, gravitational_particle_production=gpp,
               assumptions=dict(C_g=Cg, g_dof=gdof, gstar_BBN=10.75, gstar_production=106.75))
    dump('M5_OBSERVATIONAL_SCREEN.json', out)
    print(json.dumps(out, indent=1))

if __name__ == '__main__':
    main()
