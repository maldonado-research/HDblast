"""Mechanism screen, part 2: energy budget, residual vacuum, and physical scales.

For each detuning delta in {3e-4,1e-3,3e-3,1e-2,3e-2,1e-1}:
  * initial (unstable) static shell: cone near phi=-1, shell near phi=0 (frozen registered solver
    supplies the root; this script re-integrates the same radial ODE adding I=int (rho/rho_b)^2 dy),
  * final '+1' static branch (22 Sept package solver, imported read-only), same integral.
The graviton zero-mode normalisation of each static configuration gives the brane Planck mass
  1/kappa_4^2 = (2/kappa_5^2) int_0^{y_b} (rho/rho_b)^2 dy     (two Z2 copies),
compared with the local RS identification kappa_4^2 = kappa_5^2 sigma/6.
In a 4D effective description with a field-dependent Planck mass, the Einstein-frame vacuum energy
of a static configuration is V_E = 3 M_E^4 H^2 / M^2 (M_E fixed), so the fraction of the initial
vacuum energy that remains as residual vacuum is f_E = (H_vac^2/M_f^2)/(H_0^2/M_i^2) and radiation
created from the released energy can dominate the residual vacuum for at most
  N_max = (1/4) ln[(1-f_E)/f_E]  e-folds  (upper bound: instantaneous, lossless conversion).
Then physical scales: if H_vac is identified with today's dark-energy rate, the blast rate and delta
follow; and the c-tuning c* that removes the residual vacuum is computed from the series.
"""
import json, math, sys, time
import numpy as np
from scipy.integrate import solve_ivp, simpson
from common import *

sys.path.insert(0, str(PKG/'frozen')); sys.path.insert(0, str(PKG/'static_branch'))
import registered_solver as RS          # read-only import
import solve_plus_branch as PB          # read-only import

DELTAS = [3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
GAMMA_OVER_H0 = 1.657193663767          # registered instability rate (folder 152, finest eigenvalue)

def initial_shell(delta, rtol=2e-12):
    _, meta = RS.shell_background(delta, rtol=rtol)
    psi_h, yb = meta['psi_h'], meta['yb']; ph = psi_h - 1; y0 = 1e-4
    u, u1 = RS.stable_potential(psi_h)
    u2 = 4*ph*ph + 2*(ph*ph-1) - (4/3)*((ph*ph-1)**2 + 2*ph*W(ph)); k4 = u1*(u2/280 + u/630)
    ini = [y0 - u*y0**3/36, psi_h + u1*y0*y0/10 + k4*y0**4, u1*y0/5 + 4*k4*y0**3, y0**3/3]
    def rhs(y, v):
        rho, psi, py, _ = v; pot, pot1 = RS.stable_potential(psi)
        rpy = np.sqrt(1 + rho*rho*(py*py/12 - pot/6))
        return [rpy, py, pot1 - 4*rpy*py/rho, rho*rho]
    sol = solve_ivp(rhs, (y0, yb), ini, method='DOP853', rtol=rtol, atol=[1e-13, 1e-21, 1e-21, 1e-13])
    rho, psi, py, I = sol.y[:, -1]; phi_b = psi - 1
    H2 = 1/rho**2
    H2_id = sigma(phi_b, delta)**2/36 - sigma1(phi_b, delta)**2/48 + U(phi_b)/6
    return dict(delta=delta, rho_b=float(rho), phi_b=float(phi_b), y_b=float(yb), H2=float(H2),
                H2_identity_residual=float(H2 - H2_id), sigma=float(sigma(phi_b, delta)),
                I=float(I/rho**2), kappa5sq_over_kappa4sq=float(2*I/rho**2),
                RS_local_kappa5sq_over_kappa4sq=float(6/sigma(phi_b, delta)), root_residual=meta['root_residual'])

def plus_branch(delta, rtol=2e-12, y0=1e-4):
    row, prof = PB.solve(delta, rtol=rtol, y0=y0)
    ys, rho = prof[0], prof[1]
    I = simpson((rho/row['rho_b'])**2, x=ys) + y0**3/3/row['rho_b']**2
    phi_b = row['phi_b']
    return dict(delta=delta, rho_b=row['rho_b'], phi_b=phi_b, y_b=row['y_b'], H2=row['H2'],
                H2_identity_residual=row['H2_identity_residual'], sigma=float(sigma(phi_b, delta)),
                I=float(I), kappa5sq_over_kappa4sq=float(2*I),
                RS_local_kappa5sq_over_kappa4sq=float(6/sigma(phi_b, delta)), junction_residual=row['junction_residual'])

def budget(ini, fin):
    fH = fin['H2']/ini['H2']
    # careful: M^2 = kappa5sq_over_kappa4sq/kappa5^2 ; V_E ∝ H^2/M^2  =>  f_E = fH * M_i^2/M_f^2
    fE = fH*ini['kappa5sq_over_kappa4sq']/fin['kappa5sq_over_kappa4sq']
    fE_RS = fH*ini['RS_local_kappa5sq_over_kappa4sq']/fin['RS_local_kappa5sq_over_kappa4sq']
    Nmax = lambda f: 0.25*math.log((1-f)/f) if 0 < f < 0.5 else (0.0 if f >= 0.5 else None)
    return dict(H_vac2_over_H02=fH, M_i2_over_M_f2=ini['kappa5sq_over_kappa4sq']/fin['kappa5sq_over_kappa4sq'],
                f_E=fE, f_E_RS_local=fE_RS, f_brane_frame_fixedM=fH,
                radiation_to_vacuum_max_E=(1-fE)/fE, N_max_E=Nmax(fE), N_max_RS_local=Nmax(fE_RS), N_max_fixedM=Nmax(fH),
                H02_over_delta=ini['H2']/ini['delta'], Hvac2_over_delta=fin['H2']/fin['delta'])

def sudden_conversion_screening(ini, fin, eps_list=(0.01, 0.1, 0.5, 1.0)):
    """Local 5D algebra: convert a fraction eps of the tension drop into radiation at fixed
    sigma+kappa5^2 rho (nA continuous) at the final configuration.  H is continuous across the
    conversion, so the Weyl term must absorb the difference:
      Wy = H_before^2 - [(sigma_f+R)^2/36 - sigma_f'^2/48 + U_f/6].
    Returns the ratio of this Weyl term to the radiation's own contributions."""
    sf = fin['sigma']; d = fin['delta']; pf = fin['phi_b']; Hb2 = ini['H2']
    dsig = ini['sigma'] - sf; out = []
    for e in eps_list:
        R = e*dsig
        rad_terms = (sf*R/18 + R*R/36)
        Wy = Hb2 - ((sf+R)**2/36 - sigma1(pf, d)**2/48 + U(pf)/6)
        out.append(dict(eps=e, R=R, radiation_terms=rad_terms, Wy_required=Wy, Wy_over_radiation_terms=Wy/rad_terms,
                        R_over_Rcrit=R/(math.sqrt(sf*sf+36*fin['H2'])-sf)))
    return out

def physical(bud, ini, fin):
    """Map to physical units assuming (hypothesis) H_vac = observed dark-energy rate."""
    H_L = math.sqrt(OMEGA_L)*H0_TODAY_EV
    H0_blast = H_L/math.sqrt(bud['H_vac2_over_H02'])
    gamma = GAMMA_OVER_H0*H0_blast
    t_efold_Gyr = SEC_PER_EV_INV/gamma/GYR_S
    rows = []
    for ell_um in [0.1, 1.0, 10.0, 38.6]:
        Lm = ell_um*1e-6/9; invL = HBARC_EV_M/Lm
        delta_req = fin['delta']*(H_L/invL)**2/fin['H2']   # H_vac^2 = (H2/delta)*delta/L^2 (leading order linear in delta)
        rows.append(dict(ell_plus_micron=ell_um, inv_L_eV=invL, delta_required_for_Hvac_eq_HLambda=delta_req))
    return dict(H_Lambda_eV=H_L, H0_blast_eV=H0_blast, gamma_eV=gamma, efold_time_Gyr=t_efold_Gyr,
                note='delta enters H_vac^2 and H_0^2 linearly at leading order, so H0/H_vac is delta-independent; '
                     'identifying H_vac with the observed dark-energy rate fixes the blast rate independently of ell.',
                delta_required=rows)

def c_tuning(delta, H02, r):
    """Series (22 Sept package): H_vac^2 = delta(1+c)/27 + delta^2[(1+c)^2/36 - c^2/384] + O(delta^3)."""
    from scipy.optimize import brentq
    f = lambda c: delta*(1+c)/27 + delta**2*((1+c)**2/36 - c*c/384)
    cstar = brentq(f, -1.5, 0.0)
    need = {}
    for N in [1, 5, 10, 22]:
        # need f_E <= e^{-4N}/(1+e^{-4N}) ; f_E ≈ H_vac^2/H0^2 * (M_i^2/M_f^2) with ratio ~ 1/3 ; solve for |c-c*|
        target = math.exp(-4*N)/(1+math.exp(-4*N))
        slope = delta/27  # dH_vac^2/dc at leading order
        need[str(N)] = dict(max_abs_c_minus_cstar=target*H02/(r*slope))
    return dict(c_star=cstar, one_plus_c_star=1+cstar, leading_order_one_plus_cstar=27*delta/384,
                registered_c=C, required_precision_on_c_for_N_efolds=need,
                note='series-based (O(delta^3) remainder); changing c is a change of the registered model')

def reheating_temperature_window(bud, ini, fin, gstar=10.75):
    """If c were tuned so the residual vacuum is negligible, the maximum radiation density is the
    released Einstein-frame energy ~ V_E,i. Express T_max for a given (ell, delta) and the minimal delta
    giving T_max >= 5 MeV (BBN) for each ell."""
    out = []
    MP = MPL_RED_EV
    for ell_um in [0.1, 1.0, 10.0, 38.6]:
        Lm = ell_um*1e-6/9; invL = HBARC_EV_M/Lm
        # rho_max (in final-frame units) = 3 M_f^2 H0^2 (M_f^2/M_i^2)
        coef = 3*MP**2*(1/bud['M_i2_over_M_f2'])*bud['H02_over_delta']*invL**2   # times delta
        Tmin = 5e6; rho_need = math.pi**2/30*gstar*Tmin**4
        dmin = rho_need/coef
        out.append(dict(ell_plus_micron=ell_um, rho_max_per_delta_eV4=coef, delta_min_for_Tmax_5MeV=dmin,
                        H0_at_delta_min_eV=math.sqrt(bud['H02_over_delta']*dmin)*invL))
    return out

def main():
    t0 = time.time(); rows = []
    for d in DELTAS:
        ini = initial_shell(d); ini2 = initial_shell(d, rtol=1e-13)
        fin = plus_branch(d); fin2 = plus_branch(d, rtol=8e-14, y0=5e-5)
        b = budget(ini, fin); b2 = budget(ini2, fin2)
        rows.append(dict(delta=d, initial=ini, final=fin, budget=b,
                         tolerance_check=dict(f_E_rel_change=abs(b2['f_E']-b['f_E'])/b['f_E'],
                                              I_i_rel_change=abs(ini2['I']-ini['I'])/ini['I'],
                                              I_f_rel_change=abs(fin2['I']-fin['I'])/fin['I'])))
        print('delta=%g rho_i=%.6f H0^2/d=%.6f Hv^2/d=%.6f fH=%.6f Mi2/Mf2=%.6f (RS %.6f) fE=%.6f Nmax=%.4f  tol=%.1e' % (
            d, ini['rho_b'], b['H02_over_delta'], b['Hvac2_over_delta'], b['H_vac2_over_H02'], b['M_i2_over_M_f2'],
            ini['RS_local_kappa5sq_over_kappa4sq']/fin['RS_local_kappa5sq_over_kappa4sq'], b['f_E'], b['N_max_E'],
            rows[-1]['tolerance_check']['f_E_rel_change']), flush=True)
    reg = [r for r in rows if r['delta'] == 1e-3][0]
    # controls
    controls = dict(
        registered_rho_b_matches_archive=abs(reg['initial']['rho_b'] - 78.82817714224423) < 1e-6,
        registered_plus_H2_matches_package=abs(reg['final']['H2'] - 5.9240147943292764e-05) < 1e-15,
        H2_identity_initial_all=max(abs(r['initial']['H2_identity_residual']) for r in rows) < 1e-10,
    )
    # known-limit control: pure AdS l=9 dS brane at the plus-branch rho_b: I = [sinh(2k y)/(4k) - y/2]/sinh^2(k y_b) with k=1/9
    k = 1/9; yb = math.asinh(k*reg['final']['rho_b'])/k
    I_ads = (math.sinh(2*k*yb)/(4*k) - yb/2)/math.sinh(k*yb)**2
    controls['analytic_AdS_I'] = I_ads; controls['plus_branch_I_minus_analytic_rel'] = (reg['final']['I']-I_ads)/I_ads
    controls['plus_branch_matches_analytic_dS_in_AdS_I_to_1e-6'] = abs(controls['plus_branch_I_minus_analytic_rel']) < 1e-6
    I_wrong = (math.cosh(k*yb)-1)/(k*math.sinh(k*yb))      # wrong-formula control: unsquared warp factor
    controls['wrong_formula_unsquared_warp_I'] = I_wrong
    controls['wrong_formula_control_differs'] = abs(I_wrong-reg['final']['I'])/reg['final']['I'] > 0.5
    out = dict(status='numerical (floating-point static solutions) + conditional (4D effective energy budget)',
               rows=rows, controls=controls,
               sudden_conversion_screening_registered=sudden_conversion_screening(reg['initial'], reg['final']),
               physical_scales_registered=physical(reg['budget'], reg['initial'], reg['final']),
               c_tuning_registered_delta=c_tuning(1e-3, reg['initial']['H2'], reg['budget']['M_i2_over_M_f2']),
               reheating_temperature_window_if_c_tuned=reheating_temperature_window(reg['budget'], reg['initial'], reg['final']),
               observational_target=dict(
                   rho_r_over_rho_Lambda_at_T5MeV=(math.pi**2/30*10.75*(5e6)**4)/(3*MPL_RED_EV**2*OMEGA_L*H0_TODAY_EV**2),
                   N_equivalent=0.25*math.log((math.pi**2/30*10.75*(5e6)**4)/(3*MPL_RED_EV**2*OMEGA_L*H0_TODAY_EV**2))),
               runtime_s=time.time()-t0)
    dump('M2_ENERGY_BUDGET_AND_SCALES.json', out)
    print(json.dumps(dict(controls=controls, sudden=out['sudden_conversion_screening_registered'],
                          phys=out['physical_scales_registered'], ctune=out['c_tuning_registered_delta'],
                          window=out['reheating_temperature_window_if_c_tuned'], target=out['observational_target']), indent=1))

if __name__ == '__main__':
    main()
