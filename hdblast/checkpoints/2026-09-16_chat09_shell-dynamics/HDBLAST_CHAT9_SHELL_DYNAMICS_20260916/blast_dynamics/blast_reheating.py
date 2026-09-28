#!/usr/bin/env python3
"""Reheating ledger for MODIFIED detunings with V = 0 somewhere (NOT the registered tension).
 (A) V = t (1-phi)(1+b phi), b = 1+c_star  (zero at the field-space end phi_b=+1; hilltop kept at phi_b=0)
 (B) V = t (1-phi/phi_m)^2 (1+b phi), b = c_star + 2/phi_m  (interior Minkowski minimum at phi_m; hilltop kept at phi_b=0)
EFT quantities at the minimum: m^2 = d^2V_E/dTheta^2, alpha = -(1/2) dln f/dTheta (coupling of brane-frame matter),
Gamma ~ alpha^2 m^3/M_pl^2 (O(1) prefactor not computed), T_rh ~ 0.5 sqrt(Gamma M_pl)."""
import numpy as np, json
from blast_common import *
from blast_homogeneous import evolve, first
Mpl_GeV = 2.435e18; T_BBN = 4e-3/Mpl_GeV
T, C, Th_end = build_tables(); R = {}
c = c_star
# ---------- (A)
b = 1 + c
LA = Landscape(lambda p: (1-p)*(1+b*p), lambda p: -(1+b*p) + b*(1-p), T, C, Th_end)
d2 = np.gradient(LA.dv, LA.Th)
sel = (LA.phi > 0.99) & (LA.phi < 1-1e-9) & (LA.v > 0)
dTh = LA.Th[sel] - LA.Th_plus_end; slope = np.polyfit(np.log(dTh[::200]), np.log(LA.v[sel][::200]), 1)[0]
R['A'] = dict(potential="t(1-phi)(1+b phi), b=1+c", mu2_hilltop=float(3*d2[LA.i0]/LA.v0), dlnV_E_dln_DeltaTheta_near_end=float(slope),
              statement="V_E ~ DeltaTheta^18 near phi_b=+1 because 1-phi ~ e^{2y} and DeltaTheta ~ e^{y/9}: the EFT modulus mass at the end point is ZERO; no oscillation there.")
t = 1e-3; HE = np.sqrt(t*LA.v0/3); dq = HE/(2*np.pi)
sA = (-3+np.sqrt(9-4*R['A']['mu2_hilltop']))/2
o = evolve(LA, -dq, -sA*dq, 0, 60.0, dtau=5e-4)
ik = int(np.argmax(o['K'])); w = (o['K']-o['u'])/(o['K']+o['u']); N = np.log(o['a']/o['a'][0])
R['A'].update(s_linear=float(sA), N_E_at_10pct_drop=float(N[first(o['u'] < 0.9*o['u'][0])]), peak_K_over_V0=float(o['K'][ik]), N_E_at_peak_K=float(N[ik]),
              eps_H_gt_1_at_tau=(float(o['tau'][first(3*o['K'] > o['h']**2)]) if first(3*o['K'] > o['h']**2) is not None else None),
              end=dict(tau=float(o['tau'][-1]), DeltaTheta=float(o['Th'][-1]-LA.Th_plus_end), w=float(w[-1]), h=float(o['h'][-1]),
                       sinhpsi_over_sqrt_t=float(o['sinhpsi_over_sqrt_t'][-1]), reached_boundary=bool(o['Th'][-1] <= LA.Th_plus_end+1e-12)))
print(json.dumps(R['A'], indent=1))
# ---------- (B)
R['B'] = {}
for pm in (0.5, 0.9, -0.5):
    b = c + 2/pm
    Vf = lambda p, pm=pm, b=b: (1-p/pm)**2*(1+b*p)
    dVf = lambda p, pm=pm, b=b: -2/pm*(1-p/pm)*(1+b*p) + b*(1-p/pm)**2
    LB = Landscape(Vf, dVf, T, C, Th_end); d2 = np.gradient(LB.dv, LB.Th)
    im = np.argmin(np.abs(LB.phi - pm) + 10*(LB.branch > 0)); 
    m2_over_H2 = 3*d2[im]/LB.v0; alpha = -0.5*LB.fTh[im]/LB.f[im]; mu2 = 3*d2[LB.i0]/LB.v0
    # closed form check: m^2 = V''(phi_m)/(f^2 Z_E(phi_m)) , Z_E = Z/f + 1.5 (f'/f)^2 with d/dphi
    ZEphi = 1/LB.phiTh[im]**2; m2_cf = (2*(1+b*pm)/pm**2)/(LB.f[im]**2*ZEphi)/(LB.v0/3)
    rm = np.sqrt(m2_over_H2)
    r = dict(phi_m=pm, b=b, mu2_hilltop=float(mu2), Theta_m=float(LB.Th[im]), f_m=float(LB.f[im]), m_over_H_top=float(rm), m_over_H_top_closed_form=float(np.sqrt(m2_cf)),
             alpha=float(alpha), Gamma_over_Htop_per_HoverMpl_squared=float(alpha**2*rm**3), Trh_over_Mpl_per_HoverMpl_to_1p5=float(0.5*abs(alpha)*rm**1.5))
    # scales
    tab = {}
    for tt in (1e-3, 1e-6, 1e-10, 1e-20, 1e-28):
        H = np.sqrt(tt*LB.v0/3); m = rm*H; Gam = alpha**2*m**3; Trh = 0.5*np.sqrt(Gam)
        tab["t=%g" % tt] = dict(H_top_over_Mpl=float(H), H_top_GeV=float(H*Mpl_GeV), m_GeV=float(m*Mpl_GeV), Gamma_over_H_top=float(Gam/H), T_rh_GeV=float(Trh*Mpl_GeV), passes_BBN_4MeV=bool(Trh > T_BBN),
                                efolds_matter_like_oscillation=float((2/3)*np.log(H/Gam)))
    r['scales'] = tab
    Hmin = (T_BBN/(0.5*abs(alpha)*rm**1.5))**(2/3)
    r['BBN_bound'] = dict(H_top_over_Mpl_min=float(Hmin), H_top_GeV_min=float(Hmin*Mpl_GeV), m_GeV_min=float(rm*Hmin*Mpl_GeV), t_min=float(3*Hmin**2/LB.v0))
    # dynamics, k=0: e-folds, oscillation
    H = np.sqrt(1e-3*LB.v0/3); dq = H/(2*np.pi); s_ = (-3+np.sqrt(9-4*mu2))/2; sg = -1 if pm > 0 else 1
    o = evolve(LB, sg*dq, sg*s_*dq, 0, 12.0, dtau=1e-4); N = np.log(o['a']/o['a'][0])
    icr = first((o['phi']-pm)*sg*(-1) >= 0) if pm > 0 else first(o['phi'] <= pm)
    ik = int(np.argmax(o['K']))
    r['k0_dynamics'] = dict(s_linear=float(s_), N_E_first_crossing_of_minimum=float(N[icr]), tau_first_crossing=float(o['tau'][icr]), K_over_V0_at_crossing=float(o['K'][icr]),
                            fraction_of_V0_lost_to_Hubble_friction_by_then=float(o['fric'][icr]), eps_H_max=float((3*o['K']/o['h']**2).max()),
                            phi_b_range_after=[float(o['phi'][icr:].min()), float(o['phi'][icr:].max())], max_sinhpsi_t1em3=float(np.abs(o['sinhpsi_over_sqrt_t']).max()*np.sqrt(1e-3)),
                            rho_a3_late_over_V0=float(((o['K']+o['u'])*o['a']**3)[-1]), a_end=float(o['a'][-1]))
    # closed shell: numerical recollapse demo with a large displacement, and analytic a_max for the quantum displacement
    oc = evolve(LB, sg*0.08, 0.0, 1, 400.0, dtau=1e-3, stop=lambda s: s[3] < 0 and s[2] < 0.5)
    itn = first(oc['h'] < 0)
    ioc = first((oc['phi']-pm)*(1 if pm < 0 else -1) >= 0)
    rho_a3 = ((oc['K']+oc['u'])*oc['a']**3)
    r['closed_demo(dTheta0=0.08)'] = dict(a_bounce=float(oc['a'][0]), tau_max_expansion=(float(oc['tau'][itn]) if itn else None), a_max=(float(oc['a'][itn]) if itn else None),
                                         a_max_dust_estimate=float(np.mean(rho_a3[itn-2000:itn])) if itn else None, recollapses=bool(itn is not None))
    oq = evolve(LB, sg*dq, 0.0, 1, 12.0, dtau=1e-4); rq = ((oq['K']+oq['u'])*oq['a']**3)[-20000:].mean()
    r['closed_quantum_start(t=1e-3)'] = dict(a_bounce=float(oq['a'][0]), rho_a3_oscillating=float(rq), a_max_if_dust=float(rq), a_max_over_a_bounce=float(rq/oq['a'][0]),
                                            time_to_max_expansion_dust_in_inv_Htop=float(np.pi*rq/2), total_lifetime_dust_in_inv_Htop=float(np.pi*rq))
    R['B']["phi_m=%g" % pm] = r; print(json.dumps(r, indent=1))
json.dump(R, open("BLAST_REHEATING_RESULTS.json", "w"), indent=1)
