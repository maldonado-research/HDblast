#!/usr/bin/env python3
"""Homogeneous roll-off of the shell modulus in the Einstein-frame EFT (registered linear detuning).
Units: tau = H_top * t_E, H_top^2 = t v0/3 (M_pl = 1).  u = v/v0.  Equations (t-independent in these units):
   Th'' + 3 h Th' + 3 u_Th = 0,   h^2 = Th'^2/6 + u - k/a^2,   a''/a = u - Th'^2/3 .
Only the initial quantum displacement dTheta = H_top/(2 pi) and the 5D-unit diagnostics (rapidity) depend on t."""
import numpy as np, json, sys
from blast_common import *

def evolve(L, th0, dth0, k, tau_max, dtau=5e-4, a0=None, stop=None):
    v0 = L.v0; Th = L.Th; uarr = L.v/v0; duarr = L.dv/v0
    thp = L.Th_plus_end; vinf = L.v[0]/v0
    def u_du(th):
        if th < L.Th_min:                       # analytic tail toward phi_b=+1: v = vinf_true (1 + dTh^2/3)
            d = th - thp; base = vinf/(1 + L.dTh_min**2/3)
            return base*(1 + d*d/3), base*2*d/3
        return np.interp(th, Th, uarr), np.interp(th, Th, duarr)
    u0, _ = u_du(th0)
    if k == 0: a = 1.0 if a0 is None else a0; ad = a*np.sqrt(dth0**2/6 + u0)
    else:      a = 1/np.sqrt(dth0**2/6 + u0); ad = 0.0          # time-symmetric bounce, adot = 0
    def rhs(s):
        th, dth, a, ad = s; u, du = u_du(th)
        return np.array([dth, -3*(ad/a)*dth - 3*du, ad, a*(u - dth*dth/3)])
    s = np.array([th0, dth0, a, ad]); n = int(tau_max/dtau); rec = np.empty((n+1, 6)); fric = 0.0
    for i in range(n+1):
        u, du = u_du(s[0]); rec[i] = [i*dtau, s[0], s[1], s[2], s[3], fric]
        if s[0] >= Th[-1] or s[0] <= thp or not np.isfinite(s).all() or (stop and stop(s)): rec = rec[:i+1]; break
        if i == n: break
        k1 = rhs(s); k2 = rhs(s+dtau/2*k1); k3 = rhs(s+dtau/2*k2); k4 = rhs(s+dtau*k3)
        sn = s + dtau/6*(k1+2*k2+2*k3+k4)
        fric += dtau*0.5*(3*(s[3]/s[2])*s[1]**2 + 3*(sn[3]/sn[2])*sn[1]**2)   # int 3 h Th'^2 dtau; d(K+u)/dtau = -h Th'^2 -> /3 below
        s = sn
    tau, th, dth, a, ad, fr = rec.T
    inside = th >= L.Th_min
    u = np.where(inside, np.interp(th, Th, uarr), (vinf/(1+L.dTh_min**2/3))*(1+(th-thp)**2/3))
    h = ad/a
    out = dict(tau=tau, Th=th, dTh=dth, a=a, h=h, u=u, K=dth**2/6, fric=fr/3, k=k,
               cons=(h*h - (dth**2/6 + u - k/a**2)))
    f = np.where(inside, np.interp(th, Th, L.f), 9 - 1.5*(th-thp)**2)
    fTh = np.where(inside, np.interp(th, Th, L.fTh), -3*(th-thp))
    dThdy = np.where(inside, np.interp(th, Th, L.dThdy), (th-thp)/9)
    phi = np.where(inside, np.interp(th, Th, L.phi), 1.0)
    out.update(f=f, phi=phi, y=np.where(inside, np.interp(th, Th, L.y), np.nan), branch=np.interp(th, Th, L.branch),
               aJ=a/np.sqrt(f)*np.sqrt(L.f[L.i0]), hJ=np.sqrt(f/L.f[L.i0])*(h - fTh*dth/(2*f)),
               # bulk velocity of the shell relative to the phi-leaves, per sqrt(t):  sinh(psi)/sqrt(t)
               sinhpsi_over_sqrt_t=np.sqrt(f*v0/3)*dth/dThdy,
               # d phi_b / d tau_J in 5D units per sqrt(t)
               phidotJ_over_sqrt_t=np.sqrt(f*v0/3)*dth*np.where(inside, np.interp(th, Th, L.phiTh), 0.0))
    return out

def first(cond):
    i = np.where(cond)[0]; return int(i[0]) if len(i) else None

if __name__ == "__main__":
    T, C, Th_end = build_tables(ycoth_min=0.15, dyc=1e-4); L = Landscape(V_reg, dV_reg, T, C, Th_end)
    t = 1e-3; v0 = L.v0; HE = np.sqrt(t*v0/3); dq = HE/(2*np.pi)
    mu2 = -4*(3*c_star**2-4*c_star+8)/(c_star*(3*c_star+4)); s_lin = (-3+np.sqrt(9-4*mu2))/2
    R = dict(t=t, H_top_E_in_Mpl=HE, H_top_J_5Dunits=np.sqrt(L.f[L.i0])*HE, dTheta_quantum=dq, mu2=mu2, s_linear=s_lin, runs={})
    print("H_E(top)/M_pl = %.5e, dTheta = %.5e, s = %.5f" % (HE, dq, s_lin))
    saved = {}
    for side, sgn in (("plus(phi_b->+1)", -1), ("throat(phi_b->-1)", +1)):
        for k in (0, 1):
            th0 = sgn*dq; dth0 = s_lin*th0 if k == 0 else 0.0
            o = evolve(L, th0, dth0, k, 40.0 if sgn < 0 else 14.0, stop=(lambda s: s[3]/s[2] < -6.0))
            tau, th, a, h, u, K = o['tau'], o['Th'], o['a'], o['h'], o['u'], o['K']
            N = np.log(a/a[0]); r = {}
            # linear growth-rate check (k=0: pure growing mode)
            i1 = first(np.abs(th) > 1.5*dq); i2 = first(np.abs(th) > 6*dq)
            r['measured_growth_rate_s'] = float(np.log(th[i2]/th[i1])/(tau[i2]-tau[i1]))
            for name, cond in (("u_drop_1pct", u < 0.99*u[0]), ("u_drop_10pct", u < 0.9*u[0]), ("u_drop_50pct", u < 0.5*u[0]),
                               ("eps_H=1 (acceleration ends)", (o['dTh']**2/2 > h*h) ), ("|phi_b|=0.5", np.abs(o['phi']) > 0.5),
                               ("|phi_b|=0.9", np.abs(o['phi']) > 0.9)):
                i = first(cond)
                r[name] = None if i is None else dict(tau=float(tau[i]), N_E=float(N[i]), N_J=float(np.log(o['aJ'][i]/o['aJ'][0])), phi_b=float(o['phi'][i]), Theta=float(th[i]))
            ik = int(np.argmax(K)); r['peak_kinetic'] = dict(tau=float(tau[ik]), N_E=float(N[ik]), K_over_V0=float(K[ik]), K_over_rho=float(K[ik]/(K[ik]+u[ik])),
                                     u_over_V0=float(u[ik]), phi_b=float(o['phi'][ik]), eps_H=float(3*K[ik]/h[ik]**2), wmax=float(((K-u)/(K+u)).max()))
            r['max_constraint_violation'] = float(np.abs(o['cons']).max())
            # energy ledger at the end of the run: V0 - V = K + friction loss (k=0) ; all in units of V0
            r['ledger_end'] = dict(tau=float(tau[-1]), dV=float(u[0]-u[-1]), K=float(K[-1]), hubble_friction_loss=float(o['fric'][-1]), K0=float(K[0]),
                                   residual=float(u[0]+K[0]-u[-1]-K[-1]-o['fric'][-1]))
            # pulse duration: FWHM of d phi_b/d tau_J profile -> sech^2 tau_p = FWHM/1.7627  (Jordan time in units 1/H_top_J)
            tauJ = np.concatenate([[0], np.cumsum(0.5*(1/np.sqrt(o['f'][1:]/L.f[L.i0]) + 1/np.sqrt(o['f'][:-1]/L.f[L.i0]))*np.diff(tau))])
            o['tauJ'] = tauJ
            pd = np.abs(o['phidotJ_over_sqrt_t']); ip = int(np.argmax(pd))
            if sgn < 0:
                half = np.where(pd > 0.5*pd[ip])[0]; fwhm = tauJ[half[-1]] - tauJ[half[0]]
                r['phi_dot_pulse'] = dict(tauJ_peak=float(tauJ[ip]), FWHM_in_inv_HtopJ=float(fwhm), tau_p_sech2=float(fwhm/1.7627),
                                          peak_phidot_over_H_topJ=float(pd[ip]/np.sqrt(L.f[L.i0]*v0/3)))
                hh = np.abs(np.gradient(o['hJ'], tauJ)); ih = int(np.argmax(hh)); half = np.where(hh > 0.5*hh[ih])[0]
                r['Hdot_J_pulse'] = dict(tauJ_peak=float(tauJ[ih]), FWHM=float(tauJ[half[-1]]-tauJ[half[0]]), tau_p_sech2=float((tauJ[half[-1]]-tauJ[half[0]])/1.7627),
                                         peak_minus_HdotJ_over_HtopJ2=float(hh[ih]))
                r['hJ_min_over_HtopJ'] = float(o['hJ'].min()); r['hJ_final_over_HtopJ'] = float(o['hJ'][-1]); r['hJ_final_exact_RS'] = float(np.sqrt(((1+c_star)/27)/(L.f[L.i0]*v0/3)))
                r['end_state'] = dict(tau=float(tau[-1]), Theta_minus_Theta_end=float(th[-1]-L.Th_plus_end), y_b=float(9*np.log((th[-1]-L.Th_plus_end)/L.dTh_min)+L.y[0]) ,
                                      h_E_final=float(h[-1]), h_E_expected=float(np.sqrt(((1+c_star)/81)/v0)),
                                      dlnDeltaTheta_dtau_over_h=float(o['dTh'][-1]/(th[-1]-L.Th_plus_end)/h[-1]),
                                      sinhpsi_over_sqrt_t_final=float(o['sinhpsi_over_sqrt_t'][-1]), minus_l_HJ_over_sqrt_t_exact=float(-9*np.sqrt((1+c_star)/27)),
                                      max_abs_sinhpsi_over_sqrt_t=float(np.abs(o['sinhpsi_over_sqrt_t']).max()),
                                      min_DeltaTheta_reached=float((th-L.Th_plus_end).min()), reached_boundary=bool((th-L.Th_plus_end).min() <= 0))
            else:
                sp = np.abs(o['sinhpsi_over_sqrt_t'])*np.sqrt(t)
                for thr in (0.1, 0.3, 1.0):
                    i = first(sp > thr)
                    r['naive_rapidity |sinh psi|>%.1f (t=1e-3)' % thr] = None if i is None else dict(tau=float(tau[i]), phi_b=float(o['phi'][i]), one_plus_phi=float(1+o['phi'][i]), y_b=float(o['y'][i]), N_E=float(N[i]))
                ic = first(th >= L.Th_throat)
                r['crossing phi_b=-1'] = dict(tau=float(tau[ic]), N_E=float(N[ic]), N_J=float(np.log(o['aJ'][ic]/o['aJ'][0])), dTheta_dtau=float(o['dTh'][ic]), K_over_V0=float(K[ic]), u_over_V0=float(u[ic]),
                                              h_E=float(h[ic]), hJ_over_HtopJ=float(o['hJ'][ic]),
                                              phidot_J_5Dunits_t1em3=float(o['phidotJ_over_sqrt_t'][ic]*np.sqrt(t)), HJ_times_l_minus_t1em3=float(o['hJ'][ic]*np.sqrt(L.f[L.i0])*HE*9/5))
                iz = first(u <= 0); r['V_E=0 crossing (phi_b=-1/c)'] = None if iz is None else dict(tau=float(tau[iz]), N_E=float(N[iz]), phi_b=float(o['phi'][iz]))
                it = first(h <= 0); r['Einstein-frame turnaround h=0'] = None if it is None else dict(tau=float(tau[it]), N_E=float(N[it]), phi_b=float(o['phi'][it]), a_over_a0=float(a[it]/a[0]))
                itj = first(o['hJ'] <= 0); r['Jordan-frame turnaround H_J=0'] = None if itj is None else dict(tau=float(tau[itj]), phi_b=float(o['phi'][itj]))
                r['end_of_run'] = dict(tau=float(tau[-1]), phi_b=float(o['phi'][-1]), h=float(h[-1]), K_over_V0=float(K[-1]), phidotJ_over_sqrt_t=float(o["phidotJ_over_sqrt_t"][-1]), hJ_over_HtopJ=float(o["hJ"][-1]), reason="table end or h<-6")
            R['runs']["%s k=%d" % (side, k)] = r
            saved["%s_k%d" % ("plus" if sgn < 0 else "throat", k)] = o
            print(side, "k=", k); print(json.dumps(r, indent=1))
    # e-folds before roll-off vs t (k=0): classical dynamics is t-independent -> N shifts by ln(ratio)/s
    Ns = {}
    for tt in (1e-2, 1e-3, 1e-6, 1e-9, 1e-12):
        d = np.sqrt(tt*v0/3)/(2*np.pi); o = evolve(L, -d, -s_lin*d, 0, 60.0, dtau=1e-3, stop=lambda s: s[0] < -1.0)
        i = first(o['u'] < 0.9*o['u'][0]); Ns["t=%g" % tt] = dict(H_over_Mpl=float(np.sqrt(tt*v0/3)), N_E_to_10pct_drop=float(np.log(o['a'][i]/o['a'][0])), estimate_ln_ratio_over_s=float(np.log(0.2/d)/s_lin))
    R['efolds_before_rolloff_vs_t (plus side, k=0)'] = Ns; print(json.dumps(Ns, indent=1))
    json.dump(R, open("BLAST_HOMOGENEOUS_RESULTS.json", "w"), indent=1)
    for kx, o in saved.items():
        np.savez("traj_%s.npz" % kx, **{a_: b_ for a_, b_ in o.items() if isinstance(b_, np.ndarray)})
