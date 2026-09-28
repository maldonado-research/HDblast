#!/usr/bin/env python3
"""Spectator particle production along the actual roll-off trajectory, side (i) phi_b -> +1, k=0, brane (Jordan) frame.
Field chi_s minimally coupled to the induced metric; v = a^{3/2} chi:  v'' + Om^2 v = 0,
   Om^2 = k^2/a_J^2 + m_eff^2 - (9/4) H_J^2 - (3/2) dH_J/dtau_J        (units: H_topJ = 1, a_J(start) = 1)
Bogoliubov ODEs in the WKB basis:  alpha' = (W'/2W) e^{+2i th} beta,  beta' = (W'/2W) e^{-2i th} alpha,  th' = W.
Particle number = first-order adiabatic number |beta - i W'/(4W^2) e^{-2i th} alpha|^2 (in and out states in de Sitter)."""
import numpy as np, json, sys, time
from blast_common import *
from blast_homogeneous import evolve, first

def bogoliubov(tgrid, A, M, ks, eps=0.05, dtmax=0.01, t_start=None, t_end=None, trace=False):
    """tgrid uniform-ish background grid; A=1/a^2, M = m_eff^2-(9/4)H^2-(3/2)Hdot sampled on tgrid. ks: array of comoving k."""
    Ad = np.gradient(A, tgrid); Md = np.gradient(M, tgrid)
    kmax = ks.max(); t = t_start; nodes = [t]
    while t < t_end:
        Wm = np.sqrt(kmax**2*np.interp(t, tgrid, A) + np.interp(t, tgrid, M)); t = min(t + min(dtmax, eps/Wm), t_end); nodes.append(t)
    nodes = np.array(nodes); mids = 0.5*(nodes[1:]+nodes[:-1])
    def WWd(tt):
        a_, m_, ad_, md_ = (np.interp(tt, tgrid, q) for q in (A, M, Ad, Md))
        W2 = ks[None, :]**2*a_[:, None] + m_[:, None]
        if W2.min() <= 0: raise ValueError("Omega^2 <= 0: adiabatic basis undefined")
        Wv = np.sqrt(W2); return Wv, (ks[None, :]**2*ad_[:, None] + md_[:, None])/(2*Wv)
    Wn, Wdn = WWd(nodes); Wmid, Wdm = WWd(mids)
    gn = Wdn/(2*Wn); gm = Wdm/(2*Wmid)
    al = np.ones(len(ks), complex); be = 1j*Wdn[0]/(4*Wn[0]**2)*al; th = np.zeros(len(ks))
    tr = []
    for i in range(len(nodes)-1):
        h = nodes[i+1]-nodes[i]
        thm = th + h*(5*Wn[i] + 8*Wmid[i] - Wn[i+1])/24; th1 = th + h*(Wn[i] + 4*Wmid[i] + Wn[i+1])/6
        e0 = gn[i]*np.exp(2j*th); em = gm[i]*np.exp(2j*thm); e1 = gn[i+1]*np.exp(2j*th1)
        k1a = e0*be;              k1b = np.conj(e0)*al
        k2a = em*(be+h/2*k1b);    k2b = np.conj(em)*(al+h/2*k1a)
        k3a = em*(be+h/2*k2b);    k3b = np.conj(em)*(al+h/2*k2a)
        k4a = e1*(be+h*k3b);      k4b = np.conj(e1)*(al+h*k3a)
        al = al + h/6*(k1a+2*k2a+2*k3a+k4a); be = be + h/6*(k1b+2*k2b+2*k3b+k4b); th = th1
        if trace and i % 200 == 0:
            tr.append((nodes[i+1], np.abs(be - 1j*Wdn[i+1]/(4*Wn[i+1]**2)*np.exp(-2j*th)*al)**2))
    N = np.abs(be - 1j*Wdn[-1]/(4*Wn[-1]**2)*np.exp(-2j*th)*al)**2
    return N, np.abs(al)**2 - np.abs(be)**2 - 1, tr, len(nodes)

def selftest():
    """flat-space exact benchmarks: sech^2 pulse (Chat 8, eq. 6) and tanh step (Bernard-Duncan / Birrell-Davies 3.4)."""
    tg = np.linspace(-25, 25, 200001); out = {}
    for lam, kap in ((0.5, 0.5), (0.5, 1.0), (1.0, 0.7), (1.5, 0.3)):
        M = lam*(lam+1)/np.cosh(tg)**2; N, un, _, _ = bogoliubov(tg, np.ones_like(tg), M, np.array([kap]), eps=0.02, dtmax=0.005, t_start=-25, t_end=25)
        out["sech2 lam=%.1f kappa=%.1f" % (lam, kap)] = dict(numeric=float(N[0]), exact=float(np.sin(np.pi*lam)**2/np.sinh(np.pi*kap)**2), unitarity=float(un[0]))
    Ai, Bi, rho = 2.0, 1.0, 1.0; M = Ai + Bi*np.tanh(rho*tg); kk = 0.5
    N, un, _, _ = bogoliubov(tg, np.ones_like(tg), M, np.array([kk]), eps=0.02, dtmax=0.005, t_start=-25, t_end=25)
    wi, wo = np.sqrt(kk*kk+Ai-Bi), np.sqrt(kk*kk+Ai+Bi)
    out["tanh step"] = dict(numeric=float(N[0]), exact=float(np.sinh(np.pi*(wo-wi)/(2*rho))**2/(np.sinh(np.pi*wi/rho)*np.sinh(np.pi*wo/rho))), unitarity=float(un[0]))
    return out

if __name__ == "__main__":
    R = dict(selftest=selftest()); print(json.dumps(R['selftest'], indent=1)); sys.stdout.flush()
    T, C, Th_end = build_tables(); L = Landscape(V_reg, dV_reg, T, C, Th_end)
    t = 1e-3; v0 = L.v0; HE = np.sqrt(t*v0/3); dq = HE/(2*np.pi)
    mu2 = -4*(3*c_star**2-4*c_star+8)/(c_star*(3*c_star+4)); s_lin = (-3+np.sqrt(9-4*mu2))/2
    Nearly = 9.0                                         # extra e-folds of hilltop de Sitter before the registered start (same growing-mode trajectory)
    d0 = dq*np.exp(-s_lin*Nearly)
    o = evolve(L, -d0, -s_lin*d0, 0, Nearly+42.0, dtau=5e-4)
    f0 = L.f[L.i0]; tau = o['tau']; fr = o['f']/f0
    tJ = np.concatenate([[0], np.cumsum(0.5*(fr[1:]**-0.5 + fr[:-1]**-0.5)*np.diff(tau))])
    aJ = o['aJ']/o['aJ'][0]; hJ = o['hJ']; hJd = np.gradient(hJ, tJ); phi = o['phi']
    phid = np.gradient(phi, tJ); ip = int(np.argmax(np.abs(phid))); t_star = tJ[ip]; a_star = aJ[ip]
    ih = int(np.argmax(-hJd)); half = np.where(-hJd > 0.5*(-hJd[ih]))[0]; tau_p_H = (tJ[half[-1]]-tJ[half[0]])/1.7627
    half = np.where(np.abs(phid) > 0.5*np.abs(phid[ip]))[0]; tau_p_phi = (tJ[half[-1]]-tJ[half[0]])/1.7627
    R['background'] = dict(tJ_star=float(t_star), N_J_star=float(np.log(a_star)), phi_b_star=float(phi[ip]), phidot_star_over_H=float(phid[ip]), HJ_star=float(hJ[ip]),
                           tau_p_phidot=float(tau_p_phi), tau_p_Hdot=float(tau_p_H), peak_minus_Hdot=float(-hJd[ih]), HJ_final=float(hJ[-1]), tJ_end=float(tJ[-1]),
                           sech2_ceiling_rho_in_H4=dict(phidot=float(9.3775e-4/tau_p_phi**4), Hdot=float(9.3775e-4/tau_p_H**4)))
    print(json.dumps(R['background'], indent=1)); sys.stdout.flush()
    A = 1/aJ**2; Mgeo = -(9/4)*hJ**2 - 1.5*hJd
    i_dS = first(o['u'] < 0.999*o['u'][0]); t_dS = tJ[i_dS]          # last time at which the background is still hilltop de Sitter to 1e-3
    kap = np.exp(np.linspace(np.log(0.02), np.log(12.0), 40))        # k/(a_star m)
    def run(m, g=0.0, phis=0.0, label=""):
        M = m*m + g*g*(phi-phis)**2 + Mgeo
        Nk = np.empty(len(kap)); un = np.empty(len(kap)); nst = 0; t0 = time.time()
        meff0 = np.sqrt(m*m + g*g*phis**2)
        for b in range(0, len(kap), 8):
            ks = kap[b:b+8]*a_star*m
            a_s = min(ks.min()/(100*meff0), aJ[i_dS]); ts = np.interp(a_s, aJ, tJ)
            Nb, ub, _, nn = bogoliubov(tJ, A, M, ks, t_start=max(ts, tJ[5]), t_end=tJ[-5]); Nk[b:b+8] = Nb; un[b:b+8] = ub; nst += nn
        mu_i = np.sqrt(meff0**2 - 2.25); mf2 = m*m + g*g*(1-phis)**2; mu_f = np.sqrt(mf2/hJ[-1]**2 - 2.25)
        base_i = 1/(np.exp(2*np.pi*mu_i)-1); base_f = 1/(np.exp(2*np.pi*mu_f)-1)
        k = kap*a_star*m; ex = Nk - base_f
        n_star = np.trapz(k**3*ex, np.log(k))/(2*np.pi**2*a_star**3)
        om = np.sqrt(k**2/a_star**2 + mf2); rho_star = np.trapz(k**3*om*ex, np.log(k))/(2*np.pi**2*a_star**3)
        r = dict(m=m, g=g, phi_star=phis, N_k=[float(x) for x in Nk], kappa=[float(x) for x in kap], unitarity_max=float(np.abs(un).max()),
                 dS_baseline_initial=float(base_i), dS_baseline_final=float(base_f), N_IR=float(Nk[0]), N_UV=float(Nk[-1]), N_max=float(Nk.max()), kappa_at_max=float(kap[int(np.argmax(Nk))]),
                 n_excess_at_a_star_in_H3=float(n_star), rho_excess_at_a_star_in_H4=float(rho_star), steps=nst, seconds=time.time()-t0)
        # sech^2 benchmark with sin^2 = 1 (envelope) and fitted tau_p
        for nm, tp in (("phidot", tau_p_phi), ("Hdot", tau_p_H)):
            r["sech2_envelope_N_k_tau_p_"+nm] = [float(1/np.sinh(np.pi*tp*np.sqrt((q*m)**2 + mf2))**2) for q in kap[::8]]
        if g > 0 and 0 < phis < 1:
            j = first(phi > phis); pd = abs(phid[j]); aj = aJ[j]
            r['KLS_estimate'] = dict(phidot_cross=float(pd), g_phidot=float(g*pd), N_k=[float(np.exp(-np.pi*((q*m*a_star/aj)**2 + m*m - 2.25*hJ[j]**2 - 1.5*hJd[j])/(g*pd))) for q in kap[::8]],
                                     N_k_numeric=[float(x) for x in Nk[::8]])
        R[label] = r
        print(label, json.dumps({k_: v_ for k_, v_ in r.items() if k_ not in ('N_k', 'kappa')}, indent=1)); sys.stdout.flush()
    for m in (1.6, 2.0, 3.0): run(m, label="(a) minimal m=%.1f H_topJ" % m)
    run(2.0, g=10.0, phis=0.5, label="(b) m=2, g=10 H_topJ, phi_star=0.5")
    run(2.0, g=30.0, phis=0.5, label="(b) m=2, g=30 H_topJ, phi_star=0.5")
    run(2.0, g=5.0, phis=1.0, label="(b) m=2, g=5 H_topJ, phi_star=phi_inf=1 (Chat-8 convention)")
    # ratio to the total energy: rho_tot = 3 M_pl^2 H^2 -> rho_prod/rho_tot = (C/3) (H_J/M_pl,J)^2 with rho_prod = C H_J^4 ; (H/M_pl)^2 = t v0/3
    R['rho_ratio_formula'] = "rho_prod/rho_tot = (C/3) (H/M_pl)^2,  C = rho_excess_in_H4,  (H/M_pl)^2 = t*v0/3 = %.4e at t=1e-3" % (t*v0/3)
    json.dump(R, open("BLAST_PARTICLE_RESULTS.json", "w"), indent=1)
