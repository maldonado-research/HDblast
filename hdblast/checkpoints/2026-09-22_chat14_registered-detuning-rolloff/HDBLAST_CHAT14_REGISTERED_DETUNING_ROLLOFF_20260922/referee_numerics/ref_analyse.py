#!/usr/bin/env python3
"""Referee (numerics) re-analysis of the Chat 14 saved time series. Independent of analyse_v2.py.
rec columns: [t, phi_b, H_J/H0, b_b, a_t_b, dphi/dtau*rho_b, M_bulk, max|f|, max|b|, H0*tau, phi_min, phi_max, M_shell]"""
import json, glob, os, sys, math
import numpy as np
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS = sys.argv[1:] or [os.path.join(BASE, "runs")]
S_LIN = 1.6571907
out = {}
def load(fn):
    d = np.load(fn); rec = d["rec"]
    summ_fn = fn.replace("_timeseries.npz", "_summary.json")
    summ = json.load(open(summ_fn)) if os.path.exists(summ_fn) else {}
    return d, rec, summ
def interp_at(x, y, x0):
    """linear interpolation of y at x0 on a monotone x (first crossing)"""
    k = np.where((x[:-1] - x0)*(x[1:] - x0) <= 0)[0]
    if k.size == 0: return None
    i = k[0]; w = (x0 - x[i])/(x[i+1] - x[i]); return float(y[i] + w*(y[i+1] - y[i]))
for rd in RUNS:
    for fn in sorted(glob.glob(os.path.join(rd, "*_timeseries.npz"))):
        tag = os.path.basename(fn).replace("_timeseries.npz", "")
        d, rec, summ = load(fn)
        t, phib, HJ, bb, pab, M_bulk, tau, M_shell = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3], rec[:, 4], rec[:, 6], rec[:, 9], rec[:, 12]
        pb0 = summ.get("phi_b_static", phib[0]); dc = summ.get("dc", 0.0)
        r = dict(closure=summ.get("closure"), dc=dc, pts_wall=summ.get("points_across_wall"), t_end=float(t[-1]), H0tau_end=float(tau[-1]))
        dev = phib - pb0
        # (1) growth rates: sliding 1-Hubble-time windows of the log of |phi_b - phi_b0| restricted to the linear regime |dev| < 1e-3
        rates = {}
        for a in np.arange(1.0, min(t[-1], 9.0) - 0.99, 0.5):
            m = (t >= a) & (t < a + 1) & (np.abs(dev) < 1e-3) & (np.abs(dev) > 0)
            if m.sum() > 20:
                sl = np.polyfit(t[m], np.log(np.abs(dev[m])), 1)[0]; rates["%.1f-%.1f" % (a, a + 1)] = round(float(sl), 5)
        r["window_rates"] = rates
        # also a rate from the metric field b_b (independent field) and from a_t_b
        rates_b = {}
        for a in np.arange(2.0, min(t[-1], 9.0) - 0.99, 1.0):
            m = (t >= a) & (t < a + 1) & (np.abs(dev) < 1e-3) & (np.abs(bb) > 0)
            if m.sum() > 20:
                sl = np.polyfit(t[m], np.log(np.abs(bb[m])), 1)[0]; rates_b["%.1f-%.1f" % (a, a + 1)] = round(float(sl), 4)
        r["window_rates_b_b"] = rates_b
        # (2) +1 side: values at fixed proper time H0 tau = 6.9 and at 6.5, 6.7, and the last point; b_b (lapse) history
        if dc is not None and dc > 0:
            for tau0 in (6.0, 6.5, 6.7, 6.8, 6.9, 7.0):
                v = interp_at(tau, HJ, tau0); p = interp_at(tau, phib, tau0); tb = interp_at(tau, t, tau0); lb = interp_at(tau, bb, tau0)
                if v is not None: r["at_H0tau_%.1f" % tau0] = dict(HJ_over_H0=round(v, 5), phi_b=round(p, 5), t=round(tb, 3), b_b=round(lb, 3))
            # H_J vs tau slope over the last 0.2 in tau before 6.9
            m = (tau > 6.7) & (tau < 6.9)
            if m.sum() > 5: r["dHJ_dtau_6.7_6.9"] = round(float(np.polyfit(tau[m], HJ[m], 1)[0]), 4)
            j = np.where(HJ < 0)[0]
            if j.size: r["first_negative_HJ"] = dict(t=float(t[j[0]]), H0tau=float(tau[j[0]]), b_b=float(bb[j[0]]), a_t_b=float(pab[j[0]]))
            r["HJ_end"] = float(HJ[-1]); r["phib_end"] = float(phib[-1]); r["b_b_end"] = float(bb[-1]); r["a_t_b_end"] = float(pab[-1])
            r["H0tau_at_phib_0.5"] = interp_at(phib, tau, 0.5); r["t_at_phib_0.5"] = interp_at(phib, t, 0.5)
            r["H0tau_at_phib_0.98"] = interp_at(phib, tau, 0.98)
        # (3) throat side
        if dc is not None and dc < 0:
            j = np.where(HJ < 0)[0]
            if j.size:
                i = j[0]; r["turnaround"] = dict(H0tau=round(interp_at(HJ[i-1:i+1], tau[i-1:i+1], 0.0), 5), phi_b=round(interp_at(HJ[i-1:i+1], phib[i-1:i+1], 0.0), 5), t=round(float(t[i]), 4))
            r["cross_phi_-1"] = dict(H0tau=interp_at(phib, tau, -1.0), HJ_over_H0=interp_at(phib, HJ, -1.0), t=interp_at(phib, t, -1.0))
            for thr in (-5.0, -10.0, -30.0, -100.0, -300.0):
                v = interp_at(HJ, tau, thr)
                if v is not None: r["H0tau_at_HJ_%g" % thr] = round(v, 4)
            r["H0tau_at_phib_-1.67"] = interp_at(phib, tau, -1.0/0.5975949350280)
            r["HJ_end"] = float(HJ[-1]); r["phib_end"] = float(phib[-1])
        # (4) constraint monitors: max over time of shell and bulk monitors before/after key times
        ok = np.isfinite(M_shell)
        r["M_shell_max_all"] = float(np.nanmax(M_shell)) if ok.any() else None
        for tcut in (4.0, 6.0, 7.0, 8.0):
            m = ok & (t <= tcut)
            if m.any(): r["M_shell_max_t<=%g" % tcut] = float(np.max(M_shell[m])); r["M_bulk_max_t<=%g" % tcut] = float(np.max(M_bulk[m]))
        if dc is not None and dc > 0:
            m = ok & (tau <= 6.9)
            if m.any(): r["M_shell_max_H0tau<=6.9"] = float(np.max(M_shell[m])); r["M_bulk_max_H0tau<=6.9"] = float(np.max(M_bulk[m]))
        # (5) end-state profile: kink location and grid-scale noise
        if "f_end" in d and "phis" in d and "z" in d:
            z, phi = d["z"], d["phis"] + d["f_end"]; f = d["f_end"]
            k = np.where((phi[:-1] - 0.0)*(phi[1:] - 0.0) <= 0)[0]
            r["phi_zero_crossings_z"] = [round(float(z[i]), 4) for i in k[-5:]]
            r["phi_range_end"] = [float(phi.min()), float(phi.max())]
            # grid-scale noise indicator: rms of the 2nd difference relative to rms of field over the fine region and the coarse region
            dd = f[2:] - 2*f[1:-1] + f[:-2]
            fine = z[1:-1] > -0.06; coarse = z[1:-1] < -0.2
            r["noise_fine_rms_d2f_over_rms_f"] = float(np.sqrt(np.mean(dd[fine]**2))/(np.sqrt(np.mean(f[1:-1][fine]**2)) + 1e-300))
            r["noise_coarse_rms_d2f_over_rms_f"] = float(np.sqrt(np.mean(dd[coarse]**2))/(np.sqrt(np.mean(f[1:-1][coarse]**2)) + 1e-300))
            # location of largest |d2 f| (grid-scale activity)
            i = int(np.argmax(np.abs(dd))); r["max_d2f_at_z"] = float(z[i+1]); r["max_d2f"] = float(np.abs(dd[i]))
        out[tag] = r
        print(tag, json.dumps(r, indent=None))
json.dump(out, open(os.path.join(BASE, "referee_numerics", "REFEREE_ANALYSIS.json"), "w"), indent=1)
