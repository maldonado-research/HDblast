#!/usr/bin/env python3
"""Chat 14 - analysis of the v2 runs: growth-rate fits on sliding windows (compare with the 5D linear theory), end states, turnaround."""
import json, glob, os, math
import numpy as np
S_LIN = {0.001: 1.6571907, 0.03: 1.64836, 0.1: 1.62702}      # growth rates from the Chat 9 linear theory (mu^2 = -7.7178716 at t = 1e-3)
Ip = 1.0357712571566784; c = 2/Ip - 4/3
out = {}
for fn in sorted(glob.glob("runs/*_timeseries.npz")):
    tag = os.path.basename(fn).replace("_timeseries.npz", "")
    if tag.startswith("diag"): continue
    summ = json.load(open(fn.replace("_timeseries.npz", "_summary.json"))) if os.path.exists(fn.replace("_timeseries.npz", "_summary.json")) else None
    rec = np.load(fn)["rec"]; t, phib, HJ, bb, tau = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3], rec[:, 9]
    td = summ["t_det"] if summ else None; pb0 = summ["phi_b_static"] if summ else phib[0]; rb = summ["rho_b"] if summ else None
    r = dict(t_det=td, dc=summ["dc"] if summ else None, points_across_wall=summ.get("points_across_wall") if summ else None, closure=summ.get("closure") if summ else None, t_end=float(t[-1]))
    d = np.abs(phib - pb0); ok = d > 0
    # sliding-window log-slope of |phi_b - phi_b0| where the deviation is still small (< 1e-3): local growth rate
    rates = []
    for a in np.arange(1.0, t[-1] - 1.0, 0.5):
        m = (t >= a) & (t < a + 1.0) & ok & (d < 1e-3) & (d > 10*d[0] if d[0] > 0 else d > 0)
        if m.sum() > 20: rates.append((round(float(a), 1), round(float(np.polyfit(t[m], np.log(d[m]), 1)[0]), 4)))
    r["local_growth_rates"] = rates; r["linear_theory_rate"] = S_LIN.get(td)
    if rates: r["late_linear_rate"] = rates[-1][1]
    if summ and summ["dc"] > 0:
        i = np.argmax(phib); r["phi_b_max"] = float(phib[i]); r["H0tau_end"] = float(tau[-1]); r["HJ_over_H0_end"] = float(HJ[-1])
        sig1 = 2/3 + td*(1 + c); r["HJ_over_H0_RS_exact"] = float(math.sqrt(sig1**2/36 - 1/81)*rb)
        k = phib > 0.98
        if k.any(): j = np.argmax(k); r["at_phi_0.98"] = dict(H0tau=float(tau[j]), HJ_over_H0=float(HJ[j]))
    elif summ and summ["dc"] < 0:
        j = np.where(HJ < 0)[0]
        if j.size: r["turnaround"] = dict(H0tau=float(tau[j[0]]), phi_b=float(phib[j[0]]), t=float(t[j[0]]))
        k = phib < -1
        if k.any(): j1 = np.argmax(k); r["cross_phi_-1"] = dict(H0tau=float(tau[j1]), HJ_over_H0=float(HJ[j1]))
        k5 = HJ < -5
        if k5.any(): r["H0tau_at_HJ_-5"] = float(tau[np.argmax(k5)])
    out[tag] = r; print(tag, json.dumps(r))
json.dump(out, open("runs/ANALYSIS_V2.json", "w"), indent=1)
