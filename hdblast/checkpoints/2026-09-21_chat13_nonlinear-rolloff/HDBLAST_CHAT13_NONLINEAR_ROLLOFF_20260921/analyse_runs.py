#!/usr/bin/env python3
"""Chat 13 - post-processing of the nonlinear runs: linear growth rate, end states, crunch time, comparison of resolutions."""
import json, glob, math, os
import numpy as np
P = json.load(open("runs/linear_predictions.json"))
out = {}
for fn in sorted(glob.glob("runs/*_timeseries.npz")):
    tag = os.path.basename(fn).replace("_timeseries.npz", "")
    if tag.startswith(("smoke", "diag")): continue
    D = np.load(fn); rec = D["rec"]; summ = json.load(open(fn.replace("_timeseries.npz", "_summary.json")))
    t, phib, HJ, bb = rec[:, 0], rec[:, 1], rec[:, 2], rec[:, 3]
    tau = rec[:, 9] if rec.shape[1] > 9 else None
    pb0 = summ["phi_b_static"]; dphi = np.abs(phib - pb0)
    r = {"t_det": summ["t_det"], "dz": summ["dz"], "dc": summ["dc"], "stop": summ["stop_reason"]}
    m = (dphi > 5*dphi[0]) & (dphi < 5e-3) & (t > 0.5)
    if m.sum() > 20:
        sl = np.polyfit(t[m], np.log(dphi[m]), 1)[0]; r["growth_rate_fit"] = float(sl)
        key = str(summ["t_det"]); r["growth_rate_linear_theory"] = P[key]["growth_rate"] if key in P else None
    if summ["dc"] > 0:
        i = np.argmax(phib); r["phi_b_max"] = float(phib[i]); r["t_at_phi_max"] = float(t[i])
        late = t > (t[i] + 1.0)
        if late.any(): r["HJ_over_H0_plateau"] = [float(HJ[late].min()), float(HJ[late].max())]
        Ip = 1.0357712571566784; c = 2/Ip - 4/3
        sig1 = 2/3 + summ["t_det"]*(1 + c); r["HJ_over_H0_dS_brane_in_AdS_plus_exact"] = float(math.sqrt(sig1**2/36 - 1/81)*summ["rho_b"])   # exact RS brane in AdS_+ (referee correction: earlier version dropped the t^2 term)
        if tau is not None: r["H0tau_at_phi_0p9"] = float(tau[np.argmax(phib > 0.9)]) if (phib > 0.9).any() else None
    else:
        j = np.where(HJ < 0)[0]
        if j.size:
            j0 = j[0]; r["turnaround"] = {"t": float(t[j0]), "H0tau": float(tau[j0]) if tau is not None else None, "phi_b": float(phib[j0])}
            if tau is not None:
                # crunch: fit H_J ~ -alpha/(tau* - tau) on the last part
                k = HJ < -3
                if k.sum() > 5:
                    x = tau[k]; yv = -1/HJ[k]                                   # = (tau* - tau)/alpha  -> linear in tau
                    A_ = np.polyfit(x, yv, 1); r["crunch"] = {"H0_tau_star": float(-A_[1]/A_[0]), "alpha": float(-1/A_[0]), "n_points": int(k.sum()), "last_H0tau": float(tau[-1]), "last_HJ": float(HJ[-1]), "last_phi_b": float(phib[-1])}
        k1 = np.argmax(phib < -1) if (phib < -1).any() else None
        if k1: r["cross_phi_-1"] = {"t": float(t[k1]), "H0tau": float(tau[k1]) if tau is not None else None}
    out[tag] = r; print(tag, json.dumps(r))
json.dump(out, open("runs/ANALYSIS.json", "w"), indent=1)
