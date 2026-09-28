#!/usr/bin/env python3
"""Referee post-processing: compares the referee runs with the author's runs (../runs) at equal conformal time t."""
import json, os, glob, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(HERE, "..", "runs"); R = os.path.join(HERE, "runs")
def load(path):
    D = np.load(path); return D["rec"]
def at(rec, T, col):
    i = np.argmin(np.abs(rec[:, 0] - T)); return rec[i, col] if abs(rec[i, 0] - T) < 0.02 else float("nan")
out = {}
S = 1.6270213424531335            # linear growth rate at t_det = 0.1 (author's linear_predictions.json)
# ---- (b) seed sign flipped, half amplitude: expected time shift ln2/s for every event
ref = load(os.path.join(A, "t01_minus_coarse_timeseries.npz"))
if os.path.exists(os.path.join(R, "ref_minus_half_timeseries.npz")):
    h = load(os.path.join(R, "ref_minus_half_timeseries.npz"))
    def first_cross(rec, val):
        j = np.where(rec[:, 1] < val)[0]; return rec[j[0], 0] if j.size else float("nan")
    def turnaround(rec):
        j = np.where(rec[:, 2] < 0)[0]; return (rec[j[0], 0], rec[j[0], 9], rec[j[0], 1]) if j.size else (float("nan"),)*3
    res = {"expected_shift_ln2_over_s": math.log(2)/S}
    for name, val in [("cross_phi_-0.5", -0.5), ("cross_phi_-1", -1.0), ("cross_phi_-1.5", -1.5)]:
        res[name] = {"author_dc=-1e-4": first_cross(ref, val), "referee_dc=-5e-5": first_cross(h, val), "shift": first_cross(h, val) - first_cross(ref, val)}
    ta, th = turnaround(ref), turnaround(h)
    res["turnaround"] = {"author": {"t": ta[0], "H0tau": ta[1], "phi_b": ta[2]}, "referee_half_seed": {"t": th[0], "H0tau": th[1], "phi_b": th[2]}, "shift_t": th[0] - ta[0]}
    # invariance check: the trajectory in (phi_b, H_J) space should be seed-independent (to O(seed^2)); compare H_J at equal phi_b
    for pv in [-1.0, -1.5, -1.9, -3.0, -5.0]:
        ia = np.argmax(ref[:, 1] < pv) if (ref[:, 1] < pv).any() else None; ih = np.argmax(h[:, 1] < pv) if (h[:, 1] < pv).any() else None
        if ia and ih: res["HJ_at_phi_b=%g" % pv] = {"author": float(ref[ia, 2]), "referee": float(h[ih, 2])}
    out["seed_test"] = res
# ---- (c) KO strength, (pslope) closure variant, (d) L = 16: compare with author's dz = 2e-3 plus run at equal t
base = load(os.path.join(A, "t01_plus_coarse_timeseries.npz"))
cmp = {}
for tag in ["ref_plus_ko002", "ref_plus_ko01", "ref_plus_pslope", "ref_plus_L16"]:
    fn = os.path.join(R, tag + "_timeseries.npz")
    if not os.path.exists(fn): continue
    r = load(fn); d = {}
    for T in [4, 6, 8, 10, 11, 12, 13, 14, 15, 16]:
        if T > r[-1, 0] + 0.02: continue
        d["t=%g" % T] = {"phi_b": float(at(r, T, 1)), "HJ/H0": float(at(r, T, 2)), "H0tau": float(at(r, T, 9)),
                         "d_phi_b_vs_author": float(at(r, T, 1) - at(base, T, 1)), "d_HJ_vs_author": float(at(r, T, 2) - at(base, T, 2)),
                         "M_int": float(np.nanmax(r[(r[:, 0] > T - 0.5) & (r[:, 0] <= T + 0.01), 6])), "M_shell": float(np.nanmax(r[(r[:, 0] > T - 0.5) & (r[:, 0] <= T + 0.01), 12])),
                         "H_int": float(np.nanmax(r[(r[:, 0] > T - 0.5) & (r[:, 0] <= T + 0.01), 13])), "H_shell": float(np.nanmax(r[(r[:, 0] > T - 0.5) & (r[:, 0] <= T + 0.01), 14]))}
    cmp[tag] = d
out["plus_variants_vs_author_dz2e-3"] = cmp
basem = load(os.path.join(A, "t01_minus_coarse_timeseries.npz")); cm = {}
for tag in ["ref_minus_ko01", "ref_minus_pslope"]:
    fn = os.path.join(R, tag + "_timeseries.npz")
    if not os.path.exists(fn): continue
    r = load(fn); d = {}
    j = np.where(r[:, 2] < 0)[0]
    d["turnaround"] = {"t": float(r[j[0], 0]), "H0tau": float(r[j[0], 9]), "phi_b": float(r[j[0], 1])} if j.size else None
    for T in [5.5, 6.0, 6.25, 6.5, 7.0, 8.0, 9.0]:
        if T > r[-1, 0] + 0.02: continue
        d["t=%g" % T] = {"phi_b": float(at(r, T, 1)), "HJ/H0": float(at(r, T, 2)), "H0tau": float(at(r, T, 9)), "a_b": float(at(r, T, 15)), "b_b": float(at(r, T, 3)),
                         "author_phi_b": float(at(basem, T, 1)), "author_HJ": float(at(basem, T, 2)), "M_shell": float(np.nanmax(r[(r[:, 0] > T - 0.5) & (r[:, 0] <= T + 0.01), 12]))}
    cm[tag] = d
out["minus_variants_vs_author_dz2e-3"] = cm
# ---- t = 1e-3 resolution scan: growth-rate fit of |phi_b - phi_b_static| on 1 < t < 3
sc = {}
for tag in ["ref_t1e3_dz1e-3", "ref_t1e3_dz5e-4", "ref_t1e3_dz2.5e-4", "ref_t1e3_dz5e-4_pslope"]:
    fn = os.path.join(R, tag + "_timeseries.npz")
    if not os.path.exists(fn): continue
    r = load(fn); summ = json.load(open(fn.replace("_timeseries.npz", "_summary.json")))
    dphi = np.abs(r[:, 1] - summ["phi_b_static"]); t = r[:, 0]
    d = {"dz": summ["dz"], "max_M_int": summ["max_momentum_constraint"], "max_M_shell": summ.get("max_M_shell"), "max_H_shell": summ.get("max_H_shell")}
    for (t1, t2) in [(0.5, 1.5), (1.0, 2.0), (1.5, 2.5), (2.0, 3.0)]:
        m = (t > t1) & (t < t2) & (dphi > 0)
        if m.sum() > 10: d["rate_%g_%g" % (t1, t2)] = float(np.polyfit(t[m], np.log(dphi[m]), 1)[0])
    d["dphi_b(t=0,1,2,3)"] = [float(at(r, T, 1) - summ["phi_b_static"]) for T in (0, 1, 2, 3)]
    d["HJ/H0(t=1,2,3)"] = [float(at(r, T, 2)) for T in (1, 2, 3)]
    sc[tag] = d
out["t1e-3_resolution_scan(linear_theory_rate=1.657)"] = sc
json.dump(out, open(os.path.join(HERE, "REFEREE_ANALYSIS.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
