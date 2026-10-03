#!/usr/bin/env python3
"""Read-only saved-trajectory audit. No producer imports or simulations.
Requires numpy. Newly defined diagnostic; does not replace the historical rule.
Report files are new outputs only. Original inputs are read without pickle.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

DEFAULT_TAGS = ["main_Y0.5_dc1e-2_dzf5e-4","main_Y0.7_dc1e-2_dzf5e-4","main_Y1.5_dc1e-2_dzf5e-4","fine_Y2_dc1e-4_dzf2.5e-4","fine_Y2_dc1e-2_dzf2.5e-4","fine_Y3_dc1e-2_dzf1.25e-4","fine_Y5_dc1e-2_dzf1.25e-4_cfl0.25","main_Y2_dc1e-2_dzf5e-4","main_Y2_dc1e-4_dzf5e-4"]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def clean(x):
    if isinstance(x, dict): return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)): return [clean(v) for v in x]
    if isinstance(x, np.ndarray): return clean(x.tolist())
    if isinstance(x, (np.floating, float)): return float(x) if np.isfinite(x) else None
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (np.bool_,)): return bool(x)
    return x

def load_series(path, index, ancestry=()):
    path = Path(path)
    if path in ancestry: raise ValueError("Restart ancestry cycle: " + str(path))
    summary = json.loads(path.read_text())
    npzpath = path.with_name(path.name.replace("_summary.json", "_timeseries.npz"))
    with np.load(npzpath, allow_pickle=False) as z:
        columns = [str(c) for c in z["cols"]]
        rec = z["rec"]
        if rec.ndim != 2 or rec.shape[1] != len(columns): raise ValueError("Invalid cols/rec schema")
        t = {c: rec[:, j].copy() for j, c in enumerate(columns)}
    sources = [dict(summary=str(path), summary_sha256=sha(path), npz=str(npzpath), npz_sha256=sha(npzpath))]
    restart = summary.get("restart_info")
    if restart:
        parent_name = Path(restart["file"]).name.split("_state_T")[0] + "_summary.json"
        candidates = index.get(parent_name, [])
        if len(candidates) != 1:
            raise ValueError("Restart parent not uniquely available: " + parent_name)
        _, parent, parent_sources = load_series(candidates[0], index, ancestry + (path,))
        before = parent["T"] < float(restart["T"]) - 1e-9
        t = {k: np.concatenate([parent[k][before], v]) for k, v in t.items() if k in parent}
        sources = parent_sources + sources
    if len(t["T"]) < 2: raise ValueError("Too few saved records")
    return summary, t, sources

def stats(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    return dict(n=int(len(x)), max=float(x.max()) if len(x) else None,
                median=float(np.median(x)) if len(x) else None,
                p95=float(np.percentile(x, 95)) if len(x) else None)

def longest(t, gate):
    best = None
    start = None
    def candidate(a, b):
        return dict(start_index=int(a), end_index=int(b), n_saved=int(b-a+1),
                    Delta_N=float(t["ln_a"][b]-t["ln_a"][a]),
                    start_H0tau=float(t["H0tau"][a]), end_H0tau=float(t["H0tau"][b]),
                    start_T=float(t["T"][a]), end_T=float(t["T"][b]))
    for i in range(len(gate)):
        continuous = i == 0 or (t["ln_a"][i] > t["ln_a"][i-1] and t["H0tau"][i] > t["H0tau"][i-1])
        if gate[i]:
            if start is None: start = i
            elif not continuous:
                item = candidate(start, i-1)
                if best is None or item["Delta_N"] > best["Delta_N"]: best = item
                start = i
        elif start is not None:
            item = candidate(start, i-1)
            if best is None or item["Delta_N"] > best["Delta_N"]: best = item
            start = None
    if start is not None:
        item = candidate(start, len(gate)-1)
        if best is None or item["Delta_N"] > best["Delta_N"]: best = item
    if best is None: return dict(Delta_N=0.0, n_saved=0, passes_Delta_N_0p5=False)
    best["passes_Delta_N_0p5"] = best["Delta_N"] >= 0.5
    return best

def audit(path, index, savedpoints):
    s, t, sources = load_series(path, index)
    required = ["T", "H0tau", "ln_a", "H_over_H0", "R", "Wy", "vac", "rad",
                "kin", "fric", "phi_b", "v_over_H0", "Hmax_near", "Mmax_near"]
    missing = [k for k in required if k not in t]
    if missing: raise ValueError("Required archived diagnostics absent: " + ", ".join(missing))
    if np.any(np.diff(t["T"]) < -1e-8) or np.any(np.diff(t["H0tau"]) < -1e-8):
        raise ValueError("Nonmonotone time in stitched history")
    p = s["params"]
    rb = float(s["rho_b"])
    Y, delta, c, d = [float(p[k]) for k in ("Y", "delta", "c", "d")]
    if rb <= 0 or Y < 0: raise ValueError("Invalid rb or one-way friction coefficient")
    phi, vh = t["phi_b"], t["v_over_H0"]
    sigma = 2*(1-phi+phi**3/3) + delta*(1+c*phi+d*phi**2/2)
    sigma_phi = 2*(phi**2-1) + delta*(c+d*phi)
    Q = rb**2*t["R"]**2/36
    lin = rb**2*sigma*t["R"]/18
    kinetic = vh**2/12
    mixed = np.abs(sigma_phi*Y*rb*vh)/24
    source_square = Y**2*vh**2/48
    V, W = t["vac"], t["Wy"]
    den0 = lin + Q + np.abs(V) + np.abs(W)
    den_abs = den0 + kinetic + mixed + source_square
    with np.errstate(divide="ignore", invalid="ignore"):
        F0 = lin/den0
        F_abs = lin/den_abs
        r_linear = np.abs(W)/lin
        Omega_r_legacy = t["rad"]/t["H_over_H0"]**2
        r_total_legacy = np.abs(W)/t["rad"]
    checked = np.isfinite(t["Hmax_near"]) & np.isfinite(t["Mmax_near"])
    bad = np.flatnonzero((t["Hmax_near"] > 0.05) | (t["Mmax_near"] > 0.05))
    ie = int(bad[0])-1 if len(bad) else len(t["T"])-1
    if not checked.any(): raise ValueError("No finite near-shell constraint checkpoints")
    reliable = np.arange(len(t["T"])) <= ie
    finite = np.ones(len(t["T"]), dtype=bool)
    for k in required:
        if k not in ("Hmax_near", "Mmax_near"): finite &= np.isfinite(t[k])
    baseline = reliable & finite & (t["H_over_H0"] > 0) & (lin > 0)
    good_abs = baseline & (F_abs >= 0.9) & (r_linear <= 0.03)
    good0 = baseline & (F0 >= 0.9) & (r_linear <= 0.03)
    legacy = reliable & finite & (t["H_over_H0"] > 0) & (Omega_r_legacy >= 0.9) & (r_total_legacy <= 0.03)
    closure = V+t["rad"]+t["kin"]+t["fric"]+W-t["H_over_H0"]**2
    closure_scale = np.abs(V)+np.abs(t["rad"])+np.abs(t["kin"])+np.abs(t["fric"])+np.abs(W)+t["H_over_H0"]**2+1e-300
    ri = s.get("restart_info")
    restart_public = None if not ri else {"parent_state_basename": Path(ri["file"]).name, "T": ri["T"]}
    out = dict(tag=s.get("tag", Path(path).name), params=p, rho_b=rb, sources=sources,
               n_saved=int(len(t["T"])), stop_reason=s.get("stop_reason"), restart_info=restart_public,
               reliability=dict(last_index=ie, first_bad_constraint_index=int(bad[0]) if len(bad) else None,
                   last_H0tau=float(t["H0tau"][ie]) if ie >= 0 else None,
                   checked_records=int(checked.sum()), inherited_cutoff="first near-shell H or M residual > 0.05",
                   note="Residuals are archived intermittently; unchecked records inherit the historical prefix cutoff."),
               longest_F_abs=longest(t, good_abs), longest_F0_upper_bound=longest(t, good0),
               longest_legacy_point_gate=longest(t, legacy),
               F_abs_expanding_reliable=stats(F_abs[baseline]), F0_expanding_reliable=stats(F0[baseline]),
               algebra=dict(rad_reconstruction_abs=stats(np.abs((t["rad"]-lin-Q)[reliable])),
                   kin_reconstruction_abs=stats(np.abs((t["kin"]-kinetic)[reliable])),
                   fric_reconstruction_abs=stats(np.abs((t["fric"]+(sigma_phi*Y*rb*vh)/24+source_square)[reliable])),
                   closure_relative=stats((np.abs(closure)/closure_scale)[reliable]),
                   near_shell_H=stats(t["Hmax_near"][reliable]), near_shell_M=stats(t["Mmax_near"][reliable])))
    use = np.flatnonzero(reliable & finite)
    if len(use) >= 4:
        tau = t["H0tau"][use]
        keep = np.concatenate(([True], np.diff(tau) > 1e-9))
        use = use[keep]
        tau, lna = t["H0tau"][use], t["ln_a"][use]
        R, vv = t["R"][use], vh[use]
        a4 = np.exp(4*lna)
        source = (Y/rb)*vv**2*a4
        integral = float(np.sum(0.5*(source[1:]+source[:-1])*np.diff(tau)))
        delta_Ra4 = float(R[-1]*a4[-1]-R[0]*a4[0])
        dR = np.gradient(R, tau, edge_order=2)
        dilution = 4*t["H_over_H0"][use]*R
        feeding = (Y/rb)*vv**2
        ledger_relative = np.abs(dR+dilution-feeding)/(np.abs(dR)+np.abs(dilution)+np.abs(feeding)+1e-300)
        out["radiation_source_ledger"] = dict(Delta_Ra4=delta_Ra4, integral_Y_vhat2_a4_over_rb=integral,
            relative_endpoint_residual=abs(delta_Ra4-integral)/(abs(delta_Ra4)+abs(integral)+1e-300),
            endpoint_ratio=delta_Ra4/integral if integral else None,
            differential_relative_interior=stats(ledger_relative[3:-3]),
            note="Trapezoidal integration/finite-difference diagnostic of archived points, not an exact certificate.")
    out["historical_savedpoint_comparison"] = []
    tag = Path(path).name.replace("_summary.json", "")
    historical = savedpoints.get(tag, {})
    for label, value in [("plateau", historical.get("plateau")), ("at_max_Omega_r", historical.get("at_max_Omega_r"))]:
        if not value or value.get("H0tau") is None: continue
        target = float(value["H0tau"])
        i = int(np.argmin(np.abs(t["H0tau"]-target)))
        out["historical_savedpoint_comparison"].append(dict(label=label, target_H0tau=target,
            nearest_H0tau=float(t["H0tau"][i]), index=i, reliable=bool(reliable[i]),
            H_over_H0=float(t["H_over_H0"][i]), F0=float(F0[i]), F_abs=float(F_abs[i]),
            abs_W_over_linear_radiation=float(r_linear[i]), Omega_r_legacy=float(Omega_r_legacy[i]),
            conservative_point_pass=bool(good_abs[i])))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-root", type=Path, default=Path("research/HDBLAST_CHECKPOINT_20260930/B3_Y_scan"))
    ap.add_argument("--source-ref", default="not supplied")
    ap.add_argument("--tag", action="append", help="Repeat to override the nine documented default target tags")
    ap.add_argument("--output", type=Path, default=Path(__file__).resolve().parent.parent/"outputs"/"trajectory_audit.json")
    args = ap.parse_args()
    allpaths = list(args.scan_root.rglob("*_summary.json"))
    index = {}
    for p in allpaths: index.setdefault(p.name, []).append(p)
    tags = args.tag or DEFAULT_TAGS
    selected = sorted([p for p in allpaths if p.parent.name in ("main", "fine")
                       and p.name.replace("_summary.json", "") in tags])
    savedpoints = {}
    archived_results = args.scan_root/"B3_RESULTS.json"
    if archived_results.exists(): savedpoints = json.loads(archived_results.read_text()).get("runs", {})
    report = dict(status="post-hoc read-only saved-trajectory diagnostic", source_ref=args.source_ref,
        input_baseline="8f67197b730d4e1c43554b86f224c29cc72629eb",
        script_sha256=sha(__file__), target_tags=tags,
        definitions=dict(linear_radiation="rb^2 sigma(phi) R / 18",
            Q="rb^2 R^2/36", F0="linear/(linear+Q+abs(V)+abs(W)); scalar-free UPPER BOUND",
            F_abs="linear/(linear+Q+abs(V)+abs(W)+vhat^2/12+abs(sigma_phi*Y*rb*vhat)/24+Y^2*vhat^2/48)",
            gate="F_abs>=0.9, abs(W)/linear<=0.03, H>0, linear>0, finite saved diagnostics, historical reliability prefix",
            interval="longest contiguous passing saved-point block with strictly increasing ln_a and proper time",
            limits="Saved-point diagnostic only; no guarantee between samples, no independent resolution convergence; fine restarts share ancestors."),
        runs={}, errors={})
    found = {p.name.replace("_summary.json", "") for p in selected}
    for tag in tags:
        if tag not in found: report["errors"][tag] = "Requested archived summary not found"
    for p in selected:
        try: report["runs"][str(p)] = audit(p, index, savedpoints)
        except Exception as e: report["errors"][str(p)] = type(e).__name__+": "+str(e)
    result = clean(report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print("HDBLAST_TRAJECTORY_JSON="+json.dumps(result, separators=(",", ":"), allow_nan=False))
    return 1 if report["errors"] else 0

if __name__ == "__main__": raise SystemExit(main())
