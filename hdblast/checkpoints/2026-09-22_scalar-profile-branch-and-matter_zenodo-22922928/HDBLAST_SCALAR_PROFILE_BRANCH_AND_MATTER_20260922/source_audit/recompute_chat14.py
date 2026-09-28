#!/usr/bin/env python3
"""Read archived Chat14 data; never import or execute its evolution code.

NumPy only. Rerun: python recompute_chat14.py
All event locations use linear interpolation between adjacent recorded samples.
On-shell R5 uses the Einstein equation plus scalar junction condition and is NOT
an independent curvature/constraint test. No long evolution is rerun here.
"""
from pathlib import Path
import hashlib
import io
import json
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "inputs/CHAT14_COMPLETED_REFERENCE.zip"
PREFIX = "HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/"
C = 2 / 1.0357712571566784 - 4 / 3


def sha(b):
    return hashlib.sha256(b).hexdigest()


def crossing(rec, col, target):
    x = rec[:, col] - target
    hits = np.flatnonzero((x[:-1] * x[1:] <= 0) & (x[:-1] != x[1:]))
    if not len(hits):
        return None
    i = int(hits[0])
    w = -x[i] / (x[i+1] - x[i])
    # Constraint records include NaNs; don't interpolate these sparse monitors.
    row = (1-w) * rec[i] + w * rec[i+1]
    return {"coordinate_t":float(row[0]), "H0_tau":float(row[9]),
            "phi_b":float(row[1]), "H_over_H0":float(row[2]),
            "phi_proper_velocity_over_H0":float(row[5]), "b_b":float(row[3])}


def at_coordinate(rec, t):
    return crossing(rec, 0, t)


def at_proper(rec, tau):
    return crossing(rec, 9, tau)


def W(p):
    return 1-p+p**3/3


def U(p):
    return .5*(p*p-1)**2 - 2*W(p)**2/3


def json_finite(value):
    """Retain failed-run records using null for NaN/Inf, not invalid JSON."""
    if isinstance(value, dict):
        return {k:json_finite(v) for k,v in value.items()}
    if isinstance(value, list):
        return [json_finite(v) for v in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main():
    with zipfile.ZipFile(ARCHIVE) as z:
        manifest = json.loads(z.read(PREFIX+"SHA256_MANIFEST.json"))
        entries = manifest["sha256"]
        missing = [n for n in entries if PREFIX+n not in z.namelist()]
        mismatch = [n for n,h in entries.items() if PREFIX+n in z.namelist()
                    and sha(z.read(PREFIX+n)) != h]
        checks = {
            "zip_crc_ok": z.testzip() is None,
            "manifest_hashes_match": not missing and not mismatch,
            "manifest_file_count_match": len(entries) == manifest["n_files"],
            "manifest_total_bytes_match": sum(len(z.read(PREFIX+n)) for n in entries)
                                          == manifest["total_bytes"],
        }
        out = {"archive_sha256":sha(ARCHIVE.read_bytes()),
               "manifest_payload_count":len(entries), "missing":missing,
               "hash_mismatch":mismatch, "runs":{}, "checks":checks}
        recs = {}
        names = sorted(n for n in entries if n.endswith("_summary.json"))
        for name in names:
            summary = json.loads(z.read(PREFIX+name))
            stem = name.removesuffix("_summary.json")
            data = np.load(io.BytesIO(z.read(PREFIX+stem+"_timeseries.npz")))
            rec = data["rec"]
            recs[stem] = rec
            rb = summary["rho_b"]
            end = rec[-1]
            th = summary["L"] * (1-summary["taper"])
            proper_points = {str(t):at_proper(rec,t) for t in [6.,6.5,6.7,6.8,6.9,7.,7.05]}
            p = end[1]
            scalar_norm_over_H0_sq = -end[5]**2 + rb**2*(2*(p*p-1)+summary["t_det"]*C)**2/4
            # Metric junction benchmark ONLY: constant phi=1 fails sigma'(1) BC.
            sig_plus = 2*W(1) + summary["t_det"]*(1+C)
            rs_H = rb*np.sqrt((sig_plus/6)**2 - (1/9)**2)
            run = {"summary":summary,"saved_arrays":data.files,
                   "all_recorded_dynamic_columns_finite":bool(np.all(np.isfinite(rec[:,[0,1,2,3,4,5,7,8,9,10,11]]))),
                   "proper_time_strictly_increasing":bool(np.all(np.diff(rec[:,9])>0)),
                   "H_identity_max_abs_error":float(np.max(np.abs(rec[:,2]-(1+rec[:,4])*np.exp(-rec[:,3])))),
                   "events":{"phi_minus_one":crossing(rec,1,-1),
                             "subbalanced_tension":crossing(rec,1,-1/C),
                             "turnaround":crossing(rec,2,0),
                             "H_minus_5":crossing(rec,2,-5),
                             "H_minus_10":crossing(rec,2,-10),
                             "H_minus_30":crossing(rec,2,-30)},
                   "at_proper_time":proper_points,
                   "taper_earliest_continuum_arrival_t":th,
                   "at_taper_earliest_arrival":at_coordinate(rec,th),
                   "end_after_taper_arrival":bool(end[0]>th),
                   "end_after_far_boundary_arrival":bool(end[0]>-data['z'][0]),
                   "at_end":{"coordinate_t":float(end[0]),"H0_tau":float(end[9]),
                             "phi_b":float(end[1]),"H_over_H0":float(end[2]),
                             "lapse_ratio":float(np.exp(end[3])),
                             "U_phi":float(U(p)),
                             "R5_over_H0_sq_ON_SHELL_INFERENCE":float(scalar_norm_over_H0_sq+rb**2*10*U(p)/3)},
                   "metric_RS_benchmark_NOT_full_scalar_solution_H_over_H0":float(rs_H),
                   "constant_phi_one_scalar_junction_residual":float(summary["t_det"]*C/2),
                   "early_growth_fit":{}}
            turn = run["events"]["turnaround"]
            if turn:
                turn["U_phi"] = float(U(turn["phi_b"]))
                turn["detuning_term"] = float(summary["t_det"]*(1+C*turn["phi_b"]))
            for lo,hi in [(1,2),(2,3),(2.5,3.5),(3,4),(4,5)]:
                dif = np.abs(rec[:,1]-summary["phi_b_static"])
                mask = (rec[:,0]>=lo)&(rec[:,0]<=hi)&(dif>0)&(dif<1e-3)
                if mask.sum()>10:
                    run["early_growth_fit"][f"{lo}:{hi}"] = float(np.polyfit(rec[mask,0],np.log(dif[mask]),1)[0])
            out["runs"][stem] = run
        main_runs = [n for n in out["runs"] if "/v2_t1e3_" in n and out["runs"][n]["summary"]["closure"]==4]
        checks["registered_main_runs_H_identity"] = max(out["runs"][n]["H_identity_max_abs_error"] for n in main_runs)<1e-10
        # Three available registered negative-branch runs have independent event interpolation.
        negative = ["runs/v2_t1e3_minus_c4", "runs/v2_t1e3_minus_c4_fine", "runs/v2_t1e3_minus_c4_finecoarse"]
        vals = [out["runs"][n]["events"]["turnaround"]["H0_tau"] for n in negative]
        checks["three_negative_turnaround_times_within_0p001"] = max(vals)-min(vals)<.001
        vals = [out["runs"][n]["events"]["turnaround"]["phi_b"] for n in negative]
        checks["three_negative_turnaround_fields_within_0p002"] = max(vals)-min(vals)<.002
        # Compare H values at the same proper time, not merely time of threshold crossing.
        a,b = recs[negative[0]],recs[negative[1]]
        out["negative_wall_grid_same_proper_time_comparison"] = {}
        for h in [-5,-10,-22]:
            event = crossing(a,2,h)
            other = at_proper(b,event["H0_tau"])
            out["negative_wall_grid_same_proper_time_comparison"][str(h)] = {
                "coarse_event":event,"fine_at_same_time":other,
                "relative_H_difference":abs(other["H_over_H0"]-h)/abs(h) if other else None}
        out["passed"] = all(checks.values())
        (HERE/"CHAT14_RECOMPUTED_RESULTS.json").write_text(json.dumps(json_finite(out),indent=2,allow_nan=False)+"\n")
        print(json.dumps({"checks":checks,"passed":out["passed"],"run_count":len(out["runs"]),"payload_count":len(entries)},indent=2))
        assert out["passed"]


if __name__ == "__main__":
    main()
