#!/usr/bin/env python3
"""Check that the key numbers quoted in the 27 Sept 2026 checkpoint report exist
in the workstreams' machine-readable outputs (or follow from them).

Each entry: (id, description, relative path of JSON, quoted value, relative tolerance).
A quoted value "is backed" if some numeric leaf of that JSON (numbers or numeric
strings) lies within the tolerance.  Derived entries compute a value from named
leaves.  Negative controls use exact key paths and must NOT match (stale or
refuted numbers the verifiers corrected).  Also checks the analytic_structure
manifest incident reported by its verifier.

Writes synthesis/QUOTED_NUMBER_CHECKS.json and exits non-zero on any failure.
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def leaves(o, p=()):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, p + (k,))
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, p + (i,))
    elif isinstance(o, bool) or o is None:
        return
    else:
        try:
            yield p, float(o)
        except (TypeError, ValueError):
            return


_cache = {}


def js(rel):
    if rel not in _cache:
        with open(os.path.join(ROOT, rel)) as f:
            _cache[rel] = json.load(f)
    return _cache[rel]


def get(rel, path):
    o = js(rel)
    for k in path:
        o = o[k]
    return o


def backed(rel, value, rtol):
    for p, v in leaves(js(rel)):
        if abs(v - value) <= abs(value) * rtol:
            return list(p), v
    return None, None


Q = [
    # ---- stability, gauge-invariant (Method A)
    ("A1", "calibration mu^2 original shell", "stability_gauge_invariant/SPECTRUM_RESULTS.json", -7.717871625260, 1e-12),
    ("A2", "calibration growth rate", "stability_gauge_invariant/SPECTRUM_RESULTS.json", 1.657193631259, 1e-12),
    ("A3", "B on +1 branch at delta=0.1", "stability_gauge_invariant/SPECTRUM_RESULTS.json", 3.47856, 1e-5),
    ("A4", "B on +1 branch at delta=0.0003", "stability_gauge_invariant/SPECTRUM_RESULTS.json", 3.55532, 1e-5),
    ("A5", "min normalized mismatch (delta=0.1)", "stability_gauge_invariant/SPECTRUM_RESULTS.json", 0.992741, 1e-6),
    ("A6", "B on original shell", "stability_gauge_invariant/SPECTRUM_RESULTS.json", -2.68e-4, 5e-3),
    ("A7", "flipped-sigma'' control root, delta=0.001", "stability_gauge_invariant/CONTROLS_RESULTS.json", -28555.115, 1e-7),
    ("A8", "delta->0 extrapolated calibration", "stability_gauge_invariant/CONTROLS_RESULTS.json", -7.7197956, 1e-7),
    ("A9", "closed-form 4D calibration value", "stability_gauge_invariant/CONTROLS_RESULTS.json", -7.7197959, 1e-7),
    # ---- stability, time domain (Method B)
    ("B1", "time-domain pencil rate original shell h=5e-5", "stability_time_domain/results/KEY_RESULTS.json", 1.657193368, 1e-9),
    ("B2", "flipped control rate shift-invert", "stability_time_domain/results/KEY_RESULTS.json", 167.48925, 1e-7),
    ("B3", "Richardson amplitude rate +1 branch", "stability_time_domain/results/KEY_RESULTS.json", -1.50012, 1e-5),
    ("B4", "gauge eigenvalue L=10", "stability_time_domain/results/KEY_RESULTS.json", 0.0594476, 1e-5),
    ("B5", "far-boundary mode L=6", "stability_time_domain/results/KEY_RESULTS.json", 0.1040027, 1e-6),
    # ---- analytic structure
    ("S1", "H^2 at delta=1e-3 (30 digits)", "analytic_structure/runs/BVP_SCAN_reg.json", 5.9240147943288795996e-5, 1e-14),
    ("S2", "exact delta^9 coefficient of eta_b (verifier)", "analytic_structure/verify/V1_SERIES_INDEPENDENT.json", 3.249147254269e-6, 1e-11),
    ("S3", "exact delta^9 coefficient of H^2 (verifier)", "analytic_structure/verify/V1_SERIES_INDEPENDENT.json", -2.759694789e-7, 1e-9),
    ("S4", "shell position u_b threshold for decoupled scalar bound state (H/k=1/sinh u_b=16.65)", "analytic_structure/SPECTRUM_LEADING_ORDER.json", 0.060034, 1e-4),
    # ---- preheating
    ("P1", "max R at M5 cutoff, screen-passing (data-end)", "preheating/SUMMARY.json", 1.678e-4, 1e-3),
    ("P2", "max R at M5 cutoff over whole scan", "preheating/SUMMARY.json", 2.769e-3, 1e-3),
    ("P3", "max R over post-window profile (verifier)", "preheating/verify/CHECK_CLAIMS.json", 5.91e-4, 2e-3),
    ("P4", "max decay y=1 ratio, screen-passing (corrected)", "preheating/SUMMARY.json", 2.34e-4, 1e-3),
    ("P5", "min lambda_c for R=1 (full)", "preheating/SUMMARY.json", 18.13, 1e-3),
    ("P6", "min lambda_c for R=1 (post)", "preheating/SUMMARY.json", 6.275, 1e-3),
    ("P7", "backreaction range if R=1 (low end)", "preheating/SUMMARY.json", 0.1005, 1e-3),
    ("P8", "backreaction range if R=1 (high end)", "preheating/SUMMARY.json", 0.7917, 1e-3),
    ("P9", "residual-vacuum threshold kappa5^2 rho_crit/H0", "preheating/SUMMARY.json", 0.12563, 1e-4),
    # ---- evolution constraints
    ("E1", "baseline H_max(t=0.5) finest grid", "evolution_constraints/ANALYSIS.json", 6.38e-3, 5e-3),
    ("E2", "remedy A H_max(t=0.5) finest grid", "evolution_constraints/ANALYSIS.json", 9.35e-4, 5e-3),
    ("E3", "remedy A+B H_max(t=0.5) finest grid", "evolution_constraints/ANALYSIS.json", 2.28e-5, 5e-3),
    ("E4", "baseline finest-pair order t=0.5", "evolution_constraints/ANALYSIS.json", 0.24, 3e-2),
    ("E5", "remedy C finest-pair order t=0.5 (max)", "evolution_constraints/ANALYSIS.json", 3.85, 3e-3),
    # ---- mechanisms
    ("M1", "Weyl peak (finest grid)", "mechanisms/M1_WEYL_ALONG_CHAT14.json", 0.0806, 1e-3),
    ("M2", "Weyl peak (other grids)", "mechanisms/M1_WEYL_ALONG_CHAT14.json", 0.0809, 1e-3),
    ("M3", "H_vac^2/H0^2 at delta=1e-3", "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", 0.3681, 1e-3),
    ("M4", "N_max (4D budget, zero-mode Planck masses)", "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", 0.5917, 1e-3),
    ("M5", "e-folds of radiation domination needed", "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", 21.84, 1e-3),
    ("M6", "blast e-folding time if H_vac = observed dark energy (Gyr)", "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", 6.417, 1e-3),
    ("M7", "d* at delta=0.1", "mechanisms/M8_QUADRATIC_TENSION_TUNING.json", -3.1069, 1e-4),
    ("M8", "d* pilot radiation share of H^2", "mechanisms/verify/V2_PILOT_CONSISTENCY.json", 0.690, 2e-3),
    ("M9", "c* signed family: phi_b at c*", "mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json", -1.7287, 1e-4),
    ("M10", "c* pilot: phi_b at window end", "mechanisms/verify/V5_CSTAR_PILOT_SIGNED.json", 0.8002, 1e-3),
    ("M11", "Y=0 pilot reproduces Chat 14 archive |dphi|", "mechanisms/M6_PILOT_COUPLED_5D.json", 7.6e-12, 2e-2),
    # ---- literature / observations
    ("L1", "knee frequency / (1/yr)", "literature/observations_checks/knee_and_scale_checks.json", 0.9972, 2e-4),
    ("L2", "Delta N_eff of knee spectrum", "literature/observations_checks/knee_and_scale_checks.json", 7.09e-4, 2e-3),
    ("L3", "H*ell in registered model", "literature/observations_checks/knee_and_scale_checks.json", 0.069271, 1e-5),
    ("L4", "full vs thin-brane H shift (ppm)", "literature/theory_checks/rs_thin_brane_crosscheck.json", -7.86893, 1e-5),
    ("L5", "germ amplification at delta=1e-3", "literature/methods_checks/certification_pilot.json", 6.288e18, 1e-3),
    # ---- synthesis
    ("Y1", "B(delta->0) fit from Method A numerics", "synthesis/RECONCILIATION.json", 3.5555554, 1e-7),
    ("Y2", "slowest linear decay rate on +1 branch in H0 units", "synthesis/RECONCILIATION.json", -0.91008, 1e-5),
]

DERIVED = [
    ("D1", "Weyl/radiation at extended d* window end (verifier)", 1.23, 5e-3,
     lambda: [r for r in js("mechanisms/verify/V2_PILOT_CONSISTENCY.json")["runs"]
              if r["tag"] == "ext_dstar_Y1_dzf2e-3_phimax1.3"][0]["budget_at_window_end_H0sq_units"]["fractions_of_H2"],
     lambda f: f["weyl"] / f["radiation_terms"], "mechanisms/verify/V2_PILOT_CONSISTENCY.json"),
    ("D2", "Weyl/radiation at original d* stop (producer)", 1.756, 1e-3,
     lambda: [r for r in js("mechanisms/verify/V2_PILOT_CONSISTENCY.json")["runs"]
              if r["tag"] == "pilot_dstar_Y1_dc1e-2_L17_dzf2e-3"][0]["budget_at_window_end_H0sq_units"]["fractions_of_H2"],
     lambda f: f["weyl"] / f["radiation_terms"], "mechanisms/verify/V2_PILOT_CONSISTENCY.json"),
    ("D3", "flipped control: A vs B relative difference", 3.7e-10, 5e-2,
     lambda: [c for c in js("synthesis/RECONCILIATION.json")["checks"] if c["name"].startswith("flipped-sigma'' control growth")][0],
     lambda c: c["rel_diff"], "synthesis/RECONCILIATION.json"),
]

# negative controls: stale/refuted values at exact key paths must NOT equal the quoted old numbers
NEG = [
    ("N1", "stale README decay maximum 1.3e-4 is not the JSON value", "preheating/SUMMARY.json",
     ["max_decay_y1_ratio_at_cutoff_full_screen_passing"], 1.3e-4, 5e-2),
    ("N2", "old 'f_H >= 0.21 along the family' is contradicted: signed family reaches f_H<0.21",
     "mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json", ["rows", -1, "fH_series"], 0.21, 0.5),
    ("N3", "fabricated value not present: calibration mu^2 = -7.7179 (4 digits wrong at 1e-9)",
     "stability_gauge_invariant/SPECTRUM_RESULTS.json", ["calibration", "mu2"], -7.7179, 1e-9),
]

results = []
for qid, desc, rel, val, tol in Q:
    path, v = backed(rel, val, tol)
    results.append({"id": qid, "desc": desc, "file": rel, "quoted": val, "rtol": tol,
                    "found_at": path, "found_value": v, "pass": path is not None})
for qid, desc, val, tol, getter, fn, rel in DERIVED:
    v = fn(getter())
    ok = abs(v - val) <= abs(val) * tol
    results.append({"id": qid, "desc": desc, "file": rel, "quoted": val, "rtol": tol, "derived_value": v, "pass": ok})
for qid, desc, rel, path, val, tol in NEG:
    v = float(get(rel, path))
    match = abs(v - val) <= abs(val) * tol
    results.append({"id": qid, "desc": desc, "file": rel, "path": path, "value_at_path": v, "stale_value": val,
                    "rtol": tol, "control": True, "pass": not match})

# ---- analytic_structure manifest incident (reported by its verifier)
man_rel = "analytic_structure/MANIFEST.sha256.json"
man = js(man_rel)
mism, n_verify, n_orig = [], 0, 0
for k, h in man.items():
    if k.startswith("verify/"):
        n_verify += 1
        continue
    n_orig += 1
    p = os.path.join(ROOT, "analytic_structure", k)
    if not os.path.exists(p):
        mism.append([k, "missing"])
        continue
    with open(p, "rb") as f:
        if hashlib.sha256(f.read()).hexdigest() != h:
            mism.append([k, "hash differs"])
manifest_check = {"file": man_rel, "n_entries": len(man), "n_non_verify_entries": n_orig,
                  "n_verify_entries_added_by_auditor": n_verify, "non_verify_mismatches": mism,
                  "note": "Verifier reported regenerating this manifest by accident; the non-verify entries are "
                          "compared with the current files here. The manifest was left as found (not edited)."}

out = {"status": "consistency check of quoted numbers against saved JSON (no new physics)",
       "n_checks": len(results), "n_pass": sum(r["pass"] for r in results),
       "all_pass": all(r["pass"] for r in results), "results": results,
       "analytic_structure_manifest": manifest_check}
with open(os.path.join(HERE, "QUOTED_NUMBER_CHECKS.json"), "w") as f:
    json.dump(out, f, indent=1)
for r in results:
    print(("PASS " if r["pass"] else "FAIL ") + r["id"] + " " + r["desc"])
print("manifest:", n_orig, "original entries,", n_verify, "verify entries, mismatches:", mism)
print("all_pass:", out["all_pass"])
sys.exit(0 if out["all_pass"] else 1)
