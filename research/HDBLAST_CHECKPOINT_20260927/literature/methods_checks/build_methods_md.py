#!/usr/bin/env python3
"""Validate the sweep-3 (methods) records and render ../methods.md and methods_sources.json.

The build FAILS if any check fails:
  - required fields present; implication / access / provenance labels valid; ids and URLs unique;
  - provenance "sweep1": the sweep-1 record exists and its URL is byte-identical (so the URL was obtained by a
    WebSearch in this session); access must be "search-snippet-only";
  - provenance "programme": the named programme file exists (read-only) and contains the URL verbatim;
    access must be "programme-record-only";
  - LEADS carry no URL (no 'url' key and no 'http' substring anywhere);
  - at least 15 distinct planned queries; every planned-query index cited by a lead exists;
  - the pilot JSON exists, was produced by the current pilot script (SHA-256), and reports exact_ok and controls_ok.

Usage: python3 build_methods_md.py [--sources PATH] [--out DIR]   (defaults: this folder's records, ../)
"""
import argparse
import hashlib
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DBLAST = "/home/user/unified-theory-maldonado/new-files/D-Blast 3"
SWEEP1 = os.path.join(ROOT, "theory_checks", "theory_sources.json")

ap = argparse.ArgumentParser()
ap.add_argument("--sources", default=os.path.join(HERE, "methods_sources.py"))
ap.add_argument("--out", default=ROOT)
args = ap.parse_args()

spec = importlib.util.spec_from_file_location("ms", args.sources)
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)

ALLOWED_IMP = {"supports", "constrains", "method", "context", "challenges"}
ALLOWED_ACC = {"sweep1": "search-snippet-only", "programme": "programme-record-only"}
REQ = ("id", "theme", "provenance", "access", "title", "authors", "year", "venue", "url", "finding",
       "relevance", "implication", "problems")
errors = []
s1 = {s["id"]: s for s in json.load(open(SWEEP1))["sources"]}
ids, urls = set(), set()
file_cache = {}
for s in ms.SOURCES:
    for k in REQ:
        if not s.get(k):
            errors.append(f"{s.get('id')}: missing {k}")
    if s["id"] in ids:
        errors.append(f"duplicate id {s['id']}")
    ids.add(s["id"])
    if s["url"] in urls:
        errors.append(f"duplicate url {s['url']}")
    urls.add(s["url"])
    if s["implication"] not in ALLOWED_IMP:
        errors.append(f"{s['id']}: bad implication {s['implication']}")
    if s["provenance"] not in ALLOWED_ACC:
        errors.append(f"{s['id']}: bad provenance {s['provenance']}")
        continue
    if s["access"] != ALLOWED_ACC[s["provenance"]]:
        errors.append(f"{s['id']}: access {s['access']} inconsistent with provenance {s['provenance']}")
    if not s["url"].startswith("https://"):
        errors.append(f"{s['id']}: url not https")
    for p in s["problems"]:
        if p not in ms.PROBLEMS:
            errors.append(f"{s['id']}: unknown problem {p}")
    if s["provenance"] == "sweep1":
        rec = s1.get(s.get("sweep1_id"))
        if rec is None:
            errors.append(f"{s['id']}: sweep1 id {s.get('sweep1_id')} not in theory_sources.json")
        elif rec["url"] != s["url"]:
            errors.append(f"{s['id']}: url differs from sweep-1 record {rec['id']} ({rec['url']})")
    else:
        path = os.path.join(DBLAST, s.get("programme_file", ""))
        if not os.path.isfile(path):
            errors.append(f"{s['id']}: programme file missing: {s.get('programme_file')}")
        else:
            if path not in file_cache:
                file_cache[path] = open(path, encoding="utf-8", errors="replace").read()
            if s["url"] not in file_cache[path]:
                errors.append(f"{s['id']}: url not found verbatim in {s['programme_file']}")
for L in ms.LEADS:
    if "url" in L or any("http" in str(v) for v in L.values()):
        errors.append(f"lead carries a URL: {L.get('pointer', '')[:60]}")
    if not (1 <= L["planned_query"] <= len(ms.QUERIES_PLANNED)):
        errors.append(f"lead cites unknown planned query {L['planned_query']}")
if len(set(ms.QUERIES_PLANNED)) < 15:
    errors.append("fewer than 15 distinct planned queries")
pj = os.path.join(HERE, "certification_pilot.json")
if not os.path.exists(pj):
    errors.append("certification_pilot.json missing: run certification_pilot.py first")
    pilot = None
else:
    pilot = json.load(open(pj))
    cur = hashlib.sha256(open(os.path.join(HERE, "certification_pilot.py"), "rb").read()).hexdigest()
    if pilot.get("script_sha256") != cur:
        errors.append("certification_pilot.json is stale (script changed): rerun certification_pilot.py")
    if not (pilot.get("exact_ok") and pilot.get("controls_ok")):
        errors.append("pilot exact checks or controls did not behave as required")
if errors:
    print("BUILD FAILED")
    print("\n".join(errors))
    sys.exit(1)

# ---------------------------------------------------------------- render
THEMES = {"A": "A. Constraint-damped formulations, outer boundaries and error control",
          "B": "B. Thin shells, junction conditions, moving boundaries and collapse diagnostics",
          "C": "C. Gauge-invariant perturbations and negative modes of de Sitter-sliced walls",
          "D": "D. Non-perturbative particle production and backreaction",
          "E": "E. Validated numerics and computer-assisted proofs"}
ACC_TXT = {"search-snippet-only": "search snippet only (second-hand, through the sweep-1 record)",
           "programme-record-only": "programme record only (not opened, searched or re-verified in this sweep)"}
bib = []
for th in "ABCDE":
    bib.append(f"### {THEMES[th]}\n")
    for s in [x for x in ms.SOURCES if x["theme"] == th]:
        if s["provenance"] == "sweep1":
            origin = f"sweep-1 record `{s['sweep1_id']}` (`theory_checks/theory_sources.json`)"
        else:
            origin = f"programme file `D-Blast 3/{s['programme_file']}` ({s['programme_access']})"
        cav = f"\n- *Caveat:* {s['caveat']}" if s.get("caveat") else ""
        bib.append(
            f"**{s['id']}. {s['title']}** ({s['year']})\n"
            f"- Authors: {s['authors']}. {s['venue']}.\n"
            f"- URL: <{s['url']}>\n"
            f"- Access: {ACC_TXT[s['access']]}. Origin of the URL: {origin}.\n"
            f"- Finding: {s['finding']}\n"
            f"- Relevance to HDBLAST ({', '.join(s['problems'])}): {s['relevance'][0].upper() + s['relevance'][1:]}\n"
            f"- Implication: `{s['implication']}`{cav}\n")
leads = []
for th in "ABCDE":
    rows = [L for L in ms.LEADS if L["theme"] == th]
    if not rows:
        continue
    leads.append(f"*{THEMES[th]}*\n")
    for L in rows:
        leads.append(f"- {L['pointer']}. Why: {L['why']} (planned query {L['planned_query']})")
    leads.append("")
pq = "\n".join(f"{i}. `{q}`" for i, q in enumerate(ms.QUERIES_PLANNED, 1))
att = "\n".join(f"- `{a['query']}`: {a['outcome']}" for a in ms.QUERIES_ATTEMPTED)
att += "\n" + "\n".join(f"- fetch `{a['url']}`: {a['outcome']}" for a in ms.FETCH_ATTEMPTED)
pm = "\n".join(f"| {k} | {v} | {', '.join(s['id'] for s in ms.SOURCES if k in s['problems'])} |"
               for k, v in ms.PROBLEMS.items())
ptab = ["| delta | G(u_b) at checkpoint u_b | (eta_b/eta_h)/G - 1 | eta_b: LO/ckpt - 1 | eta_h: LO/ckpt - 1 | H^2: LO/ckpt - 1 | LO vs ckpt `metric_only_H2` |",
        "|---:|---:|---:|---:|---:|---:|---:|"]
for r in pilot["rows"]:
    ptab.append(f"| {r['delta']:g} | {r['G_at_checkpoint_u_b']:.6e} | {r['reldev_ratio_vs_G']:.3e} | "
                f"{r['reldev_eta_b_LO_vs_checkpoint']:.3e} | {r['reldev_eta_h_LO_vs_checkpoint']:.3e} | "
                f"{r['reldev_H2_LO_vs_checkpoint_H2']:.3e} | {r['reldev_H2_LO_vs_checkpoint_metric_only_H2']:.1e} |")
sl = pilot["scaling"]["loglog_slope_three_smallest_delta"]
co = pilot["scaling"]["reldev_over_delta_by_increasing_delta"]
pl = pilot["planning"]
n_s1 = sum(1 for s in ms.SOURCES if s["provenance"] == "sweep1")
n_pr = sum(1 for s in ms.SOURCES if s["provenance"] == "programme")
imp = {k: sum(1 for s in ms.SOURCES if s["implication"] == k) for k in sorted(ALLOWED_IMP)}
yrs = [int(str(s["year"])[:4]) for s in ms.SOURCES]
counts = (f"{len(ms.SOURCES)} sources with URLs ({n_s1} from sweep-1 search records, {n_pr} from programme records); "
          f"{len(ms.LEADS)} unverified leads without URLs; {len(ms.QUERIES_PLANNED)} planned queries (none run); "
          f"implications: " + ", ".join(f"{k} {v}" for k, v in imp.items() if v) +
          f"; published 2020-2026: {sum(1 for y in yrs if y >= 2020)}, before 2020: {sum(1 for y in yrs if y < 2020)}")
repl = {
    "{{COUNTS}}": counts, "{{ATTEMPTED}}": att, "{{PLANNED}}": pq, "{{PROBLEM_MAP}}": pm,
    "{{BIBLIOGRAPHY}}": "\n".join(bib), "{{LEADS}}": "\n".join(leads), "{{PILOT_TABLE}}": "\n".join(ptab),
    "{{SLOPE_RATIO}}": f"{sl['reldev_ratio_vs_G']:.4f}", "{{SLOPE_ETAB}}": f"{sl['reldev_eta_b_LO_vs_checkpoint']:.4f}",
    "{{SLOPE_ETAH}}": f"{sl['reldev_eta_h_LO_vs_checkpoint']:.4f}", "{{SLOPE_H2}}": f"{sl['reldev_H2_LO_vs_checkpoint_H2']:.4f}",
    "{{COEF_RATIO}}": f"{co['reldev_ratio_vs_G'][0]:.4f}", "{{COEF_ETAB}}": f"{co['reldev_eta_b_LO_vs_checkpoint'][0]:.4f}",
    "{{COEF_H2}}": f"{co['reldev_H2_LO_vs_checkpoint_H2'][0]:.5f}",
    "{{COEF_ETAH_LIST}}": ", ".join(f"{v:+.4f}" for v in co["reldev_eta_h_LO_vs_checkpoint"]),
    "{{AMP}}": f"{pl['amplification_G_at_registered_u_b']:.6e}", "{{LOGAMP}}": f"{pl['log10_amplification']:.3f}",
    "{{SUPP}}": f"{pl['unwanted_mode_suppression_log10_exp32_times_ub_minus_u0']:.1f}",
    "{{DIGITS}}": str(pl["digits_for_16_significant_digits_of_eta_h_in_direct_phi_form"]),
    "{{UB}}": f"{pl['registered_u_b']:.6f}", "{{NCONT}}": str(pilot["expectations"]["E2_n_containment_checks"]),
    "{{WRONGSIGN}}": f"{pilot['control_design_note']['wrong_sign_junction_eta_b']:.3e}",
    "{{PILOT_SHA}}": pilot["script_sha256"][:12], "{{RUNTIME}}": str(pilot["runtime_s"]),
    "{{H2EXP}}": f"{pilot['H2_defect_coefficient_expected']:.5f}",
    "{{H2OBS}}": f"{pilot['H2_defect_coefficient_observed_smallest_delta']:.5f}",
    "{{AHAT}}": f"{pl['a_hat_eta_h_times_G_over_delta_registered']:.5f}", "{{M9C64}}": f"{pl['minus_9c_over_64']:.5f}",
}
tpl = open(os.path.join(HERE, "methods_template.md")).read()
for k, v in repl.items():
    tpl = tpl.replace(k, v)
if "{{" in tpl:
    print("BUILD FAILED: unreplaced placeholder"); sys.exit(1)
os.makedirs(args.out, exist_ok=True)
open(os.path.join(args.out, "methods.md"), "w").write(tpl)
dump = {"access_note": ("No WebSearch was executed in sweep 3 (session budget exhausted); sources come from sweep-1 search "
                        "records (snippet-level) or programme records (not re-verified)."),
        "queries_attempted": ms.QUERIES_ATTEMPTED, "fetch_attempted": ms.FETCH_ATTEMPTED,
        "queries_planned_not_run": ms.QUERIES_PLANNED, "problems": ms.PROBLEMS,
        "sources": ms.SOURCES, "leads_unverified_no_url": ms.LEADS, "counts": counts}
jpath = os.path.join(HERE if args.out == ROOT else args.out, "methods_sources.json")
json.dump(dump, open(jpath, "w"), indent=1)
print("BUILD OK:", counts)
