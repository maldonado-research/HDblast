#!/usr/bin/env python3
"""Critic corrections and checks for the HDBLAST checkpoint of 27 September 2026.

A final completeness and overclaim review of the top-level report found statements stronger than the audits
support, missing caveats, one number with no JSON source, an incomplete list of unverified citations, and AI
product names.  The hand-written top-level Markdown files were edited directly.  The two files that are generated
by `synthesis/` builders (RESULTS_SUMMARY.json and CLAIM_REVIEW.md) are corrected here, reproducibly, so that the
builders themselves (outside the top level) stay untouched.

Usage (from this folder):
    python3 critic_corrections.py              # patch RESULTS_SUMMARY.json, re-render CLAIM_REVIEW.md, run checks,
                                               # write CRITIC_CHECKS.json (idempotent)
    python3 critic_corrections.py --check-only # verify the published files and numbers; writes nothing

Run order: synthesis/build_results_summary.py -> synthesis/build_claim_review.py -> this script -> make_manifest.
Reads only files inside this folder.  Exit code 0 only if every check passes.
"""
import glob
import json
import os
import re
import codecs
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
CHECK_ONLY = "--check-only" in sys.argv


def P(rel):
    return os.path.join(ROOT, rel)


def js(rel):
    with open(P(rel)) as f:
        return json.load(f)


def get(obj, path):
    for k in path:
        obj = obj[k]
    return obj


# ---------------------------------------------------------------- 1. corrections to RESULTS_SUMMARY.json
# (target, field, op, old, new): op 'replace' swaps a substring; op 'append' adds text if absent.
CRITIC_TAG = "[critic 28 Sept 2026]"
JSON_PATCHES = [
    ("author", None, "regex", r"prepared with \w+ \(\w+\)", "prepared with AI assistance"),
    ("headline", 0, "replace", "at all six tested detunings;",
     "at all six tested detunings (the second method covers only the homogeneous sector, at three detunings);"),
    ("headline", 1, "replace", "instant preheating delivers <=1.7e-4 of the residual-vacuum threshold at the M5 cutoff (conditional)",
     "instant preheating delivers <=1.7e-4 of the residual-vacuum threshold at the M5 cutoff for screen-passing "
     "couplings (data-end energy; <=5.9e-4 over the stored profile; conditional on the fixed archived trajectory)"),
    ("headline", 2, "replace", "(radiation reaches 69% of H^2 but bulk Weyl radiation stays above N_eff limits so far)",
     "(in delta=0.1 pilots radiation reaches 69% of H^2, inside a budget with large cancellations, but bulk Weyl "
     "radiation stays above N_eff limits so far)"),
    ("headline", 2, "replace", "(found by the auditor).",
     "(found by the auditor at delta=0.1; one exploratory single-grid pilot)."),
    ("GI-4", "note", "append", None,
     "Real mu^2 scanned on a 360-point grid in [-400, 2.2499]; mu^2<-400 and complex mu^2 outside the winding "
     "contours are excluded only by the conditional B>0 theorem (GI-6). Regularity (normalisability) at the cone is "
     "an assumed criterion. Auditor's margin test: an instability needs B in [-0.659, 0.0081], actual B=3.479-3.555. "
     + CRITIC_TAG),
    ("TD-3", "note", "append", None,
     "Homogeneous (FRW-symmetric) sector only, at delta=0.001 with 0.003 and 0.01 as perturbed-parameter checks. "
     "The time-domain signal barely reaches the shell (peak |f_b|~8.5e-9), so this claim rests on the spectral "
     "evidence. " + CRITIC_TAG),
    ("MX-6", "claim", "replace", "delta~1e-62..1e-67", "delta~1.1e-62..7.6e-68"),
    ("MX-8", "note", "append", None,
     "5D pilots run at delta=0.1 (100x the registered detuning). " + CRITIC_TAG),
    ("MX-10", "correction", "replace", "and the c* pilot rolls",
     "and the c* pilot (single grid, one seed sign, Y=0, delta=0.1 only) rolls"),
    ("MX-11", "correction", "replace",
     "and ~0.6 after the chart freeze, still >=6x above the N_eff limit (<=0.1); then H turns negative (gauge artefact or reversal).",
     "(about 12x the Delta N_eff-derived limit of 0.1). Values after the chart freeze are not usable: proper time "
     "advances by only ~0.2 H0 tau while the ratio falls through 0.6 and below zero and then swings back as H "
     "turns negative (gauge artefact or reversal); the earlier '~0.6 after the freeze' figure is withdrawn "
     "(CRITIC_CHECKS.json, post_freeze_dstar)."),
    ("MX-11", "note", "append", None,
     "Pilots at delta=0.1 only; d* convergence rests on one grid pair (dz_fine 1e-3 and 2e-3). " + CRITIC_TAG),
    ("LT-2", "note", "append", None,
     "The O(delta^2) term predicts -7.849 ppm of the measured -7.869 ppm shift in H; higher orders supply the rest. "
     + CRITIC_TAG),
    ("SY-1", "note", "append", None,
     "Slowest-decay agreement is within the supported +/-1e-3 of TD-4. " + CRITIC_TAG),
]

CRITIC_REVIEW_JSON = {
    "date": "2026-09-28",
    "scope": "completeness and overclaim review of the top-level report files; workstream outputs not modified",
    "script": "critic_corrections.py",
    "checks": "CRITIC_CHECKS.json",
    "summary": [
        "wording stronger than the audits support was tightened (stability scope, detection-power controls, "
        "decay-rate precision, thin-brane ppm attribution, audit re-run coverage)",
        "missing caveats added (delta=0.1 pilots, Method A scan range and cone-regularity assumption, Method B weak "
        "shell excitation, preheating fixed trajectory and 'post' cutoff maximum, projection non-convergence)",
        "the unsourced '~0.6 after the chart freeze' Weyl/radiation value was withdrawn",
        "the unverified-citation list was completed; arXiv:2609.21421 does have a search record",
        "contradiction fixed: 'each workstream audited' vs literature/program_map self-checked only",
        "AI product names removed from the top-level files",
    ],
}


def apply_json_patches():
    R = js("RESULTS_SUMMARY.json")
    by_id = {c["id"]: c for c in R["claims"]}
    log = []
    for target, field, op, old, new in JSON_PATCHES:
        if target in by_id:
            holder, key = by_id[target], field
        elif isinstance(field, int):
            holder, key = R[target], field
        else:
            holder, key = R, target
        cur = holder[key]
        if op == "regex":
            if new in cur and not re.search(old, cur):
                log.append((target, field, "already"))
                continue
            if len(re.findall(old, cur)) != 1:
                raise SystemExit(f"regex patch target not found exactly once: {target}")
            holder[key] = re.sub(old, new, cur)
            log.append((target, field, "applied"))
            continue
        if op == "append":
            if new in cur:
                log.append((target, field, "already"))
                continue
            holder[key] = (cur + " " + new).strip() if cur else new
        else:
            if new in cur and old not in cur:
                log.append((target, field, "already"))
                continue
            if cur.count(old) != 1:
                raise SystemExit(f"patch target not found exactly once: {target}.{field}: {old[:60]}")
            holder[key] = cur.replace(old, new)
        log.append((target, field, "applied"))
    R["critic_review"] = CRITIC_REVIEW_JSON
    with open(P("RESULTS_SUMMARY.json"), "w") as f:
        json.dump(R, f, indent=1, ensure_ascii=False)
        f.write("\n")
    return log


def json_patches_present():
    R = js("RESULTS_SUMMARY.json")
    by_id = {c["id"]: c for c in R["claims"]}
    bad = []
    for target, field, op, old, new in JSON_PATCHES:
        if target in by_id:
            cur = by_id[target][field]
        elif isinstance(field, int):
            cur = R[target][field]
        else:
            cur = R[target]
        if new not in cur or (op == "replace" and old in cur) or (op == "regex" and re.search(old, cur)):
            bad.append(f"{target}.{field}")
    if "critic_review" not in R:
        bad.append("critic_review")
    return bad


# ---------------------------------------------------------------- 2. citations cross-check
LIT_FILES = [p for p in glob.glob(P("literature/**/*"), recursive=True) if os.path.isfile(p)
             and not p.endswith((".npz", ".pyc"))] + [P("LITERATURE_2022_2026.md"), P("synthesis/LITERATURE_MERGED.json")]
ARXIV = re.compile(r"((?:hep-th|hep-ph|gr-qc|astro-ph)/\d{7}|(?<![\d.])\d{4}\.\d{4,5}(?![\d]))")
URL = re.compile(r"https?://(?:[^\s()\]>`]|\([^\s()]*\))+")
DOC_GLOBS = ["00_READ_FIRST.md", "CLAIM_REVIEW.md", "PUBLIC_SUMMARY.md", "NEXT_TESTS.md", "README.md", "REPRODUCE.md",
             "*/README*.md", "*/VERIFICATION.md"]


def lit_corpus():
    s = []
    for p in LIT_FILES:
        with open(p, errors="ignore") as f:
            s.append(f.read())
    return "\n".join(s)


def citation_check(corpus):
    out = {}
    for g in DOC_GLOBS:
        for p in sorted(glob.glob(P(g))):
            rel = os.path.relpath(p, ROOT)
            if rel.startswith("literature/"):
                continue
            t = open(p).read()
            if rel == "CLAIM_REVIEW.md":  # the critic list itself names the missing ids
                t = t.split("## Citations")[0]
            ids = sorted(set(ARXIV.findall(t)))
            # only count ids that appear as arXiv links or arXiv-style references, not numbers like 2026.0928
            ids = [i for i in ids if "/" in i or re.search(r"(arxiv[^\s]*|arXiv:)" + re.escape(i), t)]
            miss_ids = [i for i in ids if i not in corpus]
            urls = sorted(set(u.rstrip(".,;:'\"") for u in URL.findall(t)))
            miss_urls = []
            for u in urls:
                core = re.sub(r"^https?://", "", u)
                if core in corpus or any(i in u and i in corpus for i in ARXIV.findall(u)):
                    continue
                doi = re.search(r"10\.\d{4,}/[^\s]+", u)
                if doi and doi.group(0) in corpus:
                    continue
                miss_urls.append(u)
            if miss_ids or miss_urls:
                out[rel] = {"arxiv_ids_without_literature_record": miss_ids, "other_links_without_record": miss_urls}
    return out


# ---------------------------------------------------------------- 3. CLAIM_REVIEW.md rendering and text patches
def citations_section(cit):
    ids_by_ws = {}
    links_by_ws = {}
    for rel, v in cit.items():
        ws = rel.split("/")[0] if "/" in rel else "top level"
        ids_by_ws.setdefault(ws, set()).update(v["arxiv_ids_without_literature_record"])
        links_by_ws.setdefault(ws, set()).update(u for u in v["other_links_without_record"] if "arxiv.org" not in u)
    lines = ["## Citations", "",
             "Paper sites were blocked for downloads and the shared web-search budget ran out during the audits, so no "
             "auditor could re-verify any citation. All citations in this round are search-snippet-only or taken from "
             "earlier programme records (`LITERATURE_2022_2026.md`).", "",
             "A cross-check by `critic_corrections.py` (`CRITIC_CHECKS.json`, key `citations`) compares every arXiv "
             "identifier and link in the top-level files, workstream READMEs and audit reports with the literature files "
             "(`literature/`, `LITERATURE_2022_2026.md`, `synthesis/LITERATURE_MERGED.json`). The following have **no "
             "record** there. Their provenance is only the workstream's own statement, so treat them as unverified:", ""]
    for ws in sorted(ids_by_ws):
        ids = sorted(ids_by_ws[ws])
        links = sorted(links_by_ws.get(ws, ()))
        if not ids and not links:
            continue
        parts = []
        if ids:
            parts.append("arXiv " + ", ".join(ids))
        if links:
            parts.append(f"{len(links)} non-arXiv link(s): " + ", ".join(f"<{u}>" for u in links))
        lines.append(f"- `{ws}`: " + "; ".join(parts) + ".")
    lines += ["",
              "Corrections to the earlier version of this section: arXiv:2609.21421 *does* have a search-snippet record "
              "(theory E11, methods C9). Not every unrecorded citation is context only: the ΔN_eff input 2.990 ± 0.070 "
              "(arXiv:2603.13226) feeds the dark-radiation limit used in MX-2 and MX-11, and the ACT DR6 N_eff value is "
              "quoted differently in two sweeps (2.86 ± 0.13 and 2.89 ± 0.11). The negative conclusions of MX-2 and "
              "MX-11 hold for either value, but the thresholds must be checked against the papers before any external "
              "use. No computed (non-observational) result depends on these citations.", ""]
    return "\n".join(lines)


CRITIC_SECTION = """## Critic review (28 September 2026)

A final completeness and overclaim review of the top-level documents (`00_READ_FIRST.md`, `README.md`,
`PUBLIC_SUMMARY.md`, `NEXT_TESTS.md`, `REPRODUCE.md`, this file and `RESULTS_SUMMARY.json`) against every workstream
README, audit report and JSON output. Workstream outputs were not modified. Checks: `critic_corrections.py` →
`CRITIC_CHECKS.json`.

Corrections made:

1. **d\\* pilot after the chart freeze (MX-11).** The report quoted Weyl/radiation "≈ 0.6 after the chart freeze,
   ≥ 6× the N_eff limit". No JSON contained that value, and the mechanisms workstream itself rules post-freeze values
   out as gauge artefacts. Recomputed from the saved time series: after the freeze proper time advances by only
   ≈ 0.2 H₀τ while the ratio falls through 0.6 and below zero and then swings back as H turns negative, so there is
   no plateau. The figure is withdrawn; the supported statement is 1.23 at the end of the reliable window (≈ 12× the
   ΔN_eff-derived limit of 0.1).
2. **Stability scope.** "Two completely different methods" and "stable at δ = 0.0003–0.1" were stronger than the
   evidence: Method B covers only the homogeneous sector at three δ; Method A scans real μ² on [−400, 2.2499] and
   relies on the conditional B > 0 theorem beyond it and on cone regularity as the normalisability criterion. The
   planted-tachyon controls test only a strong instability; the auditor's B\\* margin is the relevant detection-power
   evidence. Method B's time-domain runs barely excite the shell (≈ 8.5×10⁻⁹), so its conclusion rests on the spectra.
3. **Decay-rate precision.** The reconciliation row "−1.5001, agree within 1×10⁻⁴" contradicted the audit
   (−1.500 ± 0.001, order-sensitive Richardson); corrected.
4. **Pilots at δ = 0.1.** The dissipative-coupling, d\\* and c\\* pilots all ran at δ = 0.1 (100× the registered
   detuning); the c\\* run is a single-grid, single-seed, Y = 0 exploration; d\\* convergence rests on one grid pair.
   These caveats were missing from the summary documents.
5. **Preheating bound.** The fixed-trajectory condition and the "post" cutoff maximum (≤ 4.0×10⁻³) were added next
   to the headline ≤ 1.7×10⁻⁴.
6. **Smaller numerical wording.** δ–Λ lock range is 1.1×10⁻⁶² to 7.6×10⁻⁶⁸ (not "10⁻⁶²–10⁻⁶⁷"); the O(δ²)
   thin-brane term explains −7.849 of the −7.869 ppm shift (not all of it); H_vac²/H₀² is 0.372 at δ = 0.01 and
   0.406 at δ = 0.1; the multiprecision solvers' agreement is 10⁻³⁹–10⁻⁴¹, their working precision 40–65 digits.
7. **Contradictions between documents.** README and 00_READ_FIRST said every workstream had an independent audit,
   while the literature sweeps and program map are self-checked only; REPRODUCE said every audited script was re-run
   bit-identically, while several were re-run only in subsets, some agree to 10⁻⁹–10⁻⁸, and m7 was not re-run.
8. **Reproduction.** Commands for the audit scripts that produce quoted corrections (mechanisms v1/v2/v4 and the
   extended d\\* pilot, preheating `check_claims.py`, evolution `recompute_saved.py`, time-domain
   `v4_compare_and_extract.py`) and the order of the post-build critic step were added to `REPRODUCE.md`.
9. **Citations.** The unverified-citation list was completed (section above); CMB-S4, which has no literature record,
   was removed from `NEXT_TESTS.md`; the N_eff normalisation mismatch (ρ_DR/ρ_SM at production vs ρ_dr/ρ_γ at
   recombination) is now flagged in the pre-registered decision rule of next test A1.
10. **AI product names** were removed from the top-level files.

Remaining gaps the critic could not fix (outside the top level): workstream READMEs keep their pre-audit wording
(for example the stale preheating decay column and the "13–100×" damping range); `synthesis/build_results_summary.py`
and `synthesis/build_claim_review.py` still emit the pre-critic text, which this script re-corrects;
`analytic_structure/MANIFEST.sha256.json` still carries the 25 accidental `verify/` entries; `__pycache__` folders
are present in some workstream folders.
"""


def render_claim_review(cit):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    subprocess.run([sys.executable, P("synthesis/build_claim_review.py")], check=True, env=env,
                   stdout=subprocess.DEVNULL)
    t = open(P("CLAIM_REVIEW.md")).read()
    head = r"Prepared with \w+ \(\w+\) for Ricardo Maldonado's HDBLAST program\."
    assert len(re.findall(head, t)) == 1
    t = re.sub(head, "Prepared with AI assistance for Ricardo Maldonado's HDBLAST program.", t)
    t = t.replace("The table is generated from `RESULTS_SUMMARY.json` by `synthesis/build_claim_review.py`.",
                  "The table is generated from `RESULTS_SUMMARY.json` by `synthesis/build_claim_review.py`; the "
                  "critic corrections (last section) are then applied by `critic_corrections.py`.")
    a, b = t.index("## Citations"), t.index("## Process notes")
    t = t[:a] + citations_section(cit) + "\n" + t[b:]
    t = t.rstrip("\n") + "\n\n" + CRITIC_SECTION
    with open(P("CLAIM_REVIEW.md"), "w") as f:
        f.write(t)


# ---------------------------------------------------------------- 4. numbers added or changed by the critic
NUMBERS = [
    # (id, description, quoted, rtol, file, key path)
    ("P1", "preheating: max R over stored profile, screen-passing", 5.9e-4, 0.01,
     "preheating/verify/CHECK_CLAIMS.json", ["aggregates", "max_R_over_profile_screen", "value"]),
    ("P2", "preheating: max R with 'post' cutoff, screen-passing", 4.0e-3, 0.02,
     "preheating/verify/CHECK_CLAIMS.json", ["aggregates", "max_R_post_screen", "value"]),
    ("P3", "preheating: decay maximum y=1, screen-passing", 2.34e-4, 0.005,
     "preheating/verify/CHECK_CLAIMS.json", ["aggregates", "max_decay_y1_screen", "value"]),
    ("M1", "delta-Lambda lock, ell=38.6 um", 1.1e-62, 0.03,
     "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", ["physical_scales_registered", "delta_required", 3, "delta_required_for_Hvac_eq_HLambda"]),
    ("M2", "delta-Lambda lock, ell=0.1 um", 7.6e-68, 0.01,
     "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", ["physical_scales_registered", "delta_required", 0, "delta_required_for_Hvac_eq_HLambda"]),
    ("M3", "H_vac^2/H0^2 at delta=0.01", 0.372, 0.002,
     "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", ["rows", 3, "budget", "H_vac2_over_H02"]),
    ("M4", "H_vac^2/H0^2 at delta=0.1", 0.406, 0.002,
     "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", ["rows", 5, "budget", "H_vac2_over_H02"]),
    ("M5", "N_max at delta=0.1", 0.47, 0.01,
     "mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json", ["rows", 5, "budget", "N_max_E"]),
    ("M6", "c* (signed family, delta=0.1)", -0.99307, 1e-5, "mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json", ["rows", 17, "c"]),
    ("M7", "c* shell cone value phi_h", -1.00695, 1e-5, "mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json", ["rows", 17, "phi_h"]),
    ("M8", "c* pilot H/H0 at window end", 0.14, 0.01, "mechanisms/verify/V5_CSTAR_PILOT_SIGNED.json", ["window_end", "H_over_H0"]),
    ("L1", "thin-brane O(delta^2) prediction (ppm)", -7.849, 2e-4,
     "literature/theory_checks/rs_thin_brane_crosscheck.json", ["numbers_at_registered_point", "predicted_ppm_from_c2_over_384"]),
    ("L2", "full vs thin-brane shift (ppm)", -7.869, 2e-4,
     "literature/theory_checks/rs_thin_brane_crosscheck.json", ["numbers_at_registered_point", "ppm_shift_in_H_full_vs_thin"]),
    ("T1", "Method B Richardson assuming order 2", -1.4994, 1e-4,
     "stability_time_domain/verify/v4_compare_and_extract.json", ["richardson_assuming_order_2"]),
    ("T2", "Method B Richardson assuming order 4", -1.5007, 1e-4,
     "stability_time_domain/verify/v4_compare_and_extract.json", ["richardson_assuming_order_4"]),
    ("T3", "Method B peak shell response |f_b| (finest grid)", 8.5e-9, 0.01,
     "stability_time_domain/verify/v4_compare_and_extract.json", ["independent_extraction", "plus_lin_h5e-5", "f_b_abs_max"]),
    ("G1", "Method A real mu^2 scan: number of points", 360, 0.0, "stability_gauge_invariant/SPECTRUM_RESULTS.json", None),
]


def number_checks():
    res = []
    for nid, desc, q, rt, f, path in NUMBERS:
        d = js(f)
        v = len(d["grid"]) if nid == "G1" else get(d, path)
        v = float(v) if isinstance(v, str) else v  # some JSON outputs store numbers as strings
        ok = abs(v - q) <= rt * abs(q) if rt > 0 else v == q
        res.append(dict(id=nid, desc=desc, quoted=q, rtol=rt, file=f, value=v, passed=bool(ok)))
    # Method A scan interval and B* margin over the +1 branches
    sp = js("stability_gauge_invariant/SPECTRUM_RESULTS.json")
    lo, hi = min(sp["grid"]), max(sp["grid"])
    res.append(dict(id="G2", desc="Method A scan interval [-400, 2.2499]", quoted=[-400, 2.2499], value=[lo, hi],
                    passed=bool(abs(lo + 400) < 1e-9 and abs(hi - 2.2499) < 1e-12)))
    br = js("stability_gauge_invariant/verify/indep_spectrum.json")["branches"]
    plus = [b["Bstar_control"] for b in br if b["Bstar_control"]["B_actual"] > 1]
    bmin = min(b["Bstar_min_mu2_neg"] for b in plus)
    bmax = max(b["Bstar_max_mu2_neg"] for b in plus)
    bact = [b["B_actual"] for b in plus]
    res.append(dict(id="G3", desc="B* range over +1 branches [-0.659, 0.0081]; actual B 3.479-3.555",
                    quoted=[-0.659, 0.0081, 3.479, 3.555], value=[bmin, bmax, min(bact), max(bact)],
                    passed=bool(abs(bmin + 0.659) < 5e-4 and abs(bmax - 0.0081) < 5e-5
                                and abs(min(bact) - 3.479) < 5e-4 and abs(max(bact) - 3.555) < 5e-4)))
    # damping factors of Remedy B vs baseline (t=0.5 and 1): quoted 'about 25-200'
    tb = js("evolution_constraints/verify/recompute_saved.json")["table_H_max"]
    fac = [b / k for tt in ("0.5", "1.0") for b, k in zip(tb["base"][tt]["values"], tb["k10"][tt]["values"])]
    res.append(dict(id="E1", desc="Remedy B / baseline H_max factors at t=0.5 and 1 lie in 25-200", quoted=[25, 200],
                    value=[min(fac), max(fac)], passed=bool(min(fac) >= 25 and max(fac) <= 200)))
    # d* budget at the original stop
    v2 = {r["tag"]: r for r in js("mechanisms/verify/V2_PILOT_CONSISTENCY.json")["runs"]}
    fr = v2["pilot_dstar_Y1_dc1e-2_L17_dzf2e-3"]["budget_at_window_end_H0sq_units"]["fractions_of_H2"]
    res.append(dict(id="D1", desc="d* stop budget: vacuum -117%, Weyl +121%, radiation +69%",
                    quoted=[-1.17, 1.21, 0.69],
                    value=[fr["vacuum_sigma2_minus_sigma1sq_plus_U"], fr["weyl"], fr["radiation_terms"]],
                    passed=bool(abs(fr["vacuum_sigma2_minus_sigma1sq_plus_U"] + 1.17) < 0.006
                                and abs(fr["weyl"] - 1.21) < 0.006 and abs(fr["radiation_terms"] - 0.69) < 0.006)))
    fe = v2["ext_dstar_Y1_dzf2e-3_phimax1.3"]["budget_at_window_end_H0sq_units"]["fractions_of_H2"]
    r_end = fe["weyl"] / fe["radiation_terms"]
    res.append(dict(id="D2", desc="extended d*: Weyl/radiation at reliable-window end, and multiple of the 0.1 limit",
                    quoted=[1.23, 12], value=[r_end, r_end / 0.1],
                    passed=bool(abs(r_end - 1.23) < 0.006 and 11.5 < r_end / 0.1 < 12.5)))
    return res


# ---------------------------------------------------------------- 5. the post-freeze Weyl/radiation ratio
def post_freeze_dstar():
    """Recompute Weyl/(radiation terms) along the extended d* run, with the same definitions and window rule as
    mechanisms/verify/v2_pilot_consistency.py, and characterise what happens after the chart freeze."""
    f = P("mechanisms/verify/pilot/ext_dstar_Y1_dzf2e-3_phimax1.3_summary.json")
    s = json.load(open(f))
    rec = np.load(f.replace("_summary.json", "_timeseries.npz"))["rec"]
    W = lambda p: 1 - p + p ** 3 / 3
    td, c, dq, rb = s["t_det"], s["c"], s["d_quad"], s["rho_b"]
    sig = lambda p: 2 * W(p) + td * (1 + c * p + dq * p * p / 2)
    t, ss, phi, h, R, Wy = rec[:, 0], rec[:, 8], rec[:, 1], rec[:, 2], rec[:, 6], rec[:, 7]
    lapse = np.gradient(ss, t)
    bad = np.where((lapse < 0.2) & (t > 1.0))[0]
    ie = int(bad[0]) - 1
    rad = rb * rb * (sig(phi) * R / 18 + R * R / 36)
    ratio = Wy / np.where(rad != 0, rad, np.nan)
    ineg = int(np.where((h < 0) & (t > t[ie]))[0][0])
    seg = ratio[ie:ineg]
    # control: the same diagnostic must reproduce the audited window-end value (1.2320)
    ctrl = abs(ratio[ie] - 1.232045878803561) < 1e-9
    out = dict(
        window_end_t=float(t[ie]), window_end_H0tau=float(ss[ie]), ratio_at_window_end=float(ratio[ie]),
        control_reproduces_V2_window_end_ratio=bool(ctrl),
        first_H_negative_t=float(t[ineg]), H0tau_at_first_H_negative=float(ss[ineg]),
        proper_time_advance_window_end_to_H_negative=float(ss[ineg] - ss[ie]),
        coordinate_time_advance=float(t[ineg] - t[ie]),
        ratio_min_before_H_negative=float(np.nanmin(seg)), ratio_max_before_H_negative=float(np.nanmax(seg)),
        ratio_crosses_0p6=bool(np.nanmax(seg) > 0.6 > np.nanmin(seg)),
        ratio_changes_sign=bool(np.nanmin(seg) < 0 < np.nanmax(seg)),
        samples=[[float(t[i]), float(ss[i]), float(h[i]), float(ratio[i])]
                 for i in range(ie, min(len(t), ineg + 60), max(1, (ineg - ie) // 12))],
        samples_columns=["t", "H0tau", "H/H0", "Weyl/radiation_terms"],
    )
    # 'plateau at ~0.6' would need the ratio to stay within +/-10% of 0.6 over the post-freeze positive-H segment
    out["plateau_near_0p6"] = bool(np.nanmax(np.abs(seg / 0.6 - 1)) < 0.1)
    out["passed"] = bool(ctrl and out["ratio_crosses_0p6"] and out["ratio_changes_sign"]
                         and not out["plateau_near_0p6"] and out["proper_time_advance_window_end_to_H_negative"] < 0.25)
    return out


# ---------------------------------------------------------------- 6. wording scans (AI names, withdrawn phrases)
# product names are stored rot13-encoded so that this file does not itself contain them
AI_NAME = re.compile(codecs.decode(r"Pynhqr|Naguebcvp|PungTCG|BcraNV|TCG-?\q|Trzvav|Nfgen \q|Bchf \q|Fbaarg \q|Unvxh \q", "rot13"))
TOP_TEXT = ["00_READ_FIRST.md", "CLAIM_REVIEW.md", "PUBLIC_SUMMARY.md", "NEXT_TESTS.md", "README.md", "REPRODUCE.md",
            "RESULTS_SUMMARY.json", "LITERATURE_2022_2026.md", "requirements.txt", "critic_corrections.py"]
WITHDRAWN = [
    "≈ 0.6 after the chart freeze", "~0.6 after the chart freeze, still", "≥ 6× the N_eff limit", "≈ 0.6× after the chart freeze",
    "agree within 1×10⁻⁴", "δ ≈ 10⁻⁶²–10⁻⁶⁷", "delta~1e-62..1e-67", "one independent audit per workstream",
    "Seven workstreams, each with an independent adversarial audit", "which accounts for the 7.87 ppm correction",
    "CMB-S4", "Two completely different methods", "obtained bit-identical results, apart from runtime fields",
    "Two independent 40–65-digit solvers",
]


def wording_scans():
    ai_hits, stale_hits = [], []
    for rel in TOP_TEXT:
        t = open(P(rel)).read()
        for m in AI_NAME.finditer(t):
            ai_hits.append([rel, m.group(0)])
        if rel == "critic_corrections.py":
            continue  # this file lists the withdrawn phrases on purpose
        if rel == "CLAIM_REVIEW.md":
            t = t.split("## Critic review")[0]  # the critic section quotes withdrawn phrases on purpose
        for w in WITHDRAWN:
            if w in t:
                stale_hits.append([rel, w])
    # controls: the scanners must fire on planted text and must not fire on the physics term 'anthropic predictions'
    planted = "prepared with " + codecs.decode("Pynhqr", "rot13")
    ctrl_ai = bool(AI_NAME.search(planted)) and not AI_NAME.search("a framework for anthropic predictions")
    ctrl_stale = any(w in "Weyl/radiation ≈ 0.6 after the chart freeze" for w in WITHDRAWN)
    elsewhere = []
    for p in glob.glob(P("**/*"), recursive=True):
        rel = os.path.relpath(p, ROOT)
        if os.path.isfile(p) and "/" in rel and p.endswith((".md", ".py", ".json", ".txt", ".sh")) \
                and "MANIFEST" not in rel:
            try:
                t = open(p, errors="ignore").read()
            except OSError:
                continue
            if AI_NAME.search(t):
                elsewhere.append(rel)
    return dict(ai_names_in_top_level=ai_hits, withdrawn_phrases_present=stale_hits,
                control_ai_regex=ctrl_ai, control_withdrawn_scanner=ctrl_stale,
                ai_names_outside_top_level_not_editable_here=sorted(elsewhere),
                passed=bool(not ai_hits and not stale_hits and ctrl_ai and ctrl_stale))


# ---------------------------------------------------------------- main
def main():
    corpus = lit_corpus()
    log = []
    if not CHECK_ONLY:
        log = apply_json_patches()
    cit = citation_check(corpus)
    if not CHECK_ONLY:
        render_claim_review(cit)
        cit = citation_check(corpus)
    # citation control: a fabricated identifier must be reported missing, a recorded one must not
    ctrl_cit = ("hep-th/9999999" not in corpus) and ("2609.21421" in corpus) and ("hep-th/9912118" in corpus)
    cr = open(P("CLAIM_REVIEW.md")).read()
    listed = cr.split("## Citations")[1].split("## Process notes")[0]
    all_missing = sorted({i for v in cit.values() for i in v["arxiv_ids_without_literature_record"]})
    cit_listed_ok = all(i in listed for i in all_missing)
    checks = dict(
        json_patches_missing=json_patches_present(),
        claim_review_has_critic_section="## Critic review (28 September 2026)" in cr,
        numbers=number_checks(),
        post_freeze_dstar=post_freeze_dstar(),
        wording=wording_scans(),
        citations=dict(per_file=cit, all_missing_arxiv_ids=all_missing, all_listed_in_CLAIM_REVIEW=cit_listed_ok,
                       control_fabricated_id_absent_and_known_ids_present=ctrl_cit,
                       corpus_files=[os.path.relpath(p, ROOT) for p in LIT_FILES if p.endswith((".md", ".json", ".py"))][:60]),
    )
    ok = (not checks["json_patches_missing"] and checks["claim_review_has_critic_section"]
          and all(n["passed"] for n in checks["numbers"]) and checks["post_freeze_dstar"]["passed"]
          and checks["wording"]["passed"] and cit_listed_ok and ctrl_cit)
    out = dict(status="critic checks: numbers are floating-point reads of saved JSON/NPZ outputs (tolerances per entry); "
                      "wording and citation scans are text searches",
               mode="check-only" if CHECK_ONLY else "apply+check", patch_log=log, all_pass=bool(ok), **checks)
    for n in checks["numbers"]:
        print(("PASS " if n["passed"] else "FAIL ") + n["id"], n["desc"], n["value"])
    pf = checks["post_freeze_dstar"]
    print(("PASS" if pf["passed"] else "FAIL"), "post-freeze d*: ratio at window end %.4f; min/max before H<0: %.3f/%.3f; "
          "proper-time advance %.3f H0tau" % (pf["ratio_at_window_end"], pf["ratio_min_before_H_negative"],
                                               pf["ratio_max_before_H_negative"], pf["proper_time_advance_window_end_to_H_negative"]))
    print(("PASS" if checks["wording"]["passed"] else "FAIL"), "wording scans", checks["wording"]["ai_names_in_top_level"],
          checks["wording"]["withdrawn_phrases_present"])
    print(("PASS" if cit_listed_ok and ctrl_cit else "FAIL"), "citations: %d arXiv ids without literature record, all listed"
          % len(all_missing))
    print("json patches missing:", checks["json_patches_missing"])
    print("all_pass:", ok)
    if not CHECK_ONLY:
        with open(P("CRITIC_CHECKS.json"), "w") as f:
            json.dump(out, f, indent=1, ensure_ascii=False)
            f.write("\n")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
