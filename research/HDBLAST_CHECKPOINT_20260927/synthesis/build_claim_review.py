#!/usr/bin/env python3
"""Render CLAIM_REVIEW.md from RESULTS_SUMMARY.json (single source of truth).

Run after build_results_summary.py.  The narrative sections are fixed text; the
table rows come only from the JSON so the two files cannot drift apart.
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
with open(os.path.join(ROOT, "RESULTS_SUMMARY.json")) as f:
    R = json.load(f)

WS_TITLE = {
    "stability_gauge_invariant": "Stability of the +1 branch, Method A (gauge-invariant perturbation theory)",
    "stability_time_domain": "Stability of the +1 branch, Method B (time domain and discrete operator)",
    "analytic_structure": "Exact analytic structure of the +1 branch",
    "preheating": "Instant preheating along the archived roll-off",
    "evolution_constraints": "Constraint control in 5D evolutions (numerical method)",
    "mechanisms": "Screen of mechanisms for a radiation era",
    "literature": "Literature sweeps (theory, observations, methods)",
    "program_map": "Program map",
    "synthesis": "Synthesis checks done for this report",
}


def esc(s):
    return s.replace("|", "/").replace("\n", " ")


L = []
A = L.append
A("# Claim-by-claim review — HDBLAST checkpoint, 27 September 2026")
A("")
A("Prepared with Claude (Anthropic) for Ricardo Maldonado's HDBLAST program. This is an internal review, not external "
  "peer review. The table is generated from `RESULTS_SUMMARY.json` by `synthesis/build_claim_review.py`.")
A("")
A("## How to read this table")
A("")
A("- **Producer label**: what the workstream itself called the result.")
A("- **Verifier verdict**: the independent auditor's verdict (each workstream folder has a `VERIFICATION.md` and a "
  "`verify/` folder with the auditor's own scripts). *not-independently-audited* means only the workstream's own "
  "checks exist (literature, program map, and the synthesis checks done for this report).")
A("- **Status here**: the label this checkpoint uses after applying the verifier's corrections. "
  "*corrected* means the original claim was refuted and the corrected statement is given instead.")
A("- **Established** (✔) only when the verifier confirmed the claim and it is not inconclusive or corrected.")
A("")
A("Status vocabulary: exact-verified (symbolic identity), numerical (floating point with measured tolerances, never "
  "interval-certified here), conditional (valid under stated assumptions), negative (the idea fails the test), "
  "inconclusive, corrected, context, procedural-negative.")
A("")
A(f"Totals: {R['n_claims']} claims. Verdicts: " + ", ".join(f"{k} {v}" for k, v in R["counts"]["verifier_verdict"].items())
  + ". Status: " + ", ".join(f"{k} {v}" for k, v in R["counts"]["status"].items()) + ".")
A("")
A("## Refuted or corrected claims (read these first)")
A("")
for c in R["claims"]:
    if c["verifier_verdict"] == "refuted":
        A(f"- **{c['id']}** — original: *{c['claim']}* → **Correction:** {c['correction']}")
A("")
A("## Other corrections from the audits (partially confirmed)")
A("")
for c in R["claims"]:
    if c["verifier_verdict"] == "partially-confirmed" and c["correction"]:
        A(f"- **{c['id']}**: {c['correction']}")
A("")
ws_order = []
for c in R["claims"]:
    if c["workstream"] not in ws_order:
        ws_order.append(c["workstream"])
for ws in ws_order:
    A(f"## {WS_TITLE.get(ws, ws)} (`{ws}/`)")
    A("")
    A("| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |")
    A("|---|---|---|---|---|---|---|---|")
    for c in R["claims"]:
        if c["workstream"] != ws:
            continue
        ev = "<br>".join(f"`{p}`" for p in c["evidence"])
        cn = " ".join(x for x in (c["correction"], c["note"]) if x)
        A(f"| {c['id']} | {esc(c['claim'])} | {c['producer_label']} | {c['verifier_verdict']} | {c['status']} | "
          f"{'✔' if c['established'] else ''} | {ev} | {esc(cn)} |")
    A("")
A("## Reconciliation of the two stability methods")
A("")
A("| Question | Method A (gauge-invariant) | Method B (time domain) | Agree? | Why / scope |")
A("|---|---|---|---|---|")
A("| Calibration growth of the original shell | 1.657193631259 | 1.657193368 (pencil), 1.6571936573 (eigenvalue) | yes, 1.6e-7 / 1.6e-8 | same physical tachyon; B's residual is time-step/grid error |")
A("| Artificial tachyon (sign of σ″ flipped), δ=0.001 | μ²=−28555.115 → rate 167.4892463 | rate 167.4892463 (shift-invert) | yes, 3.7e-10 | independent codes see the same wrong-model instability, so both can detect one |")
A("| Scalar bound states below the 9/4 continuum on the +1 branch | none at 6 δ | none (homogeneous sector) at δ=0.001, 0.003, 0.01 | yes | B covers only the homogeneous (FRW-symmetric) sector; every dS₄ harmonic μ² has a homogeneous representative, so both test the same eigenvalue condition |")
A("| Slowest decay | continuum edge μ²=9/4 ⇒ e^{−3H₊τ/2} | −1.500 ± 0.001 in H₊ units | yes | = −0.910 H₀ |")
A("| Tensor sector | only the massless graviton, gap 3H/2 | not tested | — | analytic_structure's leading-order spectrum {0} ∪ [9/4, ∞) agrees with A |")
A("| Positive discrete rates | none physical | gauge mode (λ→1), far-boundary artefacts, grid-scale modes | no conflict | B's positive rates have no shell support and change with L or h; they are properties of the truncated grid |")
A("| The λ≈1 mode | ℓ=1 harmonic μ²=−4 is pure gauge (no physical mode since Bφ′_b≠0) | residual gauge mode, λ=0.99999993 at L=8 | consistent | λ=1 ⇔ μ²=−λ(λ+3)=−4; mapped value agrees with −4 to <1e-6 |")
A("| δ→0 limit | B-coefficient → 3.5555554 (fit) | — | matches exact 32/9 to 3e-8 | new cross-check (`synthesis/RECONCILIATION.json`) |")
A("")
A("Overall: the two methods agree wherever they overlap, and their scopes are complementary. Neither is "
  "interval-certified; the vector sector, nonlinear stability, quantum tunnelling and dynamical attraction remain open.")
A("")
A("## Citations")
A("")
A("Paper sites were blocked for downloads and the shared web-search budget ran out during the audits, so no auditor "
  "could re-verify any citation. All citations in this round are search-snippet-only or taken from earlier programme "
  "records (`LITERATURE_2022_2026.md`). A few links used as context in workstream READMEs (Karch–Randall hep-th/0011156, "
  "Himemoto–Sasaki gr-qc/0010035, Langlois–Maartens–Sasaki–Wands hep-th/0012044, arXiv:2609.21421) have no search record "
  "in this round and are marked unverifiable. No result depends on them.")
A("")
A("## Process notes")
A("")
A("- The analytic_structure auditor accidentally regenerated that folder's `MANIFEST.sha256.json`, adding 25 `verify/` "
  "entries. `synthesis/QUOTED_NUMBER_CHECKS.json` confirms that all 29 original entries still match the current files. "
  "The file was left as found; the root `MANIFEST.sha256.json` covers every file.")
A("- No claim in this checkpoint is a discovery, a proof, a hot Big Bang mechanism, or observational support for the "
  "hypothesis. Novelty is claimed only relative to this project.")
A("")
with open(os.path.join(ROOT, "CLAIM_REVIEW.md"), "w") as f:
    f.write("\n".join(L))
print("wrote CLAIM_REVIEW.md with", R["n_claims"], "claims")
