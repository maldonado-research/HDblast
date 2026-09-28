#!/usr/bin/env python3
"""Validate the theory-sweep source records and render ../theory.md and theory_sources.json.

Checks (the build fails if any fails):
  - at least 15 distinct queries recorded;
  - every source has id, title, year, url, finding, relevance, a valid implication and access label;
  - every cited query number exists; URLs are unique; ids are unique;
  - when a venue string names an arXiv id and the URL is an arxiv.org link, the two ids agree;
  - the cross-check JSON exists and reports all_pass = true.
"""
import json, os, re, sys, importlib.util, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location('ts', os.path.join(HERE, 'theory_sources.py'))
ts = importlib.util.module_from_spec(spec); spec.loader.exec_module(ts)

ALLOWED_IMP = {"supports", "constrains", "method", "context", "challenges"}
ALLOWED_ACC = {"search-snippet-only", "full-text-read"}
errors = []
Q = ts.QUERIES
if len(set(Q)) < 15:
    errors.append("fewer than 15 distinct queries")
ids, urls = set(), set()
for s in ts.SOURCES:
    for k in ("id", "title", "year", "url", "finding", "relevance", "implication", "access", "queries"):
        if not s.get(k):
            errors.append(f"{s.get('id')}: missing {k}")
    if s["id"] in ids: errors.append(f"duplicate id {s['id']}")
    ids.add(s["id"])
    if s["url"] in urls: errors.append(f"duplicate url {s['url']}")
    urls.add(s["url"])
    if s["implication"] not in ALLOWED_IMP: errors.append(f"{s['id']}: bad implication")
    if s["access"] not in ALLOWED_ACC: errors.append(f"{s['id']}: bad access")
    for q in s["queries"]:
        if not (1 <= q <= len(Q)): errors.append(f"{s['id']}: bad query {q}")
    m_url = re.search(r'arxiv\.org/(?:abs|pdf|html)/([^\s?#]+?)(?:v\d+)?$', s["url"])
    m_ven = re.findall(r'arXiv:([\w\-/\.]+)', s.get("venue", ""))
    if m_url and m_ven:
        uid = m_url.group(1)
        if not any(uid == v.rstrip('.,;)') for v in m_ven):
            errors.append(f"{s['id']}: arXiv id in url ({uid}) not in venue {m_ven}")
    if not re.match(r'^https://', s["url"]): errors.append(f"{s['id']}: url not https")

xc_path = os.path.join(HERE, 'rs_thin_brane_crosscheck.json')
if not os.path.exists(xc_path):
    errors.append("cross-check JSON missing: run rs_thin_brane_crosscheck.py first")
    xc = None
else:
    xc = json.load(open(xc_path))
    if not xc.get("all_pass"): errors.append("cross-check all_pass is false")

if errors:
    print("BUILD FAILED"); print("\n".join(errors)); sys.exit(1)

pm_path = os.path.join(HERE, 'prior_mentions.json')
if not os.path.exists(pm_path):
    print("BUILD FAILED: run prior_mention_scan.py first"); sys.exit(1)
PMJ = json.load(open(pm_path)); PM = PMJ["results"]
if set(PM) != {s["id"] for s in ts.SOURCES}:
    print("BUILD FAILED: prior_mentions.json out of date; rerun prior_mention_scan.py"); sys.exit(1)
n_prior = sum(1 for v in PM.values() if v["n_files"] > 0)
theme_names = {"A": "A. Dark bubbles: de Sitter on a wall from decaying AdS5",
               "B": "B. Creation, nucleation, catalysis and bubbles of nothing",
               "C": "C. Brane and bubble collisions; ekpyrotic and cyclic models",
               "D": "D. Randall–Sundrum / Karch–Randall foundations and brane cosmology",
               "E": "E. Stability of de Sitter-sliced walls and branes",
               "F": "F. Einstein–scalar de Sitter-sliced flows; holographic reheating",
               "G": "G. Swampland and string-theory context"}
bib = []
for th in "ABCDEFG":
    bib.append(f"### {theme_names[th]}\n")
    for s in [x for x in ts.SOURCES if x["theme"] == th]:
        qs = ", ".join(str(q) for q in s["queries"])
        cav = f"\n- *Caveat:* {s['caveat']}" if s.get("caveat") else ""
        pm = PM.get(s["id"], {})
        prior = (f"yes, in {pm['n_files']} scanned programme file(s), e.g. `{pm['files'][0]}`" if pm.get("n_files")
                 else "not found in the scanned programme files")
        bib.append(
            f"**{s['id']}. {s['title']}** ({s['year']})\n"
            f"- Authors: {s['authors']}\n"
            f"- Venue or identifier: {s['venue']}\n"
            f"- URL: <{s['url']}>\n"
            f"- Access: {s['access']} (queries {qs}). Implication: **{s['implication']}**.\n"
            f"- *Finding (from snippet):* {s['finding']}\n"
            f"- *Relevance to HDBLAST:* {s['relevance']}{cav}\n"
            f"- *Already in the programme record:* {prior}\n")
bib_txt = "\n".join(bib)

imp = collections.Counter(s["implication"] for s in ts.SOURCES)
def period(y):
    y = int(re.match(r'\d{4}', y).group(0))
    return "2022–2026" if y >= 2022 else ("2010–2021" if y >= 2010 else "before 2010")
per = collections.Counter(period(s["year"]) for s in ts.SOURCES)
imp_txt = ", ".join(f"{k} {imp[k]}" for k in ["supports", "constrains", "challenges", "method", "context"])
per_txt = ", ".join(f"{k}: {per[k]}" for k in ["2022–2026", "2010–2021", "before 2010"])

n = xc["numbers_at_registered_point"]; c = xc["controls"]
rows = [
    ("AdS radius at φ=+1 and φ=−1 (from U=−6/ℓ²)", f"ℓ = {xc['phi=1']['ell']} and {xc['phi=-1']['ell']}", "EXACT"),
    ("Balanced tension 2W equals the critical tension 6/ℓ", f"{xc['phi=1']['balanced_tension_2W']} = {xc['phi=1']['critical_tension_6_over_ell']}; {xc['phi=-1']['balanced_tension_2W']} = {xc['phi=-1']['critical_tension_6_over_ell']}", "EXACT"),
    ("dS-sliced AdS5 with ρ=ℓ sinh(y/ℓ) is Einstein with R_AB=−(4/ℓ²)g_AB; regular cone ρ′(0)", f"{xc['dS_sliced_AdS5_is_Einstein_R_eq_minus4_over_ell2_g']}; ρ′(0) = {xc['regular_cone_rho_prime_at_0']}", "EXACT"),
    ("First-order flow φ′=W_φ, A′=−W/3 solves R_AB=φ_Aφ_B+(2/3)U g_AB", str(xc["first_order_flow_solves_Einstein_eqs_(DeWolfe_et_al_structure)"]), "EXACT"),
    ("Flat kink φ=−tanh y solves φ′=W_φ", str(xc["flat_kink_phi_eq_minus_tanh_solves_phi_prime_eq_W_phi"]), "EXACT"),
    ("U for large φ (either sign)", f"U ~ ({xc['U_large_phi_leading_coefficient_of_phi6']}) φ⁶, unbounded below", "EXACT"),
    ("Thin-brane H²=(σ(1)/6)²−1/81 minus checkpoint expansion, through O(δ²)", f"checkpoint − thin = {xc['checkpoint_minus_thin_brane_through_delta2']}", "EXACT"),
    ("Thin-brane H² vs checkpoint `metric_only_H2` at δ=10⁻³", f"{n['H2_thin_brane'][:16]} vs {n['checkpoint_metric_only_H2']}; relative difference {n['rel_diff_thin_vs_checkpoint_metric_only']}", "NUMERICAL"),
    ("H/H₀: thin-brane vs full +1 branch", f"{n['H_over_H0_thin_brane']} vs {n['H_over_H0_full_branch']}", "NUMERICAL"),
    ("Shift of H, full vs thin (ppm); second-order prediction −c²δ²/384", f"{n['ppm_shift_in_H_full_vs_thin']} vs {n['predicted_ppm_from_c2_over_384']}", "NUMERICAL"),
    ("KK continuum threshold (units of H²)", xc["KK_threshold_over_H2"], "EXACT (literature value)"),
    ("Control: one-sided junction σ/3 matches?", str(c["one_sided_junction_sigma_over_3_matches"]), "control, expected False"),
    ("Control: φ=−1 radius with φ=+1 tension matches?", str(c["wrong_vacuum_radius_matches"]), "control, expected False"),
    ("Control: wrong-sign flow A′=+W/3 solves Einstein eqs?", str(c["wrong_sign_flow_A_prime_eq_plus_W_over_3_solves"]), "control, expected False"),
    ("Control: c perturbed by 1%, relative difference from checkpoint", c["perturbed_c_rel_diff_vs_checkpoint_metric_only"], "control, expected ≫10⁻¹²"),
    ("All checks and controls", str(xc["all_pass"]), ""),
]
res_txt = "| Check | Result | Label |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {l} |" for a, b, l in rows)
qt = "| # | Query |\n|---|---|\n" + "\n".join(f"| {i+1} | {q.replace('|', '/')} |" for i, q in enumerate(Q))
exc = "\n".join(f"- Query {e['query']}: {e['reason']}" for e in ts.EXCLUDED)

tpl = open(os.path.join(HERE, 'theory_template.md')).read()
out = (tpl.replace("{{N_QUERIES}}", str(len(Q))).replace("{{N_SOURCES}}", str(len(ts.SOURCES)))
          .replace("{{IMPLICATION_COUNTS}}", imp_txt).replace("{{PERIOD_COUNTS}}", per_txt)
          .replace("{{INPUT_SHA}}", xc["input_sha256"]).replace("{{PRIOR_FILES}}", str(PMJ["files_scanned"]))
          .replace("{{N_PRIOR}}", str(n_prior)).replace("{{N_NEW}}", str(len(ts.SOURCES) - n_prior)).replace("{{RESULTS_TABLE}}", res_txt)
          .replace("{{BIBLIOGRAPHY}}", bib_txt).replace("{{EXCLUDED}}", exc).replace("{{QUERY_TABLE}}", qt))
assert "{{" not in out
for s in ts.SOURCES:
    assert f"**{s['id']}." in out
open(os.path.join(ROOT, 'theory.md'), 'w').write(out)
json.dump(dict(queries=Q, sources=ts.SOURCES, excluded=ts.EXCLUDED,
               counts=dict(queries=len(Q), sources=len(ts.SOURCES), implication=dict(imp), period=dict(per),
                           already_in_programme_record=n_prior, not_found_in_programme_record=len(ts.SOURCES) - n_prior),
               prior_mentions={k: dict(n_files=v["n_files"], files=v["files"]) for k, v in PM.items()},
               access_note="All entries search-snippet-only; arxiv/zenodo/springer/diva/uu.se fetches returned EGRESS_BLOCKED."),
          open(os.path.join(HERE, 'theory_sources.json'), 'w'), indent=1, ensure_ascii=False)
print(f"BUILD OK: {len(Q)} queries, {len(ts.SOURCES)} sources; implication {dict(imp)}; period {dict(per)}")
