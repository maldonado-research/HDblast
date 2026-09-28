#!/usr/bin/env python3
"""Merge and deduplicate the three literature sweeps of the 27 Sept 2026 round.

Inputs (read-only):
  literature/theory_checks/theory_sources.json          (sweep 1, theory)
  literature/observations_checks/observations_sources.json (sweep 2, observations)
  literature/methods_checks/methods_sources.json        (sweep 3, methods; no new searches ran)

Outputs:
  LITERATURE_2022_2026.md          (checkpoint root)
  synthesis/LITERATURE_MERGED.json (machine-readable merged list + dedup log)

Deduplication key: arXiv identifier if the URL contains one (new or old style,
version suffix stripped), otherwise the lower-cased URL without scheme, 'www.',
trailing slash or query string.  Every URL is copied verbatim from the sweep
records; nothing is added here.  A negative control checks that two distinct
arXiv ids never merge and that a version-suffixed duplicate does merge.
"""
import json
import os
import re
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
LIT = os.path.join(ROOT, "literature")

NEW = re.compile(r"(\d{4}\.\d{4,5})(v\d+)?")
OLD = re.compile(r"((?:hep-th|hep-ph|gr-qc|astro-ph|hep-ex|nucl-th|quant-ph|math-ph)/\d{7})(v\d+)?")


def key_of(url):
    u = url.strip()
    if "arxiv.org" in u:
        m = OLD.search(u) or NEW.search(u)
        if m:
            return "arxiv:" + m.group(1)
    u = re.sub(r"^https?://", "", u.lower())
    u = re.sub(r"^www\.", "", u).split("?")[0].rstrip("/")
    return "url:" + u


def year_int(y):
    try:
        return int(str(y)[:4])
    except ValueError:
        return None


def load(rel):
    with open(os.path.join(LIT, rel)) as f:
        return json.load(f)


# ---------------------------------------------------------------- controls
assert key_of("https://arxiv.org/abs/2511.21362") == key_of("https://arxiv.org/html/2511.21362v2"), "version merge"
assert key_of("https://arxiv.org/abs/hep-th/9912118") == key_of("https://arxiv.org/pdf/hep-th/9912118v3")
assert key_of("https://arxiv.org/abs/2511.21362") != key_of("https://arxiv.org/abs/2511.21363"), "distinct ids"
assert key_of("https://arxiv.org/abs/hep-th/9912118") != key_of("https://arxiv.org/abs/hep-th/9912119")

theory = load("theory_checks/theory_sources.json")["sources"]
obs = load("observations_checks/observations_sources.json")["sources"]
meth = load("methods_checks/methods_sources.json")["sources"]

THEME_T = {"A": "dark bubbles", "B": "creation / nucleation", "C": "collisions / ekpyrosis",
           "D": "RS/KR foundations", "E": "stability of dS walls", "F": "dS-sliced flows / holographic reheating",
           "G": "swampland / strings"}
THEME_M = {"A": "constraint control / boundaries", "B": "shells and junctions", "C": "gauge-invariant perturbations",
           "D": "particle production / backreaction", "E": "validated numerics"}

merged = OrderedDict()
log = []


def add(rec, sweep, topic):
    k = key_of(rec["url"])
    entry = merged.get(k)
    origin = f"{sweep}:{rec['id']}"
    if entry is None:
        merged[k] = entry = {
            "key": k, "title": rec["title"], "authors": rec.get("authors", ""), "year": rec.get("year", ""),
            "venue": rec.get("venue", ""), "url": rec["url"], "urls_seen": [rec["url"]],
            "findings": [rec["finding"]], "relevance": [], "implications": [], "topics": [],
            "access": [], "origins": [],
        }
    else:
        log.append({"merged": origin, "into": entry["origins"][0], "key": k})
        if rec["url"] not in entry["urls_seen"]:
            entry["urls_seen"].append(rec["url"])
        if rec["finding"] not in entry["findings"] and not rec["finding"].startswith("Per the sweep-1"):
            entry["findings"].append(rec["finding"])
        if not entry["authors"] and rec.get("authors"):
            entry["authors"] = rec["authors"]
    entry["origins"].append(origin)
    if rec.get("relevance") and rec["relevance"] not in entry["relevance"]:
        entry["relevance"].append(f"[{sweep}] " + rec["relevance"])
    if rec["implication"] not in entry["implications"]:
        entry["implications"].append(rec["implication"])
    if topic not in entry["topics"]:
        entry["topics"].append(topic)
    if rec["access"] not in entry["access"]:
        entry["access"].append(rec["access"])


for r in theory:
    add(r, "theory", THEME_T.get(r["theme"], r["theme"]))
for r in obs:
    add(r, "observations", r["section"])
for r in meth:
    add(r, "methods", THEME_M.get(r["theme"], r["theme"]))

n_in = len(theory) + len(obs) + len(meth)
entries = list(merged.values())
RANK = ["challenges", "constrains", "supports", "method", "context"]


def primary_impl(e):
    return min(e["implications"], key=lambda x: RANK.index(x) if x in RANK else 99)


recent = [e for e in entries if (year_int(e["year"]) or 0) >= 2022]
unknown = [e for e in entries if year_int(e["year"]) is None]
older = [e for e in entries if (year_int(e["year"]) or 9999) < 2022]

# ------------------------------------------------------------------ markdown
L = []
A = L.append
A("# Literature 2022–2026 for HDBLAST: merged, deduplicated and annotated")
A("")
A("Checkpoint of 27 September 2026 (searches run 28 September 2026). Built by `synthesis/build_literature.py` from the "
  "three sweep records in `literature/`. Every URL below was returned by a web search in this round or is written in the "
  "programme's own archive files; none was added by hand.")
A("")
A("**Read this first: evidence level.** Paper sites (arxiv.org, zenodo.org, journal sites) were blocked for downloads, "
  "so **no paper was opened**. Each entry is based on a search-result snippet (`search-snippet-only`) or on an earlier "
  "programme record (`programme-record-only`). Snippets can misattribute authors, venues or numbers. Before any external "
  "use, open the paper and check. No independent auditor re-checked the literature sweeps in this round; the stability "
  "and mechanism audits could not re-verify citations either, because the shared web-search budget ran out. The methods "
  "sweep ran **no** new searches (budget exhausted), so its entries are second-hand.")
A("")
A(f"**Counts.** {n_in} records in the three sweeps ({len(theory)} theory, {len(obs)} observations, {len(meth)} methods) "
  f"merge to **{len(entries)} distinct sources** ({len(log)} duplicates merged). Of these, **{len(recent)} are dated "
  f"2022–2026**, {len(unknown)} have no year in the snippet, and {len(older)} are earlier foundational work (listed "
  "briefly in the last section).")
A("")
A("**Implication labels** (the sweeps' own judgement, relative to HDBLAST): *challenges* = could undercut the hypothesis "
  "or the matter extension; *constrains* = sets a bound any HDBLAST version must satisfy; *supports* = consistent with, or "
  "provides a mechanism HDBLAST could use (never observational support for HDBLAST itself); *method* = a tool we can use; "
  "*context* = prior art or background.")
A("")
A("## What the 2022–2026 literature implies for HDBLAST (summary)")
A("")
A("1. **Prior art is substantial.** \"Our universe as a wall created in a 5D event\" is an active published line "
  "(dark-bubble cosmology, black-hole-catalysed nucleation, detuned-brane scenarios). HDBLAST must not claim the concept; "
  "its own content is the registered model and its computed results.")
A("2. **Every published hot-origin model adds an ingredient** (a bulk black hole, colliding walls, a holographic hot sector, "
  "or matter on the wall). A pure-tension wall produces no heat. HDBLAST's own negative results this round (preheating "
  "and mechanism screens) agree with that pattern.")
A("3. **Observations constrain, not test, the current model.** The registered 5D model has no radiation era and no "
  "perturbation spectrum, so CMB, BBN, JWST, DESI and H0 data are requirements for the future, not tests now. Dark "
  "radiation (N_eff) is the sharpest direct bound on any 5D remnant (bulk Weyl term).")
A("4. **The PTA knee is untested at likelihood level** and sits at 0.997×(1/yr), on the pulsar-position/proper-motion "
  "fitting blind spot (standard in PTA work, new to this project). 2025–2026 noise reanalyses trend towards γ=13/3, "
  "which is unfavourable to the frozen low-frequency slope but is not a rejection.")
A("5. **Most useful new tools:** the Einstein–scalar de Sitter-sliced flow programme (Kiritsis, Nitti and collaborators), "
  "dark-bubble junction/perturbation methods, and constraint-damping practice from numerical relativity.")
A("")
A("Known inconsistencies between sweep snippets (unresolved, paper not opened): the ACT DR6 N_eff value is quoted as "
  "2.86 ± 0.13 in the observations sweep and 2.89 ± 0.11 in the mechanisms workstream; dark-radiation allowances are "
  "quoted in different normalisations (ρ_dr/ρ_γ vs ρ_DR/ρ_SM at production) and are not directly comparable.")
A("")
A("## Annotated sources dated 2022–2026 (and undated), grouped by implication")
A("")
for impl in RANK:
    group = [e for e in recent + unknown if primary_impl(e) == impl]
    group.sort(key=lambda e: -(year_int(e["year"]) or 0))
    if not group:
        continue
    A(f"### {impl.capitalize()} ({len(group)})")
    A("")
    for e in group:
        au = f" — {e['authors']}" if e["authors"] else ""
        yr = e["year"] if year_int(e["year"]) else "year not given"
        A(f"- **{e['title']}**{au} ({yr}). {e['venue']}. <{e['url']}>")
        A(f"  - *Finding (snippet):* {e['findings'][0]}")
        for extra in e["findings"][1:]:
            A(f"  - *Other sweep's reading:* {extra}")
        for rel in e["relevance"]:
            A(f"  - *For HDBLAST:* {rel}")
        A(f"  - *Labels:* {', '.join(e['implications'])}; topic: {', '.join(e['topics'])}; access: "
          f"{', '.join(e['access'])}; records: {', '.join(e['origins'])}")
    A("")

A("## Earlier foundational work cited by the sweeps (before 2022)")
A("")
A("Listed for completeness; see `literature/theory.md` and `literature/methods.md` for full annotations.")
A("")
A("| Year | Source | Label | Why it matters for HDBLAST (first sentence of the sweep note) |")
A("|---|---|---|---|")
for e in sorted(older, key=lambda e: (year_int(e["year"]), e["title"])):
    rel = e["relevance"][0].split("] ", 1)[-1] if e["relevance"] else ""
    first = re.split(r"(?<=[.!?])\s", rel, maxsplit=1)[0].replace("|", "/")
    A(f"| {e['year']} | [{e['title'].replace('|', '/')}]({e['url']}) | {primary_impl(e)} | {first} |")
A("")
A("## Reproduction")
A("")
A("```bash")
A("cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927")
A("python3 synthesis/build_literature.py   # rewrites this file and synthesis/LITERATURE_MERGED.json")
A("```")
A("")
with open(os.path.join(ROOT, "LITERATURE_2022_2026.md"), "w") as f:
    f.write("\n".join(L))

out = {
    "status": "bibliographic merge (no new searches); all entries search-snippet-only or programme-record-only",
    "n_input_records": n_in, "n_distinct": len(entries), "n_recent_2022_2026": len(recent),
    "n_year_unknown": len(unknown), "n_before_2022": len(older), "duplicates_merged": log,
    "implication_counts_recent": {k: sum(1 for e in recent + unknown if primary_impl(e) == k) for k in RANK},
    "controls": {"version_suffix_merges": True, "distinct_ids_do_not_merge": True},
    "entries": entries,
}
with open(os.path.join(HERE, "LITERATURE_MERGED.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("entries", "duplicates_merged")}, indent=1))
print("duplicates merged:", len(log))
sys.exit(0)
