#!/usr/bin/env python3
"""Read-only scan: does any project text file discuss the 1/yr timing-model
(sky-position / proper-motion / parallax) degeneracy of PTA spectra?

Scans .md/.txt/.json files (< 20 MB) under new-files/D-Blast 3 and
new-files/latest-work.  Writes project_scan_1yr.json.  Positive control: the
pattern 'one cycle per year' must be found (the project page uses it).
"""
import json
import re
import time
from pathlib import Path

ROOT = Path("/home/user/unified-theory-maldonado/new-files")
DIRS = [ROOT / "D-Blast 3", ROOT / "latest-work"]
EXT = {".md", ".txt", ".json"}
PATTERNS = {
    "degeneracy": re.compile(r"sky[- ]position fit|position fitting|astrometric fit|proper[- ]motion fit|"
                             r"parallax.{0,40}(1/yr|2/yr)|(1/yr|one cycle per year).{0,80}(sky position|astrometr|timing[- ]model fit)",
                             re.I),
    "positive_control_one_cycle_per_year": re.compile(r"one cycle per year", re.I),
}
t0 = time.time()
hits = {k: [] for k in PATTERNS}
n_files = 0
for d in DIRS:
    for p in d.rglob("*"):
        if p.suffix.lower() not in EXT or not p.is_file():
            continue
        try:
            if p.stat().st_size > 20_000_000:
                continue
            txt = p.read_text(errors="ignore")
        except OSError:
            continue
        n_files += 1
        for k, rx in PATTERNS.items():
            if rx.search(txt):
                hits[k].append(str(p.relative_to(ROOT)))
        if time.time() - t0 > 600:
            break
out = {"files_scanned": n_files, "seconds": round(time.time() - t0, 1),
       "n_hits": {k: len(v) for k, v in hits.items()},
       "degeneracy_hits": hits["degeneracy"][:50],
       "positive_control_examples": hits["positive_control_one_cycle_per_year"][:5],
       "positive_control_passed": len(hits["positive_control_one_cycle_per_year"]) > 0}
Path(__file__).with_suffix(".json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ["files_scanned", "seconds", "n_hits", "positive_control_passed"]}, indent=1))
