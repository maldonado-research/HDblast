#!/usr/bin/env python3
"""Negative control for build_methods_md.py: corrupted copies of the records must be REFUSED.

Each case writes a corrupted copy of methods_sources.py to a temporary directory and runs the builder on it
with --out pointing to that temporary directory (methods.md in this folder is never touched). Pass = the
builder exits non-zero AND its message names the expected failure. The unmodified records must build.
"""
import json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, "methods_sources.py")).read()
CASES = {
    "fabricated_programme_url": ('url="https://arxiv.org/abs/2504.21644v2"', 'url="https://arxiv.org/abs/2504.21645v2"',
                                 "not found verbatim"),
    "altered_sweep1_url": ('url="https://arxiv.org/abs/2006.04999"', 'url="https://arxiv.org/abs/2006.04998"',
                           "differs from sweep-1 record"),
    "lead_with_url": ('why="Left open; not cited."', 'why="Left open; see https://example.org/gregory-ruth"',
                      "lead carries a URL"),
    "access_label_upgraded": ('id="E1", theme="E", provenance="programme",\n         programme_file="untitled folder 136/M489G_LITERATURE_AND_PUBLICATION_REVIEW_20260903.md",\n         programme_access="publisher full text and versioned preprint, read 3 Sept 2026 by an earlier agent", access="programme-record-only"',
                              'id="E1", theme="E", provenance="programme",\n         programme_file="untitled folder 136/M489G_LITERATURE_AND_PUBLICATION_REVIEW_20260903.md",\n         programme_access="publisher full text and versioned preprint, read 3 Sept 2026 by an earlier agent", access="full-text-read"',
                              "inconsistent with provenance"),
}
results = {}
with tempfile.TemporaryDirectory() as td:
    # unmodified records must build
    p0 = os.path.join(td, "ok_sources.py"); open(p0, "w").write(SRC)
    r0 = subprocess.run([sys.executable, os.path.join(HERE, "build_methods_md.py"), "--sources", p0, "--out", os.path.join(td, "ok")],
                        capture_output=True, text=True)
    results["unmodified_builds"] = r0.returncode == 0
    for name, (old, new, expect) in CASES.items():
        if old not in SRC:
            results[name] = {"passed": False, "reason": "pattern not found in records (control could not be applied)"}
            continue
        p = os.path.join(td, name + ".py"); open(p, "w").write(SRC.replace(old, new, 1))
        r = subprocess.run([sys.executable, os.path.join(HERE, "build_methods_md.py"), "--sources", p, "--out", os.path.join(td, name)],
                           capture_output=True, text=True)
        results[name] = {"passed": r.returncode != 0 and expect in r.stdout, "returncode": r.returncode,
                         "message": r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-200:]}
all_pass = results["unmodified_builds"] and all(v["passed"] for k, v in results.items() if k != "unmodified_builds")
out = {"all_pass": all_pass, "results": results}
json.dump(out, open(os.path.join(HERE, "builder_negative_control.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(0 if all_pass else 1)
