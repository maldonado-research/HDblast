"""Reject freshly manifested mathematical claim mutations with the full checker."""
import argparse
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(condition, message):
    if not condition:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    need(spec is not None and spec.loader is not None, "Control loader unavailable")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def rehash(root):
    files = {}
    for path in sorted(root.rglob("*")):
        need(not path.is_symlink(), "Control fixture symlink rejected")
        if path.is_file() and path != root / "MANIFEST.json":
            raw = path.read_bytes()
            files[path.relative_to(root).as_posix()] = {"bytes": len(raw), "sha256": sha(raw)}
    raw = (json.dumps({"schema": "HDBLAST_UNIFORM_CAUCHY_DATA_MANIFEST_V1",
                       "files": files}, sort_keys=True, indent=2) + "\n").encode("utf-8")
    (root / "MANIFEST.json").write_bytes(raw)
    return sha(raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    script = Path(__file__).absolute()
    need(script.is_file() and all(not item.is_symlink() for item in (script, *script.parents)),
         "Regular mathematical-controls path required")
    root = script.resolve(strict=True).parent.parent
    out = args.output_dir.absolute()
    need(".." not in out.parts and not out.is_relative_to(root),
         "Fresh mathematical-controls directory must be outside capsule")
    need(not out.exists() and not out.is_symlink(), "Mathematical-controls directory exists")
    need(out.parent.is_dir() and all(not item.is_symlink() for item in out.parents),
         "Regular existing mathematical-controls parent required")
    out.mkdir()

    def invoke(case_root, pin, name):
        report = out / (name + ".json")
        child = subprocess.run([sys.executable, "-I", "-B",
            str(case_root / "checks/check_checkpoint.py"), "--manifest-sha256", pin,
            "--output", str(report)], capture_output=True, text=True, timeout=300)
        return child, report

    child, report = invoke(root, args.manifest_sha256, "positive")
    need(child.returncode == 0 and report.is_file(), "Mathematical-controls baseline failed")
    baseline = json.loads(report.read_text(encoding="utf-8"))
    need(baseline["status"] == "PASS" and baseline["manifest_sha256"] == args.manifest_sha256,
         "Mathematical-controls baseline report failed")
    # The baseline closed checker has authenticated all modules before these imports.
    modules={suite:module(root/"checks"/(name+".py"),"mutations_"+suite)
             for suite,name in (("arithmetic","cauchy_arithmetic"),("chart","chart_chain"),("algebra","cauchy_algebra"))}
    data=json.loads((root/"DATA.json").read_text(encoding="utf-8"))
    mutations=[]
    for suite,loaded in modules.items():
        need(hasattr(loaded,"MUTATIONS"),"Mutation roster required: "+suite)
        mutations.extend((suite,field,value) for field,value in loaded.MUTATIONS)
    need(len({(suite, field) for suite, field, value in mutations}) == len(mutations),
         "Repeated mathematical mutation")
    cases = []
    for suite, field, value in mutations:
        name = "rehashed_" + suite + "_" + field
        fixture = out / "fixtures" / name
        shutil.copytree(root, fixture)
        mutant = copy.deepcopy(data)
        before = mutant[suite][field]
        need(before != value, "Mutation leaves claim unchanged: " + name)
        mutant[suite][field] = value
        (fixture / "DATA.json").write_text(json.dumps(mutant, sort_keys=True, indent=2) + "\n",
                                             encoding="utf-8")
        pin = rehash(fixture)
        need(pin != args.manifest_sha256, "Mathematical fixture was not freshly manifested")
        child, child_report = invoke(fixture, pin, name)
        lines = child.stderr.strip().splitlines()
        diagnostic = lines[-1] if lines else ""
        need(child.returncode != 0 and not child_report.exists()
             and diagnostic.startswith("ValueError: ")
             and ("Arithmetic failed: " + suite + ":" in diagnostic
                  or "oracle" in diagnostic or "reject" in diagnostic),
             "Mathematical fault did not reach an exact oracle: " + name + "; diagnostic=" + diagnostic)
        cases.append({"name": name, "suite": suite, "field": field,
            "before": before, "after": value, "rehashed_manifest_sha256": pin,
            "status": "REJECTED", "report_created": False, "diagnostic": diagnostic})
    receipt = {"schema": "HDBLAST_UNIFORM_CAUCHY_DATA_MATH_CONTROLS_V1",
        "status": "PASS", "manifest_sha256": args.manifest_sha256,
        "baseline_producer_executions": 1, "rejected_children": len(cases),
        "suite_counts": {suite: sum(case["suite"] == suite for case in cases)
                         for suite in ("arithmetic", "chart", "algebra")},
        "cases": cases,
        "scope": "Every mathematical claim mutation receives a new complete capsule manifest "
                 "and is rejected by a fresh full-checker process after successful source authentication. "
                 "Continuous analytic arguments are not validated by these finite controls."}
    with (out / "MATH_CONTROLS.json").open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "rejected_children": len(cases)}, sort_keys=True))


if __name__ == "__main__":
    main()
