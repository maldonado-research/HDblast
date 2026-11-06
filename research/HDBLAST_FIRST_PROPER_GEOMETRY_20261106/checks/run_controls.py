"""Six fail-closed controls for check_checkpoint.py; not PDE validation."""
import argparse
import hashlib
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


def manifest(root):
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != root / "MANIFEST.json":
            data = path.read_bytes()
            files[path.relative_to(root).as_posix()] = {
                "bytes": len(data), "sha256": sha(data)}
    value = {"schema": "HDBLAST_FIRST_PROPER_GEOMETRY_MANIFEST_V1", "files": files}
    (root / "MANIFEST.json").write_text(
        json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return sha((root / "MANIFEST.json").read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(script.is_file() and all(not p.is_symlink() for p in
         (script, *script.parents)), "Regular controls script path required")
    root = script.resolve(strict=True).parent.parent
    out = args.output_dir.absolute()
    need(".." not in out.parts and not out.is_relative_to(root),
         "Fresh controls directory must be outside capsule")
    need(not out.exists() and not out.is_symlink(), "Controls directory exists")
    need(out.parent.is_dir() and all(not p.is_symlink() for p in out.parents),
         "Regular existing controls parent required")
    out.mkdir()

    def invoke(case_root, pin, name):
        report = out / (name + ".json")
        command = [sys.executable, "-I", "-S", "-B",
                   str(case_root / "checks/check_checkpoint.py"),
                   "--manifest-sha256", pin, "--output", str(report)]
        result = subprocess.run(command, capture_output=True, text=True)
        return result, report

    # Authenticate the original source before constructing any mutation.
    result, report = invoke(root, args.manifest_sha256, "positive")
    need(result.returncode == 0 and report.is_file(), "Baseline checker failed")
    baseline = json.loads(report.read_text())
    need(baseline["status"] == "PASS", "Baseline report failed")
    cases = []

    def reject(name, case_root, pin):
        result, report = invoke(case_root, pin, name)
        need(result.returncode != 0 and not report.exists(),
             "Control was not rejected: " + name)
        cases.append({"name": name, "exit_code": result.returncode,
                      "report_created": False, "status": "REJECTED",
                      "diagnostic": result.stderr.strip().splitlines()[-1]})

    wrong_pin = "0" * 64 if args.manifest_sha256 != "0" * 64 else "1" * 64
    reject("wrong_manifest_pin", root, wrong_pin)

    def copy(name):
        target = out / "fixtures" / name
        shutil.copytree(root, target)
        return target

    changed = copy("changed_source")
    with (changed / "inputs/HIGHER_JETS_NORMAL.json").open("ab") as stream:
        stream.write(b"\n")
    reject("source_byte_change", changed, args.manifest_sha256)

    extra = copy("extra_member")
    (extra / "UNLISTED.txt").write_text("unexpected source\n")
    reject("unlisted_file", extra, args.manifest_sha256)

    mutant = copy("rehashed_math_mutant")
    evaluator_path = mutant / "evaluate_first_cut.py"
    original = evaluator_path.read_text()
    needle = "g1 = ustar*(2+2*target*active+4*active**2)"
    need(original.count(needle) == 1, "Mathematical control target changed")
    evaluator_path.write_text(original.replace(needle,
        "g1 = 10*ustar*(2+2*target*active+4*active**2)"))
    mutant_pin = manifest(mutant)
    reject("rehashed_metric_budget_error", mutant, mutant_pin)

    linked = copy("symlink_member")
    data = (linked / "inputs/HIGHER_JETS_NORMAL.json").read_bytes()
    target = out / "symlink_input_target.json"
    target.write_bytes(data)
    (linked / "inputs/HIGHER_JETS_NORMAL.json").unlink()
    (linked / "inputs/HIGHER_JETS_NORMAL.json").symlink_to(target)
    reject("symlink_source", linked, args.manifest_sha256)

    alternate = copy("alternate_script_path")
    inside = alternate / "forbidden_report.json"
    alternate_script = alternate / "checks" / ".." / "checks" / "check_checkpoint.py"
    result = subprocess.run(
        [sys.executable, "-I", "-S", "-B", str(alternate_script),
         "--manifest-sha256", args.manifest_sha256, "--output", str(inside)],
        capture_output=True, text=True)
    need(result.returncode != 0 and not inside.exists(),
         "Alternate invocation escaped output containment")
    cases.append({"name": "alternate_script_path_inside_output",
                  "exit_code": result.returncode, "report_created": False,
                  "status": "REJECTED"})

    receipt = {
        "schema": "HDBLAST_FIRST_PROPER_GEOMETRY_CONTROLS_V1", "status": "PASS",
        "manifest_sha256": args.manifest_sha256,
        "baseline_producer_executions": 1, "rejected_children": len(cases),
        "cases": cases,
        "scope": "Protocol and arithmetic failure controls, not independent "
                 "continuous proof or PDE validation.",
    }
    with (out / "CONTROLS.json").open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "rejected_children": len(cases)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
