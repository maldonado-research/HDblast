"""Five independent closed-capsule authentication rejection controls.

Mathematical mutation controls are supplied separately by the replay modules.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(script.is_file() and all(not item.is_symlink() for item in (script, *script.parents)),
         "Regular controls script path required")
    root = script.resolve(strict=True).parent.parent
    out = args.output_dir.absolute()
    need(".." not in out.parts and not out.is_relative_to(root),
         "Fresh controls directory must be outside capsule")
    need(not out.exists() and not out.is_symlink(), "Controls directory exists")
    need(out.parent.is_dir() and all(not item.is_symlink() for item in out.parents),
         "Regular existing controls parent required")
    out.mkdir()

    def invoke(case_root, pin, name):
        report = out / (name + ".json")
        result = subprocess.run([sys.executable, "-I", "-B",
            str(case_root / "checks/check_checkpoint.py"), "--manifest-sha256", pin,
            "--output", str(report)], capture_output=True, text=True, timeout=300)
        return result, report

    result, report = invoke(root, args.manifest_sha256, "positive")
    need(result.returncode == 0 and report.is_file(), "Baseline checker failed")
    baseline = json.loads(report.read_text(encoding="utf-8"))
    need(baseline["status"] == "PASS" and baseline["manifest_sha256"] == args.manifest_sha256,
         "Baseline report failed")
    cases = []

    def copy(name):
        target = out / "fixtures" / name
        shutil.copytree(root, target)
        return target

    def reject(name, fixture, pin, expected):
        child, child_report = invoke(fixture, pin, name)
        lines = child.stderr.strip().splitlines()
        diagnostic = lines[-1] if lines else ""
        need(child.returncode != 0 and not child_report.exists()
             and diagnostic.startswith("ValueError: ") and expected in diagnostic,
             "Control did not reach intended rejection: " + name + "; diagnostic=" + diagnostic)
        cases.append({"name": name, "exit_code": child.returncode,
                      "report_created": False, "status": "REJECTED", "diagnostic": diagnostic})

    source = Path("inputs/CURRENT_SEED_DATA.json")
    changed = copy("changed_source")
    with (changed / source).open("ab") as stream:
        stream.write(b"\n")
    reject("changed_source", changed, args.manifest_sha256, "Source bytes differ: " + source.as_posix())

    extra = copy("extra_file")
    (extra / "UNLISTED.txt").write_text("unlisted\n", encoding="utf-8")
    reject("extra_file", extra, args.manifest_sha256, "Closed file roster mismatch")

    missing = copy("missing_file")
    (missing / source).unlink()
    reject("missing_file", missing, args.manifest_sha256, "Closed file roster mismatch")

    linked = copy("symlink_source")
    target = out / "symlink_input_target.json"
    target.write_bytes((linked / source).read_bytes())
    (linked / source).unlink()
    (linked / source).symlink_to(target)
    reject("symlink_source", linked, args.manifest_sha256, "Capsule symlink rejected")

    stale = copy("stale_manifest_pin")
    with (stale / "MANIFEST.json").open("ab") as stream:
        stream.write(b"\n")
    reject("stale_manifest_pin", stale, args.manifest_sha256, "Manifest SHA256 mismatch")

    need(len(cases) == 5, "Five authentication rejection controls required")
    receipt = {"schema": "HDBLAST_VARIABLE_COEFFICIENT_ENERGY_AUTH_CONTROLS_V1",
        "status": "PASS", "manifest_sha256": args.manifest_sha256,
        "baseline_checker_executions": 1, "rejected_children": len(cases), "cases": cases,
        "scope": "Five authentication faults rejected before evaluator import. "
                 "These controls do not validate continuous analytic proofs or "
                 "replace separately rehashed mathematical mutation controls."}
    with (out / "AUTH_CONTROLS.json").open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "rejected_children": len(cases)}, sort_keys=True))


if __name__ == "__main__":
    main()
