"""Six CLI rejection controls for the closed capsule; not PDE validation."""
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
    value = {"schema": "HDBLAST_WEIGHTED_REGULARITY_MANIFEST_V1", "files": files}
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

    def invoke(case_root, pin, name, script_path=None, report_path=None):
        report = report_path if report_path is not None else out / (name + ".json")
        checker = script_path if script_path is not None else case_root / "checks/check_checkpoint.py"
        result = subprocess.run(
            [sys.executable, "-I", "-S", "-B", str(checker),
             "--manifest-sha256", pin, "--output", str(report)],
            capture_output=True, text=True, timeout=60)
        return result, report

    # Authenticate the original source before constructing any mutation.
    result, report = invoke(root, args.manifest_sha256, "positive")
    need(result.returncode == 0 and report.is_file(), "Baseline checker failed")
    baseline = json.loads(report.read_text(encoding="utf-8"))
    need(baseline["status"] == "PASS"
         and baseline["manifest_sha256"] == args.manifest_sha256,
         "Baseline report failed")
    cases = []

    def reject(name, case_root, pin, expected, **paths):
        result, report = invoke(case_root, pin, name, **paths)
        lines = result.stderr.strip().splitlines()
        diagnostic = lines[-1] if lines else ""
        need(result.returncode != 0 and not report.exists()
             and diagnostic.startswith("ValueError: ") and expected in diagnostic,
             "Control did not reach its intended rejection: " + name
             + "; diagnostic=" + diagnostic)
        cases.append({"name": name, "exit_code": result.returncode,
                      "report_created": False, "status": "REJECTED",
                      "diagnostic": diagnostic})

    wrong_pin = "0" * 64 if args.manifest_sha256 != "0" * 64 else "1" * 64
    reject("wrong_manifest_pin", root, wrong_pin, "Manifest SHA256 mismatch")

    # The baseline has authenticated these retained actual-input bytes.
    source = Path("inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json")

    def copy(name):
        target = out / "fixtures" / name
        shutil.copytree(root, target)
        return target

    changed = copy("changed_source")
    with (changed / source).open("ab") as stream:
        stream.write(b"\n")
    reject("source_byte_change", changed, args.manifest_sha256,
           "Source bytes differ: " + source.as_posix())

    extra = copy("extra_member")
    (extra / "UNLISTED.txt").write_text("unexpected source\n", encoding="utf-8")
    reject("unlisted_file", extra, args.manifest_sha256, "Closed file roster mismatch")

    mutant = copy("rehashed_math_mutant")
    evaluator_path = mutant / "derive_scalar_envelopes.py"
    original = evaluator_path.read_text(encoding="utf-8")
    needle = "SCALAR_LOWER_JET_CEILING = F(1, 10**6)"
    need(original.count(needle) == 1, "Mathematical control target changed")
    evaluator_path.write_text(original.replace(needle,
        "SCALAR_LOWER_JET_CEILING = F(1, 10**12)"), encoding="utf-8")
    mutant_pin = manifest(mutant)
    reject("rehashed_scalar_envelope_arithmetic_error", mutant, mutant_pin,
           "scalar normal lower coefficient jets lambda0")

    linked = copy("symlink_member")
    target = out / "symlink_input_target.json"
    target.write_bytes((linked / source).read_bytes())
    (linked / source).unlink()
    (linked / source).symlink_to(target)
    reject("symlink_source", linked, args.manifest_sha256, "Capsule symlink rejected")

    alternate = copy("alternate_script_path")
    reject("alternate_script_path_inside_output", alternate, args.manifest_sha256,
           "Output must be outside capsule",
           script_path=alternate / "checks" / ".." / "checks" / "check_checkpoint.py",
           report_path=alternate / "forbidden_report.json")

    need(len(cases) == 6, "Six rejection controls required")
    receipt = {
        "schema": "HDBLAST_WEIGHTED_REGULARITY_CONTROLS_V1", "status": "PASS",
        "manifest_sha256": args.manifest_sha256,
        "baseline_producer_executions": 1, "rejected_children": len(cases),
        "cases": cases,
        "scope": "Protocol and finite arithmetic failure controls, not independent "
                 "continuous proof, nonlinear amplitude bounds or PDE validation.",
    }
    with (out / "CONTROLS.json").open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "rejected_children": len(cases)},
                     sort_keys=True))


if __name__ == "__main__":
    main()
