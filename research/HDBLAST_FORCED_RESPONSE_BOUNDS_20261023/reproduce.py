#!/usr/bin/env python3
"""Authenticate the capsule and retained dependency before twelve replay children."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        need(key not in result, "DUPLICATE_JSON_KEY:"+key)
        result[key] = value
    return result


def decode(raw):
    return json.loads(raw, object_pairs_hook=unique_object)


def inventory(expected):
    manifest_path = ROOT/"MANIFEST.json"
    need(manifest_path.is_file() and not manifest_path.is_symlink(), "MANIFEST_NOT_REGULAR")
    raw = manifest_path.read_bytes()
    need(sha(raw) == expected, "MANIFEST_EXTERNAL_PIN")
    manifest = decode(raw)
    need(manifest["capsule"] == "HDBLAST_FORCED_RESPONSE_BOUNDS_20261023", "CAPSULE_NAME")
    rows = manifest["files"]
    expected_names = set()
    for row in rows:
        name = row["path"]
        pure = PurePosixPath(name)
        need(name == pure.as_posix() and not pure.is_absolute() and name not in (".", "MANIFEST.json")
             and ".." not in pure.parts and "\\" not in name, "NONCANONICAL_MANIFEST_PATH")
        need(name not in expected_names, "DUPLICATE_MANIFEST_PATH")
        expected_names.add(name)
        path = ROOT/name
        need(path.is_file() and not path.is_symlink(), "NONREGULAR_MEMBER:"+name)
        data = path.read_bytes()
        need(len(data) == row["bytes"] and sha(data) == row["sha256"], "MEMBER_PIN:"+name)
    actual_names = set()
    for path in ROOT.rglob("*"):
        need(not path.is_symlink(), "SOURCE_SYMLINK")
        need(path.is_file() or path.is_dir(), "SOURCE_SPECIAL_FILE")
        if path.is_file():
            actual_names.add(path.relative_to(ROOT).as_posix())
    need(actual_names == expected_names | {"MANIFEST.json"}, "SOURCE_INVENTORY_MISMATCH")
    return {"manifest_sha256": expected, "files_including_manifest": len(actual_names),
            "members": rows}


def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode,
         "USE_PYTHON_I_S_B")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    need(len(args.manifest_sha256) == 64 and all(c in "0123456789abcdef" for c in args.manifest_sha256),
         "INVALID_MANIFEST_HASH")
    output = args.output_dir.resolve()
    need(args.output_dir.is_absolute() and not output.is_relative_to(ROOT), "OUTPUT_LOCATION")
    need(not output.exists(), "OUTPUT_EXISTS")
    before = inventory(args.manifest_sha256)
    output.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location("retained_source_dependency", ROOT/"retained_dependency.py")
    need(spec is not None and spec.loader is not None, "DEPENDENCY_HELPER_LOADER")
    dependency = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dependency)
    retained, dependency_record = dependency.extract_authenticated(
        ROOT/"inputs/GEOMETRIC_SOURCE_DEPENDENCY.zip", output/"retained-source")
    retained_before = dependency.extracted_inventory(retained)
    rows = []
    fresh_numeric = {}
    jobs = [
        ("actual_cone", "compute/certify_exterior_response.py",
         "compute/ACTUAL_CONE_RESPONSE_NORMAL.json", "compute/ACTUAL_CONE_RESPONSE_OPTIMIZED.json"),
        ("independent_numeric_review", "numeric-review/check_cone_response.py",
         "numeric-review/REVIEW_FINAL_NORMAL.json", "numeric-review/REVIEW_FINAL_OPTIMIZED.json"),
        ("numeric_rejections", "numeric-review/check_rejection_controls.py",
         "numeric-review/CONTROLS_FINAL_NORMAL.json", "numeric-review/CONTROLS_FINAL_OPTIMIZED.json"),
        ("manuscript_display", "manuscript/build_causal_section.py",
         "manuscript/MANUSCRIPT_ARITHMETIC.json", "manuscript/MANUSCRIPT_ARITHMETIC.json"),
        ("exact_bounds", "check_bounds.py", "golden/BOUNDS.json", "golden/BOUNDS.json"),
        ("bounds_rejections", "check_controls.py", "golden/CONTROLS.json", "golden/CONTROLS.json"),
    ]
    for label, program, normal_golden, optimized_golden in jobs:
        collected = []
        for optimized in (False, True):
            name = label+("_optimized.json" if optimized else "_normal.json")
            child_output = output/name
            extra = []
            output_args = ["--output", str(child_output)]
            auxiliary = []
            if label == "actual_cone":
                extra = ["--source-root", str(retained), "--n", "4", "--require-positive"]
            elif label in ("independent_numeric_review", "numeric_rejections"):
                need(optimized in fresh_numeric, "FRESH_NUMERIC_REPORT_REQUIRED")
                extra = ["--source-root", str(retained), "--report", str(fresh_numeric[optimized]),
                         "--producer-script", str(ROOT/"compute/certify_exterior_response.py")]
                if label == "independent_numeric_review":
                    extra += ["--oracle-spots"]
            elif label == "exact_bounds":
                need(optimized in fresh_numeric, "FRESH_NUMERIC_REPORT_REQUIRED")
                extra = ["--numeric-report", str(fresh_numeric[optimized])]
            elif label == "manuscript_display":
                need(optimized in fresh_numeric, "FRESH_NUMERIC_REPORT_REQUIRED")
                child_folder = output/name.removesuffix(".json")
                need(not child_folder.exists(), "DISPLAY_OUTPUT_EXISTS")
                child_output = child_folder/"MANUSCRIPT_ARITHMETIC.json"
                extra = ["--report", str(fresh_numeric[optimized]), "--output-dir", str(child_folder)]
                output_args = []
            argv = [sys.executable]+(["-O"] if optimized else [])+[
                "-I", "-S", "-B", str(ROOT/program)]+extra+output_args
            completed = subprocess.run(argv, cwd=output, text=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, timeout=180, check=False)
            need(completed.returncode == 0, "CHILD_FAILED:"+name+":"+completed.stderr[-2000:])
            raw = child_output.read_bytes()
            golden = optimized_golden if optimized else normal_golden
            saved = (ROOT/golden).read_bytes()
            need(raw == saved and decode(raw) == decode(saved), "GOLDEN_MISMATCH:"+name)
            if label == "manuscript_display":
                need({p.name for p in child_folder.iterdir()} == {"CAUSAL_RESULT.tex", "MANUSCRIPT_ARITHMETIC.json"}
                     and all(p.is_file() and not p.is_symlink() for p in child_folder.iterdir()),
                     "DISPLAY_OUTPUT_INVENTORY")
                tex = (child_folder/"CAUSAL_RESULT.tex").read_bytes()
                need(tex == (ROOT/"manuscript/CAUSAL_RESULT.tex").read_bytes(), "DISPLAY_TEX_RAW_GOLDEN")
                auxiliary.append({"output": (child_folder/"CAUSAL_RESULT.tex").relative_to(output).as_posix(),
                                  "bytes": len(tex), "sha256": sha(tex), "raw_golden_equal": True})
            if label == "actual_cone":
                fresh_numeric[optimized] = child_output
            collected.append(raw)
            rows.append({"role": label, "program": program, "optimized": optimized, "exit_code": completed.returncode,
                         "output": child_output.relative_to(output).as_posix(), "bytes": len(raw), "sha256": sha(raw),
                         "additional_raw_outputs": auxiliary,
                         "typed_and_raw_golden_equal": True})
        need(collected[0] == collected[1], "EXECUTION_MODE_MISMATCH:"+program)
    retained_after = dependency.extracted_inventory(retained)
    need(retained_after == retained_before, "RETAINED_DEPENDENCY_CHANGED_DURING_REPLAY")
    after = inventory(args.manifest_sha256)
    need(after == before, "SOURCE_CHANGED_DURING_REPLAY")
    need(len(rows) == 12, "COMBINED_CHILD_COUNT")
    report = {"status": "PASS_AUTHENTICATED_FORCED_RESPONSE_REPRODUCTION",
              "manifest_sha256": args.manifest_sha256,
              "files_including_manifest": before["files_including_manifest"],
              "top_level_children": 12,
              "actual_source_integration_processes": 2,
              "independent_numerical_review_processes": 2,
              "exact_bound_arithmetic_processes": 2,
              "manuscript_display_arithmetic_processes": 2,
              "numerical_control_processes": 2, "bound_control_processes": 2,
              "targeted_in_process_controls_per_mode": {"numeric": 9, "bounds": 5, "total": 14},
              "retained_dependency": dependency_record,
              "retained_dependency_unchanged": True,
              "children": rows, "source_unchanged": True,
              "upstream_scientific_producers_replayed": False,
              "scope": "Actual causal-cone integration, independent joins/oracle spots, exact global bounds, displayed manuscript arithmetic and targeted controls. No PDF typesetting, full-domain waveform, full-domain continuous response residual or upstream full scientific replay."}
    target = output/"REPRODUCTION_RECEIPT.json"
    with target.open("x") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": report["status"], "receipt_sha256": sha(target.read_bytes()),
                      "files_including_manifest": report["files_including_manifest"]}))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, subprocess.TimeoutExpired) as error:
        sys.exit("RESPONSE_BOUND_REPRODUCTION_FAILED:"+str(error))
