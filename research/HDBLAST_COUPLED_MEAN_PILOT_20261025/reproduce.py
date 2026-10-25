#!/usr/bin/env python3
"""Authenticate, replay both modes, and retain a bounded scientific outcome."""
import argparse
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def authenticate(pin):
    raw = (HERE/"MANIFEST.json").read_bytes()
    need(digest(raw) == pin, "external capsule manifest pin")
    record = json.loads(raw)
    seen = set()
    for entry in record["files"]:
        rel = Path(entry["path"])
        need(not rel.is_absolute() and ".." not in rel.parts and str(rel) not in seen,
             "manifest member path")
        seen.add(str(rel))
        path = HERE/rel
        need(path.is_file() and not path.is_symlink(), "regular source member")
        value = path.read_bytes()
        need(len(value) == entry["bytes"] and digest(value) == entry["sha256"],
             "source bytes "+str(rel))
    for path in HERE.rglob("*"):
        need(not path.is_symlink(), "no source symlink")
        if path.is_file():
            need(str(path.relative_to(HERE)) in seen|{"MANIFEST.json"},
                 "unlisted source file")
    return {"manifest_sha256":pin,"files_including_manifest":len(seen)+1}


def child(command, expected=None):
    process = subprocess.run(command, cwd=HERE, capture_output=True, check=True)
    if expected is not None:
        need(process.stdout == expected, "exact replay reference")
    return process.stdout, process.stderr


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest-sha256", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--native", action="store_true")
    p.add_argument("--repo-root", type=Path)
    p.add_argument("--require-accuracy-pass", action="store_true")
    p.add_argument("--export-native-json", action="store_true")
    a = p.parse_args()
    output = a.output_dir.resolve()
    need(not output.exists() and HERE not in output.parents and output != HERE,
         "fresh output outside source")
    need(not a.native or a.repo_root is not None, "native needs explicit repository root")
    need(not a.require_accuracy_pass or a.native, "accuracy assertion needs native")
    before = authenticate(a.manifest_sha256)
    output.mkdir(parents=True)
    node = shutil.which("node")
    need(node is not None, "Node with BigInt support required")
    roles = []
    raw, stderr = child([node,str(HERE/"independent_kernel_checks.js")],
                        (HERE/"INDEPENDENT_EXACT_REPORT.json").read_bytes())
    (output/"INDEPENDENT_EXACT_REPORT.json").write_bytes(raw)
    roles.append({"role":"independent_polygon_kernel_controls","sha256":digest(raw)})
    fixture, stderr = child([node,str(HERE/"generate_fixtures.js")],
                            (HERE/"EXACT_FIXTURES.json").read_bytes())
    roles.append({"role":"independent_exact_fixture_regeneration","sha256":digest(fixture)})
    mode_reports = []
    for mode, flags in (("normal",[]),("optimized",["-O"])):
        command = [sys.executable,"-I","-S","-B",*flags,str(HERE/"check_kernel.py")]
        raw, stderr = child(command)
        report = json.loads(raw)
        need(report["status"] == "PASS_EXACT_PYTHON_RIDGE_AND_PROFILE_CONTROLS",
             "Python representation controls")
        (output/("KERNEL_"+mode.upper()+".json")).write_bytes(raw)
        mode_reports.append(raw)
        roles.append({"role":"python_kernel_"+mode,"sha256":digest(raw)})
    need(mode_reports[0] == mode_reports[1], "normal and optimized kernel reports differ")

    native_result = None
    native_files = None
    if a.native:
        mode_payloads = []
        for mode, flags in (("normal",[]),("optimized",["-O"])):
            destination = output/("native_"+mode)
            command = [sys.executable,"-I","-S","-B",*flags,str(HERE/"native_pilot.py"),
                       "--repo-root",str(a.repo_root.resolve()),
                       "--output-dir",str(destination)]
            raw, stderr = child(command)
            (output/("NATIVE_"+mode.upper()+".stdout.txt")).write_bytes(raw)
            (output/("NATIVE_"+mode.upper()+".stderr.txt")).write_bytes(stderr)
            evaluated = destination/"EVALUATED_CANDIDATE.json"
            eval_command = [sys.executable,"-I","-S","-B",*flags,str(HERE/"evaluate_candidate.py"),
                            "--native-output",str(destination),"--output",str(evaluated)]
            eval_raw, eval_stderr = child(eval_command)
            payload = {name:(destination/name).read_bytes()
                       for name in ("REPORT.json","PROFILES.json","EVALUATED_CANDIDATE.json")}
            mode_payloads.append(payload)
            roles.append({"role":"bounded_native_"+mode,
                          "files":{name:digest(value) for name,value in payload.items()}})
        need(mode_payloads[0] == mode_payloads[1],
             "normal and optimized native/evaluated files differ")
        native_files = mode_payloads[0]
        golden = HERE/"results"
        if golden.exists():
            for name, raw in native_files.items():
                need((golden/name).read_bytes() == raw, "frozen native result "+name)
        native_report = json.loads(native_files["REPORT.json"])
        native_result = {"status":native_report["status"],
                         "accuracy_target_passed":native_report["passes_target"],
                         "selected_mesh_denominator":native_report["levels"][-1]["mesh_denominator"],
                         "total_X_error_upper":native_report["levels"][-1]["ledger"]["total_X_error_upper"],
                         "input_hashes":native_report["input_hashes"]}
        if a.require_accuracy_pass:
            need(native_report["passes_target"], "required scientific accuracy gate failed")

    after = authenticate(a.manifest_sha256)
    need(before == after, "source changed during replay")
    receipt = {"status":"PASS_REPLAY_WITH_RECORDED_SCIENTIFIC_OUTCOME",
               "source":before,"python":platform.python_version(),
               "node":subprocess.check_output([node,"--version"],text=True).strip(),
               "roles":roles,"native":native_result,
               "limits":["Portable controls verify the exact kernel, not the retained field.",
                         "Native accuracy outcome is reported separately from software replay.",
                         "No pointwise Einstein residual, full angular solution, invariant sign, nonlinear family or physical-origin inference."]}
    (output/"REPRODUCTION_RECEIPT.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt,sort_keys=True),flush=True)
    if native_files is not None and a.export_native_json:
        export = {name:raw.decode("utf-8") for name,raw in native_files.items()}
        export["REPRODUCTION_RECEIPT.json"] = (output/"REPRODUCTION_RECEIPT.json").read_text()
        print("HDB_RIDGE_EXPORT:"+json.dumps(export,sort_keys=True,separators=(",",":")),
              flush=True)


if __name__ == "__main__":
    main()
