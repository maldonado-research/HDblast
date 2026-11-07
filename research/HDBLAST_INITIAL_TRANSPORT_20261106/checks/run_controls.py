"""Seven closed-capsule rejection controls; not continuous-proof validation."""
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


def rehash_manifest(root):
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path != root / "MANIFEST.json":
            raw = path.read_bytes()
            files[path.relative_to(root).as_posix()] = {"bytes":len(raw), "sha256":sha(raw)}
    value = {"schema":"HDBLAST_INITIAL_TRANSPORT_MANIFEST_V1", "files":files}
    (root / "MANIFEST.json").write_text(json.dumps(value,sort_keys=True,indent=2)+"\n", encoding="utf-8")
    return sha((root / "MANIFEST.json").read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    script = Path(__file__).absolute()
    need(script.is_file() and all(not p.is_symlink() for p in (script,*script.parents)),
         "Regular controls script path required")
    root = script.resolve(strict=True).parent.parent
    out = args.output_dir.absolute()
    need(".." not in out.parts and not out.is_relative_to(root), "Fresh controls directory must be outside capsule")
    need(not out.exists() and not out.is_symlink(), "Controls directory exists")
    need(out.parent.is_dir() and all(not p.is_symlink() for p in out.parents),
         "Regular existing controls parent required")
    out.mkdir()

    def invoke(case_root, pin, name):
        report = out / (name+".json")
        result = subprocess.run([sys.executable,"-I","-S","-B",
            str(case_root / "checks/check_checkpoint.py"), "--manifest-sha256",pin,
            "--output",str(report)],capture_output=True,text=True,timeout=60)
        return result, report

    result, report = invoke(root,args.manifest_sha256,"positive")
    need(result.returncode == 0 and report.is_file(), "Baseline checker failed")
    baseline = json.loads(report.read_text(encoding="utf-8"))
    need(baseline["status"] == "PASS" and baseline["manifest_sha256"] == args.manifest_sha256,
         "Baseline report failed")
    cases = []

    def copy(name):
        target = out / "fixtures" / name
        shutil.copytree(root,target)
        return target

    def reject(name, fixture, pin, expected, **details):
        result, report = invoke(fixture,pin,name)
        lines = result.stderr.strip().splitlines()
        diagnostic = lines[-1] if lines else ""
        need(result.returncode != 0 and not report.exists()
             and diagnostic.startswith("ValueError: ") and expected in diagnostic,
             "Control did not reach intended rejection: "+name+"; diagnostic="+diagnostic)
        cases.append({"name":name,"exit_code":result.returncode,"report_created":False,
                      "status":"REJECTED","diagnostic":diagnostic,**details})

    source = Path("inputs/GAUGE_CAUCHY_COMPLETION_06G.md")
    changed = copy("changed_source")
    with (changed / source).open("ab") as stream:
        stream.write(b"\n")
    reject("changed_source",changed,args.manifest_sha256,"Source bytes differ: "+source.as_posix())

    extra = copy("extra_file")
    (extra / "UNLISTED.txt").write_text("unlisted\n",encoding="utf-8")
    reject("extra_file",extra,args.manifest_sha256,"Closed file roster mismatch")

    missing = copy("missing_file")
    (missing / source).unlink()
    reject("missing_file",missing,args.manifest_sha256,"Closed file roster mismatch")

    linked = copy("symlink_source")
    target = out / "symlink_input_target.json"
    target.write_bytes((linked / source).read_bytes())
    (linked / source).unlink()
    (linked / source).symlink_to(target)
    reject("symlink_source",linked,args.manifest_sha256,"Capsule symlink rejected")

    stale = copy("stale_manifest_pin")
    with (stale / "MANIFEST.json").open("ab") as stream:
        stream.write(b"\n")
    reject("stale_manifest_pin",stale,args.manifest_sha256,"Manifest SHA256 mismatch")

    def arithmetic_mutant(name, field, before, after, expected):
        mutant = copy(name)
        data_path = mutant / "DATA.json"
        data = json.loads(data_path.read_text(encoding="utf-8"))
        need(data["symbolic"][field] == before, "Mathematical control target changed")
        data["symbolic"][field] = after
        data_path.write_text(json.dumps(data,sort_keys=True,indent=2)+"\n",encoding="utf-8")
        mutant_pin = rehash_manifest(mutant)
        reject(name,mutant,mutant_pin,expected,
               rehashed_manifest_sha256=mutant_pin,
               mutation="symbolic."+field+": "+str(before)+" -> "+str(after))

    arithmetic_mutant("rehashed_wrong_root_hessian","root_second_derivative","-1/4","0",
                      "Arithmetic failed: symbolic:root_second_derivative")
    arithmetic_mutant("rehashed_omitted_velocity_factor","velocity_cross_factor",2,1,
                      "Arithmetic failed: symbolic:root_second_velocity")

    need(len(cases) == 7, "Seven rejection controls required")
    receipt = {"schema":"HDBLAST_INITIAL_TRANSPORT_CONTROLS_V1", "status":"PASS",
        "manifest_sha256":args.manifest_sha256,"baseline_producer_executions":1,
        "rejected_children":len(cases),"cases":cases,
        "scope":"Protocol and actual finite-arithmetic failure controls; no continuous Sobolev proof, PDE, evolution or clock validation."}
    with (out / "CONTROLS.json").open("x",encoding="utf-8") as stream:
        stream.write(json.dumps(receipt,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"status":"PASS","rejected_children":len(cases)},sort_keys=True))


if __name__ == "__main__":
    main()
