"""Portable exact replay of a completed registered source/operator archive.

Only the standard library is imported before validating the externally pinned
whole-payload manifest. An archived actual execution review must pass before
either new registered route is launched. A scientific negative status is kept.
No root, interpreter or credential path is hardcoded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import types

ROUTES=("primary","independent")
MANIFEST_NAME="MANIFEST.json"
HELPER_NAME="execution/replay_checkpoint.py"


def require(condition,message):
    if not condition:raise ValueError(message)


def load_json(raw):
    def unique(items):
        output={}
        for key,value in items:
            require(key not in output,"Duplicate JSON key")
            output[key]=value
        return output
    def invalid(value):raise ValueError("Nonfinite JSON constant")
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)


def sha_file(path):
    digest=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda:handle.read(65536),b""):digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value):
    return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()


def digest(value):return hashlib.sha256(canonical_json(value)).hexdigest()


def exact_directory(path):
    path=Path(os.path.abspath(path))
    require(path.is_dir() and path.resolve()==path,"Real checkpoint directory required")
    for ancestor in (path,*path.parents):
        require(not ancestor.is_symlink(),"Symlink checkpoint ancestor forbidden")
    return path


def safe_relative_file(root,name):
    require(type(name) is str and name and "\\" not in name,"Canonical relative path required")
    parsed=PurePosixPath(name)
    require(not parsed.is_absolute() and parsed.as_posix()==name and
            not any(part in (".","..","") for part in name.split("/")),"Unsafe manifest path")
    path=root
    for part in parsed.parts:
        path=path/part
        require(not path.is_symlink(),"Manifest symlink forbidden")
    require(path.exists() and stat.S_ISREG(path.stat().st_mode),"Regular manifest file absent: "+name)
    return path


def verify_manifest(checkpoint,expected_sha256):
    root=exact_directory(checkpoint)
    require(re.fullmatch("[0-9a-f]{64}",expected_sha256 or ""),"Explicit expected manifest SHA256 required")
    manifest_path=safe_relative_file(root,MANIFEST_NAME)
    require(manifest_path.stat().st_size<=32*1024*1024,"Unexpected oversized manifest")
    raw=manifest_path.read_bytes()
    require(hashlib.sha256(raw).hexdigest()==expected_sha256,"Whole-payload manifest hash differs")
    manifest=load_json(raw)
    require(type(manifest) is dict and set(manifest)=={"schema_version","files"} and
            type(manifest["schema_version"]) is int and manifest["schema_version"]==1,
            "Unsupported exact manifest schema")
    files=manifest["files"]
    require(type(files) is dict and files and MANIFEST_NAME not in files,"Manifest must exclude itself")
    required={HELPER_NAME,"FULL_REGISTRATION.json","execution/registration_guard.py",
              "execution/execute_bounded.py","execution/review_completed.py"}
    for route in ROUTES:required|={f"results/{route}/OUTPUT.json",f"results/{route}/EXECUTION.json"}
    require(required<=set(files),"Completed replay payload inventory incomplete")
    actual=set()
    for path in root.rglob("*"):
        require(not path.is_symlink(),"Payload symlink forbidden")
        if path.is_dir():continue
        require(stat.S_ISREG(path.stat().st_mode),"Payload special file forbidden")
        name=path.relative_to(root).as_posix()
        if name!=MANIFEST_NAME:actual.add(name)
    require(actual==set(files),"Whole-payload membership differs")
    for name,item in files.items():
        require(type(item) is dict and set(item)=={"sha256","bytes"} and
                type(item["bytes"]) is int and item["bytes"]>=0 and
                re.fullmatch("[0-9a-f]{64}",item["sha256"] or ""),"Invalid manifest item")
        path=safe_relative_file(root,name)
        require(path.stat().st_size==item["bytes"] and sha_file(path)==item["sha256"],
                "Whole-payload member changed: "+name)
    require(sha_file(Path(__file__))==files[HELPER_NAME]["sha256"],
            "Executing replay helper differs from pinned archive helper")
    return root,manifest


def load_verified_module(root,relative,module_name):
    """Caller has already verified whole manifest and source hash before load."""
    path=safe_relative_file(root,relative)
    module=types.ModuleType(module_name)
    module.__file__=str(path);module.__package__=""
    sys.modules[module_name]=module
    exec(compile(path.read_bytes(),str(path),"exec"),module.__dict__)
    return module


def replay(checkpoint,expected_manifest_sha256,output_directory,python,
           receipt_sha256,registration_sha256,receipt=None,optimized=False):
    sys.dont_write_bytecode=True
    root,manifest=verify_manifest(checkpoint,expected_manifest_sha256)
    for pin in (receipt_sha256,registration_sha256):
        require(re.fullmatch("[0-9a-f]{64}",pin or ""),"Explicit receipt and registration SHA256 pins required")
    interpreter=Path(os.path.abspath(python))
    require(interpreter.is_file() and os.access(interpreter,os.X_OK),"Explicit executable Python path required")
    receipt_path=Path(os.path.abspath(receipt)) if receipt is not None else root/"FREEZE_RECEIPT.json"
    require(receipt_path.is_file() and not receipt_path.is_symlink(),"Pinned freeze receipt required")
    require(sha_file(receipt_path)==receipt_sha256,"Freeze receipt hash differs")
    output=Path(os.path.abspath(output_directory))
    require(output!=root and root not in output.parents and not output.exists() and not output.is_symlink(),
            "Fresh external output directory required")
    for ancestor in output.parents:
        require(not ancestor.is_symlink(),"Symlink output ancestor forbidden")
    # Import only authenticated standard-library guard/review modules. This
    # review rejects fabricated or incomplete archived execution receipts.
    guard=load_verified_module(root,"execution/registration_guard.py","registration_guard")
    guard.authenticate(root,receipt_path,receipt_sha256,registration_sha256)
    reviewer=load_verified_module(root,"execution/review_completed.py","replay_verified_completed_review")
    baseline_dirs={route:root/"results"/route for route in ROUTES}
    baseline_review,baseline_outputs=reviewer.review(root,receipt_path,receipt_sha256,
                                                    registration_sha256,baseline_dirs)
    verify_manifest(root,expected_manifest_sha256)
    output.mkdir(parents=True)
    baseline_record={"status":"PASS_ARCHIVED_ACTUAL_EXECUTION_REVIEW",
        "review":baseline_review,"manifest_sha256":expected_manifest_sha256,
        "registration_sha256":registration_sha256,"freeze_receipt_sha256":receipt_sha256,
        "payload_sha256":{route:digest(baseline_outputs[route]["payload"]) for route in ROUTES},
        "stable_frame_sha256":{route:digest(reviewer.stable_frame(baseline_outputs[route])) for route in ROUTES}}
    (output/"ARCHIVED_BASELINE_REVIEW.json").write_text(json.dumps(baseline_record,indent=2,sort_keys=True)+"\n")
    fresh_dirs={route:output/route for route in ROUTES};commands=[];codes={}
    try:
        for route in ROUTES:
            verify_manifest(root,expected_manifest_sha256)
            command=[str(interpreter),"-B"]+(["-O"] if optimized else [])+[
                str(root/"execution/execute_bounded.py"),"--root",str(root),
                "--receipt",str(receipt_path),"--receipt-sha256",receipt_sha256,
                "--registration-sha256",registration_sha256,"--python",str(interpreter),
                "--route",route,"--output-dir",str(fresh_dirs[route])]
            if optimized:command.append("--optimized")
            commands.append(command)
            environment=dict(os.environ)
            environment.pop("PYTHONPATH",None);environment.pop("PYTHONHOME",None)
            environment["PYTHONDONTWRITEBYTECODE"]="1"
            with (output/(route+"_LAUNCH.log")).open("xb") as log:
                completed=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,
                                         env=environment,check=False)
            codes[route]=completed.returncode
        require(all(code==0 for code in codes.values()),"A fresh bounded route failed; artifacts retained")
        fresh_review,fresh_outputs=reviewer.review(root,receipt_path,receipt_sha256,
                                                  registration_sha256,fresh_dirs)
        comparisons={}
        for route in ROUTES:
            old=baseline_outputs[route];new=fresh_outputs[route]
            old_frame=reviewer.stable_frame(old);new_frame=reviewer.stable_frame(new)
            comparisons[route]={"payload_exact":canonical_json(old["payload"])==canonical_json(new["payload"]),
                "stable_frame_exact":canonical_json(old_frame)==canonical_json(new_frame),
                "archived_payload_sha256":digest(old["payload"]),"fresh_payload_sha256":digest(new["payload"]),
                "archived_stable_frame_sha256":digest(old_frame),"fresh_stable_frame_sha256":digest(new_frame)}
        verify_manifest(root,expected_manifest_sha256)
        exact=all(value["payload_exact"] and value["stable_frame_exact"] for value in comparisons.values())
        same_science=baseline_review["status"]==fresh_review["status"]
        report={"schema_version":1,"status":"PASS_EXACT_FRESH_REPLAY" if exact and same_science else "FAIL_EXACT_FRESH_REPLAY",
            "checkpoint":str(root),"manifest_sha256":expected_manifest_sha256,
            "registration_sha256":registration_sha256,"freeze_receipt_sha256":receipt_sha256,
            "python":str(interpreter),"optimized":optimized,"commands":commands,"launcher_exit_codes":codes,
            "archived_scientific_status":baseline_review["status"],"fresh_scientific_status":fresh_review["status"],
            "scientific_status_exact":same_science,"comparisons":comparisons,
            "frozen_payload_bytes_unchanged":True,"archived_review":baseline_review,"fresh_review":fresh_review,
            "interpretation":"Exact reproduction preserves the archived scientific state, including an unresolved or negative result. It does not upgrade that state or establish a twelve-case pressure/contact certificate."}
        (output/"FRESH_REPLAY.json").write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+"\n")
        return report
    except BaseException as error:
        failure={"schema_version":1,"status":"FAIL_FRESH_REPLAY_EXECUTION","reason":type(error).__name__+": "+str(error),
                 "archived_scientific_status":baseline_review["status"],"launcher_exit_codes":codes,
                 "commands":commands,"artifacts_retained":True}
        (output/"FRESH_REPLAY_FAILURE.json").write_text(json.dumps(failure,indent=2,sort_keys=True)+"\n")
        raise


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--checkpoint",required=True)
    parser.add_argument("--expected-manifest-sha256",required=True)
    parser.add_argument("--python",required=True)
    parser.add_argument("--receipt")
    parser.add_argument("--receipt-sha256",required=True)
    parser.add_argument("--registration-sha256",required=True)
    parser.add_argument("--output-directory",required=True)
    parser.add_argument("--optimized",action="store_true")
    args=parser.parse_args()
    try:
        report=replay(args.checkpoint,args.expected_manifest_sha256,args.output_directory,args.python,
                      args.receipt_sha256,args.registration_sha256,args.receipt,args.optimized)
        print(json.dumps(report,sort_keys=True,allow_nan=False))
        return 0 if report["status"]=="PASS_EXACT_FRESH_REPLAY" else 1
    except (ValueError,OSError,KeyError,TypeError) as error:
        print("REPLAY_REJECTED: "+str(error),file=sys.stderr)
        return 1


if __name__=="__main__":raise SystemExit(main())
