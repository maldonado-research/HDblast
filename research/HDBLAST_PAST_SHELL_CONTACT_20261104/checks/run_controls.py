"""Exercise six input/CLI rejections and two source-authentication controls.

These are transport/authentication controls, not mathematical proof mutants.
Scientific producer programs are not executed by the rejected cases.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def invoke(capsule, pin, output=None, optimized=False):
    command = [sys.executable, "-I", "-S", "-B"]
    if optimized:
        command.append("-O")
    command.extend([str(capsule / "checks/run_checkpoint.py"),
                    "--manifest-sha256", pin])
    if output is None:
        command.append("--authenticate-only")
    else:
        command.extend(["--output-dir", str(output)])
    result = subprocess.run(command, capture_output=True, timeout=30, check=False)
    need(result.stderr == b"" and len(result.stdout.splitlines()) == 1,
         "CONTROL_SINGLE_CLEAN_SUMMARY")
    summary = json.loads(result.stdout)
    need(isinstance(summary, dict), "CONTROL_SUMMARY_OBJECT")
    return result.returncode, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--manifest-sha256")
    args = parser.parse_args()
    out = args.output_dir.absolute()
    need(all(not path.is_symlink() for path in (out, *out.parents)), "OUTPUT_SYMLINK")
    out = out.resolve()
    need(not out.exists() and not out.is_relative_to(ROOT), "FRESH_EXTERNAL_OUTPUT")
    manifest_bytes = (ROOT / "MANIFEST.json").read_bytes()
    pin = hashlib.sha256(manifest_bytes).hexdigest()
    if args.manifest_sha256 is not None:
        need(pin == args.manifest_sha256, "MANIFEST_PIN")
    # A positive source-authentication run validates the baseline before copies
    # or any deliberate invalid case are created. It does not run producers.
    positives = []
    for optimized in (False, True):
        code, summary = invoke(ROOT, pin, optimized=optimized)
        need(code == 0 and summary.get("status") == "PASS_CLOSED_PAST_SHELL_CONTACT_SOURCE"
             and summary.get("manifest_sha256") == pin, "BASELINE_AUTHENTICATION")
        positives.append({"optimized": optimized, "exit_code": code, "summary": summary})
    out.mkdir(parents=True)
    rejections = []
    for optimized in (False, True):
        code, summary = invoke(ROOT, "0" * 64, optimized=optimized)
        need(code == 2 and summary == {"status": "FAIL", "witness": "MANIFEST_PIN"},
             "WRONG_PIN_REJECTION")
        rejections.append({"case": "wrong_manifest_pin_" + ("optimized" if optimized else "normal"),
                           "exit_code": code, "summary": summary})
    cases = [
        ("changed_selected_input", "primary/inputs/SHELL_BOX.json"),
        ("changed_reference_report", "primary/references/NORMAL.json"),
        ("changed_producer", "primary/verify.py"),
    ]
    for label, member in cases:
        fixture = out / label
        shutil.copytree(ROOT, fixture, symlinks=False)
        path = fixture / member
        path.write_bytes(path.read_bytes() + b"\n")
        code, summary = invoke(fixture, pin)
        need(code == 2 and summary == {"status": "FAIL", "witness": "MEMBER_BYTES:" + member},
             "MEMBER_TAMPER_REJECTION:" + label)
        rejections.append({"case": label, "exit_code": code, "summary": summary})
    occupied = out / "occupied_output"
    occupied.mkdir()
    marker = occupied / "preserve.txt"
    marker.write_bytes(b"existing-output-must-remain-unchanged\n")
    code, summary = invoke(ROOT, pin, output=occupied)
    need(code == 2 and summary == {"status": "FAIL", "witness": "FRESH_OUTPUT_DIRECTORY"}
         and marker.read_bytes() == b"existing-output-must-remain-unchanged\n"
         and {path.name for path in occupied.iterdir()} == {"preserve.txt"},
         "OCCUPIED_OUTPUT_REJECTION_AND_PRESERVATION")
    rejections.append({"case": "occupied_output", "exit_code": code, "summary": summary})
    need(len(positives) == 2 and len(rejections) == 6, "CONTROL_ACCOUNTING")
    need((ROOT / "MANIFEST.json").read_bytes() == manifest_bytes, "SOURCE_MANIFEST_UNCHANGED")
    # Final baseline authentication also occurs inside every actual producer
    # replay. Here the original bytes are read without adding another child.
    original = json.loads(manifest_bytes)["files"]
    for name, row in original.items():
        data = (ROOT / name).read_bytes()
        need(len(data) == row["bytes"] and hashlib.sha256(data).hexdigest() == row["sha256"],
             "ORIGINAL_CAPSULE_UNCHANGED")
    report = {
        "schema": "HDBLAST_PAST_SHELL_CONTACT_INPUT_CONTROLS_V1",
        "status": "PASS_PAST_SHELL_CONTACT_INPUT_CONTROLS",
        "manifest_sha256": pin, "actual_rejected_children": len(rejections),
        "positive_authentication_only_children": len(positives),
        "scientific_producer_children": 0,
        "positive_controls": positives, "rejected_controls": rejections,
        "scope": "Deliberate source/reference/producer transport corruption, incorrect manifest pins and occupied output are rejected. These controls are not independent mathematical proofs, field/PDE queries or proof-mutant evaluations.",
    }
    with (out / "CONTROLS.json").open("xb") as stream:
        stream.write((json.dumps(report, sort_keys=True, indent=2) + "\n").encode())
    print(json.dumps({"status": report["status"], "actual_rejected_children": len(rejections),
                      "positive_authentication_only_children": len(positives),
                      "manifest_sha256": pin}, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, subprocess.SubprocessError) as error:
        print(json.dumps({"status": "FAIL", "witness": str(error)}, sort_keys=True))
        sys.exit(2)
