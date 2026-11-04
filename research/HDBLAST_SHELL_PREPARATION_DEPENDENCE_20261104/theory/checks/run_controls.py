#!/usr/bin/env python3
"""Portable exact replay and strict pinned-input rejection controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def need(ok, message):
    if not ok:
        raise ValueError(message)


def line(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(args.output.is_absolute() and str(args.output) == str(args.output.resolve()), "normalized absolute control output required")
    need(not args.output.exists() and not args.output.is_symlink(), "fresh control output required")
    component = Path(__file__).resolve().parents[1]
    local = json.loads((component / "LOCAL_PINS.json").read_text())
    for relative, expected in local.items():
        payload = (component / relative).read_bytes()
        need(hashlib.sha256(payload).hexdigest() == expected, "local pin: " + relative)
    source = json.loads((component / "SOURCE_PINS.json").read_text())
    reference = (component / "references/NORMAL.json").read_bytes()
    reference_sha = hashlib.sha256(reference).hexdigest()
    need(reference == (component / "references/OPTIMIZED.json").read_bytes(), "frozen mode match")
    expected_success = line({"status": "PASS_COARSE_ANALYTIC_PREPARATION_SIGN_CERTIFICATE", "output_sha256": reference_sha})
    rows = []
    with tempfile.TemporaryDirectory(prefix="preparation-coarse-controls-") as temporary:
        scratch = Path(temporary).resolve()
        for optimized in (False, True):
            output = scratch / ("optimized.json" if optimized else "normal.json")
            flags = ["-I", "-S", "-B"] + (["-O"] if optimized else [])
            command = [sys.executable, *flags, str(component / "produce.py"),
                       "--repo-root", str(args.repo_root), "--output", str(output)]
            completed = subprocess.run(command, capture_output=True, text=True)
            need(completed.returncode == 0 and completed.stdout == expected_success and completed.stderr == "",
                 "strict positive process contract")
            need(output.read_bytes() == reference, "frozen scientific reference mismatch")
            rows.append({"case_id": "optimized_replay" if optimized else "normal_replay", "pass": True,
                         "flags": flags, "returncode": completed.returncode, "stdout": completed.stdout,
                         "stderr": completed.stderr, "output_exists": True, "output_sha256": reference_sha})
        altered_root = scratch / "altered-root"
        for relative in source:
            destination = altered_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes((args.repo_root / relative).read_bytes())
        for index, relative in enumerate(source):
            destination = altered_root / relative
            original = destination.read_bytes()
            destination.write_bytes(original + b"\n")
            expected_rejection = line({"status": "REJECTED", "reason": "source hash: " + relative})
            for optimized in (False, True):
                output = scratch / ("rejected-%s-%s.json" % (index, optimized))
                flags = ["-I", "-S", "-B"] + (["-O"] if optimized else [])
                command = [sys.executable, *flags, str(component / "produce.py"),
                           "--repo-root", str(altered_root), "--output", str(output)]
                completed = subprocess.run(command, capture_output=True, text=True)
                need(completed.returncode == 2 and completed.stdout == expected_rejection and completed.stderr == "",
                     "strict source rejection process contract")
                need(not output.exists(), "source rejection must precede report creation")
                rows.append({"case_id": "source_%s_%s" % (index, "optimized" if optimized else "normal"),
                             "source": relative, "pass": True, "flags": flags,
                             "returncode": completed.returncode, "stdout": completed.stdout,
                             "stderr": completed.stderr, "output_exists": False})
            destination.write_bytes(original)
    result = {"status": "PASS_COARSE_PORTABLE_REPLAY_AND_SOURCE_CONTROLS", "controls": rows,
              "reference_sha256": reference_sha, "strict_child_process_count": 6,
              "new_oracle_calls": 0, "scientific_imports": 0}
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    with args.output.open("x") as stream:
        stream.write(payload)
    print(line({"status": result["status"], "strict_child_process_count": 6}).rstrip("\n"))


if __name__ == "__main__":
    main()
