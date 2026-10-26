#!/usr/bin/env python3
"""Replay exact positive checks in normal/-O modes and all deliberate faults."""
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
from hashlib import sha256
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def fresh_output_path(value):
    """Validate before importing the project checker or spawning work."""
    path = Path(os.path.abspath(value))
    for part in (path,) + tuple(path.parents):
        if part.is_symlink():
            raise ValueError("OUTPUT_SYMLINK")
    if not path.parent.is_dir():
        raise ValueError("OUTPUT_PARENT_REQUIRED")
    resolved = path.resolve()
    protected = [HERE.resolve()]
    protected.extend(p.resolve() for p in (HERE,) + tuple(HERE.parents)
                     if p.name == "HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026")
    if any(resolved == p or p in resolved.parents for p in protected):
        raise ValueError("OUTPUT_PROTECTED_SOURCE")
    if os.path.lexists(path):
        raise ValueError("OUTPUT_EXISTS")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = fresh_output_path(args.output)
    checker = HERE / "check_curl.py"
    spec = importlib.util.spec_from_file_location("curl_checker_controls", checker)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reports, rows = [], []
    with tempfile.TemporaryDirectory(prefix="screen-curl-controls-") as tmp:
        for optimized in (False, True):
            command = [sys.executable, "-I", "-S", "-B"] + (["-O"] if optimized else []) + [str(checker)]
            path = Path(tmp) / ("optimized.json" if optimized else "normal.json")
            proc = subprocess.run(command + ["--output", str(path)], capture_output=True, text=True)
            if proc.returncode:
                raise ValueError("POSITIVE_REPLAY_FAILED: " + proc.stderr)
            reports.append(json.loads(path.read_text()))
            for mutation, label in module.CONTROLS.items():
                proc = subprocess.run(command + ["--inject", mutation], capture_output=True, text=True)
                observed = proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else ""
                if proc.returncode == 0 or observed != "ValueError: " + label:
                    raise ValueError("UNEXPECTED_CONTROL: " + mutation + " " + proc.stderr)
                rows.append({"optimized": optimized, "mutation": mutation,
                             "expected_gate": label, "returncode": proc.returncode,
                             "observed": observed})
    if reports[0] != reports[1]:
        raise ValueError("NORMAL_OPTIMIZED_MISMATCH")
    report = {"status": "PASS_CURL_POSITIVE_AND_NEGATIVE_CONTROLS",
              "normal_optimized_equal": True, "controls": rows,
              "checker_sha256": sha256(checker.read_bytes()).hexdigest(),
              "runner_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    with output.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(report["status"], len(rows))


if __name__ == "__main__":
    main()
