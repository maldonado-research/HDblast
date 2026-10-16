#!/usr/bin/env python3
"""Replay exact normal/-O certificates and targeted wrong controls."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile


def need(condition, code):
    if not condition:
        raise ValueError(code)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), "OUTPUT_EXISTS")
    checker = Path(__file__).with_name("check_evaluation_bound.py")
    controls = runpy.run_path(str(checker))["CONTROLS"]
    normal_reports = []
    records = []
    with tempfile.TemporaryDirectory(prefix="evaluation-replay-") as directory:
        base = Path(directory)
        for optimized in [False, True]:
            mode = "optimized" if optimized else "normal"
            command = [sys.executable] + (["-O"] if optimized else []) + [str(checker)]
            target = base / (mode + ".json")
            result = subprocess.run(command + ["--output", str(target)], text=True, capture_output=True)
            need(result.returncode == 0 and target.exists(), "BASELINE_FAILED")
            report = json.loads(target.read_text())
            normal_reports.append(report)
            records.append({"mode": mode, "case": "baseline", "returncode": result.returncode,
                            "status": report["status"]})
            for mutation, expected_code in sorted(controls.items()):
                target = base / (mode + "-" + mutation + ".json")
                result = subprocess.run(command + ["--inject", mutation, "--output", str(target)],
                                        text=True, capture_output=True)
                need(result.returncode != 0 and not target.exists(), "WRONG_CONTROL_ACCEPTED")
                need(result.stderr.rstrip().endswith("ValueError: " + expected_code),
                     "WRONG_CONTROL_DIAGNOSTIC")
                records.append({"mode": mode, "case": mutation, "returncode": result.returncode,
                                "diagnostic": expected_code, "no_report_written": True})
    need(normal_reports[0] == normal_reports[1], "NORMAL_OPTIMIZED_DISAGREE")
    out = {
        "status": "PASS_UNPUBLISHED_SELECTED23_EVALUATION_REPLAY",
        "normal_optimized_identical": True,
        "baseline_count": 2,
        "rejected_wrong_control_count": 2 * len(controls),
        "checker_sha256": sha256(checker.read_bytes()).hexdigest(),
        "replay_checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "records": records,
        "scope": "Angular evaluation algebra only; no source/initial/operator closure or reduced evolution certified.",
    }
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(out["status"])


if __name__ == "__main__":
    main()
