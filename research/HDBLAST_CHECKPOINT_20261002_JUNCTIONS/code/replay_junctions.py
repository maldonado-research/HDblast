#!/usr/bin/env python3
"""Verify frozen inputs and replay algebra checks into a new external directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def read(path):
    return json.loads(path.read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    package = Path(__file__).resolve().parents[1]
    repo, output = args.repo.resolve(), args.output.resolve()
    if output == repo or repo in output.parents:
        raise ValueError("Use a fresh output directory outside the repository")
    pins = read(package / "evidence/SOURCE_PINS.json")
    for name, expected in pins["files"].items():
        if digest(repo / name) != expected:
            raise ValueError("Inherited source changed: " + name)
    manifest = read(package / "MANIFEST.json")
    actual_files = {p.relative_to(package).as_posix() for p in package.rglob("*")
                    if p.is_file() and p.name != "MANIFEST.json"}
    if actual_files != set(manifest["files"]):
        raise ValueError("Checkpoint file membership differs from the frozen manifest")
    for name, expected in manifest["files"].items():
        if digest(package / name) != expected:
            raise ValueError("Checkpoint payload changed: " + name)
    output.mkdir(parents=True, exist_ok=False)
    work = output / "checkpoint"
    shutil.copytree(package, work)
    logs = output / "logs"
    logs.mkdir()
    environment = os.environ.copy()
    environment.update(PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
                       OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    commands = [
        ("action", ["code/verify_action_junctions.py", "--output", "evidence/REPLAY_ACTION.json"]),
        ("action_optimized", ["-O", "code/verify_action_junctions.py", "--output", "evidence/REPLAY_ACTION_OPTIMIZED.json"]),
        ("independent_action", ["review/independent_junction_checks.py"]),
        ("static_bridge", ["static_bridge/verify_static_bridge.py", "--repo", str(repo)]),
        ("independent_static", ["review/independent_static_bridge_checks.py", "--repo", str(repo)]),
    ]
    receipts = []
    for name, command in commands:
        with (logs / (name + ".log")).open("w") as stream:
            run = subprocess.run([sys.executable, *command], cwd=work, env=environment,
                                 stdout=stream, stderr=subprocess.STDOUT)
        receipts.append({"check": name, "command": ["python", *command], "exit_code": run.returncode})
        (output / "EXECUTION.json").write_text(json.dumps(receipts, indent=2) + "\n")
        if run.returncode:
            raise RuntimeError(f"{name} failed; see {logs / (name + '.log')}")
        print(name + " PASS", flush=True)
    action = read(work / "evidence/REPLAY_ACTION.json")
    optimized = read(work / "evidence/REPLAY_ACTION_OPTIMIZED.json")
    independent = read(work / "review/INDEPENDENT_JUNCTION_CHECKS.json")
    bridge = read(work / "static_bridge/STATIC_BRIDGE_CHECKS.json")
    static_review = read(work / "review/INDEPENDENT_STATIC_BRIDGE_CHECKS.json")
    if action != optimized or action["positive_checks"] != 21 or action["negative_controls"] != 14:
        raise RuntimeError("Action check reports are incomplete or differ under optimization")
    if independent["checks_passed"] != 21 or independent["negative_controls_rejected"] != 9:
        raise RuntimeError("Independent action checks are incomplete")
    if (len(bridge["checks"]) != 33 or len(bridge["negative_controls"]) != 5
            or not all(row["passed"] for row in bridge["checks"])
            or not all(row["detected"] for row in bridge["negative_controls"])
            or not bridge["all_source_hashes_match_pin"]):
        raise RuntimeError("Static bridge checks are incomplete")
    if static_review["checks_passed"] != 16 or static_review["negative_controls_rejected"] != 2:
        raise RuntimeError("Independent static review checks are incomplete")
    receipt = {"status": "PASS", "scope": "algebra and archived-data reanalysis only",
               "physical_evolution_executed": False, "commands": receipts,
               "manifest_sha256": digest(package / "MANIFEST.json"),
               "inherited_source_count": len(pins["files"])}
    (output / "VALIDATION.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
