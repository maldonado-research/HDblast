#!/usr/bin/env python3
"""Trusted round selection and publication gate; no model or remote calls."""
import argparse
import hashlib
import json
import math
import re
import shutil
import stat
from pathlib import Path

TERMINAL = {"PASS", "FAIL", "BLOCKED", "INCONCLUSIVE"}
MODEL = "gpt-6.1-sol"
EFFORT = "ultra"
MAX_EVALUATIONS = 3
MAX_BYTES = 20 * 1024 * 1024
SUFFIXES = {".md", ".json", ".py", ".sh", ".csv", ".npy", ".npz", ".png", ".pdf", ".txt", ".sha256"}


def read(path):
    def finite_float(value):
        result = float(value)
        if not math.isfinite(result): raise ValueError("Nonfinite JSON number")
        return result
    return json.loads(Path(path).read_text(), parse_float=finite_float,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))


def write(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def select(queue, state, source_sha, run_id, model, effort, scheduler_sha=None):
    if model != MODEL or effort != EFFORT:
        raise ValueError("Require the catalog-backed gpt-6.1-sol / ultra configuration; no fallback")
    if not re.fullmatch(r"[0-9a-f]{40}", source_sha):
        raise ValueError("Require full Git commit SHA")
    scheduler_sha = scheduler_sha or source_sha
    if not re.fullmatch(r"[0-9a-f]{40}", scheduler_sha):
        raise ValueError("Require full scheduler Git commit SHA")
    if not re.fullmatch(r"[0-9]+-[0-9]+", run_id):
        raise ValueError("Run ID must be GITHUB_RUN_ID-GITHUB_RUN_ATTEMPT")
    if queue["schema_version"] != 1 or state["schema_version"] != 1:
        raise ValueError("Unsupported queue/state schema")
    records = state["records"]
    if any(record["run_id"] == run_id for record in records):
        return {"ready": False, "reason": "This run ID already has a terminal record"}
    completed = {record["task_id"]: record for record in records}
    for item in queue["items"]:
        if item["id"] in completed or item.get("status", "pending") != "pending":
            continue
        if any(completed.get(dep, {}).get("status") != "PASS" for dep in item["depends_on"]):
            continue
        return {
            "ready": True, "schema_version": 1, "source_sha": source_sha,
            "scheduler_sha": scheduler_sha,
            "run_id": run_id, "task": item, "model": model, "reasoning_effort": effort,
            "round_path": "research/HDBLAST_ROUND_" + run_id,
            "max_evaluations": MAX_EVALUATIONS,
            "frozen_sources": queue["frozen_sources"],
        }
    return {"ready": False, "reason": "No uncompleted task with passing prerequisites; preserve the last outcome"}


def safe_files(root):
    root = Path(root)
    if root.is_symlink() or not root.is_dir():
        raise ValueError("Artifact root must be a regular directory")
    total = 0
    files = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("Symlink in round artifact")
        if path.is_dir():
            continue
        if not stat.S_ISREG(path.lstat().st_mode):
            raise ValueError("Nonregular artifact file")
        relative = path.relative_to(root).as_posix()
        if any(part.startswith(".") for part in Path(relative).parts) or path.suffix not in SUFFIXES:
            raise ValueError("Unsupported artifact path: " + relative)
        size = path.stat().st_size
        total += size
        if size > MAX_BYTES or total > MAX_BYTES:
            raise ValueError("Round artifact exceeds 20 MiB publication bound")
        files[relative] = digest(path)
    return files


def validate(context, state, root):
    if not context.get("ready"):
        raise ValueError("Cannot publish an unselected round")
    files = safe_files(root)
    result = read(Path(root) / "ROUND_RESULT.json")
    for key, expected in {
        "schema_version": 1, "run_id": context["run_id"], "task_id": context["task"]["id"],
        "source_sha": context["source_sha"], "model": MODEL, "reasoning_effort": EFFORT,
        "scheduler_sha": context["scheduler_sha"],
    }.items():
        if result.get(key) != expected:
            raise ValueError("Incorrect result provenance: " + key)
    if result.get("status") not in TERMINAL or result.get("historical_failures_preserved") is not True:
        raise ValueError("Missing honest terminal scientific outcome")
    for key in ("summary", "next_step"):
        if not isinstance(result.get(key), str) or not result[key].strip():
            raise ValueError("Missing result field: " + key)
    for key in ("claims", "limitations"):
        if not isinstance(result.get(key), list) or not result[key] or any(not isinstance(item, str) or not item.strip() for item in result[key]):
            raise ValueError("Require a nonempty array of strings: " + key)
    for name in ("REGISTRATION.md", "RESULTS.md", "INDEPENDENT_REVIEW.md"):
        if name not in files or (Path(root) / name).stat().st_size < 40:
            raise ValueError("Missing substantive round evidence: " + name)
    artifacts = result.get("artifacts", {})
    if artifacts != {name: sha for name, sha in files.items() if name != "ROUND_RESULT.json"}:
        raise ValueError("Artifact digest manifest is incomplete or incorrect")
    evaluations = result.get("evaluations", [])
    if not isinstance(evaluations, list) or not 1 <= len(evaluations) <= MAX_EVALUATIONS:
        raise ValueError("Require 1 to 3 documented meaningful evaluations")
    keys = set()
    for evaluation in evaluations:
        if not isinstance(evaluation, dict) or not isinstance(evaluation.get("purpose"), str):
            raise ValueError("Require evaluation object with purpose")
        if not isinstance(evaluation.get("command"), list) or not evaluation["command"] or any(not isinstance(item, str) or not item for item in evaluation["command"]):
            raise ValueError("Require evaluation command as a nonempty argv array")
        key = (evaluation.get("purpose"), tuple(evaluation.get("command", [])), evaluation.get("inputs_sha256"))
        if not key[0] or not key[1] or not re.fullmatch(r"[0-9a-f]{64}", str(key[2])) or key in keys:
            raise ValueError("Missing, duplicate, or unverifiable evaluation")
        keys.add(key)
        if type(evaluation.get("exit_code")) is not int or evaluation.get("outcome") not in TERMINAL:
            raise ValueError("Missing actual evaluation exit/outcome")
        if not isinstance(evaluation.get("evidence_paths"), list) or not evaluation["evidence_paths"] or any(not isinstance(path, str) or path not in artifacts for path in evaluation["evidence_paths"]):
            raise ValueError("Evaluation does not point to retained evidence")
    new_state = json.loads(json.dumps(state))
    existing = [row for row in state["records"] if row["run_id"] == context["run_id"] or row["task_id"] == context["task"]["id"]]
    record = {"run_id": context["run_id"], "task_id": context["task"]["id"],
              "source_sha": context["source_sha"], "status": result["status"],
              "scheduler_sha": context["scheduler_sha"],
              "result_sha256": files["ROUND_RESULT.json"], "round_path": context["round_path"]}
    if existing:
        if len(existing) == 1 and existing[0] == record:
            return result, new_state, False
        raise ValueError("Refuse to overwrite an earlier task/run outcome")
    new_state["records"].append(record)
    return result, new_state, True


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    start = sub.add_parser("select")
    for name in ("queue", "state", "source-sha", "run-id", "model", "effort", "out"):
        start.add_argument("--" + name, required=True)
    start.add_argument("--scheduler-sha")
    finish = sub.add_parser("publish")
    for name in ("context", "state", "round-root", "repo-root", "state-output", "metadata-out"):
        finish.add_argument("--" + name, required=True)
    args = parser.parse_args()
    if args.command == "select":
        context = select(read(args.queue), read(args.state), args.source_sha, args.run_id, args.model, args.effort, args.scheduler_sha)
        write(args.out, context)
        print(json.dumps({"ready": context["ready"], "reason": context.get("reason"), "round_path": context.get("round_path")}))
        return
    context, state = read(args.context), read(args.state)
    result, new_state, changed = validate(context, state, args.round_root)
    repo = Path(args.repo_root).resolve()
    destination = repo / context["round_path"]
    if destination.parent != repo / "research" or not re.fullmatch(r"research/HDBLAST_ROUND_[0-9]+-[0-9]+", context["round_path"]):
        raise ValueError("Unsafe destination")
    if changed:
        if destination.exists():
            raise ValueError("Refuse to overwrite existing checkpoint")
        shutil.copytree(args.round_root, destination)
        write(args.state_output, new_state)
    write(args.metadata_out, {"changed": changed, "title": "HDBLAST round: " + context["task"]["id"],
                              "body": result["summary"] + "\n\nOutcome: " + result["status"] +
                              "\n\nHDBLAST-AUTO-ROUND task=" + context["task"]["id"] + " source=" + context["source_sha"] +
                              "\n\nsource_sha: " + context["source_sha"] +
                              "\n\nLimits: " + str(result["limitations"]) +
                              "\n\nEvidence and independent review are retained in " + context["round_path"] +
                              ". This PR does not merge itself or publish a DOI."})


if __name__ == "__main__":
    main()
