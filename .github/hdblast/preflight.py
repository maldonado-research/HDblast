#!/usr/bin/env python3
"""Read-only remote metadata adapter for the offline invocation gate."""
import argparse
import json
import os
import re
import subprocess
from pathlib import Path
from urllib.parse import urlencode

from round_control import digest, read, select, write
from runtime_gates import gate


def api(path):
    return json.loads(subprocess.check_output(["gh", "api", path], text=True))


def metadata(repo, campaign_start):
    runs = []
    for page in range(1, 101):
        query = urlencode({"per_page": 100, "page": page, "created": ">=" + campaign_start})
        response = api(f"repos/{repo}/actions/workflows/hdblast_research_round.yml/runs?{query}")
        rows = response["workflow_runs"]
        runs.extend(rows)
        if len(runs) >= response["total_count"]:
            break
        if not rows:
            raise ValueError("Incomplete workflow pagination")
    else:
        raise ValueError("Workflow pagination bound exceeded")
    pulls = []
    for page in range(1, 101):
        rows = api(f"repos/{repo}/pulls?state=all&per_page=100&page={page}")
        pulls.extend(rows)
        if len(rows) < 100:
            break
    else:
        raise ValueError("PR pagination bound exceeded")
    return runs, pulls


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "research"], required=True)
    parser.add_argument("--controls", required=True)
    parser.add_argument("--science")
    parser.add_argument("--context-out", required=True)
    args = parser.parse_args()
    mapping = {
        "enabled": "HDBLAST_RESEARCH_ENABLED", "smoke_enabled": "HDBLAST_SMOKE_ENABLED",
        "smoke_passed": "HDBLAST_MODEL_SMOKE_PASSED", "model": "HDBLAST_CODEX_MODEL",
        "effort": "HDBLAST_CODEX_REASONING_EFFORT", "start_utc": "HDBLAST_CAMPAIGN_START_UTC",
        "until_utc": "HDBLAST_RESEARCH_UNTIL_UTC", "daily_limit": "HDBLAST_MAX_PAID_ROUNDS_PER_DAY",
        "total_limit": "HDBLAST_MAX_CAMPAIGN_INVOCATIONS", "repository": "GITHUB_REPOSITORY",
        "ref_name": "GITHUB_REF_NAME", "default_branch": "HDBLAST_DEFAULT_BRANCH",
        "event_name": "GITHUB_EVENT_NAME", "run_id": "GITHUB_RUN_ID", "run_attempt": "GITHUB_RUN_ATTEMPT",
    }
    config = {key: os.environ.get(name, "") for key, name in mapping.items()}
    config["mode"] = args.mode
    if not os.environ.get("HDBLAST_API_KEY_CHECK"):
        raise ValueError("Secure API binding is required; its value is never logged")
    context = None
    run_id = config["run_id"] + "-" + config["run_attempt"]
    if args.mode == "research":
        source = os.environ.get("HDBLAST_SCIENCE_SHA", "")
        if not re.fullmatch("[0-9a-f]{40}", source) or not args.science:
            raise ValueError("Require exact science SHA and an existing science checkout")
        actual = subprocess.check_output(["git", "-C", args.science, "rev-parse", "HEAD"], text=True).strip()
        if source != actual:
            raise ValueError("Science checkout does not match requested full SHA")
        queue = read(Path(args.controls) / "queue.json")
        for relative, expected in queue.get("required_sources", {}).items():
            path = Path(args.science) / relative
            if not path.is_file() or path.is_symlink() or digest(path) != expected:
                raise ValueError("Reviewed junction/static bridge source missing or changed: " + relative)
        context = select(queue, read(Path(args.controls) / "round_state.json"),
                         source, run_id, config["model"], config["effort"], os.environ["GITHUB_SHA"])
        if not context["ready"]:
            with open(os.environ["GITHUB_OUTPUT"], "a") as output: output.write("ready=false\n")
            print(json.dumps(context))
            return
        for frozen in context["frozen_sources"]:
            if not (Path(args.science) / frozen["path"]).is_dir():
                raise ValueError("Required frozen science source is absent")
        config.update(task_id=context["task"]["id"], source_sha=source)
    runs, pulls = metadata(config["repository"], config["start_utc"])
    decision = gate(config, runs, pulls, history_complete=True, pulls_complete=True)
    print(json.dumps(decision, sort_keys=True))
    with open(os.environ["GITHUB_OUTPUT"], "a") as output:
        output.write("ready=" + str(decision["allowed"]).lower() + "\nrun_id=" + run_id + "\n")
    if decision["allowed"] and context: write(args.context_out, context)


if __name__ == "__main__":
    main()
