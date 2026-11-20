#!/usr/bin/env python3
"""Fail-closed, offline preflight for bounded paid HDBLAST invocations.

The caller must fetch complete canonical workflow history and PR metadata using
read-only credentials. Every workflow attempt, including failed/skipped attempts,
consumes the conservative invocation allowance; this is not a dollar/token cap.
No model output or credentials are accepted here.
"""
import datetime as dt
import re

REPOSITORY = "maldonado-research/HDblast"
MODEL = "gpt-6.1-sol"
EFFORT = "ultra"
MAX_DAILY = 4
MAX_TOTAL = 12
MAX_CAMPAIGN_DAYS = 7
ROUND_TITLES = {"smoke": "HDBLAST smoke round", "research": "HDBLAST research round"}
MARKER = re.compile(r"(?m)^HDBLAST-AUTO-ROUND task=([a-z0-9-]+) source=([0-9a-f]{40})$")
FAILED_INVOCATIONS = {"failure", "cancelled", "timed_out", "action_required", "startup_failure"}


def utc(value):
    if not isinstance(value, str):
        raise ValueError("Require an ISO-8601 UTC timestamp")
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None or parsed.utcoffset() != dt.timedelta(0):
        raise ValueError("Timestamp must have an explicit UTC offset")
    return parsed.astimezone(dt.timezone.utc)


def bounded_integer(value, upper):
    if type(value) is int:
        parsed = value
    elif isinstance(value, str) and re.fullmatch(r"[1-9][0-9]*", value):
        parsed = int(value)
    else:
        raise ValueError("Budget must be a positive integer")
    if not 1 <= parsed <= upper:
        raise ValueError("Budget exceeds the hard invocation bound")
    return parsed


def gate(config, runs, pulls, *, now=None, history_complete=False, pulls_complete=False):
    """Return {allowed, reason, campaign_invocations, daily_invocations}.

    Invalid/missing metadata raises ValueError before any paid action. Normal
    disabled, expired, exhausted, pending or duplicate-round states deny cleanly.
    Enabled values and smoke_passed are the exact strings 'true' from repository
    variables. Research additionally requires trusted selected task_id/source_sha.
    """
    now = now or dt.datetime.now(dt.timezone.utc)
    if now.tzinfo is None or now.utcoffset() != dt.timedelta(0):
        raise ValueError("Clock must use UTC")
    decision = {"allowed": False, "reason": "", "campaign_invocations": 0, "daily_invocations": 0}

    def deny(reason):
        decision["reason"] = reason
        return decision

    mode = config.get("mode")
    if mode not in ROUND_TITLES:
        raise ValueError("Unsupported paid invocation mode")
    if config.get("repository") != REPOSITORY:
        return deny("Paid invocation is restricted to the canonical repository")
    if not config.get("default_branch") or config.get("ref_name") != config["default_branch"]:
        return deny("Paid invocation is restricted to the current default branch")
    event = config.get("event_name")
    if event not in {"schedule", "workflow_dispatch"} or (mode == "smoke" and event != "workflow_dispatch"):
        return deny("Unsupported paid invocation event")
    if config.get("smoke_enabled" if mode == "smoke" else "enabled") != "true":
        return deny("This paid invocation mode is disabled")
    if config.get("model") != MODEL or config.get("effort") != EFFORT:
        return deny("Require exactly gpt-6.1-sol with ultra reasoning effort")
    if mode == "research" and config.get("smoke_passed") != "true":
        return deny("Exact-model smoke has not been marked verified")
    if str(config.get("run_attempt")) != "1":
        return deny("Paid reruns are forbidden; create a fresh bounded invocation")
    run_id = str(config.get("run_id", ""))
    if not re.fullmatch(r"[1-9][0-9]*", run_id):
        raise ValueError("Require the trusted GitHub workflow run ID")
    start, until = utc(config.get("start_utc")), utc(config.get("until_utc"))
    if until <= start or until - start > dt.timedelta(days=MAX_CAMPAIGN_DAYS):
        raise ValueError("Require a campaign window greater than zero and at most seven days")
    if not start <= now < until:
        return deny("Campaign has not started or has expired")
    daily_limit = bounded_integer(config.get("daily_limit"), MAX_DAILY)
    total_limit = bounded_integer(config.get("total_limit"), MAX_TOTAL)
    if history_complete is not True or pulls_complete is not True:
        raise ValueError("Require complete workflow and PR metadata; pagination/truncation must fail closed")
    if not isinstance(runs, list) or not isinstance(pulls, list):
        raise ValueError("Metadata must be lists")

    seen = set()
    current = None
    smoke_invocations = 0
    prior_research_failure = False
    for row in runs:
        row_id = str(row.get("id", ""))
        if not re.fullmatch(r"[1-9][0-9]*", row_id) or row_id in seen:
            raise ValueError("Missing or duplicate workflow run ID")
        seen.add(row_id)
        created = utc(row.get("created_at"))
        attempt = bounded_integer(row.get("run_attempt"), 1000)
        if created > now + dt.timedelta(minutes=5):
            raise ValueError("Workflow metadata contains a future run")
        if row_id == run_id:
            current = row
            if attempt != 1:
                return deny("Current API run metadata identifies a rerun")
            if not start <= created < until or created.date() != now.date():
                return deny("Stale queued invocation or run created outside the campaign")
            if row.get("display_title") != ROUND_TITLES[mode]:
                raise ValueError("Current run-name must identify smoke versus research")
        if not start <= created < until:
            continue
        decision["campaign_invocations"] += attempt
        if created.date() == now.date():
            decision["daily_invocations"] += attempt
        # An unknown dispatch title cannot safely establish that no smoke ran.
        if row.get("display_title") == ROUND_TITLES["smoke"] or (
            row.get("event") == "workflow_dispatch" and row.get("display_title") != ROUND_TITLES["research"]
        ):
            smoke_invocations += attempt
        if row_id != run_id and row.get("display_title") == ROUND_TITLES["research"] and (
            row.get("conclusion") in FAILED_INVOCATIONS or attempt > 1
        ):
            # GitHub exposes the latest rerun conclusion on a workflow row. A
            # rerun cannot establish that the first research attempt succeeded.
            prior_research_failure = True
    if current is None:
        raise ValueError("Current invocation absent from workflow history; retry later without a paid call")
    if decision["campaign_invocations"] > total_limit:
        return deny("Finite campaign invocation budget is exhausted")
    if decision["daily_invocations"] > daily_limit:
        return deny("UTC daily invocation budget is exhausted")
    if mode == "smoke" and smoke_invocations > 1:
        return deny("At most one paid smoke invocation is permitted per campaign")

    if mode == "research":
        if prior_research_failure:
            return deny("An earlier research invocation failed or was rerun; diagnose it before explicitly renewing the campaign")
        task = config.get("task_id", "")
        source = config.get("source_sha", "")
        if not re.fullmatch(r"[a-z0-9-]+", task) or not re.fullmatch(r"[0-9a-f]{40}", source):
            raise ValueError("Require a trusted selected task and pinned science-source SHA")
        for row in pulls:
            if row.get("state") not in {"open", "closed"} or not isinstance(row.get("head", {}).get("ref"), str):
                raise ValueError("Incomplete PR metadata")
            automatic = row["head"]["ref"].startswith("research/auto-")
            if automatic and row["state"] == "open":
                return deny("An automation PR still awaits review; no paid research call")
            markers = MARKER.findall(row.get("body") or "")
            if automatic and (task, source) in markers:
                return deny("This task/source pair already has a reviewable PR; preserve its outcome")

    decision.update(allowed=True, reason="Bounded first-attempt paid invocation authorized by explicit configuration")
    return decision
