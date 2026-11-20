#!/usr/bin/env python3
"""Offline tests of actual paid-invocation guards; no model/network calls."""
import copy
import datetime as dt
import unittest

from runtime_gates import gate

NOW = dt.datetime(2026, 10, 2, 12, 0, tzinfo=dt.timezone.utc)
SHA = "a" * 40


class RuntimeGuardTests(unittest.TestCase):
    def setUp(self):
        self.config = dict(
            repository="maldonado-research/HDblast", default_branch="main", ref_name="main",
            event_name="workflow_dispatch", mode="research", enabled="true", smoke_enabled="true",
            model="gpt-6.1-sol", effort="ultra", smoke_passed="true", run_id="102", run_attempt="1",
            start_utc="2026-10-02T00:00:00Z", until_utc="2026-10-09T00:00:00Z",
            daily_limit="4", total_limit="12", task_id="common-reference-construction", source_sha=SHA,
        )
        self.runs = [dict(id=102, run_attempt=1, created_at="2026-10-02T11:59:00Z",
                          display_title="HDBLAST research round", event="workflow_dispatch")]
        self.pulls = []

    def call(self, **kwargs):
        return gate(self.config, self.runs, self.pulls, now=NOW,
                    history_complete=kwargs.get("history_complete", True),
                    pulls_complete=kwargs.get("pulls_complete", True))

    def test_research_default_branch_first_attempt_allowed(self):
        result = self.call()
        self.assertTrue(result["allowed"])
        self.assertEqual(result["daily_invocations"], 1)

    def test_paid_rerun_denied_even_if_api_created_time_unchanged(self):
        self.config["run_attempt"] = "2"
        self.assertFalse(self.call()["allowed"])

    def test_api_rerun_disagreement_denied(self):
        self.runs[0]["run_attempt"] = 2
        self.assertFalse(self.call()["allowed"])

    def test_fork_other_branch_pr_and_unsupported_event_denied(self):
        for key, value in (("repository", "fork/HDblast"), ("ref_name", "research/pr"),
                           ("event_name", "pull_request_target")):
            original = self.config[key]
            self.config[key] = value
            self.assertFalse(self.call()["allowed"], (key, value))
            self.config[key] = original

    def test_disabled_missing_smoke_and_model_fallback_denied(self):
        for key, value in (("enabled", "false"), ("smoke_passed", "false"),
                           ("model", "gpt-6-sol"), ("effort", "high")):
            original = self.config[key]
            self.config[key] = value
            self.assertFalse(self.call()["allowed"], (key, value))
            self.config[key] = original

    def test_exact_expiry_is_denied(self):
        self.config["until_utc"] = NOW.isoformat()
        self.assertFalse(self.call()["allowed"])

    def test_future_campaign_is_denied(self):
        self.config["start_utc"] = "2026-10-03T00:00:00Z"
        self.assertFalse(self.call()["allowed"])

    def test_invalid_campaign_window_and_timezone_fail_closed(self):
        for value in ("2026-10-12T00:00:00Z", "2026-10-09T00:00:00", "2026-10-09T00:00:00+02:00"):
            self.config["until_utc"] = value
            with self.assertRaises(ValueError):
                self.call()

    def test_daily_budget_counts_failed_skipped_and_prior_attempts(self):
        for status, ident in (("failure", 99), ("cancelled", 100), ("skipped", 101)):
            self.runs.append(dict(self.runs[0], id=ident, conclusion=status, display_title="HDBLAST smoke round"))
        self.assertTrue(self.call()["allowed"])
        self.runs[1]["run_attempt"] = 2
        self.assertFalse(self.call()["allowed"])

    def test_campaign_total_budget_counts_prior_days(self):
        self.config["start_utc"] = "2026-10-01T00:00:00Z"
        self.config["until_utc"] = "2026-10-08T00:00:00Z"
        self.config["total_limit"] = "2"
        self.runs += [dict(self.runs[0], id=i, created_at="2026-10-01T11:00:00Z") for i in (100, 101)]
        self.assertFalse(self.call()["allowed"])

    def test_failed_research_invocation_stops_campaign_until_diagnosed_renewal(self):
        for outcome in ('failure', 'cancelled', 'timed_out', 'action_required', 'startup_failure'):
            self.runs = [self.runs[0], dict(self.runs[0], id=101, created_at="2026-10-02T01:00:00Z",
                                           conclusion=outcome)]
            self.assertFalse(self.call()["allowed"], outcome)
        self.config['start_utc'] = "2026-10-02T11:00:00Z"
        self.assertTrue(self.call()["allowed"])

    def test_prior_rerun_cannot_hide_original_failure_from_fail_stop(self):
        self.runs.append(dict(self.runs[0], id=101, run_attempt=2, conclusion='success'))
        self.assertFalse(self.call()["allowed"])

    def test_above_hard_budget_and_boolean_budget_fail_closed(self):
        for key, value in (("daily_limit", "5"), ("total_limit", "13"), ("daily_limit", True)):
            original = self.config[key]
            self.config[key] = value
            with self.assertRaises(ValueError):
                self.call()
            self.config[key] = original

    def test_incomplete_metadata_and_missing_current_run_fail_closed(self):
        for field in ("history_complete", "pulls_complete"):
            with self.assertRaises(ValueError):
                self.call(**{field: False})
        self.runs = []
        with self.assertRaises(ValueError):
            self.call()

    def test_duplicate_run_history_fails_closed(self):
        self.runs.append(copy.deepcopy(self.runs[0]))
        with self.assertRaises(ValueError):
            self.call()

    def test_stale_cross_day_current_run_denied(self):
        self.config["start_utc"] = "2026-10-01T00:00:00Z"
        self.config["until_utc"] = "2026-10-08T00:00:00Z"
        self.runs[0]["created_at"] = "2026-10-01T23:59:00Z"
        self.assertFalse(self.call()["allowed"])

    def test_pending_automation_pr_pauses_research(self):
        self.pulls = [dict(state="open", head={"ref": "research/auto-101-1"}, body="")]
        self.assertFalse(self.call()["allowed"])

    def test_closed_unmerged_task_source_pr_is_not_repeated(self):
        self.pulls = [dict(state="closed", merged_at=None, head={"ref": "research/auto-101-1"},
                           body=f"Research evidence\n\nHDBLAST-AUTO-ROUND task=common-reference-construction source={SHA}\n")]
        self.assertFalse(self.call()["allowed"])
        self.pulls[0]["body"] = self.pulls[0]["body"].replace(SHA, "b" * 40)
        self.assertTrue(self.call()["allowed"])

    def test_smoke_has_separate_enable_and_shared_budgets(self):
        self.config.update(mode="smoke", enabled="false", smoke_passed="false")
        self.runs[0]["display_title"] = "HDBLAST smoke round"
        self.assertTrue(self.call()["allowed"])
        self.config["smoke_enabled"] = "false"
        self.assertFalse(self.call()["allowed"])

    def test_second_smoke_in_campaign_is_denied(self):
        self.config["mode"] = "smoke"
        self.runs[0]["display_title"] = "HDBLAST smoke round"
        self.runs.append(dict(self.runs[0], id=101, conclusion="failure"))
        self.assertFalse(self.call()["allowed"])

    def test_smoke_unknown_dispatch_history_is_conservative(self):
        self.config["mode"] = "smoke"
        self.runs[0]["display_title"] = "HDBLAST smoke round"
        self.runs.append(dict(self.runs[0], id=101, display_title="old ambiguous workflow name"))
        self.assertFalse(self.call()["allowed"])

    def test_schedule_research_allowed_but_schedule_smoke_denied(self):
        self.config["event_name"] = "schedule"
        self.assertTrue(self.call()["allowed"])
        self.config["mode"] = "smoke"
        self.assertFalse(self.call()["allowed"])


if __name__ == "__main__":
    unittest.main()
