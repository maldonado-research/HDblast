"""Offline control-plane tests. Fixtures are not scientific results."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

import round_control as rc


class RoundControlTests(unittest.TestCase):
    def setUp(self):
        self.queue = rc.read(Path(__file__).with_name("queue.json"))
        self.state = {"schema_version": 1, "records": []}
        self.sha = "e17a01b428bb8049e919c42376ab0359e41d613c"
        self.scheduler_sha = "a" * 40
        self.context = rc.select(self.queue, self.state, self.sha, "123-1", rc.MODEL, rc.EFFORT, self.scheduler_sha)
        self.temp = tempfile.TemporaryDirectory(prefix="hdblast-control-fixture-")
        self.root = Path(self.temp.name) / "round"
        self.root.mkdir()
        for name in ("REGISTRATION.md", "RESULTS.md", "INDEPENDENT_REVIEW.md"):
            (self.root / name).write_text("Offline control-plane unit-test fixture; this is not scientific evidence.\n")
        self.result = {key: self.context[key] for key in ("schema_version", "run_id", "source_sha", "scheduler_sha", "model", "reasoning_effort")}
        self.result.update({
            "task_id": self.context["task"]["id"], "status": "FAIL", "historical_failures_preserved": True,
            "summary": "Unit-test fixture", "claims": ["Control-plane fixture only"], "limitations": ["No science executed"],
            "next_step": "Run the actual registered science independently",
            "evaluations": [{"purpose": "fixture", "command": ["python", "fixture.py"], "inputs_sha256": "b" * 64,
                             "exit_code": 1, "outcome": "FAIL", "evidence_paths": ["RESULTS.md"]}],
            "artifacts": rc.safe_files(self.root),
        })
        self.save()

    def tearDown(self):
        self.temp.cleanup()

    def save(self):
        rc.write(self.root / "ROUND_RESULT.json", self.result)

    def test_exact_requested_model_without_fallback(self):
        for model, effort in [("gpt-6.1-sol", "high"), ("other", "ultra")]:
            with self.assertRaises(ValueError): rc.select(self.queue, self.state, self.sha, "1-1", model, effort)

    def test_negative_outcome_preserved_and_idempotent(self):
        result, state, changed = rc.validate(self.context, self.state, self.root)
        self.assertTrue(changed)
        self.assertEqual(result["status"], "FAIL")
        self.assertEqual(state["records"][0]["status"], "FAIL")
        _, second, changed = rc.validate(self.context, state, self.root)
        self.assertFalse(changed)
        self.assertEqual(state, second)
        self.assertFalse(rc.select(self.queue, state, self.sha, "124-1", rc.MODEL, rc.EFFORT)["ready"])

    def test_passing_dependency_advances_exactly_one_task(self):
        self.result["status"] = "PASS"
        self.save()
        _, state, _ = rc.validate(self.context, self.state, self.root)
        following = rc.select(self.queue, state, self.sha, "124-1", rc.MODEL, rc.EFFORT)
        self.assertTrue(following["ready"])
        self.assertEqual(following["task"]["id"], "constrained-smooth-initial-data")

    def test_false_source_provenance_rejected(self):
        for field in ("source_sha", "scheduler_sha"):
            with self.subTest(field=field):
                result = copy.deepcopy(self.result)
                self.result[field] = "f" * 40
                self.save()
                with self.assertRaises(ValueError): rc.validate(self.context, self.state, self.root)
                self.result = result

    def test_tampered_evidence_hash_rejected(self):
        (self.root / "RESULTS.md").write_text("Tampered fixture\n")
        with self.assertRaises(ValueError): rc.validate(self.context, self.state, self.root)

    def test_missing_independent_review_rejected(self):
        (self.root / "INDEPENDENT_REVIEW.md").unlink()
        with self.assertRaises(ValueError): rc.validate(self.context, self.state, self.root)

    def test_duplicate_and_over_budget_evaluations_rejected(self):
        self.result["evaluations"] *= 2
        self.save()
        with self.assertRaises(ValueError): rc.validate(self.context, self.state, self.root)
        self.result["evaluations"] *= 2
        self.save()
        with self.assertRaises(ValueError): rc.validate(self.context, self.state, self.root)

    def test_symlink_and_hidden_control_file_rejected(self):
        target = self.root / "outside.md"
        target.symlink_to(Path(__file__))
        with self.assertRaises(ValueError): rc.safe_files(self.root)
        target.unlink()
        (self.root / ".github").mkdir()
        (self.root / ".github" / "workflow.json").write_text("{}")
        with self.assertRaises(ValueError): rc.safe_files(self.root)

    def test_empty_queue_stops_without_a_model_call(self):
        self.queue["items"] = []
        self.assertFalse(rc.select(self.queue, self.state, self.sha, "124-1", rc.MODEL, rc.EFFORT)["ready"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
