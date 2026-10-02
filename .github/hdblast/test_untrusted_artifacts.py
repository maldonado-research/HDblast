"""Independent negative tests; all fixtures are control-plane data, not science."""
import json
import math
import os
import subprocess
import sys
import unittest
from pathlib import Path

STAGED = Path(__file__).resolve().parent
sys.path.insert(0, str(STAGED))
import round_control as rc
from test_round_control import RoundControlTests


class UntrustedArtifactTests(RoundControlTests):
    def test_string_command_is_rejected(self):
        self.result['evaluations'][0]['command'] = 'python fixture.py'
        self.save()
        with self.assertRaises(ValueError):
            rc.validate(self.context, self.state, self.root)

    def test_claims_and_limitations_require_string_arrays(self):
        for key in ('claims', 'limitations'):
            original = self.result[key]
            self.result[key] = 'Untrusted output has the wrong schema'
            self.save()
            with self.assertRaises(ValueError):
                rc.validate(self.context, self.state, self.root)
            self.result[key] = original

    def test_symlink_artifact_root_is_rejected(self):
        alias = self.root.parent / 'alias'
        alias.symlink_to(self.root, target_is_directory=True)
        with self.assertRaises(ValueError):
            rc.safe_files(alias)

    def test_json_overflow_metadata_is_rejected(self):
        self.result['untrusted_numeric_metadata'] = math.inf
        payload = json.dumps(self.result).replace('Infinity', '1e999')
        (self.root / 'ROUND_RESULT.json').write_text(payload)
        with self.assertRaises(ValueError):
            rc.validate(self.context, self.state, self.root)

    def test_published_model_script_is_retained_without_execution(self):
        sentinel = self.root.parent / 'must-not-exist'
        (self.root / 'submitted.py').write_text(f"from pathlib import Path\nPath({str(sentinel)!r}).write_text('executed')\n")
        self.result['artifacts'] = {name: digest for name, digest in rc.safe_files(self.root).items()
                                    if name != 'ROUND_RESULT.json'}
        self.save()
        repo = self.root.parent / 'fresh-publish'
        (repo / 'research').mkdir(parents=True)
        context_file, state_file = self.root.parent / 'context.json', self.root.parent / 'state.json'
        rc.write(context_file, self.context)
        rc.write(state_file, self.state)
        completed = subprocess.run([
            sys.executable, str(STAGED / 'round_control.py'), 'publish',
            '--context', str(context_file), '--state', str(state_file), '--round-root', str(self.root),
            '--repo-root', str(repo), '--state-output', str(repo / 'round_state.json'),
            '--metadata-out', str(repo / 'metadata.json')], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertTrue((repo / self.context['round_path'] / 'submitted.py').is_file())
        self.assertFalse(sentinel.exists())


def load_tests(loader, tests, pattern):
    names = ['test_string_command_is_rejected', 'test_claims_and_limitations_require_string_arrays', 'test_symlink_artifact_root_is_rejected', 'test_json_overflow_metadata_is_rejected', 'test_published_model_script_is_retained_without_execution']
    return unittest.TestSuite(UntrustedArtifactTests(name) for name in names)


if __name__ == '__main__':
    unittest.main(verbosity=2)
