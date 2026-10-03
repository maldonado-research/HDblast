#!/usr/bin/env python3
"""Synthetic integrity, outcome, path and resource guards; zero physical inputs."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import zipfile

import build_package
import ledger_integrity as integrity
import replay_ledger
import verify_public_freeze

FREEZE = '1' * 40
TOY_ROUTE = '''import argparse, json, pathlib, sys
p=argparse.ArgumentParser()
for n in ('checkpoint-root','registration-sha256','freeze-commit','output-dir'): p.add_argument('--'+n, required=True)
a=p.parse_args(); out=pathlib.Path(a.output_dir); out.mkdir(parents=True, exist_ok=False)
(out/'progress.jsonl').write_text(json.dumps({'toy':True,'completed_cases':1})+'\\n')
{fault}
(out/'diagnostic.json').write_text(json.dumps({'toy':True,'declared_status':'completed'})+'\\n')
'''
TOY_VALIDATOR = '''import argparse, json, pathlib, sys
p=argparse.ArgumentParser()
for n in ('checkpoint-root','registration-sha256','freeze-commit','output-dir','primary','independent','execution-receipts'): p.add_argument('--'+n, required=True)
a=p.parse_args(); receipts=json.loads(pathlib.Path(a.execution_receipts).read_text())
routes=[r for r in receipts if r['physical_route']]
if len(routes)!=2 or any(r['status']!='PASS_EXECUTION' or r['exit_code']!=0 for r in routes): sys.exit(19)
out=pathlib.Path(a.output_dir); out.mkdir(parents=True, exist_ok=False)
(out/'VALIDATION.json').write_text(json.dumps({'classification':{classification!r}, 'old_metric_status':'FAIL','python_optimization':sys.flags.optimize,'toy':True})+'\\n')
'''


class GuardTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='ledger-synthetic-guards-')
        self.base = Path(self.temporary.name)
        self.root = self.base / integrity.CHECKPOINT_NAME
        self.root.mkdir()
        (self.root / 'code').mkdir()
        (self.root / 'independent').mkdir()

    def tearDown(self):
        self.temporary.cleanup()

    def fixture(self, fault='', classification='NO_GATE_SCALE_ATTRIBUTION'):
        (self.root / 'code/diagnostic_primary.py').write_text(TOY_ROUTE.replace('{fault}', fault))
        (self.root / 'independent/diagnostic_independent.py').write_text(TOY_ROUTE.replace('{fault}', ''))
        (self.root / 'code/validate_ledger.py').write_text(TOY_VALIDATOR.replace('{classification!r}', repr(classification)))
        integrity.write_json(self.root / 'EXPERIMENT.json', {'synthetic_fixture': True, 'reproduction': {
            'routes': {'primary': 'code/diagnostic_primary.py', 'independent': 'independent/diagnostic_independent.py'},
            'validator': 'code/validate_ledger.py', 'result_filename': 'diagnostic.json', 'validation_filename': 'VALIDATION.json',
            'postfreeze_output_directories': ['recorded_outputs']}})
        self.register()

    def register(self, extra=None):
        files = {p.relative_to(self.root).as_posix(): integrity.file_sha(p)
                 for p in self.root.rglob('*') if p.is_file()
                 and p.name not in ('FULL_REGISTRATION.json', 'FREEZE_RECEIPT.json', 'MANIFEST.json')}
        if extra:
            files.update(extra)
        integrity.write_json(self.root / 'FULL_REGISTRATION.json', {'schema_version': 1, 'files': files})
        self.pin = integrity.file_sha(self.root / 'FULL_REGISTRATION.json')
        integrity.write_json(self.root / 'FREEZE_RECEIPT.json',
                             {'public_freeze_commit': FREEZE, 'registration_sha256': self.pin})

    def do_replay(self, output='replay', plan_only=False):
        # Unit-only overrides permit pure stdlib toys and tests under -O.
        # Real replay CLI has no such override and enforces all pinned packages.
        with patch.object(replay_ledger, 'runtime_check'), \
             patch.object(replay_ledger.sys, 'flags', SimpleNamespace(optimize=0)):
            return replay_ledger.replay(self.root, self.base / output, self.pin, FREEZE, plan_only)

    def test_frozen_valid_bytes(self):
        self.fixture()
        found = integrity.verify_frozen(self.root, self.pin, FREEZE)
        self.assertEqual(found['registration_sha256'], self.pin)

    def test_full_registration_mutation_rejected(self):
        self.fixture()
        with (self.root / 'FULL_REGISTRATION.json').open('a') as out:
            out.write(' ')
        with self.assertRaisesRegex(RuntimeError, 'registration SHA256'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_frozen_source_mutation_rejected(self):
        self.fixture()
        (self.root / 'code/diagnostic_primary.py').write_text('modified')
        with self.assertRaisesRegex(RuntimeError, 'Frozen input changed'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_receipt_commit_and_sha_rejected(self):
        self.fixture()
        for field, value in (('public_freeze_commit', '2' * 40), ('registration_sha256', '2' * 64)):
            receipt = {'public_freeze_commit': FREEZE, 'registration_sha256': self.pin, field: value}
            integrity.write_json(self.root / 'FREEZE_RECEIPT.json', receipt)
            with self.assertRaisesRegex(RuntimeError, 'Freeze receipt'):
                integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_path_traversal_absolute_empty_and_backslash_rejected(self):
        self.fixture()
        for name in ('../external', '/tmp/external', 'code//file', './file', 'code/../file', 'x\\y', 'C:/file'):
            with self.subTest(name=name), self.assertRaises(RuntimeError):
                integrity.safe_file(self.root, name)

    def test_registered_symlink_file_and_directory_rejected(self):
        self.fixture()
        external = self.base / 'external'
        external.mkdir()
        (external / 'value').write_text('toy')
        (self.root / 'link').symlink_to(external, target_is_directory=True)
        with self.assertRaisesRegex(RuntimeError, 'Symlink'):
            integrity.safe_file(self.root, 'link/value')
        with self.assertRaisesRegex(RuntimeError, 'Symlink'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_self_reference_and_omitted_experiment_rejected(self):
        self.fixture()
        self.register({'FULL_REGISTRATION.json': '0' * 64})
        with self.assertRaisesRegex(RuntimeError, 'Self-referential'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)
        registration = integrity.read_json(self.root / 'FULL_REGISTRATION.json')
        del registration['files']['FULL_REGISTRATION.json']; del registration['files']['EXPERIMENT.json']
        integrity.write_json(self.root / 'FULL_REGISTRATION.json', registration)
        self.pin = integrity.file_sha(self.root / 'FULL_REGISTRATION.json')
        integrity.write_json(self.root / 'FREEZE_RECEIPT.json', {'public_freeze_commit': FREEZE, 'registration_sha256': self.pin})
        with self.assertRaisesRegex(RuntimeError, 'Experiment omitted'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_duplicate_json_keys_and_nan_rejected(self):
        path = self.base / 'bad.json'
        for raw in ('{"a":1,"a":2}', '{"a":NaN}'):
            path.write_text(raw)
            with self.assertRaises(RuntimeError):
                integrity.read_json(path)

    def test_public_git_freeze_binds_exact_registered_bytes(self):
        self.fixture()
        def git(*args):
            return subprocess.check_output(['git', '-C', str(self.base), *args], stderr=subprocess.DEVNULL)
        git('init', '-q')
        git('add', integrity.CHECKPOINT_NAME)
        git('-c', 'user.name=Synthetic Guard', '-c', 'user.email=synthetic@example.invalid',
            'commit', '-q', '-m', 'Synthetic byte fixture only')
        freeze = git('rev-parse', 'HEAD').decode().strip()
        integrity.write_json(self.root / 'FREEZE_RECEIPT.json', {'public_freeze_commit': freeze, 'registration_sha256': self.pin})
        report = verify_public_freeze.verify(self.root, self.base, self.pin, freeze)
        self.assertEqual(report['status'], 'PASS_LOCAL_FREEZE_GIT_BYTES')
        (self.root / 'code/diagnostic_primary.py').write_text('changed but newly re-registered')
        self.register()
        integrity.write_json(self.root / 'FREEZE_RECEIPT.json', {'public_freeze_commit': freeze, 'registration_sha256': self.pin})
        with self.assertRaisesRegex(RuntimeError, 'Public freeze Git bytes differ'):
            verify_public_freeze.verify(self.root, self.base, self.pin, freeze)

    def test_plan_only_executes_zero_children(self):
        self.fixture()
        state = self.do_replay(plan_only=True)
        self.assertEqual(state['status'], 'VERIFIED_PLAN_ONLY')
        self.assertEqual(state['attempted_commands'], 0)
        plan = integrity.read_json(self.base / 'replay/COMMAND_PLAN.json')
        self.assertEqual(len([c for c in plan if c['physical_route']]), 2)
        self.assertTrue(all(c['timeout_seconds'] == 900 and c['memory_limit_kib'] == 262144 for c in plan))
        self.assertTrue(all('-O' not in c['command'] for c in plan if c['physical_route']))

    def test_negative_and_consistency_outcomes_are_completed_without_science_pass(self):
        for index, classification in enumerate(('NO_GATE_SCALE_ATTRIBUTION', 'CONSISTENCY_FAILURE', 'LEDGER_ERROR_DEMONSTRATED')):
            self.fixture(classification=classification)
            state = self.do_replay(output='replay-' + str(index))
            self.assertEqual(state['status'], 'COMPLETED_REPLAY', state)
            self.assertEqual(state['classification'], classification)
            self.assertEqual(state['old_metric_status'], 'FAIL')
            self.assertEqual(state['physical_routes_completed'], 2)
            self.assertTrue(all(c['exit_code'] == 0 and c['peak_rss_kib'] > 0 for c in state['commands']))

    def test_nonzero_route_retains_progress_and_real_exit(self):
        self.fixture(fault="(out/'diagnostic.json').write_text('{\\\"declared_status\\\":\\\"completed\\\"}'); sys.exit(7)")
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertEqual(state['commands'][0]['exit_code'], 7)
        self.assertEqual(state['physical_routes_attempted'], 2)
        self.assertEqual(state['classification'], 'UNCOMPUTED')
        self.assertTrue((self.base / 'replay/fresh/primary/progress.jsonl').is_file())
        self.assertEqual(len(state['skipped_commands']), 2)

    def test_zero_exit_without_diagnostic_is_failure_and_second_route_runs(self):
        self.fixture(fault='sys.exit(0)')
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertEqual(state['commands'][0]['exit_code'], 0)
        self.assertEqual(state['commands'][0]['status'], 'FAIL_EXECUTION')
        self.assertEqual(state['physical_routes_attempted'], 2)

    def test_normal_optimized_validator_disagreement_rejected(self):
        self.fixture()
        path = self.root / 'code/validate_ledger.py'
        path.write_text(path.read_text().replace("'toy':True", "'toy':sys.flags.optimize==0"))
        self.register()
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertIn('validator reports differ', state['failure']['exception'])

    def test_source_mutation_during_route_rejected_and_progress_retained(self):
        self.fixture(fault="pathlib.Path(a.checkpoint_root,'EXPERIMENT.json').write_text('changed')")
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertEqual(state['source_verification_after'], 'FAIL')
        self.assertTrue((self.base / 'replay/fresh/primary/progress.jsonl').is_file())

    def test_fresh_output_and_external_output_required(self):
        self.fixture()
        with self.assertRaisesRegex(RuntimeError, 'fresh output'):
            with patch.object(replay_ledger.sys, 'flags', SimpleNamespace(optimize=0)):
                replay_ledger.replay(self.root, self.root / 'inside', self.pin, FREEZE)
        self.do_replay(plan_only=True)
        with self.assertRaisesRegex(RuntimeError, 'fresh output'):
            self.do_replay()

    def command(self, code, timeout=2, memory=262144):
        logs = self.base / 'logs'; logs.mkdir()
        return replay_ledger.run_command({'name': 'toy', 'physical_route': True,
            'command': [sys.executable, '-c', code], 'cwd': str(self.base),
            'timeout_seconds': timeout, 'memory_limit_kib': memory}, logs, os.environ.copy())

    def test_wall_timeout_is_nonzero_actual_receipt(self):
        record = self.command('import time; print("progress",flush=True); time.sleep(30)', timeout=0.15)
        self.assertEqual(record['status'], 'FAIL_EXECUTION')
        self.assertTrue(record['timed_out'])
        self.assertLess(record['exit_code'], 0)
        self.assertTrue((self.base / 'logs/toy.stdout.log').read_text().startswith('progress'))

    def test_memory_budget_is_measured_and_enforced(self):
        record = self.command('import time; x=bytearray(80*1024*1024); time.sleep(30)', memory=30000)
        self.assertEqual(record['status'], 'FAIL_EXECUTION')
        self.assertTrue(record['memory_exceeded'])
        self.assertGreater(record['peak_rss_kib'], 30000)
        self.assertLess(record['exit_code'], 0)

    def test_launch_error_cannot_be_inferred_success(self):
        logs = self.base / 'logs'; logs.mkdir()
        record = replay_ledger.run_command({'name': 'missing', 'physical_route': True,
            'command': [str(self.base / 'missing-command')], 'cwd': str(self.base),
            'timeout_seconds': 2, 'memory_limit_kib': 262144}, logs, os.environ.copy())
        self.assertEqual(record['status'], 'FAIL_EXECUTION')
        self.assertIsNone(record['exit_code'])
        self.assertIn('exception', record)

    def test_orphan_descendant_escaping_process_group_is_rejected_and_reaped(self):
        code = "import subprocess,sys; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(30)'],start_new_session=True); print(p.pid,flush=True)"
        # Unit-only bypass independently tests the redundant orphan cleanup.
        # Actual CLI has no bypass for in-kernel child/thread prohibition.
        with patch.object(replay_ledger.resource, 'setrlimit'):
            record = self.command(code)
        self.assertEqual(record['status'], 'FAIL_EXECUTION')
        self.assertEqual(record['exit_code'], 0)
        pid = int((self.base / 'logs/toy.stdout.log').read_text().strip())
        self.assertFalse(Path('/proc', str(pid)).exists())
        self.assertTrue(record['descendant_exit_records'])

    def test_kernel_prevents_descendant_creation(self):
        code = "import subprocess,sys\ntry:\n subprocess.run([sys.executable,'-c','print(1)'],check=True)\nexcept BlockingIOError:\n print('KERNEL_PROCESS_CREATION_FORBIDDEN',flush=True)\nelse:\n sys.exit(31)"
        record = self.command(code)
        self.assertEqual(record['status'], 'PASS_EXECUTION')
        self.assertEqual(record['process_creation_limit'], 0)
        self.assertIn('KERNEL_PROCESS_CREATION_FORBIDDEN', (self.base / 'logs/toy.stdout.log').read_text())

    def test_unregistered_source_shadow_and_output_module_rejected(self):
        self.fixture()
        for name in ('code/json.py', 'recorded_outputs/attack.py', 'recorded_outputs/attack.zip'):
            path = self.root / name; path.parent.mkdir(exist_ok=True)
            path.write_text('toy')
            with self.assertRaisesRegex(RuntimeError, 'Unregistered payload'):
                integrity.verify_frozen(self.root, self.pin, FREEZE)
            path.unlink()

    def test_validator_symlink_failure_clears_unaccepted_classification(self):
        self.fixture(classification='LEDGER_ERROR_DEMONSTRATED')
        path = self.root / 'code/validate_ledger.py'
        path.write_text(path.read_text() + "(out/'badlink').symlink_to(pathlib.Path(a.primary))\n")
        self.register()
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertEqual(state['classification'], 'UNCOMPUTED')
        self.assertEqual(state['unaccepted_validator_classification'], 'LEDGER_ERROR_DEMONSTRATED')

    def test_unknown_scientific_outcome_rejected(self):
        self.fixture(classification='SCIENCE_PASS')
        state = self.do_replay()
        self.assertEqual(state['status'], 'FAIL_REPLAY')
        self.assertEqual(state['classification'], 'UNCOMPUTED')

    def test_deterministic_package_roundtrip_includes_recorded_outputs(self):
        self.fixture()
        (self.root / 'recorded_outputs').mkdir()
        (self.root / 'recorded_outputs/diagnostic.json').write_text('{"toy":true}\n')
        (self.root / '__pycache__').mkdir()
        (self.root / '__pycache__/cache.pyc').write_bytes(b'ignored cache')
        (self.root / integrity.ZIP_NAME).write_bytes(b'excluded package itself')
        a = build_package.build(self.root, self.base / 'package-a', self.pin, FREEZE)
        b = build_package.build(self.root, self.base / 'package-b', self.pin, FREEZE)
        self.assertEqual(a['zip_sha256'], b['zip_sha256'])
        self.assertFalse(a['science_evaluated'])
        manifest = integrity.read_json(self.root / 'MANIFEST.json')
        self.assertIn('recorded_outputs/diagnostic.json', manifest['files'])
        self.assertIn('FREEZE_RECEIPT.json', manifest['files'])
        self.assertNotIn(integrity.ZIP_NAME, manifest['files'])
        self.assertNotIn('__pycache__/cache.pyc', manifest['files'])
        extracted = build_package.extract_verified(a['zip_path'], self.base / 'extracted', a['package_manifest_sha256'])
        integrity.verify_frozen(extracted, self.pin, FREEZE)
        with patch.object(replay_ledger, 'runtime_check'), patch.object(replay_ledger.sys, 'flags', SimpleNamespace(optimize=0)):
            state = replay_ledger.replay(extracted, self.base / 'zip-replay', self.pin, FREEZE)
        self.assertEqual(state['status'], 'COMPLETED_REPLAY')

    def test_payload_after_manifest_mutation_rejected(self):
        self.fixture()
        build_package.build(self.root, self.base / 'package', self.pin, FREEZE)
        (self.root / 'unregistered-after-manifest').write_text('toy')
        with self.assertRaisesRegex(RuntimeError, 'payload'):
            integrity.verify_frozen(self.root, self.pin, FREEZE)

    def test_zip_traversal_duplicate_and_symlink_rejected(self):
        self.fixture()
        receipt = build_package.build(self.root, self.base / 'package', self.pin, FREEZE)
        for index, attack in enumerate(('traversal', 'duplicate', 'symlink')):
            path = self.base / ('attack-' + str(index) + '.zip')
            with zipfile.ZipFile(receipt['zip_path']) as original, zipfile.ZipFile(path, 'w') as changed:
                for item in original.infolist():
                    changed.writestr(item, original.read(item))
                info = zipfile.ZipInfo(integrity.CHECKPOINT_NAME + ('/../escaped' if attack == 'traversal' else '/FREEZE_RECEIPT.json' if attack == 'duplicate' else '/link'))
                info.create_system = 3
                info.external_attr = (0o120777 if attack == 'symlink' else 0o100644) << 16
                changed.writestr(info, 'toy')
            with self.subTest(attack=attack), self.assertRaises(RuntimeError):
                build_package.extract_verified(path, self.base / ('extract-' + str(index)), receipt['package_manifest_sha256'])
            self.assertFalse((self.base / ('extract-' + str(index))).exists())


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(GuardTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    print(json.dumps({'status': 'PASS_SYNTHETIC_GUARDS' if result.wasSuccessful() else 'FAIL_SYNTHETIC_GUARDS',
                      'tests_run': result.testsRun, 'physical_evaluations': 0,
                      'python_optimization': sys.flags.optimize}, sort_keys=True))
    raise SystemExit(0 if result.wasSuccessful() else 1)
