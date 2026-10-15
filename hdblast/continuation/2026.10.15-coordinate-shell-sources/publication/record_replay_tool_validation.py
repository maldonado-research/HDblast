"""Bind actual authenticated replay artifacts, approved roster and tested support tools."""
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import importlib.util
import json
from pathlib import Path

CAPSULE = 'HDBLAST_COORDINATE_SHELL_SOURCES_20261015'


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def validate_output_inventory(configuration, replay, source, replay_dir):
    """Recover every declared complete output from the approved build config."""
    template = module(Path(__file__).resolve().parent/'runner_template.py', 'known_runner_template')
    checkers = {c['name']: c for c in configuration['checkers']}
    controls = {(c['checker'], c['name']): c for c in configuration['controls']}
    require(len(checkers) == len(configuration['checkers']), 'Unique configured checker names required')
    for step in replay['steps']:
        checker = checkers[step['checker']]
        prefix = step['mode']+'/'+step['checker']+'/'+step.get('control', 'positive')+'/'
        require(step['log']['path'] == prefix+'execution.log' and step['output']['path'] == prefix+'fresh.json',
                'Configured process log/result path differs')
        expected = checker.get('additional_outputs', []) if step['kind'] == 'positive' else []
        require(len(step['additional_outputs']) == len(expected), 'Complete configured additional output inventory required')
        if step['kind'] == 'positive':
            golden_name = checker.get('optimized_golden', checker['golden']) if step['mode'] == 'optimized' else checker['golden']
            golden_raw = (source/golden_name).read_bytes()
            result_raw = (replay_dir/step['output']['path']).read_bytes()
            require(template.json_equal(json.loads(result_raw), json.loads(golden_raw)),
                    'Complete fresh JSON differs from accepted mode golden')
            require(not checker.get('byte_equal', False) or result_raw == golden_raw, 'Complete deterministic JSON bytes differ')
        else:
            template.validate_failure_diagnostic(controls[(step['checker'], step['control'])],
                                                (replay_dir/step['log']['path']).read_text())
        for row, (fresh_name, golden_name) in zip(step['additional_outputs'], expected):
            golden_raw = (source/golden_name).read_bytes()
            require(row['output']['path'] == prefix+fresh_name and row['golden']['path'] == golden_name
                    and row['golden']['bytes'] == len(golden_raw) and row['golden']['sha256'] == sha256(golden_raw).hexdigest()
                    and row['complete_bytes_equal'] is True
                    and (replay_dir/row['output']['path']).read_bytes() == golden_raw,
                    'Complete configured additional bytes/golden identity differ')


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=task/'fresh-zip-source'/CAPSULE)
    parser.add_argument('--replay-dir', type=Path, default=task/'fresh-zip-replay')
    parser.add_argument('--output', type=Path, default=task/'publication/REPLAY_TOOL_VALIDATION.json')
    parser.add_argument('--review', action='append', type=Path, default=[])
    args = parser.parse_args()
    publication = task/'publication'
    require(not args.output.exists() and not any(p.is_symlink() for p in [args.output, *args.output.parents])
            and not args.output.resolve().is_relative_to(args.source.resolve())
            and not args.output.resolve().is_relative_to(args.replay_dir.resolve())
            and not args.output.resolve().is_relative_to((task/'frozen-source').resolve()),
            'New validation receipt must remain outside frozen source and replay artifacts')
    def record(path):
        require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents]), 'Real evidence file required')
        raw = path.read_bytes()
        return {'path': str(path), 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()}
    def authenticate(path, pin):
        actual = record(path)
        require(actual['bytes'] == pin['bytes'] and actual['sha256'] == pin['sha256'], 'Evidence bytes changed')
    replay_path = args.replay_dir/'REPRODUCTION_RECEIPT.json'
    replay = json.loads(replay_path.read_bytes())
    freeze_path = task/'frozen-source/SOURCE_FREEZE_RECEIPT.json'
    freeze = json.loads(freeze_path.read_bytes())
    manifest_raw = (args.source/'MANIFEST.json').read_bytes()
    manifest_pin = sha256(manifest_raw).hexdigest()
    require(freeze['status'] == 'PASS' and freeze['capsule'] == CAPSULE
            and freeze['manifest']['bytes'] == len(manifest_raw)
            and freeze['manifest']['sha256'] == manifest_pin, 'External source freeze pin differs')
    verified = module(publication/'replay-tools/verify_manifest.py', 'known_verifier').verify(args.source, manifest_pin)
    require(replay['authentication']['manifest_sha256'] == manifest_pin, 'Replay external manifest authentication differs')
    roster_path = args.source/'REPLAY_ROSTER.json'
    roster = json.loads(roster_path.read_bytes())
    counts = module(publication/'build_final_evidence.py', 'roster_validator').validate_roster(roster, replay)
    config_path = publication/'RUNNER_CONFIG.json'
    require(replay['embedded_configuration_sha256'] == record(config_path)['sha256'], 'Replay config differs from accepted build config')
    configuration = json.loads(config_path.read_bytes())
    require(configuration['status'] == 'APPROVED_FOR_FROZEN_REPLAY' and configuration['capsule'] == CAPSULE,
            'Approved current replay config required')
    build_path = publication/'replay-tools/RUNNER_BUILD_RECEIPT.json'
    build = json.loads(build_path.read_bytes())
    for row in build['source_input_pins']:
        authenticate(task/row['path'], row)
    for name in ('reproduce.py', 'REPLAY_ROSTER.json', 'verify_manifest.py'):
        require((args.source/name).read_bytes() == (publication/'replay-tools'/name).read_bytes(), 'Accepted runner/tool bytes differ')
    require(replay['runner']['sha256'] == record(args.source/'reproduce.py')['sha256']
            and replay['integrity_verifier']['sha256'] == record(args.source/'verify_manifest.py')['sha256'], 'Replay runner/verifier pins differ')
    validate_output_inventory(configuration, replay, args.source, args.replay_dir)
    artifacts = {}
    for step in replay['steps']:
        for row in [step['log'], step['output']]+[x['output'] for x in step['additional_outputs']]:
            name = row['path']
            require(isinstance(name, str) and name not in artifacts
                    and not Path(name).is_absolute() and '..' not in Path(name).parts and '\\' not in name,
                    'Unique canonical replay artifact required')
            authenticate(args.replay_dir/name, row)
            artifacts[name] = row
        if step['kind'] == 'external_corruption':
            failure = json.loads((args.replay_dir/step['output']['path']).read_bytes())
            log = (args.replay_dir/step['log']['path']).read_text()
            require(failure['status'] == 'FAIL' and failure['failure'] == step['expected_reason'] == step['reported_failure']
                    and failure['complete_diagnostic'] == log and failure['exit_code'] == step['exit_code'] == 1
                    and failure['passing_output_absent'] is True, 'Exact named failure artifact consistency required')
    synthetic_path = publication/'RUNNER_PLUMBING_VALIDATION.json'
    synthetic = json.loads(synthetic_path.read_bytes())
    require(synthetic['status'] == 'PASS_SYNTHETIC_REPLAY_PLUMBING'
            and synthetic['rejected_integrity_and_contract_cases']
            and all(synthetic['rejected_integrity_and_contract_cases'].values()), 'Passing synthetic corruption hardening required')
    for name, pin in synthetic['helper_sha256'].items():
        require(record(publication/name)['sha256'] == pin, 'Tested helper changed')
    result = {'status': 'PASS', 'observed_utc': datetime.now(timezone.utc).isoformat(), 'capsule': CAPSULE,
              'validation_basis': 'authenticated_fresh_zip_replay', 'source_files': verified['files_verified'],
              'manifest_sha256': manifest_pin, 'runner': record(args.source/'reproduce.py'), 'roster': record(roster_path),
              'verifier': record(args.source/'verify_manifest.py'), 'approved_build_config': record(config_path),
              'build_receipt': record(build_path), 'synthetic_plumbing_validation': record(synthetic_path),
              'independent_reviews': [record(p) for p in args.review], 'replay_receipt': record(replay_path),
              'source_freeze_receipt': record(freeze_path),
              'fresh_processes': len(replay['steps']), 'replay_counts': counts,
              'fresh_artifact_files_excluding_receipt': len(artifacts),
              'fresh_artifact_files_including_receipt': len(artifacts)+1,
              'source_input_pins_rechecked_unchanged': len(build['source_input_pins']),
              'artifact_records': [artifacts[n] for n in sorted(artifacts)],
              'scope': 'Actual authenticated fresh-ZIP replay and publication-tool validation. Counts derive from the approved roster. No remote CI, merge, live publication, full U2 response or analytic PDE proof is certified.'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'counts': counts, 'artifacts': len(artifacts)}, indent=2))


if __name__ == '__main__':
    main()
