"""Install approved completed local evidence; remote publication remains pending."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path, PurePosixPath

CAPSULE = 'HDBLAST_COORDINATE_SHELL_SOURCES_20261015'
DELIVERY = 'hdblast/continuation/2026.10.15-coordinate-shell-sources'


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def canonical(name):
    require(isinstance(name, str) and name and '\\' not in name, 'Relative string path required')
    p = PurePosixPath(name)
    require(p.parts and not p.is_absolute() and '..' not in p.parts and p.as_posix() == name,
            'Canonical relative evidence path required')
    require(not ({part.lower() for part in p.parts} &
                 {'private', 'internal', '_private', '_internal', '.git', '__pycache__'}),
            'Private, internal or cache evidence forbidden')
    return name


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worktree', required=True, type=Path)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--replay-dir', type=Path, default=task/'fresh-zip-replay')
    args = parser.parse_args()
    config_raw = args.config.read_bytes()
    config = json.loads(config_raw)
    require(config['schema'] == 'hdblast-initial-evidence-config-v1'
            and config['status'] == 'APPROVED_FOR_INITIAL_EVIDENCE'
            and config['capsule'] == CAPSULE and isinstance(config['readme'], str)
            and config['readme'].strip(), 'Explicit approved initial evidence config required')
    destination = args.worktree.absolute()/DELIVERY
    require(not destination.exists() and not any(p.is_symlink() for p in [destination, *destination.parents]),
            'Fresh real initial evidence directory required')
    replay_raw = (args.replay_dir/'REPRODUCTION_RECEIPT.json').read_bytes()
    replay = json.loads(replay_raw)
    roster = json.loads((task/'frozen-source'/CAPSULE/'REPLAY_ROSTER.json').read_bytes())
    spec = importlib.util.spec_from_file_location('evidence_roster', task/'publication/build_final_evidence.py')
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    counts = validator.validate_roster(roster, replay)
    files = {'README.md': config['readme'].encode(), 'REPRODUCTION_RECEIPT.json': replay_raw}
    def add(name, raw):
        canonical(name)
        require(name != 'INITIAL_REPLAY_MANIFEST.json' and name not in files, 'Unique nonreserved evidence path required')
        files[name] = raw
    for row in config['files']:
        origin = canonical(row['source'])
        path = task/origin
        require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents])
                and path.resolve().is_relative_to(task.resolve()), 'Real completed evidence input required')
        raw = path.read_bytes()
        require(type(row['bytes']) is int and len(raw) == row['bytes']
                and sha256(raw).hexdigest() == row['sha256'], 'Approved evidence bytes differ: '+origin)
        add(row['path'], raw)
    for step in replay['steps']:
        for row in [step['log'], step['output']]+[x['output'] for x in step['additional_outputs']]:
            name = canonical(row['path'])
            path = args.replay_dir/name
            require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents]), 'Real replay artifact required')
            raw = path.read_bytes()
            require(len(raw) == row['bytes'] and sha256(raw).hexdigest() == row['sha256'], 'Replay artifact changed')
            add('fresh-replay/'+name, raw)
    required = {'SOURCE_FREEZE_RECEIPT.json', 'FRESH_ZIP_AUTHENTICATION.json',
                'FIGURE_INTEGRATION_RECEIPT.json', 'REPLAY_TOOL_VALIDATION.json'}
    require(required <= set(files), 'Freeze/authentication/figure/replay validation receipts required')
    freeze_raw = (task/'frozen-source/SOURCE_FREEZE_RECEIPT.json').read_bytes()
    freeze = json.loads(freeze_raw)
    validation_raw = (task/'publication/REPLAY_TOOL_VALIDATION.json').read_bytes()
    validation = json.loads(validation_raw)
    authentication_raw = (task/'FRESH_ZIP_AUTHENTICATION.json').read_bytes()
    authentication = json.loads(authentication_raw)
    figure_raw = (task/'publication/FIGURE_INTEGRATION_RECEIPT.json').read_bytes()
    source = task/'frozen-source'/CAPSULE
    manifest_raw = (source/'MANIFEST.json').read_bytes()
    require(files['SOURCE_FREEZE_RECEIPT.json'] == freeze_raw
            and files['REPLAY_TOOL_VALIDATION.json'] == validation_raw
            and files['FRESH_ZIP_AUTHENTICATION.json'] == authentication_raw
            and files['FIGURE_INTEGRATION_RECEIPT.json'] == figure_raw,
            'Initial evidence must preserve the actual current support receipts')
    require(freeze['status'] == authentication['status'] == validation['status'] == 'PASS'
            and freeze['capsule'] == validation['capsule'] == CAPSULE
            and len(manifest_raw) == freeze['manifest']['bytes']
            and sha256(manifest_raw).hexdigest() == freeze['manifest']['sha256']
            == replay['authentication']['manifest_sha256'] == validation['manifest_sha256'],
            'Replay, validation and current external freeze identity differ')
    require(authentication['archive'] == freeze['archive'] and authentication['manifest'] == freeze['manifest'],
            'Fresh ZIP authentication belongs to a different freeze')
    verifier_spec = importlib.util.spec_from_file_location('known_manifest_verifier', task/'publication/replay-tools/verify_manifest.py')
    verifier = importlib.util.module_from_spec(verifier_spec)
    verifier_spec.loader.exec_module(verifier)
    verifier.verify(source, freeze['manifest']['sha256'])
    for name, key in [('reproduce.py', 'runner'), ('verify_manifest.py', 'integrity_verifier')]:
        raw = (source/name).read_bytes()
        require(replay[key]['sha256'] == sha256(raw).hexdigest(), 'Replay frozen runner/verifier differs')
    replay_pin = validation['replay_receipt']
    require(replay_pin['bytes'] == len(replay_raw) and replay_pin['sha256'] == sha256(replay_raw).hexdigest(),
            'Initial replay differs from the completely validated replay bytes')
    expected_artifacts = {r['path']: r for r in validation['artifact_records']}
    copied_artifacts = {name.removeprefix('fresh-replay/'): raw for name, raw in files.items() if name.startswith('fresh-replay/')}
    require(len(expected_artifacts) == len(validation['artifact_records'])
            and set(expected_artifacts) == set(copied_artifacts), 'Complete validated replay artifact inventory required')
    for name, raw in copied_artifacts.items():
        require(expected_artifacts[name]['bytes'] == len(raw)
                and expected_artifacts[name]['sha256'] == sha256(raw).hexdigest(), 'Validated replay artifact bytes differ')
    manifest = {'status': 'PASS', 'capsule': CAPSULE,
                'scope': 'Completed local replay evidence; remote CI, merges and live publication remain pending.',
                'approved_config_sha256': sha256(config_raw).hexdigest(), 'replay_counts': counts,
                'files': [{'path': n, 'bytes': len(b), 'sha256': sha256(b).hexdigest()} for n, b in sorted(files.items())]}
    add('CONFIG_APPROVAL.json', config_raw)
    # Include the approved config itself in the complete manifest inventory.
    manifest['files'] = [{'path': n, 'bytes': len(b), 'sha256': sha256(b).hexdigest()} for n, b in sorted(files.items())]
    files['INITIAL_REPLAY_MANIFEST.json'] = (json.dumps(manifest, indent=2)+'\n').encode()
    destination.mkdir(parents=True)
    for name, raw in files.items():
        path = destination/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
        require(path.read_bytes() == raw, 'Initial evidence write/readback differs')
    print(json.dumps({'status': 'PASS', 'delivery': DELIVERY, 'files': len(files), 'replay_counts': counts}, indent=2))


if __name__ == '__main__':
    main()
