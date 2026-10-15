"""Publish the completed local replay as a dated, nonrecursive evidence directory."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import shutil

ROOT = Path(__file__).resolve().parent.parent
DELIVERY = 'hdblast/continuation/2026.10.15-uniform-first-variation-waveform'
destination = ROOT / 'publish-HDblast' / DELIVERY
if destination.exists():
    raise RuntimeError('New evidence directory required')
replay = json.loads((ROOT / 'fresh-zip-replay/REPRODUCTION_RECEIPT.json').read_text())
if replay['status'] != 'PASS':
    raise RuntimeError('Completed fresh replay required')
destination.mkdir(parents=True)
names = ['frozen-source/SOURCE_FREEZE_RECEIPT.json', 'FRESH_ZIP_AUTHENTICATION.json',
         'publication/SOURCE_ASSEMBLY_RECEIPT.json', 'publication/CI_INSTALL_PREFLIGHT.json',
         'publication/SUMMARY_INTEGRATION_PREFLIGHT.json',
         'publication/FIGURE_INTEGRATION_RECEIPT.json',
         'publication/GIT_INVENTORY_TOOL_VALIDATION.json',
         'publication/REPLAY_TOOL_VALIDATION.json',
         'publication/MONITOR_READBACK_TOOL_VALIDATION.json']
for name in names:
    target = destination / Path(name).name
    shutil.copyfile(ROOT / name, target)
shutil.copyfile(ROOT / 'fresh-zip-replay/REPRODUCTION_RECEIPT.json',
                destination / 'REPRODUCTION_RECEIPT.json')
for step in replay['steps']:
    results = [step['log'], step['output']]
    results += [row['output'] for row in step.get('additional_outputs', [])]
    for record in results:
        name = record['path']
        relative = PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix() != name:
            raise RuntimeError('Canonical replay result path required')
        source = ROOT / 'fresh-zip-replay' / name
        raw = source.read_bytes()
        if len(raw) != record['bytes'] or sha256(raw).hexdigest() != record['sha256']:
            raise RuntimeError('Replay result changed: ' + name)
        target = destination / 'fresh-replay' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.read_bytes() != raw:
            raise RuntimeError('Conflicting replay result path')
        target.write_bytes(raw)
(destination / 'README.md').write_text('''# Uniform first-variation waveform: initial replay evidence

This dated directory records the completed local freeze and authenticated
fresh-ZIP replay for `HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015`.
The source capsule contains the proof, exact code, continuous-cell data and
independent internal reviews. Its ZIP, manifest, runner and verifier have
external SHA256 pins in the dedicated repository workflow.

The replay receipt reports actual process counts, finite checks, embedded
controls and deliberately corrupted formula runs. Result logs, fresh JSON and
generated CSV outputs are included under `fresh-replay/`; copied execution
workspaces are excluded. These are computational
checks of the conditional mathematical checkpoint, not a formal certification
of its analytic PDE hypotheses, external peer review or cosmological evidence.

At this initial evidence stage, the present branch's eventual PR merge, remote
CI and live website readback are pending. Completed publication evidence will
be added separately, preserving these initial files. No new Zenodo version or
APS submission is part of this checkpoint.
''')
records = []
for path in sorted(destination.rglob('*')):
    if path.is_file():
        raw = path.read_bytes()
        records.append({'path': path.relative_to(destination).as_posix(),
                        'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})
manifest = {'status': 'PASS', 'scope': 'Initial completed local replay evidence; remote publication pending.',
            'files': records}
(destination / 'INITIAL_REPLAY_MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps({'status': 'PASS', 'delivery_path': DELIVERY, 'files': len(records) + 1}))
