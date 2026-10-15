"""Install completed local replay evidence without claiming pending remote facts."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parent.parent
DELIVERY = 'hdblast/continuation/2026.10.15-quadratic-shell-forcing'
destination = ROOT/'publication/hdblast-worktree'/DELIVERY
if destination.exists():
    raise RuntimeError('Fresh initial evidence directory required')
replay = json.loads((ROOT/'fresh-zip-replay/REPRODUCTION_RECEIPT.json').read_bytes())
if replay['status'] != 'PASS':
    raise RuntimeError('Completed fresh ZIP replay required')
destination.mkdir(parents=True)


def copy(relative, target_name=None):
    source = ROOT/relative
    if not source.is_file() or source.is_symlink():
        raise RuntimeError('Regular completed evidence file required: '+relative)
    raw = source.read_bytes()
    target = destination/(target_name or Path(relative).name)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.read_bytes() != raw:
        raise RuntimeError('Conflicting initial evidence bytes')
    target.write_bytes(raw)


for name in ['frozen-source/SOURCE_FREEZE_RECEIPT.json', 'FRESH_ZIP_AUTHENTICATION.json',
             'publication/SOURCE_ASSEMBLY_RECEIPT.json', 'publication/PUBLICATION_SCOPE_CHECK.json',
             'publication/CI_INSTALL_PREFLIGHT.json', 'publication/SITE_INTEGRATION_REVIEW.json',
             'publication/FIGURE_INTEGRATION_RECEIPT.json', 'publication/REPLAY_TOOL_VALIDATION.json',
             'publication/MONITOR_READBACK_TOOL_VALIDATION.json', 'REPRODUCE.md']:
    copy(name)
copy('fresh-zip-replay/REPRODUCTION_RECEIPT.json')
for step in replay['steps']:
    records = [step['log'], step['output']]+[v['output'] for v in step['additional_outputs']]
    for row in records:
        name = row['path']
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or path.as_posix() != name:
            raise RuntimeError('Canonical replay output required')
        raw = (ROOT/'fresh-zip-replay'/name).read_bytes()
        if len(raw) != row['bytes'] or sha256(raw).hexdigest() != row['sha256']:
            raise RuntimeError('Completed replay artifact changed')
        copy('fresh-zip-replay/'+name, 'fresh-replay/'+name)
(destination/'README.md').write_text('''# Quadratic material tension contribution: completed local replay

This delivery records the source freeze and authenticated fresh-ZIP replay for
`HDBLAST_QUADRATIC_SHELL_FORCING_20261015`. The source supplies explicit quadratic
material tension terms, continuous-cell enclosures, angular identities and
independent internal reviews. Complete generated JSON/CSV and process logs from
normal and optimized Python are retained. Disposable executable copies and raw
literature responses are excluded. The frozen roster, rather than older counts,
defines every executed calculation, named group and deliberate corruption.

The calculation does not evaluate the complete source in bulk coordinates,
full U2, total second stress or a nonlinear remainder. The earlier conditional
complete second-coefficient construction is preserved. Finite arithmetic does
not prove the analytic PDE assumptions or replace external peer review.

Remote CI, merge and live readback are pending at this initial stage. Completed
publication facts will be added separately, preserving every initial byte.
No Zenodo version or APS submission is part of this checkpoint.
''')
records = []
for path in sorted(destination.rglob('*')):
    if path.is_file():
        raw = path.read_bytes()
        records.append({'path': path.relative_to(destination).as_posix(), 'bytes': len(raw),
                        'sha256': sha256(raw).hexdigest()})
(destination/'INITIAL_REPLAY_MANIFEST.json').write_text(json.dumps({'status': 'PASS',
    'scope': 'Completed local replay; remote publication pending.', 'files': records}, indent=2)+'\n')
print(json.dumps({'status': 'PASS', 'delivery': DELIVERY, 'files': len(records)+1}))
