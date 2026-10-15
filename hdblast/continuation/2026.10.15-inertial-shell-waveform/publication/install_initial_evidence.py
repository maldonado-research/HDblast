"""Install completed local replay and reviewed planning; no remote claims."""
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath

ROOT=Path(__file__).resolve().parent.parent
DELIVERY='hdblast/continuation/2026.10.15-inertial-shell-waveform'
destination=ROOT/'publication/hdblast-worktree'/DELIVERY
if destination.exists():
    raise RuntimeError('Fresh evidence directory required')
replay=json.loads((ROOT/'fresh-zip-replay/REPRODUCTION_RECEIPT.json').read_text())
if replay['status']!='PASS':
    raise RuntimeError('Completed fresh ZIP replay required')
destination.mkdir(parents=True)

def copy(relative,target_name=None):
    source=ROOT/relative
    if not source.is_file() or source.is_symlink():
        raise RuntimeError('Regular evidence input required: '+relative)
    target=destination/(target_name or Path(relative).name)
    raw=source.read_bytes()
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists() and target.read_bytes()!=raw:
        raise RuntimeError('Conflicting initial evidence path')
    target.write_bytes(raw)

for name in ['frozen-source/SOURCE_FREEZE_RECEIPT.json','FRESH_ZIP_AUTHENTICATION.json',
             'publication/SOURCE_ASSEMBLY_RECEIPT.json','publication/PUBLICATION_SCOPE_CHECK.json',
             'publication/CI_INSTALL_PREFLIGHT.json','publication/SITE_INTEGRATION_REVIEW.json',
             'publication/FIGURE_INTEGRATION_RECEIPT.json','publication/REPLAY_TOOL_VALIDATION.json',
             'publication/MONITOR_READBACK_TOOL_VALIDATION.json','publication/REPRODUCE.md']:
    copy(name)
copy('fresh-zip-replay/REPRODUCTION_RECEIPT.json')
for step in replay['steps']:
    records=[step['log'],step['output']]+[x['output'] for x in step.get('additional_outputs',[])]
    for record in records:
        name=record['path'];relative=PurePosixPath(name)
        if relative.is_absolute() or '..' in relative.parts or relative.as_posix()!=name:
            raise RuntimeError('Canonical replay path required')
        raw=(ROOT/'fresh-zip-replay'/name).read_bytes()
        if len(raw)!=record['bytes'] or sha256(raw).hexdigest()!=record['sha256']:
            raise RuntimeError('Completed replay bytes changed')
        copy('fresh-zip-replay/'+name,'fresh-replay/'+name)
for path in sorted((ROOT/'planning').glob('*')):
    if path.suffix in {'.md','.py','.json'}:
        copy(path.relative_to(ROOT).as_posix(),'planning/'+path.name)
(destination/'README.md').write_text('''# Inertial shell waveform: completed local replay evidence

This directory records the completed source freeze and authenticated fresh-ZIP
replay for `HDBLAST_INERTIAL_SHELL_WAVEFORM_20261015`. The source contains the
conditional linear derivation, exact arithmetic, continuous-cell enclosures
and independent internal reviews. The manifest, ZIP, runner and verifier have
external SHA256 pins in the dedicated repository workflow and receipts here.

The replay receipt records actual processes, finite check groups and deliberate
corruption controls in normal and optimized Python modes. Fresh JSON, full CSV
and process logs are retained; execution workspaces are excluded. Planning
identities are separately checked and are not included in the frozen replay's
counts. These checks do not formally certify the analytic PDE hypotheses or
constitute external peer review.

At this initial stage, remote CI, merge and live readback remain pending.
Completed publication evidence will be added separately, preserving these
initial bytes. No Zenodo version or APS submission is part of this checkpoint.
''')
records=[]
for path in sorted(destination.rglob('*')):
    if path.is_file():
        raw=path.read_bytes()
        records.append({'path':path.relative_to(destination).as_posix(),'bytes':len(raw),'sha256':sha256(raw).hexdigest()})
(destination/'INITIAL_REPLAY_MANIFEST.json').write_text(json.dumps({'status':'PASS',
    'scope':'Completed local replay and separately reviewed planning; remote publication pending.',
    'files':records},indent=2)+'\n')
print(json.dumps({'status':'PASS','delivery':DELIVERY,'files':len(records)+1}))
