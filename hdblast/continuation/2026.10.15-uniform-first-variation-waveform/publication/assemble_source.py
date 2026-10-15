"""Assemble only approved science, reviews and minimal inputs into a fresh tree."""
from hashlib import sha256
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'assembled-source'
if OUT.exists():
    raise RuntimeError('Fresh assembly directory required')
OUT.mkdir()
top = ['README.md', 'UNIFORM_FIRST_VARIATION_WAVEFORM.md', 'requirements.txt',
       'waveform_intervals.py', 'compute_waveform.py', 'check_waveform.py',
       'CHECK_WAVEFORM_NORMAL.json', 'plot_waveform.py']
records = []


def copy(path, relative):
    raw = path.read_bytes()
    target = OUT / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
    if target.read_bytes() != raw:
        raise RuntimeError('Assembly copy differs: ' + relative)
    records.append({'path': relative, 'bytes': len(raw), 'sha256': sha256(raw).hexdigest()})


for name in top:
    copy(ROOT / name, name)
for folder in ['inputs', 'archive', 'astra', 'sol', 'referee', 'computed', 'figures', 'literature']:
    for path in sorted((ROOT / folder).rglob('*')):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if set(relative.parts) & {'__pycache__', 'internal', 'final-review-replay'}:
            continue
        if path.suffix not in {'.md', '.json', '.py', '.txt', '.csv', '.png', '.svg'}:
            raise RuntimeError('Unreviewed assembly file type: ' + str(relative))
        copy(path, relative.as_posix())
for name in ['reproduce.py', 'verify_manifest.py']:
    copy(ROOT / 'publication/replay-tools' / name, name)
receipt = {'status': 'PASS', 'scope': 'Byte-preserving source assembly, excluding raw retrieved bodies and scratch workspaces.',
           'source_files': len(records), 'files': records}
(ROOT / 'publication/SOURCE_ASSEMBLY_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps({'status': 'PASS', 'source_files': len(records)}))
