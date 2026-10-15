"""Authenticate external freeze pins and extract matching bytes before scientific imports."""
from hashlib import sha256
import json
from pathlib import Path
import stat
import zipfile

ROOT = Path(__file__).resolve().parent.parent
CAP = 'HDBLAST_UNIFORM_FIRST_VARIATION_WAVEFORM_20261015'
frozen = ROOT / 'frozen-source'
receipt = json.loads((frozen / 'SOURCE_FREEZE_RECEIPT.json').read_text())
archive = frozen / (CAP + '.zip')
raw = archive.read_bytes()
if sha256(raw).hexdigest() != receipt['archive']['sha256'] or len(raw) != receipt['archive']['bytes']:
    raise RuntimeError('External archive pin differs')
source = frozen / CAP
manifest_raw = (source / 'MANIFEST.json').read_bytes()
if sha256(manifest_raw).hexdigest() != receipt['manifest']['sha256']:
    raise RuntimeError('External manifest pin differs')
manifest = json.loads(manifest_raw)
expected = {CAP + '/MANIFEST.json': manifest_raw}
for row in manifest['files']:
    path = source / row['path']
    if path.is_symlink() or not path.is_file():
        raise RuntimeError('Regular source input required')
    content = path.read_bytes()
    if len(content) != row['bytes'] or sha256(content).hexdigest() != row['sha256']:
        raise RuntimeError('Frozen manifest entry differs: ' + row['path'])
    expected[CAP + '/' + row['path']] = content
destination = ROOT / 'fresh-zip-source'
if destination.exists():
    raise RuntimeError('Fresh ZIP extraction required')
with zipfile.ZipFile(archive) as zipped:
    names = zipped.namelist()
    if len(names) != len(set(names)) or set(names) != set(expected):
        raise RuntimeError('ZIP inventory differs')
    for info in zipped.infolist():
        if stat.S_IFMT(info.external_attr >> 16) != stat.S_IFREG:
            raise RuntimeError('Regular ZIP entries required')
        if zipped.read(info.filename) != expected[info.filename]:
            raise RuntimeError('ZIP/source bytes differ')
    destination.mkdir()
    for name, content in expected.items():
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        if path.read_bytes() != content:
            raise RuntimeError('Extracted bytes differ')
pins = {name: sha256((destination / CAP / name).read_bytes()).hexdigest()
        for name in ['MANIFEST.json', 'reproduce.py', 'verify_manifest.py']}
result = {'status': 'PASS', 'archive': receipt['archive'], 'manifest': receipt['manifest'],
          'bootstrap_pins_before_execution': pins, 'extracted_members': len(expected),
          'scope': 'External freeze/ZIP/source byte authentication and extraction before local scientific execution.'}
(ROOT / 'FRESH_ZIP_AUTHENTICATION.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result, indent=2))
