"""SHA-256 manifest of every file in this folder (except the manifest itself)."""
import hashlib, json
from pathlib import Path
HERE = Path(__file__).resolve().parent
m = {str(p.relative_to(HERE)): hashlib.sha256(p.read_bytes()).hexdigest()
     for p in sorted(HERE.rglob('*')) if p.is_file() and p.name != 'MANIFEST.sha256.json' and '__pycache__' not in p.parts}
(HERE / 'MANIFEST.sha256.json').write_text(json.dumps(m, indent=1) + '\n')
print(len(m), 'files')
