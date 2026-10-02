import hashlib
import json
from pathlib import Path

meta=json.loads(Path("INPUTS.json").read_text())
for name, record in meta["inputs"].items():
    p=Path("data/A1")/name
    assert p.stat().st_size == record["bytes"], name
    assert hashlib.sha256(p.read_bytes()).hexdigest() == record["sha256"], name
print(json.dumps({"status":"PASS","files":len(meta["inputs"])}))
