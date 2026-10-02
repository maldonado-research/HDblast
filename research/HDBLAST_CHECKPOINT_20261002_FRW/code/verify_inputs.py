import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
d=json.loads((root/"INPUTS.json").read_text())
for name,meta in d["inputs"].items():
    b=(root/"data"/"A1"/name).read_bytes()
    assert len(b)==meta["bytes"],name
    assert hashlib.sha256(b).hexdigest()==meta["sha256"],name
print(json.dumps({"status":"PASS","input_hashes":len(d["inputs"])}))

