import hashlib,json
from pathlib import Path
r=Path(__file__).resolve().parents[1]
hashes={"REGISTRATION.md":"939e7f6ad178c32cdcf638d01b8694d2194e805f0aad38964416518e3b72a2b5","PROTOCOL_ADDENDUM.md":"e7a98bea92852ac47a57a720e16efadd48444283ae71bbe72c5daa99575337c4"}
for p,h in hashes.items():
    assert hashlib.sha256((r/p).read_bytes()).hexdigest()==h,p
print(json.dumps({"status":"PASS","prospective_protocol_hashes":hashes}))
