import hashlib
from pathlib import Path
checks={"REGISTRATION.md":"fc72bf086bb7747af0a51848ef997bf11d9d4451ee894aef899ea0ea8ca2888d","PROTOCOL_ADDENDUM.md":"52ed41cff38c6ee9dbdc9f86d109922238b71a29c2953f1d5f2638c10d366066"}
for name,digest in checks.items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,name
print("Protocol and pre-execution addendum hashes PASS")
