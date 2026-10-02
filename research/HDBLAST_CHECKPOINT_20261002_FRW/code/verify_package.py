"""Verify the curated checkpoint ZIP and its per-file SHA-256 manifest."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import sys
import zipfile

path = sys.argv[1] if len(sys.argv) > 1 else "HDBLAST_CHECKPOINT_20261002_FRW.zip"
root = "HDBLAST_CHECKPOINT_20261002_FRW/"
with zipfile.ZipFile(path) as archive:
    names = archive.namelist()
    if len(names) != len(set(names)):
        raise ValueError("Duplicate ZIP paths")
    for name in names:
        parts = PurePosixPath(name).parts
        if not name.startswith(root) or ".." in parts or name.startswith("/"):
            raise ValueError("Unexpected ZIP path: " + name)
    damaged = archive.testzip()
    if damaged:
        raise ValueError("CRC failure: " + damaged)
    manifest_name = root + "MANIFEST.sha256.json"
    manifest = json.loads(archive.read(manifest_name))
    expected = {root + p for p in manifest["sha256"]}
    if set(names) != expected | {manifest_name}:
        raise ValueError("ZIP inventory differs from manifest")
    for relative, expected_digest in manifest["sha256"].items():
        digest = hashlib.sha256(archive.read(root + relative)).hexdigest()
        if digest != expected_digest:
            raise ValueError("SHA-256 failure: " + relative)
        checked_out = Path(__file__).resolve().parent.parent / relative
        if not checked_out.is_file() or hashlib.sha256(checked_out.read_bytes()).hexdigest() != expected_digest:
            raise ValueError("Checked-out payload differs from ZIP: " + relative)
print(json.dumps({"status": "PASS", "files_checked": len(expected),
                  "zip": path, "algorithm": "sha256"}))
