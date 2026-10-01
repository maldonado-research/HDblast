"""Verify the curated checkpoint ZIP and its per-file SHA-256 manifest."""
import hashlib
import json
from pathlib import PurePosixPath
import sys
import zipfile

path = sys.argv[1] if len(sys.argv) > 1 else "HDBLAST_CHECKPOINT_20261001.zip"
root = "HDBLAST_CHECKPOINT_20261001/"
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
print(json.dumps({"status": "PASS", "files_checked": len(expected),
                  "zip": path, "algorithm": "sha256"}))
