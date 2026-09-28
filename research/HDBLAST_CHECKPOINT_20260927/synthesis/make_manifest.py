#!/usr/bin/env python3
"""SHA-256 manifest of every file under the checkpoint folder (run last).

  python3 synthesis/make_manifest.py          # writes MANIFEST.sha256.json at the checkpoint root
  python3 synthesis/make_manifest.py --check  # re-hashes and reports added / missing / changed files

The manifest file itself is excluded. Symlinks are recorded by target path, not
followed. A negative control in --check mode confirms that a one-byte change of
an in-memory copy produces a different hash.
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NAME = "MANIFEST.sha256.json"


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def scan():
    files, links = {}, {}
    for dirpath, dirnames, filenames in os.walk(ROOT, followlinks=False):
        dirnames.sort()
        for d in list(dirnames):
            p = os.path.join(dirpath, d)
            if os.path.islink(p):
                links[os.path.relpath(p, ROOT)] = os.readlink(p)
                dirnames.remove(d)
        for fn in sorted(filenames):
            p = os.path.join(dirpath, fn)
            rel = os.path.relpath(p, ROOT)
            if rel == NAME:
                continue
            if os.path.islink(p):
                links[rel] = os.readlink(p)
                continue
            files[rel] = digest(p)
    return files, links


files, links = scan()
if "--check" in sys.argv:
    with open(os.path.join(ROOT, NAME)) as f:
        man = json.load(f)
    old = man["sha256"]
    added = sorted(set(files) - set(old))
    missing = sorted(set(old) - set(files))
    changed = sorted(k for k in set(files) & set(old) if files[k] != old[k])
    sample = next(iter(sorted(files)))
    with open(os.path.join(ROOT, sample), "rb") as f:
        data = bytearray(f.read())
    data.append(0)
    control_ok = hashlib.sha256(bytes(data)).hexdigest() != files[sample]
    ok = not (added or missing or changed) and control_ok
    print(json.dumps({"n_files": len(files), "added": added, "missing": missing, "changed": changed,
                      "negative_control_detects_change": control_ok, "ok": ok}, indent=1))
    sys.exit(0 if ok else 1)

out = {"root": "research/HDBLAST_CHECKPOINT_20260927", "algorithm": "sha256",
       "note": "Every regular file under the checkpoint folder except this manifest; symlinks listed separately.",
       "n_files": len(files), "sha256": files, "symlinks": links}
with open(os.path.join(ROOT, NAME), "w") as f:
    json.dump(out, f, indent=1, sort_keys=False)
print("files:", len(files), "symlinks:", len(links))
