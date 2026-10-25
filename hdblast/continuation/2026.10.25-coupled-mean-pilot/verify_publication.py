#!/usr/bin/env python3
"""Verify the published source/ZIP against external pins before native replay."""
import argparse
import hashlib
import json
import subprocess
import sys
import zipfile
from pathlib import Path

def need(ok,message):
    if not ok:
        raise ValueError(message)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def relative(root, text):
    path=Path(text)
    need(not path.is_absolute() and ".." not in path.parts,"relative publication path")
    return root/path

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root",type=Path,required=True)
    p.add_argument("--pins",type=Path,required=True)
    p.add_argument("--pins-sha256",required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    repo=a.repo_root.resolve()
    raw=a.pins.read_bytes()
    need(sha(raw)==a.pins_sha256,"external publication pins")
    pins=json.loads(raw)
    source=relative(repo,pins["capsule"]["path"])
    archive=relative(repo,pins["archive"]["path"])
    manifest=(source/"MANIFEST.json").read_bytes()
    need(sha(manifest)==pins["capsule"]["manifest_sha256"],"source manifest")
    need(sha((source/"reproduce.py").read_bytes())==pins["runner_sha256"],"source runner")
    frozen=archive.read_bytes()
    need(len(frozen)==pins["archive"]["bytes"] and sha(frozen)==pins["archive"]["sha256"],
         "source archive bytes")
    entries=json.loads(manifest)["files"]
    expected={source.name+"/"+e["path"]:e for e in entries}
    expected[source.name+"/MANIFEST.json"]={"path":"MANIFEST.json","bytes":len(manifest),
                                          "sha256":sha(manifest)}
    with zipfile.ZipFile(archive) as z:
        infos=z.infolist()
        need(len(infos)==len(expected) and {i.filename for i in infos}==set(expected),
             "complete unique archive closure")
        need(z.testzip() is None,"archive CRC")
        for info in infos:
            need(not info.is_dir() and not (info.flag_bits&1)
                 and (info.external_attr>>16)&0o170000 != 0o120000,
                 "regular unencrypted archive member")
            row=expected[info.filename]
            member=z.read(info)
            need(len(member)==row["bytes"] and sha(member)==row["sha256"],"archive member hash")
            need(member==(source/row["path"]).read_bytes(),"archive/source equality")
    package={"status":"PASS_EXACT_SOURCE_ARCHIVE_PUBLICATION_PINS",
             "files_including_manifest":len(expected),
             "manifest_sha256":sha(manifest),"archive_sha256":sha(frozen),
             "archive_bytes":len(frozen)}
    print(json.dumps(package,sort_keys=True),flush=True)
    subprocess.run([sys.executable,"-I","-S","-B",str(source/"reproduce.py"),
                    "--manifest-sha256",pins["capsule"]["manifest_sha256"],
                    "--output-dir",str(a.output_dir),"--native","--repo-root",str(repo),
                    "--require-accuracy-pass"],check=True)
    need(sha(archive.read_bytes())==pins["archive"]["sha256"],"archive changed after replay")
    need(sha(a.pins.read_bytes())==a.pins_sha256,"pins changed after replay")
    (a.output_dir/"PUBLICATION_PACKAGE_RECEIPT.json").write_text(json.dumps(package,indent=2,sort_keys=True)+"\n")

if __name__=="__main__":
    main()
