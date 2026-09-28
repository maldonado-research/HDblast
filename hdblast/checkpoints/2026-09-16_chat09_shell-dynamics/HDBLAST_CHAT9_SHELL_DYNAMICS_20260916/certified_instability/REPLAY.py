#!/usr/bin/env python3
"""Replay of the tachyon certificate (Python 3.9 standard library only; about 6-7 minutes).
  1. exact containment self-test of the interval / Taylor-model core;
  2. re-run of both certification configurations into *_REPLAY.json files;
  3. comparison of results, brackets and gates with the stored certificates (the computation is deterministic)."""
import json, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, HDB_CERT_SUFFIX="_REPLAY")
def run(args):
    r = subprocess.run([sys.executable] + args, cwd=here, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    print(r.stdout.strip().splitlines()[-1])
    if r.returncode != 0:
        raise SystemExit("FAILED: " + " ".join(args))
run(["selftest_core.py"])
ok = True
for args, name in (([], "TACHYON_CERTIFICATE"), (["-7.7180", "-7.7177", "-7.71788", "-7.71786"], "TACHYON_CERTIFICATE_TIGHT")):
    run(["certify_tachyon.py"] + args)
    a = json.load(open(os.path.join(here, name + ".json"))); b = json.load(open(os.path.join(here, name + "_REPLAY.json")))
    same = all(a[k] == b[k] for k in ("result", "certified_bracket", "gates", "shell_point_enclosures", "source_sha256"))
    print(name, "replay identical:", same, " bracket:", b["certified_bracket"]["float"], " gates:", len(b["gates"]), "all pass:", all(g["pass"] for g in b["gates"]))
    ok = ok and same
print("REPLAY", "PASS" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
