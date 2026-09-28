#!/usr/bin/env python3
"""Replay of the S_8/5 existence certificate (stdlib only, ~20 s).
  1. re-runs certify_existence.py into EXISTENCE_CERTIFICATE_S85_REPLAY.json and compares it with the archived
     certificate (everything except runtime_seconds must be identical: the arithmetic is deterministic integer arithmetic);
  2. negative control: the same script with the p-box shifted by 4 radii off the root must FAIL a Miranda gate.
CENTER_GUESS.json is non-rigorous guidance (regenerate with  python3 find_center.py 5  if desired, ~35 s)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, HDB_CERT_SUFFIX="_REPLAY")
r = subprocess.run([sys.executable, os.path.join(HERE, "certify_existence.py")], cwd=HERE, env=env, capture_output=True, text=True)
print(r.stdout[-600:])
if r.returncode != 0:
    raise SystemExit("REPLAY FAILED: " + r.stderr[-500:])
a = json.load(open(os.path.join(HERE, "EXISTENCE_CERTIFICATE_S85.json")))
b = json.load(open(os.path.join(HERE, "EXISTENCE_CERTIFICATE_S85_REPLAY.json")))
a.pop("runtime_seconds"); b.pop("runtime_seconds")
if a != b:
    raise SystemExit("REPLAY MISMATCH with the archived certificate")
print("replay identical to the archived certificate; %d gates" % len(a["gates"]))
env = dict(os.environ, HDB_CERT_SUFFIX="_REPLAY", HDB_NEG_CONTROL_SHIFT="4")
r = subprocess.run([sys.executable, os.path.join(HERE, "certify_existence.py")], cwd=HERE, env=env, capture_output=True, text=True)
msg = (r.stderr.strip().splitlines() or [""])[-1]
if r.returncode == 0 or "GATE FAILED" not in msg:
    raise SystemExit("NEGATIVE CONTROL DID NOT FAIL AS EXPECTED")
print("negative control (box shifted off the root) correctly rejected:", msg[:160])
print("REPLAY PASS")
