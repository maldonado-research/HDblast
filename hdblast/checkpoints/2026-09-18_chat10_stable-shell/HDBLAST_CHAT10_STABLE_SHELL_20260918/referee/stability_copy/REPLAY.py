#!/usr/bin/env python3
"""Deterministic replay of the S_8/5 stability certificate (stdlib only).
  1. re-runs certify_stability.py into STABILITY_CERTIFICATE_S85_REPLAY.json and checks that it equals the delivered
     STABILITY_CERTIFICATE_S85.json in everything except runtime_seconds;
  2. negative control: the same run with d replaced by 13/10 in the shell condition (B = ... + t d/2; for d = 13/10 a scalar
     bound state exists at mu^2 = 1.2949 in floating point) MUST fail a junction-mismatch gate.
Usage: python3 REPLAY.py            (~5 min)"""
import os, sys, json, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
env = dict(os.environ, HDB_CERT_SUFFIX="_REPLAY")
r = subprocess.run([sys.executable, os.path.join(HERE, "certify_stability.py")], env=env, cwd=HERE)
if r.returncode != 0:
    raise SystemExit("REPLAY FAILED: certify_stability.py exited with status %d" % r.returncode)
a = json.load(open(os.path.join(HERE, "STABILITY_CERTIFICATE_S85.json"))); b = json.load(open(os.path.join(HERE, "STABILITY_CERTIFICATE_S85_REPLAY.json")))
a.pop("runtime_seconds"); b.pop("runtime_seconds")
print("replay identical to delivered certificate (except runtime):", a == b)
if a != b:
    raise SystemExit("REPLAY MISMATCH")
env = dict(os.environ, HDB_NEG_CONTROL_D="13/10", HDB_CERT_SUFFIX="_X")
r = subprocess.run([sys.executable, os.path.join(HERE, "certify_stability.py")], env=env, cwd=HERE, capture_output=True, text=True)
last = (r.stdout + r.stderr).strip().splitlines()[-1]
print("negative control (d = 13/10 in the shell condition): exit status %d; last line: %s" % (r.returncode, last))
if r.returncode == 0 or "GATE FAILED" not in last or "junction mismatch" not in last:
    raise SystemExit("NEGATIVE CONTROL DID NOT FAIL AS EXPECTED")
print("REPLAY OK: certificate reproduced, negative control correctly rejected.")
