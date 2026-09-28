#!/usr/bin/env python3
"""Replay the Chat 12 symbolic sector checks.  Needs SymPy:
    "/Users/ricardomaldonado/Documents/D-Blast 3/.hdblast_venv/bin/python" REPLAY.py          (about 1 minute)
The floating-point spot check tensor_float_spotcheck.py is separate (system python3 + numpy, ~6 minutes, needs folder 146)."""
import json, subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__)); ok = True; out = {}
for script in ("verify_tensor_sector.py", "verify_vector_sector.py", "verify_special_harmonics.py"):
    r = subprocess.run([sys.executable, script], cwd=here, capture_output=True, text=True)
    n_pass = sum(1 for l in r.stdout.splitlines() if l.startswith("[PASS]")); n_fail = sum(1 for l in r.stdout.splitlines() if l.startswith("[FAIL]"))
    out[script] = dict(exit_code=r.returncode, passed=n_pass, failed=n_fail); ok &= (r.returncode == 0 and n_fail == 0)
    print("%-32s exit=%d PASS=%d FAIL=%d" % (script, r.returncode, n_pass, n_fail), flush=True)
json.dump(dict(all_pass=ok, scripts=out), open(os.path.join(here, "REPLAY_RESULT.json"), "w"), indent=1)
print("REPLAY:", "ALL PASS" if ok else "FAILED"); sys.exit(0 if ok else 1)
