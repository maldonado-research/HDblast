#!/usr/bin/env python3
"""Replay all Chat 11 symbolic verifications.  Needs SymPy; use the isolated environment created for the project:
    "/Users/ricardomaldonado/Documents/D-Blast 3/.hdblast_venv/bin/python" REPLAY.py
(or any Python >= 3.9 with sympy >= 1.12).  Exits non-zero on any failure.  Runtime about 30 seconds."""
import json, subprocess, sys, os
here = os.path.dirname(os.path.abspath(__file__)); ok = True; out = {}
for script in ("verify_linearised_equations.py", "verify_junction_logic.py", "verify_certificate_odes.py", "negative_controls.py"):
    r = subprocess.run([sys.executable, script], cwd=here, capture_output=True, text=True)
    last = r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-300:]
    n_pass = sum(1 for l in r.stdout.splitlines() if l.startswith("[PASS]")); n_fail = sum(1 for l in r.stdout.splitlines() if l.startswith("[FAIL]"))
    out[script] = dict(exit_code=r.returncode, passed=n_pass, failed=n_fail, last_line=last); ok &= (r.returncode == 0)
    print("%-36s exit=%d  PASS=%d FAIL=%d  | %s" % (script, r.returncode, n_pass, n_fail, last), flush=True)
json.dump(dict(all_pass=ok, scripts=out), open(os.path.join(here, "REPLAY_RESULT.json"), "w"), indent=1)
print("REPLAY:", "ALL PASS" if ok else "FAILED"); sys.exit(0 if ok else 1)
