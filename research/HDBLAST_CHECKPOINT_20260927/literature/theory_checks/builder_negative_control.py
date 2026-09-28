#!/usr/bin/env python3
"""Negative control for build_theory_md.py: corrupt one record (wrong arXiv id, nonexistent query
number) in a temporary copy and confirm the builder refuses to build. Writes builder_negative_control.json."""
import json, os, shutil, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
dst = os.path.join(tmp, 'lit', 'theory_checks'); shutil.copytree(HERE, dst)
p = os.path.join(dst, 'theory_sources.py'); s = open(p).read()
s = s.replace('url="https://arxiv.org/abs/1807.01570"', 'url="https://arxiv.org/abs/1807.01571"', 1)
s = s.replace('queries=[1]),', 'queries=[999]),', 1)
open(p, 'w').write(s)
r = subprocess.run([sys.executable, os.path.join(dst, 'build_theory_md.py')], capture_output=True, text=True)
res = dict(returncode=r.returncode, stdout=r.stdout.strip().splitlines(),
           expected="nonzero return code and both corruptions reported",
           passed=(r.returncode != 0 and 'bad query 999' in r.stdout and '1807.01571' in r.stdout),
           wrote_into_real_folder=False)
shutil.rmtree(tmp)
json.dump(res, open(os.path.join(HERE, 'builder_negative_control.json'), 'w'), indent=1)
print(json.dumps(res, indent=1)); sys.exit(0 if res['passed'] else 1)
