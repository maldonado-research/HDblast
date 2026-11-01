#!/usr/bin/env python3
"""Read-only closed authentication helper; imports no scientific dependency."""
import argparse,importlib.util,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;CAP=HERE.parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',required=True,type=Path);p.add_argument('--manifest-sha256');a=p.parse_args()
 s=importlib.util.spec_from_file_location('checkpoint_standard_library_guards',HERE/'run_checkpoint.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
 own=g.authenticate(CAP,a.manifest_sha256);plan=g.read(HERE/'REPLAY_PLAN.json');before=g.sources(a.repo_root.resolve(strict=True),plan)
 g.need(own==g.authenticate(CAP,own['manifest_sha256'])and before==g.sources(a.repo_root.resolve(strict=True),plan),'AUTHENTICATION_POSTCHECK')
 print(json.dumps({'status':'PASS','witness':'CLOSED_SOURCE_AUTHENTICATION','old_capsules':len(before)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Exception as exc:
  print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
