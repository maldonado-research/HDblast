#!/usr/bin/env python3
"""Read-only closed authentication helper; imports no scientific dependency."""
import argparse,importlib.util,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;CAP=HERE.parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',required=True,type=Path);p.add_argument('--manifest-sha256');a=p.parse_args()
 s=importlib.util.spec_from_file_location('checkpoint_standard_library_guards',HERE/'run_checkpoint.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g)
 repo=a.repo_root.absolute();g.need(all(not p.is_symlink() for p in [repo,*repo.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True)
 own=g.authenticate(CAP,a.manifest_sha256);plan=g.read(HERE/'REPLAY_PLAN.json');g.need(plan['schema']=='HDBLAST_SHELL_PREPARATION_DEPENDENCE_ROOT_REPLAY_V1','REPLAY_PLAN_SCHEMA');before=g.sources(repo,plan);g.need(len(before)==27,'COMPLETE_27_INHERITED_CLOSED_CAPSULES')
 g.need(own==g.authenticate(CAP,own['manifest_sha256'])and before==g.sources(repo,plan),'AUTHENTICATION_POSTCHECK')
 print(json.dumps({'status':'PASS','witness':'CLOSED_SOURCE_AUTHENTICATION','old_capsules':len(before)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Exception as exc:
  print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
