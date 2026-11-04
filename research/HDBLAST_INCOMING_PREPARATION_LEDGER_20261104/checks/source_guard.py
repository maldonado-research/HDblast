#!/usr/bin/env python3
"""Standalone standard-library closed-source authentication; no project imports."""
import argparse,hashlib,json,sys
from pathlib import Path,PurePosixPath
HERE=Path(__file__).resolve().parent
CAP=HERE.parent
class Failure(Exception):pass
def need(ok,witness):
    if not ok:raise Failure(witness)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)
def relative(x):
    need(isinstance(x,str) and bool(x),'SAFE_RELATIVE_PATH');p=PurePosixPath(x)
    need(not p.is_absolute() and '..' not in p.parts and p.as_posix()==x and x!='.','SAFE_RELATIVE_PATH');return Path(*p.parts)
def authenticate(root,pin=None,old=False):
    prefix='SOURCE' if old else 'ROOT'
    need(root.is_dir() and not root.is_symlink(),prefix+'_CAPSULE_ROOT')
    manifest=root/'MANIFEST.json';need(regular(manifest),prefix+'_MANIFEST_REGULAR');actual_pin=sha(manifest)
    need(pin is None or actual_pin==pin,prefix+'_MANIFEST_PIN')
    rows=read(manifest)['files']
    if isinstance(rows,dict):rows=[{'path':k,**v} for k,v in rows.items()]
    need(isinstance(rows,list),prefix+'_MANIFEST_ROWS');seen={'MANIFEST.json'}
    for row in rows:
        rel=relative(row['path']);need(str(rel) not in seen,prefix+'_DUPLICATE_MEMBER');seen.add(str(rel));p=root/rel
        need(regular(p) and type(row['bytes']) is int and row['bytes']>=0 and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],prefix+'_MEMBER_BYTES')
    actual=set()
    for p in root.rglob('*'):
        need(not p.is_symlink(),prefix+'_SYMLINK')
        if p.is_file():actual.add(str(p.relative_to(root)))
        else:need(p.is_dir(),prefix+'_SPECIAL_FILE')
    need(actual==seen,prefix+'_CLOSED_INVENTORY')
    return {'manifest_sha256':actual_pin,'members_including_manifest':len(seen)}
def sources(repo,plan):
    result={name:authenticate(repo/'research'/relative(name),pin,True) for name,pin in sorted(plan['source_capsules'].items())}
    for rel,row in sorted(plan['selected_source_files'].items()):
        p=repo/relative(rel);need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SELECTED_OLD_SOURCE')
    for rel,row in sorted(plan['critical_new_files'].items()):
        p=CAP/relative(rel);need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SELECTED_NEW_SOURCE')
    return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--manifest-sha256');a=ap.parse_args()
 repo=a.repo_root.absolute();need(all(not p.is_symlink() for p in [repo,*repo.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True)
 own=authenticate(CAP,a.manifest_sha256);plan=read(HERE/'REPLAY_PLAN.json');need(plan['schema']=='HDBLAST_INCOMING_PREPARATION_LEDGER_ROOT_REPLAY_V1','REPLAY_PLAN_SCHEMA');need(plan.get('assembly_status')=='FROZEN_READY','ASSEMBLY_NOT_FROZEN');before=sources(repo,plan);need(len(before)==30,'COMPLETE_30_INHERITED_CLOSED_CAPSULES')
 need(own==authenticate(CAP,own['manifest_sha256']) and before==sources(repo,plan),'AUTHENTICATION_POSTCHECK')
 print(json.dumps({'status':'PASS','witness':'CLOSED_SOURCE_AUTHENTICATION','old_capsules':len(before)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Exception as error:
  print(json.dumps({'status':'FAIL','witness':str(error)},sort_keys=True));sys.exit(2)
