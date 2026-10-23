#!/usr/bin/env python3
"""Authenticate every declared byte in each selected capsule before import."""
import hashlib,json
from pathlib import Path
ROOT=Path('/workspace/hdblast-research-work')
OUT=Path(__file__).resolve().parent
SPECS={
 'P':('common-gauge-forcing-20261016/capsule-source','8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6','common-gauge-forcing-20261016/publication/hdblast-worktree/hdblast/continuation/2026.10.16-common-gauge-forcing/PUBLICATION_PINS.json'),
 'N':('numerical-response-20261016/capsule-source','928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44','numerical-response-20261016/publication/research-worktree/hdblast/continuation/2026.10.16-numerical-response-inputs/PUBLICATION_PINS.json'),
 'G':('geometric-bulk-source-20261016/capsule-source','3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2','boundary-response-20261023/followon/NEXT_SCIENTIFIC_TARGET.md'),
 'C':('forced-response-20261023/package-work/attempt-2/source','72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76','boundary-response-20261023/followon/NEXT_SCIENTIFIC_TARGET.md')}
def sha(raw):return hashlib.sha256(raw).hexdigest()
def run(dependency_map=None):
 report={'status':'PASS_COMPLETE_SELECTED_CAPSULE_CLOSURES','closures':{}}
 for label,(rel,expected,external) in SPECS.items():
  root=Path(dependency_map[label]).resolve() if dependency_map is not None else ROOT/rel;raw=(root/'MANIFEST.json').read_bytes()
  if sha(raw)!=expected:raise ValueError('manifest pin '+label)
  if dependency_map is None:
   external_source=str(ROOT/external);eraw=(ROOT/external).read_bytes()
   if expected not in eraw.decode():raise ValueError('external manifest pin '+label)
  else:
   external_source='embedded independently retained manifest SHA256 in authenticate_inputs.py';eraw=Path(__file__).read_bytes()
   if expected not in eraw.decode():raise ValueError('embedded external manifest pin '+label)
  mf=json.loads(raw);seen=set();nbytes=0;entries=[]
  for e in mf['files']:
   relpath=e['path'];p=Path(relpath)
   if p.is_absolute() or '..' in p.parts or relpath in seen:raise ValueError('unsafe/duplicate path')
   seen.add(relpath);f=root/p
   if f.is_symlink():raise ValueError('unexpected symbolic input')
   content=f.read_bytes()
   if len(content)!=e['bytes'] or sha(content)!=e['sha256']:raise ValueError('closure mismatch '+str(f))
   nbytes+=len(content);entries.append(e)
  undeclared=sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and str(p.relative_to(root)) not in seen|{'MANIFEST.json'} and '__pycache__' not in p.parts)
  if undeclared:raise ValueError('undeclared closure files '+str(undeclared))
  report['closures'][label]={'root':str(root),'manifest_sha256':expected,'external_pin_source':external_source,'external_pin_source_sha256':sha(eraw),'files_verified':len(entries),'declared_bytes_verified':nbytes,'undeclared_regular_files':undeclared,'entries':entries}
 return report
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
 r=run();args.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:{'files':v['files_verified'],'bytes':v['declared_bytes_verified']} for k,v in r['closures'].items()}))
