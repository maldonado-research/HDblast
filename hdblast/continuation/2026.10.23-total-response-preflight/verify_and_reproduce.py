#!/usr/bin/env python3
"""Authenticate exact source/ZIP closure, then run only the portable replay."""
import argparse,hashlib,json,subprocess,sys,zipfile
from pathlib import Path
def need(ok,label):
 if not ok:raise ValueError(label)
def sha(b):return hashlib.sha256(b).hexdigest()
p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--pins',type=Path,required=True);p.add_argument('--pins-sha256',required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
raw=a.pins.read_bytes();need(sha(raw)==a.pins_sha256,'external publication pins');pins=json.loads(raw)
repo=a.repo_root.resolve();source=repo/pins['source_path'];archive=repo/pins['archive_path'];out=a.output_dir.resolve()
need(not out.exists() and repo!=out and repo not in out.parents,'fresh output outside repository');out.mkdir(parents=True)
for rel,pin in pins['bootstrap_files'].items():
 f=repo/rel;need(not f.is_symlink() and sha(f.read_bytes())==pin,'bootstrap file hash')
need(sha((source/'MANIFEST.json').read_bytes())==pins['manifest_sha256'],'manifest')
need(sha(archive.read_bytes())==pins['archive_sha256'],'archive')
mf=json.loads((source/'MANIFEST.json').read_bytes());entries={r['path']:r for r in mf['files']};need(len(entries)==len(mf['files']),'unique manifest')
expected=set(entries)|{'MANIFEST.json'};actual={str(f.relative_to(source)) for f in source.rglob('*') if f.is_file()};need(actual==expected,'complete source inventory')
for n in expected:
 f=source/n;need(not f.is_symlink(),'no source symlink')
 if n!='MANIFEST.json':need(len(f.read_bytes())==entries[n]['bytes'] and sha(f.read_bytes())==entries[n]['sha256'],'source member')
with zipfile.ZipFile(archive) as z:
 names=z.namelist();wanted={pins['capsule']+'/'+n for n in expected};need(len(names)==len(wanted) and set(names)==wanted,'exact ZIP members')
 for n in expected:need(z.read(pins['capsule']+'/'+n)==(source/n).read_bytes(),'ZIP/source bytes')
subprocess.run([sys.executable,'-I','-S','-B',str(source/'reproduce.py'),'--manifest-sha256',pins['manifest_sha256'],'--output-dir',str(out/'portable')],check=True)
(out/'AUTHENTICATION.json').write_text(json.dumps({'status':'PASS_SOURCE_ZIP_BOOTSTRAP_AUTHENTICATION','capsule':pins['capsule'],'source_files':len(expected),'manifest_sha256':pins['manifest_sha256'],'archive_sha256':pins['archive_sha256']},indent=2)+'\n')
print('PASS_SOURCE_ZIP_AND_PORTABLE_REPLAY')
