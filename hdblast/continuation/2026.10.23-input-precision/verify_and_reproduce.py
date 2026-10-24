#!/usr/bin/env python3
"""Authenticate exact source and ZIP before any selected numerical replay."""
import argparse,hashlib,json,subprocess,sys,zipfile
from pathlib import Path
def need(ok,why):
 if not ok:raise ValueError(why)
def sha(b):return hashlib.sha256(b).hexdigest()
ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--pins-sha256',required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--dependency-map',type=Path);a=ap.parse_args()
raw=a.pins.read_bytes();need(sha(raw)==a.pins_sha256,'independent publication pins');pins=json.loads(raw);repo=a.repo_root.resolve();source=repo/pins['source_path'];archive=repo/pins['archive_path'];out=a.output_dir.resolve();need(not out.exists() and repo!=out and repo not in out.parents,'fresh output outside repo');out.mkdir(parents=True)
for rel,pin in pins['bootstrap_files'].items():
 p=repo/rel;need(not p.is_symlink() and sha(p.read_bytes())==pin,'bootstrap member')
need(sha((source/'MANIFEST.json').read_bytes())==pins['manifest_sha256'],'source manifest');need(sha(archive.read_bytes())==pins['archive_sha256'],'exact ZIP')
m=json.loads((source/'MANIFEST.json').read_text());entries={e['path']:e for e in m['files']};expected=set(entries)|{'MANIFEST.json'};need(len(entries)==len(m['files']),'unique manifest')
actual=set()
for p in source.rglob('*'):
 need(not p.is_symlink(),'no symbolic source path')
 if p.is_file():actual.add(str(p.relative_to(source)))
need(actual==expected,'complete source inventory')
for n,e in entries.items():
 p=Path(n);need(not p.is_absolute() and '..' not in p.parts,'safe member');b=(source/n).read_bytes();need(len(b)==e['bytes'] and sha(b)==e['sha256'],'source member')
with zipfile.ZipFile(archive) as z:
 wanted={pins['capsule']+'/'+n for n in expected};need(len(z.namelist())==len(wanted) and set(z.namelist())==wanted,'complete exact ZIP inventory')
 for n in expected:need(z.read(pins['capsule']+'/'+n)==(source/n).read_bytes(),'ZIP/source agreement')
cmd=[sys.executable,'-I','-S','-B',str(source/'reproduce.py'),'--manifest-sha256',pins['manifest_sha256'],'--output-dir',str(out/'replay')]
if a.dependency_map:cmd+=['--native-dependency-map',str(a.dependency_map.resolve())]
subprocess.run(cmd,check=True)
(out/'AUTHENTICATION.json').write_text(json.dumps({'status':'PASS_EXACT_SOURCE_ZIP_BOOTSTRAP_AUTHENTICATION','capsule':pins['capsule'],'source_files':len(expected),'manifest_sha256':pins['manifest_sha256'],'archive_sha256':pins['archive_sha256'],'native_requested':bool(a.dependency_map)},indent=2)+chr(10))
print('PASS_EXACT_SOURCE_ZIP_AND_SELECTED_REPLAY')
