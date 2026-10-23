#!/usr/bin/env python3
"""Authenticate the selected capsule and replay its portable checks only."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path

def need(ok,label):
 if not ok:raise ValueError(label)
def unique(rows):
 out={}
 for k,v in rows:
  need(k not in out,'duplicate JSON key');out[k]=v
 return out
def sha(raw):return hashlib.sha256(raw).hexdigest()
def auth(root,expected):
 raw=(root/'MANIFEST.json').read_bytes();need(sha(raw)==expected,'external manifest hash')
 mf=json.loads(raw,object_pairs_hook=unique);seen={'MANIFEST.json'}
 for row in mf['files']:
  rel=row['path'];p=Path(rel);need(not p.is_absolute() and '..' not in p.parts and rel not in seen,'safe unique member')
  seen.add(rel);f=root/p;need(f.is_file() and not f.is_symlink(),'regular member')
  b=f.read_bytes();need(len(b)==row['bytes'] and sha(b)==row['sha256'],'member mismatch: '+rel)
 actual=set()
 for f in root.rglob('*'):
  need(not f.is_symlink(),'no symlink');
  if f.is_file():actual.add(str(f.relative_to(root)))
 need(actual==seen,'exact declared source closure');return len(seen)
def main():
 p=argparse.ArgumentParser();p.add_argument('--manifest-sha256',required=True);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 root=Path(__file__).resolve().parent;out=a.output_dir.resolve();need(out!=root and root not in out.parents,'output outside frozen source');need(not out.exists(),'fresh output');out.mkdir(parents=True)
 count=auth(root,a.manifest_sha256)
 jobs=[('conditional-arithmetic','analytic/check_portable_local_total.py',[],'analytic/PORTABLE_LOCAL_TOTAL_NORMAL.json'),
 ('initial-formula','initial/check_initial_formula_identities.py',[],'initial/INITIAL_FORMULA_REVIEW_NORMAL.json'),
 ('initial-ledger','initial/check_initial_gate_ledger.py',['--report',str(root/'initial/INITIAL_GATE_NORMAL.json')],'initial/LEDGER_NORMAL.json'),
 ('initial-ledger-controls','initial/run_gate_controls.py',['--report',str(root/'initial/INITIAL_GATE_NORMAL.json')],'initial/CONTROLS_NORMAL.json'),
 ('observable-algebra','observable/check_observable_algebra.py',[],'observable/ALGEBRA_REVIEWED_NORMAL.json')]
 results=[]
 for mode in ('normal','optimized'):
  for role,script,args,golden in jobs:
   target=out/(role+'_'+mode+'.json')
   cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+[str(root/script),*args,'--output',str(target)]
   r=subprocess.run(cmd,cwd=out,capture_output=True,text=True);need(r.returncode==0,role+' failed '+r.stderr[-1000:])
   need(target.read_bytes()==(root/golden).read_bytes(),role+' golden mismatch');results.append({'role':role,'mode':mode,'sha256':sha(target.read_bytes()),'matches_frozen_golden':True})
 need((root/'initial/INITIAL_GATE_NORMAL.json').read_bytes()==(root/'initial/INITIAL_GATE_OPTIMIZED.json').read_bytes(),'native pilot normal/O reports differ')
 need(auth(root,a.manifest_sha256)==count,'postclosure')
 receipt={'status':'PASS_FRESH_PORTABLE_REPLAY','manifest_sha256':a.manifest_sha256,'source_files_including_manifest':count,'children':results,'scope':'Portable exact arithmetic, symbolic polynomial identities and initial-report ledger only. The optional native initial producer has separate complete upstream closure authentication; upstream PDE premises and large source campaigns are not replayed here. No total response or observable value is evaluated.'}
 (out/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'children':len(results),'source_files':count}))
if __name__=='__main__':main()
