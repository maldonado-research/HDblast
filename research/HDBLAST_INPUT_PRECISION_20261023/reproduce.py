#!/usr/bin/env python3
"""Authenticate this complete capsule, then replay its declared finite roles."""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path
def need(ok,why):
 if not ok:raise ValueError(why)
def sha(b):return hashlib.sha256(b).hexdigest()
def unique(items):
 d={}
 for k,v in items:
  need(k not in d,'duplicate JSON key');d[k]=v
 return d
def authenticate(root,pin):
 raw=(root/'MANIFEST.json').read_bytes();need(sha(raw)==pin,'independent manifest pin')
 m=json.loads(raw,object_pairs_hook=unique);expected={'MANIFEST.json'}
 for e in m['files']:
  n=e['path'];q=Path(n);need(n not in expected and not q.is_absolute() and '..' not in q.parts,'safe unique member')
  expected.add(n);p=root/q;need(p.is_file() and not p.is_symlink(),'regular source member')
  b=p.read_bytes();need(len(b)==e['bytes'] and sha(b)==e['sha256'],'source member '+n)
 actual=set()
 for p in root.rglob('*'):
  need(not p.is_symlink(),'no symbolic source path')
  if p.is_file():actual.add(str(p.relative_to(root)))
 need(actual==expected,'exact complete source inventory');return len(expected)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--manifest-sha256',required=True);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--native-dependency-map',type=Path);a=ap.parse_args()
 root=Path(__file__).resolve().parent;out=a.output_dir.resolve();need(root!=out and root not in out.parents and not out.exists(),'fresh output outside source');out.mkdir(parents=True)
 count=authenticate(root,a.manifest_sha256);jobs=[]
 def run(role,script,args,mode='normal',golden=None):
  target=out/(role+'_'+mode+'.json');cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+[str(root/script),*args,'--output',str(target)]
  r=subprocess.run(cmd,cwd=out,capture_output=True,text=True);need(r.returncode==0,role+' '+r.stderr[-1500:]);need(target.is_file(),'child receipt')
  if golden:need(target.read_bytes()==(root/golden).read_bytes(),'exact golden '+role)
  jobs.append({'role':role,'mode':mode,'sha256':sha(target.read_bytes()),'frozen_golden_compared':bool(golden)})
  return target
 for mode in ('normal','optimized'):
  run('branch','branch/portable/replay_branch.py',[],mode,'branch/portable/PORTABLE_NORMAL.json')
  run('initial-ledger','initial/check_difference_ledger.py',['--report',str(root/'initial/RECENTERED_CANDIDATE_NORMAL.json')],mode,'initial/RECENTERED_LEDGER.json')
  run('fixed-candidate-ledger','initial/check_difference_ledger.py',['--report',str(root/'initial/FIXED_CANDIDATE_VALIDATED_BETA.json')],mode,'initial/FIXED_CANDIDATE_LEDGER.json')
  run('independent-initial-ledger','review/verify_initial_ledger.py',['--report',str(root/'initial/RECENTERED_CANDIDATE_NORMAL.json'),'--branch',str(root/'branch/BRANCH_COMPARISON_REPORT.json')],mode,'review/INDEPENDENT_INITIAL_LEDGER_NORMAL.json')
  run('model-identities','review/check_independent_model_identities.py',[],mode,'review/INDEPENDENT_MODEL_IDENTITIES_NORMAL.json')
  # The independent spot control uses no project inputs and is copied into
  # a writable directory; the frozen source remains read-only throughout.
  work=out/('paired-control-'+mode);work.mkdir();script=work/'check_difference_controls.py';shutil.copyfile(root/'initial/check_difference_controls.py',script)
  cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+[str(script)]
  r=subprocess.run(cmd,cwd=work,capture_output=True,text=True);need(r.returncode==0,'paired algebra '+r.stderr[-1000:]);target=work/'DIFFERENCE_CONTROLS.json'
  need(target.read_bytes()==(root/'initial/DIFFERENCE_CONTROLS.json').read_bytes(),'paired control golden')
  jobs.append({'role':'paired-algebra-and-kernel-spots','mode':mode,'sha256':sha(target.read_bytes()),'frozen_golden_compared':True})
 run('branch-relocation-controls','branch/portable/verify_portable.py',[],golden='branch/portable/PORTABLE_CHECKS.json')
 native=[]
 if a.native_dependency_map:
  # Both modes execute the full 96-cell formulas against five complete,
  # fixed-hash upstream capsules, not only their saved ledger reports.
  processes=[]
  for mode in ('normal','optimized'):
   target=out/('native_'+mode+'.json');log=out/('native_'+mode+'.log');stream=log.open('w')
   cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='optimized' else [])+[str(root/'initial/direct_difference_pilot.py'),'--case','branch','--dependency-map',str(a.native_dependency_map.resolve()),'--recenter','--output',str(target)]
   processes.append((mode,target,stream,subprocess.Popen(cmd,cwd=out,stdout=stream,stderr=subprocess.STDOUT)))
  for mode,target,stream,process in processes:
   code=process.wait();stream.close();need(code==0,'native '+mode+' failed; inspect external log');native.append({'mode':mode,'sha256':sha(target.read_bytes())})
  run('native-comparison','initial/compare_difference_reports.py',[str(root/'initial/RECENTERED_CANDIDATE_NORMAL.json'),str(out/'native_normal.json'),str(out/'native_optimized.json')])
  run('native-ledger','initial/check_difference_ledger.py',['--report',str(out/'native_normal.json')],golden='initial/RECENTERED_LEDGER.json')
  run('native-authentication-controls','initial/check_authentication_controls.py',['--dependency-map',str(a.native_dependency_map.resolve()),'--scratch-dir',str(out/'owned-authentication-copies')])
 need(authenticate(root,a.manifest_sha256)==count,'complete post-replay authentication')
 receipt={'status':'PASS_COMPLETE_INPUT_PRECISION_REPLAY','manifest_sha256':a.manifest_sha256,'source_files_including_manifest':count,'portable_and_ledger_roles':jobs,'native_interval_evaluations':native,'scope':'Finite exact branch comparison with inherited theorem premises; complete initial report/coverage arithmetic and controls. Native entries, when present, additionally reexecute all96 interval cells with full dependency authentication. No bulk iteration, total-response waveform, invariant observable sign or nonlinear family is evaluated.'}
 (out/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'roles':len(jobs),'native_modes':len(native),'source_files':count}))
if __name__=='__main__':main()
