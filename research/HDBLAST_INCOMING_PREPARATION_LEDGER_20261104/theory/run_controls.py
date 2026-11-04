#!/usr/bin/env python3
"""One canonical pair plus exact real CLI rejection controls; stdlib only."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
class Failure(Exception):pass
def need(v,w):
 if not v:raise Failure(w)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()
def regular(p):return p.is_file() and not p.is_symlink() and all(not q.is_symlink()for q in p.absolute().parents)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();repo=a.repo_root.absolute();out=a.output.absolute()
 need(all(not p.is_symlink()for p in [repo,*repo.parents,out,*out.parents]),'CONTROL_PATH_SYMLINK');repo=repo.resolve(strict=True);out=out.resolve();need(not out.exists(),'CONTROL_FRESH_OUTPUT');need(out.parent.is_dir(),'CONTROL_OUTPUT_PARENT');need(all(b!=out and b not in out.parents for b in [repo,HERE]),'CONTROL_OUTPUT_OUTSIDE_SOURCE')
 artifacts=out.with_suffix('.artifacts');need(not artifacts.exists(),'CONTROL_FRESH_ARTIFACTS');artifacts.mkdir();sourcepins=json.loads((HERE/'SOURCE_PINS.json').read_bytes());rows=[];positives=[]
 cases={'source_bytes':'AUTHENTICATED_SOURCE_BYTES','output_exists':'FRESH_OUTPUT_FILE','output_inside_repo':'OUTPUT_OUTSIDE_CODE_AND_SOURCE','design_bytes':'FROZEN_DESIGN_BYTES'}
 def run(repo_path,target,opt,design=None):
  flags=['-I','-S','-B',*(['-O']if opt else [])];cmd=[sys.executable,*flags,str(HERE/'verify.py'),'--repo-root',str(repo_path),'--output',str(target)]
  if design is not None:cmd+=['--design',str(design)]
  r=subprocess.run(cmd,capture_output=True,text=True,timeout=60);need(not r.stderr and len(r.stdout.splitlines())==1,'EXACT_CHILD_PROTOCOL');return r,json.loads(r.stdout),flags
 for opt in [False,True]:
  target=artifacts/('OPTIMIZED.json'if opt else 'NORMAL.json');r,document,flags=run(repo,target,opt);need(r.returncode==0 and document=={'status':'PASS','witness':'INCOMING_PREPARATION_THEORY_IDENTITIES','finite_checks':32,'output_sha256':sha(target)},'CANONICAL_CHILD_WITNESS');report=json.loads(target.read_bytes());need(report['status']=='PASS_CONDITIONAL_INCOMING_PREPARATION_IDENTITIES','CANONICAL_REPORT_STATUS');positives.append({'optimized':opt,'flags':flags,'returncode':r.returncode,'stdout':document,'stderr_empty':True,'report_sha256':sha(target)})
  for case,reason in cases.items():
   with tempfile.TemporaryDirectory(prefix='hdblast-incoming-theory-control-',dir=out.parent)as td:
    tmp=Path(td);fixture=tmp/'repo';fixture.mkdir()
    for name in sourcepins['files']:
     p=fixture/name;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(repo/name,p)
    output=tmp/'result.json';design=None;before=None
    if case=='source_bytes':
     p=fixture/next(iter(sourcepins['files']));p.write_bytes(p.read_bytes()+b'\n')
    elif case=='output_exists':output.write_bytes(b'preexisting-bytes-must-remain\n');before=output.read_bytes()
    elif case=='output_inside_repo':output=fixture/'never-written.json'
    elif case=='design_bytes':design=tmp/'wrong-design.json';design.write_bytes((HERE/'FROZEN_DESIGN.json').read_bytes()+b'\n')
    r,document,flags=run(fixture,output,opt,design);expected={'status':'FAIL','witness':reason};need(r.returncode==2 and document==expected,'EXACT_REJECTION_WITNESS');preserved=output.read_bytes()==before if before is not None else not output.exists();need(preserved,'REJECTED_OUTPUT_ARTIFACT_CONTRACT');rows.append({'case':case,'optimized':opt,'flags':flags,'returncode':r.returncode,'stdout':document,'stderr_empty':True,'output_absent':before is None and not output.exists(),'preexisting_bytes_preserved':before is not None and preserved})
 need(positives[0]['report_sha256']==positives[1]['report_sha256'],'CANONICAL_NORMAL_OPTIMIZED_BYTES')
 report={'status':'PASS_STRICT_INCOMING_THEORY_CONTROLS','actual_rejected_subprocesses':len(rows),'distinct_rejection_cases':len(cases),'positive_subprocesses':len(positives),'negative_cases':rows,'canonical_children':positives,'canonical_pair_sha256':positives[0]['report_sha256'],'new_PDE_solves':0,'new_background_or_evolved_field_queries':0,'new_bump_queries':0,'scope':'Each negative is a real isolated child; source authentication/path/output protocol controls do not prove analytic PDE premises.'}
 with out.open('xb')as f:f.write(canonical(report))
 print(json.dumps({'status':'PASS','witness':'STRICT_INCOMING_THEORY_CONTROLS','actual_rejected_subprocesses':len(rows),'positive_subprocesses':len(positives),'output_sha256':sha(out)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except (Failure,ArithmeticError,ValueError,OSError,KeyError,TypeError)as e:
  print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));sys.exit(2)
