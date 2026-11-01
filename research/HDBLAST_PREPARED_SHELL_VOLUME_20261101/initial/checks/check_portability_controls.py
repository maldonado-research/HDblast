#!/usr/bin/env python3
"""Bounded import/authentication controls, with optional full relocated replay."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
CAP=Path('research/HDBLAST_PREPARED_SHELL_VOLUME_20261101')

def need(ok,msg):
 if not ok:raise ValueError(msg)

def copy_source(source,dest):
 for rel in json.loads((ROOT/'SOURCE_PINS.json').read_text()):
  p=dest/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/rel,p)

def invoke(script,repo,gauge,out,optimized):
 p=subprocess.run([sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(script),'--repo-root',str(repo),'--gauge-api',str(gauge),'--output-dir',str(out)],text=True,capture_output=True)
 lines=p.stdout.splitlines();need(len(lines)==1 and bool(lines[0]),'stdout must be exactly one JSON line');need(p.stderr=='','stderr must be empty')
 witness=json.loads(lines[0]);need(isinstance(witness,dict) and set(witness)=={'FINAL'},'not an exact FINAL witness')
 return p,witness['FINAL']

def main():
 a=argparse.ArgumentParser();a.add_argument('--repo-root',type=Path,required=True);a.add_argument('--gauge-api',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a.add_argument('--relocated-positive-output-dir',type=Path);v=a.parse_args();need(not v.output.exists(),'refuse overwrite')
 cases=[]
 with tempfile.TemporaryDirectory(prefix='hdblast-initial-portable-controls-') as td:
  tmp=Path(td);repo=tmp/'repo';copy_source(v.repo_root,repo);gauge=tmp/'gauge';gauge.mkdir()
  for name in ('mean_initial.py','conditioned_mean_initial.py'):shutil.copy2(v.gauge_api/name,gauge/name)
  for optimized in (False,True):
   for kind,expected in [('source','SOURCE_AUTHENTICATION'),('gauge','GAUGE_AUTHENTICATION'),('local','LOCAL_AUTHENTICATION'),('reference','REFERENCE_AUTHENTICATION'),('existing','OUTPUT_ALREADY_EXISTS'),('inside','OUTPUT_INSIDE_SOURCE')]:
    case=tmp/(kind+str(optimized));case.mkdir();script=ROOT/'reproduce.py';r=repo;g=gauge;out=case/'result';changed=None
    if kind=='source':
     r=case/'repo';copy_source(v.repo_root,r);changed=r/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/interval_core.py';changed.write_text(changed.read_text()+'\nraise RuntimeError("scientific import must be blocked")\n')
    if kind=='gauge':
     g=case/'gauge';shutil.copytree(gauge,g);changed=g/'conditioned_mean_initial.py';changed.write_text(changed.read_text()+'\nraise RuntimeError("gauge import must be blocked")\n')
    if kind in ('local','reference'):
     pkg=case/'initial';shutil.copytree(ROOT,pkg,ignore=shutil.ignore_patterns('__pycache__'));script=pkg/'reproduce.py';changed=pkg/('scientific.py' if kind=='local' else 'references/RAW_BOXES.json');changed.write_text(changed.read_text()+'\n')
     if kind=='reference':
      # Satisfy the outer local-byte layer, leaving the independently pinned
      # scientific-reference layer to reject before importing science.
      pins=json.loads((pkg/'LOCAL_PINS.json').read_text());raw=changed.read_bytes();pins['references/RAW_BOXES.json']={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()};(pkg/'LOCAL_PINS.json').write_text(json.dumps(pins,sort_keys=True,indent=2)+'\n')
    if kind=='existing':out.mkdir()
    if kind=='inside':out=repo/'must-not-be-created'
    p,witness=invoke(script,r,g,out,optimized);need(p.returncode==2,'wrong rejection exit');need(witness=={'status':'REJECTED','witness':expected},'wrong exact FINAL witness');need((out.exists() and not any(out.iterdir())) if kind=='existing' else not out.exists(),'rejection created scientific output')
    cases.append({'case':kind,'optimized':optimized,'child_flags':['-I','-S','-B']+(['-O'] if optimized else []),'exit_code':p.returncode,'stdout_exact':p.stdout,'stderr_exact':p.stderr,'stdout_one_json':True,'stderr_empty':True,'FINAL':witness,'new_scientific_output_absent':True,'actual_initial_background_calls_executed':0,'evolved_actual_field_oracle_calls_executed':0})
  positives=[]
  if v.relocated_positive_output_dir:
   out=v.relocated_positive_output_dir.resolve();need(not out.exists(),'positive relocation output exists');pkg=repo/CAP/'initial';shutil.copytree(ROOT,pkg,ignore=shutil.ignore_patterns('__pycache__'));gg=repo/CAP/'gauge';shutil.copytree(gauge,gg)
   p,w=invoke(pkg/'reproduce.py',repo,gg,out,False);expected={'status':'PASS','witness':'FROZEN_INITIAL_SCIENCE_MATCH','initial_background_calls':1990,'evolved_actual_field_calls':0,'scientific_outputs':['RAW_BOXES.json','CONDITIONED_RAW_BOXES.json']};need(p.returncode==0 and w==expected,'relocated positive replay');positives.append({'case':'fresh_source_code_and_gauge_relocation','optimized':False,'child_flags':['-I','-S','-B'],'exit_code':0,'stdout_exact':p.stdout,'stderr_exact':p.stderr,'stdout_one_json':True,'stderr_empty':True,'FINAL':w,'receipt_sha256':hashlib.sha256((out/'RECEIPT.json').read_bytes()).hexdigest()})
 report={'status':'PASS','bounded_negative_cases':len(cases),'fresh_relocated_positive_generations':len(positives),'negative_controls':cases,'relocated_positive':positives,'counting':'Each negative exits before scientific imports and performs zero initial-background/evolved-field oracle calls. Each full positive generation performs 1990 actual-initial-background calls and zero evolved actual-field oracle calls. Positive relocation is counted separately from negative controls.'}
 v.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'FINAL':{'status':'PASS','negative_cases':len(cases),'relocated_positive':len(positives)}}))
if __name__=='__main__':main()
