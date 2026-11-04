#!/usr/bin/env python3
"""Run isolated independent normal/O reconstructions and ten strict CLI controls."""
import argparse,ast,hashlib,json
from pathlib import Path
import shutil,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
RECOMPUTE_SHA='6e730528b134423efc965464f828d68c94ef06b65d44f50afea13c84d4b617ad'
class Rejected(Exception):pass
class Parser(argparse.ArgumentParser):
 def error(self,message):raise Rejected('ARGUMENTS')
def need(x,w):
 if not x:raise Rejected(w)
def canonical(v):return (json.dumps(v,sort_keys=True,indent=2)+'\n').encode()
def sha(raw):return hashlib.sha256(raw).hexdigest()
def safe(p):return p.is_absolute() and not any(q.is_symlink() for q in [p,*p.parents])
def main():
 ap=Parser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--candidate',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
 need(safe(a.repo_root) and a.repo_root.is_dir(),'REPO_PATH');repo=a.repo_root.resolve();need(safe(a.candidate) and a.candidate.is_file(),'CANDIDATE_PATH');candidate=a.candidate.resolve()
 need(safe(a.output) and a.output.parent.is_dir(),'OUTPUT_PATH');out=a.output.resolve();artifacts=out.with_suffix('.artifacts');need(not out.exists() and not artifacts.exists(),'OUTPUT_EXISTS');need(all(x not in out.parents and x!=out for x in [repo,HERE]),'OUTPUT_LOCATION')
 source=HERE/'recompute.py';need(source.is_file() and not source.is_symlink() and sha(source.read_bytes())==RECOMPUTE_SHA,'RECOMPUTE_PIN')
 # Parse the fixed literal source roster; never import or execute a project module.
 tree=ast.parse(source.read_text());pins=None
 for node in tree.body:
  if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PINS' for t in node.targets):pins=ast.literal_eval(node.value)
 need(isinstance(pins,dict) and set(pins)=={'premises','static','base','witness'},'SOURCE_ROSTER')
 for alias,(rel,pin) in pins.items():
  p=repo/rel;need(p.is_file() and not p.is_symlink() and sha(p.read_bytes())==pin,'INPUT_PIN:'+alias)
 artifacts.mkdir();positive=[];negative=[];modes=[]
 with tempfile.TemporaryDirectory(prefix='hdblast-independent-h1-controls-') as tmp:
  fixture=Path(tmp)/'bad-repo';fixture.mkdir()
  for rel,pin in pins.values():
   dest=fixture/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(repo/rel,dest)
  badsource=fixture/pins['premises'][0];before=badsource.read_text();need('1/10000' in before,'SOURCE_MUTATION_TARGET');badsource.write_text(before.replace('1/10000','9/10000',1))
  badprimary=Path(tmp)/'bad-primary.json';wrong=json.loads(candidate.read_bytes());wrong['class']['radius_closed']=str(2*10**31+1);badprimary.write_bytes(canonical(wrong))
  unknownprimary=Path(tmp)/'unknown-primary.json';wrong=json.loads(candidate.read_bytes());wrong['unregistered_extra_field']=True;unknownprimary.write_bytes(canonical(wrong))
  badfraction=Path(tmp)/'bad-fraction.json';wrong=json.loads(candidate.read_bytes());wrong['energy']['initial_energy_per_R_squared']='1/0';badfraction.write_bytes(canonical(wrong))
  def invoke(opt,rep,prim,dest):
   flags=['-I','-S','-B']+(['-O'] if opt else []);cmd=[sys.executable,*flags,str(source),'--repo-root',str(rep),'--candidate',str(prim),'--output',str(dest)]
   proc=subprocess.run(cmd,capture_output=True,text=True,timeout=60)
   need(not proc.stderr and len(proc.stdout.splitlines())==1,'CHILD_EXACT_PROTOCOL')
   try:value=json.loads(proc.stdout)
   except ValueError:raise Rejected('CHILD_EXACT_PROTOCOL') from None
   return proc,value,flags
  for opt,label in [(False,'NORMAL'),(True,'OPTIMIZED')]:
   dest=artifacts/(label+'.json');proc,value,flags=invoke(opt,repo,candidate,dest)
   need(proc.returncode==0 and dest.is_file(),'POSITIVE_CHILD')
   payload=dest.read_bytes();need(value=={'status':'PASS','witness':'INDEPENDENT_H1_RECONSTRUCTION','report_sha256':sha(payload)},'POSITIVE_WITNESS')
   positive.append({'mode':label,'flags':flags,'returncode':0,'stdout':value,'stderr_empty':True,'artifact':dest.name,'artifact_sha256':sha(payload)});modes.append(payload)
   cases=[('source-pin',fixture,candidate,'INPUT_PIN:premises',False),('shared-primary-radius',repo,badprimary,'PRIMARY_COMPARISON',False),('unknown-primary-field',repo,unknownprimary,'PRIMARY_SCHEMA',False),('invalid-primary-fraction',repo,badfraction,'INVALID_INPUT_STRUCTURE',False),('occupied-output',repo,candidate,'OUTPUT_EXISTS',True)]
   for case,rep,prim,witness,occupied in cases:
    rejected=Path(tmp)/(label+'-'+case+'.json');sentinel=b'PRESERVE_EXISTING_OUTPUT\n'
    if occupied:rejected.write_bytes(sentinel)
    proc,value,flags=invoke(opt,rep,prim,rejected)
    need(proc.returncode==2 and value=={'status':'FAIL','witness':witness},'NEGATIVE_EXACT_WITNESS')
    if occupied:need(rejected.read_bytes()==sentinel,'EXISTING_OUTPUT_PRESERVED')
    else:need(not rejected.exists(),'REJECTED_OUTPUT_ABSENT')
    negative.append({'case':case,'mode':label,'flags':flags,'returncode':2,'stdout':value,'stderr_empty':True,'output_absent':not occupied,'existing_output_preserved':occupied})
 need(modes[0]==modes[1],'NORMAL_OPTIMIZED_CANONICAL_IDENTITY')
 report={'schema':'HDBLAST_INDEPENDENT_H1_STRICT_CONTROLS_V1','status':'PASS_INDEPENDENT_H1_STRICT_CONTROLS','positive_children':positive,'negative_children':negative,'canonical_normal_optimized_byte_identical':True,'actual_positive_subprocesses':len(positive),'actual_rejected_subprocesses':len(negative),'inprocess_mutant_executions':0,'activity':{'PDE_solves':0,'background_or_evolved_field_queries':0,'bump_queries':0,'analytic_target_bump_replays':0},'scope':'One final wrapper execution runs two canonical mathematical children and ten exact rejection children. Source authentication, not inherited PDE-proof replay.'}
 raw=canonical(report);out.write_bytes(raw);print(json.dumps({'status':'PASS','witness':'INDEPENDENT_H1_STRICT_CONTROLS','positive_children':len(positive),'rejected_children':len(negative),'report_sha256':sha(raw)},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Rejected as e:print(json.dumps({'status':'FAIL','witness':str(e)},sort_keys=True));sys.exit(2)
 except (OSError,ValueError,KeyError,TypeError,ArithmeticError) as e:print(json.dumps({'status':'FAIL','witness':'INVALID_INPUT_STRUCTURE'},sort_keys=True));sys.exit(2)
