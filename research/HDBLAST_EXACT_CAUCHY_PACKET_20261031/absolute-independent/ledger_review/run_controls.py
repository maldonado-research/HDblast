#!/usr/bin/env python3
"""Real subprocess controls for the separate exact Fraction review."""
import argparse,copy,hashlib,json,shutil,subprocess,sys
from pathlib import Path

ASSETS=Path(__file__).resolve().parent
BASE=None
WORK=None
SOURCE=None
OWN=None
RESPONSE=None
MAP_SHA='255d3c3372d97a6a5a7ffa5663fe55363c904129adf62102d781014df4124272'
ABS_SHA='b7b0511e899b71b022d56e02dcbef0d116cf3fede9fc4e1c7bf7f73c77514275'
RESP_SHA='6a0adb48d79b210a60322920e70ad6a60f2862b1ade3222c67a285419198b74e'
RESP_SRC_SHA='b4d6cd3651cf71a987bb1d32d2a746eb3bc3c64a88f245caf8d956c433e2a417'

def need(x,why):
 if not x: raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
def invoke(checker,root,pins,absolute,support,response,response_source,out,optimized=False,**over):
 argv=[sys.executable]+(['-O'] if optimized else [])+[str(checker),
 '--repo-root',str(root),'--source-pins',str(pins),'--source-pins-sha256',over.get('mapsha',MAP_SHA),
 '--absolute',str(absolute),'--absolute-sha256',over.get('abssha',ABS_SHA),
 '--absolute-support-root',str(support),'--producer-response',str(response),
 '--producer-response-sha256',over.get('respsha',RESP_SHA),'--producer-source',str(response_source),
 '--producer-source-sha256',RESP_SRC_SHA,'--output',str(out)]
 return subprocess.run(argv,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)

def main():
 global BASE,WORK,SOURCE,OWN,RESPONSE
 ap=argparse.ArgumentParser()
 ap.add_argument('--source-root',required=True,type=Path)
 ap.add_argument('--absolute-root',required=True,type=Path)
 ap.add_argument('--response-root',required=True,type=Path)
 ap.add_argument('--output-dir',required=True,type=Path)
 a=ap.parse_args()
 SOURCE,OWN,RESPONSE=a.source_root,a.absolute_root,a.response_root
 BASE=a.output_dir;WORK=BASE/'controls'
 need(not BASE.exists(),'fresh control output directory')
 BASE.mkdir(parents=True)
 need(not WORK.exists(),'fresh controls directory')
 WORK.mkdir(parents=True)
 checker=ASSETS/'check_refined_response_independent.py'
 before={checker:sha(checker),ASSETS/'source-pins.json':sha(ASSETS/'source-pins.json')}
 paths=(checker,SOURCE,ASSETS/'source-pins.json',OWN/'reference/INDEPENDENT_ABSOLUTE.json',OWN,
 RESPONSE/'reference/RESPONSE.json',RESPONSE/'checks/response_adapter.py')
 outputs=[]
 for opt in (False,True):
  out=BASE/('INDEPENDENT_RESPONSE_OPTIMIZED.json' if opt else 'INDEPENDENT_RESPONSE.json')
  need(not out.exists(),'fresh normal optimized result')
  got=invoke(*paths,out,opt)
  need(got.returncode==0,got.stderr)
  outputs.append(out)
 need(outputs[0].read_bytes()==outputs[1].read_bytes(),'normal optimized canonical equality')
 # A real sparse relocation of every selected source plus all executed helpers.
 stage=WORK/'relocated';root=stage/'source';support=stage/'own';checks=stage/'checks';prod=stage/'producer'
 for d in (root,support,checks,prod):d.mkdir(parents=True)
 fullpins=json.loads((OWN/'source-pins.json').read_text())
 stagepins=dict(fullpins);stagepins.update(json.loads((ASSETS/'source-pins.json').read_text()))
 for rel in stagepins:
  p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(SOURCE/rel,p)
 absolute=support/'INDEPENDENT_ABSOLUTE.json';shutil.copyfile(paths[3],absolute)
 ownreport=json.loads(absolute.read_text())
 for rel in [*ownreport['own_sha256'],'FREEZE.json','source-pins.json']:
  p=support/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(OWN/rel,p)
 relocated_checker=checks/checker.name;shutil.copyfile(checker,relocated_checker)
 pins=checks/'source-pins.json';shutil.copyfile(ASSETS/'source-pins.json',pins)
 response=prod/'RESPONSE.json';shutil.copyfile(paths[5],response)
 response_source=prod/'response_adapter.py';shutil.copyfile(paths[6],response_source)
 shutil.copyfile(RESPONSE/'checks/FREEZE.json',prod/'FREEZE.json')
 relocated_paths=(relocated_checker,root,pins,absolute,support,response,response_source)
 for opt in (False,True):
  out=WORK/('PORTABLE_OPTIMIZED.json' if opt else 'PORTABLE_NORMAL.json')
  got=invoke(*relocated_paths,out,opt)
  need(got.returncode==0,got.stderr)
  need(out.read_bytes()==outputs[0].read_bytes(),'relocated canonical equality')
 records=[]
 def rejected(name,expected,paths=relocated_paths,**over):
  for opt in (False,True):
   out=WORK/(name+('-O' if opt else '')+'.json')
   got=invoke(*paths,out,opt,**over)
   need(got.returncode!=0 and not out.exists() and expected in got.stderr,'effective rejection:'+name+':'+got.stderr)
   records.append(dict(case=name,optimized=opt,nonzero_exit=True,no_success_artifact=True,expected_rejection=expected))
 original_response=response.read_bytes()
 for name,change,marker in [
  ('response_curvature_leaf',lambda d:d['ledger']['second_coordinate_error_H_caps'].__setitem__('rr',str(__import__('fractions').Fraction(d['ledger']['second_coordinate_error_H_caps']['rr'])+1)),'every response ledger node independently equal'),
  ('response_frozen_error',lambda d:d['frozen_initial_error_H_caps'].__setitem__('v0_r','0'),'producer retains exact actual-error caps'),
  ('response_original_Fx',lambda d:d['constants'].__setitem__('Fx','0'),'every producer operator source contraction commutator constant'),
  ('response_absolute_jet',lambda d:d['candidate_initial_absolute_H_caps'].__setitem__('u_rr',str(__import__('fractions').Fraction(d['candidate_initial_absolute_H_caps']['u_rr'])+1)),'every response ledger node independently equal')]:
  d=json.loads(original_response);change(d);write(response,d)
  rejected(name,marker,respsha=sha(response));response.write_bytes(original_response)
 original_absolute=absolute.read_bytes();d=json.loads(original_absolute);d['source_postauth_equal']=False;write(absolute,d)
 rejected('absolute_authentication_flag','independent absolute report role',abssha=sha(absolute));absolute.write_bytes(original_absolute)
 own_helper=support/'scripts/pinned_independent_exact.py';raw=own_helper.read_bytes();own_helper.write_bytes(raw+b'\n# deliberate disk corruption\n')
 rejected('absolute_arithmetic_source_disk','selected member hash: scripts/pinned_independent_exact.py');own_helper.write_bytes(raw)
 raw=response_source.read_bytes();response_source.write_bytes(raw+b'\n# deliberate disk corruption\n')
 rejected('producer_response_source_disk','producer response source hash');response_source.write_bytes(raw)
 fx=root/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json';raw=fx.read_bytes();fx.write_bytes(raw+b'\n')
 rejected('original_Fx_disk','selected member hash: research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json');fx.write_bytes(raw)
 raw=pins.read_bytes();pins.write_bytes(raw+b'\n');rejected('selected_pin_map_disk','selected pin map hash before any scientific read');pins.write_bytes(raw)
 alias=support/'alias-absolute.json';alias.symlink_to(absolute.name)
 alt=list(relocated_paths);alt[3]=alias
 rejected('absolute_report_symlink','regular file required: alias-absolute.json',tuple(alt))
 scripts=support/'scripts';real_scripts=support/'scripts-real';scripts.rename(real_scripts);scripts.symlink_to(real_scripts.name,target_is_directory=True)
 rejected('absolute_source_ancestor_symlink','selected member symlink: scripts/pinned_independent_exact.py')
 scripts.unlink();real_scripts.rename(scripts)
 need(all(sha(p)==h for p,h in before.items()),'executed independent checker/source map unchanged')
 receipt=dict(status='PASS_REAL_SUBPROCESS_INDEPENDENT_RESPONSE_CONTROLS',
  checker_sha256=before[checker],canonical_report_sha256=sha(outputs[0]),
  normal_optimized_byte_identical=True,relocated_normal_optimized_byte_identical=True,
  relocated_selected_scientific_source_count=len(stagepins),
  report_corruptions_rehashed_before_invocation=True,controls=records,
  all_rejections_have_nonzero_exit_and_no_success_artifact=True,
  scope='Separate exact rational ledger controls; independent interval expression and coverage controls remain separate.')
 write(BASE/'CONTROLS.json',receipt)
 print(json.dumps({'status':receipt['status'],'corruption_subprocess_count':len(records),'canonical_report_sha256':receipt['canonical_report_sha256']}))

if __name__=='__main__':main()
