#!/usr/bin/env python3
"""Exact failure witnesses and absent artifacts for response mutations."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent;RESPONSE=HERE.parent;CAP=RESPONSE.parent
class Failure(ValueError):pass
def need(ok,w):
 if not ok:raise Failure(w)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return (json.dumps(x,indent=2,sort_keys=True)+'\n').encode()
CASES={'drop_conditioned_scalar':('CONDITIONED_ADAPTER_INCLUSION',0),
       'reverse_conditioned_scalar':('CONDITIONED_ADAPTER_INCLUSION',0),
       'reverse_volume_source':('REGISTERED_SIGN_GATES',1),
       'omit_proper_mean_correction':('FROZEN_SCIENTIFIC_REFERENCE',1),
       'drop_actual_source_history':('FROZEN_SCIENTIFIC_REFERENCE',1),
       'old_total_cap_substitution':('FROZEN_SCIENTIFIC_REFERENCE',1)}
def invoke(producer,repo,out,mode,case):
 command=[sys.executable,'-I','-S','-B']+(['-O']if mode=='optimized'else [])+[str(producer),'--repo-root',str(repo),'--output-dir',str(out)]
 if case:command+=['--inject',case]
 proc=subprocess.run(command,capture_output=True,text=True,timeout=180)
 lines=proc.stdout.splitlines();need(len(lines)==1 and not proc.stderr,'EXACT_SINGLE_FINAL_JSON')
 final=json.loads(lines[0]);positive=case is None
 expected={'status':'PASS','witness':'NEW_TOTAL_MEAN_ENCLOSURES','actual_first_oracle_calls':1}if positive else {'status':'FAIL','witness':CASES[case][0],'actual_first_oracle_calls':CASES[case][1]}
 need(final==expected,'EXACT_FINAL_WITNESS');need(proc.returncode==(0 if positive else 2),'EXACT_EXIT_CODE')
 need(out.exists()==positive and (out/'REPORT.json').exists()==positive,'EXACT_ARTIFACT_CONTRACT')
 if positive:need((out/'REPORT.json').read_bytes()==(RESPONSE/'references/FINAL_NORMAL.json').read_bytes(),'POSITIVE_FROZEN_REPORT_BYTES')
 return {'mode':mode,'subprocess_exit':proc.returncode,'final_witness':final,'output_directory_absent':not out.exists(),'report_absent':not(out/'REPORT.json').exists(),'stderr_empty':True}
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 repo=a.repo_root.resolve(strict=True);output=a.output.absolute();need(not output.exists()and all(not x.is_symlink()for x in [output,*output.parents]),'FRESH_OUTPUT')
 need(all(output!=x and x not in output.parents for x in [repo,CAP]),'OUTPUT_OUTSIDE_SOURCES')
 own=sha(Path(__file__));rows=[];positives=[]
 with tempfile.TemporaryDirectory(prefix='prepared-response-controls-')as temporary:
  scratch=Path(temporary);relocated=scratch/'relocated';shutil.copytree(CAP,relocated)
  for location,producer in [('original',RESPONSE/'producer.py'),('relocated',relocated/'response/producer.py')]:
   for mode in ('normal','optimized'):
    positives.append({'location':location,**invoke(producer,repo,scratch/(location+'-'+mode),mode,None)})
  for case in CASES:
   for mode in ('normal','optimized'):
    rows.append({'case':case,'expected_final_witness':{'status':'FAIL','witness':CASES[case][0],'actual_first_oracle_calls':CASES[case][1]},**invoke(RESPONSE/'producer.py',repo,scratch/(case+'-'+mode),mode,case)})
 need(sha(Path(__file__))==own,'CONTROL_SOURCE_UNCHANGED')
 oracle=sum(r['final_witness']['actual_first_oracle_calls']for r in rows+positives)
 result={'status':'PASS_STRICT_RESPONSE_CONTROLS','source_sha256':own,'distinct_negative_cases':len(CASES),
  'actual_rejected_subprocesses':len(rows),'positive_subprocesses':len(positives),'new_actual_oracle_calls':oracle,
  'exact_final_json_witnesses_required':True,'positive_relocation_separate_from_negatives':True,'rejections':rows,'positive_checks':positives}
 output.write_bytes(canonical(result));print(json.dumps({'status':'PASS','witness':'STRICT_RESPONSE_CONTROLS','actual_rejected_subprocesses':len(rows),'positive_subprocesses':len(positives),'new_actual_oracle_calls':oracle},sort_keys=True))
if __name__=='__main__':
 try:main()
 except (Failure,ValueError,OSError,KeyError)as exc:
  print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
