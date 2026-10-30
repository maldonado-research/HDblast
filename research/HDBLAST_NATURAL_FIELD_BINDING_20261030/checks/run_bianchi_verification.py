#!/usr/bin/env python3
"""Authenticated isolated normal/optimized replays and corrupted subprocesses."""
import argparse, ast, concurrent.futures, hashlib, json, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
PINS={'check_natural_bulk_bianchi.py': '8f8e7ec7ee3955d0f115e245d40a59b37e3173e5793dccaea002a7b2ce0a8175', 'codazzi_laurent.py': '178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2', 'angular_generic.py': 'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9', 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json': '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'}
CONTROLS=('reduced_instead_of_natural','forcing_instead_of_natural','half_second_variation','omit_deltaB_E1','omit_S1_deta','positive_time_contraction','missing_trace_half','missing_radial_connection','missing_FRW_connection','wrong_potential_chain','frozen_rho_radial','wrong_scalar_source_sign','premature_harmonic_scalar')
def need(ok,label):
 if not ok:raise ValueError(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def auth(root):
 out={}
 for rel,pin in PINS.items():
  p=root/rel if rel.startswith('research/') else HERE/rel
  need(p.is_file() and not p.is_symlink() and sha(p)==pin,'SOURCE_AUTHENTICATION:'+rel)
  out[rel]={'sha256':pin,'bytes':p.stat().st_size}
  if p.suffix=='.py':
   tree=ast.parse(p.read_text());need(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'ASSERT_FORBIDDEN:'+rel)
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
 root=a.repo_root.resolve();out=a.output_dir.resolve();need(not out.exists(),'REFUSE_OVERWRITE');out.mkdir(parents=True)
 before=auth(root);driver_pin=sha(Path(__file__))
 def run(mode,mutation):
  auth(root);tag=mode+('_'+mutation if mutation else '')
  report=out/(tag+'.json');cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='OPTIMIZED' else [])+[str(HERE/'check_natural_bulk_bianchi.py'),'--repo-root',str(root),'--output',str(report)]+(['--mutation',mutation] if mutation else [])
  started=time.monotonic();result=subprocess.run(cmd,text=True,capture_output=True,timeout=180);elapsed=time.monotonic()-started
  evidence={'command':result.args,'elapsed_seconds':round(elapsed,6),'returncode':result.returncode}
  auth(root)
  if mutation:
   marker='ValueError: NONZERO:';need(result.returncode!=0 and marker in result.stderr and not report.exists(),'CORRUPTION_NOT_REJECTED:'+tag+':'+result.stderr[-2000:])
   witness=result.stderr.split(marker,1)[1].splitlines()[0]
   return {**evidence,'mode':mode,'mutation':mutation,'exact_nonzero_witness':witness,'successful_report_absent':True}
  need(result.returncode==0 and report.is_file(),'BASELINE_FAILURE:'+tag+':'+result.stderr)
  value=json.loads(report.read_text());need(value['status']=='PASS_EXACT_RESIDUAL_CORRECTED_NATURAL_BULK_BIANCHI' and len(value['rows'])==8,'BASELINE_STATUS:'+tag)
  return {**evidence,'mode':mode,'sha256':sha(report),'rows':len(value['rows']),'ordinary_checks':len(value['checks']),'angular_checks':len(value['angular_coefficient_checks'])}
 baselines=[run(mode,None) for mode in ('NORMAL','OPTIMIZED')]
 need(baselines[0]['sha256']==baselines[1]['sha256'],'NORMAL_OPTIMIZED_CANONICAL_MISMATCH')
 jobs=[(mode,m) for mode in ('NORMAL','OPTIMIZED') for m in CONTROLS]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:controls=list(pool.map(lambda job:run(*job),jobs))
 for m in CONTROLS:
  rows=[x for x in controls if x['mutation']==m];need(rows[0]['exact_nonzero_witness']==rows[1]['exact_nonzero_witness'],'NORMAL_OPTIMIZED_WITNESS_MISMATCH:'+m)
 need(auth(root)==before,'SOURCE_CLOSURE_CHANGED');need(sha(Path(__file__))==driver_pin,'DRIVER_CHANGED_DURING_EXECUTION')
 report={'status':'PASS_PRIVATE_NATURAL_BULK_BIANCHI_VERIFICATION','baselines':baselines,'normal_optimized_reports_byte_identical':True,'fresh_corrupted_subprocesses':controls,'closure_pre_post_equal':True,'source_authentication':before,'driver_sha256':driver_pin,'scope':'Internal review. No external specialist review, publishing, source edits, actual first-equation membership or unconditional homogeneous identity claimed.'}
 (out/'VERIFICATION.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':report['status'],'baselines':len(baselines),'rejected_corruptions':len(controls)}))
if __name__=='__main__':main()
