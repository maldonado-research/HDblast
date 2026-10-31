#!/usr/bin/env python3
"""Authenticated isolated normal/optimized replays and corrupted subprocesses."""
import argparse, ast, concurrent.futures, hashlib, json, subprocess, sys, time
from pathlib import Path
HERE=Path(__file__).resolve().parent
PINS={'covector_wave_rows.py':'cb025648deda6eec4cfd77b3bb96c6cf6e0713f2b55987b3fdfff22dbb066bbf','check_reduced_source_binding.py': '8863101fa6733836f50ea06d9067dc64f4ee639f5629c0928a4456df32fec017', 'codazzi_laurent.py': '178b0a4e58db63eae6c24112ed9d91bf430f70dca9f6f8126c5e352e429f8fc2', 'angular_generic.py': 'fd14d8cc3ac6f58c206b784c5c57d15430f63f38ca091ac03c7769f9bfc89df9', 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json': '34b10b62d0b36a4506a35407a591a9fdff319971994f25e0b5f031f48eb4f615'}
CONTROLS=('premature_QC_background','q_missing_two','q_positive_time','PC_missing_two','PC_positive_time','drop_PC1_role','wrong_gauge_gradient_sign','unnormalized_symgrad','omit_TT_PC1','drop_BPC1','QE_instead_general','QEred_instead_general','forcing_instead_general','premature_Qphi','omit_deltaB_E1','omit_S1_deta','wrong_wave_sign','missing_half_wave','frozen_rho_radial','wrong_potential_chain')
def need(ok,label):
 if not ok:raise ValueError(label)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def auth(root):
 out={}
 cap=root/'research/HDBLAST_NATURAL_FIELD_BINDING_20261030';manifest=cap/'MANIFEST.json';manifest_pin='5121e4304fad7ab5856977b38c6ebd9fdc205900b543bc1c3d97384d246f378f'
 need(manifest.is_file() and not manifest.is_symlink() and sha(manifest)==manifest_pin,'OLD_MANIFEST_PIN')
 data=json.loads(manifest.read_text());need(len(data['files'])==data['file_count_excluding_this_manifest'],'OLD_MANIFEST_COUNT')
 for row in data['files']:
  rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts,'OLD_MANIFEST_PATH');f=cap/rel
  need(f.is_file() and not f.is_symlink() and f.stat().st_size==row['bytes'] and sha(f)==row['sha256'],'OLD_CAPSULE_SOURCE:'+str(rel))
 out['old_natural_manifest']={'sha256':manifest_pin,'file_count':len(data['files'])}

 for rel,pin in PINS.items():
  p=root/rel if rel.startswith('research/') else HERE/rel
  need(p.is_file() and not p.is_symlink() and sha(p)==pin,'SOURCE_AUTHENTICATION:'+rel)
  out[rel]={'sha256':pin,'bytes':p.stat().st_size}
  if p.suffix=='.py':
   tree=ast.parse(p.read_text());need(not any(isinstance(n,ast.Assert) for n in ast.walk(tree)),'ASSERT_FORBIDDEN:'+rel)
   allowed={'argparse','ast','hashlib','importlib','json','sys','fractions','itertools','math','pathlib','functools'}
   for node in ast.walk(tree):
    if isinstance(node,ast.Import):need(all(x.name.split('.')[0] in allowed for x in node.names),'IMPORT_CLOSURE:'+rel)
    if isinstance(node,ast.ImportFrom):need(node.module and node.module.split('.')[0] in allowed,'IMPORT_CLOSURE:'+rel)
    if isinstance(node,ast.Call) and isinstance(node.func,ast.Name):need(node.func.id not in {'eval','exec','compile','__import__'},'DYNAMIC_CODE_FORBIDDEN:'+rel)
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);ap.add_argument('--workers',type=int,default=2);a=ap.parse_args()
 root=a.repo_root.resolve();out=a.output_dir.resolve();need(not out.exists(),'REFUSE_OVERWRITE');out.mkdir(parents=True)
 before=auth(root);driver_pin=sha(Path(__file__))
 def run(mode,mutation):
  need(auth(root)==before,'CLOSURE_CHANGED_BEFORE');tag=mode+('_'+mutation if mutation else '')
  report=out/(tag+'.json');need(not report.exists(),'OUTPUT_PREEXISTS:'+tag);cmd=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='OPTIMIZED' else [])+[str(HERE/'check_reduced_source_binding.py'),'--repo-root',str(root),'--output',str(report)]+(['--mutation',mutation] if mutation else [])
  started=time.monotonic();result=subprocess.run(cmd,text=True,capture_output=True,timeout=600);elapsed=time.monotonic()-started
  evidence={'command':result.args,'elapsed_seconds':round(elapsed,6),'returncode':result.returncode}
  need(auth(root)==before,'CLOSURE_CHANGED_AFTER:'+tag)
  if mutation:
   marker='ValueError: NONZERO:';need(result.returncode!=0 and marker in result.stderr and not report.exists(),'CORRUPTION_NOT_REJECTED:'+tag+':'+result.stderr[-2000:])
   witness=result.stderr.split(marker,1)[1].splitlines()[0]
   return {**evidence,'mode':mode,'mutation':mutation,'exact_nonzero_witness':witness,'successful_report_absent':True}
  need(result.returncode==0 and report.is_file(),'BASELINE_FAILURE:'+tag+':'+result.stderr)
  value=json.loads(report.read_text());need(value['status']=='PASS_EXACT_REDUCED_SOURCE_BINDING' and len(value['rows'])==8,'BASELINE_STATUS:'+tag)
  return {**evidence,'mode':mode,'sha256':sha(report),'rows':len(value['rows']),'ordinary_checks':len(value['checks']),'angular_checks':len(value['angular_coefficient_checks'])}
 baselines=[run(mode,None) for mode in ('NORMAL','OPTIMIZED')]
 need(baselines[0]['sha256']==baselines[1]['sha256'],'NORMAL_OPTIMIZED_CANONICAL_MISMATCH')
 jobs=[(mode,m) for mode in ('NORMAL','OPTIMIZED') for m in CONTROLS]
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:controls=list(pool.map(lambda job:run(*job),jobs))
 for m in CONTROLS:
  rows=[x for x in controls if x['mutation']==m];need(rows[0]['exact_nonzero_witness']==rows[1]['exact_nonzero_witness'],'NORMAL_OPTIMIZED_WITNESS_MISMATCH:'+m)
 need(auth(root)==before,'SOURCE_CLOSURE_CHANGED');need(sha(Path(__file__))==driver_pin,'DRIVER_CHANGED_DURING_EXECUTION')
 report={'status':'PASS_PRIVATE_REDUCED_SOURCE_VERIFICATION','baselines':baselines,'normal_optimized_reports_byte_identical':True,'fresh_corrupted_subprocesses':controls,'closure_pre_post_equal':True,'source_authentication':before,'driver_sha256':driver_pin,'scope':'Internal exact source/operator review; first-natural and harmonic specialization remains conditional. No actual numerical third-jet adapter or publication performed.'}
 (out/'VERIFICATION.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':report['status'],'baselines':len(baselines),'rejected_corruptions':len(controls)}))
if __name__=='__main__':main()
