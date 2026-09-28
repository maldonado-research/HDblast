#!/usr/bin/env python3
"""Hash verification and bounded independent checks in a disposable copy."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile,time

ROOT=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--hashes-only',action='store_true');p.add_argument('--receipt',type=Path)
a=p.parse_args()
if a.receipt and (a.receipt.resolve()==ROOT or ROOT in a.receipt.resolve().parents):
    raise SystemExit('Write replay receipt outside the frozen package.')
manifest=json.loads((ROOT/'MANIFEST.sha256.json').read_text())
for rel,digest in manifest['files'].items():
    f=ROOT/rel
    if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=digest:raise SystemExit('Hash failure: '+rel)
out=dict(status='PASS_HASHES',payload_files_checked=len(manifest['files']),stages=[],
         scope='Hash checks, saved-data analyses, exact identities, independent static IVPs, initial-data reconstruction and a ten-step matched-time replay. No complete PDE trajectory is rerun.')
jobs=[('external_source_audit','source_audit/recompute_chat14.py'),
      ('static_branch_independent_IVPs','static_branch/independent_branch_checks.py'),
      ('matter_exact_identities','matter/verify_matter_extension.py'),
      ('crossing_screen','matter/analyze_crossing_eligibility.py'),
      ('sampling_difference_identities','evolution_review/verify_sampling_and_difference.py'),
      ('initial_data_reconstruction','evolution_review/check_wrapper_initialization.py'),
      ('inverse_sampling_tolerance','evolution_review/check_inverse_tolerance.py'),
      ('saved_full_constraint_analysis','evolution_review/analyze_evolution_residuals.py'),
      ('matched_time_refinement','evolution_review/matched_time_refinement.py'),
      ('result_accounting','build_results_summary.py')]
if not a.hashes_only:
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['OPENBLAS_NUM_THREADS']='1'
    if env.get('PYTHONPATH'):
        env['PYTHONPATH']=os.pathsep.join(str(Path(x).resolve()) for x in env['PYTHONPATH'].split(os.pathsep) if x)
    with tempfile.TemporaryDirectory(prefix='hdblast-branch-matter-replay-') as temp:
        copied=Path(temp)/'checkpoint';shutil.copytree(ROOT,copied,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name,script in jobs:
            if not (copied/script).is_file():raise SystemExit('Missing replay script: '+script)
            start=time.monotonic()
            run=subprocess.run([sys.executable,'-B',script],cwd=copied,env=env,text=True,
                               stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=180)
            row=dict(name=name,script=script,returncode=run.returncode,seconds=time.monotonic()-start,
                     stdout=run.stdout,stderr=run.stderr);out['stages'].append(row)
            if run.returncode:
                out['status']='FAIL_BOUNDED_REPLAY'
                if a.receipt:a.receipt.write_text(json.dumps(out,indent=2)+'\n')
                print(json.dumps(row,indent=2));raise SystemExit(1)
            print(name+': PASS',flush=True)
    out['status']='PASS_HASHES_AND_BOUNDED_REPLAY'
if a.receipt:a.receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='stages'}))
