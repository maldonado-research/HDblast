#!/usr/bin/env python3
"""Verify hashes and replay bounded checks in a temporary copy; leave package intact."""
from pathlib import Path
import argparse, hashlib, json, os, shutil, subprocess, sys, tempfile, time

root=Path(__file__).resolve().parent
p=argparse.ArgumentParser()
p.add_argument('--symbolic-python',default=sys.executable)
p.add_argument('--hashes-only',action='store_true')
p.add_argument('--receipt',type=Path,help='Optional receipt path outside the frozen package')
a=p.parse_args()
manifest=root/'MANIFEST.sha256.json'
if not manifest.exists():
    raise SystemExit('MANIFEST.sha256.json is required; use the final package.')
expected=json.loads(manifest.read_text())['files']
failures=[]
for rel,digest in expected.items():
    f=root/rel
    if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=digest:
        failures.append(rel)
if failures:raise SystemExit('Hash failures: '+repr(failures))
out=dict(status='PASS_HASHES',payload_files_checked=len(expected),stages=[],
         scope='Hash verification and bounded check replay. No spectral solve, PDE evolution or initial-data shooting is repeated.')
if not a.hashes_only:
    jobs=[
      ('saved_run_analysis','analyze_controls.py',[],False),
      ('proper_clock','literature/verify_proper_clock.py',[],False),
      ('archived_controls','numerics_review/check_archived_numerics.py',['inputs/CHAT13_REFERENCE.zip'],False),
      ('replacement_solver_controls','numerics_review/review_new_solver.py',[],False),
      ('constraint_transport','numerics_review/verify_constraint_transport.py',[],True),
      ('initial_data_algebra','initial_data/verify_constraint_seed_algebra.py',[],True),
      ('independent_initial_data_review','numerics_review/verify_seed_compatibility_independent.py',[],True),
      ('sampled_initial_constraints','initial_data/check_sampled_constraints.py',[],False),
    ]
    # This also works after extraction elsewhere. Historic local paths in input
    # provenance are metadata only and are never opened by these replay stages.
    with tempfile.TemporaryDirectory(prefix='hdblast-sep22-replay-') as tmp:
        copied=Path(tmp)/'checkpoint'
        shutil.copytree(root,copied,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for name,script,args,symbolic in jobs:
            executable=a.symbolic_python if symbolic else sys.executable
            env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['OPENBLAS_NUM_THREADS']='1'
            if symbolic and executable!=sys.executable:env.pop('PYTHONPATH',None)
            # Turn relative PYTHONPATH entries into absolute paths before changing cwd.
            elif env.get('PYTHONPATH'):
                env['PYTHONPATH']=os.pathsep.join(str(Path(x).resolve()) for x in env['PYTHONPATH'].split(os.pathsep) if x)
            start=time.monotonic()
            result=subprocess.run([executable,'-B',script,*args],cwd=copied,env=env,text=True,
                                  stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=120)
            row=dict(name=name,returncode=result.returncode,seconds=time.monotonic()-start,
                     stdout=result.stdout,stderr=result.stderr)
            out['stages'].append(row)
            if result.returncode:
                out['status']='FAIL_CHECK_REPLAY'
                print(json.dumps(row,indent=2))
                if a.receipt:a.receipt.write_text(json.dumps(out,indent=2)+'\n')
                raise SystemExit(1)
            print(name+': PASS',flush=True)
        out['status']='PASS_HASHES_AND_BOUNDED_CHECK_REPLAY'
if a.receipt:
    if root in a.receipt.resolve().parents:
        raise SystemExit('Receipt must be outside the frozen package.')
    a.receipt.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='stages'}))
