#!/usr/bin/env python3
"""Run real checker-source and closed-model corruptions in isolated copies."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);args=ap.parse_args()
    output_dir=args.output_dir;need(not output_dir.exists() and not output_dir.is_symlink(),'FRESH_OUTPUT_DIR');output_dir.mkdir(parents=True)
    scratch=output_dir/'controls';scratch.mkdir()
    base=(HERE/'independent_numeric.py').read_text();rows=[]
    changes={
        'wrong_boundary_lift_sign':('LIFT_SIGN = -1','LIFT_SIGN = 1','LIFT_SIGN_WITNESS:got=7/3;expected=-7/3'),
        'wrong_reflected_endpoint_sign':('REFLECTED_LOWER_SIGN = 1','REFLECTED_LOWER_SIGN = -1','REFLECTED_KERNEL_WITNESS:got=0;expected=5/4'),
        'missing_fixed_gauge_quadratic':('SCREEN_GAUGE_FACTOR = 1','SCREEN_GAUGE_FACTOR = 0','GAUGE_QUADRATIC_WITNESS:got=0;expected=-1/2'),
        'missing_scalar_feedback_tail':('FEEDBACK_KEEP_BD = 1','FEEDBACK_KEEP_BD = 0','SCALAR_FEEDBACK_WITNESS:got=1/7;expected=2969/5544')}
    for name,(old,new,witness) in changes.items():
        d=scratch/name;d.mkdir();checker=d/'independent_numeric.py'
        need(base.count(old)==1,'EXACT_ONE_MUTATION');checker.write_text(base.replace(old,new))
        for optimized in (False,True):
            output=d/('OPTIMIZED.json' if optimized else 'NORMAL.json')
            command=[sys.executable]+(['-O'] if optimized else [])+[str(checker),'--source-root',str(args.source_root),'--output',str(output),'--guards-only']
            run=subprocess.run(command,capture_output=True,text=True)
            # Store the exact reported witness, not the whole traceback.
            lines=[s for s in run.stderr.splitlines() if s.startswith('ValueError: ')]
            actual=lines[-1].split('ValueError: ',1)[1] if lines else run.stderr.strip()
            need(run.returncode!=0 and not output.exists(),'NEGATIVE_CONTROL_REJECT:'+name)
            need(actual==witness,'EXACT_FAILED_WITNESS:'+name)
            rows.append({'control':name,'optimized':optimized,'exit_code':run.returncode,'no_success_report':not output.exists(),'checker_sha256':sha(checker),'failed_witness':actual})
    # Authenticate before imports: copy only the first capsule needed by the
    # ordered pin map, then alter a listed file or replace it with a symlink.
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text());first=next(k for k in pins if not k.startswith('research/'))
    for name,kind in (('closed_source_bytes','bytes'),('closed_source_symlink','symlink')):
        d=scratch/name;d.mkdir();checker=d/'independent_numeric.py';checker.write_text(base);shutil.copyfile(HERE/'SOURCE_PINS.json',d/'SOURCE_PINS.json')
        repo=d/'source';target=repo/'research'/first;target.parent.mkdir(parents=True)
        shutil.copytree(args.source_root/'research'/first,target)
        manifest=json.loads((target/'MANIFEST.json').read_text())
        rel=manifest['files'][0]['path'];member=target/rel
        if kind=='bytes':member.write_bytes(member.read_bytes()+b'\nCORRUPTION\n')
        else:
            member.unlink();member.symlink_to(args.source_root/'research'/first/rel)
        for optimized in (False,True):
            output=d/('OPTIMIZED.json' if optimized else 'NORMAL.json')
            command=[sys.executable]+(['-O'] if optimized else [])+[str(checker),'--source-root',str(repo),'--output',str(output)]
            run=subprocess.run(command,capture_output=True,text=True)
            lines=[s for s in run.stderr.splitlines() if s.startswith('ValueError: ')]
            actual=lines[-1].split('ValueError: ',1)[1] if lines else run.stderr.strip()
            need(run.returncode!=0 and not output.exists() and actual=='CLOSED_MEMBER:'+first+'/'+rel,'CLOSED_SOURCE_REJECT:'+name)
            rows.append({'control':name,'optimized':optimized,'exit_code':run.returncode,'no_success_report':not output.exists(),'failed_witness':actual,'mutated_member':first+'/'+rel})
    report={'status':'PASS_TWELVE_REAL_NEGATIVE_CONTROL_SUBPROCESSES','controls':rows,'checker_source_sha256':sha(HERE/'independent_numeric.py'),'scope':'Four independent computational source changes and two closed-source corruptions; each rejected before a success artifact in normal and optimized modes.'}
    (output_dir/'CONTROLS.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(report['status'])

if __name__=='__main__':main()
