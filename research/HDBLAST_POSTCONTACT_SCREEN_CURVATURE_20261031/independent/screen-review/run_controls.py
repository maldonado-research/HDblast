#!/usr/bin/env python3
"""Real source and independent-equation corruptions in isolated scratch copies."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--source-root',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    here=Path(__file__).resolve().parent
    need(not args.output.exists(),'FRESH_CONTROLS_REPORT')
    scratch=here/'controls-scratch'
    need(not scratch.exists(),'FRESH_CONTROLS_SCRATCH')
    scratch.mkdir()
    source=here/'check_screen.py'
    original=source.read_text()
    cases=[('omit_epsilon', 'EPSILON_TERM_ENABLED = True', 'EPSILON_TERM_ENABLED = False',
            'SYMBOLIC_FOUR_CONTRIBUTION_IDENTITY'),
           ('wrong_B_sign', 'COMPACT_B_SIGN = 1', 'COMPACT_B_SIGN = -1',
            'SYMBOLIC_FOUR_CONTRIBUTION_IDENTITY'),
           ('extra_half_parameter', 'LINEAR_PARAMETER_FACTOR = Q(1)',
            'LINEAR_PARAMETER_FACTOR = Q(1,2)', 'SECOND_PARAMETER_AND_TIME_UNITS'),
           ('wrong_y_connection', '+a*jr/(ell*rho)', '-a*jr/(ell*rho)',
            'PHYSICAL_YY_CONNECTION')]
    results=[]
    for name,old,new,diagnostic in cases:
        need(original.count(old)==1,'UNIQUE_MUTATION:'+name)
        mutant=scratch/(name+'.py')
        mutant.write_text(original.replace(old,new))
        for optimized in (False,True):
            output=scratch/(name+('_optimized' if optimized else '_normal')+'.json')
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(mutant),
                 '--source-root',str(args.source_root),'--output',str(output)]
            run=subprocess.run(cmd,capture_output=True,text=True)
            need(run.returncode!=0 and not output.exists() and diagnostic in run.stderr,
                 'MUTANT_NOT_REJECTED:'+name)
            results.append({'control':name,'optimized':optimized,'exit_code':run.returncode,
                            'no_success_report':True,'diagnostic':diagnostic})
    # Copy only the four independently pinned files. Mutation is caught before
    # the old screen module can execute, in both normal and optimized Python.
    pins=json.loads((here/'SCREEN_NORMAL.json').read_text())['source_pins']
    clone=scratch/'source-auth-root'
    for rel in pins:
        dest=clone/rel;dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(args.source_root/rel,dest)
    target=clone/'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py'
    target.write_text(target.read_text()+'\n# actual corrupted source bytes\n')
    for optimized in (False,True):
        output=scratch/('source_auth'+('_optimized' if optimized else '_normal')+'.json')
        run=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(source),
                           '--source-root',str(clone),'--output',str(output)],capture_output=True,text=True)
        need(run.returncode!=0 and not output.exists() and 'SOURCE_HASH:' in run.stderr,
             'SOURCE_CORRUPTION_NOT_REJECTED')
        results.append({'control':'source_byte_corruption','optimized':optimized,
                        'exit_code':run.returncode,'no_success_report':True,
                        'diagnostic':'SOURCE_HASH before project import'})
    target.unlink()
    target.symlink_to(args.source_root/'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py')
    for optimized in (False,True):
        output=scratch/('source_symlink'+('_optimized' if optimized else '_normal')+'.json')
        run=subprocess.run([sys.executable]+(['-O'] if optimized else [])+[str(source),
                           '--source-root',str(clone),'--output',str(output)],capture_output=True,text=True)
        need(run.returncode!=0 and not output.exists() and 'SOURCE_TYPE:' in run.stderr,
             'SOURCE_SYMLINK_NOT_REJECTED')
        results.append({'control':'source_symlink','optimized':optimized,
                        'exit_code':run.returncode,'no_success_report':True,
                        'diagnostic':'SOURCE_TYPE before project import'})
    target.unlink()
    normal=(here/'SCREEN_NORMAL.json').read_bytes()
    optimized=(here/'SCREEN_OPTIMIZED.json').read_bytes()
    need(normal==optimized,'NORMAL_OPTIMIZED_REPORT_IDENTITY')
    report={'status':'PASS_PRIVATE_SCREEN_REAL_NEGATIVE_CONTROLS',
            'normal_optimized_byte_identical':True,
            'report_sha256':hashlib.sha256(normal).hexdigest(),
            'rejected_subprocesses':results,
            'scope':'Independent formula/unit corruption and selected-source authentication; no PDE campaign.'}
    with args.output.open('x') as f:
        f.write(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(report['status'])


if __name__=='__main__':
    main()
