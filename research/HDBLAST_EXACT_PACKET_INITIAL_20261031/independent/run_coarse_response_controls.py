#!/usr/bin/env python3
"""Execute meaningful on-disk corruptions of the independent response ledger."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--checker',type=Path,required=True)
    ap.add_argument('--source-root',type=Path,required=True)
    ap.add_argument('--initial-report',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output_dir.exists(),'OUTPUT_DIRECTORY_MUST_BE_ABSENT')
    checker_sha=sha(args.checker);own_sha=sha(Path(__file__))
    initial_sha=sha(args.initial_report)
    original=args.checker.read_text()
    mutations=[
        ('omit_intermediate_initial_acceleration',
         'Dx=P*(f00+T*f0x+act[\'M\']*Bs+Gx)',
         'Dx=P*(T*f0x+act[\'M\']*Bs+Gx)',
         'REPRODUCE_RETAINED_COMPLETE_DIFFERENTIATED_LEDGER'),
        ('omit_correction_initial_acceleration',
         'Cx=P*(r00+T*r0x)',
         'Cx=P*(T*r0x)',
         'REPRODUCE_RETAINED_COMPLETE_DIFFERENTIATED_LEDGER'),
        ('omit_original_Fx',
         "+(diff['D']+comm)*Bs+Fx",
         "+(diff['D']+comm)*Bs",
         'REPRODUCE_RETAINED_COMPLETE_DIFFERENTIATED_LEDGER'),
        ('omit_differentiated_boundary_source',
         "Dx=P*(f00+T*f0x+act['M']*Bs+Gx)",
         "Dx=P*(f00+T*f0x)",
         'REPRODUCE_RETAINED_COMPLETE_DIFFERENTIATED_LEDGER'),
        ('interchange_absolute_and_error_curvature',
         "errors.update(v0=olderr['v0'],v0_r=olderr['v0_r'])",
         "errors.update(v0=olderr['v0'],v0_r=olderr['v0_r'])\n    absolute['u_rr']=errors['e0_rr']",
         'EXPLICIT_NONZERO_OMISSION_WITNESSES'),
        ('omit_unchanged_velocity_gradient_error',
         "errors.update(v0=olderr['v0'],v0_r=olderr['v0_r'])",
         "errors.update(v0=olderr['v0'],v0_r=Q(0))",
         'EXPLICIT_NONZERO_OMISSION_WITNESSES'),
    ]
    args.output_dir.mkdir(parents=True)
    results=[]
    for name,before,after,expected in mutations:
        need(original.count(before)==1,'MUTATION_TARGET_MUST_BE_UNIQUE:'+name)
        corrupt=args.output_dir/(name+'.py')
        corrupt.write_text(original.replace(before,after))
        artifact=args.output_dir/(name+'.json')
        command=[sys.executable,'-I','-S','-B']
        if sys.flags.optimize:
            command.append('-O')
        command += [str(corrupt),'--source-root',str(args.source_root),
                    '--initial-report',str(args.initial_report),'--output',str(artifact)]
        proc=subprocess.run(command,capture_output=True,text=True)
        tail=proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else ''
        need(proc.returncode!=0 and tail=='RuntimeError: '+expected and not artifact.exists(),
             'CORRUPTION_NOT_EFFECTIVELY_REJECTED:'+name)
        results.append({'name':name,'corrupt_source_sha256':sha(corrupt),
                        'nonzero_exit':True,'rejection':expected,'successful_artifact_absent':True})
    need(sha(args.checker)==checker_sha and sha(args.initial_report)==initial_sha and
         sha(Path(__file__))==own_sha,'CONTROL_INPUTS_UNCHANGED')
    report={'status':'PASS_EXECUTED_INDEPENDENT_RESPONSE_FORMULA_CORRUPTIONS',
            'subprocess_isolated_no_site_no_bytecode':True,
            'baseline_checker_sha256_before_and_after':checker_sha,
            'initial_report_sha256_before_and_after':initial_sha,
            'runner_sha256_before_and_after':own_sha,'controls':results,
            'scope':'Actual copied scripts executed; each nonzero formula mutation rejected with no successful JSON artifact. Retained-ledger reproduction and explicit nonzero witness predicates remain active under optimization.'}
    (args.output_dir/'CONTROLS.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(report['status'])


if __name__=='__main__':
    main()
