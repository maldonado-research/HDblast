#!/usr/bin/env python3
"""Run actual altered response-formula executions in both Python modes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

FAULTS=('absolute_error_jet_interchange','omit_frozen_initial_acceleration',
 'omit_Robin_free_wave_boundary_source','lose_W0_minus_Bs_cancellation',
 'omit_bulk_initial_acceleration','omit_original_Fx','filtered_Fx_substitution',
 'omit_complete_commutator','omit_velocity_gradient_error',
 'filtered_initial_error_forcing','omit_full_feedback_inverse')

def need(ok,message):
    if not ok: raise RuntimeError(message)

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--candidate-report',type=Path,required=True)
    ap.add_argument('--candidate-sha256',required=True)
    ap.add_argument('--oracle',type=Path,required=True)
    ap.add_argument('--oracle-sha256',required=True)
    ap.add_argument('--normal-report',type=Path,required=True)
    ap.add_argument('--optimized-report',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output_dir.exists(),'fresh absent output dir')
    me=Path(__file__).resolve()
    before={str(me):sha(me),str(args.oracle.resolve()):args.oracle_sha256,
      str(args.candidate_report.resolve()):args.candidate_sha256}
    need(sha(args.oracle)==args.oracle_sha256,'oracle authentication')
    need(sha(args.candidate_report)==args.candidate_sha256,'candidate authentication')
    normal=args.normal_report.read_bytes();optimized=args.optimized_report.read_bytes()
    need(normal==optimized,'complete normal/optimized response byte identity')
    doc=json.loads(normal)
    portable=dict(doc['input_sha256'])
    portable['run_response_controls.py']=before[str(me)]
    for role,digest in doc['input_sha256'].items():
        if role=='response_ledger_oracle.py': path=args.oracle
        elif role=='candidate_initial_report': path=args.candidate_report
        else:
            need(role.startswith('research/') and '..' not in Path(role).parts,
                 'portable source path role')
            path=args.repo_root/role
        before[str(path.resolve())]=digest
    need(all(sha(p)==h for p,h in before.items()),'all input pins before controls')
    args.output_dir.mkdir(parents=True)
    rows=[]
    for mode,flags in [('normal',[]),('optimized',['-O'])]:
        for fault in FAULTS:
            output=args.output_dir/(mode+'-'+fault+'.json')
            command=[sys.executable,*flags,'-I','-S','-B',str(args.oracle),
              '--repo-root',str(args.repo_root),'--candidate-report',str(args.candidate_report),
              '--candidate-sha256',args.candidate_sha256,'--output',str(output),'--fault',fault]
            p=subprocess.run(command,capture_output=True,text=True,check=False)
            token='CORRUPTION_REJECTED:'+fault+':nonzero_exact_residual='
            need(p.returncode!=0 and token in p.stderr and not output.exists(),
                 'effective failed execution without success artifact: '+mode+':'+fault)
            line=next(line for line in p.stderr.splitlines() if token in line)
            witness=line.split(token,1)[1]
            need(witness not in ('0','0/1'),'nonzero exact corruption witness')
            rows.append({'mode':mode,'fault':fault,'returncode':p.returncode,
              'explicit_failure_token':token,'nonzero_exact_residual':witness,
              'stdout_empty':p.stdout=='','successful_output_absent':True})
    need(all(sha(p)==h for p,h in before.items()),'all input pins unchanged on readback')
    out={'status':'PASS_RESPONSE_REPRODUCIBILITY_AND_22_REAL_FORMULA_CORRUPTIONS',
      'input_sha256':portable,'normal_optimized_response_sha256':hashlib.sha256(normal).hexdigest(),
      'normal_optimized_bytes_identical':True,'actual_child_executions':rows,
      'source_readback_unchanged':True,
      'scope':'Checks canonical response arithmetic identities and authenticated roles; interval formula/coverage controls remain in separate producer and independent evaluator.'}
    (args.output_dir/'CONTROLS.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')

if __name__=='__main__': main()
