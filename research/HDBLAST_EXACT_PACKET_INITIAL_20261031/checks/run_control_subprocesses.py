#!/usr/bin/env python3
"""Runtime receipts for explicit rejection in normal and optimized Python."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True

CHECK_PIN='6630a8cdde8af9821a99a8d17b10916125b637537cd375b1d76e52f05488044f'
INITIAL_PIN='692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3'
CASES=('restored-cutoff','restored-convolution','midpoint-exponential','wrong-radial-scale',
       'zeroed-mean-override','omitted-background-channel','omitted-remainder-channel',
       'velocity-absolute-error-interchange','omitted-velocity-charge','global-denominator-at-local-point')

def need(ok,reason):
    if not ok: raise RuntimeError(reason)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--initial-report',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    need(not args.output_dir.exists(),'runtime receipt directory must be fresh')
    script=Path(__file__).resolve().parent/'check_controls.py'
    need(sha(script)==CHECK_PIN and sha(args.initial_report)==INITIAL_PIN,'authenticate executed checker and canonical report before subprocesses')
    self_before=sha(Path(__file__))
    receipts=[]
    for optimized in (False,True):
        for case in CASES:
            output=args.output_dir/('optimized' if optimized else 'normal')/(case+'.json')
            cmd=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(script),'--repo-root',str(args.repo_root),
                '--initial-report',str(args.initial_report),'--output',str(output),'--inject-error',case]
            result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=False)
            marker='CONTROL_REJECTED '+case+' '
            need(result.returncode!=0 and marker in result.stderr and not output.exists(),
                 'effective failure/no artifact required for '+case)
            witness_line=next(line.split(marker,1)[1] for line in result.stderr.splitlines() if marker in line)
            receipts.append({'case':case,'optimized':optimized,'returncode':result.returncode,
                'successful_output_exists':False,'explicit_witness':json.loads(witness_line),'stdout':result.stdout})
    need(sha(script)==CHECK_PIN and sha(args.initial_report)==INITIAL_PIN and sha(Path(__file__))==self_before,
         'checker report and runtime driver unchanged after subprocesses')
    args.output_dir.mkdir(parents=True)
    result={'status':'PASS_20_EFFECTIVE_REJECTION_SUBPROCESSES','checker_sha256':CHECK_PIN,
        'initial_sha256':INITIAL_PIN,'driver_sha256':self_before,'runtime_receipts':receipts}
    (args.output_dir/'RECEIPTS.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'subprocesses':len(receipts)},sort_keys=True))

if __name__=='__main__':main()
