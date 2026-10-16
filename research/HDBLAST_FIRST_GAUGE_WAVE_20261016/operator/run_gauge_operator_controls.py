#!/usr/bin/env python3
"""Exercise real CLI rejections in normal and optimized Python; stdlib only."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


CASES = {
    'source_sign':'longitudinal_radial_defect',
    'radial_drift':'inward_radial_drift',
    'tangential_drift':'inward_tangential_drift',
    'angular_curvature':'desitter_angular_component',
    'missing_tangential_coupling':'warped_tangential_component',
    'reflection_drift':'reflection_B_parity',
    'omit_master_time_factor':'reconstruction_source',
    'contraction_underbound':'contractive_cap',
}


def require(condition,message):
    if not condition:
        raise ValueError(message)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    require(not args.output.exists(),'output path already exists')
    root=Path(__file__).resolve().parent
    checker=root/'check_gauge_operator.py'
    results=[]
    with tempfile.TemporaryDirectory(prefix='operator-controls-',dir=args.output.resolve().parent) as tmp:
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            base=[sys.executable,'-I','-B']+flags+[str(checker)]
            for mutation,identity in CASES.items():
                target=Path(tmp)/(mode+'-'+mutation+'.json')
                proc=subprocess.run(base+['--output',str(target),'--inject',mutation],
                                    text=True,capture_output=True,timeout=30)
                require(proc.returncode==2,'wrong failure code:'+mutation)
                require(proc.stderr.strip()=='IDENTITY_FAILED:'+identity,
                        'wrong failure reason:'+mutation+':'+proc.stderr)
                require(not target.exists(),'mutation wrote output:'+mutation)
                results.append({'mode':mode,'control':mutation,'status':'EXPECTED_NAMED_FAILURE_NO_OUTPUT'})
            target=Path(tmp)/(mode+'-existing.json')
            target.write_text('KEEP EXISTING OUTPUT\n')
            proc=subprocess.run(base+['--output',str(target)],text=True,capture_output=True,timeout=30)
            require(proc.returncode==2 and proc.stderr.strip()=='output path already exists',
                    'existing-output rejection failed')
            require(target.read_text()=='KEEP EXISTING OUTPUT\n','existing output changed')
            results.append({'mode':mode,'control':'existing_output','status':'EXPECTED_REJECTION_PRESERVED_OUTPUT'})
    result={'status':'PASS_ALL_NORMAL_AND_OPTIMIZED_CLI_CONTROLS','control_count':len(results),
            'checker_sha256':hashlib.sha256(checker.read_bytes()).hexdigest(),
            'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'controls':results}
    with args.output.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True)
        handle.write('\n')
    print(result['status'])


if __name__=='__main__':
    try:
        main()
    except (ValueError,OSError,subprocess.TimeoutExpired) as exc:
        print(str(exc),file=sys.stderr)
        sys.exit(2)
