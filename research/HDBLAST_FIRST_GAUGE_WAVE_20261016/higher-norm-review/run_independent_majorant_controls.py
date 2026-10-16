#!/usr/bin/env python3
"""Real mutation rejections, both interpreter modes, with no output on failure."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--coefficients',required=True,type=Path)
    parser.add_argument('--candidate',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    a=parser.parse_args()
    if a.output.exists():raise ValueError('output already exists')
    script=Path(__file__).resolve().with_name('check_independent_higher_majorants.py')
    cases={'initial_velocity':'WRONG_SECOND_TIME_INITIAL_VELOCITY',
           'commutator':'WRONG_SECOND_COMMUTATOR_MULTIPLICITY',
           'radial_recovery':'WRONG_THIRD_RADIAL_RECOVERY'}
    results=[]
    with tempfile.TemporaryDirectory(prefix='majorant-controls-',dir=a.output.resolve().parent) as tmp:
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            base=[sys.executable,'-I','-B']+flags+[str(script),'--coefficients',str(a.coefficients),'--candidate',str(a.candidate)]
            for mutation,failure in cases.items():
                dest=Path(tmp)/(mode+'-'+mutation+'.json')
                proc=subprocess.run(base+['--inject',mutation,'--output',str(dest)],text=True,capture_output=True,timeout=60)
                if proc.returncode!=2 or proc.stderr.strip()!=failure or dest.exists():
                    raise ValueError('control failed:'+mode+':'+mutation+':'+proc.stderr)
                results.append({'mode':mode,'mutation':mutation,'status':'EXPECTED_FAILURE_NO_OUTPUT'})
    report={'status':'PASS_INDEPENDENT_MAJORANT_CLI_CONTROLS','count':len(results),
            'checker_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),
            'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'controls':results}
    with a.output.open('x') as f:json.dump(report,f,indent=2,sort_keys=True);f.write('\n')
    print(report['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.TimeoutExpired) as e:
        print(str(e),file=sys.stderr);sys.exit(2)
