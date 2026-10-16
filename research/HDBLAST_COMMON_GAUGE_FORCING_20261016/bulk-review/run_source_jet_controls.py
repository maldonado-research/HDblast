#!/usr/bin/env python3
"""Run specific source-formula mutations in normal and optimized interpreters."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

CASES={'connection_factor':'CONNECTION_HESSIAN',
       'background_connection_derivative':'COVARIANT_CONNECTION_DERIVATIVE',
       'ricci_product_factor':'NATURAL_EINSTEIN_HESSIAN',
       'stress_factor':'NATURAL_EINSTEIN_HESSIAN',
       'potential_cross_factor':'NATURAL_EINSTEIN_HESSIAN',
       'scalar_connection':'SCALAR_HESSIAN',
       'harmonic_quadratic_contraction':'QUADRATIC_GAUGE_VECTOR',
       'reduction_symmetrization':'REDUCED_EINSTEIN_HESSIAN',
       'reduction_lowering_cross':'REDUCED_EINSTEIN_HESSIAN'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args()
    if a.output.exists():raise ValueError('output already exists')
    checker=Path(__file__).resolve().with_name('check_common_gauge_sources.py')
    rows=[]
    with tempfile.TemporaryDirectory(prefix='source-controls-',dir=a.output.resolve().parent) as tmp:
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            base=[sys.executable,'-I','-B']+flags+[str(checker)]
            for mutation,reason in CASES.items():
                dest=Path(tmp)/(mode+'-'+mutation+'.json')
                p=subprocess.run(base+['--inject',mutation,'--output',str(dest)],text=True,capture_output=True,timeout=60)
                if p.returncode!=2 or p.stderr.strip()!=reason or dest.exists():
                    raise ValueError('control failed:'+mutation+':'+p.stderr)
                rows.append({'mode':mode,'mutation':mutation,'status':'NAMED_FAILURE_NO_OUTPUT'})
    result={'status':'PASS_SOURCE_JET_NEGATIVE_CONTROLS','count':len(rows),'controls':rows,
            'checker_sha256':sha256(checker.read_bytes()).hexdigest(),
            'runner_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    with a.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(result['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.TimeoutExpired) as e:
        print(str(e),file=sys.stderr);sys.exit(2)
