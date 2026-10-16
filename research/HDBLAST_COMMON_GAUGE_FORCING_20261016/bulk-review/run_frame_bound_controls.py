#!/usr/bin/env python3
"""Independent frame and quantitative-cap controls in both Python modes."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    p=argparse.ArgumentParser();p.add_argument('--inputs',required=True,type=Path)
    p.add_argument('--candidate',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    if a.output.exists():raise ValueError('output already exists')
    root=Path(__file__).resolve().parent
    cases=[('frame','constant_normal_connection','SPHERE_CONNECTION_DERIVATIVE_REQUIRED'),
           ('frame','missing_derivative_slot','METRIC_SECOND_JET_DERIVATIVE_SLOT_REQUIRED'),
           ('bound','metric_cap','INDEPENDENT_METRIC_CAP'),
           ('bound','scalar_cap','INDEPENDENT_SCALAR_CAP'),
           ('bound','source_cap','INDEPENDENT_SOURCE_CAP:QE')]
    rows=[]
    with tempfile.TemporaryDirectory(prefix='frame-bound-controls-',dir=a.output.resolve().parent) as tmp:
        for mode,flags in [('normal',[]),('optimized',['-O'])]:
            for family,mutation,reason in cases:
                script=root/('check_frame_connection_jets.py' if family=='frame' else 'check_independent_bulk_bounds.py')
                cmd=[sys.executable,'-I','-B']+flags+[str(script)]
                if family=='bound':cmd+=['--inputs',str(a.inputs),'--candidate',str(a.candidate)]
                out=Path(tmp)/(mode+'-'+mutation+'.json')
                proc=subprocess.run(cmd+['--inject',mutation,'--output',str(out)],text=True,capture_output=True,timeout=60)
                if proc.returncode!=2 or proc.stderr.strip()!=reason or out.exists():
                    raise ValueError('bad control:'+mode+':'+mutation+':'+proc.stderr)
                rows.append({'mode':mode,'mutation':mutation,'status':'NAMED_FAILURE_NO_OUTPUT'})
    scripts=['check_frame_connection_jets.py','check_independent_bulk_bounds.py','run_frame_bound_controls.py']
    result={'status':'PASS_FRAME_AND_BOUND_NEGATIVE_CONTROLS','count':len(rows),'controls':rows,
            'file_sha256':{name:sha256((root/name).read_bytes()).hexdigest() for name in scripts}}
    with a.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(result['status'])


if __name__=='__main__':
    try:main()
    except (ValueError,OSError,subprocess.TimeoutExpired) as e:print(str(e),file=sys.stderr);sys.exit(2)
