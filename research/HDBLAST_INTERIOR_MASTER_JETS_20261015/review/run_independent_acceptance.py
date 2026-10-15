#!/usr/bin/env python3
"""Portable replay of reviewed scientific programs and independent checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok,why):
    if not ok:raise ValueError(why)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    root=args.root.resolve();out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    own_checker=(Path(__file__).resolve().parent/'check_independent_algebra.py').relative_to(root).as_posix()
    # Relative script and report paths are serialized; absolute paths are only
    # execution selectors inside the temporary/replay workspace.
    jobs=[
        ('POTENTIAL_CONTROLS',False,'potential/run_potential_controls.py',[]),
        ('JETS_NORMAL',False,'jets/check_higher_jets.py',['--root',str(root)]),
        ('JETS_OPTIMIZED',True,'jets/check_higher_jets.py',['--root',str(root)]),
        ('JETS_CONTROLS',False,'jets/check_higher_jets.py',['--run-controls']),
        ('BULK_NORMAL',False,'bulk/check_longitudinal_sources.py',['--run-controls']),
        ('BULK_OPTIMIZED',True,'bulk/check_longitudinal_sources.py',['--run-controls']),
        ('INDEPENDENT_NORMAL',False,own_checker,['--root',str(root)]),
        ('INDEPENDENT_OPTIMIZED',True,own_checker,['--root',str(root)]),
    ]
    def execute(job):
        name,optimized,script,extra=job
        target=out/(name+'.json')
        command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(root/script),*extra,'--output',str(target)]
        result=subprocess.run(command,cwd=root,text=True,capture_output=True,timeout=180)
        need(result.returncode==0,name+' failed: '+result.stderr)
        need(target.is_file(),name+' did not create receipt')
        return {'name':name,'script':script,'optimized':optimized,'returncode':result.returncode,
                'output':target.name,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
    with ThreadPoolExecutor(max_workers=4) as pool:processes=list(pool.map(execute,jobs))
    def read(name):return json.loads((out/(name+'.json')).read_text())
    golden=[('JETS_NORMAL','jets/HIGHER_JETS_NORMAL.json'),('JETS_OPTIMIZED','jets/HIGHER_JETS_OPTIMIZED.json'),
            ('BULK_NORMAL','bulk/ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS.json'),
            ('BULK_OPTIMIZED','bulk/ACTUAL_LONGITUDINAL_BULK_SOURCE_BOUNDS_OPTIMIZED.json'),
            ('JETS_CONTROLS','jets/HIGHER_JETS_CONTROLS.json')]
    for fresh,relative in golden:
        need((out/(fresh+'.json')).read_bytes()==(root/relative).read_bytes(),'fresh golden differs: '+fresh)
    for prefix in ('JETS','BULK','INDEPENDENT'):
        a=read(prefix+'_NORMAL');b=read(prefix+'_OPTIMIZED')
        if prefix=='BULK':a['controls'].pop('optimized');b['controls'].pop('optimized')
        else:a.pop('optimized');b.pop('optimized')
        need(a==b,'mode science differs: '+prefix)
    reasons={'half_ricci':'longitudinal Einstein Hessian (0, 0): exact residual ',
             'drop_mixed_radial':'longitudinal Einstein Hessian (0, 1): exact residual ',
             'wrong_scalar_wave_factor':'on-shell scalar cancellation: exact residual ',
             'omit_physical_scale':'physical order-1 derivative scale: exact residual '}
    for mode in ('NORMAL','OPTIMIZED'):
        controls=read('BULK_'+mode)['controls']['cases']
        need(len(controls)==4,'four bulk controls')
        for row in controls:
            need(row['returncode']==1 and row['rejection'].startswith(reasons[row['mutation']]),'wrong bulk control rejection reason')
    pot=read('POTENTIAL_CONTROLS')
    need(pot['normal_optimized_byte_identity'] is True,'potential mode identity')
    need(sum(x.get('status')=='REJECTED_AS_REQUIRED' for x in pot['cases'].values())==8,'eight potential controls')
    need(read('JETS_CONTROLS')['control_runs']==22,'twenty-two jet controls')
    report={'status':'PASS_REVIEWED_SOURCES_AND_INDEPENDENT_REPLAY','processes':processes,
            'scientific_positive_runs_by_family':{'potential':2,'jets':2,'bulk':2,'independent_referee':2},
            'expected_failure_subprocesses_by_family':{'potential':8,'jets':22,'bulk':8},
            'golden_byte_comparisons':len(golden),'science_equal_normal_optimized':True,
            'scope':'Separate process and witness counts; no aggregate mathematical theorem count is inferred.'}
    (out/'REPLAY_ACCEPTANCE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'scientific_positive_runs_by_family':report['scientific_positive_runs_by_family'],
                      'expected_failure_subprocesses_by_family':report['expected_failure_subprocesses_by_family']}))


if __name__=='__main__':main()
