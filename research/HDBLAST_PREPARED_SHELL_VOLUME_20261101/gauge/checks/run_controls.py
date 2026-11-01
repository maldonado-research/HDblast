#!/usr/bin/env python3
"""Real-process fail-closed controls; exact expected witnesses, never substring tests."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
from authenticate import authenticate, output_path, COMPONENT, failure_payload


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repo-root',required=True,type=Path)
    parser.add_argument('--output',required=True);args=parser.parse_args()
    try:
        destination=output_path(args.output,args.repo_root)
        before,pins=authenticate(args.repo_root)
        records=[];calls=0
        with tempfile.TemporaryDirectory(prefix='prepared-gauge-controls-') as name:
            temp=Path(name);component=temp/'component';fixture=temp/'repo'
            for relative in before['local_files']+['gauge/ESSENTIAL_PINS.json']:
                target=component/relative;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(COMPONENT/relative,target)
            for row in pins['sources']:
                target=fixture/row['path'];target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(args.repo_root/row['path'],target)
            def invoke(kind,label,mutation=None,expected=None):
                nonlocal calls
                out=temp/(label+'.json')
                command=[sys.executable,'-I','-B']
                if sys.flags.optimize:command.append('-O')
                command += [str(component/'gauge/checks'/('check_'+kind+'.py')),
                            '--repo-root',str(fixture),'--output',str(out)]
                if mutation:command += ['--mutation',mutation]
                result=subprocess.run(command,capture_output=True,text=True);calls+=1
                if expected is None:
                    if result.returncode!=0 or result.stdout!='' or result.stderr!='' or not out.is_file():
                        raise RuntimeError('baseline failure '+label)
                    report=json.loads(out.read_text())
                    if report['status']!='PASS':raise RuntimeError('baseline status '+label)
                    records.append({'case':label,'expected':'PASS','observed':'PASS'})
                else:
                    witness={'returncode':result.returncode,'stdout':result.stdout,
                             'stderr':result.stderr,'output_exists':out.exists()}
                    wanted={'returncode':2,'stdout':json.dumps({'status':'FAIL','witness':expected},sort_keys=True)+'\n',
                            'stderr':'','output_exists':False}
                    if witness!=wanted:raise RuntimeError('exact failure witness mismatch '+label+': '+repr(witness))
                    records.append({'case':label,'expected_witness':wanted,'observed_witness':witness})
            for kind in ('preparation','initial_mean','conditioned'):invoke(kind,'baseline_'+kind)
            for mutation,expected in (
                ('omit_radial_second_profile','INITIAL_GAUGE_RESIDUAL_OMIT_RADIAL_SECOND_PROFILE'),
                ('omit_trace_inverse_cross','INITIAL_GAUGE_RESIDUAL_OMIT_TRACE_INVERSE_CROSS'),
                ('omit_metric_mixed_cross','INITIAL_GAUGE_RESIDUAL_OMIT_METRIC_MIXED_CROSS'),
                ('omit_lapse_square_value','INITIAL_GAUGE_RESIDUAL_OMIT_LAPSE_SQUARE_VALUE'),
                ('halve_second_gauge_velocity','INITIAL_GAUGE_RESIDUAL_HALVE_SECOND_GAUGE_VELOCITY')):
                invoke('initial_mean',mutation,mutation,expected)
            for mutation,expected in [('omit_time_jacobian','CONDITIONED_LAPSE_VALUE'),
                                      ('omit_vector_potential','CONDITIONED_MIXED_VELOCITY'),
                                      ('omit_scalar_compensation','CONDITIONED_PULLBACK')]:
                invoke('conditioned',mutation,mutation,expected)
            invoke('preparation','wrong_momentum','wrong_momentum','PREPARATION_MOMENTUM')
            source=fixture/pins['sources'][0]['path'];saved=source.read_bytes();source.write_bytes(saved+b'\n')
            invoke('initial_mean','changed_upstream_proof',expected='SOURCE_AUTHENTICATION')
            source.write_bytes(saved)
            local=component/'gauge/conditioned_mean_initial.py';saved=local.read_bytes();local.write_bytes(saved+b'\n')
            invoke('conditioned','changed_local_api',expected='LOCAL_AUTHENTICATION')
            local.write_bytes(saved)
        after,_=authenticate(args.repo_root)
        if before!=after:raise RuntimeError('authentication changed during controls')
        report={'status':'PASS','authentication':{'before':before,'after':after,'equal':True},
                'real_subprocess_count':calls,'positive_subprocess_count':3,'negative_subprocess_count':calls-3,
                'in_process_mutant_count':0,'exact_witness_checks':records}
        payload=json.dumps(report,indent=2)+'\n'
        if destination is None:sys.stdout.write(payload)
        else:
            with destination.open('x') as stream:stream.write(payload)
    except Exception as exc:
        sys.stdout.write(failure_payload(exc));return 2
    return 0


if __name__=='__main__':raise SystemExit(main())
