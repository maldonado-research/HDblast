#!/usr/bin/env python3
"""Run semantic higher-jet corruption controls in ordinary and optimized Python."""
import argparse
from pathlib import Path
import subprocess
import sys
import json

ROOT=Path(__file__).resolve().parent
CONTROLS={
 'use_Ax0_initial_velocity':'second-time initial velocity uses A(0), not A_x(0)',
 'drop_second_commutator_two':'second-time commutator factor two',
 'drop_field_value_factor':'weighted field value factor is X/2=3/2',
 'omit_source_radial_derivative':'third radial recovery retains source radial derivative',
 'wrong_physical_scaling':'unscaled derivatives require ell to negative total order',
}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();destination=args.output.resolve()
    if destination==ROOT or ROOT in destination.parents or destination.exists():
        raise ValueError('control output must be new and outside higher-norm input tree')
    runs=[]
    for optimized in (False,True):
        for name,reason in CONTROLS.items():
            cmd=[sys.executable]+(['-O'] if optimized else [])+[
                str(ROOT/'check_gauge_higher_norms.py'),'--inject-error',name,
                '--output',str(destination.parent/('rejected_'+name+str(optimized)+'.json'))]
            result=subprocess.run(cmd,text=True,capture_output=True)
            if result.returncode!=1 or reason not in result.stderr:
                raise ValueError('corruption not rejected for intended reason: '+name+'\n'+result.stderr)
            runs.append({'control':name,'optimized':optimized,'exit_code':result.returncode,
                         'expected_exception':reason,'status':'REJECTED_FOR_INTENDED_REASON'})
    report={'status':'PASS_HIGHER_NORM_SEMANTIC_CONTROLS_NORMAL_AND_OPTIMIZED',
            'runs':runs,'scope':'Initial A(0), commutator2, field-value3/2, radial source derivative and ell inverse powers.'}
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'runs':len(runs),'output':str(destination)},indent=2))


if __name__=='__main__':main()
