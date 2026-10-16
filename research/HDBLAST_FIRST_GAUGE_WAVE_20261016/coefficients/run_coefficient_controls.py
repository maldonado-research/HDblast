#!/usr/bin/env python3
"""Run all explicit semantic corruption controls in normal and optimized modes."""
import argparse
from pathlib import Path
import subprocess
import sys
import json

ROOT=Path(__file__).resolve().parent
CONTROLS={
 'wrong_jet_factorial':'formal normalized jet identity',
 'wrong_source_sign':'proper-y/scaled-r source identity',
 'annulus_alpha_reuse':'alpha integral must cover the whole collar from zero',
 'wrong_alpha_origin':'alpha(0)=1 lies in global enclosure',
 'omit_ell_A':'weighted C1 Volterra matrix contraction',
}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    destination=args.output.resolve()
    if destination==ROOT or ROOT in destination.parents or destination.exists():
        raise ValueError('control output must be new and outside coefficient input tree')
    runs=[]
    for optimized in (False,True):
        for name,message in CONTROLS.items():
            cmd=[sys.executable]+(['-O'] if optimized else [])+[
                str(ROOT/'check_global_gauge_coefficients.py'),'--inject-error',name,
                '--output',str(destination.parent/('CONTROL_SHOULD_NOT_EXIST_'+name+str(optimized)+'.json'))]
            result=subprocess.run(cmd,text=True,capture_output=True)
            if result.returncode==0 or message not in result.stderr:
                raise ValueError('control was not rejected for intended reason: '+name+'\n'+result.stderr)
            runs.append({'control':name,'optimized':optimized,'exit_code':result.returncode,
                         'expected_exception':message,'status':'REJECTED_FOR_INTENDED_REASON'})
    report={'status':'PASS_ALL_SEMANTIC_CONTROLS_REJECTED_NORMAL_AND_OPTIMIZED',
            'runs':runs,'scope':'Wrong normalization, sign, domain, origin and scaling claims fail with active exception guards.'}
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'runs':len(runs),'output':str(args.output)},indent=2))


if __name__=='__main__':main()
