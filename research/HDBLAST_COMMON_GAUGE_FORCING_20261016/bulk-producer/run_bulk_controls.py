#!/usr/bin/env python3
"""Reason-specific common-gauge source controls in ordinary and optimized Python."""
import argparse
from pathlib import Path
import subprocess
import sys
import json

ROOT=Path(__file__).resolve().parent
CONTROLS={
 'freeze_frame_scales':'radial frame scale derivative retained',
 'drop_derivative_slot':'second covariant jet includes derivative-slot connection',
 'omit_scalar_pullback':'metric and scalar use the same pullback plus pf',
 'undersize_sphere_fourth':'actual quadrupole fourth sphere derivative cap',
 'drop_gauge_double_sum':'quadratic gauge row retains both contracted indices',
 'drop_scalar_hessian_connection':'scalar Hessian retains its frame connection',
 'half_scalar_stress_hessian':'second parameter derivative retains scalar stress factor two',
}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path)
    args=p.parse_args();out=args.output.resolve()
    if out==ROOT or ROOT in out.parents or out.exists():raise ValueError('new external output required')
    rows=[]
    for opt in (False,True):
        for name,reason in CONTROLS.items():
            result=subprocess.run([sys.executable]+(['-O'] if opt else [])+[
                str(ROOT/'check_common_gauge_bulk.py'),'--inject-error',name,
                '--output',str(out.parent/('REJECTED_'+name+str(opt)+'.json'))],text=True,capture_output=True)
            if result.returncode!=1 or reason not in result.stderr:
                raise ValueError('wrong rejection: '+name+'\n'+result.stderr)
            rows.append({'mutation':name,'optimized':opt,'exit_code':result.returncode,
                         'expected_exception':reason,'report_written':False})
    report={'status':'PASS_ALL_COMMON_GAUGE_SEMANTIC_CONTROLS_NORMAL_AND_OPTIMIZED','runs':rows}
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'runs':len(rows),'output':str(out)},indent=2))


if __name__=='__main__':main()
