#!/usr/bin/env python3
"""Replay normal/optimized exact audit and require purposeful failures."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent


def require(ok, text):
    if not ok:
        raise ValueError(text)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'POTENTIAL_CONTROL_RESULTS.json')
    args=parser.parse_args()
    cases={
        'input_identity':'input authentication failed',
        'drop_Wzz_cross_factor':'formal identity failed: Wzz_total_derivative',
        'drop_Qzz_reflection_cross':'formal identity failed: Qzz_total_derivative',
        'undersized_Wzz':'claimed strict bound failed: absolute_W_zz',
    }
    results={}
    positive=[]
    with tempfile.TemporaryDirectory(prefix='potential-controls-') as tmp:
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal'
            base=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])
            output=Path(tmp)/(mode+'.json')
            process=subprocess.run(base+[str(ROOT/'check_actual_potential.py'),'--output',str(output)],capture_output=True,text=True)
            require(process.returncode==0,mode+' positive run failed: '+process.stderr)
            positive.append(output.read_bytes())
            results[mode]={'status':'PASS','sha256':sha256(positive[-1]).hexdigest()}
            for name,expected in cases.items():
                target=Path(tmp)/(mode+'-'+name+'.json')
                bad=subprocess.run(base+[str(ROOT/'check_actual_potential.py'),'--corrupt',name,'--output',str(target)],capture_output=True,text=True)
                require(bad.returncode!=0 and expected in bad.stderr and not target.exists(),mode+' failed to reject '+name)
                results[mode+'_'+name]={'status':'REJECTED_AS_REQUIRED','exception':expected,'exit_code':bad.returncode}
    require(positive[0]==positive[1],'normal/optimized report bytes differ')
    require(positive[0]==(ROOT/'ACTUAL_POTENTIAL_JET_BOUNDS.json').read_bytes(),'retained report differs from fresh controls')
    report={'status':'PASS_NORMAL_OPTIMIZED_AND_EIGHT_NEGATIVE_CONTROLS','normal_optimized_byte_identity':True,'cases':results}
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
