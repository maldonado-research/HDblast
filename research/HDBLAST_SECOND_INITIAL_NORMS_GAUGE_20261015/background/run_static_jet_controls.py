#!/usr/bin/env python3
"""Normal/-O static-jet replay and ten intended-error controls."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent


def require(ok,text):
 if not ok:raise ValueError(text)


def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--output',type=Path,default=ROOT/'STATIC_JET_CONTROL_RESULTS.json')
 args=ap.parse_args()
 cases={'input_tamper':'actual retained input pin mismatch',
        'wrong_jet_normalization':'formal jet identity failed: rho order 2',
        'flip_rho_ODE':'formal jet identity failed: rho order 1',
        'halve_B_log_derivative':'formal jet identity failed: B order 1',
        'unsafe_p_denominator':'actual annulus denominator signs failed'}
 results={};positive=[]
 with tempfile.TemporaryDirectory(prefix='actual-static-jet-controls-') as tmp:
  for optimized in (False,True):
   mode='optimized' if optimized else 'normal'
   command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])
   output=Path(tmp)/(mode+'.json')
   p=subprocess.run(command+[str(ROOT/'bound_static_coefficient_jets.py'),'--output',str(output)],capture_output=True,text=True,timeout=180)
   require(p.returncode==0,mode+' positive failed: '+p.stderr)
   positive.append(output.read_bytes())
   results[mode]={'status':'PASS','sha256':sha256(positive[-1]).hexdigest()}
   for mutation,reason in cases.items():
    target=Path(tmp)/(mode+'-'+mutation+'.json')
    bad=subprocess.run(command+[str(ROOT/'bound_static_coefficient_jets.py'),'--inject-error',mutation,'--output',str(target)],capture_output=True,text=True,timeout=180)
    require(bad.returncode!=0 and ('ValueError: '+reason) in bad.stderr and not target.exists(),mode+' control failed intended rejection: '+mutation)
    results[mode+'_'+mutation]={'status':'REJECTED_AS_REQUIRED','exception':reason,'exit_code':bad.returncode,'no_output_report':True}
 require(positive[0]==positive[1],'normal/optimized exact report bytes differ')
 require(positive[0]==(ROOT/'ACTUAL_ANNULUS_STATIC_JETS.json').read_bytes(),'fresh replay differs from retained exact report')
 report={'status':'PASS_STATIC_JET_NORMAL_OPTIMIZED_AND_TEN_NEGATIVE_CONTROLS',
         'normal_optimized_byte_identity':True,'cases':results}
 args.output.write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(report,indent=2))


if __name__=='__main__':main()
