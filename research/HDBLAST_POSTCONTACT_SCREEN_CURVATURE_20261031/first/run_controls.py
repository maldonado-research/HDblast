#!/usr/bin/env python3
"""Normal/optimized replay and real in-memory formula mutation controls."""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ORACLE=HERE/'reflected_first_oracle.py'
EXPECTED={
 'wrong_reflection':'REFLECTED_ODD_DERIVATIVE_PARITY',
 'wrong_branch':'REFLECTED_COMPARATOR_POLYNOMIAL_JET',
 'wrong_shell_sign':'PHYSICAL_TT_SHELL_SOURCE_SIGN_AND_NORMALIZATION',
 'omit_G':'PHYSICAL_TT_SHELL_SOURCE_SIGN_AND_NORMALIZATION',
 'omit_dynamic_lift':'dynamic boundary third-jet retained',
}

def need(ok,message):
 if not ok:raise ValueError(message)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 pins={p.name:sha(p) for p in (Path(__file__),ORACLE)}
 rows=[]
 for optimized in (False,True):
  for mutation,reason in EXPECTED.items():
   tag=('OPTIMIZED' if optimized else 'NORMAL')+'_'+mutation
   out=HERE/('CONTROL_'+tag+'.json')
   need(not out.exists(),'FRESH_CONTROL_OUTPUT')
   cmd=[sys.executable,'-B']+(['-O'] if optimized else [])
   cmd += [str(ORACLE),'--inject-error',mutation,'--output',str(out)]
   proc=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=90)
   need(proc.returncode!=0 and reason in proc.stderr and not out.exists(),
        'CONTROL_REJECTION:'+tag)
   rows.append({'mutation':mutation,'optimized':optimized,'returncode':proc.returncode,
                'expected_stderr':reason,'success_certificate_absent':True,
                'stdout_sha256':hashlib.sha256(proc.stdout.encode()).hexdigest(),
                'stderr_sha256':hashlib.sha256(proc.stderr.encode()).hexdigest()})
 normal=HERE/'FIRST_NORMAL.json';optimized=HERE/'FIRST_OPTIMIZED.json'
 need(normal.is_file() and optimized.is_file() and normal.read_bytes()==optimized.read_bytes(),
      'BYTE_IDENTICAL_NORMAL_OPTIMIZED_REPORTS')
 need(pins=={p.name:sha(p) for p in (Path(__file__),ORACLE)},'LOCAL_SOURCE_CHANGED')
 result={'status':'PASS_REFLECTED_ACTUAL_FIRST_ORACLE_CONTROLS',
         'normal_optimized_report_sha256':sha(normal),'byte_identical':True,
         'formula_mutant_subprocesses':rows,'source_sha256':pins,
         'in_memory_mutations':'Only authenticated source text compiled in isolated control subprocesses is altered; frozen source files are untouched.'}
 out=HERE/'CONTROLS.json';need(not out.exists(),'FRESH_CONTROL_REPORT')
 with out.open('x') as handle:json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
 print(result['status'])

if __name__=='__main__':main()
