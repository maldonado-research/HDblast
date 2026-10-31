#!/usr/bin/env python3
"""Real record-corruption rejection probes for independent report replay."""
from concurrent.futures import ThreadPoolExecutor
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
INPUT=Path('/workspace/hdblast-tt-postcontact-31d/POSTCONTACT_NORMAL.json')


def need(ok,message):
    if not ok:raise ValueError(message)


def main():
    dest=HERE/'report-controls';need(not dest.exists(),'FRESH_CONTROLS_DIRECTORY');dest.mkdir()
    original=json.loads(INPUT.read_text())
    cases=[]
    def add(name,tag,change):
        d=copy.deepcopy(original);change(d)
        path=dest/(name+'.json');path.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
        cases.append((name,tag,path,hashlib.sha256(path.read_bytes()).hexdigest()))
    add('omitted_nonzero_G','GENERIC_NONZERO_ENDPOINT_HISTORY',lambda d:d.update(shell_endpoint_G=['0','0']))
    add('reversed_reflection_radius','REFLECTED_CHARACTERISTIC_GEOMETRY',
        lambda d:d['characteristic_segments'][0].update(reflected_r=['-9/8','-1']))
    add('removed_full_feedback_tail','COMPLETE_FEEDBACK_TAIL',lambda d:d.update(complete_feedback_Ux_error='0'))
    add('removed_fixed_quadratic','FULL_FIRST_FIELD_QUADRATIC_COEFFICIENT',lambda d:d.update(screen_quadratic=['0','0']))
    add('dropped_uniform_source_cell','BULK_COVER_COUNT',lambda d:d['bulk_cap_cover'].pop())
    add('wrong_physical_clock','PHYSICAL_CLOCK_ELL',lambda d:d['apex_background'].update(ell='1/1000'))
    add('claimed_positive_sign','EXACT_SIGN_ACCEPTANCE',lambda d:d.update(sign='positive'))
    jobs=[]
    for name,tag,path,pin in cases:
        for optimized in (False,True):
            out=dest/(name+('_optimized' if optimized else '_normal')+'_report.json')
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_report.py'),
              '--report',str(path),'--report-sha256',pin,'--output',str(out)]
            jobs.append((name,tag,optimized,cmd,out))
    def execute(job):
        name,tag,optimized,cmd,out=job
        run=subprocess.run(cmd,capture_output=True,text=True,timeout=20)
        need(run.returncode!=0 and not out.exists() and tag in run.stderr,
          'CONTROL_ACCEPTED_OR_WRONG_FAILURE:'+name+':'+run.stderr[-400:])
        return {'case':name,'optimized':optimized,'exit_code':run.returncode,
          'expected_rejection':tag,'success_report_created':out.exists(),
          'stderr_last_line':run.stderr.strip().splitlines()[-1]}
    with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(execute,jobs))
    need((HERE/'REPORT_NORMAL.json').read_bytes()==(HERE/'REPORT_OPTIMIZED.json').read_bytes(),
         'POSITIVE_NORMAL_OPTIMIZED_MISMATCH')
    result={'status':'ALL_REAL_POSTCONTACT_REPORT_CORRUPTIONS_REJECTED','subprocesses':len(rows),
      'cases':rows,'positive_normal_optimized_byte_identical':True,
      'positive_report_sha256':hashlib.sha256((HERE/'REPORT_NORMAL.json').read_bytes()).hexdigest(),
      'scope':'Mutated report inputs are deliberately repinned to reach their arithmetic/geometry guard; the independent checker source is unchanged.'}
    (dest/'CONTROLS.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])


if __name__=='__main__':main()
