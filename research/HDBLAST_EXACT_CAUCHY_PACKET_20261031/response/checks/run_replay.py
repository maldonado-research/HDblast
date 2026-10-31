#!/usr/bin/env python3
"""Rebuild canonical rational reports with normal/-O Python and disk controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--absolute-dir',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    base=Path(__file__).resolve().parent.parent
    out=args.output_dir
    need(not out.exists(),'fresh replay output directory')
    out.mkdir(parents=True)
    reference=json.loads((base/'reference/RESPONSE.json').read_text())
    pins=reference['selected_input_sha256_before']
    adapter=base/'checks/response_adapter.py'; freeze=base/'checks/FREEZE.json'
    need(sha(adapter)==pins['response_adapter.py'],'canonical adapter source pin')
    need(sha(freeze)==pins['FREEZE.json'],'canonical freeze pin')
    cap=args.absolute_dir/'REFINEMENT.json'
    need(sha(cap)==pins['refined_caps_report'],'canonical absolute report pin')
    common=[str(adapter),'--repo-root',str(args.repo_root),'--freeze-sha256',pins['FREEZE.json']]
    refined=['--caps',str(cap),'--caps-sha256',pins['refined_caps_report']]
    for filename in ('refine.py','check_refinement.py','auth.py','FREEZE.json'):
        p=args.absolute_dir/filename; expected=pins['proof-input:'+filename]
        need(sha(p)==expected,'absolute proof source pin: '+filename)
        refined += ['--proof-input',str(p)+'='+expected]
    reports={}
    for label,flags in [('BASELINE',[]),('RESPONSE',refined)]:
        outputs=[]
        for mode,opts in [('normal',[]),('optimized',['-O'])]:
            p=out/(label+'_'+mode+'.json')
            subprocess.run([sys.executable]+opts+common+flags+['--output',str(p)],check=True)
            outputs.append(p)
        need(outputs[0].read_bytes()==outputs[1].read_bytes(),'normal/optimized byte identity: '+label)
        need(outputs[0].read_bytes()==(base/'reference'/(label+'.json')).read_bytes(),'canonical reference exact bytes: '+label)
        reports[label]=sha(outputs[0])
    subprocess.run([sys.executable,str(base/'checks/run_controls.py'),'--repo-root',str(args.repo_root),
        '--caps',str(cap),'--response',str(out/'RESPONSE_normal.json'),
        '--output-dir',str(out/'controls')],check=True)
    controls=json.loads((out/'controls/CONTROLS.json').read_text())
    need(controls['count']==62,'complete frozen response corruption count')
    result={'status':'PASS_PORTABLE_SAME_FUNCTION_RESPONSE_REFINEMENT',
        'canonical_normal_optimized_bytes_identical':True,'reports':reports,
        'rejected_corrupted_subprocesses':62,'runner_sha256':sha(__file__),
        'same_function_true_error_unchanged':True,'conditional_physical_side_estimates':True}
    (out/'REPLAY.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result))

if __name__=='__main__':
    main()
