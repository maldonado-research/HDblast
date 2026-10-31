#!/usr/bin/env python3
"""Portable normal/optimized absolute-cap refinement replay and fault controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args();root=args.repo_root.resolve();out=args.output_dir.resolve()
    need(not out.exists(),'fresh absent output directory required');out.mkdir(parents=True)
    pairs={}
    def run(command,name):
        p=subprocess.run(command,capture_output=True,timeout=600)
        (out/(name+'.stdout')).write_bytes(p.stdout);(out/(name+'.stderr')).write_bytes(p.stderr)
        need(p.returncode==0,'replay '+name+' failed: '+p.stderr.decode(errors='replace')[-2000:])
    for optimized in (False,True):
        mode='optimized' if optimized else 'normal';flags=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])
        report=out/(mode+'-REFINEMENT.json');review=out/(mode+'-CHECK.json')
        run(flags+[str(HERE/'refine.py'),'--repo-root',str(root),'--output',str(report)],mode+'-producer')
        need(report.read_bytes()==(HERE/'REFINEMENT.json').read_bytes(),'canonical refinement byte identity')
        run(flags+[str(HERE/'check_refinement.py'),'--repo-root',str(root),'--report',str(report),'--output',str(review)],mode+'-checker')
        for key,path in [('REFINEMENT',report),('CHECK',review)]:
            sha=hashlib.sha256(path.read_bytes()).hexdigest();need(key not in pairs or pairs[key]==sha,'normal/optimized identity '+key);pairs[key]=sha
    run([sys.executable,'-I','-S','-B',str(HERE/'run_controls.py'),'--repo-root',str(root),'--report',str(out/'normal-REFINEMENT.json'),'--output-dir',str(out/'controls')],'controls')
    controls=json.loads((out/'controls/CONTROLS.json').read_bytes());need(controls['case_count']==14 and controls['execution_count']==28,'controls complete')
    result={'status':'PASS_PORTABLE_ABSOLUTE_REFINEMENT','canonical_sha256':pairs,'normal_optimized_byte_identity':True,'real_corruption_cases':14,'real_corruption_executions':28,'same_candidate_and_frozen_initial_error_ledger':True}
    (out/'VERIFICATION.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
