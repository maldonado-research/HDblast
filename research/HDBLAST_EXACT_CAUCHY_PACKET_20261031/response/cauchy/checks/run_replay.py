#!/usr/bin/env python3
"""Rebuild the distinct Cauchy response and rerun all effective controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(ok,message):
    if not ok:raise RuntimeError(message)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    for role in ('same-response','absolute','velocity','independent-velocity'):
        ap.add_argument('--'+role+'-dir',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args();base=Path(__file__).resolve().parent.parent;out=args.output_dir
    need(not out.exists(),'fresh replay directory');out.mkdir(parents=True)
    reference=base/'reference/RESPONSE.json';r=json.loads(reference.read_text())
    adapter=base/'checks/response_cauchy.py';freeze=base/'checks/FREEZE.json'
    need(sha(adapter)==r['selected_input_sha256_before']['response_cauchy.py'],'canonical arithmetic recipe pin')
    need(sha(freeze)==r['selected_input_sha256_before']['FREEZE.json'],'canonical input freeze pin')
    common=['--repo-root',str(args.repo_root)]
    for role in ('same_response','absolute','velocity','independent_velocity'):
        common += ['--'+role.replace('_','-')+'-dir',str(getattr(args,role+'_dir'))]
    flags=['--freeze-sha256',sha(freeze)]
    paths=[]
    for mode,opts in [('normal',[]),('optimized',['-O'])]:
        p=out/('RESPONSE_'+mode+'.json');paths.append(p)
        subprocess.run([sys.executable]+opts+[str(adapter)]+common+flags+['--verify-report',str(reference),'--output',str(p)],check=True)
    need(paths[0].read_bytes()==paths[1].read_bytes()==reference.read_bytes(),'normal/optimized canonical exact byte identity')
    subprocess.run([sys.executable,str(base/'checks/run_controls.py')]+common+['--response',str(reference),'--output-dir',str(out/'controls')],check=True)
    controls=json.loads((out/'controls/CONTROLS.json').read_text());need(controls['count']==60,'complete Cauchy corruption count')
    receipt={'status':'PASS_PORTABLE_SEPARATE_EXACT_CAUCHY_RESPONSE','response_sha256':sha(reference),
        'normal_optimized_canonical_identity':True,'rejected_corrupted_subprocesses':60,
        'full_response_gates_against_PR109':r['full_response_gates_against_PR109'],
        'runner_sha256':sha(__file__),'conditional_physical_side_estimates_only':True,
        'distinct_candidate_and_fresh_velocity_charges':True}
    (out/'REPLAY.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':main()
