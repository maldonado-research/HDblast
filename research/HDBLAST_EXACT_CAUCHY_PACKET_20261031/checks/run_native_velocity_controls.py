#!/usr/bin/env python3
"""Portable execution of the frozen native velocity corruption roster."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

def need(ok,why):
    if not ok:raise RuntimeError(why)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--velocity-dir',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();need(not a.output_dir.exists(),'fresh native velocity controls output');a.output_dir.mkdir(parents=True)
    canonical=json.loads((a.velocity_dir/'CONTROLS_NORMAL.json').read_bytes());names=sorted(canonical['cases']);need(len(names)==14,'fourteen frozen semantic cases')
    receipts=[]
    def run(name,script,arguments,expected):
        for optimized in (False,True):
            mode='optimized' if optimized else 'normal';out=a.output_dir/(name+'-'+mode+'.json')
            command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(script)]+arguments+['--output',str(out)]
            p=subprocess.run(command,capture_output=True,timeout=180);stderr=p.stderr.decode(errors='replace')
            need(p.returncode!=0 and not out.exists() and expected in stderr,'ineffective native velocity rejection '+name+' '+mode)
            (a.output_dir/(name+'-'+mode+'.stderr')).write_bytes(p.stderr)
            receipts.append({'case':name,'mode':mode,'returncode':p.returncode,'success_artifact_absent':True,'expected_rejection':expected,'expected_rejection_observed':True})
    for name in names:run(name,a.velocity_dir/'check_controls.py',['--repo-root',str(a.repo_root),'--inject-error',name],'CONTROL_REJECTED:'+name)
    corrupt=a.output_dir/'bad_auth';corrupt.mkdir();shutil.copyfile(a.velocity_dir/'producer.py',corrupt/'producer.py')
    frozen=json.loads((a.velocity_dir/'FREEZE.json').read_bytes());rel='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/checks/producer.py';frozen['input_sha256'][rel]='0'*64
    (corrupt/'FREEZE.json').write_text(json.dumps(frozen,sort_keys=True,indent=2)+'\n')
    run('source-hash-corruption-before-import',corrupt/'producer.py',['--repo-root',str(a.repo_root),'--local-only'],'AUTH_SOURCE:'+rel)
    need(len(receipts)==30,'thirty actual rejected processes')
    result={'status':'PASS_PORTABLE_NATIVE_VELOCITY_CONTROLS','rejected_process_count':len(receipts),'semantic_cases':14,'source_authentication_cases':1,'receipts':receipts}
    (a.output_dir/'CONTROLS.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'count':len(receipts)}))
if __name__=='__main__':main()
