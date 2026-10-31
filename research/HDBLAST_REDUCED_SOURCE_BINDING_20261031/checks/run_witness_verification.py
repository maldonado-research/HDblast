#!/usr/bin/env python3
"""Fresh normal/optimized witness baselines and deliberate premise misuse."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
PIN='0bf27b73a4bee19b6caec83f7841335ef0680443027ba067a87cdb7a7e28aecf'
CONTROLS=('exchange_QC_under_wave','point_E1_implies_prolongation')
def need(ok,why):
    if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def auth(root):
    source=HERE/'coherent_witness.py'
    need(source.is_file() and not source.is_symlink() and sha(source)==PIN,'WITNESS_SOURCE_PIN')
    for p in (source,Path(__file__).resolve()):
        tree=ast.parse(p.read_text());need(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'ASSERT_FORBIDDEN')
        allowed={'argparse','ast','fractions','hashlib','importlib','itertools','json','math','pathlib','subprocess','sys','time'}
        for n in ast.walk(tree):
            if isinstance(n,ast.Import):need(all(x.name.split('.')[0] in allowed for x in n.names),'IMPORT_CLOSURE')
            if isinstance(n,ast.ImportFrom):need(n.module and n.module.split('.')[0] in allowed,'IMPORT_CLOSURE')
    spec=importlib.util.spec_from_file_location('coherent_witness_auth',source);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m.authenticate(root)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);a=ap.parse_args()
    root=a.repo_root.resolve();out=a.output_dir.resolve();need(not out.exists(),'REFUSE_OVERWRITE');out.mkdir(parents=True)
    before=auth(root);own_sha=sha(Path(__file__).resolve());baselines=[];controls=[]
    for mode in ('NORMAL','OPTIMIZED'):
        for mutation in (None,)+CONTROLS:
            tag=mode+('_'+mutation if mutation else '');output=out/(tag+'.json')
            need(not output.exists(),'OUTPUT_PREEXISTS:'+tag);need(auth(root)==before,'INPUT_CHANGED_BEFORE:'+tag)
            command=[sys.executable,'-I','-S','-B']+(['-O'] if mode=='OPTIMIZED' else [])+[str(HERE/'coherent_witness.py'),'--repo-root',str(root),'--output',str(output)]+(['--mutation',mutation] if mutation else [])
            started=time.monotonic();result=subprocess.run(command,capture_output=True,text=True,timeout=180)
            row={'mode':mode,'command':command,'returncode':result.returncode,'elapsed_seconds':round(time.monotonic()-started,6)}
            need(auth(root)==before,'INPUT_CHANGED_AFTER:'+tag)
            if mutation:
                marker='ValueError: NONZERO:'
                need(result.returncode!=0 and marker in result.stderr and not output.exists(),'CONTROL_NOT_REJECTED:'+tag+':'+result.stderr)
                row.update(mutation=mutation,exact_nonzero_witness=result.stderr.split(marker,1)[1].splitlines()[0],successful_report_absent=True);controls.append(row)
            else:
                need(result.returncode==0 and output.exists(),'BASELINE_FAILURE:'+tag+':'+result.stderr)
                value=json.loads(output.read_text());need(value['status']=='PASS_EXACT_COHERENT_POINTWISE_COUNTEREXAMPLES','BASELINE_STATUS:'+tag)
                row['sha256']=sha(output);baselines.append(row)
    need(baselines[0]['sha256']==baselines[1]['sha256'],'NORMAL_OPTIMIZED_MISMATCH')
    for mutation in CONTROLS:
        rows=[row for row in controls if row['mutation']==mutation];need(rows[0]['exact_nonzero_witness']==rows[1]['exact_nonzero_witness'],'CONTROL_WITNESS_MISMATCH:'+mutation)
    need(auth(root)==before,'AUTHENTICATION_CHANGED');need(sha(Path(__file__).resolve())==own_sha,'DRIVER_CHANGED')
    report={'status':'PASS_COHERENT_WITNESS_VERIFICATION','authentication':before,'driver_sha256':own_sha,'baseline_reports_byte_identical':True,'fresh_positive_subprocesses':baselines,'fresh_negative_subprocesses':controls}
    (out/'VERIFICATION.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(report['status'])
if __name__=='__main__':main()
