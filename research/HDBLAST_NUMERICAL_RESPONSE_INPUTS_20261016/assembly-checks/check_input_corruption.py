#!/usr/bin/env python3
"""Fresh changed-master-byte checks in isolated normal and optimized children."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def need(ok,reason):
    if not ok:raise ValueError(reason)


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();root=Path(__file__).resolve().parents[1];numeric=root/'numeric'
    need(not args.output.exists(),'OUTPUT_MUST_BE_NEW')
    need(root not in args.output.resolve().parents,'OUTPUT_MUST_BE_EXTERNAL')
    files=[numeric/'check_approximation_errors.py',numeric/'INPUT_PINS.json',*sorted((numeric/'inputs').iterdir())]
    before={p.relative_to(numeric).as_posix():digest(p) for p in files}
    rows=[]
    with tempfile.TemporaryDirectory(prefix='hdblast-isolated-input-corruption-') as td:
        fixture=Path(td)/'fixture';fixture.mkdir()
        for name in ('check_approximation_errors.py','INPUT_PINS.json'):shutil.copy2(numeric/name,fixture/name)
        shutil.copytree(numeric/'inputs',fixture/'inputs')
        pins=json.loads((fixture/'INPUT_PINS.json').read_text());master=fixture/pins['master']['path']
        master.write_bytes(master.read_bytes()+b' ')
        for mode in ('normal','optimized'):
            target=Path(td)/(mode+'.json')
            command=[sys.executable,'-I','-S',*(['-O'] if mode=='optimized' else []),'-B',
                     str(fixture/'check_approximation_errors.py'),'--output',str(target)]
            run=subprocess.run(command,cwd=fixture,capture_output=True,timeout=120)
            reason='input authentication: master'
            need(run.returncode!=0 and reason.encode() in run.stderr and not target.exists(),
                 'CHANGED_MASTER_BYTE_NOT_REJECTED_FOR_AUTHENTICATION')
            rows.append({'mode':mode,'returncode':run.returncode,'specific_reason':reason,'report_absent':True,
                         'isolated_flags':['-I','-S','-B']+(['-O'] if mode=='optimized' else [])})
    need(before=={p.relative_to(numeric).as_posix():digest(p) for p in files},'ORIGINAL_NUMERIC_INPUTS_CHANGED')
    result={'status':'PASS_FRESH_ISOLATED_CHANGED_MASTER_BYTE_REJECTIONS','nested_rejection_processes':2,
            'nested_rejections':rows,'original_numeric_source_unchanged':True,
            'checker_sha256':digest(Path(__file__)),
            'scope':'One registered positive harness, two isolated changed-input rejection children. Not two extra direct roster rejections and not a new PDE proof.'}
    with args.output.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')


if __name__=='__main__':main()
