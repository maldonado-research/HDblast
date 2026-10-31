#!/usr/bin/env python3
"""Execute this first experiment with relocated code and four copied capsules."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
CAPS=['HDBLAST_GEOMETRIC_BULK_SOURCES_20261016','HDBLAST_FULL_ANGULAR_RESPONSE_20261026','HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026','HDBLAST_EXACT_PACKET_INITIAL_20261031']
FILES=['producer.py','FREEZE.json','STRUCTURE_NORMAL.json','check_tt_block.py','SOURCE_PINS.json','REFLECTED_TT_PROOF.md']
def need(ok,message):
    if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    need(not a.output.exists(),'fresh receipt');expected=(HERE/'POSTCONTACT_NORMAL.json').read_bytes()
    with tempfile.TemporaryDirectory(prefix='relocated-experiment-',dir=HERE) as td:
        root=Path(td);repo=root/'inputs';code=root/'unrelated-code-location';code.mkdir();(repo/'research').mkdir(parents=True)
        for cap in CAPS:shutil.copytree(a.repo_root/'research'/cap,repo/'research'/cap)
        for name in FILES:shutil.copyfile(HERE/name,code/name)
        out=root/'report.json';command=[sys.executable,'-O',str(code/'producer.py'),'--repo-root',str(repo),'--output',str(out)]
        p=subprocess.run(command,capture_output=True,text=True)
        need(p.returncode==0 and out.exists(),'relocated successful execution:'+p.stderr)
        need(out.read_bytes()==expected,'relocated exact byte identity')
        report={'status':'PASS_REAL_RELOCATED_CODE_AND_FOUR_CAPSULE_EXECUTION','optimized':True,'report_sha256':sha(out),'canonical_report_sha256':sha(HERE/'POSTCONTACT_NORMAL.json'),'report_byte_identity':True,'copied_capsules':CAPS,'copied_local_files':{name:sha(code/name) for name in FILES}}
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(report['status'])
if __name__=='__main__':main()
