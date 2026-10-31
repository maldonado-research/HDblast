#!/usr/bin/env python3
"""Real copied-source corruptions must fail before any successful report."""
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
G='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/'
S='research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/'

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def mutate(root,name):
    if name=='drop_TT_bulk_term':
        p=root/G/'projection/GEOMETRIC_COEFFICIENT_MAPS.json';d=json.loads(p.read_bytes());d['maps']['metric_forcing']['ab'][5].pop();p.write_text(json.dumps(d))
    elif name=='reverse_shell_units':
        p=root/G/'projection/inputs/linear_operator.py';v=p.read_text();old="'ab':tuple(2*ell*rho*x for x in sources['QKab'])";need(old in v,'mutation target');p.write_text(v.replace(old,"'ab':tuple(-2*ell*rho*x for x in sources['QKab'])",1))
    elif name=='omit_fixed_quadratic_term':
        p=root/S/'theory/screen_coefficient.py';v=p.read_text();old="-m['E']*(m['D']+m['L'])/(2*R)";need(old in v,'mutation target');p.write_text(v.replace(old,'-0',1))
    elif name=='omit_reflected_actual_errors':
        p=root/G/'numeric/actual_first_jets.py';v=p.read_text();old="if mode=='actual':";need(old in v,'mutation target');p.write_text(v.replace(old,"if False:",1))
    else:raise ValueError('unknown fault')
    return p

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    need(not a.output.exists(),'fresh output');rows=[]
    with tempfile.TemporaryDirectory(prefix='source-controls-',dir=HERE) as tmp:
        base=Path(tmp);copy=base/'repo';(copy/'research').mkdir(parents=True)
        for cap in CAPS:shutil.copytree(a.repo_root/'research'/cap,copy/'research'/cap)
        for name in ['drop_TT_bulk_term','reverse_shell_units','omit_fixed_quadratic_term','omit_reflected_actual_errors']:
            if name=='omit_fixed_quadratic_term':rel=S+'theory/screen_coefficient.py'
            elif name=='omit_reflected_actual_errors':rel=G+'numeric/actual_first_jets.py'
            elif name=='reverse_shell_units':rel=G+'projection/inputs/linear_operator.py'
            else:rel=G+'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
            original=(copy/rel).read_bytes();oldhash=sha(copy/rel);changed=mutate(copy,name)
            need(sha(changed)!=oldhash,'real bytes changed')
            for optimized in [False,True]:
                out=base/(name+('-optimized' if optimized else '-normal')+'.json')
                cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'producer.py'),'--repo-root',str(copy),'--output',str(out)]
                p=subprocess.run(cmd,capture_output=True,text=True)
                need(p.returncode!=0 and not out.exists(),'EFFECTIVE_SOURCE_REJECTION:'+name)
                need('MEMBER_BYTES:' in p.stderr,'expected authentication failure')
                rows.append({'fault':name,'optimized':optimized,'relative_source':rel,'original_sha256':oldhash,'mutated_sha256':sha(changed),'exit_code':p.returncode,'no_report':not out.exists(),'failure_reason':p.stderr.splitlines()[-1]})
            (copy/rel).write_bytes(original)
    result={'status':'PASS_REAL_COPIED_SOURCE_CORRUPTIONS_NORMAL_AND_OPTIMIZED','producer_sha256':sha(HERE/'producer.py'),'cases':rows}
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(result['status'])

if __name__=='__main__':main()
