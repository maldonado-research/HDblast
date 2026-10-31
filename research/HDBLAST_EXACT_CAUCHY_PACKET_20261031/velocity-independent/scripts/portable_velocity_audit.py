#!/usr/bin/env python3
"""Run the sealed arithmetic audit and write canonical portable input roles."""
import argparse,hashlib,json,subprocess,sys,tempfile
from pathlib import Path
sys.dont_write_bytecode=True
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser()
    for key in ('source','own','native','producer-root','output'):ap.add_argument('--'+key,required=True,type=Path)
    args=ap.parse_args();need(not args.output.exists(),'FRESH_PORTABLE_OUTPUT');helper=Path(__file__).with_name('check_velocity_reports.py');before={p:sha(p) for p in (Path(__file__),helper,args.own,args.native)}
    with tempfile.TemporaryDirectory(prefix='velocity-audit-') as d:
        raw=Path(d)/'raw.json';cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(helper),'--source',str(args.source),'--own',str(args.own),'--native',str(args.native),'--producer-root',str(args.producer_root),'--out',str(raw)]
        subprocess.run(cmd,check=True,capture_output=True,text=True);x=json.loads(raw.read_bytes())
    mapped=x.pop('input_and_self_sha256_before_after')
    need(mapped=={str(helper):sha(helper),str(args.own):sha(args.own),str(args.native):sha(args.native)},'RAW_AUDIT_EXACT_INPUT_MAP')
    x['portable_input_roles_sha256_before_after']={'independent_checker_source':sha(helper),'independent_initial_report':sha(args.own),'native_initial_report':sha(args.native)}
    x['portable_wrapper_sha256_before_after']=sha(Path(__file__))
    x['projection_scope']='Only provenance metadata path keys are replaced by portable input roles; all scientific/audit fields retain their exact values.'
    need(all(sha(p)==h for p,h in before.items()),'PORTABLE_AUDIT_SOURCE_INPUT_CLOSURE_AFTER');data=json.dumps(x,sort_keys=True,indent=2)+'\n';need('/workspace/' not in data,'NO_PRIVATE_PATHS');args.output.write_text(data);print(x['status'])
if __name__=='__main__':main()
