#!/usr/bin/env python3
"""Run independent response equations and canonicalize provenance input roles."""
import argparse,hashlib,json,subprocess,sys,tempfile
from pathlib import Path
sys.dont_write_bytecode=True
BASE='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def arguments(args,helper,output):
    roles={'own_initial':args.independent_dir/'INITIAL_NORMAL.json','own_regularity':args.independent_dir/'REGULARITY_FINAL_NORMAL.json','own_audit':args.independent_dir/'AUDIT_PORTABLE_FINAL_NORMAL.json','absolute':args.absolute_dir/'REFINEMENT.json','baseline':args.source/(BASE+'RESPONSE.json'),'baseline_initial':args.source/(BASE+'INITIAL.json'),'native_response':args.native_response}
    cmd=[sys.executable]+(['-O'] if sys.flags.optimize else [])+[str(helper),'--source',str(args.source),'--freeze',str(args.freeze),'--output',str(output)]
    for role,path in roles.items():cmd+=['--input',role+'='+str(path)]
    return cmd,roles
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--absolute-dir',required=True,type=Path);ap.add_argument('--native-response',required=True,type=Path);ap.add_argument('--independent-dir',type=Path,default=Path(__file__).resolve().parent.parent);ap.add_argument('--freeze',type=Path,default=Path(__file__).resolve().parent.parent/'CONSERVATIVE_RESPONSE_PORTABLE_FREEZE.json');ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    need(not args.output.exists(),'FRESH_PORTABLE_RESPONSE_OUTPUT');helper=Path(__file__).with_name('check_conservative_response.py');before={p:sha(p) for p in (Path(__file__),helper,args.freeze)}
    with tempfile.TemporaryDirectory(prefix='conservative-response-') as d:
        raw=Path(d)/'raw.json';cmd,roles=arguments(args,helper,raw);subprocess.run(cmd,check=True,capture_output=True,text=True);x=json.loads(raw.read_bytes())
    old=x.pop('authenticated_files_before_after');mapped={}
    exact_roles={str(helper):'independent/scripts/check_conservative_response.py',str(args.freeze):'independent/CONSERVATIVE_RESPONSE_PORTABLE_FREEZE.json'}
    for role,path in roles.items():
        name=('research/'+str(path).split('/research/',1)[1]) if role in ('baseline','baseline_initial') else {'own_initial':'independent/INITIAL_NORMAL.json','own_regularity':'independent/REGULARITY_FINAL_NORMAL.json','own_audit':'independent/AUDIT_PORTABLE_FINAL_NORMAL.json','absolute':'absolute/REFINEMENT.json','native_response':'cauchy/reference/RESPONSE.json'}[role]
        exact_roles[str(path)]=name
    for path,h in old.items():
        if path in exact_roles:role=exact_roles[path]
        else:
            p=Path(path);need(p.is_relative_to(args.source),'EXPECTED_SOURCE_PROVENANCE_ROLE');role=str(p.relative_to(args.source))
        need(role not in mapped,'UNIQUE_PORTABLE_PROVENANCE_ROLE');mapped[role]=h
    x['portable_authenticated_files_sha256_before_after']=mapped;x['portable_wrapper_sha256_before_after']=sha(Path(__file__));x['projection_scope']='Only provenance path keys are replaced by portable source/input roles. Exact rational equations, caps, terms and gates are unchanged.'
    need(all(sha(p)==h for p,h in before.items()),'PORTABLE_RESPONSE_SOURCE_CLOSURE_AFTER');out=json.dumps(x,sort_keys=True,indent=2)+'\n';need('/workspace/' not in out,'NO_PRIVATE_PATHS');args.output.write_text(out);print(x['status'])
if __name__=='__main__':main()
