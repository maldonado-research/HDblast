#!/usr/bin/env python3
"""Execute report/disk controls against the separate independent tree auditor."""
import argparse
import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def need(x,why):
    if not x:raise RuntimeError(why)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path);ap.add_argument('--producer-root',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    here=Path(__file__).resolve().parent.parent;work=here/'audit-controls';need(not work.exists() and not args.output.exists(),'FRESH_CONTROL_OUTPUT');work.mkdir()
    auditor=here/'scripts/audit_native_structure.py';before={str(p.relative_to(here)):sha(p) for p in (auditor,here/'source-pins.json',here/'producer-source-pins.json')}
    original=json.loads((args.producer_root/'REFINEMENT.json').read_text());cases=[]
    def mutate(name,fun,reason):
        doc=copy.deepcopy(original);fun(doc);p=work/(name+'.json');p.write_text(json.dumps(doc,sort_keys=True,indent=2)+'\n');cases.append((name,p,args.source_root,args.producer_root,reason))
    mutate('coverage_gap',lambda r:r['nodes'][0]['r'].__setitem__(1,'999/1000'),'PARENT_CHILD_DOMAIN')
    mutate('coverage_overlap',lambda r:r['nodes'][1]['r'].__setitem__(0,'999/1000'),'PARENT_CHILD_DOMAIN')
    mutate('dropped_component',lambda r:r['nodes'][0]['raw_component_caps'][0].pop(),'COMPLETE_3_BY_23')
    mutate('wrong_gram',lambda r:r['Gram'].__setitem__(0,'2'),'FULL23_LABELS_GRAM')
    mutate('understated_leaf_H',lambda r:r['nodes'][0]['absolute_H_caps'].__setitem__(0,'1'),'FULL23_H_AGGREGATION')
    mutate('understated_global_H',lambda r:r['candidate_initial_absolute_H_caps'].__setitem__('u_rr','0'),'ALL_LEAF_AGGREGATION')
    mutate('changed_error',lambda r:r['initial_error_H_caps'].__setitem__('e0_rr','0'),'FROZEN_ERROR_LEDGER')
    mutate('changed_velocity',lambda r:r['velocity']['absolute_H_caps'].__setitem__('v','0'),'FROZEN_VELOCITY')
    def parent(r):
        n=next(n for n in r['nodes'] if n['parent'] is not None);n['effective_component_caps'][0][12]='0'
    mutate('invalid_parent_min',parent,'PARENT_MIN_INTERSECTION')
    mutate('omitted_tail',lambda r:r['leaves'].remove('s1'),'COMPLETE_LEAF_ROSTER')
    mutate('scope_changed',lambda r:r.__setitem__('same_candidate_unknown_true_error_unchanged',False),'SCOPE_FLAGS')
    mutate('old_absolute_changed',lambda r:r['previous_candidate_initial_absolute_H_caps'].__setitem__('u_rr','1'),'FROZEN_OLD_ABSOLUTE')
    # Disk failures must happen while authenticating, before any formulas or
    # report certificate gates can run. Copy only selected paths, never source.
    pins=json.loads((here/'source-pins.json').read_text());source_copy=work/'source-copy'
    for rel in pins:
        to=source_copy/rel;to.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.source_root/rel,to)
    rel='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/checks/producer.py';bad=source_copy/rel;bad.write_text(bad.read_text()+'\n# disk control\n')
    cases.append(('source_formula_tamper',args.producer_root/'REFINEMENT.json',source_copy,args.producer_root,'SOURCE_HASH:'+rel))
    producer_copy=work/'producer-copy';producer_copy.mkdir()
    for rel in json.loads((here/'producer-source-pins.json').read_text()):shutil.copyfile(args.producer_root/rel,producer_copy/rel)
    (producer_copy/'refine.py').write_text((producer_copy/'refine.py').read_text()+'\n# disk control\n')
    cases.append(('producer_formula_tamper',args.producer_root/'REFINEMENT.json',args.source_root,producer_copy,'SOURCE_HASH:refine.py'))
    # A selected source file reached through a symbolic root is refused too.
    linked=work/'source-symlink';linked.symlink_to(args.source_root,target_is_directory=True)
    cases.append(('source_root_symlink',args.producer_root/'REFINEMENT.json',linked,args.producer_root,'SOURCE_ROOT_TYPE'))
    rows=[]
    for name,report,source,producer,reason in cases:
        for optimized in (False,True):
            out=work/(name+('-optimized' if optimized else '-normal')+'.result.json')
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(auditor),'--source-root',str(source),'--producer-root',str(producer),'--report',str(report),'--source-pins',str(here/'source-pins.json'),'--source-pins-sha256',sha(here/'source-pins.json'),'--producer-pins',str(here/'producer-source-pins.json'),'--producer-pins-sha256',sha(here/'producer-source-pins.json'),'--output',str(out)]
            p=subprocess.run(cmd,capture_output=True,text=True)
            need(p.returncode!=0 and reason in p.stderr and not out.exists(),'CONTROL_FAILED:'+name)
            rows.append({'control':name,'optimized':optimized,'returncode':p.returncode,'expected_rejection':reason,'expected_rejection_observed':True,'no_success_output':not out.exists()})
    need(before=={rel:sha(here/rel) for rel in before},'AUDITOR_AND_MAPS_UNCHANGED')
    report={'status':'PASS_EFFECTIVE_INDEPENDENT_NATIVE_AUDIT_CONTROLS','controls':len(cases),'actual_failing_subprocesses':len(rows),'every_case_rejected_in_normal_and_optimized':True,'rows':rows,'auditor_and_maps_sha256':before,'limitations':'Structural auditor does not claim independent numerical replay of native raw tables. Separate native numerical controls validate formula-table understatement.'}
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n');print(report['status'])


if __name__=='__main__':main()
