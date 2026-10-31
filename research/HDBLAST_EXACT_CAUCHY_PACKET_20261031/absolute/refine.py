#!/usr/bin/env python3
"""Tighten complete absolute caps for the unchanged exact-packet function."""
import argparse
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('refinement_auth',HERE/'auth.py')
auth=importlib.util.module_from_spec(spec);spec.loader.exec_module(auth)
need=auth.need

def make_node(identifier,seed,depth,interval,raw,e,parent=None):
    effective=raw if parent is None else [[min(v,w) for v,w in zip(row,oldrow)] for row,oldrow in zip(raw,parent['effective_component_caps'])]
    return {'id':identifier,'seed':seed,'depth':depth,'r':interval,'parent':None if parent is None else parent['id'],
            'children':[],'raw_component_caps':raw,'effective_component_caps':effective,
            'absolute_H_caps':[e.norm(row) for row in effective],
            'raw_source':'authenticated_published_seed' if parent is None else 'fresh_completed_ExactPacket.evaluate'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();need(not args.output.exists(),'fresh output required')
    f,cap,receipt,p=auth.load_producer(args.repo_root,HERE/'FREEZE.json')
    initial=auth.json_read(cap/'reference/INITIAL.json');old=initial['campaigns'][0]
    need(len(initial['campaigns'])==1 and old['bits']==128,'single published128 campaign')
    need(old['seeds']==f['seed_partition'],'same seed partition')
    need(old['initial_error_H_caps']==f['frozen_complete_initial_error_H_caps'],'frozen error ledger')
    e=p.ExactPacket(args.repo_root,128);den=e.denominators()
    nodes=[];leaves=[];splits=[]
    target=Q(f['old_absolute_H_caps']['u_rr'])/2
    for seed,row in enumerate(old['leaves']):
        need(row['seed']==seed and row['depth']==0 and row['r']==f['seed_partition'][seed],'published seed leaf identity')
        raw=[[Q(x) for x in rr] for rr in row['absolute_component_caps']]
        node=make_node('s'+str(seed),seed,0,list(map(Q,row['r'])),raw,e)
        need(node['absolute_H_caps']==list(map(Q,row['absolute_H_caps'])),'published all23 seed aggregation')
        stack=[node]
        while stack:
            n=stack.pop();nodes.append(n)
            if n['absolute_H_caps'][2]>target and n['depth']<f['max_depth_relative_to_seed'] and len(old['leaves'])+len(splits)<f['max_total_leaves']:
                lo,hi=n['r'];mid=(lo+hi)/2;children=[]
                for suffix,aa,bb in [('L',lo,mid),('R',mid,hi)]:
                    _,nominal,_=e.evaluate(aa,bb)
                    raw=[[v.mag() for v in rr] for rr in nominal]
                    children.append(make_node(n['id']+suffix,seed,n['depth']+1,[aa,bb],raw,e,n))
                n['children']=[c['id'] for c in children];splits.append(n['id']);stack.extend(reversed(children))
            else:
                leaves.append(n)
    leaves.sort(key=lambda n:n['r'][0]);need(p.cover([n['r'] for n in leaves],Q(0),Q(7,4)),'exact complete leaf coverage')
    absolute=[max(n['absolute_H_caps'][k] for n in leaves) for k in range(3)]
    accepted=all(absolute[k]<=Q(f['old_absolute_H_caps'][key]) for k,key in enumerate(('u','u_r','u_rr'))) and absolute[2]<Q(f['old_absolute_H_caps']['u_rr'])
    desired=absolute[2]<=target
    _,_,after=auth.authenticate(args.repo_root,HERE/'FREEZE.json');need(after==receipt,'source closure unchanged')
    caps=dict(zip(('u','u_r','u_rr'),absolute));caps.update({k:Q(old['candidate_initial_absolute_H_caps'][k]) for k in ('v','v_r')})
    out={'schema':'HDBLAST_SAME_CANDIDATE_ABSOLUTE_REFINEMENT_V1','status':'PASS_DESIRED_ABSOLUTE_RR_REDUCTION' if accepted and desired else ('PASS_STRICT_ABSOLUTE_RR_REDUCTION' if accepted else 'INCONCLUSIVE'),
        'candidate':'W1_exact_packet','bits':128,'freeze_sha256':auth.digest(HERE/'FREEZE.json'),
        'producer_sha256':auth.digest(Path(__file__)),'authentication':receipt,'authentication_after_matches':True,
        'source_initial_sha256':f['initial_report_sha256'],'labels':old['labels'],'Gram':old['Gram'],'denominators':den,
        'seeds':old['seeds'],'nodes':nodes,'leaves':[n['id'] for n in leaves],'splits':splits,
        'leaf_count':len(leaves),'split_count':len(splits),'max_depth':max(n['depth'] for n in leaves),
        'candidate_initial_absolute_H_caps':caps,'previous_candidate_initial_absolute_H_caps':old['candidate_initial_absolute_H_caps'],
        'initial_error_H_caps':old['initial_error_H_caps'],'velocity':old['velocity'],
        'complete_initial_errors_frozen':True,'all23_completed_components_in_each_node':True,
        'strict_rr_reduction':accepted,'desired_factor_two_rr_reduction':desired,
        'absolute_rr_certificate_reduction_factor':Q(f['old_absolute_H_caps']['u_rr'])/absolute[2],
        'precision192_used':False,'exact_support_collar':[Q(0),Q(1)],'analytic_tail_retained':[Q(1),Q(65,64)],
        'definition_and_scope':f['scope'],
        'same_candidate_unknown_true_error_unchanged':True,'response_evolution_performed':False,
        'endpoint_continuation':'Same completed-component quadratic extension matching exact value/r/rr at7/4; initial physical footprint ends at7/4; no exterior samples needed for this initial cap.'}
    args.output.write_text(json.dumps(p.encode(out),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'leaves':len(leaves),'splits':len(splits),'caps':{k:float(v) for k,v in caps.items()},'rr_factor':float(out['absolute_rr_certificate_reduction_factor'])},sort_keys=True))
if __name__=='__main__':main()
