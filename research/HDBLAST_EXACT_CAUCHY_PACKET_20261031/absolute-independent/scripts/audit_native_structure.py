#!/usr/bin/env python3
"""Independent report semantics and exact complete tree/Gram audit, no imports."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path


def need(x, why):
    if not x: raise RuntimeError(why)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def authenticate(root, pins):
    need(root.is_dir() and not root.is_symlink(),'SOURCE_ROOT_TYPE')
    for rel, pin in pins.items():
        name=Path(rel);need(not name.is_absolute() and '..' not in name.parts,'SAFE_SOURCE_PATH')
        p=root/name;need(p.is_file(),'SOURCE_FILE:'+rel)
        for pp in (p,*p.parents):
            if pp==root: break
            need(not pp.is_symlink(),'SOURCE_SYMLINK:'+rel)
        need(sha(p)==pin,'SOURCE_HASH:'+rel)


def norm(row, weights):
    square=sum((w*x*x for w,x in zip(weights,row)),Q(0));grid=1<<128
    k=math.isqrt((square*grid*grid).__floor__())
    if Q(k*k,grid*grid)<square: k+=1
    return Q(k,grid)


def table(x):
    need(isinstance(x,list) and len(x)==3 and all(isinstance(r,list) and len(r)==23 for r in x),'COMPLETE_3_BY_23')
    out=[list(map(Q,r)) for r in x];need(all(v>=0 for r in out for v in r),'NONNEGATIVE_CAPS');return out


def main():
    ap=argparse.ArgumentParser()
    for arg in ('source-root','producer-root','report','source-pins','producer-pins','output'):
        ap.add_argument('--'+arg,required=True,type=Path)
    ap.add_argument('--source-pins-sha256',required=True);ap.add_argument('--producer-pins-sha256',required=True)
    args=ap.parse_args();selfpin=sha(Path(__file__))
    need(not args.output.exists(),'FRESH_OUTPUT')
    need(sha(args.source_pins)==args.source_pins_sha256 and sha(args.producer_pins)==args.producer_pins_sha256,'PIN_MAP_AUTH')
    pins=json.loads(args.source_pins.read_text());ppins=json.loads(args.producer_pins.read_text())
    authenticate(args.source_root,pins);authenticate(args.producer_root,ppins)
    root=args.source_root/'research/HDBLAST_EXACT_PACKET_INITIAL_20261031'
    old=json.loads((root/'reference/INITIAL.json').read_text())['campaigns'][0]
    f=json.loads((args.producer_root/'FREEZE.json').read_text());r=json.loads(args.report.read_text())
    need(r['candidate']=='W1_exact_packet' and r['bits']==128 and r['status']=='PASS_DESIRED_ABSOLUTE_RR_REDUCTION','IDENTITY')
    need(r['labels']==old['labels'] and len(set(r['labels']))==23 and r['Gram']==old['Gram'],'FULL23_LABELS_GRAM')
    weights=list(map(Q,r['Gram']));need(min(weights)>0,'POSITIVE_GRAM')
    need(r['initial_error_H_caps']==old['initial_error_H_caps']==f['frozen_complete_initial_error_H_caps'],'FROZEN_ERROR_LEDGER')
    need(r['velocity']==old['velocity']==f['unchanged_velocity'],'FROZEN_VELOCITY')
    need(r['previous_candidate_initial_absolute_H_caps']==old['candidate_initial_absolute_H_caps'],'FROZEN_OLD_ABSOLUTE')
    need(r['seeds']==old['seeds']==f['seed_partition'],'SEED_IDENTITY')
    need(r['freeze_sha256']==sha(args.producer_root/'FREEZE.json') and r['producer_sha256']==sha(args.producer_root/'refine.py'),'PRODUCER_AND_FREEZE_IDENTITY')
    need(r['source_initial_sha256']==sha(root/'reference/INITIAL.json'),'INITIAL_IDENTITY')
    nodes={n['id']:n for n in r['nodes']};need(len(nodes)==len(r['nodes']),'UNIQUE_NODES')
    visited=set();leaf=[];splits=[]
    def walk(name, seed, depth, bounds, parent):
        need(name in nodes and name not in visited,'REACHABLE_UNIQUE_NODE');visited.add(name)
        n=nodes[name];raw=table(n['raw_component_caps']);eff=table(n['effective_component_caps'])
        need(n['seed']==seed and n['depth']==depth and 0<=depth<=16 and list(map(Q,n['r']))==bounds and n['parent']==parent,'PARENT_CHILD_DOMAIN')
        if parent is None:
            need(n['raw_component_caps']==old['leaves'][seed]['absolute_component_caps'],'AUTHENTICATED_ROOT_RAW')
            expected=raw
        else:
            par=table(nodes[parent]['effective_component_caps']);expected=[[min(a,b) for a,b in zip(x,y)] for x,y in zip(raw,par)]
        need(eff==expected,'PARENT_MIN_INTERSECTION')
        need(list(map(Q,n['absolute_H_caps']))==[norm(v,weights) for v in eff],'FULL23_H_AGGREGATION')
        if n['children']:
            need(n['children']==[name+'L',name+'R'],'BINARY_CHILDREN');splits.append(name)
            lo,hi=bounds;mid=(lo+hi)/2
            walk(name+'L',seed,depth+1,[lo,mid],name);walk(name+'R',seed,depth+1,[mid,hi],name)
        else: leaf.append(n)
    for seed,bounds in enumerate(r['seeds']): walk('s'+str(seed),seed,0,list(map(Q,bounds)),None)
    need(visited==set(nodes),'FULL_NODE_REACHABILITY')
    leaf.sort(key=lambda n:Q(n['r'][0]));need(r['leaves']==[n['id'] for n in leaf],'COMPLETE_LEAF_ROSTER')
    need(r['splits']==splits,'COMPLETE_SPLIT_ROSTER')
    bounds=[list(map(Q,n['r'])) for n in leaf]
    need(bounds[0][0]==0 and bounds[-1][1]==Q(7,4) and all(a[1]==b[0] for a,b in zip(bounds,bounds[1:])),'EXACT_CLOSED_COVERAGE')
    need(len(leaf)==r['leaf_count']==len(r['seeds'])+len(splits)<=4096 and r['split_count']==len(splits),'FROZEN_LEAF_BUDGET')
    need(r['max_depth']==max(n['depth'] for n in leaf),'DEPTH_SUMMARY')
    caps={k:max(Q(n['absolute_H_caps'][j]) for n in leaf) for j,k in enumerate(('u','u_r','u_rr'))}
    caps.update({k:Q(old['candidate_initial_absolute_H_caps'][k]) for k in ('v','v_r')})
    need(caps=={k:Q(v) for k,v in r['candidate_initial_absolute_H_caps'].items()},'ALL_LEAF_AGGREGATION')
    need(all(caps[k]<=Q(old['candidate_initial_absolute_H_caps'][k]) for k in ('u','u_r','u_rr')) and caps['u_rr']<=Q(old['candidate_initial_absolute_H_caps']['u_rr'])/2,'ABSOLUTE_CERTIFICATE_GATES')
    need(r['exact_support_collar']==['0','1'] and r['analytic_tail_retained']==['1','65/64'],'ZERO_TAIL_COVERAGE')
    need(r['complete_initial_errors_frozen'] is True and r['same_candidate_unknown_true_error_unchanged'] is True and r['response_evolution_performed'] is False,'SCOPE_FLAGS')
    need(Q(r['absolute_rr_certificate_reduction_factor'])==Q(old['candidate_initial_absolute_H_caps']['u_rr'])/caps['u_rr'],'EXACT_RR_FACTOR')
    authenticate(args.source_root,pins);authenticate(args.producer_root,ppins)
    need(sha(Path(__file__))==selfpin and sha(args.source_pins)==args.source_pins_sha256 and sha(args.producer_pins)==args.producer_pins_sha256,'POSTAUTH')
    out={'status':'PASS_INDEPENDENT_NATIVE_TREE_AND_LEDGER_AUDIT','report_sha256':sha(args.report),'auditor_sha256':selfpin,
         'source_pins_sha256':args.source_pins_sha256,'producer_pins_sha256':args.producer_pins_sha256,
         'pre_post_all_selected_source_authentication':True,'no_producer_or_core_imports':True,
         'nodes':len(nodes),'leaves':len(leaf),'splits':len(splits),'all23_Gram_aggregations':3*len(nodes),
         'exact_footprint':['0','7/4'],'tail':['1','65/64'],'frozen_error_velocity_verified':True,
         'native_raw_formula_replay_role':'Separate native assertion-free checker replays every raw node. This auditor independently checks tree, inheritance, all23 aggregation, complete footprint and frozen inputs.',
         'scope':'Same exact candidate; stronger absolute upper certificates only; actual error unchanged.'}
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
    print(out['status'])


if __name__=='__main__':main()
