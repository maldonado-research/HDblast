#!/usr/bin/env python3
"""Assertion-free complete replay of absolute-cap refinement and all23 aggregation."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('refinement_auth',HERE/'auth.py');a=importlib.util.module_from_spec(s);s.loader.exec_module(a)
need=a.need

def norm(caps,weights,bits):
    square=sum((w*x*x for w,x in zip(weights,caps)),Q(0));den=1<<bits
    n=math.isqrt(square.numerator*den*den//square.denominator)
    if Q(n*n,den*den)<square:n+=1
    return Q(n,den)

def table(rows):
    need(isinstance(rows,list) and len(rows)==3 and all(isinstance(r,list) and len(r)==23 for r in rows),'all23 complete3-jet table')
    out=[[Q(x) for x in row] for row in rows]
    need(all(x>=0 for row in out for x in row),'nonnegative component caps')
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();need(not args.output.exists(),'fresh checker output required')
    f,cap,receipt,p=a.load_producer(args.repo_root,HERE/'FREEZE.json')
    r=a.json_read(args.report);old=a.json_read(cap/'reference/INITIAL.json')['campaigns'][0]
    need(r['schema']=='HDBLAST_SAME_CANDIDATE_ABSOLUTE_REFINEMENT_V1' and r['candidate']=='W1_exact_packet','same candidate schema')
    need(r['freeze_sha256']==a.digest(HERE/'FREEZE.json') and r['source_initial_sha256']==f['initial_report_sha256'],'freeze and source report pins')
    need(r['producer_sha256']==a.digest(HERE/'refine.py'),'refinement producer pin')
    need(r['authentication']==receipt and r['authentication_after_matches'] is True,'source authentication receipt')
    need(r['bits']==128 and r['precision192_used'] is False,'first precision campaign only')
    need(r['labels']==old['labels'] and len(r['labels'])==23 and len(set(r['labels']))==23,'exact all23 labels')
    need(r['Gram']==old['Gram'],'exact all23 Gram weights')
    weights=list(map(Q,r['Gram']));need(all(w>0 for w in weights),'positive Gram')
    need(r['seeds']==f['seed_partition']==old['seeds'],'frozen exact seeds')
    need(r['initial_error_H_caps']==old['initial_error_H_caps']==f['frozen_complete_initial_error_H_caps'],'frozen complete initial errors')
    need(r['velocity']==old['velocity']==f['unchanged_velocity'],'unchanged complete velocity')
    need(r['previous_candidate_initial_absolute_H_caps']==old['candidate_initial_absolute_H_caps'],'previous absolute caps identity')
    need(r['complete_initial_errors_frozen'] is True and r['all23_completed_components_in_each_node'] is True and r['same_candidate_unknown_true_error_unchanged'] is True and r['response_evolution_performed'] is False,'same-function scope flags')
    need(r['definition_and_scope']==f['scope'],'scope identity')
    need(r['exact_support_collar']==['0','1'] and r['analytic_tail_retained']==['1','65/64'],'collar and tail retained')
    need(r['endpoint_continuation']=='Same completed-component quadratic extension matching exact value/r/rr at7/4; initial physical footprint ends at7/4; no exterior samples needed for this initial cap.','same completed exterior continuation')
    nodes={n['id']:n for n in r['nodes']};need(len(nodes)==len(r['nodes']),'unique node IDs')
    need(len(nodes)==len(old['seeds'])+2*r['split_count'],'full binary refinement node count')
    parsed={}
    for n in r['nodes']:
        raw=table(n['raw_component_caps']);effective=table(n['effective_component_caps']);bounds=list(map(Q,n['r']))
        need(len(bounds)==2 and bounds[0]<bounds[1],'positive leaf interval')
        need(isinstance(n['depth'],int) and 0<=n['depth']<=f['max_depth_relative_to_seed'],'depth budget')
        parsed[n['id']]=(raw,effective,bounds)
    visited=[]
    def walk(identifier,seed,depth,interval,parent):
        need(identifier in nodes and identifier not in visited,'reachable unique tree node')
        visited.append(identifier);n=nodes[identifier];raw,effective,bounds=parsed[identifier]
        need(n['seed']==seed and n['depth']==depth and bounds==interval and n['parent']==parent,'exact parent-child subset identity')
        if parent is None:
            need(n['raw_source']=='authenticated_published_seed','root raw source')
            need(n['raw_component_caps']==old['leaves'][seed]['absolute_component_caps'],'authenticated root component caps')
            expected=raw
        else:
            need(n['raw_source']=='fresh_completed_ExactPacket.evaluate','child raw source')
            previous=parsed[parent][1]
            expected=[[min(x,y) for x,y in zip(row,oldrow)] for row,oldrow in zip(raw,previous)]
        need(effective==expected,'valid parent subset component minimum')
        need(list(map(Q,n['absolute_H_caps']))==[norm(row,weights,128) for row in effective],'all23 node Gram aggregation')
        if n['children']:
            need(n['children']==[identifier+'L',identifier+'R'],'binary child IDs')
            lo,hi=interval;mid=(lo+hi)/2
            walk(identifier+'L',seed,depth+1,[lo,mid],identifier)
            walk(identifier+'R',seed,depth+1,[mid,hi],identifier)
    for seed,interval in enumerate(f['seed_partition']):walk('s'+str(seed),seed,0,list(map(Q,interval)),None)
    need(set(visited)==set(nodes),'no disconnected nodes')
    leaves=sorted((n for n in nodes.values() if not n['children']),key=lambda n:Q(n['r'][0]))
    need(r['leaves']==[n['id'] for n in leaves],'complete ordered leaf roster')
    intervals=[list(map(Q,n['r'])) for n in leaves]
    need(intervals[0][0]==0 and intervals[-1][1]==Q(7,4) and all(x[1]==y[0] for x,y in zip(intervals,intervals[1:])),'exact complete leaf coverage')
    need(r['splits']==[n['id'] for n in r['nodes'] if n['children']],'split roster')
    need(r['leaf_count']==len(leaves)==len(old['seeds'])+r['split_count'] and len(leaves)<=f['max_total_leaves'],'leaf budget and counts')
    need(r['max_depth']==max(n['depth'] for n in leaves),'depth summary')
    caps={key:max(Q(n['absolute_H_caps'][k]) for n in leaves) for k,key in enumerate(('u','u_r','u_rr'))}
    caps.update({k:Q(old['candidate_initial_absolute_H_caps'][k]) for k in ('v','v_r')})
    need({k:Q(v) for k,v in r['candidate_initial_absolute_H_caps'].items()}==caps,'complete all-leaf absolute aggregation')
    oldcaps={k:Q(v) for k,v in old['candidate_initial_absolute_H_caps'].items()}
    accepted=all(caps[k]<=oldcaps[k] for k in ('u','u_r','u_rr')) and caps['u_rr']<oldcaps['u_rr']
    desired=caps['u_rr']<=oldcaps['u_rr']/2
    need(r['strict_rr_reduction']==accepted and r['desired_factor_two_rr_reduction']==desired,'exact reduction gates')
    need(accepted and desired and r['status']=='PASS_DESIRED_ABSOLUTE_RR_REDUCTION','successful frozen absolute gates')
    need(Q(r['absolute_rr_certificate_reduction_factor'])==oldcaps['u_rr']/caps['u_rr'],'exact rr factor')
    # This numerical replay evaluates every root and every child. In particular,
    # an attacker cannot lower a child raw table and recompute a matching norm.
    e=p.ExactPacket(args.repo_root,128)
    need(p.encode(e.denominators())==r['denominators'],'actual and nominal denominator separation')
    for n in r['nodes']:
        lo,hi=map(Q,n['r']);_,nominal,_=e.evaluate(lo,hi)
        raw=[[v.mag() for v in row] for row in nominal]
        need(raw==parsed[n['id']][0],'completed exact numerical raw cap '+n['id'])
    _,_,after=a.authenticate(args.repo_root,HERE/'FREEZE.json');need(after==receipt,'source authentication after replay')
    result={'status':'PASS_COMPLETE_ABSOLUTE_REFINEMENT_REPLAY','report_sha256':a.digest(args.report),
      'checker_sha256':a.digest(Path(__file__)),'freeze_sha256':a.digest(HERE/'FREEZE.json'),
      'source_initial_sha256':f['initial_report_sha256'],'complete_numerical_nodes_replayed':len(nodes),
      'complete_all23_H_aggregations_checked':3*len(nodes),'leaf_count':len(leaves),'split_count':r['split_count'],
      'exact_covered_interval':['0','7/4'],'tail_exact_interval':['1','65/64'],'max_depth':r['max_depth'],
      'candidate_initial_absolute_H_caps':r['candidate_initial_absolute_H_caps'],'initial_error_H_caps':old['initial_error_H_caps'],
      'same_candidate_unknown_true_error_unchanged':True,'authentication_before_and_after':True,'precision_bits':128}
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps({'status':result['status'],'nodes':len(nodes),'leaves':len(leaves)},sort_keys=True))
if __name__=='__main__':main()
