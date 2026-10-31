#!/usr/bin/env python3
"""Stdlib-only independent semantic audit of the own absolute interval table."""
import argparse,hashlib,json,math
from fractions import Fraction as Q
from pathlib import Path

PACK='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/'
BASELINE_SHA='692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3'
OWN_REPORT_SHA='b7b0511e899b71b022d56e02dcbef0d116cf3fede9fc4e1c7bf7f73c77514275'
LABELS=[f'{f}:{l}' for f in ('K','L','N','J','w') for l in (0,2,4)]+[f'{f}:{l}' for f in ('v','u') for l in (2,4)]+['H:2','H:4','T:+','T:-']
GRAM=list(map(Q,['1','1/12','19/720','2','1/6','19/360','1','1/12','19/720','3','1/4','19/240','1','1/12','19/720','4/3','19/15','4/3','19/15','20/9','133/15','1/6','1/6']))
SEEDS=[(Q(0),Q(1)),(Q(1),Q(65,64))]+[(Q(k,64),Q(k+1,64)) for k in range(65,112)]

def need(ok,why):
 if not ok:raise RuntimeError(why)
def sha(p):
 p=Path(p);need(p.is_file() and not p.is_symlink() and all(not q.is_symlink() for q in p.parents),'REGULAR_FILE:'+p.name)
 return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,Q):return str(x)
 if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [enc(v) for v in x]
 return x
def closure(root,pins):
 need(root.is_dir() and not root.is_symlink(),'SOURCE_ROOT')
 for rel,pin in pins.items():
  name=Path(rel);need(not name.is_absolute() and '..' not in name.parts,'SOURCE_PATH')
  need(sha(root/name)==pin,'SOURCE_HASH:'+rel)
 return pins
def matrix(raw):
 need(isinstance(raw,list) and len(raw)==3 and all(isinstance(row,list) and len(row)==23 for row in raw),'FULL_3_BY_23')
 out=[[Q(v) for v in row] for row in raw];need(all(v>=0 for row in out for v in row),'NONNEGATIVE_COMPONENT_CAPS')
 need(all(row[i]==0 for row in out for i in [*range(6),*range(15,23)]),'STRUCTURAL_ZERO_COMPONENTS')
 return out
def minimum(a,b):return [[min(v,w) for v,w in zip(x,y)] for x,y in zip(a,b)]
def norm(row,bits):
 if all(v==0 for v in row):return Q(0)
 square=sum((g*v*v for g,v in zip(GRAM,row)),Q(0));grid=1<<bits
 return Q(math.isqrt((square*grid*grid).__ceil__())+1,grid)
def interval(raw):
 need(len(raw)==2,'INTERVAL_PAIR');lo,hi=map(Q,raw);need(lo<=hi,'INTERVAL_ORDER');return lo,hi

def tail_certificate(cert):
 def plus(a,b):return [(a[j] if j<len(a) else 0)+(b[j] if j<len(b) else 0) for j in range(max(len(a),len(b)))]
 def mul(a,b):
  out=[Q(0)]*(len(a)+len(b)-1)
  for j,x in enumerate(a):
   for k,y in enumerate(b):out[j+k]+=x*y
  return out
 ps=[[Q(1)]]
 for n in range(5):
  p=ps[-1];dp=[k*p[k] for k in range(1,len(p))]
  q=mul([0,0,1] if n%2==0 else [0,-4,1],plus(p,[-x for x in dp]))
  if n%2:q=plus(q,[-2*x for x in p])
  ps.append(q)
 need([[Q(v) for v in p] for p in cert['polynomials_ascending']]==ps,'TAIL_EXACT_POLYNOMIAL_RECURRENCE')
 s0=Q(cert['s_min']);need(s0==Q(4096,63),'TAIL_EXACT_DOMAIN')
 checks=[]
 for n,p in enumerate(ps):
  for sign in (-1,1):
   pol=plus([Q(0)]*(2*n)+[Q(1)],[sign*x for x in p])
   shifted=[sum((math.comb(k,j)*v*s0**(k-j) for k,v in enumerate(pol) if k>=j),Q(0)) for j in range(len(pol))]
   need(all(v>=0 for v in shifted) and s0>=2*n,'TAIL_DOMINANCE_AND_MONOTONICITY')
   checks.append(dict(order=n,sign=sign,shifted_coefficients=shifted))
 need(cert['shifted_nonnegative_coefficients']==enc(checks),'TAIL_SHIFTED_COEFFICIENTS')
 need('Every derivative tends to zero' in cert['flat_endpoint'],'TAIL_FLAT_ENDPOINT')
 return dict(exact_recurrence_orders=list(range(6)),polynomial_dominance_checks=len(checks),s_min=s0,flat_endpoint=True)

def audit_campaign(c,base,bg):
 bits=c['bits'];need(bits==128,'FROZEN_ACCEPTED_PRECISION')
 need(c['labels']==LABELS==base['labels'] and list(map(Q,c['Gram']))==GRAM==list(map(Q,base['Gram'])),'SAME_LABEL_ORDER_AND_POSITIVE_GRAM')
 seeds=[tuple(map(Q,v)) for v in c['seeds']];need(seeds==SEEDS and len(seeds)==49,'EXACT_49_SEEDS')
 leaves,splits=c['leaves'],c['splits']
 need(c['leaf_count']==len(leaves)==135 and c['split_count']==len(splits)==86 and len(leaves)==len(seeds)+len(splits),'EXACT_FINITE_COUNTS')
 need(len(leaves)<=4096,'FROZEN_LEAF_BUDGET')
 limits=[Q(base['candidate_initial_absolute_H_caps'][s]) for s in ('u','u_r','u_rr')];limits[2]/=2
 need(c['baseline_absolute_H_caps']=={s:base['candidate_initial_absolute_H_caps'][s] for s in ('u','u_r','u_rr')},'EXACT_BASELINE_CAPS')
 children={i:[] for i in range(len(splits))};roots=[];accepted=[];split_caps=[]
 def geometry(node,parent,kind,index):
  lo,hi=map(Q,node['r']);depth=node['depth'];seed=node['seed']
  need(type(depth)==int and 0<=depth<=16 and type(seed)==int and 0<=seed<len(seeds),'BOUNDED_DYADIC_DEPTH_AND_SEED')
  a,b=seeds[seed];need(a<=lo<hi<=b and hi-lo==(b-a)/(1<<depth),'DYADIC_DOMAIN_INCLUSION')
  need(((lo-a)/(hi-lo)).denominator==1,'DYADIC_ENDPOINT_GRID')
  if parent is None:
   need(depth==0 and (lo,hi)==seeds[seed],'SEED_ROOT_DOMAIN');roots.append(seed)
  else:
   need(type(parent)==int and 0<=parent<len(splits) and (kind!='split' or parent<index),'PARENT_SPLIT_REFERENCE')
   p=splits[parent];need(seed==p['seed'] and depth==p['depth']+1,'PARENT_DEPTH_AND_SEED')
   children[parent].append((lo,hi))
  return lo,hi
 for i,node in enumerate(splits):
  parent=node['parent_split'];lo,hi=geometry(node,parent,'split',i)
  need(Q(node['midpoint'])==(lo+hi)/2,'EXACT_SPLIT_MIDPOINT')
  raw=matrix(node['raw_component_caps']);cap=matrix(node['inherited_component_caps'])
  expected=raw if parent is None else minimum(raw,split_caps[parent])
  need(cap==expected,'EXACT_PARENT_COMPONENT_INTERSECTION');split_caps.append(cap)
 cursor=Q(0)
 for i,node in enumerate(leaves):
  parent=node['parent_split'];lo,hi=geometry(node,parent,'leaf',i)
  need(lo==cursor,'COMPLETE_LEAF_COVER_NO_GAP_OR_OVERLAP');cursor=hi
  raw=matrix(node['raw_component_caps']);cap=matrix(node['component_caps'])
  expected=raw if parent is None else minimum(raw,split_caps[parent])
  need(cap==expected,'EXACT_PARENT_COMPONENT_INTERSECTION')
  hh=[norm(row,bits) for row in cap]
  need(list(map(Q,node['absolute_H_caps']))==hh,'EXACT_LEAF_POSITIVE_GRAM_SQRT')
  gate=all(v<=threshold for v,threshold in zip(hh,limits))
  need(node['threshold_satisfied'] is gate and gate,'EVERY_LEAF_ACCEPTED');accepted.append(hh)
 need(cursor==Q(7,4),'COMPLETE_CLOSED_PHYSICAL_FOOTPRINT')
 need(sorted(roots)==list(range(49)),'EVERY_SEED_ROOT_EXACTLY_ONCE')
 for i,node in enumerate(splits):
  lo,hi=map(Q,node['r']);mid=(lo+hi)/2
  need(sorted(children[i])==[(lo,mid),(mid,hi)],'EXACT_TWO_CHILD_PARTITION')
 caps=[max(row[j] for row in accepted) for j in range(3)]
 need(c['complete_absolute_H_caps']==enc(dict(zip(('u','u_r','u_rr'),caps))),'EXACT_MAX_LEAF_GLOBAL_CAPS')
 need(all(v<=Q(base['candidate_initial_absolute_H_caps'][s]) for s,v in zip(('u','u_r','u_rr'),caps)) and caps[2]<=limits[2],'NO_WORSE_AND_HALF_CURVATURE')
 need(c['strict_no_worse'] is True and c['rr_factor_two_gate'] is True and c['all_leaf_thresholds_satisfied'] is True,'EXACT_GATE_FLAGS')
 need(c['exact_zero_region']==['0','1'] and c['flat_endpoint_all_jets_zero'] is True and all(Q(v)==0 for row in leaves[0]['component_caps'] for v in row),'EXACT_ZERO_FOOTPRINT')
 need(sum(node['seed']==1 for node in leaves)==1 and leaves[1]['r']==['1','65/64'],'COMPLETE_ANALYTIC_TAIL_LEAF')
 tail=tail_certificate(c['tail_certificate'])
 need('each completed component' in c['outer_extension'] and 'no C3 claim' in c['outer_extension'],'COMPLETED_QUADRATIC_EXTENSION_DEFINITION')
 need(len(c['endpoint_rows'])==23,'COMPLETE_ENDPOINT_COMPONENTS')
 for index,(label,g,row) in enumerate(zip(LABELS,GRAM,c['endpoint_rows'])):
  need(row['label']==label and Q(row['Gram'])==g,'ENDPOINT_LABEL_AND_GRAM')
  need(len(row['interior_closed_jet'])==3 and row['quadratic_exterior_minus_interior_trace']==[0,0,0],'EXACT_QUADRATIC_MATCHING_TRACES')
  for order,raw in enumerate(row['interior_closed_jet']):
   lo,hi=interval(raw);need(max(abs(lo),abs(hi))<=Q(leaves[-1]['component_caps'][order][index]),'ENDPOINT_IN_CLOSED_LAST_LEAF_CAP')
   if index<6 or index>=15:need(lo==hi==0,'ENDPOINT_STRUCTURAL_ZERO')
 # Exact rational Horner interval for the nominal denominator polynomials.
 z=(Q(0),Q(7,4*10**6))
 for name in ('rho','p'):
  out=(Q(0),Q(0))
  for v in reversed(list(map(Q,bg['representation_background_coefficients_ascending'][name]))):
   prods=[a*b for a in out for b in z];out=(min(prods)+v,max(prods)+v)
  lo,hi=interval(c['nominal_denominators'][name]);need(lo<=out[0]<=out[1]<=hi and (lo>0 or hi<0),'NOMINAL_DENOMINATOR_INCLUSION_AND_NONZERO')
 return dict(bits=bits,leaf_count=len(leaves),split_count=len(splits),seed_count=len(seeds),max_depth=max(v['depth'] for v in leaves),
   exact_parent_min_intersections=True,exact_Gram_aggregation=True,complete_closed_footprint=['0','7/4'],
   complete_absolute_H_caps=dict(zip(('u','u_r','u_rr'),caps)),endpoint_component_count=23,tail_proof=tail)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path);ap.add_argument('--support-root',required=True,type=Path)
 ap.add_argument('--report',required=True,type=Path);ap.add_argument('--report-sha256',required=True);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
 need(not args.output.exists(),'FRESH_OUTPUT')
 before={Path(__file__).resolve():sha(__file__)}
 need(sha(args.report)==args.report_sha256,'REPORT_HASH');before[args.report]=args.report_sha256
 report=json.loads(args.report.read_text());need(report['schema']=='HDBLAST_SAME_CANDIDATE_INDEPENDENT_ABSOLUTE_V1' and report['candidate']=='W1_exact_packet' and report['status']=='PASS_INDEPENDENT_ABSOLUTE_REFINEMENT','REPORT_ROLE')
 need(report['initial_error_velocity_target_source_operator_changed'] is False and report['source_preauth_before_helper_import'] is True and report['source_postauth_equal'] is True,'UNCHANGED_ROLE_AND_AUTHENTICATION')
 support_pins=dict(report['own_sha256']);support_pins.update({'FREEZE.json':report['freeze_sha256'],'source-pins.json':report['selected_source_map_sha256']})
 closure(args.support_root,support_pins)
 pins=json.loads((args.support_root/'source-pins.json').read_text());need(report['selected_source_count']==len(pins)==1456,'COMPLETE_SELECTED_SOURCE_COUNT');closure(args.source_root,pins)
 need(pins[PACK+'reference/INITIAL.json']==BASELINE_SHA,'PINNED_BASELINE_INITIAL')
 freeze=json.loads((args.support_root/'FREEZE.json').read_text());need(freeze['frozen_before_numerical_run'] is True and freeze['source_map_sha256']==report['selected_source_map_sha256'] and freeze['own_sha256']==report['own_sha256'],'FROZEN_PREDECLARED_SOURCE_BUDGET')
 need([(v['bits'],v['max_leaves'],v['max_dyadic_depth']) for v in freeze['campaigns']]==[(128,4096,16),(192,4096,16)],'EXACT_FROZEN_BOUNDED_CAMPAIGNS')
 baseline=json.loads((args.source_root/PACK/'reference/INITIAL.json').read_text());need(len(baseline['campaigns'])==1,'ONE_BASELINE_CAMPAIGN')
 bg=json.loads((args.source_root/'research/HDBLAST_INPUT_PRECISION_20261023/initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
 need(len(report['campaigns'])==1,'SUCCESSFUL_128_NEEDS_NO_192_FALLBACK')
 result=[audit_campaign(c,baseline['campaigns'][0],bg) for c in report['campaigns']]
 closure(args.source_root,pins);closure(args.support_root,support_pins);need(all(sha(p)==h for p,h in before.items()),'SOURCE_REPORT_CHECKER_POSTAUTH')
 output=dict(status='PASS_INDEPENDENT_ABSOLUTE_REPORT_SEMANTICS',candidate='W1_exact_packet',assertion_free=True,no_scientific_imports=True,
  checker_sha256=before[Path(__file__).resolve()],report_sha256=args.report_sha256,selected_source_map_sha256=report['selected_source_map_sha256'],selected_source_count=len(pins),
  supporting_source_sha256=support_pins,all_selected_sources_authenticated_before_after=True,campaigns=result,
  scope='Semantic source-authentication, complete tree, exact positive-Gram aggregation and endpoint/tail certificate audit. Full nominal expression interval replay was independently executed separately. Same function and true error retained; only upper allowances refined.')
 args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(enc(output),sort_keys=True,indent=2)+'\n')

if __name__=='__main__':main()
