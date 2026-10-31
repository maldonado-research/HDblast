#!/usr/bin/env python3
"""Authenticate and audit the independent evaluator's exact report ledger."""
import argparse, hashlib, importlib.util, json, math, sys
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
def need(ok,why):
    if not ok:raise RuntimeError(why)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--helper',type=Path,default=Path(__file__).with_name('check_independent_exact.py'));ap.add_argument('--source-pins',type=Path,default=Path(__file__).resolve().parent.parent/'source-pins.json');a=ap.parse_args()
    need(not a.output.exists(),'FRESH_OUTPUT_REQUIRED');own=digest(Path(__file__));raw=a.report.read_bytes();d=json.loads(raw);hp=digest(a.helper)
    need(hp==d['executed_helper_sha256'],'INDEPENDENT_HELPER_HASH')
    spec=importlib.util.spec_from_file_location('independent_scientific_oracle',a.helper);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.BITS=d['bits']
    ev=m.Evaluator(a.source,a.source_pins);need(d['selected_source_map_sha256']==m.SOURCE_MAP_PIN,'SOURCE_MAP_ROLE')
    need(d['label_order']==m.LABELS and list(map(Q,d['Gram']))==m.GRAM,'COMPLETE_POSITIVE_FULL23_GRAM')
    need(d['local_point']==['0','65/64'],'REGISTERED_LOCAL_POINT')
    (u,v),bg=ev.evaluate(Q(65,64));target=[];errors=[]
    for i,(label,row,z) in enumerate(zip(m.LABELS,d['local_components'],u)):
        need(row['label']==label and Q(row['Gram'])==m.GRAM[i],'LOCAL_COMPONENT_ORDER')
        expected=2*m.ELL**2*z.a.c[2];er=sum((2*m.ELL**2*z.d[k].c[2] for k in range(3)),m.I(0));target.append(expected);errors.append(er)
        need(list(map(Q,row['actual_rr']))==expected.row(),'LOCAL_ACTUAL_RECONSTRUCTION:'+label)
        need(list(map(Q,row['paired_error_rr']))==er.row(),'LOCAL_ERROR_RECONSTRUCTION:'+label+':expected='+str(er.row())+':reported='+str(row['paired_error_rr']))
        need(Q(row['weighted_error_square_upper'])==m.GRAM[i]*er.mag()**2,'LOCAL_ALL23_SQUARE_CONTRIBUTION')
        need(all(Q(x)==0 for x in row['channels'][1]),'SHARED_Q_E_DIFFERENCE_ZERO')
    den=m.sqrt_bounds(target);err=m.sqrt_bounds(errors)
    need(list(map(Q,d['actual_target_H_interval']))==den,'LOCAL_DENOMINATOR_FROM_LOCAL_ALL23_COMPONENTS')
    need(list(map(Q,d['new_error_H_interval']))==err,'LOCAL_ERROR_ALL23_NORM')
    need(den[0]>0 and err[1]<=den[0]/10 and err[1]<Q(d['old_error_H_strict_lower'])/2,'LOCAL_TRUE_ERROR_GATES')
    need(Q(d['new_error_over_target_upper'])==err[1]/den[0],'LOCAL_RATIO_ARITHMETIC')
    need(Q(d['old_relative_lower'])>Q(999,1000),'RETAINED_OLD_LOCAL_WITNESS')
    c=d['independent_full_displacement_coverage'];leaves=c['leaves'];need(c['unchanged_exact_zero_region']==['0','1'],'EXACT_ZERO_REGION')
    intervals=[list(map(Q,z['r'])) for z in leaves]
    need(bool(intervals) and intervals[0][0]==1 and intervals[-1][1]==Q(7,4),'FOOTPRINT_ENDPOINTS')
    for left,right in intervals:need(left<right,'POSITIVE_LEAF_WIDTH')
    for aa,bb in zip(intervals,intervals[1:]):need(aa[1]==bb[0],'FULL_FOOTPRINT_COVERAGE:'+str(aa[1])+'!='+str(bb[0]))
    need(len(leaves)==c['leaves_count']<=4096 and all(0<=z['depth']<=16 for z in leaves),'CAMPAIGN_BUDGET')
    need(leaves[0]['tail'] and intervals[0]==[Q(1),Q(65,64)],'FULL_FLAT_TAIL_CELL')
    totals=[[[Q(0)]*23 for _ in range(3)] for _ in range(3)];absolute=[[Q(0)]*23 for _ in range(3)]
    for row in leaves:
        ch=[[[Q(v) for v in rr] for rr in cc] for cc in row['channel_component_caps']];ab=[[Q(v) for v in rr] for rr in row['candidate_absolute_component_caps']]
        need(len(ch)==3 and all(len(cc)==3 and all(len(rr)==23 and min(rr)>=0 for rr in cc) for cc in ch),'COMPLETE_PAIRED_CHANNELS')
        need(all(x==0 for rr in ch[1] for x in rr),'SHARED_Q_E_ZERO_ON_FULL_FOOTPRINT')
        need(len(ab)==3 and all(len(rr)==23 and min(rr)>=0 for rr in ab),'COMPLETE_ABSOLUTE_JETS')
        for k in range(3):
            for j in range(3):totals[k][j]=[max(x,y) for x,y in zip(totals[k][j],ch[k][j])]
        for j in range(3):absolute[j]=[max(x,y) for x,y in zip(absolute[j],ab[j])]
    need(enc(totals)==c['channel_component_suprema'],'RECOMPUTE_ALL_COMPONENT_SUPREMA')
    need(enc(absolute)==c['candidate_absolute_component_suprema'],'RECOMPUTE_ABSOLUTE_COMPONENT_SUPREMA')
    def h(caps):return m.sqrt_bounds([m.I(x) for x in caps])[1]
    caps=[sum((h(totals[k][j]) for k in range(3)),Q(0)) for j in range(3)]
    need(list(map(Q,c['error_H_caps']))==caps,'RECOMPUTE_COMPLETE_PAIRED_ERROR_CAPS')
    need(list(map(Q,c['candidate_absolute_H_caps']))==[h(rr) for rr in absolute],'RECOMPUTE_COMPLETE_ABSOLUTE_CAPS')
    refs=a.source/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference';old=json.loads((refs/'C2.json').read_text());combined=json.loads((refs/'COMBINED.json').read_text())
    limits=[Q(old['initial_error_H_caps']['e0']),Q(old['initial_error_H_caps']['e0_r']),Q(combined['refined_uniform_initial_u_rr_H_upper'])/2]
    need(list(map(Q,c['initial_gate_thresholds']))==limits,'FROZEN_INITIAL_EXACT_THRESHOLDS')
    need(all(x<=y for x,y in zip(caps,limits)) and c['initial_gate_proved_by_this_mesh']==[True]*3,'COMPLETE_INITIAL_GATES')
    need(d['tail_certificate']==enc(m.tail_proof()),'FORMAL_TAIL_POLYNOMIAL_DOMINANCE_PROOF')
    regions=[]
    for name,lo,hi in [('flat_edge_to_old_cutoff',Q(1),Q(65,64)),('former_cutoff_band',Q(65,64),Q(66,64)),('active_interior',Q(66,64),Q(7,4))]:
        subset=[z for z,rr in zip(leaves,intervals) if lo<=rr[0] and rr[1]<=hi];need(subset,'REGION_NONEMPTY')
        rcaps=[[max(Q(z['channel_component_caps'][k][2][i]) for z in subset) for i in range(23)] for k in range(3)]
        regions.append({'region':name,'r':[lo,hi],'rr_H_absolute_upper':sum((h(rr) for rr in rcaps),Q(0)),'relative_ratio':'undefined at the flat support edge; no uniform positive target denominator certified elsewhere'})
    need(m.closure(a.source,ev.pins)==ev.before and a.source_pins.read_bytes()==ev.pinraw,'SOURCE_CLOSURE_MAP_UNCHANGED_AFTER')
    need(a.report.read_bytes()==raw and digest(a.helper)==hp and digest(Path(__file__))==own,'EXECUTED_HELPERS_AND_REPORT_UNCHANGED_AFTER')
    result={'status':'PASS_INDEPENDENT_INITIAL_ARTIFACT_AUDIT','local_ratio_upper':err[1]/den[0],'error_H_caps':caps,'region_rr_absolute_caps':regions,'source_map_sha256':m.SOURCE_MAP_PIN,'report_sha256':hashlib.sha256(raw).hexdigest(),'executed_helper_sha256':hp,'checker_sha256':own,'independent_scientific_local_reconstruction':True,'full_coverage_and_all23_ledger_rebuilt':True,'response_gate_claimed':False}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(enc(result),indent=2,sort_keys=True)+'\n');print(result['status'])
if __name__=='__main__':main()
