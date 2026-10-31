#!/usr/bin/env python3
"""Independent arithmetic/provenance audit of fresh velocity certificates."""
import argparse, hashlib, importlib.util, json, math, sys
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load_helper():
    s=importlib.util.spec_from_file_location('own_independent_velocity',HERE/'check_velocity_independent.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sqnorm(row,g):return sum((w*Q(x)**2 for w,x in zip(g,row)),Q(0))
def exact_upper(z,bits=128):
    grid=1<<bits;k=math.isqrt((z*grid*grid).__floor__())
    if Q(k*k,grid*grid)<z:k+=1
    return Q(k,grid)
def check_own(x,m,root):
    need(x['candidate']=='W1_exact_Cauchy_packet' and x['own_direct_all23_formula'] and not x['producer_scientific_imports'],'OWN_CANDIDATE_INDEPENDENCE')
    need(x['executed_helper_sha256']==sha(HERE/'check_velocity_independent.py'),'OWN_EXECUTED_HELPER_PIN')
    pins=json.loads((HERE.parent/'source-pins.json').read_bytes())
    need(x['source_pins_sha256']==m.SOURCE_MAP_PIN and x['source_sha256_before_and_after']==pins,'OWN_AUTHENTICATED_SOURCE_MAP')
    m.closure(root,pins)
    loc=x['local'];need(len(loc['rows'])==23,'OWN_LOCAL_ALL23')
    profilepath=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_PROFILES.json'
    profile=json.loads(profilepath.read_bytes());nodes=list(map(Q,profile['nodes']));table=[[Q(z) for z in row] for row in profile['profiles']['v']]
    need(loc['old_all_nodes_and_full_velocity_table_sha256']==sha(profilepath),'OWN_OLD_FULL_TABLE_HASH')
    delta=Q(loc['kernel_width']);point=Q(65,64)
    need(delta==Q(1,2048) and nodes[0]<point-delta<point<point+delta<nodes[2],'OWN_OLD_KERNEL_ISOLATION')
    kp=[Q(35,32)*z for z in [1,0,-3,0,3,0,-1]];mass=sum((z/Q(i+1) for i,z in enumerate(kp)),Q(0))
    need(list(map(Q,loc['kernel_polynomial']))==kp and Q(loc['kernel_half_mass'])==mass==Q(1,2),'OWN_KERNEL_EXACT_HALF')
    targets=[];errors=[];old_errors=[]
    for i,(row,label,g) in enumerate(zip(loc['rows'],m.LABELS,m.GRAM)):
        need(row['label']==label and Q(row['Gram'])==g,'OWN_LOCAL_LABEL_AND_GRAM')
        a=m.I(*row['actual_v_r']);d=m.I(*row['paired_error_v_r']);channels=[m.I(*z) for z in row['paired_channels']]
        need(len(channels)==3 and channels[1].lo==channels[1].hi==0,'OWN_ZERO_Q_E_CHANNEL')
        expected=sum(channels,m.I(0));need(d.row()==expected.row(),'OWN_PAIRED_CHANNEL_SUM')
        l=(table[1][i]-table[0][i])/(nodes[1]-nodes[0]);r=(table[2][i]-table[1][i])/(nodes[2]-nodes[1]);old=(l+r)/2
        need(Q(row['old_left_secant'])==l and Q(row['old_right_secant'])==r and Q(row['old_smoothed_v_r'])==old,'OWN_EXACT_ADJACENT_SECANTS')
        oe=a-old;need(oe.row()==m.I(*row['actual_minus_old_smoothed_v_r']).row(),'OWN_SIGNED_OLD_ERROR')
        targets.append(a);errors.append(d);old_errors.append(oe)
    den=m.sqrt_bounds(targets);err=m.sqrt_bounds(errors);olderr=m.sqrt_bounds(old_errors)
    need(list(map(Q,loc['actual_v_r_H_interval']))==den and list(map(Q,loc['new_error_v_r_H_interval']))==err and list(map(Q,loc['old_error_v_r_H_interval']))==olderr,'OWN_LOCAL_NORMS_RECONSTRUCT')
    need(den[0]>0 and err[1]<=den[0]/10 and err[1]<olderr[0]/2,'OWN_LOCAL_GATES')
    need(Q(loc['new_over_target_upper'])==err[1]/den[0] and Q(loc['new_upper_over_old_lower'])==err[1]/olderr[0],'OWN_LOCAL_RATIOS')
    c=x['coverage'];leaves=c['leaves'];need(len(leaves)==c['leaves_count']<=4096 and max(l['depth'] for l in leaves)<=16,'OWN_BUDGET')
    need(list(map(Q,c['exact_zero_region']))==[Q(0),Q(1)],'OWN_FLAT_ZERO_REGION')
    spans=[list(map(Q,l['r'])) for l in leaves]
    need(spans[0][0]==1 and spans[-1][1]==Q(7,4) and all(a<b for a,b in spans) and all(a[1]==b[0] for a,b in zip(spans,spans[1:])),'OWN_EXACT_COVER')
    for l in leaves:
        need(len(l['channel_component_caps'])==3 and all(len(row)==23 for ch in l['channel_component_caps'] for row in ch),'OWN_COVER_ALL23')
        need(all(Q(z)==0 for row in l['channel_component_caps'][1] for z in row),'OWN_GLOBAL_ZERO_Q_E')
        need(all(Q(z)>=0 for ch in l['channel_component_caps'] for row in ch for z in row),'OWN_NONNEGATIVE_COMPONENT_CAPS')
        for key in ('rho','p'):
            z=m.I(*l['actual_background_denominators'][key]);need(z.lo>0 or z.hi<0,'OWN_DENOMINATOR_SEPARATION')
    channel=[[[max(Q(l['channel_component_caps'][k][j][i]) for l in leaves) for i in range(23)] for j in range(2)] for k in range(3)]
    absolute=[[max(Q(l['nominal_absolute_component_caps'][j][i]) for l in leaves) for i in range(23)] for j in range(2)]
    hh=[[m.sqrt_bounds([m.I(z) for z in row])[1] for row in ch] for ch in channel]
    errors=[sum((hh[k][j] for k in range(3)),Q(0)) for j in range(2)]
    aa=[m.sqrt_bounds([m.I(z) for z in row])[1] for row in absolute]
    need(channel==[[list(map(Q,row)) for row in ch] for ch in c['channel_component_suprema']],'OWN_CHANNEL_SUPREMA')
    need(absolute==[list(map(Q,row)) for row in c['nominal_absolute_component_suprema']],'OWN_ABSOLUTE_SUPREMA')
    need(errors==list(map(Q,c['error_H_caps'])) and aa==list(map(Q,c['nominal_absolute_H_caps'])),'OWN_GLOBAL_NORM_RECONSTRUCT')
    old=json.loads((root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json').read_bytes())
    limits=[Q(old['initial_error_H_caps'][k])/2 for k in ('v0','v0_r')]
    need(all(a<=b for a,b in zip(errors,limits)),'OWN_STRONG_FACTOR_TWO_INITIAL_GATES')
    return {'coverage_leaves':len(leaves),'factor_two_initial_gates':True,'absolute_caps':[str(z) for z in aa],'error_caps':[str(z) for z in errors]}
def check_native(x,m,root,producer_root):
    pins=json.loads((HERE.parent/'producer-pins.json').read_bytes())
    for rel,pin in pins.items():need(sha(producer_root/rel)==pin,'NATIVE_PRODUCER_PIN:'+rel)
    need(x['freeze_sha256']==pins['FREEZE.json'] and x['producer_sha256']==pins['producer.py'],'NATIVE_EXECUTED_FILES')
    m.closure(root,x['source_authentication_before_after'])
    c=x['campaigns'][-1];need(x['candidate']=='W1_exact_Cauchy_packet' and c['bits']==128,'NATIVE_CANDIDATE_AND_PRECISION')
    need(c['labels']==m.LABELS and list(map(Q,c['Gram']))==m.GRAM,'NATIVE_ALL23_GRAM')
    leaves=c['leaves'];spans=[list(map(Q,l['r'])) for l in leaves]
    need(len(leaves)==c['leaf_count']<=4096 and max(l['depth'] for l in leaves)<=16,'NATIVE_BUDGET')
    need(spans[0][0]==0 and spans[-1][1]==Q(7,4) and all(a<b for a,b in spans) and all(a[1]==b[0] for a,b in zip(spans,spans[1:])),'NATIVE_EXACT_COVER')
    for l in leaves:
        need(all(len(row)==23 for row in l['absolute_component_caps']),'NATIVE_ABSOLUTE_ALL23')
        need(all(Q(v)==0 for row in l['paired_channel_component_caps'][1] for v in row),'NATIVE_EXACT_SHARED_Q_E_ZERO')
        for row,cap in zip(l['absolute_component_caps'],l['absolute_H_caps']):need(Q(cap)==exact_upper(sqnorm(row,m.GRAM)),'NATIVE_ABSOLUTE_LEAF_H')
        for rows,caps in zip(l['paired_channel_component_caps'],l['paired_channel_H_caps']):
            for row,cap in zip(rows,caps):need(Q(cap)==exact_upper(sqnorm(row,m.GRAM)),'NATIVE_CHANNEL_LEAF_H')
    abs_caps=[max(Q(l['absolute_H_caps'][j]) for l in leaves) for j in range(2)]
    channels=[[max(Q(l['paired_channel_H_caps'][k][j]) for l in leaves) for j in range(2)] for k in range(3)]
    errors=[sum((channels[k][j] for k in range(3)),Q(0)) for j in range(2)]
    need(abs_caps==[Q(c['candidate_initial_absolute_H_caps'][k]) for k in ('v','v_r')],'NATIVE_ABSOLUTE_AGGREGATION')
    need(errors==[Q(c['initial_error_H_caps'][k]) for k in ('v0','v0_r')],'NATIVE_ERROR_AGGREGATION')
    old=json.loads((root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json').read_bytes())
    need(all(e<=Q(old['initial_error_H_caps'][k])/2 for e,k in zip(errors,('v0','v0_r'))),'NATIVE_INITIAL_GATES')
    rows=c['endpoint_componentwise_affine_matching'];need(len(rows)==23 and [r['label'] for r in rows]==m.LABELS and all(r['exterior_minus_interior_exact_jet']==[0,0] for r in rows),'NATIVE_COMPLETED_COMPONENT_TRACE_MATCHING')
    return {'leaves':len(leaves),'source_members':len(x['source_authentication_before_after']),'aggregation_and_factor_two_gates':True,'scientific_status':'Native interval algorithm source reviewed; own direct all23 independent campaign supplies separately valid larger caps.'}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--own',required=True,type=Path);ap.add_argument('--native',required=True,type=Path);ap.add_argument('--producer-root',required=True,type=Path);ap.add_argument('--out',required=True,type=Path);args=ap.parse_args()
    need(not args.out.exists(),'FRESH_OUTPUT');before={str(p):sha(p) for p in (Path(__file__),args.own,args.native)};m=load_helper()
    own=check_own(json.loads(args.own.read_bytes()),m,args.source);native=check_native(json.loads(args.native.read_bytes()),m,args.source,args.producer_root)
    need(all(sha(Path(p))==h for p,h in before.items()),'AUDIT_INPUT_SELF_UNCHANGED_AFTER')
    report={'status':'PASS_INDEPENDENT_VELOCITY_REPORT_AND_PROVENANCE_AUDIT','own':own,'native':native,'input_and_self_sha256_before_after':before}
    args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(report['status'])
if __name__=='__main__':main()
