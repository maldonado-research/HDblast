#!/usr/bin/env python3
"""Fresh bounded complete exact analytic velocity, preserving every old input."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
SOURCE = 'research/HDBLAST_EXACT_PACKET_INITIAL_20261031'

def need(ok, message):
    if not ok: raise RuntimeError(message)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def encode(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,dict): return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [encode(v) for v in x]
    return x

def authenticate(root, frozen):
    checked = {}
    for rel,pin in sorted(frozen['input_sha256'].items()):
        path=root/rel
        need(path.is_file() and not path.is_symlink() and digest(path)==pin,'AUTH_SOURCE:'+rel)
        checked[rel]=pin
        if rel.endswith('/MANIFEST.json'):
            for item in json.loads(path.read_bytes())['files']:
                name=Path(item['path'])
                need(not name.is_absolute() and '..' not in name.parts,'AUTH_PATH')
                p=path.parent/name
                need(p.is_file() and not p.is_symlink() and p.stat().st_size==item['bytes'] and digest(p)==item['sha256'],'AUTH_MEMBER:'+str(p))
                checked[str(p.relative_to(root))]=item['sha256']
    return checked

def load(root):
    spec=importlib.util.spec_from_file_location('pinned_exact_packet',root/SOURCE/'checks/producer.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def evaluate(e,lo,hi):
    need(0<=lo<=hi<=Q(7,4),'DOMAIN')
    bg=e.h.background(e.I(lo,hi));q=e.packet(lo,hi)
    ev=e.expbox(e.ell*e.I(lo,hi))
    E=e.J([(-1)**k*ev/math.factorial(k) for k in range(6)])
    zero=e.h.constant(0)
    _,v=e.c._registered(bg,e.P(E,E,[zero,zero,zero]),e.P(q,q,[zero,zero,zero]))
    need(len(v)==len(e.weights)==23 and all(z.n>=1 for z in v),'ALL23_VALUE_AND_R')
    actual=[[scale*z.a.c[k] for z in v] for k,scale in enumerate((Q(1),e.ell))]
    nominal=[[scale*z.b.c[k] for z in v] for k,scale in enumerate((Q(1),e.ell))]
    channels=[[[scale*z.d[ch].c[k] for z in v] for k,scale in enumerate((Q(1),e.ell))] for ch in range(3)]
    need(all(v.lo==v.hi==0 for row in channels[1] for v in row),'SHARED_Q_E_ZERO')
    return actual,nominal,channels

def local(e,root,m):
    point=Q(65,64);delta=Q(1,2048)
    profile=json.loads((root/m.PROFILES).read_bytes())
    need(profile['labels']==list(e.c.op.LABELS) and list(map(Q,profile['gram']))==e.weights,'LOCAL_ALL23_GRAM')
    nodes=list(map(Q,profile['nodes']));j=nodes.index(point)
    need(0<j<len(nodes)-1 and nodes[j-1]<point-delta and point+delta<nodes[j+1],'KERNEL_WITHIN_ADJACENT_P1_SEGMENTS')
    vt=[list(map(Q,row)) for row in profile['profiles']['v']]
    left=[(vt[j][i]-vt[j-1][i])/(point-nodes[j-1]) for i in range(23)]
    right=[(vt[j+1][i]-vt[j][i])/(nodes[j+1]-point) for i in range(23)]
    old=[(a+b)/2 for a,b in zip(left,right)]
    actual,nominal,channels=evaluate(e,point,point)
    # Fresh direct actual analytic reconstruction at this point, with neither
    # the P1 profile nor paired errors used as an actual target evaluator.
    bg=e.h.background(e.I(point));q=e.packet(point,point);ev=e.expbox(e.ell*e.I(point))
    E=e.J([(-1)**k*ev/math.factorial(k) for k in range(6)])
    _,direct=e.c._registered({k:z.a for k,z in bg.items()},E,q)
    direct=[e.ell*z.c[1] for z in direct]
    need(all(x.lo==y.lo and x.hi==y.hi for x,y in zip(direct,actual[1])),'DIRECT_ACTUAL_RECONSTRUCTION')
    rows=[];denlower2=denupper2=error2=oldlower2=Q(0)
    for i,(label,w) in enumerate(zip(e.c.op.LABELS,e.weights)):
        a=direct[i];n=nominal[1][i];d=channels[0][1][i]+channels[2][1][i];old_error=a-old[i]
        lo=lambda x: Q(0) if x.lo<=0<=x.hi else min(abs(x.lo),abs(x.hi))
        denlower2+=w*lo(a)**2;denupper2+=w*a.mag()**2;error2+=w*d.mag()**2;oldlower2+=w*lo(old_error)**2
        rows.append({'label':label,'Gram':w,'actual_v_r':a.row(),'candidate_v_r':n.row(),'paired_error_v_r':d.row(),
            'old_left_slope':left[i],'old_right_slope':right[i],'old_smoothed_v_r':old[i],'actual_minus_old_v_r':old_error.row(),
            'actual_H_squared_lower_contribution':w*lo(a)**2,'old_error_H_squared_lower_contribution':w*lo(old_error)**2})
    D=m.sqrt_bound(denlower2,e.bits,False);DU=m.sqrt_bound(denupper2,e.bits);U=m.sqrt_bound(error2,e.bits);L=m.sqrt_bound(oldlower2,e.bits,False)
    need(D>0,'LOCAL_DENOMINATOR_ZERO_STOP')
    return {'point':[Q(0),point],'observable':'initial_velocity_r_full23_H','rows':rows,'actual_H_norm_interval':[D,DU],
        'new_error_H_upper':U,'old_error_H_lower':L,'new_relative_error_H_upper':U/D,
        'old_relative_error_H_lower':L/DU,'relative_gate':U<=D/10,'true_error_improvement_gate':U<L/2,
        'pass':U<=D/10 and U<L/2,'old_derivative_identity':'even normalized kernel; each adjacent constant P1 slope receives exactly mass 1/2',
        'kernel_delta':delta,'adjacent_nodes':nodes[j-1:j+2],'fresh_direct_actual_reconstruction':True,'old_point_table_reused':False}

def campaign(root,bits,m,old,initial,local_only=False):
    e=m.ExactPacket(root,bits);den=e.denominators();loc=local(e,root,m)
    if local_only:return {'bits':bits,'denominators':den,'local':loc,'status':'LOCAL_PASS' if loc['pass'] else 'LOCAL_INCONCLUSIVE'}
    need(loc['pass'],'LOCAL_GATE_UNPROVED_STOP')
    limits=[Q(old['initial_error_H_caps'][k])/2 for k in ('v0','v0_r')]
    seeds=[(Q(0),Q(1)),(Q(1),Q(65,64))]+[(Q(k,64),Q(k+1,64)) for k in range(65,112)]
    need(m.cover(seeds,Q(0),Q(7,4)),'SEED_COVER')
    stack=[(lo,hi,0,i,None) for i,(lo,hi) in reversed(list(enumerate(seeds)))];leaves=[];splits=[]
    while stack:
        lo,hi,depth,seed,parent=stack.pop();actual,nominal,channels=evaluate(e,lo,hi)
        absolute=[[v.mag() for v in row] for row in nominal]
        raw=[[[v.mag() for v in row] for row in ch] for ch in channels]
        caps=raw if parent is None else [[[min(a,b) for a,b in zip(row,oldrow)] for row,oldrow in zip(ch,oldch)] for ch,oldch in zip(raw,parent)]
        norms=[[e.norm(row) for row in ch] for ch in caps]
        missed=any(norms[ch][k]>limits[k]/2 for ch in (0,2) for k in range(2))
        if missed and depth<16 and len(seeds)+len(splits)<4096:
            mid=(lo+hi)/2;splits.append({'r':[lo,hi],'depth':depth,'seed':seed,'midpoint':mid,'component_caps':caps})
            stack.append((mid,hi,depth+1,seed,caps));stack.append((lo,mid,depth+1,seed,caps))
        else:
            leaves.append({'r':[lo,hi],'depth':depth,'seed':seed,'absolute_component_caps':absolute,'absolute_H_caps':[e.norm(row) for row in absolute],
                'paired_channel_component_caps':caps,'paired_channel_H_caps':norms,'threshold_satisfied':not missed,
                'stop_reason':'threshold' if not missed else ('depth_budget' if depth==16 else 'leaf_budget')})
    need(m.cover([x['r'] for x in leaves],Q(0),Q(7,4)),'COMPLETE_EXACT_COVER')
    need(len(leaves)==len(seeds)+len(splits)<=4096 and max(x['depth'] for x in leaves)<=16,'BUDGET')
    chcaps=[[max(z['paired_channel_H_caps'][ch][k] for z in leaves) for k in range(2)] for ch in range(3)]
    errors=[sum(chcaps[ch][k] for ch in range(3)) for k in range(2)]
    abscaps=[max(z['absolute_H_caps'][k] for z in leaves) for k in range(2)]
    _,nominal,channels=evaluate(e,Q(7,4),Q(7,4))
    endpoint=[{'label':label,'candidate_value_r_trace_enclosures':[nominal[k][i].row() for k in range(2)],
        'exterior_minus_interior_exact_jet':[0,0],'paired_error_enclosures':[(channels[0][k][i]+channels[2][k][i]).row() for k in range(2)]}
        for i,label in enumerate(e.c.op.LABELS)]
    ini=initial['campaigns'][-1];inherited_errors={k:Q(v) for k,v in ini['initial_error_H_caps'].items() if k.startswith('e')}
    inherited_absolute={k:Q(v) for k,v in ini['candidate_initial_absolute_H_caps'].items() if k.startswith('u')}
    initialpass=all(a<=b for a,b in zip(errors,limits)) and ini['initial_gate_pass']
    e.h.authenticate()
    return {'bits':bits,'status':'INITIAL_AND_LOCAL_GATES_PASS' if initialpass else 'INCONCLUSIVE_INITIAL_GATES',
        'labels':list(e.c.op.LABELS),'Gram':e.weights,'denominators':den,'local':loc,
        'initial_error_H_caps':inherited_errors|dict(zip(('v0','v0_r'),errors)),
        'candidate_initial_absolute_H_caps':inherited_absolute|dict(zip(('v','v_r'),abscaps)),
        'velocity_initial_gate_limits':dict(zip(('v0','v0_r'),limits)),'initial_gate_pass':initialpass,
        'displacement_initial_gates_inherited_exactly':True,'paired_channel_names':['background_coefficient_difference','exact_shared_q_E_zero','actual_ODE_remainder'],
        'paired_channel_H_caps':chcaps,'endpoint_componentwise_affine_matching':endpoint,'seeds':seeds,
        'leaves':leaves,'splits':splits,'leaf_count':len(leaves),'split_count':len(splits),
        'unresolved_cells':[z['r'] for z in leaves if not z['threshold_satisfied']],
        'budget':{'bits':bits,'max_leaves':4096,'max_depth':16},'response_evolution_performed':False}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--local-only',action='store_true');args=ap.parse_args()
    need(not args.output.exists(),'FRESH_OUTPUT_REQUIRED')
    freeze_path=HERE/'FREEZE.json';freeze_before=digest(freeze_path);self_before=digest(Path(__file__))
    frozen=json.loads(freeze_path.read_bytes());before=authenticate(args.repo_root,frozen)
    m=load(args.repo_root)
    old=json.loads((args.repo_root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/C2.json').read_bytes())
    initial=json.loads((args.repo_root/SOURCE/'reference/INITIAL.json').read_bytes())
    reports=[]
    for bits in (128,192):
        reports.append(campaign(args.repo_root,bits,m,old,initial,args.local_only))
        if reports[-1]['status'] in ('LOCAL_PASS','INITIAL_AND_LOCAL_GATES_PASS'):break
    need(before==authenticate(args.repo_root,frozen),'POST_SOURCE_AUTH')
    need(digest(freeze_path)==freeze_before and digest(Path(__file__))==self_before,'SELF_AND_FREEZE_AUTH')
    out={'schema':'HDBLAST_EXACT_CAUCHY_PACKET_INITIAL_V1','candidate':'W1_exact_Cauchy_packet','status':reports[-1]['status'],
        'freeze_sha256':freeze_before,'producer_sha256':self_before,'source_authentication_before_after':before,
        'campaigns':reports,'target_F_Fx_G_M_H_domain_unchanged':True,'old_candidate_functions_unchanged':True,
        'no_response_evolution_performed':True,'scope':frozen['scope']}
    args.output.write_text(json.dumps(encode(out),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'bits':[r['bits'] for r in reports],'leaf_counts':[r.get('leaf_count') for r in reports],
        'local_D':float(reports[-1]['local']['actual_H_norm_interval'][0]),'local_new_U':float(reports[-1]['local']['new_error_H_upper']),
        'local_old_L':float(reports[-1]['local']['old_error_H_lower']),
        'initial_error_H_caps':{k:float(v) for k,v in reports[-1].get('initial_error_H_caps',{}).items()}},sort_keys=True))

if __name__=='__main__':main()
