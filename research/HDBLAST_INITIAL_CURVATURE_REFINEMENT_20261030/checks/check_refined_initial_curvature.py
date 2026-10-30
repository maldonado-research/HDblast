#!/usr/bin/env python3
"""Unchanged W1_quadratic: refine only the dominant paired cutoff cell.

Every subdivision is rational and covers the same initial footprint.  The
nominal displacement, smoothing widths, outer quadratic continuation, P1
velocity, actual target, and all source coefficients are fixed.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True

CORE_PIN='f2d1bfef82220fc1e385a3bdb6de35487b1ef28c4caee37c8dddb6a2b68cb782'
INITIAL_PIN='bd8ec21a66f59736c13aaf2b5543c7724e04186a8a37e4293d8cba286ec57207'
REPAIR_PIN='4a6e4dcc02c8e0ece0d7452b45b4b2f1e37e8221303d20b6fbe8ec6496004002'
checks=[];controls=[]

def need(ok,label):
    if not ok: raise RuntimeError(label)
    checks.append(label)

def pin(path,digest):
    raw=Path(path).read_bytes()
    need(hashlib.sha256(raw).hexdigest()==digest,'authenticate '+Path(path).name)
    return raw

def enc(x):
    if isinstance(x,Q):return str(x)
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    return x

def valid_cover(rows,left,right):
    return bool(rows) and rows[0][0]==left and rows[-1][1]==right and all(a<b for a,b in rows) and all(a[1]==b[0] for a,b in zip(rows,rows[1:]))

def reject(label,ok):
    need(not ok,'reject '+label);controls.append(label)

ap=argparse.ArgumentParser()
ap.add_argument('--repo-root',required=True,type=Path)
ap.add_argument('--higher-capsule',type=Path)
ap.add_argument('--core-file',type=Path)
args=ap.parse_args()
hi=args.higher_capsule or args.repo_root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030'
core=args.core_file or hi/'checks/initial_second_core.py'
raw=pin(core,CORE_PIN)
old=json.loads(pin(hi/'reference/INITIAL.json',INITIAL_PIN))
repair=json.loads(pin(hi/'reference/REPAIR.json',REPAIR_PIN))
s=importlib.util.spec_from_file_location('authenticated_refinement_core',core)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
h=m.HigherInitial(args.repo_root);c=h.ctx;I=h.I;ell=h.ell
labels=list(c.op.LABELS);weights=c.weights
need(labels==old['labels']==repair['labels'] and weights==list(map(Q,old['gram']))==list(map(Q,repair['gram'])),'same all23 positive physical Gram')
need(len(labels)==23 and min(weights)>0,'positive complete H')
need(old['status']==repair['status']=='PASS' and repair['velocity_and_its_initial_error_bound_unchanged'] is True,'successful unchanged quadratic candidate')
need(old['initial_data_unchanged'] is True and repair['underlying_actual_target_unchanged'] is True,'same underlying actual target')
delta=Q(old['fixed_displacement_kernel_width']);R=Q(7,4)
need(delta==Q(1,2**30) and Q(old['fixed_velocity_kernel_width'])==Q(1,2048),'same fixed kernels')
left,right=Q(65,64),Q(66,64)
partition=old['continuous']['partition']
need(len(partition)==48 and valid_cover([tuple(map(Q,row['r'])) for row in partition],Q(1),R),'authenticated old complete48 cell cover')
replaced=[row for row in partition if list(map(Q,row['r']))==[left,right]]
need(len(replaced)==1,'unique dominant cutoff cell')
retained=[row for row in partition if row is not replaced[0]]
need(len(retained)==47,'all47 unaffected cells retained')
old_cutoff=[list(map(Q,row)) for row in replaced[0]['channel_H_caps']]
need(all(max(Q(row['channel_H_caps'][1][0]) for row in retained)<old_cutoff[1][0] for _ in [0]),'cutoff cell controls old paired channel')

# Stop once every paired cutoff-displacement channel norm is <=1e28.
# The entire leaf tree is retained in the report; no point sampling is used.
threshold=Q('1e28'); leaves=[];splits=[];evaluated=0

def evaluate(a,b,depth):
    global evaluated
    u,v=h.paired(I(a,b))
    caps=[[[ (2*ell**2*z.d[ch].c[2]).mag() for z in u],
           [ (ell*z.d[ch].c[1]).mag() for z in v]] for ch in range(3)]
    norms=[[c.norm(row) for row in channel] for channel in caps]
    evaluated+=1
    return {'r':[a,b],'depth':depth,'paired_channel_component_caps':caps,'channel_H_caps':norms}

stack=[(left,right,0)]
while stack:
    a,b,depth=stack.pop()
    row=evaluate(a,b,depth)
    if max(z[0] for z in row['channel_H_caps'])>threshold:
        need(depth<10,'bounded exact dyadic refinement depth')
        mid=(a+b)/2
        splits.append({'r':[a,b],'midpoint':mid,'depth':depth})
        stack.append((mid,b,depth+1));stack.append((a,mid,depth+1))
    else:leaves.append(row)
need(valid_cover([tuple(row['r']) for row in leaves],left,right),'gap-free disjoint exact refined cutoff cover')
need(all(max(z[0] for z in row['channel_H_caps'])<=threshold for row in leaves),'every refined cutoff leaf passes threshold')
need(evaluated==len(leaves)+len(splits) and len(leaves)==len(splits)+1,'complete binary refinement tree')
whole=[{'r':list(map(Q,row['r'])),'channel_H_caps':[list(map(Q,z)) for z in row['channel_H_caps']],'role':'retained'} for row in retained]
whole +=[{'r':row['r'],'channel_H_caps':row['channel_H_caps'],'role':'refined'} for row in leaves]
whole.sort(key=lambda row:row['r'][0])
need(valid_cover([tuple(row['r']) for row in whole],Q(1),R),'complete active-domain cover after replacement')

# Outside [1,R] the retained actual active data and approximation vanish on
# the inner collar. Smoothing includes that collar and the same quadratic
# exterior. The original proof controls it by the same global M3 and jump
# charges, so we do not assume the paired comparison extends past R.
jumps=old['curvature_jumps']
need([Q(row['r']) for row in jumps]==[left,right,R],'all old joins accounted for')
need(all(Q(z)==0 for z in jumps[0]['jump_abs_component_caps']),'first join curvature jump exactly zero')
inner=[max(Q(row['jump_abs_component_caps'][i]) for row in jumps[:2]) for i in range(23)]
M3=list(map(Q,old['piecewise_nominal_M3u']))
smoothing=[z*Q(35,128)*delta+j/2 for z,j in zip(M3,inner)]
smoothingH=c.norm(smoothing)
old_channel=[list(map(Q,row)) for row in old['continuous']['paired_channel_H_caps']]
need(sum(z[0] for z in old_channel)+smoothingH==Q(repair['new_uniform_initial_u_rr_H_upper']),'same PR100 repaired global smoothing charge')
# A channel H norm is bounded uniformly by the maximum of its cell H caps.
# The three telescoping channels are then added by the triangle inequality.
channelH=[[max(row['channel_H_caps'][ch][j] for row in whole) for j in range(2)] for ch in range(3)]
pairedH=sum(row[0] for row in channelH)
new_u2=pairedH+smoothingH
old_u2=Q(repair['new_uniform_initial_u_rr_H_upper'])
need(new_u2<old_u2,'strictly tighter certificate for unchanged candidate')
need(new_u2<Q('1.5e28'),'readable strict refined initial-curvature bound')
need(old_u2>Q(400000)*new_u2,'more than400000fold reduction of upper bound only')
# v_r already has a much larger retained P1/convolution charge. Preserve its
# published complete cap rather than claim the local paired refinement as
# an improvement of the actual velocity error or any absolute data norm.
vr=Q(old['continuous']['actual_minus_smoothed_P1_v_r_H_upper'])
need(vr<Q('9.672e24'),'retained full velocity-gradient error bound')
cover=[tuple(row['r']) for row in leaves]
reject('missing first leaf',valid_cover(cover[1:],left,right))
reject('overlapping duplicate leaf',valid_cover([cover[0]]+cover,left,right))
reject('reversed leaf',valid_cover([(cover[0][1],cover[0][0])]+cover[1:],left,right))
reject('retain no smoothing allowance',new_u2==pairedH)
reject('maximum instead of sum of telescoping channels',new_u2==max(z[0] for z in channelH)+smoothingH)
reject('old unresolved cutoff cell claimed refined',old_cutoff[1][0]<=threshold)
h.authenticate();pin(core,CORE_PIN);pin(hi/'reference/INITIAL.json',INITIAL_PIN);pin(hi/'reference/REPAIR.json',REPAIR_PIN)
print(json.dumps(enc({'status':'PASS_UNCHANGED_QUADRATIC_INITIAL_CURVATURE_REFINEMENT',
 'candidate':'W1_quadratic (W1_dagger), unchanged PR100 original-source finite candidate',
 'target_unchanged':True,'initial_data_unchanged':True,'F_G_and_coefficients_unchanged':True,
 'same_candidate_W1_quadratic':True,'candidate_changed':False,
 'refined_u_rr_H_cap_is_upper_bound':True,
 'bound_quantity':'sup_H_norm_actual_minus_smoothed_initial_u_rr',
 'labels':labels,'gram':weights,'Gram':weights,'ell':ell,'active_footprint':[Q(1),R],
 'fixed_displacement_kernel_width':delta,'fixed_velocity_kernel_width':Q(1,2048),
 'input_sha256':{'core':CORE_PIN,'higher_initial':INITIAL_PIN,'initial_repair':REPAIR_PIN},
 'refinement':{'replaced_old_cell':[left,right],'retained_old_cells':len(retained),'leaf_threshold_H':threshold,
  'number_of_leaves':len(leaves),'number_of_splits':len(splits),'number_of_evaluations':evaluated,
  'splits':splits,'leaves':leaves},
 'complete_cover':whole,'paired_channel_H_caps':channelH,
 'paired_u_rr_H_upper':pairedH,'unchanged_nominal_M3u_component_caps':M3,
 'unchanged_inner_jump_component_caps':inner,'unchanged_curvature_smoothing_component_caps':smoothing,
 'unchanged_curvature_smoothing_H_cap':smoothingH,
 'old_uniform_initial_u_rr_H_upper':old_u2,
 'new_uniform_initial_u_rr_H_upper':new_u2,
 'refined_uniform_initial_u_rr_H_upper':new_u2,
 'new_uniform_initial_u_rr_H_readable_strict_upper':'1.5e28',
 'upper_certificate_reduction_factor':old_u2/new_u2,
 'actual_error_reduction_claimed':False,
 'unchanged_full_initial_v_r_error_H_upper':vr,
 'scope':'Tighter continuous second initial-data error upper certificate for unchanged W1_quadratic. No new candidate, lower actual error claim, C2 transfer by this report alone, relative physical accuracy, nonlinear spacetime or cosmological result.',
 'named_checks':checks,'wrong_formula_controls_rejected':controls}),sort_keys=True))
