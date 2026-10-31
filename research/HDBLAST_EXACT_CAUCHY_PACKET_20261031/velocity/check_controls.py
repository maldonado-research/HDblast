#!/usr/bin/env python3
"""Effective mathematical/ledger rejection witnesses, active under python -O."""
import argparse
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
INITIAL_PIN='219d9fe0aeb161a31da564014dc85d4712ecae339e46485b16242bb724219c89'

def need(ok,message):
    if not ok:raise RuntimeError(message)

def module(path):
    spec=importlib.util.spec_from_file_location('velocity_control_module',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--inject-error');args=ap.parse_args()
    need(not args.output.exists(),'FRESH_OUTPUT')
    import hashlib
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt=json.loads((HERE/'FREEZE_RECEIPT.json').read_bytes())
    need(sha(HERE/'producer.py')==receipt['producer_sha256'],'PRODUCER_AUTH_BEFORE_IMPORT')
    need(sha(HERE/'FREEZE.json')==receipt['freeze_sha256'],'FREEZE_AUTH')
    need(sha(HERE/'INITIAL_NORMAL.json')==INITIAL_PIN,'INITIAL_AUTH')
    before_self=sha(Path(__file__));m=module(HERE/'producer.py')
    frozen=json.loads((HERE/'FREEZE.json').read_bytes());inputs=m.authenticate(args.repo_root,frozen)
    source=m.load(args.repo_root);e=source.ExactPacket(args.repo_root,128)
    report=json.loads((HERE/'INITIAL_NORMAL.json').read_bytes())['campaigns'][-1]
    loc=report['local'];D,DU=map(Q,loc['actual_H_norm_interval']);U=Q(loc['new_error_H_upper']);L=Q(loc['old_error_H_lower'])
    cases={}
    def register(name,predicate,witness,kind):
        need(not predicate,'INEFFECTIVE_CONTROL:'+name)
        cases[name]={'corrupted_predicate':predicate,'witness':witness,'kind':kind}
    register('restored-cutoff',D<=D/10,{'mutated_error_lower':D,'allowed_upper':D/10},'failed local true-error inequality')
    register('restored-smoothed-P1-velocity',L<=D/10,{'old_error_lower':L,'allowed_upper':D/10},'failed local true-error inequality')
    row=loc['rows'][3];need(row['label']=='L:0','MEAN_LABEL')
    n=list(map(Q,row['candidate_v_r']));factor=1/e.ell-1
    mismatch=[factor*x for x in n]
    register('wrong-radial-scale',mismatch[0]<=0<=mismatch[1],{'wrong_minus_correct_L0_v_r':mismatch,'required_scale':e.ell,'corrupt_scale':1},'nonzero formula witness')
    actual=list(map(Q,row['actual_v_r']));low=min(abs(x) for x in actual)
    need(actual[1]<0,'NONZERO_MEAN_TARGET')
    register('zeroed-L0-mean',low<=D/10,{'mutated_component_error_lower':low,'allowed_H_upper':D/10,'Gram':row['Gram']},'failed local true-error inequality')
    s0,s1,old=map(Q,(row['old_left_slope'],row['old_right_slope'],row['old_smoothed_v_r']))
    register('wrong-kernel-half-mass',s0+s1==old,{'correct_half_sum':old,'corrupt_full_sum':s0+s1,'nonzero_difference':(s0+s1)-old},'failed exact convolution identity')
    for name,ch in [('omitted-background-channel',0),('omitted-ODE-remainder-channel',2)]:
        charge=Q(report['paired_channel_H_caps'][ch][1]);total=Q(report['initial_error_H_caps']['v0_r'])
        register(name,total-charge>=total,{'complete_error_ledger':total,'corrupt_ledger':total-charge,'omitted_positive_charge':charge,
            'scope':'Rejects omission of a required proof charge, without claiming a true-error lower bound.'},'failed exact ledger inequality')
    absolute=Q(report['candidate_initial_absolute_H_caps']['v_r']);error=Q(report['initial_error_H_caps']['v0_r'])
    register('absolute-error-interchange',absolute==error,{'absolute_candidate_bound':absolute,'paired_target_error_bound':error,'nonzero_difference':absolute-error},'failed ledger type equality')
    shared=[Q(x) for x in report['paired_channel_H_caps'][1]]
    register('nonzero-shared-qE-difference',shared==[0,1],{'required_shared_channel':shared,'corrupt_shared_channel':[0,1]},'failed exact functional identity')
    register('incomplete-cover',source.cover([list(map(Q,z['r'])) for z in report['leaves'][1:]],Q(0),Q(7,4)),{'deleted_leaf':report['leaves'][0]['r']},'failed complete footprint coverage')
    R=Q(7,4);bg=e.h.background(e.I(R));q=e.packet(R,R);ev=e.expbox(e.ell*e.I(R))
    import math
    E=e.J([(-1)**k*ev/math.factorial(k) for k in range(6)])
    _,defined=e.c._registered({k:v.b for k,v in bg.items()},E,q)
    raw_affine=e.J(q.c[:2]+[0]*4)
    _,wrong=e.c._registered({k:v.b for k,v in bg.items()},E,raw_affine)
    diff=wrong[3].c[0]-defined[3].c[0]
    register('raw-q-affine-exterior',diff.lo<=0<=diff.hi,{'component':'L:0','radius':R,'raw_q_affine_right_minus_completed_left_velocity':diff.row()},'nonzero completed trace mismatch')
    e192=source.ExactPacket(args.repo_root,192);argument=e.ell*Q(65,64)
    low=e.expminus(argument);fine=e192.expminus(argument);mid=(low.lo+low.hi)/2
    register('midpoint-exponential',fine.lo<=mid<=fine.hi,{'argument':argument,'midpoint128':mid,'outward192':fine.row()},'nonzero exact function witness')
    for name,ch in [('unmatched-velocity-value',0),('unmatched-velocity-slope',1)]:
        required=report['endpoint_componentwise_affine_matching'][3]['exterior_minus_interior_exact_jet'][ch]
        register(name,required==1,{'component':'L:0','required_trace_difference':required,'corrupt_trace_difference':1},'failed exact C1 trace identity')
    e.h.authenticate();e192.h.authenticate();need(inputs==m.authenticate(args.repo_root,frozen),'POST_AUTH')
    need(sha(Path(__file__))==before_self and sha(HERE/'producer.py')==receipt['producer_sha256'] and sha(HERE/'INITIAL_NORMAL.json')==INITIAL_PIN,'POST_SELF_AUTH')
    if args.inject_error:
        need(args.inject_error in cases,'UNKNOWN_CONTROL')
        need(cases[args.inject_error]['corrupted_predicate'],'CONTROL_REJECTED:'+args.inject_error)
    out={'status':'PASS_EFFECTIVE_VELOCITY_CONTROLS','initial_sha256':INITIAL_PIN,'producer_sha256':receipt['producer_sha256'],
        'freeze_sha256':receipt['freeze_sha256'],'checker_sha256':before_self,'cases':cases,
        'scope':'Mathematical and proof-ledger rejection witnesses, separate from independent recomputation.'}
    args.output.write_text(json.dumps(m.encode(out),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':out['status'],'control_count':len(cases)},sort_keys=True))

if __name__=='__main__':main()
