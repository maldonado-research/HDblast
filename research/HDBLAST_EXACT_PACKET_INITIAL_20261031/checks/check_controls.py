#!/usr/bin/env python3
"""Executable exact failed-inequality controls for the frozen experiment.

Each --inject-error process must exit nonzero before creating its output.
Ledger controls reject unproved deletions, not infer unknown true errors.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
PRODUCER_PIN = '1d7f535ecd5a380cca98408531841ab207b336fd56a249262e7671d06d591896'
INITIAL_PIN = '692ac13f72354c526a0108ddae32800a9f22efb1a8afe4e987cbf403b5a296c3'
CASES = ('restored-cutoff','restored-convolution','midpoint-exponential','wrong-radial-scale',
         'zeroed-mean-override','omitted-background-channel','omitted-remainder-channel',
         'velocity-absolute-error-interchange','omitted-velocity-charge','global-denominator-at-local-point')

def need(ok, message):
    if not ok: raise RuntimeError(message)

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def encode(x):
    if isinstance(x,Q): return str(x)
    if isinstance(x,dict): return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [encode(v) for v in x]
    return x

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--initial-report',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--inject-error',choices=CASES)
    args=ap.parse_args()
    need(not args.output.exists(),'fresh absent control output required')
    before=digest(Path(__file__))
    producer=HERE/'producer.py'
    need(digest(producer)==PRODUCER_PIN,'authenticate producer before import')
    need(digest(args.initial_report)==INITIAL_PIN,'authenticate complete initial report')
    spec=importlib.util.spec_from_file_location('control_authenticated_exact_packet',producer)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    need(digest(HERE/'FREEZE.json')==mod.FREEZE_PIN,'authenticate frozen formula map')
    frozen=json.loads((HERE/'FREEZE.json').read_bytes())
    inputs=mod.authenticate(args.repo_root,frozen)
    data=json.loads(args.initial_report.read_bytes());report=data['campaigns'][0]
    need(report['status']=='INITIAL_AND_LOCAL_GATES_PASS','successful frozen initial input')
    local=report['local'];Dlo,Dup=map(Q,local['actual_H_norm_interval'])
    U=Q(local['new_error_H_upper']);ell=Q(1,10**6)
    cases={}
    def register(name,predicate,witness,kind):
        need(not predicate,'ineffective corruption '+name)
        cases[name]={'corrupted_predicate':predicate,'exact_witness':witness,'kind':kind}

    # A restored zero cutoff at P* has error equal to the target, whose
    # independently enclosed H norm has positive lower endpoint Dlo.
    register('restored-cutoff',Dlo<=Dlo/10,
        {'mutated_error_H_lower':Dlo,'allowed_error_H_upper':Dlo/10,'difference':9*Dlo/10},'failed local true-error inequality')

    # At R+2 delta the entire convolution window is in the completed
    # quadratic exterior. The exact even-kernel increment is kappa delta²/18.
    kappa=list(map(Q,report['endpoint_componentwise_quadratic_matching'][12]['candidate_trace_enclosures'][2]))
    delta=Q(1,2**30);shift=[x*delta*delta/18 for x in kappa]
    need(shift[1]<0,'nonzero completed scalar exterior curvature')
    register('restored-convolution',shift[0]<=0<=shift[1],
        {'point_r':Q(7,4)+2*delta,'kernel_window_entirely_in_completed_exterior':True,'component':'w:0',
         'convolved_minus_defined_candidate_value':shift,'positive_distance_from_zero':-shift[1]},'nonzero formula witness')

    # A midpoint is a rational guess, not the exact exponential. A fresh
    # 192-bit enclosure excludes the midpoint of the 128-bit enclosure.
    e128=mod.ExactPacket(args.repo_root,128)
    e192=mod.ExactPacket(args.repo_root,192)
    argument=ell*Q(65,64)
    low=e128.expminus(argument);fine=e192.expminus(argument)
    midpoint=(low.lo+low.hi)/2
    distance=max(fine.lo-midpoint,midpoint-fine.hi)
    register('midpoint-exponential',fine.lo<=midpoint<=fine.hi,
        {'argument':argument,'outward128':low.row(),'outward192':fine.row(),'claimed_exact_midpoint':midpoint,
         'certified_distance_from_exact_exponential_lower':distance},'nonzero exact-function witness')

    nominal=list(map(Q,local['rows'][12]['candidate_u_rr']))
    factor=1/(ell*ell)-1
    mismatch=[factor*x for x in nominal]
    need(mismatch[0]>0,'nonzero scalar radial-scale witness')
    register('wrong-radial-scale',mismatch[0]<=0<=mismatch[1],
        {'correct_second_Taylor_scale':2*ell*ell,'corrupt_scale':2,'wrong_minus_correct_scalar_rr':mismatch},'nonzero formula witness')

    actual=list(map(Q,local['rows'][12]['actual_u_rr']))
    need(actual[0]>0,'positive scalar mean target')
    register('zeroed-mean-override',actual[0]<=Dlo/10,
        {'label':'w:0','mutated_scalar_error_lower':actual[0],'local_allowed_H_error':Dlo/10,
         'note':'Simple removal of the override loop is algebraically neutral because all five packed means equal their overrides; zeroing the scalar mean is an effective corruption.'},'failed local true-error inequality')

    channels=[[Q(v) for v in row] for row in report['paired_channel_H_caps']]
    total=sum(row[2] for row in channels)
    for case,ch in [('omitted-background-channel',0),('omitted-remainder-channel',2)]:
        mutated=total-channels[ch][2]
        register(case,mutated>=total,{'complete_rr_ledger':total,'mutated_rr_ledger':mutated,
            'omitted_positive_charge':channels[ch][2],'scope':'Rejects a missing proof charge; does not claim this lower ledger value falsifies the unknown true error.'},'failed exact ledger inequality')

    old_velocity=Q(report['initial_error_H_caps']['v0_r'])
    absolute_velocity=Q(report['candidate_initial_absolute_H_caps']['v_r'])
    register('velocity-absolute-error-interchange',old_velocity==absolute_velocity,
        {'required_unchanged_error':old_velocity,'wrong_absolute_candidate_value':absolute_velocity,
         'nonzero_difference':old_velocity-absolute_velocity},'failed unchanged-ledger equality')
    old_velocity0=Q(report['initial_error_H_caps']['v0'])
    register('omitted-velocity-charge',old_velocity0==0,
        {'required_unchanged_velocity_value_error':old_velocity0,'corrupted_charge':0},'failed unchanged-ledger equality')

    normalization=json.loads((args.repo_root/'research/HDBLAST_REFINED_ERROR_TRANSFER_20261030/reference/NORMALIZATION.json').read_bytes())
    global_scale=Q(normalization['target_scale_strict_lower'])
    register('global-denominator-at-local-point',global_scale<=Dup,
        {'wrong_global_scale':global_scale,'actual_local_H_upper':Dup,'strict_difference':global_scale-Dup},'failed local denominator inclusion')
    need(set(cases)==set(CASES),'all controls exercised')
    e128.h.authenticate();e192.h.authenticate()
    need(inputs==mod.authenticate(args.repo_root,frozen),'all source files authenticated after controls')
    need(digest(Path(__file__))==before and digest(producer)==PRODUCER_PIN and digest(args.initial_report)==INITIAL_PIN,
         'checker helper report unchanged after evaluation')
    if args.inject_error:
        selected=cases[args.inject_error]
        need(selected['corrupted_predicate'], 'CONTROL_REJECTED '+args.inject_error+' '+json.dumps(encode(selected),sort_keys=True))
    result={'status':'PASS_EXACT_PACKET_EFFECTIVE_CONTROL_WITNESSES','producer_sha256':PRODUCER_PIN,
        'initial_sha256':INITIAL_PIN,'checker_sha256':before,'cases':cases,
        'mean_override_structural_identity':{'Y2_mean':Q(1,6),'mean_indices':[0,3,6,9,12],
            'identities':['0=0','0=0','-2*u^2*(1/6)=-u^2/3','fs[0]=fs[0]','-3*fy[0]/(2*p)=-3*fy[0]/(2*p)'],
            'simple_loop_deletion_is_effective_mathematical_corruption':False,
            'registered_loop_bytes_still_authenticated':True},
        'scope':'Exact mathematical or proof-ledger rejection witnesses. These controls do not replace independent evaluation.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(encode(result),sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'controls':len(cases)},sort_keys=True))

if __name__=='__main__':main()
