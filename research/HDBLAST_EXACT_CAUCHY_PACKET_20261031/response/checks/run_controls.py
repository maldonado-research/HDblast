#!/usr/bin/env python3
"""Fresh subprocess controls for the frozen exact rational response adapter."""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def need(ok,message):
    if not ok:
        raise RuntimeError(message)

def dump(p,data):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--caps',type=Path,required=True)
    ap.add_argument('--response',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    out=args.output_dir; need(not out.exists(),'fresh controls dir'); out.mkdir(parents=True)
    base=Path(__file__).resolve().parent
    adapter=base/'response_adapter.py'; freeze=base/'FREEZE.json'
    f=json.loads(freeze.read_text()); response=json.loads(args.response.read_text())
    cap_report=json.loads(args.caps.read_text())
    freeze_sha=sha(freeze); cap_sha=response['refined_cap_metadata']['sha256']
    need(sha(args.caps)==cap_sha,'canonical cap input pin')
    common=['--repo-root',str(args.repo_root),'--freeze',str(freeze),
            '--freeze-sha256',freeze_sha,'--caps',str(args.caps),'--caps-sha256',cap_sha]
    records=[]
    def run(name, flags=None, adapter_path=adapter, expect=None):
        for mode,extra in [('normal',[]),('optimized',['-O'])]:
            destination=out/(name+'-'+mode+'.json')
            cmd=[sys.executable]+extra+[str(adapter_path)]+(common if flags is None else flags)+['--output',str(destination)]
            r=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
            need(r.returncode!=0 and not destination.exists(),'corruption accepted: '+name+' '+mode)
            need(expect is None or expect in r.stderr,'wrong rejection: '+name+' '+r.stderr[-1000:])
            portable_stderr=r.stderr
            for path,role in [(out,'<controls>'),(args.repo_root,'<repository>'),(args.caps.parent,'<absolute-proof>'),(base.parent,'<response>')]:
                portable_stderr=portable_stderr.replace(str(path),role)
            records.append({'control':name,'mode':mode,'returncode':r.returncode,
                'success_artifact_absent':True,'expected_rejection':expect,
                'stderr_sha256':hashlib.sha256(portable_stderr.encode()).hexdigest(),
                'last_error':portable_stderr.strip().splitlines()[-1]})
    faults=('absolute_error_curvature_interchange','omit_full_Robin_source','omit_f00',
        'omit_r00','omit_original_Fx','omit_feedback_inverse',
        'omit_derivative_commutator','lose_r1_initial_cancellation')
    for fault in faults:
        run(fault,common+['--fault',fault],expect='CORRUPTION_REJECTED:'+fault)
    def cap_mutation(name, mutate, rehash=True, expected=None):
        data=deepcopy(cap_report); mutate(data); p=out/'fixtures'/(name+'.json'); dump(p,data)
        flags=list(common); flags[flags.index('--caps')+1]=str(p)
        flags[flags.index('--caps-sha256')+1]=sha(p) if rehash else cap_sha
        run(name,flags,expect=expected)
    for key in ('e0','e0_r','e0_rr','v0','v0_r'):
        cap_mutation('changed_frozen_error_'+key,
            lambda r,key=key:r['initial_error_H_caps'].__setitem__(key,str(Q(r['initial_error_H_caps'][key])+1)),
            expected='cannot change target-error caps')
    for key in ('v','v_r'):
        cap_mutation('changed_absolute_velocity_'+key,
            lambda r,key=key:r['candidate_initial_absolute_H_caps'].__setitem__(key,str(Q(r['candidate_initial_absolute_H_caps'][key])-1)),
            expected='unchanged absolute velocity caps')
    for key in ('u','u_r','u_rr'):
        cap_mutation('worse_absolute_'+key,
            lambda r,key=key:r['candidate_initial_absolute_H_caps'].__setitem__(key,str(Q(response['baseline_absolute_H_caps'][key])+1)),
            expected='no worse absolute displacement charges')
    cap_mutation('wrong_candidate',lambda r:r.__setitem__('candidate','W1_exact_Cauchy_packet'),expected='refined same candidate identity')
    cap_mutation('missing_absolute_mode',lambda r:r['candidate_initial_absolute_H_caps'].pop('u_rr'),expected='exact cap schema')
    cap_mutation('changed_table_dominant_rr_row',
        lambda r:r['nodes'][8]['raw_component_caps'][2].__setitem__(12,'0'),
        rehash=False,expected='refined caps report pin')
    cap_mutation('skipped_table_leaf',lambda r:r['leaves'].pop(),rehash=False,expected='refined caps report pin')
    for name,key in [('report_wrong_rr','rr'),('report_wrong_xr','xr')]:
        r=deepcopy(response); r['ledger']['second_coordinate_error_H_caps'][key]=str(Q(r['ledger']['second_coordinate_error_H_caps'][key])+1)
        p=out/'fixtures'/(name+'.json'); dump(p,r)
        run(name,common+['--verify-report',str(p)],expect='verified report every exact node')
    r=deepcopy(response); r['same_function_true_error_unchanged']=False
    p=out/'fixtures/report_wrong_true_error_claim.json'; dump(p,r)
    run('report_wrong_true_error_claim',common+['--verify-report',str(p)],expect='verified report same true error')
    # Minimal exact-byte source closure; no scientific input tree is edited.
    mirror=out/'fixtures/source-root'
    for relative in f['source_pins']:
        p=mirror/relative; p.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(args.repo_root/relative,p)
    for name,relative,mutate in [
        ('source_actual_Fx','research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json',lambda r:r.__setitem__('genuine_scaled_F_x_H_cap','0')),
        ('source_baseline_C2','research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/RESPONSE.json',lambda r:r['ledger']['second_coordinate_error_H_caps'].__setitem__('rr','0')),
        ('source_initial_error','research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/INITIAL.json',lambda r:r['campaigns'][0]['initial_error_H_caps'].__setitem__('e0_rr','0')),
    ]:
        p=mirror/relative; original=p.read_bytes(); data=json.loads(original); mutate(data); dump(p,data)
        flags=list(common); flags[flags.index('--repo-root')+1]=str(mirror)
        run(name,flags,expect='selected source pin: '+relative)
        p.write_bytes(original)
    changed_adapter=out/'fixtures/changed_adapter.py'
    changed_adapter.write_text(adapter.read_text().replace("k['Fx']", "Q(0)", 1))
    run('changed_recipe_source',adapter_path=changed_adapter,expect='frozen ledger recipe source')
    changed_freeze=out/'fixtures/changed_FREEZE.json'; altered=deepcopy(f); altered['fixed_inputs'].remove('Fx'); dump(changed_freeze,altered)
    flags=list(common); flags[flags.index('--freeze')+1]=str(changed_freeze)
    run('changed_freeze',flags,expect='frozen recipe hash')
    # A valid hash at a symlinked ancestor still fails before source use.
    link=out/'fixtures/symlinked-source'; link.symlink_to(mirror,target_is_directory=True)
    flags=list(common); flags[flags.index('--repo-root')+1]=str(link)
    run('symlinked_source_ancestor',flags,expect='symlink-free input ancestry')
    report={'status':'PASS_ALL_SUBPROCESS_CORRUPTIONS_REJECTED',
        'count':len(records),'normal_and_optimized':True,
        'formula_controls':len(faults)*2,
        'source_report_table_recipe_controls':len(records)-len(faults)*2,
        'all_have_nonzero_exit_and_absent_success_artifact':True,
        'response_sha256':sha(args.response),'caps_sha256':cap_sha,
        'freeze_sha256':freeze_sha,'adapter_sha256':sha(adapter),
        'runner_sha256':sha(__file__),'records':records}
    dump(out/'CONTROLS.json',report)
    print(json.dumps({'status':report['status'],'count':len(records)}))

if __name__=='__main__':
    main()
