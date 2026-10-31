#!/usr/bin/env python3
"""Persist actual corrupted inputs and require normal/-O audit rejection."""
import argparse, copy, hashlib, json, shutil, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=None
PRODUCER=None
OUT=HERE.parent/'controls'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    global ROOT,PRODUCER,OUT
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--producer-root',type=Path,required=True);ap.add_argument('--output-dir',type=Path,default=OUT);ap.add_argument('--receipt',type=Path,default=HERE.parent/'CONTROLS.json');args=ap.parse_args()
    ROOT,PRODUCER,OUT=args.source,args.producer_root,args.output_dir
    need(not OUT.exists(),'FRESH_CONTROL_DIRECTORY');OUT.mkdir()
    own=json.loads((HERE.parent/'INITIAL_NORMAL.json').read_bytes())
    native=json.loads((PRODUCER/'INITIAL_NORMAL.json').read_bytes())
    scenarios=[]
    def add(name,which,fn,group):
        x=copy.deepcopy(own if which=='own' else native);fn(x)
        p=OUT/(name+'.json');p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
        scenarios.append((name,p if which=='own' else HERE.parent/'INITIAL_NORMAL.json',p if which=='native' else PRODUCER/'INITIAL_NORMAL.json',PRODUCER,group))
    add('wrong-kernel-half-mass','own',lambda x:x['local'].__setitem__('kernel_half_mass','1'),'OWN_KERNEL_EXACT_HALF')
    add('lost-all23-mean-component','own',lambda x:x['local']['rows'].pop(0),'OWN_LOCAL_ALL23')
    add('nonzero-shared-packet-channel','own',lambda x:x['local']['rows'][0]['paired_channels'].__setitem__(1,['1','1']),'OWN_ZERO_Q_E_CHANNEL')
    add('unit-adjacent-secant-weight','own',lambda x:x['local']['rows'][3].__setitem__('old_smoothed_v_r','1'),'OWN_EXACT_ADJACENT_SECANTS')
    add('wrong-old-error-sign','own',lambda x:x['local']['rows'][3].__setitem__('actual_minus_old_smoothed_v_r',['0','0']),'OWN_SIGNED_OLD_ERROR')
    add('zero-target-denominator','own',lambda x:x['local'].__setitem__('actual_v_r_H_interval',['0','0']),'OWN_LOCAL_NORMS_RECONSTRUCT')
    add('dropped-background-remainder-channel','own',lambda x:x['coverage']['leaves'][-1]['channel_component_caps'][2].__setitem__(1,['0']*23),'OWN_CHANNEL_SUPREMA')
    add('coverage-gap','own',lambda x:x['coverage']['leaves'][2]['r'].__setitem__(0,'67/64'),'OWN_EXACT_COVER')
    add('unbounded-depth','own',lambda x:x['coverage']['leaves'][0].__setitem__('depth',17),'OWN_BUDGET')
    add('native-zero-absolute-cap','native',lambda x:x['campaigns'][-1]['leaves'][-1]['absolute_H_caps'].__setitem__(0,'0'),'NATIVE_ABSOLUTE_LEAF_H')
    add('native-zero-error-cap','native',lambda x:x['campaigns'][-1]['initial_error_H_caps'].__setitem__('v0_r','0'),'NATIVE_ERROR_AGGREGATION')
    add('native-missing-affine-jet','native',lambda x:x['campaigns'][-1]['endpoint_componentwise_affine_matching'][0].__setitem__('exterior_minus_interior_exact_jet',[0,1]),'NATIVE_COMPLETED_COMPONENT_TRACE_MATCHING')
    corrupt_producer=OUT/'changed-producer';corrupt_producer.mkdir()
    for name in ['producer.py','FREEZE.json']:shutil.copyfile(PRODUCER/name,corrupt_producer/name)
    p=corrupt_producer/'producer.py';p.write_text(p.read_text().replace('old=[(a+b)/2 for a,b in zip(left,right)]','old=[(a+b) for a,b in zip(left,right)]'))
    need(sha(p)!=sha(PRODUCER/'producer.py'),'SOURCE_CORRUPTION_EFFECTIVE')
    scenarios.append(('actual-source-half-mass-corruption',HERE.parent/'INITIAL_NORMAL.json',PRODUCER/'INITIAL_NORMAL.json',corrupt_producer,'NATIVE_PRODUCER_PIN:producer.py'))
    receipts=[]
    for name,ownpath,nativepath,prod,group in scenarios:
        for optimized in (False,True):
            output=OUT/(name+('-O' if optimized else '-normal')+'-must-not-exist.json')
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_velocity_reports.py'),'--source',str(ROOT),'--own',str(ownpath),'--native',str(nativepath),'--producer-root',str(prod),'--out',str(output)]
            r=subprocess.run(cmd,text=True,capture_output=True)
            need(r.returncode!=0 and group in r.stderr and not output.exists(),'CONTROL_NOT_REJECTED:'+name+':'+r.stderr)
            receipts.append({'case':name,'optimized':optimized,'exit_code':r.returncode,'expected_group':group,'observed_failure':r.stderr.strip().splitlines()[-1],'success_artifact_absent':not output.exists(),'own_input_sha256':sha(ownpath),'native_input_sha256':sha(nativepath)})
    # A real modified nominal denominator source must be rejected before use.
    regroot=OUT/'changed-regularity-source';regpins=json.loads((HERE.parent/'regularity-source-pins.json').read_bytes())
    for rel in regpins:
        p=regroot/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,p)
    rel='research/HDBLAST_INPUT_PRECISION_20261023/initial/RECENTERED_CANDIDATE_NORMAL.json';p=regroot/rel;x=json.loads(p.read_bytes());x['representation_background_coefficients_ascending']['rho'][0]='0';p.write_text(json.dumps(x,sort_keys=True)+'\n')
    for optimized in (False,True):
        output=OUT/('source-rho-'+str(optimized)+'-must-not-exist.json')
        cmd=[sys.executable]+(['-O'] if optimized else [])+[str(HERE/'check_velocity_regularity.py'),'--source-root',str(regroot),'--producer-root',str(PRODUCER),'--pins',str(HERE.parent/'regularity-source-pins.json'),'--output',str(output)]
        r=subprocess.run(cmd,text=True,capture_output=True);group='SOURCE_HASH:'+rel
        need(r.returncode!=0 and group in r.stderr and not output.exists(),'SOURCE_CONTROL_NOT_REJECTED')
        receipts.append({'case':'actual-source-nominal-rho-corruption','optimized':optimized,'exit_code':r.returncode,'expected_group':group,'observed_failure':r.stderr.strip().splitlines()[-1],'success_artifact_absent':not output.exists(),'corrupted_source_sha256':sha(p)})
    result={'status':'PASS_EFFECTIVE_INDEPENDENT_VELOCITY_CONTROLS','cases':len(receipts)//2,'subprocess_rejections':len(receipts),'normal_optimized_same_expected_groups':True,'receipts':receipts,'runner_sha256':sha(Path(__file__))}
    args.receipt.parent.mkdir(parents=True,exist_ok=True);args.receipt.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({k:result[k] for k in ('status','cases','subprocess_rejections')}))
if __name__=='__main__':main()
