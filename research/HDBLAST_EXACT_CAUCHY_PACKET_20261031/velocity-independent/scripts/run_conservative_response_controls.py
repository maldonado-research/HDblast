#!/usr/bin/env python3
"""Real source omissions and altered premise controls for independent response."""
import argparse,hashlib,json,shutil,subprocess,sys
from pathlib import Path
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
BASE='research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/'
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',required=True,type=Path);ap.add_argument('--absolute-dir',required=True,type=Path);ap.add_argument('--native-response',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);ap.add_argument('--receipt',required=True,type=Path);args=ap.parse_args()
    need(not args.output_dir.exists() and not args.receipt.exists(),'FRESH_CONTROLS');args.output_dir.mkdir(parents=True)
    helper=HERE/'check_conservative_response.py';original=helper.read_text();freeze=HERE.parent/'CONSERVATIVE_RESPONSE_PORTABLE_FREEZE.json'
    replacements=[
      ('omit-f00','initial_acceleration=P*fzero','initial_acceleration=F(0)*fzero','BASELINE_EXACT_NODE:W0_minus_Bs_time_derivative_X_cap'),
      ('omit-r00','initial_acceleration=P*rzero','initial_acceleration=F(0)*rzero','BASELINE_EXACT_NODE:correction_time_derivative_X_cap'),
      ('omit-robin','Robin_free_wave_derivative=P*M*Bs','Robin_free_wave_derivative=F(0)*M*Bs','BASELINE_EXACT_NODE:W0_minus_Bs_time_derivative_X_cap'),
      ('omit-original-Fx',"original_source_derivative=k['Fx']","original_source_derivative=F(0)*k['Fx']",'BASELINE_EXACT_NODE:r0_time_derivative_H_cap'),
      ('omit-commutator',"coefficient_commutator=Qinv*T*k['comm']*E","coefficient_commutator=F(0)*Qinv*T*k['comm']*E",'BASELINE_EXACT_NODE:second_coordinate_error_H_caps'),
      ('omit-feedback-inverse',"T,P,Qinv = k['T'],k['P'],1/(1-k['q'])","T,P,Qinv = k['T'],k['P'],F(1)",'BASELINE_EXACT_NODE:complete_C1_error_X_cap'),
      ('lose-common-first-Cauchy-jet-cancellation','initial=Qinv*j, differentiated_defect','initial=Qinv*(j+rzero), differentiated_defect','BASELINE_EXACT_NODE:second_coordinate_error_H_caps'),
      ('replace-absolute-curvature-with-error',"displacement_curvature=u['u_rr'], velocity_gradient","displacement_curvature=e['e0_rr'], velocity_gradient",'BASELINE_EXACT_NODE:free_wave_second_jet_H_cap')]
    scenarios=[]
    for name,a,b,group in replacements:
        need(original.count(a)==1,'UNIQUE_REAL_SCIENTIFIC_REPLACEMENT:'+name);p=args.output_dir/(name+'.py');p.write_text(original.replace(a,b));scenarios.append((name,p,freeze,{},group))
    badfreeze=args.output_dir/'wrong-frozen-own-initial-hash.json';x=json.loads(freeze.read_bytes());x['inputs']['own_initial']['sha256']='0'*64;badfreeze.write_text(json.dumps(x,sort_keys=True)+'\n');scenarios.append(('wrong-frozen-own-initial-hash',helper,badfreeze,{},'FROZEN_INPUT:own_initial'))
    badreg=args.output_dir/'bad-completed-velocity-trace.json';x=json.loads((HERE.parent/'REGULARITY_FINAL_NORMAL.json').read_bytes());x['complete_component_velocity_affine_extension_jet_differences_through1']['K:0']=['0','1'];badreg.write_text(json.dumps(x,sort_keys=True)+'\n')
    badregfreeze=args.output_dir/'bad-trace-premise-freeze.json';x=json.loads(freeze.read_bytes());x['inputs']['own_regularity']['sha256']=sha(badreg);badregfreeze.write_text(json.dumps(x,sort_keys=True)+'\n');scenarios.append(('invalid-completed-velocity-C1-premise',helper,badregfreeze,{'own_regularity':badreg},'COMPLETED_C1_VELOCITY'))
    roles={'own_initial':HERE.parent/'INITIAL_NORMAL.json','own_regularity':HERE.parent/'REGULARITY_FINAL_NORMAL.json','own_audit':HERE.parent/'AUDIT_PORTABLE_FINAL_NORMAL.json','absolute':args.absolute_dir/'REFINEMENT.json','baseline':args.source/(BASE+'RESPONSE.json'),'baseline_initial':args.source/(BASE+'INITIAL.json'),'native_response':args.native_response}
    receipts=[]
    for name,source,fp,overrides,group in scenarios:
        for optimized in (False,True):
            output=args.output_dir/(name+('-O' if optimized else '-normal')+'-must-not-exist.json');cmd=[sys.executable]+(['-O'] if optimized else [])+[str(source),'--source',str(args.source),'--freeze',str(fp),'--output',str(output)]
            for role,path in (roles|overrides).items():cmd+=['--input',role+'='+str(path)]
            r=subprocess.run(cmd,text=True,capture_output=True);need(r.returncode!=0 and group in r.stderr and not output.exists(),'CONTROL_REJECTED_AT_REQUIRED_GROUP:'+name+':'+r.stderr)
            receipts.append({'case':name,'optimized':optimized,'exit_code':r.returncode,'expected_group':group,'observed_failure':r.stderr.strip().splitlines()[-1],'success_artifact_absent':not output.exists(),'executed_scientific_source_sha256':sha(source),'freeze_sha256':sha(fp)})
    out={'status':'PASS_EFFECTIVE_INDEPENDENT_CONSERVATIVE_RESPONSE_CONTROLS','cases':len(scenarios),'subprocess_rejections':len(receipts),'normal_optimized_same_expected_groups':True,'runner_sha256':sha(Path(__file__)),'receipts':receipts}
    args.receipt.write_text(json.dumps(out,sort_keys=True,indent=2)+'\n');print(json.dumps({k:out[k] for k in ('status','cases','subprocess_rejections')}))
if __name__=='__main__':main()
