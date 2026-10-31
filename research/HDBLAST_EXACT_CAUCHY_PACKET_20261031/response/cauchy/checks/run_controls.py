#!/usr/bin/env python3
"""Fresh normal/-O response corruptions for the separate Cauchy candidate."""
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(ok,message):
    if not ok: raise RuntimeError(message)
def dump(p,o):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(o,indent=2,sort_keys=True)+'\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    for role in ('same-response','absolute','velocity','independent-velocity'):
        ap.add_argument('--'+role+'-dir',type=Path,required=True)
    ap.add_argument('--response',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args(); out=args.output_dir
    need(not out.exists(),'fresh controls directory');out.mkdir(parents=True)
    base=Path(__file__).resolve().parent; adapter=base/'response_cauchy.py';freeze=base/'FREEZE.json'
    f=json.loads(freeze.read_text()); response=json.loads(args.response.read_text())
    roots={'repo_root':args.repo_root,**{role:getattr(args,role+'_dir') for role in ('same_response','absolute','velocity','independent_velocity')}}
    def flags_for(mapping):
        flags=['--repo-root',str(mapping['repo_root']),'--freeze',str(freeze),'--freeze-sha256',sha(freeze)]
        for role in ('same_response','absolute','velocity','independent_velocity'):
            flags += ['--'+role.replace('_','-')+'-dir',str(mapping[role])]
        return flags
    common=flags_for(roots);records=[]
    def run(name,flags=None,script=adapter,expect=None):
        for mode,opts in [('normal',[]),('optimized',['-O'])]:
            p=out/(name+'_'+mode+'.json')
            r=subprocess.run([sys.executable]+opts+[str(script)]+(common if flags is None else flags)+['--output',str(p)],capture_output=True,text=True)
            need(r.returncode!=0 and not p.exists(),'corruption accepted '+name+' '+mode)
            need(expect is None or expect in r.stderr,'unexpected rejection '+name+' '+r.stderr[-700:])
            portable_stderr=r.stderr.replace(str(out),'<controls>')
            for role,root in roots.items():
                portable_stderr=portable_stderr.replace(str(root),'<'+role+'>')
            portable_stderr=portable_stderr.replace(str(base.parent),'<cauchy-response>')
            records.append({'control':name,'mode':mode,'returncode':r.returncode,
                'success_artifact_absent':True,'expected_rejection':expect,
                'last_error':portable_stderr.strip().splitlines()[-1],
                'stderr_sha256':hashlib.sha256(portable_stderr.encode()).hexdigest()})
    faults=('reuse_old_velocity_absolute','reuse_old_velocity_errors','omit_f00','omit_r00','omit_full_Robin_source','omit_original_Fx','omit_derivative_commutator','lose_r1_initial_cancellation','omit_feedback_inverse')
    for fault in faults:run(fault,common+['--fault',fault],expect='CORRUPTION_REJECTED:'+fault)
    for label,mutate,expected in [
        ('report_rr',lambda r:r['ledger']['second_coordinate_error_H_caps'].__setitem__('rr',str(Q(r['ledger']['second_coordinate_error_H_caps']['rr'])+1)),'verified Cauchy report exact nodes'),
        ('report_C1',lambda r:r['ledger'].__setitem__('complete_C1_error_X_cap','0'),'verified Cauchy report exact nodes'),
        ('report_old_candidate_role',lambda r:r.__setitem__('same_function_as_W1_exact_packet',True),'verified distinct candidate role'),
        ('report_old_velocity_error',lambda r:r['initial_error_H_caps'].__setitem__('v0_r','0'),'verified fresh Cauchy errors')]:
        r=deepcopy(response);mutate(r);p=out/'fixtures'/(label+'.json');dump(p,r)
        run(label,common+['--verify-report',str(p)],expect=expected)
    for label,mutate,expected in [
        ('report_false_absolute_velocity',lambda r:r['candidate_initial_absolute_H_caps'].__setitem__('v_r','0'),'verified Cauchy exact absolute caps'),
        ('report_false_original_Fx_constant',lambda r:r['constants'].__setitem__('Fx','0'),'verified Cauchy report unchanged constants'),
        ('report_false_field_role',lambda r:r['field_roles'].__setitem__('v/v_r','old P1 absolute caps'),'verified Cauchy complete canonical report identity'),
        ('report_false_scope_flag',lambda r:r.__setitem__('old_P1_velocity_caps_not_used_for_new_candidate',False),'verified Cauchy complete canonical report identity')]:
        r=deepcopy(response);mutate(r);p=out/'fixtures'/(label+'.json');dump(p,r)
        run(label,common+['--verify-report',str(p)],expect=expected)
    mirror=out/'fixtures/inputs'; mapping={role:mirror/role for role in roots}
    for role,pins in f['input_pins'].items():
        for rel in pins:
            p=mapping[role]/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(roots[role]/rel,p)
    # Sources are copied from the complete exact report pin closure, not an
    # entire repository. Never mutate the original research tree.
    original_pins={k:v for k,v in response['selected_input_sha256_before'].items() if k.startswith('research/')}
    for rel in original_pins:
        p=mapping['repo_root']/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(args.repo_root/rel,p)
    mapped=flags_for(mapping)
    fixtures=[
      ('native_fresh_velocity_error','velocity','INITIAL_NORMAL.json',lambda r:r['campaigns'][0]['initial_error_H_caps'].__setitem__('v0_r','0')),
      ('native_stale_velocity_absolute','velocity','INITIAL_NORMAL.json',lambda r:r['campaigns'][0]['candidate_initial_absolute_H_caps'].__setitem__('v_r','0')),
      ('native_velocity_trace','velocity','INITIAL_NORMAL.json',lambda r:r['campaigns'][0]['endpoint_componentwise_affine_matching'][0].__setitem__('exterior_minus_interior_exact_jet',[0,1])),
      ('native_skipped_velocity_leaf','velocity','INITIAL_NORMAL.json',lambda r:r['campaigns'][0]['leaves'].pop()),
      ('independent_local_gate','independent_velocity','INITIAL_NORMAL.json',lambda r:r['local'].__setitem__('gates_pass',False)),
      ('independent_outer_C1_trace','independent_velocity','REGULARITY_FINAL_NORMAL.json',lambda r:r['complete_component_velocity_affine_extension_jet_differences_through1'].__setitem__('L:0',['0','1'])),
      ('absolute_skipped_refinement_leaf','absolute','REFINEMENT.json',lambda r:r['leaves'].pop()),
      ('original_Fx','repo_root','research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json',lambda r:r.__setitem__('genuine_scaled_F_x_H_cap','0')),
      ('original_published_baseline','repo_root','research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/RESPONSE.json',lambda r:r['ledger']['second_coordinate_error_H_caps'].__setitem__('rr','0')),
    ]
    for label,role,rel,mutate in fixtures:
        p=mapping[role]/rel;original=p.read_bytes();r=json.loads(original);mutate(r);dump(p,r)
        run(label,mapped,expect='selected source/report/table pin:');p.write_bytes(original)
    helper=mapping['same_response']/'checks/response_adapter.py'; original=helper.read_bytes();helper.write_text(helper.read_text().replace("k['Fx']",'Q(0)',1))
    run('changed_imported_Fraction_recipe',mapped,expect='selected source/report/table pin: same_response/checks/response_adapter.py');helper.write_bytes(original)
    changed=out/'fixtures/changed_cauchy_recipe.py';changed.write_text(adapter.read_text()+'\n# tampered recipe\n')
    run('changed_Cauchy_recipe_source',script=changed,expect='frozen Cauchy recipe source')
    altered=out/'fixtures/changed_FREEZE.json';data=deepcopy(f);data['required_gates'].pop();dump(altered,data)
    flags=list(common);flags[flags.index('--freeze')+1]=str(altered)
    run('changed_frozen_gate_recipe',flags,expect='frozen Cauchy response inputs/recipe')
    link=out/'fixtures/velocity_link';link.symlink_to(mapping['velocity'],target_is_directory=True)
    linked=dict(mapping,velocity=link)
    run('symlinked_velocity_source',flags_for(linked),expect='regular symlink-free source:')
    report={'status':'PASS_SEPARATE_CAUCHY_RESPONSE_CORRUPTIONS_REJECTED','count':len(records),
        'formula_and_stale_velocity_cases':len(faults),'normal_and_optimized':True,
        'all_nonzero_exit_and_success_artifact_absent':True,'response_sha256':sha(args.response),
        'adapter_sha256':sha(adapter),'freeze_sha256':sha(freeze),'runner_sha256':sha(__file__),
        'records':records}
    dump(out/'CONTROLS.json',report)
    print(json.dumps({'status':report['status'],'count':len(records)}))

if __name__=='__main__':main()
