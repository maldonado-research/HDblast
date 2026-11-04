#!/usr/bin/env python3
"""One portable replay: three canonical mode pairs and eight real rejections."""
from hashlib import sha256
from pathlib import Path
import argparse
import json
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
SCRIPT_PINS={
 'recompute.py':'747ba7cbb2cab70911e0db1d23beef416b53611b917ee214f4a5f1d1a0f6140b',
 'compare_candidate.py':'3393d39a0c6d4f762c2cc629c6f3694e397fd062947d62047d38557e282d00aa',
 'witness.py':'606ffd50571c4f902e3c88db9a54928a76f18c91dd18b20cc11be968d831e0dc'}
LOCAL_PINS={
 'inputs/FROZEN_SCIENTIFIC_BASELINE.json':'e89d21d78395601cbf2af98b1fc127c10d1c6768a5cca489353d8774bcfd405a',
 'inputs/PORTABLE_PRODUCER_PROVENANCE.json':'a220391f72e8eaf98cbe15cd319f320556fc24ff8d70bbb255669b50ccc0f6c3'}
SOURCE_PINS={
 'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/inputs/ACTUAL_ANNULUS_STATIC_JETS.json':'3aa2c8d2e9ad652b72a34b4d8a6e6a6f14f72c32b8dcd1fbdc67556093265c22',
 'research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json':'1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3',
 'research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/RAW_WITNESS.json':'412c4840503166a1a1893b8f16b9a3be35a8af98643c461dd55d820c1dc79d76'}
RESPONSE_REL='research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json'
RAW_REL='research/HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/numeric/references/RAW_WITNESS.json'
CANDIDATE_SHA='16efc39b531ea3aea14aa6157f3f804df5d335a34294aa7cae7fb3fb9a6d00d5'
EXPECTED_REJECTION_REASONS={
 'missing_required_output':'invalid CLI arguments',
 'relative_output':'absolute normalized non-symlink output required',
 'existing_output':'output already exists',
 'tampered_static_byte_pin':'source hash: static',
 'witness_symlink_output':'output_already_exists',
 'omit_scalar_coefficient_derivative':'conditioned_component_caps_per_R/abs_k_r',
 'extra_half_in_proper_volume':'observable_operator_caps_per_N/volume',
 'drop_Gram_angular_weight':'Gram_norm_caps_per_R/D0'}

class Rejected(Exception): pass

def need(ok,why):
    if not ok: raise Rejected(why)

class Parser(argparse.ArgumentParser):
    def error(self,message): raise Rejected('invalid CLI arguments')

def absolute(path,label):
    need(path.is_absolute() and str(path)==str(path.resolve()),'absolute normalized non-symlink '+label+' required')

def pinned(path,digest,label):
    data=path.read_bytes()
    need(sha256(data).hexdigest()==digest,'byte pin: '+label)
    return data

def full_json(payload):
    need(payload.endswith('\n') and payload.count('\n')==1 and '\r' not in payload,'stdout must be exactly one complete JSON line')
    value=json.loads(payload)
    need(isinstance(value,dict),'stdout must be one complete JSON object')
    return value

def child(script,flags,root,candidate,output,extra=()):
    command=[sys.executable,*flags,str(HERE/script),'--repo-root',str(root),
             '--candidate',str(candidate),'--output',str(output),*extra]
    result=subprocess.run(command,capture_output=True,text=True,check=False)
    need(result.stderr=='','child stderr must be empty')
    parsed=full_json(result.stdout)
    return result,parsed

def main():
    need(sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode,'required interpreter flags: -I -S -B')
    parser=Parser(add_help=False,allow_abbrev=False)
    for flag in ['--repo-root','--candidate','--output']: parser.add_argument(flag,type=Path,required=True)
    argv=sys.argv[1:]
    for flag in ['--repo-root','--candidate','--output']:
        need(sum(x==flag or x.startswith(flag+'=') for x in argv)<=1,'duplicate CLI argument')
    args=parser.parse_args(argv)
    for path,label in [(args.repo_root,'repo-root'),(args.candidate,'candidate'),(args.output,'output')]: absolute(path,label)
    need(args.repo_root.is_dir() and args.candidate.is_file(),'source paths must exist')
    need(args.output.parent.is_dir() and not args.output.exists() and not args.output.is_symlink(),'fresh output with existing parent required')
    artifact_dir=args.output.with_suffix('.artifacts')
    need(not artifact_dir.exists() and not artifact_dir.is_symlink(),'fresh canonical artifact directory required')
    # No project imports. Every local/source byte pin precedes child execution.
    for name,digest in {**SCRIPT_PINS,**LOCAL_PINS}.items(): pinned(HERE/name,digest,name)
    for name,digest in SOURCE_PINS.items(): pinned(args.repo_root/name,digest,name)
    pinned(args.candidate,CANDIDATE_SHA,'portable producer candidate')
    artifact_dir.mkdir()
    canonical=[]
    reports={}
    try:
        for role,script,candidate,extra in [
          ('BOUND','recompute.py',args.repo_root/RESPONSE_REL,('--controls',)),
          ('COMPARISON','compare_candidate.py',args.candidate,()),
          ('WITNESS','witness.py',args.repo_root/RAW_REL,())]:
            payloads=[]
            for mode,flags in [('NORMAL',['-I','-S','-B']),('OPTIMIZED',['-I','-S','-B','-O'])]:
                output=artifact_dir/(role+'_'+mode+'.json')
                result,stdout=child(script,flags,args.repo_root,candidate,output,extra)
                need(result.returncode==0 and output.is_file(),'canonical child must succeed')
                payload=output.read_bytes();report=json.loads(payload)
                digest=sha256(payload).hexdigest()
                need(stdout.get('output_sha256',stdout.get('report_sha256'))==digest,'stdout report digest')
                need(stdout['status']==report['status'] and report['status'].startswith('PASS'),'canonical status')
                canonical.append({'artifact':output.name,'script':script,'mode':mode,
                    'interpreter_flags':flags,'sha256':digest,'exit_code':0,
                    'stdout_complete_JSON':stdout,'stderr_empty':True})
                payloads.append(payload);reports[role]=report
            need(payloads[0]==payloads[1],'normal optimized canonical bytes differ: '+role)
        need(reports['BOUND']['negative_controls']['in_process_mutant_count']==6,'six parsed-source mutants')
        need(reports['COMPARISON']['in_process_report_mutants']==9,'nine candidate mutants')
        need(reports['COMPARISON']['accepted_check_count']==30,'thirty independent grouped checks')
        need(reports['COMPARISON']['full_scientific_projection_matches_frozen_baseline'] is True,'full baseline comparison')
        controls=[]
        with tempfile.TemporaryDirectory(prefix='independent-controls-',dir=args.output.parent) as folder:
            scratch=Path(folder)
            response=args.repo_root/RESPONSE_REL
            def reject(label,script,root,candidate,output,extra=(),existing=None):
                before=existing.read_bytes() if existing else None
                result,stdout=child(script,['-I','-S','-B'],root,candidate,output,extra)
                need(result.returncode==2 and stdout=={'status':'REJECTED','reason':EXPECTED_REJECTION_REASONS[label]},'exact real CLI rejection: '+label)
                if existing:
                    need(existing.read_bytes()==before,'preexisting output changed')
                else:
                    need(not output.exists() and not output.is_symlink(),'rejected CLI created artifact')
                controls.append({'label':label,'script':script,'interpreter_flags':['-I','-S','-B'],
                    'exit_code':2,'stdout_complete_JSON':stdout,'stderr_empty':True,'no_new_output_artifact':True})
            # Missing-output is an actual CLI invocation, with the same isolation flags.
            command=[sys.executable,'-I','-S','-B',str(HERE/'recompute.py'),'--repo-root',str(args.repo_root),'--candidate',str(response)]
            result=subprocess.run(command,capture_output=True,text=True,check=False)
            stdout=full_json(result.stdout)
            need(result.returncode==2 and result.stderr=='' and stdout=={'status':'REJECTED','reason':EXPECTED_REJECTION_REASONS['missing_required_output']},'exact missing output rejection')
            controls.append({'label':'missing_required_output','script':'recompute.py','interpreter_flags':['-I','-S','-B'],
                'exit_code':2,'stdout_complete_JSON':stdout,'stderr_empty':True,'no_new_output_artifact':True})
            reject('relative_output','recompute.py',args.repo_root,response,Path('relative.json'))
            existing=scratch/'existing.json';existing.write_bytes(b'fixed sentinel\n')
            reject('existing_output','recompute.py',args.repo_root,response,existing,existing=existing)
            badroot=scratch/'bad-root'
            for relative in [next(iter(SOURCE_PINS)),RESPONSE_REL]:
                destination=badroot/relative;destination.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(args.repo_root/relative,destination)
            with (badroot/next(iter(SOURCE_PINS))).open('ab') as stream:stream.write(b' ')
            reject('tampered_static_byte_pin','recompute.py',badroot,badroot/RESPONSE_REL,scratch/'bad-pin.json')
            target=scratch/'witness-target.json';link=scratch/'witness-link.json';link.symlink_to(target)
            result,stdout=child('witness.py',['-I','-S','-B'],args.repo_root,args.repo_root/RAW_REL,link)
            need(result.returncode==2 and stdout=={'status':'REJECTED','reason':EXPECTED_REJECTION_REASONS['witness_symlink_output']} and not target.exists() and link.is_symlink(),'exact witness symlink rejection')
            controls.append({'label':'witness_symlink_output','script':'witness.py','interpreter_flags':['-I','-S','-B'],
                'exit_code':2,'stdout_complete_JSON':stdout,'stderr_empty':True,'no_new_output_artifact':True,'preexisting_symlink_preserved':True})
            candidate=json.loads(args.candidate.read_bytes())
            for label,section,key in [
              ('omit_scalar_coefficient_derivative','conditioned_component_caps_per_R','abs_k_r'),
              ('extra_half_in_proper_volume','observable_operator_caps_per_N','volume'),
              ('drop_Gram_angular_weight','Gram_norm_caps_per_R','D0')]:
                mutant=json.loads(json.dumps(candidate));mutant[section][key]['exact']='0'
                path=scratch/(label+'.json');path.write_text(json.dumps(mutant))
                reject(label,'compare_candidate.py',args.repo_root,path,scratch/(label+'-output.json'))
        need(len(controls)==8,'eight actual CLI negatives')
        report={'schema':'HDBLAST_INDEPENDENT_PORTABLE_VALIDATION_V1','status':'PASS_INDEPENDENT_PORTABLE_VALIDATION',
          'canonical_artifact_directory':artifact_dir.name,'canonical_children':canonical,
          'canonical_pairs_byte_identical':True,'canonical_child_processes':6,
          'canonical_normal_processes':3,'canonical_optimized_processes':3,
          'source_mutants_unique_cases':6,'source_mutant_evaluations_across_canonical_modes':12,
          'candidate_report_mutants_unique_cases':9,'candidate_report_mutant_evaluations_across_canonical_modes':18,
          'actual_CLI_negative_control_count':8,'actual_CLI_negative_controls':controls,'total_child_processes':14,
          'stdout_exactly_one_complete_JSON_line_checked':True,'all_eight_exact_rejection_witnesses_checked':True,
          'scope':{'root_old_source_closure_replays':0,'new_oracle_queries':0,'new_bump_queries':0,'new_PDE_solves':0,
            'fresh_analytic_target_bump_replay_executions':2,'new_directions':0,
            'count_note':'Unique in-process mutant cases, their evaluations, canonical child processes and actual CLI negative processes are separate counts.'}}
        payload=(json.dumps(report,indent=2,sort_keys=True)+'\n').encode()
        with args.output.open('xb') as stream:stream.write(payload)
        print(json.dumps({'status':report['status'],'output_sha256':sha256(payload).hexdigest(),'actual_CLI_negative_control_count':8},sort_keys=True))
    except Exception:
        shutil.rmtree(artifact_dir)
        raise

if __name__=='__main__':
    try:main()
    except Exception as error:
        print(json.dumps({'status':'REJECTED','reason':str(error)},sort_keys=True))
        sys.exit(2)
