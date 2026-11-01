#!/usr/bin/env python3
"""Eight strict root authentication mutations in both interpreter modes."""
import argparse,hashlib,json,shutil,subprocess,sys,tempfile
from pathlib import Path
HERE=Path(__file__).resolve().parent
CAP=HERE.parent
class Failure(Exception):pass
def need(ok,witness):
    if not ok:raise Failure(witness)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def manifest_rows(root):
    rows=read(root/'MANIFEST.json')['files']
    return [{'path':k,**v} for k,v in rows.items()] if isinstance(rows,dict) else rows
def repin_member(root,rel):
    document=read(root/'MANIFEST.json');rows=document['files'];p=root/rel
    if isinstance(rows,dict):rows[rel]={'bytes':p.stat().st_size,'sha256':sha(p)}
    else:
        matches=[row for row in rows if row['path']==rel];need(len(matches)==1,'ONE_MANIFEST_ROW');matches[0].update(bytes=p.stat().st_size,sha256=sha(p))
    (root/'MANIFEST.json').write_bytes(canonical(document));return sha(root/'MANIFEST.json')
def invoke(runner,repo,out,pin,optimized,witness,positive):
    command=[sys.executable,'-I','-S','-B']+(['-O'] if optimized else [])+[str(runner),'--repo-root',str(repo),'--output-dir',str(out),'--manifest-sha256',pin,'--authenticate-only']
    proc=subprocess.run(command,capture_output=True,text=True,timeout=120);lines=proc.stdout.splitlines()
    need(len(lines)==1 and not proc.stderr,'EXACT_SINGLE_FINAL_JSON_WITNESS');final=json.loads(lines[0])
    need(proc.returncode==(0 if positive else 2),'EXACT_CONTROL_RETURN_CODE')
    if positive:need(final.get('status')=='PASS' and final.get('witness')==witness,'EXACT_POSITIVE_FINAL_WITNESS')
    else:need(final=={'status':'FAIL','witness':witness},'EXACT_NEGATIVE_FINAL_WITNESS')
    need(not out.exists(),'AUTH_ONLY_CREATED_NO_OUTPUT')
    return {'subprocess_exit':proc.returncode,'final_witness':final,'stderr_empty':True,'output_directory_absent':True,'authentication_only':True}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    repo=a.repo_root.resolve(strict=True);output=a.output.absolute();need(not output.exists() and all(not x.is_symlink() for x in [output,*output.parents]),'FRESH_OUTPUT');output=output.resolve()
    need(all(x!=output and x not in output.parents for x in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    own=sha(Path(__file__));manifestpin=sha(CAP/'MANIFEST.json');plan=read(HERE/'REPLAY_PLAN.json')
    oldname=sorted(plan['source_capsules'])[0];oldroot=repo/'research'/oldname;oldpin=sha(oldroot/'MANIFEST.json')
    need(oldpin==plan['source_capsules'][oldname],'BASELINE_OLD_MANIFEST_PIN');oldmember=sorted(row['path'] for row in manifest_rows(oldroot))[0]
    newmember=sorted(x for x in plan['critical_new_files'] if x.endswith('.py') and not x.startswith('checks/'))[0]
    selectedold=sorted(plan['selected_source_files'])[0];positives=[];rejections=[]
    with tempfile.TemporaryDirectory(prefix='prepared-shell-volume-root-controls-') as temporary:
        scratch=Path(temporary)
        relocated=scratch/'relocated-component';shutil.copytree(CAP,relocated)
        for label,runner in [('original',HERE/'run_checkpoint.py'),('relocated',relocated/'checks/run_checkpoint.py')]:
            for optimized in [False,True]:
                result=invoke(runner,repo,scratch/(label+str(optimized)),manifestpin,optimized,'CLOSED_SOURCE_AUTHENTICATION',True)
                positives.append({'location':label,'optimized':optimized,**result})
        newcases=[('new_member_bytes','ROOT_MEMBER_BYTES'),('new_unlisted_member','ROOT_CLOSED_INVENTORY'),
                  ('new_manifest_bytes','ROOT_MANIFEST_PIN'),('selected_new_source_bytes','SELECTED_NEW_SOURCE'),
                  ('selected_old_source_pin','SELECTED_OLD_SOURCE')]
        for case,witness in newcases:
            root=scratch/case/'component';shutil.copytree(CAP,root);pin=manifestpin
            if case in ['new_member_bytes','selected_new_source_bytes']:
                q=root/newmember;q.write_bytes(q.read_bytes()+b'\nraise RuntimeError("ROOT_CONTROL_SENTINEL_MUST_NOT_RUN")\n')
                if case=='selected_new_source_bytes':pin=repin_member(root,newmember)
            elif case=='new_unlisted_member':(root/'ROOT_CONTROL_UNLISTED.txt').write_text('Deliberate closed-inventory mutation.\n')
            elif case=='new_manifest_bytes':(root/'MANIFEST.json').write_bytes((root/'MANIFEST.json').read_bytes()+b' ')
            elif case=='selected_old_source_pin':
                q=root/'checks/REPLAY_PLAN.json';document=read(q);document['selected_source_files'][selectedold]['sha256']='0'*64;q.write_bytes(canonical(document));pin=repin_member(root,'checks/REPLAY_PLAN.json')
            for optimized in [False,True]:
                result=invoke(root/'checks/run_checkpoint.py',repo,scratch/(case+str(optimized)),pin,optimized,witness,False)
                rejections.append({'case':case,'optimized':optimized,'expected_final_witness':{'status':'FAIL','witness':witness},**result})
        for case,witness in [('old_manifest_bytes','SOURCE_MANIFEST_PIN'),('old_member_bytes','SOURCE_MEMBER_BYTES'),('old_unlisted_member','SOURCE_CLOSED_INVENTORY')]:
            copiedrepo=scratch/case/'source';copied=copiedrepo/'research'/oldname;copied.parent.mkdir(parents=True);shutil.copytree(oldroot,copied)
            if case=='old_manifest_bytes':q=copied/'MANIFEST.json';q.write_bytes(q.read_bytes()+b' ')
            elif case=='old_member_bytes':q=copied/oldmember;q.write_bytes(q.read_bytes()+b' ')
            else:(copied/'ROOT_CONTROL_UNLISTED.txt').write_text('Deliberate old closed-inventory mutation.\n')
            for optimized in [False,True]:
                result=invoke(HERE/'run_checkpoint.py',copiedrepo,scratch/(case+str(optimized)),manifestpin,optimized,witness,False)
                rejections.append({'case':case,'optimized':optimized,'expected_final_witness':{'status':'FAIL','witness':witness},**result})
    need(sha(Path(__file__))==own and sha(CAP/'MANIFEST.json')==manifestpin and sha(oldroot/'MANIFEST.json')==oldpin,'ROOT_CONTROL_INPUTS_UNCHANGED')
    result={'schema':'HDBLAST_PREPARED_SHELL_VOLUME_ROOT_STRICT_AUTH_CONTROLS_V1','status':'PASS_STRICT_ROOT_AUTHENTICATION_CONTROLS',
        'controls_source_sha256':own,'capsule_manifest_sha256':manifestpin,'distinct_negative_cases':8,'actual_rejected_subprocesses':len(rejections),
        'positive_subprocesses':len(positives),'positive_relocation_separate_from_negatives':True,
        'exact_final_json_witnesses_required':True,'all_failure_exit_codes':2,'authentication_only_scientific_imports':0,
        'new_actual_oracle_calls':0,'positive_checks':positives,'rejections':rejections,
        'scope':'Root source authentication only. Neither new scientific computations nor inherited proofs execute in these probes.'}
    with output.open('xb') as handle:handle.write(canonical(result))
    print(json.dumps({'status':'PASS','witness':'STRICT_ROOT_AUTHENTICATION_CONTROLS','rejected_subprocesses':len(rejections)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
