#!/usr/bin/env python3
"""Closed-capsule exact Cauchy checkpoint, normal/-O and real controls."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent

def need(ok,why):
    if not ok:raise RuntimeError(why)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):return json.loads(path.read_bytes())

def authenticate(cap):
    raw=(cap/'MANIFEST.json').read_bytes();seen={'MANIFEST.json'}
    for row in json.loads(raw)['files']:
        rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'safe unique capsule member')
        seen.add(str(rel));path=cap/rel
        need(path.is_file() and not path.is_symlink() and all(not p.is_symlink() for p in path.parents),'regular capsule source')
        data=path.read_bytes();need(len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],'capsule member '+str(rel))
    need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen) for p in cap.rglob('*')),'closed capsule file set')
    return hashlib.sha256(raw).hexdigest()

def resolve(value,bindings):return value.format_map(bindings)

def control_check(document,contract):
    rows=document[contract['rows']];need(isinstance(rows,list) and len(rows)==contract['count'],'actual complete control row count')
    if contract.get('count_field'):need(document[contract['count_field']]==len(rows),'reported control count')
    for row in rows:
        if contract.get('returncode'):need(isinstance(row[contract['returncode']],int) and row[contract['returncode']]!=0,'nonzero corrupted subprocess')
        for key in contract.get('true',[]):need(row[key] is True,'actual rejection flag '+key)
        for key in contract.get('nonempty',[]):need(bool(row[key]),'explicit rejection witness '+key)
    paired={}
    for row in rows:
        key=row[contract['case_field']];mode=row[contract['mode_field']]
        paired.setdefault(key,[]).append(mode)
    need(all(len(modes)==2 and set(modes)==set(contract['mode_values']) for modes in paired.values()),'each corrupted case executed in both modes')
    need(str(document['status']).startswith('PASS'),'control driver status')
    return len(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--output-dir',required=True,type=Path);ap.add_argument('--jobs',type=int,choices=(1,2,3),default=3)
    a=ap.parse_args();root=a.repo_root.resolve();cap=HERE.parent;out=a.output_dir.resolve()
    closure=authenticate(cap);plan=read(HERE/'REPLAY_PLAN.json')
    need(plan['schema']=='HDBLAST_EXACT_CAUCHY_ROOT_REPLAY_PLAN_V1','replay plan schema')
    need(not out.exists() and cap not in out.parents and root not in out.parents,'fresh output outside source tree and capsule')
    out.mkdir(parents=True);(out/'logs').mkdir()
    bindings={'CAP':str(cap),'SOURCE':str(root),'OUT':str(out)}
    patches=[]
    for stage in plan['stages']:
        source=cap/stage['component'];dest=out/'staged'/stage['name'];shutil.copytree(source,dest)
        bindings[stage['placeholder']]=str(dest)
        for rel in stage.get('isolate_wrappers',[]):
            path=dest/rel;before=path.read_bytes();text=before.decode();n=text.count('[sys.executable]');m=text.count('[sys.executable,str(')
            text=text.replace('[sys.executable]',"[sys.executable,'-I','-S','-B']").replace('[sys.executable,str(',"[sys.executable,'-I','-S','-B',str(")
            need(n+m>0,'wrapper isolation substitution anchor '+stage['component']+'/'+rel)
            path.write_text(text)
            patches.append({'component':stage['component'],'wrapper':rel,'changes':'subprocess interpreter flags only: -I -S -B','replacement_count':n+m,'source_sha256':hashlib.sha256(before).hexdigest(),'runtime_copy_sha256':sha(path)})
    def execute_group(group):
        canonical={};modes={};controls={};dest=out/group['name'];dest.mkdir()
        for step in group['steps']:
            if step['kind']=='canonical':
                output=Path(resolve(step['output'],bindings));output.parent.mkdir(parents=True,exist_ok=True)
            command=[sys.executable,'-I','-S','-B']+(['-O'] if step.get('optimized',False) else [])+[resolve(x,bindings) for x in step['argv']]
            env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1'
            p=subprocess.run(command,capture_output=True,timeout=step.get('timeout_seconds',1800),env=env,cwd=root)
            log=out/'logs'/(group['name']+'-'+step['name']);log.with_suffix('.stdout').write_bytes(p.stdout);log.with_suffix('.stderr').write_bytes(p.stderr)
            need(p.returncode==0,'checkpoint step '+group['name']+'/'+step['name']+' failed: '+p.stderr.decode(errors='replace')[-2500:])
            if step['kind']=='canonical':
                expected=cap/step['reference'];need(output.read_bytes()==expected.read_bytes(),'complete canonical byte identity '+step['name'])
                digest=sha(output);key=step['canonical_key'];need(key not in canonical or canonical[key]==digest,'normal/optimized exact identity '+key)
                canonical[key]=digest;modes.setdefault(key,[]).append(step.get('optimized',False))
            elif step['kind']=='controls':
                receipt=read(Path(resolve(step['receipt'],bindings)));controls[step['name']]=control_check(receipt,step['contract'])
            else:raise RuntimeError('unknown replay step kind')
        need(all(len(v)==2 and set(v)=={False,True} for v in modes.values()),'every canonical report replayed in both modes')
        print(json.dumps({'completed_group':group['name'],'canonical_reports':len(canonical),'rejected_subprocesses':sum(controls.values())}),flush=True)
        return {'group':group['name'],'canonical_report_sha256':canonical,'actual_corrupted_subprocesses_rejected':controls}
    need(1<=len(plan['groups'])<=3,'at most three independent groups')
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        jobs=[pool.submit(execute_group,group) for group in plan['groups']];results=[job.result() for job in jobs]
    need(authenticate(cap)==closure,'closed source capsule unchanged after replay')
    total=sum(sum(group['actual_corrupted_subprocesses_rejected'].values()) for group in results)
    count=sum(len(group['canonical_report_sha256']) for group in results)
    need(total==plan['expected_rejected_subprocesses'] and count==plan['expected_canonical_reports'],'complete checkpoint report/control roster')
    result={'status':'PASS_EXACT_CAUCHY_PACKET_AND_SAME_CANDIDATE_REFINEMENT','capsule_manifest_sha256':closure,
      'canonical_reports_reproduced_in_normal_and_optimized':count,'actual_corrupted_subprocesses_rejected':total,
      'groups':results,'runtime_wrapper_isolation_patches':patches,'closed_capsule_authenticated_before_and_after':True,
      'no_canonical_mathematical_source_or_report_modified':True,'scope':plan['scope']}
    (out/'VERIFICATION.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
