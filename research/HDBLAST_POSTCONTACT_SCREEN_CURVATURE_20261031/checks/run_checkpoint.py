#!/usr/bin/env python3
"""Authenticate and replay the bounded postcontact screen checkpoint."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent

def need(ok,message):
    if not ok:raise RuntimeError(message)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_bytes())
def regular(path):
    return path.is_file() and not path.is_symlink() and all(not q.is_symlink() for q in path.absolute().parents)
def authenticate(root,pin=None):
    need(root.is_dir() and not root.is_symlink(),'CAPSULE_ROOT')
    manifest=root/'MANIFEST.json';need(regular(manifest),'MANIFEST_TYPE')
    digest=sha(manifest)
    if pin is not None:need(digest==pin,'MANIFEST_PIN:'+root.name)
    seen={'MANIFEST.json'}
    for row in read(manifest)['files']:
        rel=Path(row['path']);need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'SAFE_UNIQUE_MEMBER')
        seen.add(str(rel));path=root/rel
        need(regular(path) and path.stat().st_size==row['bytes'] and sha(path)==row['sha256'],'CLOSED_MEMBER:'+root.name+'/'+str(rel))
    for path in root.rglob('*'):
        need(not path.is_symlink(),'CLOSED_SYMLINK')
        if path.is_file():need(str(path.relative_to(root)) in seen,'CLOSED_UNLISTED_MEMBER')
    return digest
def sources(root,plan):
    result={label:authenticate(root/'research'/label,pin) for label,pin in plan['source_capsules'].items()}
    for rel,pin in plan['selected_source_files'].items():
        path=root/rel;need(regular(path) and sha(path)==pin,'SELECTED_SOURCE:'+rel)
    return result
def resolve(value,bindings):return value.format_map(bindings)
def controls(document,contract):
    need(document['status']==contract['status'],'CONTROL_STATUS')
    rows=document[contract['rows']];need(len(rows)==contract['count'],'COMPLETE_CONTROL_ROW_COUNT')
    pairs={}
    for row in rows:
        need(isinstance(row[contract['returncode']],int) and row[contract['returncode']]!=0,'NONZERO_CONTROL_EXIT')
        for key in contract.get('true',[]):need(row[key] is True,'CONTROL_TRUE:'+key)
        for key in contract.get('false',[]):need(row[key] is False,'CONTROL_FALSE:'+key)
        for key in contract.get('nonempty',[]):need(bool(row[key]),'CONTROL_WITNESS:'+key)
        pairs.setdefault(row[contract['case']],[]).append(row[contract['mode']])
    need(all(len(v)==2 and set(v)=={False,True} for v in pairs.values()),'CONTROL_BOTH_PYTHON_MODES')
    return len(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output-dir',required=True,type=Path)
    ap.add_argument('--jobs',type=int,choices=(1,2,3),default=3)
    a=ap.parse_args();repo=a.repo_root.resolve();cap=HERE.parent;out=a.output_dir.absolute()
    need(not a.output_dir.exists() and not a.output_dir.is_symlink(),'FRESH_OUTPUT_DIRECTORY')
    need(all(not p.is_symlink() for p in out.parents),'OUTPUT_ANCESTOR_SYMLINK')
    out=out.resolve();need(repo!=out and repo not in out.parents and cap!=out and cap not in out.parents,'OUTPUT_OUTSIDE_FROZEN_SOURCE')
    own=authenticate(cap);plan=read(HERE/'REPLAY_PLAN.json')
    need(plan['schema']=='HDBLAST_POSTCONTACT_SCREEN_ROOT_REPLAY_V1','PLAN_SCHEMA')
    before=sources(repo,plan)
    # No scientific or third-party import occurs before both closures are authenticated.
    env=dict(os.environ);env['PYTHONDONTWRITEBYTECODE']='1';env['HDBLAST_REPLAY_SOURCE_ROOT']=str(repo)
    version_command=[sys.executable,'-I','-B','-c','import json,sympy,mpmath; print(json.dumps({"sympy":sympy.__version__,"mpmath":mpmath.__version__},sort_keys=True))']
    version_process=subprocess.run(version_command,capture_output=True,text=True,timeout=30,env=env)
    need(version_process.returncode==0,'PINNED_SYMBOLIC_REQUIREMENTS_MISSING')
    versions=json.loads(version_process.stdout)
    need(versions==plan['requirements'],'PINNED_SYMBOLIC_REQUIREMENTS_VERSION')
    out.mkdir(parents=True);(out/'logs').mkdir()
    bindings={'CAP':str(cap),'SOURCE':str(repo),'OUT':str(out)}
    patches=[];staged_original={}
    for stage in plan['stages']:
        destination=out/'staged'/stage['name'];destination.mkdir(parents=True)
        bindings[stage['placeholder']]=str(destination)
        for rel in stage['files']:
            source=cap/stage['component']/rel;target=destination/rel;target.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(source,target);staged_original[str(target)]=sha(target)
    for stage in plan['stages']:
        destination=Path(bindings[stage['placeholder']])
        for patch in stage.get('runtime_wrapper_patches',[]):
            path=destination/patch['path'];old=path.read_bytes();text=old.decode()
            anchor=patch['anchor'];replacement=resolve(patch['replacement'],bindings)
            need(text.count(anchor)==patch['count'],'UNIQUE_RUNTIME_WRAPPER_ANCHOR:'+patch['path'])
            text=text.replace(anchor,replacement);path.write_text(text);staged_original[str(path)]=sha(path)
            patches.append({'component':stage['component'],'wrapper':patch['path'],'purpose':patch['purpose'],'original_sha256':hashlib.sha256(old).hexdigest(),'runtime_sha256':sha(path),'substitutions':patch['count']})
    def group_run(group):
        canonical={};modes={};counts={};timings={};extras={}
        for step in group['steps']:
            argv=[resolve(x,bindings) for x in step['argv']]
            flags=['-I','-B'] if step.get('requires_sympy') else ['-I','-S','-B']
            command=[sys.executable]+flags+(['-O'] if step.get('optimized') else [])+argv
            if 'output' in step:Path(resolve(step['output'],bindings)).parent.mkdir(parents=True,exist_ok=True)
            started=time.monotonic();p=subprocess.run(command,capture_output=True,timeout=step.get('timeout_seconds',240),env=env,cwd=repo)
            timings[step['name']]=round(time.monotonic()-started,6)
            log=out/'logs'/(group['name']+'-'+step['name']);log.with_suffix('.stdout').write_bytes(p.stdout);log.with_suffix('.stderr').write_bytes(p.stderr)
            need(p.returncode==0,'REPLAY_STEP:'+group['name']+'/'+step['name']+':'+p.stderr.decode(errors='replace')[-3000:])
            if step['kind']=='canonical':
                produced=Path(resolve(step['output'],bindings));expected=cap/step['reference']
                need(sha(expected)==step['expected_sha256'],'REFERENCE_EXPECTED_PIN:'+step['name'])
                need(produced.read_bytes()==expected.read_bytes(),'CANONICAL_BYTE_IDENTITY:'+step['name'])
                key=step['canonical_key'];digest=sha(produced)
                need(key not in canonical or canonical[key]==digest,'NORMAL_OPTIMIZED_IDENTITY')
                canonical[key]=digest;modes.setdefault(key,[]).append(bool(step.get('optimized')))
                d=read(produced)
                for key2,value in step.get('require_fields',{}).items():need(d[key2]==value,'CANONICAL_REQUIRED_FIELD:'+key2)
                if step.get('report_only_controls'):
                    rows=d['negative_report_data_controls'];need(len(rows)==5 and all(v['rejected'] is True and v['first_failed_witness'] for v in rows.values()),'REPORT_ONLY_CONTROLS')
                    extras[step['name']+'_report_only_control_evaluations']=5
            elif step['kind']=='controls':
                counts[step['name']]=controls(read(Path(resolve(step['receipt'],bindings))),step['contract'])
            elif step['kind']=='relocation':
                d=read(Path(resolve(step['receipt'],bindings)))
                need(d['status']=='PASS_REAL_RELOCATED_CODE_AND_FOUR_CAPSULE_EXECUTION' and d['report_byte_identity'] is True,'REAL_RELOCATION')
                need(d['report_sha256']==step['expected_report_sha256']==d['canonical_report_sha256'],'RELOCATION_REPORT_PIN')
                extras['actual_relocated_full_native_campaigns']=1
            else:raise RuntimeError('UNKNOWN_STEP_KIND')
        need(all(len(v)==2 and set(v)=={False,True} for v in modes.values()),'EACH_CANONICAL_BOTH_MODES')
        print(json.dumps({'completed_group':group['name'],'canonical_report_pairs':len(canonical),'real_rejected_subprocesses':sum(counts.values())}),flush=True)
        return {'group':group['name'],'canonical_sha256':canonical,'real_rejected_subprocesses':counts,'extra_checks':extras,'runtime_seconds':timings}
    started=time.monotonic()
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:results=list(pool.map(group_run,plan['groups']))
    need(authenticate(cap)==own and sources(repo,plan)==before,'ALL_FROZEN_SOURCES_UNCHANGED')
    need(all(regular(Path(path)) and sha(Path(path))==pin for path,pin in staged_original.items()),'STAGED_SCIENTIFIC_AND_WRAPPER_SOURCES_UNCHANGED')
    total=sum(sum(g['real_rejected_subprocesses'].values()) for g in results)
    reports=sum(len(g['canonical_sha256']) for g in results)
    need(total==plan['expected_rejected_subprocesses']==68 and reports==plan['expected_canonical_reports']==9,'COMPLETE_REPLAY_ROSTER')
    summary={'status':'PASS_COMPLETE_POSTCONTACT_SCREEN_CHECKPOINT','capsule_manifest_sha256':own,'source_capsule_manifest_sha256':before,
        'all_frozen_sources_authenticated_before_scientific_imports_and_after':True,'symbolic_versions':versions,
        'canonical_reports_reproduced_normal_and_optimized':reports,'real_negative_subprocesses':total,
        'report_only_controls_distinct':5,'report_only_controls_mode_evaluations':10,'focused_code_checks_per_mode':4016,
        'runtime_wrapper_patches':patches,'scientific_sources_modified':False,'runtime_path_global_overrides':plan['runtime_global_overrides'],
        'groups':results,'total_replay_seconds':round(time.monotonic()-started,6),'scope':plan['scope']}
    (out/'VERIFICATION.json').write_text(json.dumps(summary,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':summary['status'],'canonical_reports':reports,'real_negative_subprocesses':total,'runtime_seconds':summary['total_replay_seconds']}),flush=True)

if __name__=='__main__':main()
