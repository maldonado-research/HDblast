#!/usr/bin/env python3
"""Authenticate externally pinned source/ZIP closure, then run a fresh replay.

No placeholder pin values are supplied. A separately reviewed concrete pin
file is mandatory; no scientific program runs until every pin and inventory
has been validated. Standard library only; no network or Git mutation.
"""
import argparse
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import subprocess
import sys
import zipfile

CAPSULE='HDBLAST_FIRST_GAUGE_WAVE_20261016'


def need(ok,message):
    if not ok:raise ValueError(message)


def sha(raw):return hashlib.sha256(raw).hexdigest()


def canonical(name):
    need(isinstance(name,str),'source path must be a string')
    p=PurePosixPath(name)
    need(bool(name) and name!='.' and '\\' not in name and not p.is_absolute() and '..' not in p.parts and p.as_posix()==name,'noncanonical source path')
    return p


def authenticate(root,pins):
    manifest_raw=(root/'MANIFEST.json').read_bytes()
    need(sha(manifest_raw)==pins['manifest_sha256'],'external manifest pin differs')
    manifest=json.loads(manifest_raw);listed=set()
    need(manifest.get('capsule')==CAPSULE,'source capsule identity differs')
    for row in manifest['files']:
        name=row['path'];canonical(name)
        need(name not in listed and name!='MANIFEST.json','duplicate/reserved source member')
        need(type(row['bytes']) is int and row['bytes']>=0 and re.fullmatch('[0-9a-f]{64}',row['sha256']) is not None,'invalid source byte/hash entry')
        p=root/name
        need(p.is_file() and not any(q.is_symlink() for q in (p,*p.parents)),'regular source file required')
        raw=p.read_bytes()
        need(len(raw)==row['bytes'] and sha(raw)==row['sha256'],'source member differs: '+name)
        listed.add(name)
    walked=list(root.rglob('*'))
    need(all(not p.is_symlink() and (p.is_file() or p.is_dir()) for p in walked),'symlink or special source inventory entry')
    actual={p.relative_to(root).as_posix() for p in walked if p.is_file()}
    need(actual==listed|{'MANIFEST.json'},'complete source inventory differs')
    need(len(actual)==pins['source_files_including_manifest'],'external source count differs')
    need(sha((root/'reproduce.py').read_bytes())==pins['runner_sha256'],'external replay runner pin differs')
    roster=json.loads((root/'REPLAY_ROSTER.json').read_text());labels=set()
    for job in roster['jobs']:
        need(isinstance(job.get('name'),str) and re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',job['name']) is not None,'safe replay name required')
        need(isinstance(job.get('arguments',[]),list) and all(isinstance(v,str) for v in job.get('arguments',[])),'replay arguments must be strings')
        if 'timeout_seconds' in job:
            need(type(job['timeout_seconds']) is int and job['timeout_seconds']>0,'positive job timeout required')
        script=job['script'];canonical(script)
        need(script in listed,'replay script outside frozen inventory')
        need(isinstance(job.get('modes',['normal','optimized']),list),'replay modes must be an explicit list')
        for mode in job.get('modes',['normal','optimized']):
            need(mode in ('normal','optimized'),'unknown replay mode')
            label=job['name']+'.'+mode
            need(label not in labels,'duplicate replay label');labels.add(label)
            golden=job['goldens'][mode];canonical(golden)
            need(golden in listed,'replay golden outside frozen inventory')
    need(len(labels)==pins['expected_replay_processes'],'external replay process count differs')
    return actual


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,default=Path(__file__).resolve().parents[3])
    ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    pins_raw=args.pins.read_bytes();pins=json.loads(pins_raw)
    need(pins.get('capsule')==CAPSULE,'pin-file capsule identity differs')
    for key in ('manifest_sha256','archive_sha256','runner_sha256'):
        value=pins.get(key)
        need(isinstance(value,str) and re.fullmatch('[0-9a-f]{64}',value) is not None and value!='0'*64,'concrete complete external pin required: '+key)
    for key in ('archive_bytes','source_files_including_manifest','expected_replay_processes'):
        need(type(pins.get(key)) is int and pins[key]>0,'positive concrete external count required: '+key)
    if 'replay_timeout_seconds' in pins:
        need(type(pins['replay_timeout_seconds']) is int and pins['replay_timeout_seconds']>0,'positive replay timeout required')
    repo=args.repo_root.resolve();source=repo/'research'/CAPSULE;archive=repo/'research'/(CAPSULE+'.zip')
    need(source.is_dir() and not any(p.is_symlink() for p in (source,*source.parents)),'real source directory required')
    members=authenticate(source,pins)
    archive_raw=archive.read_bytes()
    need(not archive.is_symlink() and len(archive_raw)==pins['archive_bytes'] and sha(archive_raw)==pins['archive_sha256'],'external archive identity differs')
    with zipfile.ZipFile(archive) as z:
        expected={CAPSULE+'/'+p for p in members}
        need(len(z.namelist())==len(expected) and set(z.namelist())==expected,'complete ZIP inventory differs')
        for info in z.infolist():
            canonical(info.filename)
            need(((info.external_attr>>16)&0o170000)!=0o120000,'symlink ZIP member')
            relative=PurePosixPath(info.filename).relative_to(CAPSULE).as_posix()
            need(z.read(info.filename)==(source/relative).read_bytes(),'ZIP/source member bytes differ')
    out=args.output_dir.resolve()
    need(out!=repo and repo not in out.parents and not out.exists(),'new output directory outside repository required')
    out.mkdir(parents=True,exist_ok=False)
    extracted=out/'source'/CAPSULE;extracted.mkdir(parents=True)
    with zipfile.ZipFile(archive) as z:
        for name in sorted(members):
            target=extracted/name;target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes(z.read(CAPSULE+'/'+name))
    authenticate(extracted,pins)
    command=[sys.executable,'-I','-S','-B',str(extracted/'reproduce.py'),'--manifest-sha256',pins['manifest_sha256'],'--output-dir',str(out/'replay')]
    result=subprocess.run(command,cwd=extracted,capture_output=True,timeout=pins.get('replay_timeout_seconds',900))
    (out/'reproduce.stdout').write_bytes(result.stdout);(out/'reproduce.stderr').write_bytes(result.stderr)
    need(result.returncode==0,'authenticated fresh replay failed; inspect external logs')
    replay_raw=(out/'replay/REPRODUCTION_RECORD.json').read_bytes();replay=json.loads(replay_raw)
    need(replay.get('status')=='PASS_COMPLETE_REPORT_AGREEMENT' and replay.get('source_unchanged') is True,'fresh replay did not report complete agreement and source preservation')
    need(replay.get('manifest_sha256')==pins['manifest_sha256'] and replay.get('processes')==pins['expected_replay_processes'],'fresh replay identity/process count differs')
    authenticate(source,pins);authenticate(extracted,pins)
    need(sha(archive.read_bytes())==pins['archive_sha256'],'original ZIP changed during replay')
    report={'status':'PASS_PINNED_SOURCE_ZIP_AND_FRESH_REPLAY','capsule':CAPSULE,
            'publication_pins_sha256':sha(pins_raw),'verifier_sha256':sha(Path(__file__).read_bytes()),
            'manifest_sha256':pins['manifest_sha256'],'archive_sha256':pins['archive_sha256'],
            'runner_sha256':pins['runner_sha256'],'source_files_including_manifest':len(members),
            'matching_ZIP_members':len(members),'fresh_replay_processes':replay['processes'],
            'fresh_reproduction_record_sha256':sha(replay_raw),'original_and_extracted_source_unchanged':True,
            'scope':'Byte/source closure and the explicitly recorded scientific roster; no additional mathematical theorem or external peer-review claim.'}
    (out/'VERIFICATION_RECORD.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
