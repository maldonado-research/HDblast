#!/usr/bin/env python3
"""Deterministically package a completed metric checkpoint and its freeze receipt.

This tool never calls a physical source, mode producer, quadrature, or replay.
FULL_REGISTRATION protects prospective inputs; MANIFEST protects the complete
curated package including its separately added post-freeze receipt. Final ZIP
replay/publication receipts belong outside this checkpoint to avoid cycles.
"""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import zipfile

CHECKPOINT_NAME='HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_FOLLOWUP'
REQUIRED_RECEIPT_FIELDS=('public_freeze_commit','registration_sha256','independent_manifest_sha256')


def sha(data):return hashlib.sha256(data).hexdigest()


def require(condition,message):
    if not condition:raise RuntimeError(message)


def safe_file(root,name):
    relative=Path(name)
    require(not relative.is_absolute() and '..' not in relative.parts,'Unsafe payload path: '+name)
    target=root/relative
    require(target.resolve().is_relative_to(root),'Payload path escapes checkpoint: '+name)
    require(target.is_file() and not target.is_symlink(),'Missing or symlink file: '+name)
    for parent in target.parents:
        if parent==root:break
        require(not parent.is_symlink(),'Symlink directory: '+name)
    return target


def verify_inputs(root):
    registration_bytes=safe_file(root,'FULL_REGISTRATION.json').read_bytes()
    registration=json.loads(registration_bytes)
    receipt_bytes=safe_file(root,'FREEZE_RECEIPT.json').read_bytes()
    receipt=json.loads(receipt_bytes)
    require(all(field in receipt for field in REQUIRED_RECEIPT_FIELDS),'Incomplete post-freeze receipt')
    commit=receipt['public_freeze_commit']
    require(isinstance(commit,str) and len(commit)==40 and all(c in '0123456789abcdef' for c in commit),'Invalid public freeze commit')
    require(receipt['registration_sha256']==sha(registration_bytes),'Receipt registration hash mismatch')
    require(isinstance(registration.get('files'),dict) and registration['files'],'Missing registered file table')
    for name,digest in registration['files'].items():
        require(isinstance(digest,str) and len(digest)==64 and all(c in '0123456789abcdef' for c in digest),'Invalid registered hash: '+name)
        require(sha(safe_file(root,name).read_bytes())==digest,'Frozen input changed: '+name)
    independent=safe_file(root,'independent/MANIFEST.json').read_bytes()
    require(sha(independent)==receipt['independent_manifest_sha256']==registration['independent_manifest_sha256'],'Independent manifest pin mismatch')
    require(registration['files'].get('independent/MANIFEST.json')==sha(independent),'Independent manifest omitted from registration')
    items=json.loads(independent)['files']
    if isinstance(items,dict):items=[{'path':name,'sha256':digest} for name,digest in items.items()]
    require(isinstance(items,list) and items,'Missing independent source manifest')
    seen=set()
    for item in items:
        relative=item['path'];name='independent/'+relative
        require(relative not in seen,'Duplicate independent source path')
        seen.add(relative)
        require(sha(safe_file(root,name).read_bytes())==item['sha256'],'Independent frozen input changed: '+relative)
        require(registration['files'].get(name)==item['sha256'],'Independent input differs from full registration: '+relative)
    return {'public_freeze_commit':commit,'registration_sha256':sha(registration_bytes),
            'independent_manifest_sha256':sha(independent),'freeze_receipt_sha256':sha(receipt_bytes)}


def build(root):
    root=root.resolve()
    require(root.name==CHECKPOINT_NAME,'Unexpected checkpoint name')
    inputs=verify_inputs(root)
    archive_name=CHECKPOINT_NAME+'.zip'
    excluded={'MANIFEST.json',archive_name,archive_name+'.sha256'}
    files={}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(),'Symlink forbidden: '+str(path))
        if path.is_file():
            name=path.relative_to(root).as_posix()
            require('__pycache__' not in path.parts and path.suffix!='.pyc','Cache forbidden: '+name)
            if name not in excluded:files[name]=sha(safe_file(root,name).read_bytes())
    require('FREEZE_RECEIPT.json' in files and 'FREEZE_RECEIPT.json' not in excluded,'Package must protect freeze receipt')
    manifest={'schema_version':1,'files':files,'excluded':sorted(excluded),
              'freeze_receipt_sha256':inputs['freeze_receipt_sha256'],
              'scope':'Exact complete curated metric payload, including post-freeze receipt; only package bookkeeping excluded.'}
    manifest_bytes=(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    (root/'MANIFEST.json').write_bytes(manifest_bytes)
    names=sorted([*files,'MANIFEST.json'])
    def construct():
        stream=io.BytesIO()
        with zipfile.ZipFile(stream,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
            for name in names:
                item=zipfile.ZipInfo(CHECKPOINT_NAME+'/'+name,date_time=(2026,10,2,0,0,0))
                item.compress_type=zipfile.ZIP_DEFLATED
                item.create_system=3
                item.external_attr=0o100644<<16
                archive.writestr(item,safe_file(root,name).read_bytes(),compresslevel=9)
        return stream.getvalue()
    raw=construct()
    require(raw==construct(),'Repeated deterministic archive construction differed')
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        require(archive.namelist()==[CHECKPOINT_NAME+'/'+name for name in names],'Archive membership mismatch')
        require(archive.testzip() is None,'Archive CRC check failed')
        for name in names:
            require(archive.read(CHECKPOINT_NAME+'/'+name)==safe_file(root,name).read_bytes(),'Archive source bytes differ: '+name)
    require(verify_inputs(root)==inputs,'Frozen inputs/receipt changed during packaging')
    current={p.relative_to(root).as_posix():sha(safe_file(root,p.relative_to(root).as_posix()).read_bytes())
             for p in sorted(root.rglob('*')) if p.is_file() and p.relative_to(root).as_posix() not in excluded}
    require(current==files,'Payload membership/bytes changed during packaging')
    require((root/'MANIFEST.json').read_bytes()==manifest_bytes,'Package manifest changed during packaging')
    (root/archive_name).write_bytes(raw)
    (root/(archive_name+'.sha256')).write_text(sha(raw)+'  '+archive_name+'\n')
    return {'status':'PASS','payload_files':len(files),'zip_members':len(names),'package_files':len(files)+3,
            'zip_bytes':len(raw),'zip_sha256':sha(raw),'repeat_build_identical':True,
            'all_members_equal_source_bytes':True,'archive_crc_verified':True,
            'freeze_receipt_included_and_hashed':True,'frozen_inputs_unchanged':True,**inputs}


def main():
    print(json.dumps(build(Path(__file__).resolve().parents[1]),sort_keys=True))

if __name__=='__main__':main()
