"""Authenticate source capsule and every selected inherited member before imports."""
from pathlib import Path
import hashlib
import importlib.util
import json

FREEZE_SHA256 = '4ab1875803ba4bf3aa09881be3bd645b0355bb721fb08d6d469e75dffa7fe1d5'

def need(value, why):
    if not value:
        raise RuntimeError(why)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def json_read(path):
    return json.loads(path.read_bytes())

def manifest(cap, pin=None, closed=False):
    path=cap/'MANIFEST.json'
    if pin is not None:
        need(digest(path)==pin,'manifest pin '+str(path))
    seen={'MANIFEST.json'}
    for row in json_read(path)['files']:
        rel=Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,'safe unique member')
        seen.add(str(rel)); member=cap/rel
        need(member.is_file() and not member.is_symlink(),'regular manifest member '+str(member))
        need(member.stat().st_size==row['bytes'] and digest(member)==row['sha256'],'member digest '+str(member))
    if closed:
        need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen) for p in cap.rglob('*')),'closed source capsule')
    return {'manifest_sha256':digest(path),'member_count':len(seen)-1}

def authenticate(root, freeze_path):
    # Only the standard library has been imported at this point.
    f=json_read(freeze_path)
    need(digest(freeze_path)==FREEZE_SHA256,'refinement freeze pin')
    cap=root/f['capsule_relative']
    receipt={'capsule':manifest(cap,f['capsule_manifest_sha256'],True)}
    for rel,key in [('reference/INITIAL.json','initial_report_sha256'),('checks/producer.py','producer_sha256'),('checks/FREEZE.json','initial_freeze_sha256')]:
        need(digest(cap/rel)==f[key],'selected initial source '+rel)
    original=json_read(cap/'checks/FREEZE.json')
    selected={}
    for rel,pin in original['input_sha256'].items():
        path=root/rel
        need(path.is_file() and not path.is_symlink() and digest(path)==pin,'selected original source '+rel)
        selected[rel]=pin
        if rel.endswith('/MANIFEST.json'):
            manifest(path.parent,pin)
    profile=root/'research/HDBLAST_FULL_ANGULAR_RESPONSE_20261026/initial/results/INITIAL_PROFILES.json'
    need(digest(profile)==f['unchanged_velocity']['profile_sha256'],'velocity profile pin')
    receipt['selected_original']=selected
    receipt['profile_sha256']=digest(profile)
    receipt['closed_source_authenticated_before_import']=True
    return f,cap,receipt

def load_producer(root, freeze_path):
    f,cap,receipt=authenticate(root,freeze_path)
    spec=importlib.util.spec_from_file_location('authenticated_exact_packet_producer',cap/'checks/producer.py')
    p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
    p.authenticate(root,json_read(cap/'checks/FREEZE.json'))
    return f,cap,receipt,p
