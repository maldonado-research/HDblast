#!/usr/bin/env python3
"""Synthetic package guards; zero physical sources or child producers."""
import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import zipfile

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('metric_package_builder',HERE/'build_package.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p,d):p.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')


def fixture(parent):
    r=parent/B.CHECKPOINT_NAME;r.mkdir(parents=True)
    (r/'code').mkdir();(r/'independent').mkdir()
    shutil.copyfile(HERE/'build_package.py',r/'code/build_package.py')
    (r/'independent/forced_metric.py').write_text('# Synthetic inert fixture; never executed.\n')
    independent={'files':[{'path':'forced_metric.py','sha256':sha(r/'independent/forced_metric.py')}]}
    write(r/'independent/MANIFEST.json',independent)
    runtime={'python_major_minor':[3,12],'physical_producers_optimized':False,'numpy':'2.2.6','scipy':'1.15.3','sympy':'1.14.0','mpmath':'1.3.0','matplotlib':'3.10.1'}
    experiment={'gates':{'synthetic_guard':1},'runtime':runtime}
    write(r/'EXPERIMENT.json',experiment)
    reg={'schema_version':1,'independent_manifest_sha256':sha(r/'independent/MANIFEST.json'),'frozen_configuration':experiment,'frozen_gates':experiment['gates'],'files':{p.relative_to(r).as_posix():sha(p) for p in sorted(r.rglob('*')) if p.is_file()}}
    write(r/'FULL_REGISTRATION.json',reg)
    write(r/'FREEZE_RECEIPT.json',{'public_freeze_commit':'1'*40,'registration_sha256':sha(r/'FULL_REGISTRATION.json'),'independent_manifest_sha256':sha(r/'independent/MANIFEST.json')})
    (r/'README.md').write_text('Synthetic package fixture only; not scientific results.\n')
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--replay-driver',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise RuntimeError('Refusing overwrite')
    args.output.mkdir(parents=True)
    checks=[]
    def record(name,condition):
        if not condition:raise RuntimeError('Package guard failed: '+name)
        checks.append(name)
    r=fixture(args.output/'determinism')
    original={p.relative_to(r).as_posix():sha(p) for p in r.rglob('*') if p.is_file()}
    first=B.build(r);z=r/(r.name+'.zip');before=z.read_bytes();second=B.build(r)
    record('two full builds byte-identical',before==z.read_bytes() and first==second)
    manifest=json.loads((r/'MANIFEST.json').read_text())
    record('freeze receipt included and hashed','FREEZE_RECEIPT.json' in manifest['files'] and manifest['files']['FREEZE_RECEIPT.json']==sha(r/'FREEZE_RECEIPT.json'))
    record('only three package bookkeeping paths excluded',set(manifest['excluded'])=={'MANIFEST.json',r.name+'.zip',r.name+'.zip.sha256'})
    record('all original inputs and receipt bytes unchanged',all(sha(r/name)==pin for name,pin in original.items()))
    with zipfile.ZipFile(z) as archive:
        record('exact membership and every source byte retained',set(archive.namelist())=={r.name+'/'+name for name in [*manifest['files'],'MANIFEST.json']} and all(archive.read(r.name+'/'+name)==(r/name).read_bytes() for name in [*manifest['files'],'MANIFEST.json']))
    replay_spec=importlib.util.spec_from_file_location('metric_replay_input_verifier',args.replay_driver)
    replay=importlib.util.module_from_spec(replay_spec);replay_spec.loader.exec_module(replay)
    verified=replay.verify_inputs(r)
    record('portable replay accepts authenticated synthetic package',verified['package_manifest_sha256']==sha(r/'MANIFEST.json') and verified['freeze_receipt_sha256']==sha(r/'FREEZE_RECEIPT.json'))
    def reject(name,mutation):
        root=fixture(args.output/name);mutation(root)
        try:B.build(root)
        except (RuntimeError,KeyError,ValueError,TypeError):checks.append(name+' rejected')
        else:raise RuntimeError('Package guard survived: '+name)
    reject('missing_freeze_receipt',lambda root:(root/'FREEZE_RECEIPT.json').unlink())
    reject('mutated_registered_input',lambda root:(root/'independent/forced_metric.py').write_text('changed'))
    reject('unregistered_symlink_file',lambda root:(root/'symlink').symlink_to(root/'README.md'))
    reject('unregistered_symlink_directory',lambda root:(root/'symdir').symlink_to(root/'independent',target_is_directory=True))
    reject('bytecode_cache',lambda root:(root/'cache.pyc').write_bytes(b'cache'))
    def receipt_mutation(root,field,value):
        d=json.loads((root/'FREEZE_RECEIPT.json').read_text());d[field]=value;write(root/'FREEZE_RECEIPT.json',d)
    reject('wrong_registration_receipt_hash',lambda root:receipt_mutation(root,'registration_sha256','0'*64))
    reject('wrong_independent_receipt_hash',lambda root:receipt_mutation(root,'independent_manifest_sha256','0'*64))
    reject('invalid_freeze_commit',lambda root:receipt_mutation(root,'public_freeze_commit','main'))
    def unsafe_registration(root):
        d=json.loads((root/'FULL_REGISTRATION.json').read_text());d['files']['../escape.txt']='0'*64;write(root/'FULL_REGISTRATION.json',d);receipt_mutation(root,'registration_sha256',sha(root/'FULL_REGISTRATION.json'))
    reject('registered_path_escape',unsafe_registration)
    def omitted_independent(root):
        d=json.loads((root/'FULL_REGISTRATION.json').read_text());d['files'].pop('independent/forced_metric.py');write(root/'FULL_REGISTRATION.json',d);receipt_mutation(root,'registration_sha256',sha(root/'FULL_REGISTRATION.json'))
    reject('independent_input_omitted_from_registration',omitted_independent)
    report={'status':'PASS','check_count':len(checks),'checks':checks,'physical_evaluations':0,'child_producers_executed':0,'python_optimization':sys.flags.optimize,'builder_sha256':sha(HERE/'build_package.py'),'test_source_sha256':sha(Path(__file__)),'replay_driver_sha256':sha(args.replay_driver),'scope':'Synthetic inert files only; deterministic archive, membership, receipt binding, no input mutation, and rejection guards.'}
    write(args.output/'CHECKS.json',report);print(json.dumps({k:report[k] for k in ('status','check_count','physical_evaluations','child_producers_executed')}))

if __name__=='__main__':main()
