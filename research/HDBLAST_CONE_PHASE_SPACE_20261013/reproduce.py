#!/usr/bin/env python3
"""Authenticate and replay the cone-current research capsule in two Python modes.

The manifest is bound externally by the immutable Git revision and CI pin.
--manifest-sha256 optionally enforces that external pin before any replay.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
PRIOR_MANIFEST_SHA='d1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(root):
    entries={}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(),'Symbolic links are not part of this capsule')
        if path.is_file():entries[path.relative_to(root).as_posix()]=sha(path)
    return entries


def verify(root,manifest_pin=None):
    manifest_path=root/'MANIFEST.json'
    if manifest_pin:
        require(sha(manifest_path)==manifest_pin,'Manifest differs from external pin')
    manifest=json.loads(manifest_path.read_text())
    require(manifest['status']=='FROZEN_REVIEWED_CONE_CURRENT_CAPSULE','Unexpected manifest status')
    actual=snapshot(root);names=set();total=0
    for row in manifest['files']:
        name=row['path'];part=Path(name)
        require(not part.is_absolute() and '..' not in part.parts and name not in names,'Unsafe/duplicate manifest path')
        names.add(name);path=root/part
        require(actual.get(name)==row['sha256'] and path.stat().st_size==row['bytes'],'Capsule file differs: '+name)
        total+=row['bytes']
    require(names==set(actual)-{'MANIFEST.json'},'Manifest coverage differs')
    require(len(names)==manifest['file_count_excluding_this_manifest'] and total==manifest['bytes_excluding_this_manifest'],'Manifest totals differ')
    old=root/'basis'/'PRIOR_MANIFEST.json'
    require(sha(old)==PRIOR_MANIFEST_SHA,'Prior manifest differs')
    source=json.loads((root/'basis'/'SOURCE_INPUTS.json').read_text())
    old_entries={x['path']:x for x in json.loads(old.read_text())['files']}
    for row in source['copied_files']:
        require(row==old_entries[row['path']],'Prior input declaration differs')
        path=root/'basis'/'prior-checkpoint'/row['path']
        require(sha(path)==row['sha256'] and path.stat().st_size==row['bytes'],'Prior input differs')
    return actual


def run(argv,log,expected_failure=False):
    process=subprocess.run(argv,cwd=log.parent,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True)
    log.write_text(process.stdout+process.stderr)
    require((process.returncode!=0) if expected_failure else (process.returncode==0),'Unexpected exit; inspect '+str(log))
    return process


def command(script,*arguments,optimized=False):
    return [sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(script),*map(str,arguments)]


def counts(name,data):
    if name=='covariant':return data['identities'],data['wrong_formula_controls']
    if name=='domain':
        return sum(x['result']=='pass' for x in data['checks']),sum(x['result']=='rejected as intended' for x in data['checks'])
    if name=='tensor':return data['positive_count'],data['wrong_formula_count']
    if name=='independent':return data['positive_check_count'],data['wrong_formula_control_count']
    raise RuntimeError('Unknown check contract')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args()
    before=verify(ROOT,args.manifest_sha256)
    if args.verify_only:
        print('PASS_CAPSULE_AND_PRIOR_INPUT_INTEGRITY');return
    require(args.output is not None,'Provide a fresh --output directory')
    target=args.output.resolve()
    require(target!=ROOT and ROOT not in target.parents,'Output must be outside capsule')
    require(not target.exists(),'Output exists; choose a fresh directory')
    import sympy,mpmath
    require(platform.python_version()=='3.12.14','Recorded replay uses Python 3.12.14')
    require(sympy.__version__=='1.14.0' and mpmath.__version__=='1.3.0','Install the recorded requirements')
    target.mkdir(parents=True)
    work=target/'relocated-capsule';shutil.copytree(ROOT,work)
    require(verify(work,args.manifest_sha256)==before,'Relocated capsule differs')
    contract=json.loads((work/'CHECK_CONTRACTS.json').read_text())
    results=[]
    for mode in ('normal','optimized'):
        stage=target/mode;stage.mkdir()
        for item in contract['checks']:
            out=stage/(item['name']+'.json')
            run(command(work/item['script'],'--output',out,optimized=mode=='optimized'),stage/(item['name']+'.log'))
            require(out.read_bytes()==(work/item['expected_receipt']).read_bytes(),'Reviewed receipt differs: '+item['name']+'/'+mode)
            actual_counts=counts(item['name'],json.loads(out.read_text()))
            require(actual_counts==(item['positive_checks'],item['wrong_formula_controls']),'Check count differs')
            results.append({'mode':mode,'name':item['name'],'positive_checks':actual_counts[0],
                            'wrong_formula_controls':actual_counts[1],'sha256':sha(out)})
    for item in contract['checks']:
        name=item['name']+'.json'
        require((target/'normal'/name).read_bytes()==(target/'optimized'/name).read_bytes(),'Modes differ')
    # The copied prior theorem inputs must remain unchanged, even though these
    # scripts use them only to document conventions and conditional scope.
    corrupt=target/'corrupted-copy';shutil.copytree(work,corrupt)
    input_file=corrupt/'basis'/'prior-checkpoint'/'RESULTS.md'
    input_file.write_bytes(input_file.read_bytes()+b'\nDeliberate input-integrity control.\n')
    for optimized in (False,True):
        mode='optimized' if optimized else 'normal'
        run(command(corrupt/'reproduce.py','--verify-only',optimized=optimized),target/('corrupted-input-'+mode+'.log'),True)
    changed_manifest=target/'altered-manifest-copy';shutil.copytree(work,changed_manifest)
    original_pin=sha(changed_manifest/'MANIFEST.json')
    with (changed_manifest/'MANIFEST.json').open('ab') as f:f.write(b'\n')
    run(command(changed_manifest/'reproduce.py','--verify-only','--manifest-sha256',original_pin),target/'altered-manifest.log',True)
    for optimized in (False,True):
        mode='optimized' if optimized else 'normal';out=target/('MUST_NOT_EXIST_'+mode+'.json')
        run(command(work/'tensor'/'verify_tensor_cone_domain.py','--inject-failure','--output',out,optimized=optimized),target/('injected-failure-'+mode+'.log'),True)
        require(not out.exists(),'Injected failure wrote a success receipt')
    require(snapshot(ROOT)==before and snapshot(work)==before,'Replay changed a capsule')
    receipt={
        'status':'PASS_REVIEWED_CONE_CURRENT_REPLAY',
        'scope':'Exact local algebra, prototype scaling, source integrity and wrong-formula controls. Analytic domain and asymptotic hypotheses remain in the reviewed notes; no global registered radiative-solution existence or phase-space completeness proof.',
        'python':platform.python_version(),'sympy':sympy.__version__,'mpmath':mpmath.__version__,
        'source_manifest_sha256':sha(ROOT/'MANIFEST.json'),
        'positive_checks_per_mode':sum(x['positive_checks'] for x in contract['checks']),
        'wrong_formula_controls_per_mode':sum(x['wrong_formula_controls'] for x in contract['checks']),
        'normal_and_optimized_receipts_identical':True,
        'corrupted_prior_input_rejected_both_modes':True,
        'externally_pinned_altered_manifest_rejected':True,
        'injected_failure_rejected_without_success_receipt_both_modes':True,
        'source_and_relocated_capsules_unchanged':True,
        'prior_scientific_capsule_replayed':False,
        'checks':results,
    }
    (target/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':receipt['status'],'positive_checks_per_mode':receipt['positive_checks_per_mode'],
                      'wrong_formula_controls_per_mode':receipt['wrong_formula_controls_per_mode']},sort_keys=True))


if __name__=='__main__':main()
