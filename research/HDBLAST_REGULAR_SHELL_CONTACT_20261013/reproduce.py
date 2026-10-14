#!/usr/bin/env python3
"""Authenticate and replay the reviewed regular-shell-contact capsule.

Exact identities and controls support the reviewed proofs; they do not
mechanically certify analytic PDE existence. Prior science is not rerun.
--verify-only requires only the standard library and writes no files.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import stat
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
PRIOR_PINS = {
    'basis/PRIOR_MANIFEST.json':'a873132efe11755217ddb865fd0156041b8c1ad21ae1c45310183da4fb6140ae',
    'basis/SPECTRAL_MANIFEST.json':'d1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f',
}
STATUS = 'FROZEN_REVIEWED_REGULAR_SHELL_CONTACT_CAPSULE'


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def unique_object(pairs):
    result = {}
    for key,value in pairs:
        require(key not in result,'Duplicate JSON key: '+key)
        result[key] = value
    return result


def read_json(path):
    return json.loads(path.read_text(),object_pairs_hook=unique_object)


def rel(name):
    require(isinstance(name,str) and name and '\\' not in name,'Portable relative path required')
    part = PurePosixPath(name)
    require(not part.is_absolute() and '..' not in part.parts and part.as_posix() == name
            and name != '.','Unsafe/noncanonical path: '+name)
    return Path(*part.parts)


def snapshot(root):
    require(root.is_dir() and not root.is_symlink(),'Real capsule directory required')
    rows = {}
    for path in sorted(root.rglob('*')):
        mode = path.lstat().st_mode
        require(stat.S_ISDIR(mode) or stat.S_ISREG(mode),'Only regular capsule files/directories allowed')
        if stat.S_ISREG(mode):
            rows[path.relative_to(root).as_posix()] = sha(path)
    return rows


def inventories(data,item):
    kind = item['inventory_schema']
    if kind == 'detailed':
        require(data['passed'] is True and all(row['passed'] is True for row in data['checks']),
                'Detailed checker did not pass: '+item['name'])
        positive = [row['name'] for row in data['checks'] if 'wrong_residual' not in row]
        controls = [row['name'] for row in data['checks'] if 'wrong_residual' in row]
        require(all(row['residual'] == '0' for row in data['checks'] if 'wrong_residual' not in row),
                'Positive exact residual is not zero')
        require(all(row['nonzero_witness'] not in ('','0') for row in data['checks'] if 'wrong_residual' in row),
                'Wrong-formula control has no nonzero witness')
    elif kind == 'split_list':
        require(data['status'] == 'pass','Split checker did not pass: '+item['name'])
        positive = data['checks']
        controls = [row['name'] for row in data['wrong_controls']]
        require(data['check_count'] == len(positive)
                and all(row['nonzero_residual'] not in ('','0') for row in data['wrong_controls']),
                'Split checker count/control residual differs')
    else:
        raise RuntimeError('Unknown explicit inventory schema: '+kind)
    require(len(positive) == len(set(positive)) and len(controls) == len(set(controls))
            and not set(positive)&set(controls),'Duplicate/overlapping identity inventory')
    require(sorted(positive) == item['positive_names'] and sorted(controls) == item['control_names'],
            'Exact named inventory differs: '+item['name'])
    require((len(positive),len(controls)) == (item['positive_checks'],item['wrong_formula_controls']),
            'Identity/control counts differ')
    return len(positive),len(controls)


def contract(root):
    result = read_json(root/'CHECK_CONTRACTS.json')
    items = result['checks']
    require(bool(items) and len({item['name'] for item in items}) == len(items),'Explicit unique checker roster required')
    for item in items:
        require(re.fullmatch(r'[a-z][a-z0-9_]*',item['name']) is not None,'Portable checker name required')
        rel(item['script']);rel(item['expected_receipt'])
        require(sha(root/item['script']) == item['script_sha256'],'Checker source pin differs')
        require(item['positive_names'] == sorted(set(item['positive_names']))
                and item['control_names'] == sorted(set(item['control_names'])),'Sorted unique inventories required')
        require(bool(item['failure_controls']),'At least one meaningful failure control per checker required')
        inventories(read_json(root/item['expected_receipt']),item)
    require(result['positive_checks_per_mode'] == sum(item['positive_checks'] for item in items)
            and result['wrong_formula_controls_per_mode'] == sum(item['wrong_formula_controls'] for item in items),
            'Explicit roster totals differ')
    return result


def verify(root,pin=None):
    manifest_path = root/'MANIFEST.json'
    if pin is not None:
        require(re.fullmatch(r'[0-9a-f]{64}',pin) is not None and sha(manifest_path) == pin,
                'Manifest differs from external pin')
    manifest = read_json(manifest_path)
    require(manifest['status'] == STATUS,'Unexpected manifest status')
    before = snapshot(root)
    names = set()
    total = 0
    for row in manifest['files']:
        name = row['path'];rel(name)
        require(name not in names and name != 'MANIFEST.json','Duplicate/self manifest path')
        names.add(name)
        require(before.get(name) == row['sha256'] and (root/name).stat().st_size == row['bytes'],
                'Capsule file differs: '+name)
        total += row['bytes']
    require(names == set(before)-{'MANIFEST.json'},'Comprehensive manifest coverage differs')
    require(len(names) == manifest['file_count_excluding_this_manifest']
            and total == manifest['bytes_excluding_this_manifest'],'Manifest totals differ')
    entries = {}
    for name,expected in PRIOR_PINS.items():
        require(sha(root/name) == expected,'Pinned prior manifest differs: '+name)
        rows = read_json(root/name)['files']
        require(len({row['path'] for row in rows}) == len(rows),'Duplicate prior manifest rows')
        entries[name] = {row['path']:row for row in rows}
    source = read_json(root/'basis/SOURCE_INPUTS.json')
    require(source['prior_manifest_pins'] == PRIOR_PINS,'Prior input manifest registration differs')
    copied = set()
    for row in source['copied_files']:
        target = row['capsule_path'];rel(target);rel(row['source_path'])
        require(target.startswith('basis/prior-checkpoint/') and target not in copied,'Duplicate/outside prior copy')
        copied.add(target)
        original = entries[row['source_manifest']][row['source_path']]
        require(row['bytes'] == original['bytes'] and row['sha256'] == original['sha256']
                and sha(root/target) == row['sha256'] and (root/target).stat().st_size == row['bytes'],
                'Copied proof input differs: '+target)
    require(bool(copied) and copied == {'basis/prior-checkpoint/'+name
                                      for name in snapshot(root/'basis/prior-checkpoint')},
            'Copied proof-input coverage differs')
    require(bool(source['required_proof_paths'])
            and set(source['required_proof_paths']) <= copied,'Required prior proofs not copied/authenticated')
    contract(root)
    return before


def command(script,*arguments,optimized=False):
    return [sys.executable,'-I','-B']+(['-O'] if optimized else [])+[str(script),*map(str,arguments)]


def run(argv,log,fail=False):
    process = subprocess.run(argv,cwd=log.parent,capture_output=True,text=True,timeout=300,
                             env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
    log.write_text(process.stdout+process.stderr)
    require(process.returncode != 0 if fail else process.returncode == 0,'Unexpected exit; inspect '+str(log))
    return process


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--manifest-sha256')
    parser.add_argument('--verify-only',action='store_true')
    args = parser.parse_args()
    before = verify(ROOT,args.manifest_sha256)
    if args.verify_only:
        print('PASS_CAPSULE_AND_COPIED_PROOF_INTEGRITY');return
    require(args.output is not None,'Provide fresh outside-capsule --output')
    output = args.output.resolve()
    require(output != ROOT and ROOT not in output.parents,'Output must be outside capsule')
    require(not output.exists(),'Output exists; choose a fresh directory')
    import sympy,mpmath
    require(platform.python_version() == '3.12.14' and sympy.__version__ == '1.14.0'
            and mpmath.__version__ == '1.3.0','Use pinned Python/SymPy/mpmath versions')
    output.mkdir(parents=True)
    work = output/'relocated-capsule';shutil.copytree(ROOT,work)
    require(verify(work,args.manifest_sha256) == before,'Relocated capsule differs')
    roster = contract(work);results = [];failure_results = []
    for optimized in (False,True):
        mode = 'optimized' if optimized else 'normal'
        stage = output/mode;stage.mkdir()
        for item in roster['checks']:
            result = stage/(item['name']+'.json')
            run(command(work/item['script'],'--output',result,optimized=optimized),stage/(item['name']+'.log'))
            require(result.read_bytes() == (work/item['expected_receipt']).read_bytes(),
                    'Reviewed receipt bytes differ: '+item['name']+'/'+mode)
            positive,negative = inventories(read_json(result),item)
            results.append({'mode':mode,'name':item['name'],'positive_checks':positive,
                            'wrong_formula_controls':negative,'sha256':sha(result)})
            for control in item['failure_controls']:
                label = item['name']+'-'+control['name']+'-'+mode
                failed = output/('failed-'+label+'.json')
                process = run(command(work/item['script'],*control['arguments'],'--output',failed,
                                      optimized=optimized),output/(label+'.log'),True)
                if failed.exists():
                    data = read_json(failed)
                    require(data.get('passed') is False,'Failure control wrote a success/ambiguous report')
                    require(any(row['name'] == control['failed_check'] and row['passed'] is False
                                for row in data['checks']),'Expected coefficient mutation did not fail')
                else:
                    require(control['diagnostic'] in process.stdout+process.stderr,
                            'Expected failure reason missing')
                failure_results.append({'mode':mode,'checker':item['name'],'control':control['name'],
                                        'exit_nonzero':True,'success_receipt_written':False,
                                        'explicit_failed_report_written':failed.exists()})
    for item in roster['checks']:
        require((output/'normal'/(item['name']+'.json')).read_bytes()
                == (output/'optimized'/(item['name']+'.json')).read_bytes(),'Normal and optimized modes differ')
    corrupt = output/'corrupted-input-copy';shutil.copytree(work,corrupt)
    target = read_json(corrupt/'basis/SOURCE_INPUTS.json')['copied_files'][0]['capsule_path']
    with (corrupt/target).open('ab') as handle:handle.write(b'\nDeliberate proof-input corruption.\n')
    altered = output/'altered-manifest-copy';shutil.copytree(work,altered)
    original_pin = sha(altered/'MANIFEST.json')
    with (altered/'MANIFEST.json').open('ab') as handle:handle.write(b'\n')
    for optimized in (False,True):
        mode = 'optimized' if optimized else 'normal'
        process = run(command(corrupt/'reproduce.py','--verify-only',optimized=optimized),
                      output/('corrupted-input-'+mode+'.log'),True)
        require('Capsule file differs: '+target in process.stderr,'Proof corruption rejected for wrong reason')
        process = run(command(altered/'reproduce.py','--verify-only','--manifest-sha256',original_pin,
                              optimized=optimized),output/('altered-manifest-'+mode+'.log'),True)
        require('Manifest differs from external pin' in process.stderr,'Manifest mutation rejected for wrong reason')
    # Independently exercise each fixed prior pin after re-sealing only the
    # outer negative-test manifest. Neither authentic source copy is changed.
    for prior_name in PRIOR_PINS:
        altered_prior = output/('altered-prior-'+Path(prior_name).stem)
        shutil.copytree(work,altered_prior)
        with (altered_prior/prior_name).open('ab') as handle:handle.write(b'\n')
        manifest = read_json(altered_prior/'MANIFEST.json')
        for row in manifest['files']:
            if row['path'] == prior_name:
                row['bytes'] = (altered_prior/prior_name).stat().st_size
                row['sha256'] = sha(altered_prior/prior_name)
        manifest['bytes_excluding_this_manifest'] = sum(row['bytes'] for row in manifest['files'])
        (altered_prior/'MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
        for optimized in (False,True):
            mode = 'optimized' if optimized else 'normal'
            process = run(command(altered_prior/'reproduce.py','--verify-only',optimized=optimized),
                          output/('altered-prior-'+Path(prior_name).stem+'-'+mode+'.log'),True)
            require('Pinned prior manifest differs: '+prior_name in process.stderr,
                    'Prior manifest mutation rejected for wrong reason')
    require(snapshot(ROOT) == before and snapshot(work) == before,'Replay changed a source/relocated capsule')
    receipt = {'status':'PASS_REVIEWED_REGULAR_SHELL_CONTACT_REPLAY',
               'source_manifest_sha256':sha(ROOT/'MANIFEST.json'),
               'python':platform.python_version(),'sympy':sympy.__version__,'mpmath':mpmath.__version__,
               'positive_checks_per_mode':roster['positive_checks_per_mode'],
               'wrong_formula_controls_per_mode':roster['wrong_formula_controls_per_mode'],
               'normal_and_optimized_receipts_identical':True,'exact_named_inventories_verified':True,
               'proof_input_corruption_rejected_both_modes':True,
               'externally_pinned_manifest_corruption_rejected_both_modes':True,
               'both_fixed_prior_manifest_mutations_rejected_after_outer_reseal_both_modes':True,
               'checker_failures_rejected_without_success_receipts_both_modes':True,
               'source_and_relocated_capsules_unchanged':True,'prior_science_replayed':False,
               'checks':results,'failure_controls':failure_results,
               'scope':'Exact local identities, explicit wrong-formula controls, and provenance/replay integrity. Analytic shell-contact PDE existence, smoothness, energy estimates, physical pairing and theorem hypotheses remain in the reviewed mathematical proof; no machine PDE certificate or full-state-space claim.'}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({key:receipt[key] for key in ('status','positive_checks_per_mode','wrong_formula_controls_per_mode')},sort_keys=True))


if __name__ == '__main__':
    main()
